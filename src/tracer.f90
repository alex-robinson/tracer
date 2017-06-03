
module tracer 

    use tracer_precision
    use tracer_interp 
    use bspline_module, only : bspline_3d 
    use ncio   
    use nml 

    implicit none 

    real(prec), parameter :: z_scale_in  = 1e0
    real(prec), parameter :: z_scale_out = 1.0/z_scale_in

    type tracer_par_class 
        integer :: n, n_active, n_max_dep, id_max 
        logical :: is_sigma                     ! Is the defined z-axis in sigma coords
        real(prec_time) :: time_now, time_old
        real(prec_time) :: time_dep, time_write 
        real(prec_time) :: dt
        real(prec) :: thk_min                   ! Minimum thickness of tracer (m)
        real(prec) :: H_min                     ! Minimum ice thickness to track (m)
        real(prec) :: depth_max                 ! Maximum depth of tracer (fraction)
        real(prec) :: U_max                     ! Maximum horizontal velocity of tracer to track (m/a)
        real(prec) :: H_min_dep                 ! Minimum ice thickness for tracer deposition (m)
        real(prec) :: alpha                     ! Slope of probability function
        character(len=56) :: weight             ! Weighting function for generating prob. distribution
        logical    :: noise                     ! Add noise to gridded deposition location
        real(prec) :: dens_z_lim                ! Distance from surface to count density
        integer    :: dens_max                  ! Max allowed density of particles at surface
        
        character(len=56) :: interp_method  
    end type 

    type tracer_state_class 
        integer, allocatable :: active(:), id(:)
        real(prec), allocatable :: x(:), y(:), z(:), sigma(:)
        real(prec), allocatable :: ux(:), uy(:), uz(:)
        real(prec), allocatable :: ax(:), ay(:), az(:)
        real(prec), allocatable :: dpth(:), z_srf(:)
        real(prec), allocatable :: thk(:)            ! Tracer thickness (for compression)
        real(prec), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
        real(prec), allocatable :: H(:)

    end type 

    type tracer_stats_class
        ! All stats variable at precision of writing (prec_wrt), since
        ! this should not need high precision output 

        real(prec_wrt), allocatable :: x(:), y(:)
        real(prec_wrt), allocatable :: depth_norm(:)
        real(prec_wrt), allocatable :: age_iso(:) 

        real(prec_wrt), allocatable :: depth_iso(:,:,:)
        real(prec_wrt), allocatable :: depth_iso_err(:,:,:)
        real(prec_wrt), allocatable :: dep_z_iso(:,:,:)
        integer,    allocatable :: density_iso(:,:,:)
        
        real(prec_wrt), allocatable :: ice_age(:,:,:)
        real(prec_wrt), allocatable :: ice_age_err(:,:,:)
        integer,    allocatable :: density(:,:,:)
           
    end type

    type tracer_dep_class 
        ! Standard deposition information (time and place)
        real(prec), allocatable :: time(:) 
        real(prec), allocatable :: H(:) 
        real(prec), allocatable :: x(:), y(:), z(:)
        real(prec), allocatable :: lon(:), lat(:) 

        ! Additional tracer deposition information (climate, isotopes, etc)
        real(prec), allocatable :: t2m(:,:)     ! 4 seasons 
        real(prec), allocatable :: pr(:,:)      ! 4 seasons 
        real(prec), allocatable :: t2m_prann(:) ! Precip-weighted temp

!         real(prec), allocatable :: d18O(:)
!         real(prec), allocatable :: dD(:)

    end type 

    type tracer_class 
        type(tracer_par_class)   :: par 
        type(tracer_state_class) :: now 
        type(tracer_dep_class)   :: dep
        type(tracer_stats_class) :: stats 

    end type 

    type(lin3_interp_par_type) :: par_lin 
    type(bspline_3d)           :: bspline3d_ux, bspline3d_uy, bspline3d_uz 

    private 

    ! For other tracer modules
    public :: tracer_par_class
    public :: tracer_state_class
    public :: tracer_dep_class
    public :: tracer_stats_class

    public :: tracer_reshape1D_vec
    public :: tracer_reshape2D_field 
    public :: tracer_reshape3D_field

    ! General public 
    public :: tracer_class 
    public :: tracer_init 
    public :: tracer_update 
    public :: tracer_end 
    public :: tracer_write_init, tracer_write 
    public :: tracer_write_stats 

    ! Conversion constants
    public :: z_scale_in 
    public :: z_scale_out 

contains 

    subroutine tracer_init(trc,filename,time,x,y,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(prec), intent(IN) :: x(:), y(:)
        logical,    intent(IN) :: is_sigma  
        real(prec_time), intent(IN) :: time 

        ! Local variables 
        integer :: i 

        ! Load the parameters
        call tracer_par_load(trc%par,filename,is_sigma)

        ! Allocate the state variables 
        call tracer_allocate(trc%now,trc%dep,n=trc%par%n)
        call tracer_allocate_stats(trc%stats,x,y)

        ! ===== Initialize stats depth axes ===============

        trc%stats%age_iso = [11.7,29.0,57.0,115.0,130.0]

        do i = 1, size(trc%stats%depth_norm)
            trc%stats%depth_norm(i) = 0.04*real(i)
        end do 

        ! =================================================

        ! Initialize state 
        trc%now%active    = 0 

        trc%now%id        = mv 
        trc%now%x         = mv 
        trc%now%y         = mv 
        trc%now%z         = mv 
        trc%now%sigma     = mv 
        trc%now%z_srf     = mv 
        trc%now%dpth      = mv 
        trc%now%ux        = mv 
        trc%now%uy        = mv 
        trc%now%uz        = mv 
        trc%now%ax        = mv 
        trc%now%ay        = mv 
        trc%now%az        = mv 
        trc%now%thk       = mv 
        trc%now%T         = mv 
        trc%now%H         = mv 

        trc%dep%time      = mv 
        trc%dep%H         = mv 
        trc%dep%x         = mv 
        trc%dep%y         = mv 
        trc%dep%z         = mv 

        trc%par%id_max    = 0 

        ! Initialize the time (to one older than now)
        trc%par%time_now   = time - 1000.0_dp
        trc%par%time_dep   = time - 1000.0 
        trc%par%time_write = time - 1000.0 

        ! Initialize random number generator 
        call random_seed() 

        return 

    end subroutine tracer_init

    subroutine tracer_update(trc,time,x,y,z,z_srf,H,ux,uy,uz,dep_now,stats_now,order,sigma_srf)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        real(prec), intent(IN) :: z_srf(:,:), H(:,:)
        real(prec), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
        logical, intent(IN) :: dep_now, stats_now  
        character(len=*), intent(IN), optional :: order 
        real(prec), intent(IN), optional :: sigma_srf     ! Value at surface by default (1 or 0?)

        ! Local variables  
        character(len=3) :: idx_order 
        integer    :: i, j, k, nx, ny, nz
        logical    :: rev_z 
        real(prec), allocatable :: x1(:), y1(:), z1(:)
        real(prec), allocatable :: zc(:)   ! Actual cartesian z-axis after applying sigma*H 
        real(prec), allocatable :: z_srf1(:,:), H1(:,:)
        real(prec), allocatable :: ux1(:,:,:), uy1(:,:,:), uz1(:,:,:)
        real(prec), allocatable :: usig1(:,:,:)
        real(prec) :: ux0, uy0, uz0 
        real(prec) :: dt 

        ! Update current time and time step
        trc%par%time_old = trc%par%time_now 
        trc%par%time_now = time 
        trc%par%dt       = dble(trc%par%time_now) - dble(trc%par%time_old)

        ! Update record of last deposition time if dep_now
        if (dep_now) trc%par%time_dep = trc%par%time_now 

        ! Determine order of indices (default ijk)
        idx_order = "ijk"
        if (present(order)) idx_order = trim(order)

        ! Allocate helper z-axis variable 
        if (allocated(zc)) deallocate(zc)
        allocate(zc(size(z)))

        zc = z 

        if (trc%par%is_sigma) then 
            ! Ensure z-axis is properly bounded
            where (abs(zc) .lt. 1e-5) zc = 0.0 

            if (minval(zc) .lt. 0.0 .or. maxval(zc) .gt. 1.0) then 
                write(0,*) "tracer:: error: sigma axis not bounded between zero and one."
                write(0,*) "z = ", zc 
                stop 
            end if

        end if 

        ! Note: GRISLI (nx,ny,nz): sigma goes from 1 to 0, so sigma(1)=1 [surface], sigma(nz)=0 [base]
        !       SICO (nz,ny,nx): sigma(1) = 0, sigma(nz) = 1
        ! reshape routines ensure ascending z-axis (nx,ny,nz) with sigma(nz)=1 [surface]
        
        ! Correct the sigma values if necessary,
        ! so that sigma==0 [base]; sigma==1 [surface]
        if (trc%par%is_sigma .and. present(sigma_srf)) then 
            if (sigma_srf .eq. 0.0) then 
                ! Adjust sigma values 
                zc = 1.0 - zc 
            end if 
        end if 

        ! Also determine whether z-axis is initially ascending or descending 
        rev_z = (zc(1) .gt. zc(size(zc)))

        call tracer_reshape1D_vec(x, x1,rev=.FALSE.)
        call tracer_reshape1D_vec(y, y1,rev=.FALSE.)
        call tracer_reshape1D_vec(real(zc,kind=prec),z1,rev=rev_z)
        call tracer_reshape2D_field(idx_order,z_srf,z_srf1)
        call tracer_reshape2D_field(idx_order,H,H1)
        call tracer_reshape3D_field(idx_order,ux,ux1,rev_z=rev_z)
        call tracer_reshape3D_field(idx_order,uy,uy1,rev_z=rev_z)
        call tracer_reshape3D_field(idx_order,uz,uz1,rev_z=rev_z)
        
        ! Get axis sizes (if ny==5, this is a 2D profile domain) 
        nx = size(x1,1)
        ny = size(y1,1)
        nz = size(z1,1)

        ! Scale z-axis as desired (for numerics)
        if (z_scale_in .ne. 1.0) then 
            if (.not. trc%par%is_sigma) z1 = z1*z_scale_in 
            z_srf1 = z_srf1*z_scale_in 
            H1     = H1*z_scale_in
            uz1    = uz1*z_scale_in 
        end if 

        if (trim(trc%par%interp_method) .eq. "spline") then

            ! Allocate z-velocity field in sigma coordinates 
            if (allocated(usig1)) deallocate(usig1)
            allocate(usig1(nx,ny,nz))

            usig1 = 0.0
            do k = 1, nz 
                where (H1 .gt. 0.0) usig1(:,:,k) = uz1(:,:,k) / H1
            end do 

            call interp_bspline3D_weights(bspline3d_ux,x1,y1,z1,ux1)
            call interp_bspline3D_weights(bspline3d_uy,x1,y1,z1,uy1)
            call interp_bspline3D_weights(bspline3d_uz,x1,y1,z1,usig1)
            
            write(*,*) "spline weights calculated."
        end if 


        ! Interpolate to the get the right elevation and other deposition quantities
        do i = 1, trc%par%n 

            if (trc%now%active(i) .eq. 2) then 

                ! Temporarily store velocity of this time step (for accelaration calculation)
                ux0 = trc%now%ux(i)
                uy0 = trc%now%uy(i)
                uz0 = trc%now%uz(i) 

                ! Linear interpolation used for surface position 
                par_lin = interp_bilinear_weights(x1,y1,xout=trc%now%x(i),yout=trc%now%y(i))
                trc%now%H(i)     = interp_bilinear(par_lin,H1)
                trc%now%z_srf(i) = interp_bilinear(par_lin,z_srf1)
                trc%now%z(i)     = trc%now%z_srf(i) - trc%now%dpth(i)

                ! Calculate zc-axis for the current point
                ! (z_bedrock + ice thickness)
                ! Note: equivalent to (z_srf - depth) = trc%now%z_srf(i) - (1.0-z1)*trc%now%H(i)
                zc = (trc%now%z_srf(i)-trc%now%H(i)) + z1*trc%now%H(i)

                if (trim(trc%par%interp_method) .eq. "linear") then 
                    ! Trilinear interpolation 

                    ! Note: currently we redundantly obtain bilinear (horizontal) weights, because
                    ! they are needed to calculate zc. In the future, this could be improved. 

                    par_lin = interp_trilinear_weights(x1,y1,zc,xout=trc%now%x(i),yout=trc%now%y(i),zout=trc%now%z(i))

                    trc%now%ux(i)  = interp_trilinear(par_lin,ux1)
                    trc%now%uy(i)  = interp_trilinear(par_lin,uy1)
                    trc%now%uz(i)  = interp_trilinear(par_lin,uz1)

                else
                    ! Spline interpolation 
                    trc%now%ux(i) = interp_bspline3D(bspline3d_ux,trc%now%x(i),trc%now%y(i),trc%now%z(i)/trc%now%H(i))
                    trc%now%uy(i) = interp_bspline3D(bspline3d_uy,trc%now%x(i),trc%now%y(i),trc%now%z(i)/trc%now%H(i))
                    trc%now%uz(i) = interp_bspline3D(bspline3d_uz,trc%now%x(i),trc%now%y(i),trc%now%z(i)/trc%now%H(i)) *trc%now%H(i)  ! sigma => m

                end if 

                ! Update acceleration term 
                trc%now%ax(i) = (trc%now%ux(i) - ux0) / trc%par%dt
                trc%now%ay(i) = (trc%now%uy(i) - uy0) / trc%par%dt
                trc%now%az(i) = (trc%now%uz(i) - uz0) / trc%par%dt

                ! Filler values of the tracer state, in the future these should
                ! equal the surface temperature and the accumulation rate at the time of
                ! deposition and be calculated otherwise
                trc%now%T(i)   = 260.0 
                trc%now%thk(i) = 0.3 

            end if 

        end do 

        ! Update the tracer thickness, then destroy points that are too thin 
        ! == TO DO == 

        ! Update the tracer positions 
        call calc_position(trc%now%x,trc%now%y,trc%now%z,trc%now%ux,trc%now%uy,trc%now%uz, &
                           trc%now%ax,trc%now%ay,trc%now%az,trc%par%dt,trc%now%active)

!         call calc_position(trc%now%x,trc%now%y,trc%now%dpth,trc%now%ux,trc%now%uy,-trc%now%uz,trc%par%dt,trc%now%active)
        trc%now%dpth = max(trc%now%z_srf - trc%now%z, 0.0) 

        ! Destroy points that moved outside the valid region 
        call tracer_deactivate(trc%par,trc%now,x1,y1,maxval(H1))

        ! Activate new tracers if desired
        if (dep_now) call tracer_activate(trc%par,trc%now,x1,y1,H=H1,nmax=trc%par%n_max_dep)

        ! Finish activation for necessary points 
        do i = 1, trc%par%n 

            if (trc%now%active(i) .eq. 1) then 
                ! Point became active now, further initializations needed below

                ! Point is at the surface, so only bilinear interpolation is needed
                par_lin = interp_bilinear_weights(x1,y1,xout=trc%now%x(i),yout=trc%now%y(i))

                ! Apply interpolation weights to variables
                trc%now%dpth(i)  = 0.001*z_scale_in   ! Always deposit just below the surface (eg 1 mm) to avoid zero z-velocity
                trc%now%z_srf(i) = interp_bilinear(par_lin,z_srf1)
                trc%now%z(i)     = trc%now%z_srf(i)-trc%now%dpth(i)
                
                trc%now%H(i)     = interp_bilinear(par_lin,H1)
                trc%now%ux(i)    = interp_bilinear(par_lin,ux1(:,:,nz))
                trc%now%uy(i)    = interp_bilinear(par_lin,uy1(:,:,nz))
                trc%now%uz(i)    = interp_bilinear(par_lin,uz1(:,:,nz)) 
                trc%now%ax(i)    = 0.0 
                trc%now%ay(i)    = 0.0 
                trc%now%az(i)    = 0.0
                
                ! Initialize state variables
                trc%now%T(i)   = 260.0 
                trc%now%thk(i) = 0.3 

                ! Define deposition values 
                trc%dep%time(i) = trc%par%time_now 
                trc%dep%H(i)    = trc%now%H(i)
                trc%dep%x(i)    = trc%now%x(i)
                trc%dep%y(i)    = trc%now%y(i)
                trc%dep%z(i)    = trc%now%z(i) 

                trc%now%active(i) = 2 

            end if 

        end do 


        ! == TO DO == 
        ! - Attach whatever information we want to trace (age, deposition elevation and location, climate, isotopes, etc)
        ! - Potentially attach this information via a separate subroutine, using a flag to see if
        !   it was just deposited, then in main program calling eg, tracer_add_dep_variable(trc,"T",T),
        !   where the argument "T" should match a list of available variables, and T should be the variable
        !   to be stored from the main program. 

        ! Update summary statistics 
        trc%par%n_active = count(trc%now%active.gt.0)

        ! Calculate some summary information on eulerian grid if desired 
        if (stats_now) call calc_tracer_stats(trc,x,y,z,z_srf,H)

        return 

    end subroutine tracer_update

    subroutine tracer_end(trc)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        
        ! Allocate the state variables 
        call tracer_deallocate(trc%now,trc%dep)

        write(*,*) "tracer:: tracer object deallocated."
        
        return 

    end subroutine tracer_end

    ! ================================================
    !
    ! tracer management 
    !
    ! ================================================
    
    subroutine tracer_activate(par,now,x,y,H,nmax)
        ! Use this to activate individual or multiple tracers (not more than nmax)
        ! Only determine x/y position here, later interpolate z_srf and deposition
        ! information 

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(prec), intent(IN) :: x(:), y(:)
        real(prec), intent(IN) :: H(:,:)
        integer, intent(IN) :: nmax  

        integer :: ntot  
        real(prec) :: p(size(H,1),size(H,2)), p_init(size(H,1),size(H,2))
        integer :: i, j, k, ij(2)
        real(prec), allocatable :: jit(:,:), dens(:,:)
        real(prec) :: xmin, ymin, xmax, ymax 

        ! How many points can be activated?
        ntot = min(nmax,count(now%active == 0))

        ! Determine initial desired distribution of points on low resolution grid
        p_init = gen_distribution(H,H_min=par%H_min_dep*z_scale_in,alpha=par%alpha,dist=par%weight)
        p = p_init  

        ! Generate random numbers to populate points 
        allocate(jit(2,ntot))

        if (par%noise) then 
            call random_number(jit)
            jit = (jit - 0.5)
            jit(1,:) = jit(1,:)*(x(2)-x(1)) 

            if (size(y,1) .gt. 2) then 
                jit(2,:) = jit(2,:)*(y(2)-y(1)) 
            else   ! Profile
                jit(2,:) = 0.0 
            end if 
        else 
            jit = 0.0 
        end if 

!         write(*,*) "range jit: ", minval(jit), maxval(jit)
!         write(*,*) "npts: ", count(now%active == 0)
!         write(*,*) "ntot: ", ntot 
!         stop 
    
        ! Calculate domain boundaries to be able to apply limits 
        xmin = minval(x) 
        xmax = maxval(x) 
        ymin = minval(y) 
        ymax = maxval(y) 

        if (maxval(p_init) .gt. 0.0) then 
            ! Activate points in locations with non-zero probability
            ! This if-statement ensures some valid points currently exist in the domain
            
            k = 0 

            do j = 1, par%n 

                if (now%active(j)==0) then 

                    now%active(j) = 1
                    k = k + 1
                    par%id_max = par%id_max+1 
                    now%id(j)  = par%id_max 

                    ij = maxloc(p,mask=p.gt.0.0)
                    now%x(j) = x(ij(1)) + jit(1,k)
                    now%y(j) = y(ij(2)) + jit(2,k)
                    
                    if (now%x(j) .lt. xmin) now%x(j) = xmin 
                    if (now%x(j) .gt. xmax) now%x(j) = xmax 
                    if (now%y(j) .lt. ymin) now%y(j) = ymin 
                    if (now%y(j) .gt. ymax) now%y(j) = ymax 
                    
                    p(ij(1),ij(2)) = 0.0 
                     
                end if 

                ! Stop when all points have been allocated
                if (k .ge. ntot) exit 

                ! If there are no more points with non-zero probability, reset probability
                if (maxval(p) .eq. 0.0) p = p_init 

            end do 
          
        end if

        ! Summary 
!         write(*,*) "tracer_activate:: ", count(now%active == 0), count(now%active .eq. 1), count(now%active .eq. 2)

        return 

    end subroutine tracer_activate 

    subroutine tracer_deactivate(par,now,x,y,Hmax)
        ! Use this to deactivate individual or multiple tracers
        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(prec), intent(IN) :: x(:), y(:) 
        real(prec), intent(IN) :: Hmax 

        ! Deactivate points where:
        !  - Thickness of ice sheet at point's location is below threshold
        !  - Point is above maximum ice thickness Hmax (interp error)
        !  - Point is past maximum depth into the ice sheet 
        !  - Velocity of point is higher than maximum threshold 
        !  - x/y position is out of boundaries of the domain 
        where (now%active .gt. 0 .and. &
              ( now%H .lt. par%H_min*z_scale_in                .or. &
                now%H .gt. Hmax                                .or. &
                now%dpth/now%H .ge. par%depth_max              .or. &
                sqrt(now%ux**2 + now%uy**2) .gt. par%U_max     .or. &
                now%x .lt. minval(x) .or. now%x .gt. maxval(x) .or. &
                now%y .lt. minval(y) .or. now%y .gt. maxval(y) ) ) 

            now%active    = 0 

            now%id        = mv 
            now%x         = mv 
            now%y         = mv 
            now%z         = mv 
            now%sigma     = mv 
            now%z_srf     = mv 
            now%dpth      = mv 
            now%ux        = mv 
            now%uy        = mv 
            now%uz        = mv 
            now%thk       = mv 
            now%T         = mv 
            now%H         = mv 

        end where 

        return 

    end subroutine tracer_deactivate 

    ! ================================================
    !
    ! tracer physics / stats
    !
    ! ================================================
    
    elemental subroutine calc_position(x,y,z,ux,uy,uz,ax,ay,az,dt,active)

        implicit none 

        real(prec),   intent(INOUT) :: x, y, z 
        real(prec),   intent(IN)    :: ux, uy, uz 
        real(prec),   intent(IN)    :: ax, ay, az 
        real(prec_time), intent(IN)    :: dt 
        integer,         intent(IN)    :: active 

        if (active .gt. 0) then 

            x = x + ux*dt + 0.5*ax*dt**2 
            y = y + uy*dt + 0.5*ay*dt**2 
            z = z + uz*dt + 0.5*az*dt**2
            
        end if 

        return 

    end subroutine calc_position

    function gen_distribution(H,H_min,alpha,dist) result(p)

        implicit none 

        real(prec), intent(IN) :: H(:,:)
        real(prec), intent(IN) :: H_min, alpha 
        character(len=*), intent(IN) :: dist 
        real(prec) :: p(size(H,1),size(H,2))

        ! Local variables
        integer    :: k, ij(2)
        real(prec) :: p_sum 

        select case(trim(dist)) 

            case("linear")

                p = (alpha * max(H-H_min,0.0) / (maxval(H)-H_min))

            case("quadratic")

                p = (alpha * max(H-H_min,0.0) / (maxval(H)-H_min))**2

            case DEFAULT   ! "rand"

                ! Random even distribution (all points equally likely)
                call random_number(p)
                where (H .lt. H_min) p = 0.0 


        end select 

        ! Normalize probability sum to one 
        p_sum = sum(p)
        if (p_sum .gt. 0.0) p = p / p_sum

        return 

    end function gen_distribution

    subroutine calc_tracer_stats(trc,x,y,z,z_srf,H)
        ! Convert tracer information to isochrone format matching
        ! Macgregor et al. (2015)
        ! Note: this should only be called at time t=0 ka BP, 
        ! since age is defined assuming that. 

        implicit none
        
        type(tracer_class), intent(INOUT) :: trc
        real(prec), intent(IN) :: x(:), y(:), z(:), z_srf(:,:), H(:,:)

        ! Local variables 
        integer :: i, j, k, q
        integer :: nx, ny, nz, nq
        real(prec), allocatable :: dx(:), dy(:) 
        real(prec) :: dt, dz, dz_now  
        real(prec) :: zc(size(z))
        integer       :: id(trc%par%n)
        real(prec)    :: dist(trc%par%n)
        integer, allocatable :: inds(:)
        integer :: n_ind 

        nx = size(x)
        ny = size(y)
        
        allocate(dx(nx+1),dy(ny+1))
        dx(2:nx) = (x(2:nx)-x(1:nx-1))/2.0
        dx(1)    = dx(2)
        dx(nx+1) = dx(nx)
        dy(2:ny) = (y(2:ny)-y(1:ny-1))/2.0
        dy(1)    = dy(2)
        dy(ny+1) = dy(ny)
        
        dt = 1.0 ! Isochrone uncertainty of 1 ka  
        dz = (trc%stats%depth_norm(2) - trc%stats%depth_norm(1))/2.0   ! depth_norm is equally spaced

        ! Loop over grid and fill in information
        do j = 1, ny 
        do i = 1, nx 

            ! Calculate the isochrones
            nq = size(trc%stats%age_iso)
            do q = 1, nq 

                ! Filter for active particles within the grid box and age of interest
                call which (trc%now%active == 2 .and. &
                            trc%now%x .gt. x(i)-dx(i) .and. trc%now%x .le. x(i)+dx(i+1) .and. &
                            trc%now%y .gt. y(j)-dy(j) .and. trc%now%y .le. y(j)+dy(j+1) .and. &
                            (0.0 - trc%dep%time)*1e-3 .ge. trc%stats%age_iso(q)-dt .and. &
                            (0.0 - trc%dep%time)*1e-3 .le. trc%stats%age_iso(q)+dt, inds, n_ind) 

                ! Calculate range mean/sd depth for given age range
                if (n_ind .gt. 0) then
                    write(*,*) "isochrones: ", i, j, q, n_ind 
 
                    trc%stats%depth_iso(i,j,q)     = calc_mean(real(trc%now%dpth(inds)*z_scale_out,prec_wrt))
                    trc%stats%depth_iso_err(i,j,q) = calc_sd(real(trc%now%dpth(inds)*z_scale_out,prec_wrt), &
                                                             trc%stats%depth_iso(i,j,q))
                    trc%stats%density_iso(i,j,q)   = n_ind 

                    trc%stats%dep_z_iso(i,j,q)     = calc_mean(real(trc%dep%z(inds)*z_scale_out,prec_wrt))
                else 
                    trc%stats%depth_iso(i,j,q)     = MV 
                    trc%stats%depth_iso_err(i,j,q) = MV
                    trc%stats%density_iso(i,j,q)   = MV
                    
                    trc%stats%dep_z_iso(i,j,q)     = MV 

                end if 

            end do 
            
            ! Calculate the ages of each depth layer 
            nq = size(trc%stats%depth_norm)
            do q = 1, nq 

                ! Use dz_now to ensure that the first depth (0.04) includes all depths to the surface
                dz_now = dz 
                if (q .eq. 1) dz_now = dz*5.0 

                ! Filter for active particles within the grid box and age of interest
                call which (trc%now%active == 2 .and. &
                            trc%now%x .gt. x(i)-dx(i) .and. trc%now%x .le. x(i)+dx(i+1) .and. &
                            trc%now%y .gt. y(j)-dy(j) .and. trc%now%y .le. y(j)+dy(j+1) .and. &
                            trc%now%dpth/trc%now%H .gt. trc%stats%depth_norm(q)-dz_now  .and. &
                            trc%now%dpth/trc%now%H .le. trc%stats%depth_norm(q)+dz, inds, n_ind) 

                ! Calculate range mean/sd age for given depth range
                if (n_ind .gt. 0) then
                    write(*,*) "ice_ages: ", i, j, q, n_ind 
 
                    trc%stats%ice_age(i,j,q)     = calc_mean(real(0.0-trc%dep%time(inds),prec_wrt))*1e-3
                    trc%stats%ice_age_err(i,j,q) = calc_sd(real(0.0-trc%dep%time(inds),prec_wrt),trc%stats%ice_age(i,j,q))*1e-3
                    trc%stats%density(i,j,q)     = n_ind 
                else 
                    trc%stats%ice_age(i,j,q)     = MV 
                    trc%stats%ice_age_err(i,j,q) = MV 
                    trc%stats%density(i,j,q)     = MV 
                end if  

            end do 

        end do 
        end do  
        
        return

    end subroutine calc_tracer_stats

    function calc_mean(x) result(mean)

        implicit none 

        real(prec_wrt), intent(IN) :: x(:) 
        real(prec_wrt) :: mean 
        integer :: n 

        n = count(x.ne.MV)

        if (n .gt. 0) then 
            mean = sum(x,mask=x.ne.MV) / real(n)
        else 
            mean = MV 
        end if 
        
        return 

    end function calc_mean 

    function calc_sd(x,mean) result(stdev)

        implicit none 

        real(prec_wrt), intent(IN) :: x(:) 
        real(prec_wrt) :: mean 
        real(prec_wrt) :: stdev 
        integer :: n 

        n = count(x.ne.MV)

        if (n .gt. 0) then 
            stdev = sqrt( sum((x - mean)**2) / real(n) )
        else 
            stdev = MV 
        end if 

        return 

    end function calc_sd 

    ! ================================================
    !
    ! Initialization routines 
    !
    ! ================================================

    subroutine tracer_par_load(par,filename,is_sigma)

        implicit none 

        type(tracer_par_class), intent(OUT) :: par 
        character(len=*),       intent(IN)  :: filename 
        logical, intent(IN) :: is_sigma 

!         par%n         = 5000
!         par%n_max_dep = 100
        
!         par%thk_min   = 1e-2   ! m (minimum thickness of tracer 'layer')
!         par%H_min     = 1500.0 ! m 
!         par%depth_max = 0.99   ! fraction of thickness
!         par%U_max     = 200.0  ! m/a 

!         par%H_min_dep = 1000.0 ! m 
!         par%alpha     = 1.0 
!         par%weight    = "linear"
        
!         par%dens_z_lim = 50.0 ! m
!         par%dens_max   = 10   ! Number of points
        
        call nml_read(filename,"tracer_par","n",            par%n)
        call nml_read(filename,"tracer_par","n_max_dep",    par%n_max_dep)
        call nml_read(filename,"tracer_par","thk_min",      par%thk_min)
        call nml_read(filename,"tracer_par","H_min",        par%H_min)
        call nml_read(filename,"tracer_par","depth_max",    par%depth_max)
        call nml_read(filename,"tracer_par","U_max",        par%U_max)
        call nml_read(filename,"tracer_par","H_min_dep",    par%H_min_dep)
        call nml_read(filename,"tracer_par","alpha",        par%alpha)
        call nml_read(filename,"tracer_par","weight",       par%weight)
        call nml_read(filename,"tracer_par","noise",        par%noise)
        call nml_read(filename,"tracer_par","dens_z_lim",   par%dens_z_lim)
        call nml_read(filename,"tracer_par","dens_max",     par%dens_max)
        call nml_read(filename,"tracer_par","interp_method",par%interp_method)
        
        ! Define additional parameter values
        par%is_sigma  = is_sigma 
        par%n_active  = 0 
        par%time_now  = 0.0             ! year
        par%time_old  = 0.0             ! year
        par%dt        = 0.0             ! year

        ! Consistency checks 
        if (trim(par%interp_method) .ne. "linear" .and. &
            trim(par%interp_method) .ne. "spline" ) then 
            write(0,*) "tracer_init:: error: interp_method must be 'linear' &
            &or 'spline': "//trim(par%interp_method)
            stop 
        end if 

        return 

    end subroutine tracer_par_load
    

    subroutine tracer_allocate(now,dep,n)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        integer, intent(IN) :: n

        ! Make object is deallocated
        call tracer_deallocate(now,dep)

        ! Allocate tracer 
        allocate(now%active(n))
        allocate(now%id(n))
        allocate(now%x(n),now%y(n),now%z(n),now%sigma(n))
        allocate(now%z_srf(n),now%dpth(n))
        allocate(now%ux(n),now%uy(n),now%uz(n))
        allocate(now%ax(n),now%ay(n),now%az(n))
        allocate(now%thk(n))
        allocate(now%T(n))
        allocate(now%H(n))

        ! Allocate deposition properties 
        
        allocate(dep%time(n), dep%H(n))
        allocate(dep%x(n), dep%y(n), dep%z(n), dep%lon(n), dep%lat(n))
        allocate(dep%t2m(4,n), dep%pr(4,n), dep%t2m_prann(n))
        
        return

    end subroutine tracer_allocate

    subroutine tracer_deallocate(now,dep)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        
        ! Deallocate state objects
        if (allocated(now%active))    deallocate(now%active)
        if (allocated(now%x))         deallocate(now%x)
        if (allocated(now%y))         deallocate(now%y)
        if (allocated(now%z))         deallocate(now%z)
        if (allocated(now%sigma))     deallocate(now%sigma)
        if (allocated(now%z_srf))     deallocate(now%z_srf)
        if (allocated(now%dpth))      deallocate(now%dpth)
        if (allocated(now%ux))        deallocate(now%ux)
        if (allocated(now%uy))        deallocate(now%uy)
        if (allocated(now%uz))        deallocate(now%uz)
        if (allocated(now%ax))        deallocate(now%ax)
        if (allocated(now%ay))        deallocate(now%ay)
        if (allocated(now%az))        deallocate(now%az)
        if (allocated(now%thk))       deallocate(now%thk)
        if (allocated(now%T))         deallocate(now%T)
        if (allocated(now%H))         deallocate(now%H)

        ! Deallocate deposition objects
        if (allocated(dep%time))      deallocate(dep%time)
        if (allocated(dep%z))         deallocate(dep%z)
        if (allocated(dep%H))         deallocate(dep%H)
        if (allocated(dep%x))         deallocate(dep%x)
        if (allocated(dep%y))         deallocate(dep%y)
        if (allocated(dep%lon))       deallocate(dep%lon)
        if (allocated(dep%lat))       deallocate(dep%lat)
        if (allocated(dep%t2m))       deallocate(dep%t2m)
        if (allocated(dep%pr))        deallocate(dep%pr)
        if (allocated(dep%t2m_prann)) deallocate(dep%t2m_prann)
        
        return

    end subroutine tracer_deallocate

    subroutine tracer_allocate_stats(stats,x,y)

        implicit none 

        type(tracer_stats_class), intent(INOUT) :: stats 
        real(prec), intent(IN) :: x(:), y(:)
        
        ! Make surce object is deallocated
        call tracer_deallocate_stats(stats)

        ! Allocate tracer stats axes
        allocate(stats%x(size(x)))
        allocate(stats%y(size(y)))

        ! Allocate tracer stats objects
        allocate(stats%depth_norm(25))  ! To match Macgregor et al. (2015)
        allocate(stats%age_iso(5))      ! To match Macgregor et al. (2015)
        allocate(stats%depth_iso(size(x),size(y),size(stats%age_iso)))
        allocate(stats%depth_iso_err(size(x),size(y),size(stats%age_iso)))
        allocate(stats%dep_z_iso(size(x),size(y),size(stats%age_iso)))
        allocate(stats%density_iso(size(x),size(y),size(stats%age_iso)))
        allocate(stats%ice_age(size(x),size(y),size(stats%depth_norm)))
        allocate(stats%ice_age_err(size(x),size(y),size(stats%depth_norm)))
        allocate(stats%density(size(x),size(y),size(stats%depth_norm)))
        
        ! Also store axis information directly
        stats%x = x 
        stats%y = y 

        ! Initialize arrays to zeros 
        stats%depth_iso     = 0.0 
        stats%depth_iso_err = 0.0
        stats%dep_z_iso     = MV  
        stats%density_iso   = MV 
        stats%ice_age       = 0.0 
        stats%ice_age_err   = 0.0 
        stats%density       = MV 
        
        return

    end subroutine tracer_allocate_stats

    subroutine tracer_deallocate_stats(stats)

        implicit none 

        type(tracer_stats_class), intent(INOUT) :: stats 

        ! Deallocate stats objects
        if (allocated(stats%x))             deallocate(stats%x)
        if (allocated(stats%y))             deallocate(stats%y)
        if (allocated(stats%depth_norm))    deallocate(stats%depth_norm)
        if (allocated(stats%age_iso))       deallocate(stats%age_iso)
        if (allocated(stats%depth_iso))     deallocate(stats%depth_iso)
        if (allocated(stats%depth_iso_err)) deallocate(stats%depth_iso_err)
        if (allocated(stats%dep_z_iso))     deallocate(stats%dep_z_iso)
        if (allocated(stats%density_iso))   deallocate(stats%density_iso)
        if (allocated(stats%ice_age))       deallocate(stats%ice_age)
        if (allocated(stats%ice_age_err))   deallocate(stats%ice_age_err)
        if (allocated(stats%density))       deallocate(stats%density)
        
        return

    end subroutine tracer_deallocate_stats

    subroutine tracer_reshape1D_vec(var,var1,rev)

        implicit none 
     
        real(prec),    intent(IN) :: var(:)
        real(prec), intent(INOUT), allocatable :: var1(:)
        logical,    intent(IN) :: rev 

        integer :: i, nx

        nx = size(var,1)
        if (allocated(var1)) deallocate(var1)
        allocate(var1(nx))

        if (rev) then 
            do i = 1, nx
                var1(i) = var(nx-i+1)
            end do 
        else 
            var1 = var 
        end if 

        return 

    end subroutine tracer_reshape1D_vec

    subroutine tracer_reshape2D_field(idx_order,var,var1)

        implicit none 

        character(len=3), intent(IN) :: idx_order 
        real(prec),    intent(IN) :: var(:,:)
        real(prec), intent(INOUT), allocatable :: var1(:,:)
        integer :: i, j
        integer :: nx, ny

        select case(trim(idx_order))

            case("ijk")
                ! x, y, z array order 

                nx = size(var,1)
                ny = size(var,2)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny))

                var1 = var 

            case("kji")
                ! z, y, x array order 

                nx = size(var,2)
                ny = size(var,1)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny))

                do i = 1, nx 
                do j = 1, ny 
                    var1(i,j)  = var(j,i)
                end do 
                end do 

            case DEFAULT 

                write(0,*) "tracer_reshape2D_field:: error: unrecognized array order: ",trim(idx_order)
                write(0,*) "    Possible choices are: ijk, kji"
                stop  

        end select 

        return 

    end subroutine tracer_reshape2D_field

    subroutine tracer_reshape3D_field(idx_order,var,var1,rev_z)

        implicit none 

        character(len=3), intent(IN) :: idx_order 
        real(prec),    intent(IN) :: var(:,:,:)
        real(prec), intent(INOUT), allocatable :: var1(:,:,:)
        logical,    intent(IN) :: rev_z   ! Reverse the z-axis? 
        integer :: i, j, k
        integer :: nx, ny, nz 

        select case(trim(idx_order))

            case("ijk")
                ! x, y, z array order 

                nx = size(var,1)
                ny = size(var,2)
                nz = size(var,3)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny,nz))

                if (rev_z) then 
                    do i = 1, nx 
                    do j = 1, ny  
                    do k = 1, nz 
                        var1(i,j,k)  = var(i,j,nz-k+1) 
                    end do 
                    end do 
                    end do 
                else 
                    var1 = var 
                end if 

            case("kji")
                ! z, y, x array order 

                nx = size(var,3)
                ny = size(var,2)
                nz = size(var,1)

                if (allocated(var1)) deallocate(var1)
                allocate(var1(nx,ny,nz))

                if (rev_z) then 
                    do i = 1, nx 
                    do j = 1, ny 
                    do k = 1, nz 
                        var1(i,j,k)  = var(nz-k+1,j,i)
                    end do 
                    end do 
                    end do 
                else 
                    do i = 1, nx 
                    do j = 1, ny 
                        var1(i,j,k)  = var(k,j,i) 
                    end do 
                    end do 
                end if 

            case DEFAULT 

                write(0,*) "tracer_reshape3D_field:: error: unrecognized array order: ",trim(idx_order)
                write(0,*) "    Possible choices are: ijk, kji"
                stop  

        end select 

        return 

    end subroutine tracer_reshape3D_field

    ! ================================================
    !
    ! I/O routines 
    !
    ! ================================================

    subroutine tracer_write_init(trc,fldr,filename)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        character(len=*), intent(IN)   :: fldr, filename 

        ! Local variables 
        integer :: nt 
        character(len=512) :: path_out 

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"pt",x=1,dx=1,nx=trc%par%n)
        call nc_write_dim(path_out,"time",x=real(mv,prec_wrt),unlimited=.TRUE.)

        return 

    end subroutine tracer_write_init 

    subroutine tracer_write(trc,time,fldr,filename)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec_time) :: time 
        character(len=*), intent(IN) :: fldr, filename 

        ! Local variables 
        integer :: nt
        integer, allocatable :: dims(:)
        real(prec_wrt) :: time_in, mv_wrt   
        real(prec_wrt) :: tmp(size(trc%now%x))
        character(len=512) :: path_out 

        trc%par%time_write = time 

        path_out = trim(fldr)//"/"//trim(filename)

        mv_wrt = MV 

        ! Determine which timestep this is
        call nc_dims(path_out,"time",dims=dims)
        nt = dims(1)
        call nc_read(path_out,"time",time_in,start=[nt],count=[1])
        if (time_in .ne. MV .and. abs(time-time_in).gt.1e-2) nt = nt+1 

        call nc_write(path_out,"time",real(time,prec_wrt), dim1="time",start=[nt],count=[1],missing_value=mv_wrt)
        call nc_write(path_out,"n_active",trc%par%n_active,dim1="time",start=[nt],count=[1],missing_value=int(mv_wrt))
        
        tmp = trc%now%x
        where(trc%now%x .ne. mv_wrt) tmp = trc%now%x*1e-3
        call nc_write(path_out,"x",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="km")
        tmp = trc%now%y
        where(trc%now%y .ne. mv_wrt) tmp = trc%now%y*1e-3
        call nc_write(path_out,"y",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="km")
        call nc_write(path_out,"z",real(trc%now%z*z_scale_out,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"dpth",real(trc%now%dpth*z_scale_out,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"z_srf",real(trc%now%z_srf*z_scale_out,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"ux",real(trc%now%ux,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m/a")
        call nc_write(path_out,"uy",real(trc%now%uy,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m/a")
        call nc_write(path_out,"uz",real(trc%now%uz*z_scale_out,kind=prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m/a")
        call nc_write(path_out,"thk",real(trc%now%thk,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        call nc_write(path_out,"T",real(trc%now%T,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"H",real(trc%now%H*z_scale_out,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")

        call nc_write(path_out,"id",trc%now%id,dim1="pt",dim2="time", missing_value=int(mv_wrt), &
                        start=[1,nt],count=[trc%par%n ,1])

        tmp = mv_wrt
        where(trc%dep%time .ne. mv_wrt) tmp = time-trc%dep%time
        call nc_write(path_out,"age",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="a")

        ! Write deposition information
        call nc_write(path_out,"dep_time",real(trc%dep%time,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="years")
        call nc_write(path_out,"dep_H",real(trc%dep%H*z_scale_out,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")
        tmp = trc%dep%x
        where(trc%dep%x .ne. mv_wrt) tmp = trc%dep%x*1e-3
        call nc_write(path_out,"dep_x",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="km")
        tmp = trc%dep%y
        where(trc%dep%y .ne. mv_wrt) tmp = trc%dep%y*1e-3
        call nc_write(path_out,"dep_y",tmp,dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="km")
        call nc_write(path_out,"dep_z",real(trc%dep%z*z_scale_out,prec_wrt),dim1="pt",dim2="time", missing_value=mv_wrt, &
                        start=[1,nt],count=[trc%par%n ,1],units="m")

        return 

    end subroutine tracer_write 

    subroutine tracer_write_stats(trc,time,fldr,filename) !,z_srf,H)
        ! Write various meta-tracer information (ie, lagrangian => eulerian)
        ! This output belongs to a specific time slice, usually at time = 0 ka BP. 

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec_time) :: time
        character(len=*),   intent(IN) :: fldr, filename 
!         real(prec),         intent(IN) :: z_srf(:,:), H(:,:) 

        ! Local variables 
        character(len=512) :: path_out 
        real(prec_wrt) :: mv_wrt 

        path_out = trim(fldr)//"/"//trim(filename)

        mv_wrt = MV 

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",        x=trc%stats%x*1e-3,     units="km")
        call nc_write_dim(path_out,"yc",        x=trc%stats%y*1e-3,     units="km")
        call nc_write_dim(path_out,"depth_norm",x=trc%stats%depth_norm, units="1")
        call nc_write_dim(path_out,"age_iso",   x=trc%stats%age_iso,    units="ka")
        call nc_write_dim(path_out,"time",      x=time,unlimited=.TRUE.,units="ka")
        
!         call nc_write(path_out,"z_srf",z_srf,dim1="xc",dim2="yc",missing_value=mv_wrt, &
!                       units="m",long_name="Surface elevation")
!         call nc_write(path_out,"H",H,dim1="xc",dim2="yc",missing_value=mv_wrt, &
!                       units="m",long_name="Ice thickness")

        call nc_write(path_out,"ice_age",trc%stats%ice_age,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age")
        call nc_write(path_out,"ice_age_err",trc%stats%ice_age_err,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=mv_wrt, &
                      units="ka",long_name="Layer age - error")
        call nc_write(path_out,"density",trc%stats%density,dim1="xc",dim2="yc",dim3="depth_norm",missing_value=int(mv_wrt), &
                      units="1",long_name="Tracer density")

        call nc_write(path_out,"depth_iso",trc%stats%depth_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone depth")
        call nc_write(path_out,"depth_iso_err",trc%stats%depth_iso_err,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone depth - error")
        call nc_write(path_out,"dep_z_iso",trc%stats%dep_z_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=mv_wrt, &
                      units="m",long_name="Isochrone deposition elevation")
        call nc_write(path_out,"density_iso",trc%stats%density_iso,dim1="xc",dim2="yc",dim3="age_iso",missing_value=int(mv_wrt), &
                      units="1",long_name="Tracer density (for isochrones)")

        
        return 

    end subroutine tracer_write_stats

    subroutine tracer_read()

        implicit none 



        return 

    end subroutine tracer_read

    subroutine tracer_align(trc_new,trc_ref,trc,dist_max)
        ! Interpolate tracer ages from trc to those of trc_ref 
        ! Note: for this to work well, trc should be sufficiently high
        ! resolution to minimize interpolation errors 

        implicit none 

        type(tracer_class), intent(OUT) :: trc_new 
        type(tracer_class), intent(IN) :: trc_ref, trc  
        real(prec), intent(IN) :: dist_max 

        ! Local variables 
        integer :: i, k  
        real(prec) :: dist(trc%par%n)

        ! Store reference tracer information in new object 
        trc_new = trc_ref 

        ! Make sure to set tagged info to missing, since
        ! this will not be valid for trc_new 
        trc_new%dep%time = MV 
        trc_new%dep%H    = MV 
        trc_new%dep%x    = MV 
        trc_new%dep%y    = MV 
        trc_new%dep%z    = MV 
        
        do i = 1, trc_new%par%n 

            if (trc_new%now%active(i) .eq. 2) then 
                ! Only treat active locations 

                dist = MV 
                where (trc%now%active .eq. 2)
                    dist = sqrt( (trc_new%now%x(i)-trc%now%x)**2 &
                               + (trc_new%now%y(i)-trc%now%y)**2)
                end where 

                k = minloc(dist,mask=dist.ne.MV,dim=1)

                if (dist(k) .le. dist_max) then 
                    trc_new%dep%time(i) = trc%dep%time(i)
                else 
                    trc_new%now%active(i) = 0
                    trc_new%now%H(i)      = MV 
                    trc_new%now%z_srf(i)  = MV
                    trc_new%now%x(i)      = MV 
                    trc_new%now%y(i)      = MV 
                    trc_new%now%z(i)      = MV 
         
                end if 

            end if 

        end do 

        return 

    end subroutine tracer_align

    subroutine which(x,ind,stat)
        ! Analagous to R::which function
        ! Returns indices that match condition x==.TRUE.

        implicit none 

        logical :: x(:)
        integer, allocatable :: tmp(:), ind(:)
        integer, optional :: stat  
        integer :: n, i  

        n = count(x)
        allocate(tmp(n))
        tmp = 0 

        n = 0
        do i = 1, size(x) 
            if (x(i)) then 
                n = n+1
                tmp(n) = i 
            end if
        end do 

        if (present(stat)) stat = n 

        if (allocated(ind)) deallocate(ind)

        if (n .eq. 0) then 
            allocate(ind(1))
            ind = -1 
        else
            allocate(ind(n))
            ind = tmp(1:n)
        end if 
        
        return 

    end subroutine which

end module tracer 


