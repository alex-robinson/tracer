
module tracer 

    use tracer_precision
    use tracer_interp 
    use bspline_module 
    use ncio   
    use nml 

    implicit none 

    
    type tracer_par_class 
        integer :: n, n_active, n_max_dep, id_max 
        logical :: is_sigma                     ! Is the defined z-axis in sigma coords
        real(prec) :: time_now, time_old, dt 
        real(prec) :: thk_min                   ! Minimum thickness of tracer (m)
        real(prec) :: H_min                     ! Minimum ice thickness to track (m)
        real(prec) :: depth_max                 ! Maximum depth of tracer (fraction)
        real(prec) :: U_max                     ! Maximum horizontal velocity of tracer to track (m/a)
        real(prec) :: H_min_dep                 ! Minimum ice thickness for tracer deposition (m)
        real(prec) :: alpha                     ! Slope of probability function
        character(len=56) :: dist               ! Distribution for generating prob. function
        real(prec) :: dens_z_lim                ! Distance from surface to count density
        integer    :: dens_max                  ! Max allowed density of particles at surface
        
        character(len=56) :: interp_method 
    end type 

    type tracer_state_class 
        integer,    allocatable :: active(:), id(:)
        real(prec), allocatable :: x(:), y(:), z(:), sigma(:)
        real(prec), allocatable :: ux(:), uy(:), uz(:)
        real(prec), allocatable :: thk(:)            ! Tracer thickness (for compression)
        real(prec), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
        real(prec), allocatable :: H(:)

    end type 

    type tracer_stats_class
        real(prec), allocatable :: x(:), y(:), z(:) 
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
    integer                    :: bspline_flag 

    private 

    ! For tracer2D module
    public :: tracer_par_class
    public :: tracer_state_class
    public :: tracer_dep_class
    public :: tracer_stats_class

    ! General public 
    public :: tracer_class 
    public :: tracer_init 
    public :: tracer_update 
    public :: tracer_end 
    public :: tracer_write_init, tracer_write 
    public :: tracer_write_stats 

contains 

    subroutine tracer_init(trc,time,x,y,z,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        logical,    intent(IN) :: is_sigma 

        real(4) :: time 

        ! Load the parameters
        call tracer_par_load(trc%par,is_sigma)

        ! Allocate the state variables 
        call tracer_allocate(trc%now,trc%dep,n=trc%par%n)
        call tracer_allocate_stats(trc%stats,x,y,z)

        ! Initialize state 
        trc%now%active    = 0 

        trc%now%id        = mv 
        trc%now%x         = mv 
        trc%now%y         = mv 
        trc%now%z         = mv 
        trc%now%ux        = mv 
        trc%now%uy        = mv 
        trc%now%uz        = mv 
        trc%now%thk       = mv 
        trc%now%T         = mv 
        trc%now%H         = mv 

        trc%par%id_max    = 0 

        ! Initialize the time 
        trc%par%time_now  = time 


        ! Method 
        trc%par%interp_method = "spline"    ! "linear" or "spline"

        ! Initialize random number generator 
        call random_seed() 


        ! Consistency checks 
        if (trim(trc%par%interp_method) .ne. "linear" .or. &
            trim(trc%par%interp_method) .ne. "spline" ) then 
            write(*,*) "tracer_init:: error: interp_method must be 'linear' &
            &or 'spline': "//trim(trc%par%interp_method)
        end if 

        return 

    end subroutine tracer_init

    subroutine tracer_update(par,now,dep,stats,time,x,y,z,z_srf,H,ux,uy,uz,order)

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        type(tracer_stats_class), intent(INOUT) :: stats
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        real(prec), intent(IN) :: z_srf(:,:), H(:,:)
        real(prec), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
        character(len=*), optional :: order 

        ! Local variables  
        character(len=3) :: idx_order 
        real(prec) :: zc(size(z))   ! Actual cartesian z-axis after applying sigma*H
        real(prec) :: z_srf_now  
        integer    :: i, j, k, nx, ny, nz, ksrf, kbase  
        real(prec), allocatable :: x1(:), y1(:), z1(:)
        real(prec), allocatable :: ux1(:,:,:), uy1(:,:,:), uz1(:,:,:)
        logical :: z_descending 

        real(prec), allocatable :: usig1(:,:,:)

        ! Update current time 
        par%time_old = par%time_now 
        par%time_now = time 
        par%dt       = par%time_now - par%time_old 


        ! Determine order of indices (default ijk)
        idx_order = "ijk"
        if (present(order)) idx_order = trim(order)

        ! Also determine whether z-axis is ascending or descending 
        call tracer_reshape3D(idx_order,x,y,z,ux,uy,uz,x1,y1,z1,ux1,uy1,uz1)

        ! Note: GRISLI (nx,ny,nz): sigma goes from 1 to 0, so sigma(1)=1 [surface], sigma(nz)=0 [base]
        !       SICO (nz,ny,nx): sigma(1) = 0, sigma(nz) = 1
        ! Using tracer_reshape3D homogenizes them to ascending z-axis (nx,ny,nz)
        
        nz = size(ux1,3)
        z_descending = .FALSE. 
        if (z1(1) .gt. z1(nz)) z_descending = .TRUE. 

        ! Determine size of y-axis, if it is one, this is a 2D profile domain
        ny = size(y1,1)

        ! Get nx for completeness 
        nx = size(x1,1)

        kbase = 1
        ksrf = nz 
        if (z_descending) then 
            kbase = nz 
            ksrf  = 1 
        end if

        if (trim(par%interp_method) .eq. "spline") then  
            call bspline3d_ux%initialize(dble(x1),dble(y1),dble(z1),dble(ux1),kx=4,ky=4,kz=4,iflag=bspline_flag)
            call bspline3d_uy%initialize(dble(x1),dble(y1),dble(z1),dble(uy1),kx=4,ky=4,kz=4,iflag=bspline_flag)
            call bspline3d_uz%initialize(dble(x1),dble(y1),dble(z1),dble(uz1),kx=4,ky=4,kz=4,iflag=bspline_flag)
            
            ! Allocate z-velocity field in sigma coordinates 
            if (allocated(usig1)) deallocate(usig1)
            allocate(usig1(nx,ny,nz))

            usig1 = 0.0
            do k = 1, nz 
                where (H .gt. 0.0) usig1(:,:,k) = uz1(:,:,k) / H
            end do 
        end if 

        stop 

        ! Interpolate to the get the right elevation and other deposition quantities
        do i = 1, par%n 

            if (now%active(i) .eq. 2) then 

                ! == TO DO: We need 3D interpolation here!!!

                ! 1. Calc bilin weights at x/y locations
                ! 2. Get H and z-axis at that location
                ! 3. Calculate lin weights for z values 
                ! 4. Perform trilinear interpolation 

                ! Linear interpolation used for surface position 
                par_lin    = interp_bilinear_weights(x1,y1,xout=now%x(i),yout=now%y(i))
                now%H(i)   = interp_bilinear(par_lin,H)
                z_srf_now  = interp_bilinear(par_lin,z_srf)

                ! Calculate zc-axis for the current point
                ! (z_bedrock + ice thickness)
                zc = (z_srf_now-now%H(i)) + z1*now%H(i)

                if (trim(par%interp_method) .eq. "linear") then 
                    ! Trilinear interpolation 

                    par_lin = interp_trilinear_weights(x1,y1,zc,xout=now%x(i),yout=now%y(i),zout=now%z(i))

                    now%ux(i)  = interp_trilinear(par_lin,ux1)
                    now%uy(i)  = interp_trilinear(par_lin,uy1)
                    now%uz(i)  = interp_trilinear(par_lin,uz1)

    !                 now%ux(i)  = interp_bilinear(par_lin,ux1(:,:,nz))
    !                 now%uy(i)  = interp_bilinear(par_lin,uy1(:,:,nz))
    !                 now%uz(i)  = interp_bilinear(par_lin,uz1(:,:,nz))
    !                 ! Until trilinear interp is ready, maintain z-position at surface
    !                 now%z(i)   = interp_bilinear(par_lin,z_srf)

                else
                    ! Spline interpolation 
                    
                end if 

                now%T(i)   = 260.0 
                now%thk(i) = 0.3 

            end if 

        end do 

        ! Update the tracer thickness, then destroy points that are too thin 
        ! == TO DO == 

        ! Update the tracer positions 
        call calc_position(now%x,now%y,now%z,now%ux,now%uy,now%uz,par%dt,now%active)

        ! Destroy points that moved outside the valid region 
        call tracer_deactivate(par,now,x1,y1)

        ! Activate new tracers
        call tracer_activate(par,now,x1,y1,H=H,nmax=par%n_max_dep)

        ! Finish activation for necessary points 
        do i = 1, par%n 

            if (now%active(i) .eq. 1) then 
                ! Point became active now, further initializations needed below

                ! Point is at the surface, so only bilinear interpolation is needed
                par_lin = interp_bilinear_weights(x1,y1,xout=now%x(i),yout=now%y(i))

                ! Apply interpolation weights to variables
                now%z(i)   = interp_bilinear(par_lin,z_srf)
                now%H(i)   = interp_bilinear(par_lin,H)
                now%ux(i)  = interp_bilinear(par_lin,ux1(:,:,ksrf))
                now%uy(i)  = interp_bilinear(par_lin,uy1(:,:,ksrf))
                now%uz(i)  = interp_bilinear(par_lin,uz1(:,:,ksrf)) 

                ! Initialize state variables
                now%T(i)   = 260.0 
                now%thk(i) = 0.3 

                ! Define deposition values 
                dep%time(i) = par%time_now 

                now%active(i) = 2 

            end if 

        end do 


        ! == TO DO == 
        ! - Generate position based on a random algorithm plus weighting
        ! - Weighting function based on H and U for location, input par deposition frequency 
        ! - Attach whatever information we want to trace (age, deposition elevation and location, climate, isotopes, etc)


        ! Update summary statistics 
        par%n_active = count(now%active.gt.0)

        ! Calculate density 
        stats%density = 0 
        stats%density(:,:,nz) = calc_tracer_density(par,now,x1,y1,z_srf)

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
        real(prec) :: p(size(H,1),size(H,2))
        integer :: i, j, k, ij(2)
        real(prec), allocatable :: jit(:,:), dens(:,:)
        real(prec) :: xmin, ymin, xmax, ymax 

        ! How many points can be activated?
        ntot = min(nmax,count(now%active == 0))

        ! Determine desired distribution of points on low resolution grid
        p = gen_distribution(H,H_min=par%H_min_dep,ntot=ntot,alpha=par%alpha,dist=par%dist)

        ! Generate random numbers to populate points 
        allocate(jit(2,ntot))
        call random_number(jit)
        jit = (jit - 0.5)
        jit(1,:) = jit(1,:)*(x(2)-x(1)) 

        if (size(y,1) .gt. 2) then 
            jit(2,:) = jit(2,:)*(y(2)-y(1)) 
        else   ! Profile
            jit(2,:) = 0.0 
        end if 

!         write(*,*) "range jit: ", minval(jit), maxval(jit)
!         write(*,*) "npts: ", count(npts .gt. 0) 
!         write(*,*) "ntot: ", ntot 
!         stop 
    
        ! Calculate domain boundaries to apply limits 
        xmin = minval(x) 
        xmax = maxval(x) 
        ymin = minval(y) 
        ymax = maxval(y) 

        do i = 1, count(p .gt. 0.0)
            ij = maxloc(p,mask=p.gt.0.0)

            k = 0 
            do j = 1, par%n 

                if (now%active(j)==0) then 

                    now%active(j) = 1
                    k = k + 1
                    par%id_max = par%id_max+1 
                    now%id(j)  = par%id_max 
                    now%x(j) = x(ij(1)) + jit(1,k)
                    now%y(j) = y(ij(2)) + jit(2,k)
                    
                    if (now%x(j) .lt. xmin) now%x(j) = xmin 
                    if (now%x(j) .gt. xmax) now%x(j) = xmax 
                    if (now%y(j) .lt. ymin) now%y(j) = ymin 
                    if (now%y(j) .gt. ymax) now%y(j) = ymax 
                    
                    p(ij(1),ij(2)) = 0.0 

                    exit 
                end if 

            end do 

            ! Stop when all points have been allocated
            if (k .ge. ntot) exit 

        end do 

        return 

    end subroutine tracer_activate 

    subroutine tracer_deactivate(par,now,x,y)
        ! Use this to deactivate individual or multiple tracers
        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(prec), intent(IN) :: x(:), y(:) 

        ! Deactivate points where:
        !  - Thickness of ice sheet at point's location is below threshold
        !  - Point is past maximum depth into the ice sheet 
        !  - Velocity of point is higher than maximum threshold 
        !  - x/y position is out of boundaries of the domain 
        where (now%active .gt. 0 .and. &
                (now%H .lt. par%H_min .or. &
                (now%H - now%z)/now%H .ge. par%depth_max       .or. &
                sqrt(now%ux**2 + now%uy**2) .gt. par%U_max     .or. &
                now%x .lt. minval(x) .or. now%x .gt. maxval(x) .or. &
                now%y .lt. minval(y) .or. now%y .gt. maxval(y)) ) 

            now%active    = 0 

            now%id        = mv 
            now%x         = mv 
            now%y         = mv 
            now%z         = mv 
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
    
    elemental subroutine calc_position(x,y,z,ux,uy,uz,dt,active)

        implicit none 

        real(prec), intent(INOUT) :: x, y, z 
        real(prec), intent(IN)    :: ux, uy, uz 
        real(prec), intent(IN)    :: dt 
        integer,    intent(IN)    :: active 

        if (active .gt. 0) then 
            x = x + ux*dt 
            y = y + uy*dt 
            z = z + uz*dt 
        end if 

        return 

    end subroutine calc_position

    function gen_distribution(H,H_min,ntot,alpha,dist) result(p)

        implicit none 

        real(prec), intent(IN) :: H(:,:)
        real(prec), intent(IN) :: H_min, alpha 
        integer, intent(IN)    :: ntot           ! Total points to allocate
        character(len=*), intent(IN) :: dist 
        real(prec) :: p(size(H,1),size(H,2))

        ! Local variables
        integer    :: k, ij(2)


        select case(trim(dist)) 

            case("linear")

                p = (alpha * max(H-H_min,0.0) / (maxval(H)-H_min))

            case("quadratic")

                p = (alpha * max(H-H_min,0.0) / (maxval(H)-H_min))**2

            case DEFAULT

                ! Even distribution (all points equally likely)
                p = 1.0
                where (H .lt. H_min) p = 0.0 


        end select 

        ! Now make sure the cumulative probabilities sum to number of points to allocate 
        p = p / sum(p)

        return 

    end function gen_distribution

    function calc_tracer_density(par,now,x,y,z_srf) result(dens)
        ! Check tracer density near the surface
        ! (to avoid depositing too many particles in the same place)
        
        implicit none
        
        type(tracer_par_class) :: par
        type(tracer_state_class), intent(IN) :: now
        real(prec), intent(IN) :: x(:), y(:), z_srf(:,:)
        integer :: dens(size(x),size(y))
        
        integer :: i, j, k
        integer :: nx, ny
        
        nx = size(x)
        ny = size(y)
        
        ! Set initial density to zero everywhere
        dens = 0
        
        ! Loop over active points, count surface density 
        ! at each grid point
        do k = 1, par%n_active
            
            do i = 1, nx
                if (x(i) .gt. now%x(k)) exit
            end do
            if (i .gt. nx) then 
                i = nx
            else if (i .gt. 1) then
                if ( (now%x(k)-x(i-1))/(x(i)-x(i-1)) .lt. 0.5 ) i = i-1
            end if

            do j = 1, ny
                if (y(j) .gt. now%y(k)) exit
            end do
            if (j .gt. ny) then
                j = ny 
            else if (i .gt. 1) then
                if ( (now%y(k)-y(j-1))/(y(j)-y(j-1)) .lt. 0.5 ) j = j-1
            end if

            if (abs(z_srf(i,j)-now%z(k)) .lt. par%dens_z_lim) then
                ! If point is near surface, add it to density
               dens(i,j) = dens(i,j)+1
            end if
            
            
        end do
          
        
                
                
        return

    end function calc_tracer_density



    ! ================================================
    !
    ! Initialization routines 
    !
    ! ================================================

    subroutine tracer_par_load(par,is_sigma)

        implicit none 

        type(tracer_par_class), intent(OUT) :: par 
        logical, intent(IN) :: is_sigma 

        par%n         = 10000
        par%n_max_dep = 500
        par%n_active  = 0 

        par%is_sigma  = is_sigma 

        par%time_now = 0.0 
        par%time_old = 0.0 
        par%dt       = 0.0 

        par%thk_min   = 1e-2   ! m (minimum thickness of tracer 'layer')
        par%H_min     = 1500.0 ! m 
        par%depth_max = 0.99   ! fraction of thickness
        par%U_max     = 200.0  ! m/a 

        par%H_min_dep = 2000.0 ! m 
        par%alpha     = 1.0 
        par%dist      = "linear"
        
        par%dens_z_lim = 50.0 ! m
        par%dens_max   = 10   ! Number of points

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
        allocate(now%x(n),now%y(n),now%z(n), now%sigma(n))
        allocate(now%ux(n),now%uy(n),now%uz(n))
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
        if (allocated(now%active)) deallocate(now%active)
        if (allocated(now%x))         deallocate(now%x)
        if (allocated(now%y))         deallocate(now%y)
        if (allocated(now%z))         deallocate(now%z)
        if (allocated(now%sigma))     deallocate(now%sigma)
        if (allocated(now%ux))        deallocate(now%ux)
        if (allocated(now%uy))        deallocate(now%uy)
        if (allocated(now%uz))        deallocate(now%uz)
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

    subroutine tracer_allocate_stats(stats,x,y,z)

        implicit none 

        type(tracer_stats_class), intent(INOUT) :: stats 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        
        ! Make surce object is deallocated
        call tracer_deallocate_stats(stats)

        ! Allocate tracer stats axes
        allocate(stats%x(size(x)))
        allocate(stats%y(size(y)))
        allocate(stats%z(size(z)))

        ! Allocate tracer stats objects
        allocate(stats%density(size(x),size(y),size(z)))

        ! Also store axis information directly
        stats%x = x 
        stats%y = y 
        stats%z = z 

        return

    end subroutine tracer_allocate_stats

    subroutine tracer_deallocate_stats(stats)

        implicit none 

        type(tracer_stats_class), intent(INOUT) :: stats 

        ! Deallocate stats objects
        if (allocated(stats%density)) deallocate(stats%density)

        return

    end subroutine tracer_deallocate_stats

    
    subroutine tracer_reshape3D(idx_order,x,y,z,ux,uy,uz,x1,y1,z1,ux1,uy1,uz1)

        implicit none 

        character(len=3), intent(IN) :: idx_order 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        real(prec), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)

        real(prec), intent(INOUT), allocatable :: x1(:), y1(:), z1(:)
        real(prec), intent(INOUT), allocatable :: ux1(:,:,:), uy1(:,:,:), uz1(:,:,:)
        integer :: i, j, k
        integer :: nx, ny, nz 

        nx = size(x)
        ny = size(y)
        nz = size(z) 

        if (allocated(x1))  deallocate(x1)
        if (allocated(y1))  deallocate(y1)
        if (allocated(z1))  deallocate(z1)
        if (allocated(ux1)) deallocate(ux1)
        if (allocated(uy1)) deallocate(uy1)
        if (allocated(uz1)) deallocate(uz1)

        allocate(x1(nx))
        allocate(y1(ny))
        allocate(z1(nz))
        allocate(ux1(nx,ny,nz))
        allocate(uy1(nx,ny,nz))
        allocate(uz1(nx,ny,nz))
        
        x1 = x 
        y1 = y 

        select case(trim(idx_order))

            case("ijk")
                ! x, y, z array order 

!                 if (z(1) .lt. z(nz)) then 
!                     ! Already ascending z-axis

                    z1  = z  
                    ux1 = ux 
                    uy1 = uy 
                    uz1 = uz 

!                 else   
!                     ! Reversed z-axis 
!                     do k = 1, nz 
!                         z1(k)      = z(nz-k+1)
!                         ux1(:,:,k) = ux(:,:,nz-k+1)
!                         uy1(:,:,k) = uy(:,:,nz-k+1)
!                         uz1(:,:,k) = uz(:,:,nz-k+1)
!                     end do 
!                 end if 

            case("kji")
                ! z, y, x array order 

!                 if (z(1) .lt. z(nz)) then 
!                     ! Already ascending z-axis

                    z1  = z  
                    do i = 1, nx 
                    do j = 1, ny 
                        ux1(i,j,:) = ux(:,j,i)
                        uy1(i,j,:) = uy(:,j,i)
                        uz1(i,j,:) = uz(:,j,i)
                    end do 
                    end do 

!                 else   
!                     ! Reversed z-axis 
!                     z1  = z(nz:1)
!                     do i = 1, nx 
!                     do j = 1, ny
!                         do k = 1, nz 
!                             ux1(i,j,k) = ux(nz-k+1,j,i)
!                             uy1(i,j,k) = uy(nz-k+1,j,i)
!                             uz1(i,j,k) = uz(nz-k+1,j,i)
!                         end do 
!                     end do 
!                     end do 

!                 end if 

            case DEFAULT 

                write(*,*) "tracer_reshape3D:: error: unrecognized array order: ",trim(idx_order)
                write(*,*) "    Possible choices are: ijk, kji"
                stop  

        end select 

        return 

    end subroutine tracer_reshape3D


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
        call nc_write_dim(path_out,"time",x=mv,unlimited=.TRUE.)

        return 

    end subroutine tracer_write_init 

    subroutine tracer_write(trc,time,fldr,filename)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec) :: time 
        character(len=*), intent(IN)   :: fldr, filename 

        ! Local variables 
        integer :: nt
        integer, allocatable :: dims(:)
        real(prec) :: time_in  
        character(len=512) :: path_out 

        real(prec) :: tmp(size(trc%now%x))

        path_out = trim(fldr)//"/"//trim(filename)

        ! Determine which timestep this is
        call nc_dims(path_out,"time",dims=dims)
        nt = dims(1)
        call nc_read(path_out,"time",time_in,start=[nt],count=[1])
        if (time_in .ne. MV .and. abs(time-time_in).gt.1e-2) nt = nt+1 

        call nc_write(path_out,"time",time,dim1="time",start=[nt],count=[1],missing_value=MV)
        
        tmp = trc%now%x
        where(trc%now%x .ne. MV) tmp = trc%now%x*1e-3
        call nc_write(path_out,"x",tmp,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        tmp = trc%now%y
        where(trc%now%y .ne. MV) tmp = trc%now%y*1e-3
        call nc_write(path_out,"y",tmp,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"z",trc%now%z,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"ux",trc%now%ux,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"uy",trc%now%uy,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"uz",trc%now%uz,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"thk",trc%now%thk,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"T",trc%now%T,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"H",trc%now%H,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])

        call nc_write(path_out,"id",trc%now%id,dim1="pt",dim2="time", missing_value=int(MV), &
                        start=[1,nt],count=[trc%par%n ,1])

        ! Write deposition information
        call nc_write(path_out,"dep_time",trc%dep%time,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])

        return 

    end subroutine tracer_write 

    subroutine tracer_write_stats(trc,time,fldr,filename)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec) :: time
        character(len=*), intent(IN)   :: fldr, filename 

        ! Local variables 
        integer :: nt, nz 
        character(len=512) :: path_out 

        nz = size(trc%stats%density,3)

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",   x=trc%stats%x*1e-3)
        call nc_write_dim(path_out,"yc",   x=trc%stats%y*1e-3)
        call nc_write_dim(path_out,"sigma",x=trc%stats%z)
        call nc_write_dim(path_out,"time",x=time,unlimited=.TRUE.)

!         call nc_write(path_out,"density",trc%stats%density,dim1="xc",dim2="yc",dim3="sigma",missing_value=int(MV), &
!                       units="1",long_name="Tracer density (surface)")
        call nc_write(path_out,"dens_srf",trc%stats%density(:,:,nz),dim1="xc",dim2="yc",missing_value=int(MV), &
                      units="1",long_name="Tracer density (surface)")


        return 

    end subroutine tracer_write_stats

    subroutine tracer_read()

        implicit none 



        return 

    end subroutine tracer_read


end module tracer 


