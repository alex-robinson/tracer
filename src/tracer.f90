
module tracer 

    use tracer_precision
    use tracer_interp 
    use ncio   
    use nml 

    implicit none 

    
    type tracer_par_class 
        integer :: n, n_active, n_max_dep, id_max 
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
        
    end type 

    type tracer_state_class 
        integer,    allocatable :: active(:)
        integer,    allocatable :: id(:)
        real(prec), allocatable :: x(:), y(:), z(:)
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
        integer,    allocatable :: id(:)
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

    type(bilin_par_type) :: par_bilin 

    private 
    public :: tracer_class 
    public :: tracer_init 
    public :: tracer_update 
    public :: tracer_end 
    public :: tracer_write_init, tracer_write 
    public :: tracer_write_stats 

contains 

    subroutine tracer_init(trc,time,x,y,z)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        real(4) :: time 

        ! Load the parameters
        call tracer_par_load(trc%par)

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

        ! Initialize random number generator 
        call random_seed() 

        return 

    end subroutine tracer_init

    subroutine tracer_update(par,now,dep,stats,time,x,y,z,z_srf,H,ux,uy,uz)

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        type(tracer_stats_class), intent(INOUT) :: stats
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        real(prec), intent(IN) :: z_srf(:,:), H(:,:)
        real(prec), intent(IN) :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
        
        ! Local variables  
        integer :: i, nz 

        ! Update current time 
        par%time_old = par%time_now 
        par%time_now = time 
        par%dt       = par%time_now - par%time_old 

        nz = size(ux,3)


        ! Interpolate velocities to active point locations 

        ! Interpolate to the get the right elevation and other deposition quantities
        do i = 1, par%n 

            if (now%active(i) .eq. 2) then 
                par_bilin = interp_bilinear_weights(x,y,xout=now%x(i),yout=now%y(i))

                now%H(i)   = interp_bilinear(par_bilin,H)
                now%ux(i)  = interp_bilinear(par_bilin,ux(:,:,nz))
                now%uy(i)  = interp_bilinear(par_bilin,uy(:,:,nz))
!                 now%uz(i) = interp_bilinear(par_bilin,uz(:,:,nz))
                now%uz(i)  = 0.0 
                
                now%z(i)   = interp_bilinear(par_bilin,z_srf)
                now%T(i)   = 260.0 
                now%thk(i) = 0.3 

            end if 

        end do 


        ! == TO DO == 

        ! Update the tracer thickness, then destroy points that are too thin 

        ! Update the tracer positions 
        call calc_position(now%x,now%y,now%z,now%ux,now%uy,now%uz,par%dt,now%active)

        ! Destroy points that moved outside the valid region 
        call tracer_deactivate(par,now)

        
        ! Activate new tracers
        call tracer_activate(par,now,x,y,H=H,nmax=par%n_max_dep)

        ! Finish activation for necessary points 
        do i = 1, par%n 

            if (now%active(i) .eq. 1) then 
                ! Point became active now, further initializations needed below

                par_bilin = interp_bilinear_weights(x,y,xout=now%x(i),yout=now%y(i))

                now%z(i)   = interp_bilinear(par_bilin,z_srf)
                now%H(i)   = interp_bilinear(par_bilin,H)
                now%ux(i)  = interp_bilinear(par_bilin,ux(:,:,nz))
                now%uy(i)  = interp_bilinear(par_bilin,uy(:,:,nz))
!                 now%uz(i) = interp_bilinear(par_bilin,uz(:,:,nz))
                now%uz(i)  = 0.0 
                
                now%T(i)   = 260.0 
                now%thk(i) = 0.3 

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
        stats%density(:,:,1) = calc_tracer_density(par,now,x,y,z_srf)

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

        ! How many points can be activated?
        ntot = min(nmax,count(now%active == 0))

        ! Determine desired distribution of points on low resolution grid
        p = gen_distribution(H,H_min=par%H_min_dep,ntot=ntot,alpha=par%alpha,dist=par%dist)

        ! Generate random numbers to populate points 
        allocate(jit(2,ntot))
        call random_number(jit)
        jit = (jit - 0.5)
        jit(1,:) = jit(1,:)*(x(2)-x(1)) 
        jit(2,:) = jit(2,:)*(y(2)-y(1)) 

!         write(*,*) "range jit: ", minval(jit), maxval(jit)
!         write(*,*) "npts: ", count(npts .gt. 0) 
!         write(*,*) "ntot: ", ntot 
!         stop 

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
                    
                    p(ij(1),ij(2)) = 0.0 

                    exit 
                end if 

            end do 

            ! Stop when all points have been allocated
            if (k .ge. ntot) exit 

        end do 

        return 

    end subroutine tracer_activate 

    subroutine tracer_deactivate(par,now)
        ! Use this to deactivate individual or multiple tracers
        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 

        where (now%active .gt. 0 .and. &
                (now%H .lt. par%H_min .or. &
                 sqrt(now%ux**2 + now%uy**2) .gt. par%U_max) ) 

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
            if ( (now%x(k)-x(i-1))/(x(i)-x(i-1)) .lt. 0.5 ) i = i-1
            
            do j = 1, ny
                if (y(j) .gt. now%y(k)) exit
            end do
            if ( (now%y(k)-y(j-1))/(y(j)-y(j-1)) .lt. 0.5 ) j = j-1
            
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

    subroutine tracer_par_load(par)

        implicit none 

        type(tracer_par_class), intent(OUT) :: par 

        par%n         = 5000
        par%n_max_dep = 500
        par%n_active  = 0 

        par%time_now = 0.0 
        par%time_old = 0.0 
        par%dt       = 0.0 

        par%thk_min   = 1e-2   ! m 
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
        allocate(now%x(n),now%y(n),now%z(n))
        allocate(now%ux(n),now%uy(n),now%uz(n))
        allocate(now%thk(n))
        allocate(now%T(n))
        allocate(now%H(n))

        ! Allocate deposition properties 
        allocate(dep%id(n))
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

        path_out = trim(fldr)//"/"//trim(filename)

        ! Determine which timestep this is
        call nc_dims(path_out,"time",dims=dims)
        nt = dims(1)
        call nc_read(path_out,"time",time_in,start=[nt],count=[1])
        if (time_in .ne. MV .and. abs(time-time_in).gt.1e-2) nt = nt+1 

        call nc_write(path_out,"time",time,dim1="time",start=[nt],count=[1],missing_value=MV)
        call nc_write(path_out,"x",trc%now%x*1e-3,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"y",trc%now%y*1e-3,dim1="pt",dim2="time", missing_value=MV, &
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


        return 

    end subroutine tracer_write 

    subroutine tracer_write_stats(trc,time,fldr,filename)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec) :: time
        character(len=*), intent(IN)   :: fldr, filename 

        ! Local variables 
        integer :: nt 
        character(len=512) :: path_out 

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",   x=trc%stats%x*1e-3)
        call nc_write_dim(path_out,"yc",   x=trc%stats%y*1e-3)
        call nc_write_dim(path_out,"sigma",x=trc%stats%z)
        call nc_write_dim(path_out,"time",x=time,unlimited=.TRUE.)

!         call nc_write(path_out,"density",trc%stats%density,dim1="xc",dim2="yc",dim3="sigma",missing_value=int(MV), &
!                       units="1",long_name="Tracer density (surface)")
        call nc_write(path_out,"dens_srf",trc%stats%density(:,:,1),dim1="xc",dim2="yc",missing_value=int(MV), &
                      units="1",long_name="Tracer density (surface)")


        return 

    end subroutine tracer_write_stats

    subroutine tracer_read()

        implicit none 



        return 

    end subroutine tracer_read


end module tracer 


