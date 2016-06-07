
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
    end type 

    type(bilin_par_type) :: par_bilin 

    private 
    public :: tracer_class 
    public :: tracer_init 
    public :: tracer_update 
    public :: tracer_end 
    public :: tracer_write_init, tracer_write 

contains 

    subroutine tracer_init(trc,time)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        real(4) :: time 

        ! Load the parameters
        call tracer_par_load(trc%par)

        ! Allocate the state variables 
        call tracer_allocate(trc%now,trc%dep,n=trc%par%n)

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

    subroutine tracer_update(par,now,dep,time,x,y,z,z_srf,H,ux,uy,uz)

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
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
        real(prec), allocatable :: tmp(:,:)

        ! How many points can be activated?
        ntot = min(nmax,count(now%active == 0))

        ! Determine desired distribution of points on low resolution grid
        p = gen_distribution(H,H_min=par%H_min_dep,ntot=ntot,alpha=par%alpha,dist=par%dist)

        ! Generate random numbers to populate points 
        allocate(tmp(2,ntot))
        call random_number(tmp)
        tmp = (tmp - 0.5)
        tmp(1,:) = tmp(1,:)*(x(2)-x(1)) 
        tmp(2,:) = tmp(2,:)*(y(2)-y(1)) 

!         write(*,*) "range tmp: ", minval(tmp), maxval(tmp)
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
                    now%x(j) = x(ij(1)) + tmp(1,k)
                    now%y(j) = y(ij(2)) + tmp(2,k)
                    
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
        call nc_write(path_out,"x",trc%now%x,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"y",trc%now%y,dim1="pt",dim2="time", missing_value=MV, &
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

    subroutine tracer_read()

        implicit none 



        return 

    end subroutine tracer_read


end module tracer 


