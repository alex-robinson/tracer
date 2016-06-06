
module tracer 

    use tracer_precision
    use tracer_interp 
    use ncio   
    use nml 

    implicit none 

    
    type tracer_par_class 
        integer :: n, n_active
        real(prec) :: time_now, time_old, dt 
        real(prec) :: thk_lim                   ! Minimum thickness of tracer (m)
        real(prec) :: z_lim                     ! Minimum depth of tracer (fraction)
    end type 

    type tracer_state_class 
        logical, allocatable :: is_active(:) 
        real(prec), allocatable :: x(:), y(:), z(:)
        real(prec), allocatable :: u(:), v(:), w(:)
        real(prec), allocatable :: thk(:)            ! Tracer thickness (for compression)
        real(prec), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
    end type 

    type tracer_dep_class 
        ! Standard deposition information (time and place)
        real(prec), allocatable :: time(:) 
        real(prec), allocatable :: elev(:) 
        real(prec), allocatable :: x(:), y(:)
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

    subroutine tracer_init(trc)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        
        ! Load the parameters
        call tracer_par_load(trc%par)

        ! Allocate the state variables 
        call tracer_allocate(trc%now,trc%dep,n=trc%par%n)

        ! Initialize state 
        trc%now%is_active = .FALSE. 

        trc%now%x         = mv 
        trc%now%y         = mv 
        trc%now%z         = mv 
        trc%now%u         = mv 
        trc%now%v         = mv 
        trc%now%w         = mv 
        trc%now%thk       = mv 
        trc%now%T         = mv 

        ! Initialize random number generator 
        call random_seed() 

        return 

    end subroutine tracer_init

    subroutine tracer_update(par,now,dep,time,x,y,z,z_srf,H) !,u,v,w)

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), y(:), z(:)
        real(prec), intent(IN) :: z_srf(:,:), H(:,:)
!         real(prec), intent(IN) :: u(:,:,:), v(:,:,:), w(:,:,:)
        
        ! Local variables 
        real(prec) :: xlim(2), ylim(2) 
        integer :: i 

        ! Update current time 
        par%time_old = par%time_now 
        par%time_now = time 
        par%dt       = par%time_now - par%time_old 

        ! Interpolate velocities to active point locations 

        ! == TO DO == 

        ! Update the tracer thickness, then destroy points that are too thin 

        ! == TO DO == 

        ! Update the tracer positions 
        call calc_position(now%x,now%y,now%z,now%u,now%v,now%w,par%dt,now%is_active)

        ! Destroy points that moved outside the valid region 

        ! == TO DO == 

        ! Finally, deposit new tracer points
        xlim = [minval(x),maxval(x)]
        ylim = [minval(y),maxval(y)]
        
        call tracer_activate(par,now,xlim,ylim,H=H,nmax=20,method="random")

        ! Interpolate to the get the right elevation and other deposition quantities
        do i = 1, par%n 
            if (now%is_active(i)) then 
                par_bilin = interp_bilinear_weights(x,y,xout=now%x(i),yout=now%y(i))

                now%z(i) = interp_bilinear(par_bilin,z_srf)

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
    
    subroutine tracer_activate(par,now,xlim,ylim,H,nmax,method)
        ! Use this to activate individual or multiple tracers (not more than nmax)
        ! Only determine x/y position here, later interpolate z_srf and deposition
        ! information 

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(prec), intent(IN) :: xlim(2), ylim(2) 
        real(prec), intent(IN) :: H(:,:)
        integer, intent(IN) :: nmax 
        character(len=*), intent(IN) :: method 

        integer :: i, k 
        real(prec) :: tmp(2) 

        k = 0 

        if (trim(method) .eq. "random") then 

            do i = 1, par%n 

                if (.not. now%is_active(i)) then 

                    call random_number(tmp)

                    now%x(i) = tmp(1)*(xlim(2)-xlim(1)) + xlim(1)
                    now%y(i) = tmp(2)*(ylim(2)-ylim(1)) + ylim(1)
                    
                    ! Mark this point as active
                    now%is_active(i) = .TRUE. 
                    k = k+1 

                end if 

                ! Only activate the number of desired points
                if (k .ge. nmax) exit 

            end do 

        end if 


        return 

    end subroutine tracer_activate 

    subroutine tracer_deactivate()
        ! Use this to deactivate individual or multiple tracers
        implicit none 


        return 

    end subroutine tracer_deactivate 

    ! ================================================
    !
    ! tracer physics 
    !
    ! ================================================
    
    elemental subroutine calc_position(x,y,z,u,v,w,dt,is_active)

        implicit none 

        real(prec), intent(INOUT) :: x, y, z 
        real(prec), intent(IN)    :: u, v, w 
        real(prec), intent(IN)    :: dt 
        logical,    intent(IN)    :: is_active 

        if (is_active) then 
            x = x + u*dt 
            y = y + v*dt 
            z = z + w*dt 
        end if 

        return 

    end subroutine calc_position

    ! ================================================
    !
    ! Initialization routines 
    !
    ! ================================================

    subroutine tracer_par_load(par)

        implicit none 

        type(tracer_par_class), intent(OUT) :: par 

        par%n        = 1000
        par%n_active = 0 

        par%time_now = 0.0 
        par%time_old = 0.0 
        par%dt       = 0.0 

        par%thk_lim  = 1e-2  ! m 
        par%z_lim    = 0.99  ! fraction of thickness


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
        allocate(now%is_active(n))
        allocate(now%x(n),now%y(n),now%z(n))
        allocate(now%u(n),now%v(n),now%w(n))
        allocate(now%thk(n))
        allocate(now%T(n))

        ! Allocate deposition properties 
        allocate(dep%time(n), dep%elev(n))
        allocate(dep%x(n), dep%y(n), dep%lon(n), dep%lat(n))
        allocate(dep%t2m(4,n), dep%pr(4,n), dep%t2m_prann(n))
        
        return

    end subroutine tracer_allocate

    subroutine tracer_deallocate(now,dep)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        
        ! Deallocate state objects
        if (allocated(now%is_active)) deallocate(now%is_active)
        if (allocated(now%x))         deallocate(now%x)
        if (allocated(now%y))         deallocate(now%y)
        if (allocated(now%z))         deallocate(now%z)
        if (allocated(now%u))         deallocate(now%u)
        if (allocated(now%v))         deallocate(now%v)
        if (allocated(now%w))         deallocate(now%w)
        if (allocated(now%thk))       deallocate(now%thk)
        if (allocated(now%T))         deallocate(now%T)

        ! Deallocate deposition objects
        if (allocated(dep%time))      deallocate(dep%time)
        if (allocated(dep%elev))      deallocate(dep%elev)
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

!         type tracer_state_class 
!         logical, allocatable :: is_active(:) 
!         real(prec), allocatable :: x(:), y(:), z(:)
!         real(prec), allocatable :: u(:), v(:), w(:)
!         real(prec), allocatable :: thk(:)            ! Tracer thickness (for compression)
!         real(prec), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
!         end type 

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

!         type tracer_state_class 
!         logical, allocatable :: is_active(:) 
!         real(prec), allocatable :: x(:), y(:), z(:)
!         real(prec), allocatable :: u(:), v(:), w(:)
!         real(prec), allocatable :: thk(:)            ! Tracer thickness (for compression)
!         real(prec), allocatable :: T(:)              ! Current temperature of the tracer (for borehole comparison, internal melting...)
!         end type 
        
        ! Determine which timestep this is
        call nc_dims(path_out,"time",dims=dims)
        nt = dims(1)
        call nc_read(path_out,"time",time_in,start=[nt],count=[1])
        if (time_in .ne. MV .and. time .gt. time_in) nt = nt+1 
        
        call nc_write(path_out,"time",time,dim1="time",start=[nt],count=[1],missing_value=MV)
        call nc_write(path_out,"x",trc%now%x,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"y",trc%now%y,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"z",trc%now%z,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"u",trc%now%u,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"v",trc%now%v,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"w",trc%now%w,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"thk",trc%now%thk,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"T",trc%now%T,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])


        return 

    end subroutine tracer_write 

    subroutine tracer_read()

        implicit none 



        return 

    end subroutine tracer_read


end module tracer 


