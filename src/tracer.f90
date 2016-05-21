
module tracer 

    use ncio 
    use nml 

    implicit none 

    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 

    ! Precision used here
    integer,  parameter :: prec = sp 

    type tracer_par_class 
        integer :: n, n_active
        real(prec) :: time_now, time_old, dt  
    end type 

    type tracer_state_class 
        logical, allocatable :: is_active(:) 
        real(prec), allocatable :: x(:), y(:), z(:)
        real(prec), allocatable :: u(:), v(:), w(:)
        real(prec), allocatable :: time_dep(:) 

    end type 

    type tracer_class 
        type(tracer_par_class)   :: par 
        type(tracer_state_class) :: now 

    end type 

    private 
    public :: tracer_class 
    public :: tracer_init 
    public :: tracer_update 

contains 

    subroutine tracer_init(par,now)

        implicit none 

        type(tracer_par_class),   intent(OUT) :: par 
        type(tracer_state_class), intent(OUT) :: now 

        ! Load the parameters
        call tracer_par_load(par)

        ! Allocate the state variables 
        call tracer_allocate(now,n=par%n)

        ! Initialize state 
        now%is_active = .FALSE. 

        now%x         = 0.0 
        now%y         = 0.0 
        now%z         = 0.0 
        now%u         = 0.0 
        now%v         = 0.0 
        now%w         = 0.0 
        now%time_dep  = 0.0 

        return 

    end subroutine tracer_init

    subroutine tracer_update(par,now,time,x0,y0,z0,u0,v0,w0)

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x0(:,:,:), y0(:,:,:), z0(:,:,:)
        real(prec), intent(IN) :: u0(:,:,:), v0(:,:,:), w0(:,:,:)
        
        ! Update current time 
        par%time_old = par%time_now 
        par%time_now = time 
        par%dt       = par%time_now - par%time_old 

        ! Deposit new tracer points

        ! == TO DO == 
        ! - Function based on H and U for location, input par deposition frequency 
        ! - Attach whatever information we want to trace (age, climate, isotopes, etc)

        ! Interpolate velocities to active point locations 

        ! == TO DO == 

        ! Update the tracer thickness, then destroy points that are too thin 

        ! == TO DO == 

        ! Update the tracer positions 
        call calc_position(now%x,now%y,now%z,now%u,now%v,now%w,par%dt,now%is_active)

        ! Destroy points that moved outside the valid region 

        ! == TO DO == 



        return 

    end subroutine tracer_update

    ! ================================================
    !
    ! tracer management 
    !
    ! ================================================
    
    subroutine tracer_activate()
        ! Use this to activate individual or multiple tracers
        implicit none 


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

        return 

    end subroutine tracer_par_load

    subroutine tracer_allocate(now,n)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: n

        ! Make object is deallocated
        call tracer_deallocate(now)

        ! Allocate tracer 
        allocate(now%is_active(n))
        allocate(now%x(n),now%y(n),now%z(n))
        allocate(now%u(n),now%v(n),now%w(n))

        ! Allocate tracer properties 
        allocate(now%time_dep(n))

        return

    end subroutine tracer_allocate

    subroutine tracer_deallocate(now)

        implicit none 

        type(tracer_state_class), intent(INOUT) :: now 

        ! Allocate state objects
        if (allocated(now%is_active)) deallocate(now%is_active)
        if (allocated(now%x))         deallocate(now%x)
        if (allocated(now%y))         deallocate(now%y)
        if (allocated(now%z))         deallocate(now%z)
        if (allocated(now%u))         deallocate(now%u)
        if (allocated(now%v))         deallocate(now%v)
        if (allocated(now%w))         deallocate(now%w)

        if (allocated(now%time_dep))  deallocate(now%time_dep)

        return

    end subroutine tracer_deallocate

end module tracer 


