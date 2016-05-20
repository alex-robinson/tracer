
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

    end type 

    type tracer_state_class 
        logical, allocatable :: is_active(:) 
        real(prec), allocatable :: x(:), y(:), z(:)
        real(prec), allocatable :: u(:), v(:), w(:)
        
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


        return 

    end subroutine tracer_init


    subroutine tracer_update(par,now,time)

        implicit none 

        type(tracer_par_class), intent(IN)      :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        real(prec), intent(IN) :: time 

        
        return 

    end subroutine tracer_update


    ! ================================================
    !
    ! tracer physics 
    !
    ! ================================================
    
    elemental subroutine calc_position(x,y,z,u,v,w,dt)

        implicit none 

        real(prec), intent(INOUT) :: x, y, z 
        real(prec), intent(IN)    :: u, v, w 
        real(prec), intent(IN)    :: dt 

        x = x + u*dt 
        y = y + v*dt 
        z = z + w*dt 

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

        return

    end subroutine tracer_deallocate

end module tracer 


