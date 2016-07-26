
program tracertest 

    use ncio
    use tracer 
    use tracer2D 
    use tracer_precision 

    implicit none 

    type(tracer_class) :: trc1
    character(len=128) :: fldr, filename, filename_stats 

    type profile_class 
        integer :: nx, nz
        integer, allocatable :: dims(:) 
        real(4), allocatable :: xc(:), sigma(:)
        real(4), allocatable :: zs(:), zb(:), H(:), dHdx(:), ux(:,:), uz(:,:)
        real(4), allocatable :: age(:)

    end type 

    type(profile_class) :: prof1 
    
    integer :: k, kmax, q 
    real(4) :: time, time_end, dt  
    integer :: i0, i1 
    real(4) :: dt_write, dt_write_now  
    logical             :: dep_now 
    real(prec)          :: dt_dep 

    
    ! Limit x-direction of domain to close to ice divide (x=0) to minimize 
    ! computations for checking comparison
    i1 = 51 

    call calc_profile_RH2003(prof1)
    call profile_write(prof1,fldr="output",filename="profile_RH2003.nc")

    fldr     = "output"
    filename       = "profile_RH2003_trc1.nc"
    filename_stats = "profile_RH2003_trc1-stats.nc"

    ! Test tracer_update
    time     = -500000.0 
    time_end = 0.0
    dt       = 5.0 

    ! Initialize tracer and output file 
    call tracer2D_init(trc1,"RH2003.nml",time=time,x=prof1%xc,z=prof1%sigma,is_sigma=.TRUE.)
    call tracer2D_write_init(trc1,fldr,filename)

    dt_write = 10000.0 
    dt_dep   = 200.0 

    dt_write_now = 0.0 

    do k = 1, int((time_end-time)/dt)+1 

        if (k .gt. 1) time = time + dt 

        dep_now  = .FALSE.
        if (mod(time,dt_dep) .eq. 0.0) dep_now = .TRUE. 

        write(*,*) "time = ", time, trc1%par%n_active


        call tracer2D_update(trc1%par,trc1%now,trc1%dep,trc1%stats,time=time, &
                             x=prof1%xc,z=prof1%sigma,z_srf=prof1%H,H=prof1%H, &
                             ux=prof1%ux*0.0,uz=prof1%uz,dep_now=dep_now)

        if (k .gt. 1) dt_write_now = dt_write_now+dt 
        if (dt_write_now .eq. 0.0 .or. dt_write_now .ge. dt_write) then 
            call tracer2D_write(trc1,time,fldr,filename)
            dt_write_now = 0.0 
        end if 

    end do 

    ! Write stats 
    call tracer2D_write_stats(trc1,time,fldr,filename_stats)

contains 

    subroutine calc_profile_RH2003(prof,i1)
        ! Define a 2D profile (x-z) following 
        ! Rybak and Huybrechts (2003, Annals of Glaciology)

        implicit none 

        type(profile_class), intent(INOUT) :: prof 
        integer, intent(IN), optional :: i1 

        ! Local parameters 
        integer, parameter :: nx = 51 
        integer, parameter :: nz = 101 
        integer, parameter :: ng = 3          ! exponent
        real(4), parameter :: rho = 910.0     ! kg/m^3
        real(4), parameter :: g = 9.81        ! m/s
        real(4), parameter :: A = 10.0**(-16) ! Pa^3/a
        real(4), parameter :: M = 0.1         ! m/a
        real(4), parameter :: L = 10.0**6     ! m 

        ! Local variables 
        integer :: i, j 
        real(4) :: H0, dHdx, H 
        real(4) :: GG, B 

        GG = M 
        B  = 0 

        prof%nz = nz
        prof%nx = nx    ! Right-hand side symmetrical of domain
!         prof%nx = nx+50 ! Whole domain 

!         if (present(i1)) prof%nx = i1

        allocate(prof%xc(prof%nx))
        allocate(prof%sigma(prof%nz))
        allocate(prof%zs(prof%nx))
        allocate(prof%zb(prof%nx))
        allocate(prof%H(prof%nx))
        allocate(prof%dHdx(prof%nx)) 

        allocate(prof%ux(prof%nx,prof%nz))
        allocate(prof%uz(prof%nx,prof%nz))

        allocate(prof%age(prof%nz)) 

        if (prof%nx .lt. 100) then
            ! Define x-dimension (0-1000 km) 
            do i = 1, prof%nx 
                prof%xc(i) = 0.0 + (i-1)/real(nx-1) * 1000.0 
            end do
        else
            ! Define x-dimension (1000-1000 km)
            do i = 1, prof%nx 
                prof%xc(i) = -1000.0 + (i-1)/real(nx-1) * 1000.0 
            end do 
        end if 

        prof%xc = prof%xc*1e3   ! [km] => [m] 

        ! Define z-dimension (0-1 sigma)
        do i = 1, prof%nz 
            prof%sigma(i) = 0.0 + (i-1)/real(prof%nz-1) * 1.0 
        end do 
        prof%sigma(1) = 1e-6   ! To avoid singularities 

        ! Calculate H0 (should be H0=3598.4 m)
        H0 = (20.0*M/A)**(1.0/(2.0*(ng+1)))*(1/(rho*g))**(ng/(2.0*(ng+1)))*L**(1.0/2.0)
        
        ! Calculate H(x) 
        prof%H = H0 * (1.0-(abs(prof%xc)/L)**((ng+1.0)/ng))**(real(ng)/(2.0*(ng+1.0)))

        prof%ux   = 0.0 
        prof%uz   = 0.0 
        prof%age  = 0.0 
        prof%dHdx = 0.0

        ! Calculate velocities 
        do i = 1, prof%nx-1
            dHdx = (prof%H(i+1)-prof%H(i))/(prof%xc(i+1)-prof%xc(i))
!             dHdx = (prof%H(i)-prof%H(i-1))/(prof%xc(i)-prof%xc(i-1))
            prof%dHdx(i) = dHdx

            H    = prof%H(i) 
            do j = 1, prof%nz

                prof%ux(i,j) = -(2.0*A)/(ng+1.0)*(rho*g)**ng * dHdx**(ng-1.0) &
                  * dHdx * (H**(ng+1.0)-(H-prof%sigma(j)*H)**(ng+1.0))

                prof%uz(i,j) = prof%sigma(j)*(-GG+B+prof%ux(i,j)*dHdx) - B

                
            end do 
        end do 

        ! Calculate analytical age at the divide
        prof%age = (H0/GG)*log(prof%sigma)
        prof%age(1) = prof%age(2)
        
!         ! Write summary 
!         write(*,"(a,500f8.2)") "xc: ",    prof%xc*1e-3 
!         write(*,"(a,500f8.2)") "sigma: ", prof%sigma 
!         write(*,"(a,f10.2)") "H0 = ", H0 
!         write(*,"(a,2f10.2)") "range(ux): ", minval(prof%ux), maxval(prof%ux)
!         write(*,"(a,2f10.2)") "range(uz): ", minval(prof%uz), maxval(prof%uz) 

        return 

    end subroutine calc_profile_RH2003

    subroutine profile_write(prof,fldr,filename)

        implicit none 

        type(profile_class), intent(IN) :: prof 
        character(len=*), intent(IN)    :: fldr, filename 

        ! Local variables 
        character(len=512) :: path_out 

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",x=prof%xc*1e-3,units="kilometers")
        call nc_write_dim(path_out,"sigma",x=prof%sigma,units="1")

        call nc_write(path_out,"H",   prof%H,   dim1="xc",missing_value=MV)
        call nc_write(path_out,"dHdx",prof%dHdx,dim1="xc",missing_value=MV)
        call nc_write(path_out,"ux",prof%ux,dim1="xc",dim2="sigma",missing_value=MV)
        call nc_write(path_out,"uz",prof%uz,dim1="xc",dim2="sigma",missing_value=MV)
        call nc_write(path_out,"umag",sqrt(prof%ux**2+prof%uz**2),dim1="xc",dim2="sigma",missing_value=MV)
        call nc_write(path_out,"age",prof%age,dim1="sigma",missing_value=MV)
        
        return 

    end subroutine profile_write 

end program tracertest 

