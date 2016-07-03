
module tracer2D 
    ! Module to wrap a 2D (profile) version of the tracer model
    ! Makes calls to main tracer code by reshaping profile into
    ! 3D array with y-dimension thickness of 1. 

    use tracer
    use ncio   
    use nml 

    implicit none 

    private 
    public :: tracer2D_init 
    public :: tracer2D_update 
    public :: tracer2D_end 

    public :: tracer2D_write_init
    public :: tracer2D_write
    public :: tracer2D_write_stats

contains


    subroutine tracer2D_init(trc,time,x,z)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        real(prec), intent(IN) :: x(:), z(:)
        real(4) :: time 

        real(prec) :: y(1) 

        ! Define the ghost y-dimension
        y(1) = 0.0 

        ! Call 3D tracer_init
        call tracer_init(trc,time,x,y,z)

        return 

    end subroutine tracer2D_init


    subroutine tracer2D_update(par,now,dep,stats,time,x,z,z_srf,H,ux,uz)

        implicit none 

        type(tracer_par_class),   intent(INOUT) :: par 
        type(tracer_state_class), intent(INOUT) :: now 
        type(tracer_dep_class),   intent(INOUT) :: dep
        type(tracer_stats_class), intent(INOUT) :: stats
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), z(:)
        real(prec), intent(IN) :: z_srf(:), H(:)
        real(prec), intent(IN) :: ux(:,:), uz(:,:)
        
        ! Local variables
        real(prec) :: y(1) 
        real(prec), allocatable :: z_srf_2D(:,:), H_2D(:,:)
        real(prec), allocatable :: ux_3D(:,:,:), uy_3D(:,:,:), uz_3D(:,:,:)

        ! Define ghost dimension and data 
        allocate(z_srf_2D(size(x,1),1))
        allocate(H_2D(size(x,1),1))
        allocate(ux_3D(size(ux,1),1,size(ux,2)))
        allocate(uy_3D(size(ux,1),1,size(ux,2)))
        allocate(uz_3D(size(ux,1),1,size(ux,2)))

        ! Set y-dimension to one value of zero
        y(1) = 0.0 

        ! Reshape input data with a ghost y-dimension of length one
        z_srf_2D(:,1) = z_srf 
        H_2D(:,1)     = H 
        ux_3D(:,1,:)  = ux 
        uy            = 0.0 
        uz_3D(:,1,:)  = uz 

        ! Now update tracers using 3D call 
        call tracer_update(par,now,dep,stats,time,x,y,z,z_srf,H,ux,uy,uz)

        return 

    end subroutine tracer2D_update

    subroutine tracer2D_end(trc)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        
        ! Call normal tracer_end subroutine 
        call tracer_end(trc) 

        return 

    end subroutine tracer2D_end



    ! ================================================
    !
    ! I/O routines 
    !
    ! ================================================

    subroutine tracer2D_write_init(trc,fldr,filename)

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

    end subroutine tracer2D_write_init 

    subroutine tracer2D_write(trc,time,fldr,filename)

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
        call nc_write(path_out,"z",trc%now%z,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"ux",trc%now%ux,dim1="pt",dim2="time", missing_value=MV, &
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

    end subroutine tracer2D_write 

    subroutine tracer2D_write_stats(trc,time,fldr,filename)

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
        call nc_write_dim(path_out,"sigma",x=trc%stats%z)
        call nc_write_dim(path_out,"time",x=time,unlimited=.TRUE.)

!         call nc_write(path_out,"density",trc%stats%density,dim1="xc",dim2="yc",dim3="sigma",missing_value=int(MV), &
!                       units="1",long_name="Tracer density (surface)")
        call nc_write(path_out,"dens_srf",trc%stats%density(:,1,1),dim1="xc",missing_value=int(MV), &
                      units="1",long_name="Tracer density (surface)")

        return 

    end subroutine tracer2D_write_stats





end module tracer2D

