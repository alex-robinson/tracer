
module tracer2D 
    ! Module to wrap a 2D (profile) version of the tracer model
    ! Makes calls to main tracer code by reshaping profile into
    ! 3D array with y-dimension thickness of 1. 

    use tracer_precision
    use tracer
    use ncio    

    implicit none 

    private 
    public :: tracer2D_init 
    public :: tracer2D_update 
    public :: tracer2D_end 

    public :: tracer2D_write_init
    public :: tracer2D_write
    public :: tracer2D_write_stats

contains


    subroutine tracer2D_init(trc,filename,time,x,is_sigma)

        implicit none 

        type(tracer_class),   intent(OUT) :: trc 
        character(len=*),     intent(IN)  :: filename 
        real(prec), intent(IN) :: x(:)
        logical,    intent(IN) :: is_sigma 
        real(4) :: time 

        real(prec) :: y(2) 

        ! Define the ghost y-dimension
        y(1:2) = [0.0,1.0] 

        ! Call 3D tracer_init
        call tracer_init(trc,filename,time,x,y,is_sigma)

        return 

    end subroutine tracer2D_init


    subroutine tracer2D_update(trc,time,x,z,z_srf,H,ux,uz,dep_now,stats_now)

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), z(:)
        real(prec), intent(IN) :: z_srf(:), H(:)
        real(prec), intent(IN) :: ux(:,:), uz(:,:)
        logical,    intent(IN) :: dep_now, stats_now 
         
        ! Local variables
        real(prec) :: y(2) 
        real(prec), allocatable :: z_srf_2D(:,:), H_2D(:,:)
        real(prec), allocatable :: ux_3D(:,:,:), uy_3D(:,:,:), uz_3D(:,:,:)
        integer :: j 

        ! Define ghost dimension and data 
        allocate(z_srf_2D(size(x,1),2))
        allocate(H_2D(size(x,1),2))
        allocate(ux_3D(size(ux,1),2,size(ux,2)))
        allocate(uy_3D(size(ux,1),2,size(ux,2)))
        allocate(uz_3D(size(ux,1),2,size(ux,2)))

        ! Set y-dimension to one value of zero
        y(1:2) = [0.0,1.0] 

        ! Reshape input data with a ghost y-dimension of length two
        do j = 1, 2
            z_srf_2D(:,j) = z_srf 
            H_2D(:,j)     = H 
            ux_3D(:,j,:)  = ux 
            uy_3D         = 0.0 
            uz_3D(:,j,:)  = uz 
        end do 

        ! Now update tracers using 3D call 
        call tracer_update(trc,time,x,y,z,z_srf_2D,H_2D,ux_3D,uy_3D,uz_3D,dep_now,stats_now)

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

        call nc_write(path_out,"id",trc%now%id,dim1="pt",dim2="time", missing_value=int(MV), &
                        start=[1,nt],count=[trc%par%n ,1])

        ! Write deposition information
        call nc_write(path_out,"dep_time",trc%dep%time,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])
        call nc_write(path_out,"dep_H",trc%dep%H,dim1="pt",dim2="time", missing_value=MV, &
                        start=[1,nt],count=[trc%par%n ,1])

        return 

    end subroutine tracer2D_write 

    subroutine tracer2D_write_stats(trc,time,fldr,filename) !,z_srf,H)

        implicit none 

        type(tracer_class), intent(IN) :: trc 
        real(prec) :: time
        character(len=*), intent(IN)   :: fldr, filename 
!         real(prec),         intent(IN) :: z_srf(:), H(:) 

        ! Local variables 
        integer :: nt 
        character(len=512) :: path_out 

        path_out = trim(fldr)//"/"//trim(filename)

        ! Create output file 
        call nc_create(path_out)
        call nc_write_dim(path_out,"xc",        x=trc%stats%x*1e-3,     units="km")
        call nc_write_dim(path_out,"depth_norm",x=trc%stats%depth_norm, units="1")
        call nc_write_dim(path_out,"age_iso",   x=trc%stats%age_iso,    units="ka")
        call nc_write_dim(path_out,"time",      x=time,unlimited=.TRUE.,units="ka")
        
!         call nc_write(path_out,"z_srf",z_srf,dim1="xc",missing_value=MV, &
!                       units="m",long_name="Surface elevation")
!         call nc_write(path_out,"H",H,dim1="xc",missing_value=MV, &
!                       units="m",long_name="Ice thickness")

        call nc_write(path_out,"ice_age",trc%stats%ice_age(:,1,:),dim1="xc",dim2="depth_norm",missing_value=MV, &
                      units="ka",long_name="Layer age")
        call nc_write(path_out,"ice_age_err",trc%stats%ice_age_err(:,1,:),dim1="xc",dim2="depth_norm",missing_value=MV, &
                      units="ka",long_name="Layer age - error")

        call nc_write(path_out,"depth_iso",trc%stats%depth_iso(:,1,:),dim1="xc",dim2="age_iso",missing_value=MV, &
                      units="ka",long_name="Isochrone depth")
        call nc_write(path_out,"depth_iso_err",trc%stats%depth_iso_err(:,1,:),dim1="xc",dim2="age_iso",missing_value=MV, &
                      units="ka",long_name="Isochrone depth - error")
        
        call nc_write(path_out,"density",trc%stats%density(:,1,:),dim1="xc",dim2="depth_norm",missing_value=int(MV), &
                      units="1",long_name="Tracer density")
        
        return 

    end subroutine tracer2D_write_stats





end module tracer2D

