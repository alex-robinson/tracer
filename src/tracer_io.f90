
module tracer_io 

    use tracer_precision
    use tracer_interp 
    use tracer3D 

    implicit none 





contains 

    subroutine tracer_read(trc,filename,time)

        implicit none 

        type(tracer_class), intent(OUT) :: trc
        character(len=*),   intent(IN)  :: filename 
        real(prec_time),    intent(IN)  :: time 

        ! Local variables 
        integer :: i, k 



        return 

    end subroutine tracer_read

    subroutine tracer_align(trc_new,trc_ref,trc,dxy_max,dz_max)
        ! Interpolate tracer ages from trc to those of trc_ref 
        ! Note: for this to work well, trc should be sufficiently high
        ! resolution to minimize interpolation errors 

        implicit none 

        type(tracer_class), intent(OUT) :: trc_new 
        type(tracer_class), intent(IN) :: trc_ref, trc  
        real(prec), intent(IN) :: dxy_max, dz_max  

        ! Local variables 
        integer :: i, k  
        real(prec) :: dist_xy(trc%par%n), dist_z(trc%par%n)

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

                dist_xy = MV 
                dist_z  = MV
                where (trc%now%active .eq. 2)
                    dist_xy = sqrt( (trc_new%now%x(i)-trc%now%x)**2 &
                                  + (trc_new%now%y(i)-trc%now%y)**2)
                    dist_z  = abs(trc_new%now%z(i)-trc%now%z)
                end where 

                k = minloc(dist_xy,mask=dist_xy.ne.MV.and.dist_z.le.dz_max,dim=1)

                if (dist_xy(k) .le. dxy_max) then 
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

    subroutine tracer_import_eulerian(trc,time,x,y,z,age,z_srf,H,is_sigma,order,sigma_srf)
        ! Given a 3D field of Eulerian ages on a grid, 
        ! convert to tracer format. 

        implicit none 

        type(tracer_class), intent(INOUT) :: trc 
        real(prec), intent(IN) :: time 
        real(prec), intent(IN) :: x(:), y(:), z(:) 
        real(prec), intent(IN) :: z_srf(:,:), H(:,:)
        real(prec), intent(IN) :: age(:,:,:) 
        logical,    intent(IN) :: is_sigma 
        character(len=*), intent(IN), optional :: order 
        real(prec), intent(IN), optional :: sigma_srf     ! Value at surface by default (1 or 0?)

        ! Local variables  
        character(len=3) :: idx_order 
        integer :: nx, ny, nz 
        real(prec), allocatable :: x1(:), y1(:), z1(:)
        real(prec), allocatable :: z_srf1(:,:), H1(:,:)
        real(prec), allocatable :: age1(:,:,:) 
        real(prec) :: zc(size(z))
        logical :: rev_z 

        ! Determine order of indices (default ijk)
        idx_order = "ijk"
        if (present(order)) idx_order = trim(order)

        ! Correct the sigma values if necessary,
        ! so that sigma==0 [base]; sigma==1 [surface]
        zc = z 
        if (trc%par%is_sigma .and. present(sigma_srf)) then 
            if (sigma_srf .eq. 0.0) then 
                ! Adjust sigma values 
                zc = 1.0 - z 
            end if 
        end if 

        ! Also determine whether z-axis is initially ascending or descending 
        rev_z = (zc(1) .gt. zc(size(zc)))

        call tracer_reshape1D_vec(x, x1,rev=.FALSE.)
        call tracer_reshape1D_vec(y, y1,rev=.FALSE.)
        call tracer_reshape1D_vec(real(zc,kind=prec),z1,rev=rev_z)
        call tracer_reshape2D_field(idx_order,z_srf,z_srf1)
        call tracer_reshape2D_field(idx_order,H,H1)
        call tracer_reshape3D_field(idx_order,age,age1,rev_z=rev_z)
        
        ! Get axis sizes (if ny==2, this is a 2D profile domain) 
        nx = size(x1,1)
        ny = size(y1,1)
        nz = size(z1,1)





        return 

    end subroutine tracer_import_eulerian



end module tracer_io
