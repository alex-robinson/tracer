
module tracer_import 

    use tracer_precision
    use tracer_interp 
    use tracer 

    implicit none 





contains 


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
        real(prec_hi), allocatable :: x1(:), y1(:), z1(:)
        real(prec_hi), allocatable :: z_srf1(:,:), H1(:,:)
        real(prec_hi), allocatable :: age1(:,:,:) 
        real(prec_hi) :: zc(size(z))
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



end module tracer_import
