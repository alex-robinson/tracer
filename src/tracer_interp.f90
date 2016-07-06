
module tracer_interp 

    use tracer_precision 

    implicit none 

    type lin_interp_par_type 
        integer    :: i
        real(prec) :: alpha
    end type 

    type lin3_interp_par_type 
        integer    :: i, j, k 
        real(prec) :: alpha_x, alpha_y, alpha_z
    end type 

    private

    public :: lin_interp_par_type
    public :: lin3_interp_par_type
    public :: interp_bilinear_weights, interp_bilinear 
    public :: interp_trilinear_weights, interp_trilinear  

contains 

    subroutine calc_interp_linear_weights(idx,alpha,x,xout)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        implicit none 

        integer,    intent(OUT) :: idx 
        real(prec), intent(OUT) :: alpha 
        real(prec), intent(IN)  :: x(:)
        real(prec), intent(IN)  :: xout

        integer :: i, nx

        nx = size(x,1)

        ! By default, set index out of range and weight to zero
        idx   = -1 
        alpha = 0.0 

        ! Get x-index corresponding to nearest neighbor
        ! greater-than-equal-to x-value of interest (inside range)
        if (xout .gt. x(1) .and. xout .le. x(nx)) then 

            do i = 1, nx 
                if (x(i) .ge. xout) exit 
            end do 

            idx   = i 
            alpha = (xout - x(i-1)) / (x(i)-x(i-1))

        else if (xout .eq. x(1)) then 

            idx   = 2 
            alpha = 0.0 

        end if

        return 

    end subroutine calc_interp_linear_weights

    function interp_bilinear_weights(x,y,xout,yout) result(par)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        implicit none 

        real(prec), intent(IN) :: x(:), y(:) 
        real(prec), intent(IN) :: xout, yout
        type(lin3_interp_par_type)   :: par 

        call calc_interp_linear_weights(par%i,par%alpha_x,x,xout)
        call calc_interp_linear_weights(par%j,par%alpha_y,y,yout)
        
        return 

    end function interp_bilinear_weights

    function interp_trilinear_weights(x,y,z,xout,yout,zout) result(par)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        implicit none 

        real(prec), intent(IN) :: x(:), y(:), z(:) 
        real(prec), intent(IN) :: xout, yout, zout
        type(lin3_interp_par_type)   :: par 

        call calc_interp_linear_weights(par%i,par%alpha_x,x,xout)
        call calc_interp_linear_weights(par%j,par%alpha_y,y,yout)
        call calc_interp_linear_weights(par%k,par%alpha_z,z,zout)
        
        return 

    end function interp_trilinear_weights

    function interp_bilinear(par,var) result (varout)

        implicit none 

        type(lin3_interp_par_type), intent(IN) :: par
        real(prec), intent(IN) :: var(:,:) 
        real(prec) :: varout 
        real(prec) :: p1, p2
        integer    :: i, j

        varout = MV 

        if (par%i .gt. 0 .and. par%j .gt. 0) then

            i = par%i 
            j = par%j 

            ! Lower z-plane 
            p1 = var(i-1,j-1) + par%alpha_x*(var(i,j-1)-var(i-1,j-1))
            p2 = var(i-1,j)   + par%alpha_x*(var(i,j)-var(i-1,j))
            varout = p1 + par%alpha_y*(p2-p1)

        end if 

        return 

    end function interp_bilinear

    function interp_trilinear(par,var) result (varout)

        implicit none 

        type(lin3_interp_par_type), intent(IN) :: par
        real(prec), intent(IN) :: var(:,:,:) 
        real(prec) :: varout 
        real(prec) :: p1, p2, v1, v2 
        integer    :: i, j, k 

        varout = MV 

        if (par%i .gt. 0 .and. par%j .gt. 0 .and. par%k .gt. 0) then

            i = par%i 
            j = par%j 
            k = par%k 

            ! Lower z-plane 
            p1 = var(i-1,j-1,k-1) + par%alpha_x*(var(i,j-1,k-1)-var(i-1,j-1,k-1))
            p2 = var(i-1,j,k-1)   + par%alpha_x*(var(i,j,k-1)-var(i-1,j,k-1))
            v1 = p1 + par%alpha_y*(p2-p1)

            ! Upper z-plane 
            p1 = var(i-1,j-1,k) + par%alpha_x*(var(i,j-1,k)-var(i-1,j-1,k))
            p2 = var(i-1,j,k)   + par%alpha_x*(var(i,j,k)-var(i-1,j,k))
            v2 = p1 + par%alpha_y*(p2-p1)

            ! Linear z-interpolation 
            varout = v1 + par%alpha_z*(v2-v1)

        end if 

        return 

    end function interp_trilinear


end module tracer_interp 
