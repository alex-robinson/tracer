
module tracer_interp 

    use tracer_precision 

    implicit none 


    type bilin_par_type 
        integer    :: i, j  
        real(prec) :: alpha1, alpha2 
    end type 

    real(prec) :: missing_value = MISSING_VALUE_DEFAULT 
    real(prec) :: mv = MISSING_VALUE_DEFAULT 

    private
    public :: interp_bilinear_weights 


contains 


    function interp_bilinear_weights(x,y,xout,yout) result(par)

        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Save the indices and weights for interpolation of variables later

        implicit none 

        real(prec), intent(IN) :: x(:), y(:) 
        real(prec), intent(IN) :: xout, yout
        type(bilin_par_type)   :: par 

        integer :: x_idx, y_idx
        integer :: i, j, i1, j1  
        integer :: nx, ny, nx1, ny1 
        real(prec) :: alpha1, alpha2, p0, p1 

        nx = size(x,1)
        ny = size(y,1)

        ! Get x-index corresponding to nearest neighbor
        ! greater-than-equal-to x-value of interest
        if (xout .le. x(1)) then 
            x_idx = -1

        else if (xout .ge. (x(nx))) then 
            x_idx = -2 
        else 

            do i = 1, nx 
                if (x(i) .ge. xout) exit 
            end do 

            x_idx = i
        end if 


        ! Get y-index corresponding to nearest neighbor
        ! greater-than-equal-to y-value of interest
        if (yout .le. y(1)) then 
            y_idx = -1

        else if (yout .ge. (y(ny))) then 
            y_idx = -2 
        else 

            do j = 1, ny 
                if (y(j) .ge. yout) exit 
            end do 

            y_idx = j 
        end if 

        ! Now get bilinear interpolation weights
        i = x_idx
        j = y_idx 

        ! Only interpolate points inside the original grid
        if (i .gt. 0 .and. i-1 .gt. 0 .and. j .gt. 0 .and. j-1 .gt. 0) then 
            
            alpha1 = (xout - x(i-1)) / (x(i)-x(i-1))
            alpha2 = (yout - y(j-1)) / (y(j)-y(j-1))
            
        else 
            alpha1 = 0.0 
            alpha2 = 0.0 

        end if 

        par%i      = x_idx 
        par%j      = y_idx 
        par%alpha1 = alpha1 
        par%alpha2 = alpha2 

        return 

    end function interp_bilinear_weights

    function interp_bilinear(par,z) result (zout)

        implicit none 

        type(bilin_par_type), intent(IN) :: par
        real(prec), intent(IN) :: z(:,:) 
        real(prec) :: zout 
        real(prec) :: p0, p1 

        if (par%i .gt. 0) then 
            p0 = z(par%i-1,par%j-1) + par%alpha1*(z(par%i,par%j-1)-z(par%i-1,par%j-1))
            p1 = z(par%i-1,par%j)   + par%alpha1*(z(par%i,par%j)-z(par%i-1,par%j))
            zout = p0 + par%alpha2*(p1-p0)
        else 
            zout = missing_value 
        end if 

        return 

    end function interp_bilinear

end module tracer_interp 
