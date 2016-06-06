
program tracertest 

    use ncio 
    use tracer 
    use coordinates 

    implicit none 

    type(tracer_class) :: trc1
    type(grid_class)   :: grid 
    character(len=128) :: file0, fldr, filename 

    integer :: nz 
    real(4), allocatable :: zs(:,:), H(:,:), u(:,:,:), v(:,:,:), w(:,:,:)

    integer :: k, kmax 
    real(4) :: time 

    ! Initialize a grid for the tracer domain
    call grid_init(grid,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                   lon180=.TRUE.,dx=20.d0,nx=90,dy=20.d0,ny=150, &
                   lambda=-40.d0,phi=72.d0,alpha=8.4d0)


    ! Generate boundary variables 
    call grid_allocate(grid,zs)
    call grid_allocate(grid,H)

    nz = 10 
    allocate(u(size(H,1),size(H,2),nz))
    allocate(v(size(H,1),size(H,2),nz))
    allocate(w(size(H,1),size(H,2),nz))

    ! Load topography data 
    file0 = "data/"//trim(grid%name)//"_TOPO-B13_gl0.05.nc"
    call nc_read(file0,"zs",zs)
    call nc_read(file0,"H",H)



    fldr = "output"
    filename = trim(grid%name)//"_trc1.nc"

    ! Initialize tracer and output file 
    call tracer_init(trc1)
    call tracer_write_init(trc1,fldr,filename)

    
    ! Test tracer_update
    time = 0.0 
    kmax = 100
    do k = 1, kmax 

        time = time + real(k)

        call tracer_update(trc1%par,trc1%now,trc1%dep,time=0.0, &
                            x=real(grid%G%x),y=real(grid%G%y),z=real(grid%G%x*0.0), &
                            z_srf=real(zs),H=real(H))

        call tracer_write(trc1,time,fldr,filename)

    end do 

contains 



end program tracertest 

