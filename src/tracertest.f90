
program tracertest 

    use ncio 
    use tracer 

    implicit none 

    type(tracer_class) :: trc1
    character(len=128) :: file0, fldr, filename 

    integer :: nx, ny, nz
    integer, allocatable :: dims(:) 
    real(4), allocatable :: sigma(:), xc(:), yc(:)
    real(4), allocatable :: zs(:,:), zb(:,:), H(:,:), ux(:,:,:), uy(:,:,:), uz(:,:,:)

    integer :: k, kmax, q 
    real(4) :: time, time_end, dt  

    file0 = "data/trace1_relax/Grisli15_3D.nc"
    call nc_dims(file0,"Ux",dims=dims)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

    allocate(xc(nx))
    allocate(yc(ny))
    allocate(sigma(nz))
    allocate(zs(nx,ny))
    allocate(zb(nx,ny))
    allocate(H(nx,ny))
    
    allocate(ux(nx,ny,nz))
    allocate(uy(nx,ny,nz))
    allocate(uz(nx,ny,nz))

    call nc_read(file0,"xc",xc)
    call nc_read(file0,"yc",yc)
    call nc_read(file0,"sigma",sigma)

    xc = xc*1e3 
    yc = yc*1e3 

    call nc_read(file0,"z_srf",zs)
    call nc_read(file0,"z_bed",zb)
    call nc_read(file0,"H",H)

    call nc_read(file0,"Ux",ux)
    call nc_read(file0,"Uy",uy)
    call nc_read(file0,"Uz",uz)

    fldr     = "output"
    filename = "GRL-20KM_trc1.nc"

    ! Test tracer_update
    time     = 0.0 
    time_end = 1002.0
    dt       = 1.0 

    ! Initialize tracer and output file 
    call tracer_init(trc1,time=time)
    call tracer_write_init(trc1,fldr,filename)

    q = 9 

    do k = 1, int(time_end/dt), int(dt) 

        if (k .gt. 1) time = time + dt 
        write(*,*) "time = ", time, trc1%par%n_active

        call tracer_update(trc1%par,trc1%now,trc1%dep,time=time, &
                            x=xc,y=yc,z=sigma,z_srf=zs,H=H,ux=ux,uy=uy,uz=uz)

        q = q+1 
        if (q==10) then 
            call tracer_write(trc1,time,fldr,filename)
            q = 0 
        end if 

    end do 

contains 



end program tracertest 

