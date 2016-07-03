
program tracertest 

    use ncio
    use tracer 
    use tracer2D 

    implicit none 

    type(tracer_class) :: trc1
    character(len=128) :: file0, fldr, filename, filename_stats 

    integer :: nx, nz
    integer, allocatable :: dims(:) 
    real(4), allocatable :: xc(:), sigma(:)
    real(4), allocatable :: zs(:), zb(:), H(:), ux(:,:), uz(:,:)

    integer :: k, kmax, q 
    real(4) :: time, time_end, dt  

    file0 = "data/trace1_relax/Grisli15_3D.nc"
    call nc_dims(file0,"Ux",dims=dims)
    nx = dims(1)
    nz = dims(3)

    allocate(xc(nx))
    allocate(sigma(nz))
    allocate(zs(nx))
    allocate(zb(nx))
    allocate(H(nx))
    
    allocate(ux(nx,nz))
    allocate(uz(nx,nz))

    call nc_read(file0,"xc",xc)
    call nc_read(file0,"sigma",sigma)

    xc = xc*1e3 

    call nc_read(file0,"z_srf",zs)
    call nc_read(file0,"z_bed",zb)
    call nc_read(file0,"H",H)

    call nc_read(file0,"Ux",ux)
    call nc_read(file0,"Uz",uz)

    fldr     = "output"
    filename       = "GRL-20KM_trc1.nc"
    filename_stats = "GRL-20KM_trc1-stats.nc"

    ! Test tracer_update
    time     = 0.0 
    time_end = 1002.0
    dt       = 1.0 

    ! Initialize tracer and output file 
    call tracer2D_init(trc1,time=time,x=xc,z=sigma)
    call tracer2D_write_init(trc1,fldr,filename)

    q = 9 

    do k = 1, int(time_end/dt), int(dt) 

        if (k .gt. 1) time = time + dt 
        write(*,*) "time = ", time, trc1%par%n_active

        call tracer2D_update(trc1%par,trc1%now,trc1%dep,trc1%stats,time=time, &
                             x=xc,z=sigma,z_srf=zs,H=H,ux=ux,uz=uz)

        q = q+1 
        if (q==10) then 
            call tracer2D_write(trc1,time,fldr,filename)
            q = 0 
        end if 

    end do 

    ! Write stats 
    call tracer2D_write_stats(trc1,time,fldr,filename_stats)

contains 



end program tracertest 

