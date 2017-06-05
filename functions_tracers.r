filt_tracers = function(trc)
{
    
    kk = which(!is.na(trc$id))
    dat = data.frame(id=unique(as.vector(trc$id[kk])),time=NA,x=NA,z=NA,dpth=NA,z_srf=NA,H=NA,
                     ux=NA,uz=NA,dep_time=NA,dep_H=NA,dep_x=NA,dep_z=NA)

    time = t(array(trc$time,dim=rev(dim(trc$id))))

    nms = names(dat)

    for (i in 1:length(dat$id)) { 
        kk = which(trc$id==dat$id[i])
        k  = kk[which.max(time[kk])]

        for (q in 1:length(nms)) {
            nm = nms[q]
            dat[[nm]][i] = trc[[nm]][k]
        }

    }

    return(dat)
}

calc_rh_age = function(sigma,H0,G)
{
    # Calculate analytical age at the divide
    age    = -(H0/G)*log(sigma)
    age[!is.finite(age)] = NA 
    return(age)
}

calc_profile_RH2003 = function(x=NULL,z=NULL)
{
        # Define a 2D profile (x-z) following 
        # Rybak and Huybrechts (2003, Annals of Glaciology)

        if (is.null(x)) x = seq(0,1000,length.out=51)
        if (is.null(z)) z = seq(0,1,length.out=101)

        # Local parameters 
        nx  = length(x)
        nz  = length(z)
        ng  = 3           # exponent
        rho = 910.0       # kg/m^3
        gg  = 9.81        # m/s
        A   = 10.0^(-16)  # Pa^3/a
        L   = 10.0^6      # m 
        G   = 0.10        # m/a
        B   = 0.00 

        M   = G - B 

        xc = x*1e3   # [km] => [m] 
        sigma = z 

        # Calculate H0 (should be H0=3598.4 m)
        H0 = (20.0*M/A)^(1.0/(2.0*(ng+1)))*(1/(rho*gg))^(ng/(2.0*(ng+1)))*L^(1.0/2.0)
        
        prof = list(x=x,z=z,xc=xc,sigma=sigma,H0=H0,G=G)

        # Calculate H(x) 
        prof$H = H0 * (1.0-(abs(prof$xc)/L)^((ng+1.0)/ng))^(ng/(2.0*(ng+1.0)))

        prof$ux   = array(NA,dim=c(nx,nz)) 
        prof$uz   = array(NA,dim=c(nx,nz))
        prof$age  = rep(NA,nz)
        prof$dHdx = rep(0.0,nx)
        prof$xx   = array(0.0,dim=c(nx,nz))
        prof$zz   = array(0.0,dim=c(nx,nz))

        # Calculate slope (central differences)
        prof$dHdx[1] = (prof$H[i+1]-prof$H[i])/(prof$xc[i+2]-prof$xc[i])
        for (i in 2:(nx-1)) {
            prof$dHdx[i] = (prof$H[i+1]-prof$H[i-1])/(prof$xc[i+1]-prof$xc[i-1])
        }

        # Calculate velocities 
        for (i in 1:nx) {

            H    = prof$H[i] 
            dHdx = prof$dHdx[i]

            for (j in 1:nz) {

                if (H > 0) {
                    prof$ux[i,j] = -(2.0*A)/(ng+1.0)*(rho*gg)^ng * abs(dHdx)^(ng-1.0) *
                      dHdx * (H^(ng+1.0)-(H-prof$sigma[j]*H)^(ng+1.0))

                    prof$uz[i,j] = prof$sigma[j]*(-G+B+prof$ux[i,j]*dHdx) - B
                }

                prof$xx[i,j] = prof$xc[i]*1e-3
                prof$zz[i,j] = prof$sigma[j]*H

            }
        } 

        prof$umag = sqrt((prof$ux)^2 + prof$uz^2)

        # Calculate analytical age at the divide
        prof$age = calc_rh_age(prof$sigma,H0,G)*1e-3

        return(prof)
}

load_tracer_profile = function(filename,rh=NULL)
{
    # Load data 
    trc = my.read.nc(filename)

    # Get sigma axis 
    trc$sigma = trc$z / trc$H

    # Convert time to ka 
    trc$time     = trc$time*1e-3
    trc$dep_time = trc$dep_time*1e-3
    trc$age      = trc$age*1e-3 

    # Calculate analytical age and error
    if (!is.null(rh)) {
        trc$rh_age    = calc_rh_age(sigma=trc$sigma,H0=rh$H0,G=rh$G)*1e-3
        trc$age_err   = (trc$age - trc$rh_age)
        trc$age_err_p = 100* (trc$age - trc$rh_age) / trc$rh_age
        trc$age_err_p[abs(trc$age)<1e-1] = NA 
    }

    trc$umag = sqrt(trc$ux^2+trc$uz^2)

    return(trc)
}

calc_tracer_err = function(trc,ref)
{



    return(dat)
}