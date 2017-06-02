
library(myr)

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

calc_rh_age = function(sigma,H0,GG)
{
    # Calculate analytical age at the divide
    age    = -(H0/GG)*log(sigma)
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
        g   = 9.81        # m/s
        A   = 10.0^(-16)  # Pa^3/a
        M   = 0.1         # m/a
        L   = 10.0^6      # m 

        GG = M 
        B  = 0 

        xc = x*1e-3   # [km] => [m] 
        sigma = z 

        # Calculate H0 (should be H0=3598.4 m)
        H0 = (20.0*M/A)^(1.0/(2.0*(ng+1)))*(1/(rho*g))^(ng/(2.0*(ng+1)))*L^(1.0/2.0)
        
        prof = list(x=x,z=z,xc=xc,sigma=sigma,H0=H0,GG=GG)

        # Calculate H(x) 
        prof$H = H0 * (1.0-(abs(prof$xc)/L)^((ng+1.0)/ng))^(ng/(2.0*(ng+1.0)))

        prof$ux   = array(0.0,dim=c(nx,nz)) 
        prof$uz   = array(0.0,dim=c(nx,nz))
        prof$age  = rep(0.0,nz)
        prof$dHdx = array(0.0,dim=c(nx,nz))

        # Calculate velocities 
        for (i in 1:(nx-1)) {
            dHdx = (prof$H[i+1]-prof$H[i])/(prof$xc[i+1]-prof$xc[i])
            prof$dHdx[i] = dHdx

            H    = prof$H[i] 
            for (j in 1:nz) {

                prof$ux[i,j] = -(2.0*A)/(ng+1.0)*(rho*g)^ng * dHdx^(ng-1.0) *
                  dHdx * (H^(ng+1.0)-(H-prof$sigma[j]*H)^(ng+1.0))

                prof$uz[i,j] = prof$sigma[j]*(-GG+B+prof$ux[i,j]*dHdx) - B

                
            }
        } 

        # Calculate analytical age at the divide
        prof$age = calc_rh_age(prof$sigma,H0,GG)

        return(prof)
}

load_tracer_profile = function(filename,rh)
{
    # Load data 
    trc = my.read.nc(filename)

    # Get sigma axis 
    trc$sigma = trc$z / max(rh$H)

    # Convert time to ka 
    trc$time     = trc$time*1e-3
    trc$dep_time = trc$dep_time*1e-3
    trc$age      = trc$age*1e-3 

    # Calculate analytical age and error
    trc$rh_age    = calc_rh_age(sigma=trc$sigma,H0=rh$H0,GG=rh$GG)*1e-3
    trc$age_err   = (trc$age - trc$rh_age)
    trc$age_err_p = 100* (trc$age - trc$rh_age) / trc$rh_age
    trc$age_err_p[abs(trc$age)<1e-1] = NA 

    return(trc)
}

# Load data 
if (TRUE) {

    fldr = "output/RH2003_analytic"

    rh = calc_profile_RH2003()
    rh$age = rh$age*1e-3 

    rh0     = my.read.nc(file.path(fldr,"profile_RH2003.nc"))
    rh0$age = rh0$age*1e-3 

    trc_dp = load_tracer_profile(file.path(fldr,"RH2003_51_2_0.10_0.00_dp.nc"),rh=rh)
    trc_sp = load_tracer_profile(file.path(fldr,"RH2003_51_2_0.10_0.00_sp.nc"),rh=rh)
    
}

ptype = "png"
colax = "grey40"

# Plot comparison at dome
if (TRUE) {

    ylim  = c(-0.01,1.02)
    y.at  = seq(0,1,by=0.1)

    xlim  = c(-164,5)
    x.at  = seq(-160,0,by=40)
    x.lab = abs(x.at)

    xlim1 = c(-1e-4,1e-4)*1e3
    x.at1 = pretty(xlim1)
    xlim2 = c(1e-6,0.1)
    x.at2 = c(1e-5,1e-4,0.001,0.01,0.1)
    
    
    myfigure("plots","age_ice-divide",asp=1.9,pointsize=8,type=ptype)
    par(col=colax,col.axis=colax,col.lab=colax)

    par(plt=c(0.1,0.40,0.1,0.9),xaxs="i",yaxs="i")
    colnow = "magenta"
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    mtext(side=1,line=1.4,las=0,"Age (ka)")
    mtext(side=2,line=1.9,las=0,"Ice divide thickness (dimensionless)")
    axis(1,at=x.at,lab=x.lab)
    axis(2,at=y.at)
    abline(v=x.at,h=y.at,lwd=1,lty=2,col="lightgrey")

    lines(-trc_dp$rh_age,trc_dp$sigma,col="grey20",lwd=5)
    lines(-trc_dp$age,   trc_dp$sigma,col=colnow,lwd=1.5,lty=1)
    
    legend("topleft",bty="n",inset=0.01,col=c("grey20",colnow),lwd=c(5,1.5),lty=c(1,1),
        c("Analytical solution","Lagrangian tracer"),cex=0.8)
    box()  
    
    par(new=TRUE,plt=c(0.48,0.95,0.1,0.9),xaxs="i",yaxs="i")
    plot(xlim2,ylim,type="n",ann=FALSE,axes=FALSE,log="x")
    mtext(side=1,line=1.4,las=0,"Age err (%)")
    axis(1,at=x.at2)
    axis(2,at=y.at)
    abline(v=x.at2,h=y.at,lwd=1,lty=2,col="lightgrey")

    lines(abs(trc_dp$age_err_p),trc_dp$sigma,col=1,lwd=1.5)
    lines(abs(trc_sp$age_err_p),trc_sp$sigma,col=alpha(1,50),lwd=1.5)

    text(xlim2[2],mean(ylim)+diff(ylim)*0.01,pos=2,cex=0.8,col=1,          "Double prec.")
    text(xlim2[2],mean(ylim)-diff(ylim)*0.02,pos=2,cex=0.8,col=alpha(1,50),"Single prec.")
    
    box() 

    # par(new=TRUE,plt=c(0.50,0.70,0.1,0.9),xaxs="i",yaxs="i")
    # colnow = "slateblue"
    # plot(xlim1,ylim,type="n",ann=FALSE,axes=FALSE)
    # mtext(side=3,line=1.4,las=0,"Age err (a)",cex=0.8,col=colnow)
    # axis(3,at=x.at1,cex.axis=0.8,cex.lab=0.8,col=colnow,col.lab=colnow,col.axis=colnow)
    
    # lines(trc_dp$age_err*1e3,trc_dp$sigma,col=colnow,lwd=2)
    # # lines(trc_sp$age_err*1e3,trc_sp$sigma,col=alpha(1,50),lwd=2)

    graphics.off()

}