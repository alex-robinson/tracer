
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



# Load data 
if (TRUE) {

    rh     = my.read.nc("output/profile_RH2003.nc")
    # rh$x   = rh$x*1e-3 
    rh$age = rh$age*1e-3 

    trc = my.read.nc("output/profile_RH2003_trc1.nc")
    trc$time     = trc$time*1e-3
    trc$dep_time = trc$dep_time*1e-3

    trcmax = filt_tracers(trc)

    trcmax$rh_age  = approx(x=rh$sigma*rh$H[1],y=rh$age,xout=trcmax$z)$y
    trcmax$age_err = 100* (trcmax$dep_time - trcmax$rh_age) / trcmax$rh_age

}

ptype = "png"
colax = "grey40"

# Plot comparison at dome
if (TRUE) {

    dx = 100 
    it = which(trc$time>-160)

    
    xlim  = c(-205,5)
    ylim  = c(-0.1,3.9)
    xlim1 = range(abs(trcmax$age_err))
    xlim1 = c(1e-4,20)

    myfigure("plots","age_ice-divide",asp=1.6,pointsize=16,type=ptype)
    par(col=colax,col.axis=colax,col.lab=colax)

    par(mfrow=c(1,2))

    par(plt=c(0.12,0.95,0.07,0.95),xaxs="i",yaxs="i")
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    mtext(side=1,line=1.4,las=0,"Age (ka)")
    mtext(side=2,line=1.9,las=0,"Ice thickness (km)")
    axis(1,at=seq(-160,0,by=40))
    axis(2,at=seq(0,3.5,by=0.5))
    grid()

    lines(rh$age,rh$sigma*rh$H[1]*1e-3,lwd=4,col=1)
    points(trcmax$dep_time,trcmax$z*1e-3,col=2,cex=1.2)

    par(plt=c(0.12,0.95,0.07,0.95),xaxs="i",yaxs="i")
    plot(xlim1,ylim,type="n",ann=FALSE,axes=FALSE,log="x")
    mtext(side=1,line=1.4,las=0,"Age err (%)")
    mtext(side=2,line=1.9,las=0,"Ice thickness (km)")
    axis(1,at=c(1e-4,1e-3,1e-2,1e-1,1,10,20,100))
    axis(2,at=seq(0,3.5,by=0.5))
    grid()

    abline(v=0,lwd=4,col=1)
    points(abs(trcmax$age_err),trcmax$z*1e-3,col=2,cex=1.2)

    graphics.off()

}