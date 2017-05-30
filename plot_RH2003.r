
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

    trcmax$rh_age  = approx(x=rh$sigma*max(rh$H),y=rh$age,xout=trcmax$z)$y
    trcmax$age_err = 100* (trcmax$dep_time - trcmax$rh_age) / trcmax$rh_age

}

ptype = "png"
colax = "grey40"

# Plot comparison at dome
if (TRUE) {

    xlim  = c(-165,5)
    ylim  = c(-0.01,1.02)
    xlim1 = range(abs(trcmax$age_err))
    xlim1 = c(1e-2,10)

    x.at  = seq(-160,0,by=40)
    y.at  = seq(0,1,by=0.1)
    x.at1 = c(0.001,0.01,0.1,1,10)

    myfigure("plots","age_ice-divide",asp=1.4,pointsize=12,type=ptype)
    par(col=colax,col.axis=colax,col.lab=colax)

    par(plt=c(0.1,0.55,0.1,0.95),xaxs="i",yaxs="i")
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    mtext(side=1,line=1.4,las=0,"Age (ka)")
    mtext(side=2,line=1.9,las=0,"Ice divide thickness (dimensionless)")
    axis(1,at=x.at)
    axis(2,at=y.at)
    abline(v=x.at,h=y.at,lwd=1,lty=2,col="lightgrey")

    lines(rh$age,rh$sigma,lwd=4,col=1)
    points(trcmax$dep_time,trcmax$z/max(rh$H),col=2,cex=0.8)

    box()

    par(new=TRUE,plt=c(0.65,0.95,0.1,0.95),xaxs="i",yaxs="i")
    plot(xlim1,ylim,type="n",ann=FALSE,axes=FALSE,log="x")
    mtext(side=1,line=1.4,las=0,"Age err (%)")
    mtext(side=2,line=1.9,las=0,"Ice divide thickness (dimensionless)")
    axis(1,at=x.at1)
    axis(2,at=y.at)
    abline(v=x.at1,h=y.at,lwd=1,lty=2,col="lightgrey")

    abline(v=0,lwd=4,col=1)
    points(abs(trcmax$age_err),trcmax$z/max(rh$H),col=2,cex=0.8)

    box() 

    graphics.off()

}