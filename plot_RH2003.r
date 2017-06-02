
library(myr)
source("functions_tracers.r")

# Load data 
if (FALSE) {

    fldr = "output/RH2003"

    rh = calc_profile_RH2003()

    rh0     = my.read.nc(file.path(fldr,"profile_RH2003.nc"))
    rh0$age = rh0$age*1e-3 

    trc_51_101_0.10_0.00_sp_10.0 = load_tracer_profile(file.path(fldr,"RH2003_51_101_0.10_0.00_sp_10.0.nc"))

}

if (TRUE) {

    kk = which(!is.na(trc_51_101_0.10_0.00_sp_10.0$age))

    quilt.plot(x=trc_51_101_0.10_0.00_sp_10.0$x[kk],
               y=trc_51_101_0.10_0.00_sp_10.0$z[kk]/trc_51_101_0.10_0.00_sp_10.0$H[kk],
               z=trc_51_101_0.10_0.00_sp_10.0$age[kk],nx=101,ny=101,zlim=c(0,100))
    # points(x=trc_51_101_0.10_0.00_sp_10.0$x[kk],
    #            y=trc_51_101_0.10_0.00_sp_10.0$z[kk]*1e-3,cex=0.2)
    

    # quilt.plot(x=as.vector(rh$xx),y=as.vector(rh$zz),z=as.vector(rh$umag))

}

ptype = "png"
colax = "grey40"

# Plot comparison at dome
if (FALSE) {

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