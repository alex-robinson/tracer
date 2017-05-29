
library(myr)




# Load data 
if (TRUE) {

    rh     = my.read.nc("output/profile_RH2003.nc")
    rh$x   = rh$x*1e-3 
    rh$age = rh$age*1e-3 

    trc = my.read.nc("output/profile_RH2003_trc1.nc")
    trc$time     = trc$time*1e-3
    trc$dep_time = trc$dep_time*1e-3

}

ptype = "png"
colax = "grey40"

# Plot comparison at dome
if (TRUE) {

    dx = 100 
    it = which(trc$time>-160)

    xlim = c(-205,5)
    ylim = c(-0.1,3.9)

    myfigure("plots","age_ice-divide",asp=0.8,pointsize=16,type=ptype)
    par(col=colax,col.axis=colax,col.lab=colax)

    par(plt=c(0.12,0.95,0.07,0.95),xaxs="i",yaxs="i")
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    mtext(side=1,line=1.4,las=0,"Age (ka)")
    mtext(side=2,line=1.9,las=0,"Ice thickness (km)")
    axis(1,at=seq(-160,0,by=40))
    axis(2,at=seq(0,3.5,by=0.5))
    grid()

    lines(rh$age,rh$sigma*rh$H[1]*1e-3,lwd=4,col=1)
    points(trc$dep_time,trc$z*1e-3,col=2,cex=1.2)

    graphics.off()

}