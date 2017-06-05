
library(myr)
source("functions_tracers.r")

# Load data 
if (FALSE) {

    fldr = "output/RH2003"

    rh = calc_profile_RH2003()

    rh0     = my.read.nc(file.path(fldr,"profile_RH2003.nc"))
    rh0$age = rh0$age*1e-3 

    trc_ref = load_tracer_profile(file.path(fldr,"RH2003_501_201_0.10_0.00_dp_10.0.nc"))
    trc_1   = load_tracer_profile(file.path(fldr,"RH2003_51_101_0.10_0.00_sp_10.0.nc"))

    # trc_ref = load_tracer_profile(file.path(fldr,"RH2003_501_201_0.10_0.00_sp_10.0.nc"))
    # # trc_1   = load_tracer_profile(file.path(fldr,"RH2003_51_51_0.10_0.00_sp_10.0.nc"))
    
    # Generate gridded results for comparisons
    grd = list(x=seq(0,1000,length.out=101),y=seq(0,1,length.out=101))

    kk = which(!is.na(trc_ref$age))
    qp0 = as.image(Z=trc_ref$age[kk],x=cbind(trc_ref$x[kk],trc_ref$z[kk]/trc_ref$H[kk]),grid=grd)

    kk = which(!is.na(trc_1$age))
    qp1 = as.image(Z=trc_1$age[kk],x=cbind(trc_1$x[kk],trc_1$z[kk]/trc_1$H[kk]),grid=grd)

    grd = list(x=seq(0,1000,length.out=101),y=seq(0,rh$H0,length.out=101))
    rh_stream = as.image(Z=as.vector(rh$umag),x=cbind(as.vector(rh$xx),as.vector(rh$zz)),grid=grd)

}

ptype = "png"
colax = "grey40"

# PLOT age estimate with error comparison with high resolution reference sim
if (FALSE) {

    myfigure("plots","age-err_resolution",asp=2.2,pointsize=12,type=ptype)
    par(col=colax,col.axis=colax,col.lab=colax)

    par(mfrow=c(1,3),xaxs="i",yaxs="i")
    
    xlim = range(grd$x)
    ylim = range(grd$y)

    # Age estimate

    z = qp0$z 

    breaks = c(0,2,5,10,20,30,40,50,80,100,120,150)
    col    = colorRampPalette(jet.colors)(length(breaks)-1)

    z[z>max(breaks)] = max(breaks)
    z[z<min(breaks)] = min(breaks)

    par(plt=c(0.1,0.95,0.20,0.92)) 
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(grd$x,grd$y,z,add=TRUE,breaks=breaks,col=col)
    title(main="Age (ka)",line=0.5)
    box()

    par(new=TRUE,plt=c(0.2,0.8,0.1,0.14)) 
    mylegend(breaks,col,vertical=FALSE)
    
    # Age error

    z = (qp1$z-qp0$z)
    
    breaks = pretty(c(-2,2),15)
    col    = colorRampPalette(jet.colors)(length(breaks)-1)

    z[z>max(breaks)] = max(breaks)
    z[z<min(breaks)] = min(breaks)
    
    par(plt=c(0.1,0.95,0.20,0.92)) 
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(grd$x,grd$y,z,add=TRUE,breaks=breaks,col=col)
    title(main="Age error (ka)",line=0.5)
    box()

    par(new=TRUE,plt=c(0.2,0.8,0.1,0.14)) 
    mylegend(breaks,col,vertical=FALSE)

    # Age error (percent)

    z = 100*abs(qp1$z-qp0$z)/qp0$z
    
    breaks = pretty(c(0,10),15)
    col    = colorRampPalette(jet.colors,bias=2)(length(breaks)-1)

    z[z>max(breaks)] = max(breaks)
    z[z<min(breaks)] = min(breaks)
    
    par(plt=c(0.1,0.95,0.20,0.92)) 
    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(grd$x,grd$y,z,add=TRUE,breaks=breaks,col=col)
    title(main="Age error (%)",line=0.5)
    box()

    par(new=TRUE,plt=c(0.2,0.8,0.1,0.14)) 
    mylegend(breaks,col,vertical=FALSE)
    
    graphics.off()

}

# Plot streamlines to compare with RH2003
if (TRUE) {

    zlim = c(0,30)
    tmp = rh_stream 
    tmp$z[tmp$z>zlim[2]] = zlim[2]
    tmp$z[tmp$z<zlim[1]] = zlim[1]
    image.plot(tmp,zlim=zlim)

}