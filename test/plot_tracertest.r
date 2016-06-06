

library(myr)
library(RNetCDF)

load_nc = function(filename)
{
    nc = open.nc(filename)
    dat = read.nc(nc)
    close.nc(nc)

    return(dat)
}

get_col = function (x,breaks,col) 
{   
    col_x = rep("black",length(x))

    for (i in 1:length(col)) {
        jj = which(x > breaks[i] & x <= breaks[i+1])
        cat(breaks[i],breaks[i+1],length(jj),"\n")
        if (length(jj) > 0) col_x[jj] = col[i]
    }

    # End cases 
    jj = which(x < breaks[1])
    if (length(jj) > 0) col_x[jj] = col[1]

    jj = which(x > breaks[length(breaks)])
    if (length(jj) > 0) col_x[jj] = col[length(col)]

    return(col_x)
}

if (FALSE) {
    
    # Load boundary conditions 
    topo = load_nc("../data/GRL-20KM_TOPO-B13_gl0.05.nc")


    # Load experiment 
    trc = load_nc("../output/GRL-20KM_trc1.nc")
    nt = length(trc$time)

}

## Options 
outfldr = "plots"   # const/opt 
colax   = "grey20"
ptype   = "png"



# Plot points at initial elevations
if (TRUE) {

    xlim = c(-890,890)
    ylim = c(-1490,1490)

    breaks = c(-1,seq(0,3000,length.out=50))
    col    = c("grey90",terrain.colors(49))


    myfigure(outfldr,"trc_init_z",asp=0.8,pointsize=12,type=ptype)
    par(plt=c(0.05,0.95,0.05,0.95),xaxs="i",yaxs="i")

    plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
    image(topo$xc,topo$yc,topo$zs,add=TRUE,breaks=breaks,col=col)

    # Overlay points
    col_z = get_col(trc$z[,nt],breaks=breaks,col=col)
    points(trc$x[,nt],trc$y[,nt],pch=21,bg=col_z,lwd=1)




    graphics.off()

}



