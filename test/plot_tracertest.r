

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
        if (length(jj) > 0) col_x[jj] = col[i]
    }

    # End cases 
    jj = which(x < breaks[1])
    if (length(jj) > 0) col_x[jj] = col[1]

    jj = which(x > breaks[length(breaks)])
    if (length(jj) > 0) col_x[jj] = col[length(col)]

    return(col_x)
}

if (TRUE) {
    
    # Load boundary conditions 
    topo = load_nc("../data/trace1_relax/Grisli15_3D.nc")
    topo$z_srf[topo$z_srf<=2] = 0 
    topo$U = sqrt(topo$Ux^2 + topo$Uy^2)

    topo$ratio = topo$Uz[,,nz] / topo$U[,,nz]
    topo$ratio[!is.finite(topo$ratio)] = 0 

    # Load experiment 
    trc = load_nc("../output/GRL-20KM_trc1.nc")
    nt = length(trc$time)
    trc$U = sqrt(trc$ux^2+trc$uy^2)

    # time strings
    trc$time_str = paste(trc$time)
    kk = trc$time<100
    trc$time_str[kk] = paste0("0",trc$time[kk])
    kk = trc$time<10
    trc$time_str[kk] = paste0("00",trc$time[kk])
    # kk = trc$time >= 10 & trc$time<100
    # trc$time_str[kk] = paste0("0",trc$time[kk])
       
}

## Options 
outfldr = "plots/test1"   # const/opt 
colax   = "grey20"
ptype   = "png"



# Plot points at initial elevations
if (TRUE) {

    xlim = c(-890,890)
    ylim = c(-1490,1490)

    breaks = c(-1,seq(0,3500,length.out=100))
    col    = c("grey90",colorRampPalette(c("grey30","grey85","grey95"))(99))

    breaks_U = c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500)
    col_U    = colorRampPalette(jet.colors)(length(breaks_U)-1)
    for (q in 1:length(col_U)) col_U[q] = alpha(col_U[q],max(50,100*(q/length(col_U))))
    breaks_U_labs = paste(breaks_U)
    breaks_U_labs[seq(2,length(breaks_U),by=2)] = ""

    for (k in 1:nt) {
        fnm = paste0("trc_zfix_time",trc$time_str[k])
        myfigure(outfldr,fnm,asp=0.8,pointsize=8,width=120,type=ptype)
        par(plt=c(0.05,0.95,0.05,0.95),xaxs="i",yaxs="i")

        plot(xlim,ylim,type="n",ann=FALSE,axes=FALSE)
        image(topo$xc,topo$yc,topo$z_srf,add=TRUE,breaks=breaks,col=col)
        contour(topo$xc,topo$yc,topo$z_srf,add=TRUE,levels=seq(500,3500,by=500),lwd=1,col="grey30")

        # Overlay points
        col_z = get_col(trc$U[,k],breaks=breaks_U,col=col_U)
        points(trc$x[,k],trc$y[,k],pch=21,bg=col_z,lwd=1)
        
        # contour(topo$xc,topo$yc,topo$z_srf,add=TRUE,levels=2000,lwd=3,col=1)

        par(plt=c(0.75,0.8,0.1,0.3),new=TRUE)
        mylegend(breaks=log(breaks_U),col=col_U,vertical=TRUE,at=log(breaks_U),labels=breaks_U_labs)
        title("Horizontal vel.",line=0.5,cex=0.6)
        mtext(side=4,line=2.3,las=1,"[m/a]",cex=1)

        graphics.off()

    }


}



