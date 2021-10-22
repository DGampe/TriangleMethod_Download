ndvilst.plot <-
  function(NDVI,LST, t.min = c(), t.max = c(), classify = T, write=F, outpath = getwd(), outname = "ndvi-ts-plot",
           plot.unit = "C", n.highest = 2,
           plot.ylim = NULL, plot.xlim = c((-1.0),1.0), plot.main = "",
           plot.xlab = "NDVI [-]", plot.ylab = "Surface Temperature [°C]",
           plot.wetcol = "blue", plot.drycol = "black", plot.obscol="red", plot.ndvi0 = "forestgreen",
           textwet.x = -0.9, textwet.lab = "wet edge", textwet.col = "blue",
           textwet.font = 2, textwet.cex = 0.9, textobs.x = 0.5, textobs.lab = "observed dry edge",
           textobs.col = "red", textobs.font = 2, textobs.cex = 0.9, textdry.x = -0.6, textdry.lab = "true dry edge",
           textdry.col = "black", textdry.font = 2, textdry.cex = 0.9, textwater.x = -0.5, textwater.lab = "water",
           textwater.col = "black", textwater.font = 2, textwater.cex = 0.9, arrow.x0 = -0.5, arrow.y0 = 10, arrow.x1 = -0.23,
           arrow.angle = 30, arrow.code = 2, arrow.length = 0.15, arrow.lwd = 2, arrow.col = "black",
           pointmin.x = 0, pointmin.col = "black",	pointmin.type = "p", pointmin.lwd = 3, pointmin.pch = 16,
           textmin.x = 0.13, textmin.lab = "phi min", textmin.col = "black", textmin.font = 2, textmin.cex = 0.7,
           textndvi.x = 0.15, textndvi.lab = "NDVI = 0", textndvi.col = "forestgreen", textndvi.font = 2,textndvi.cex = 0.7,
           pointmax.col = "black", pointmax.type = "p", pointmax.lwd = 3, pointmax.pch = 16, textmax.x = 0.93, textmax.lab = "phi max",
           textmax.col = "black", textmax.font = 2, textmax.cex = 0.7, textwet.y = NULL, textobs.y = NULL, textdry.y = NULL,
           textndvi.y = NULL, textwater.y = NULL, arrow.y1 = NULL, pointmin.y = NULL, textmin.y = NULL, pointmax.x = NULL, pointmax.y = NULL,
           textmax.y = NULL, ...){

    if(class(NDVI)[[1]]!="RasterLayer"){
      NDVI = raster(NDVI)
    }
    if(class(LST)[[1]]!="RasterLayer"){
      LST = raster(LST)
    }
    if(is.null(plot.ylim)){plot.ylim=range(LST[],na.rm=T)+0.3*c(-1,1)*diff(range(LST[],na.rm=T))}

    if(classify==T){
      a = seq((-1),0.95,by=0.05)
      b = seq((-0.95),1,by=0.05)
      c = seq((-0.95),1,by=0.05)
      abc = cbind(a,b,c)
      class.NDVI = reclassify(NDVI,abc)
    }else{
      class.NDVI = NDVI
    }

    y1=sort(LST[], decreasing = T)[n.highest]   # Second highest temperature
    x1=mean(NDVI[which(LST[]==sort(LST[], decreasing = T)[n.highest])]) # NDVI that corresponds to 2nd highest temp
    y2= min(LST[which(class.NDVI[]==sort(class.NDVI[], decreasing = T)[n.highest])])# LST of 2nd highest NDVI
    x2=sort(NDVI[], decreasing = T)[n.highest]
    m=(y2-y1)/(x2-x1)
    t=y1-m*x1
    x=1
    y = y2

    if(!is.null(t.min)){ y = t.min}
    if(!is.null(t.max)){ t = t.max}

    m= (y-t)/(x2)

    text.tmin = paste("t.min:",round(y,2))
    text.tmax = paste("t.max:",round(t,2))

    stats.out=c()
    stats.out$t.min = round(y,2)
    stats.out$t.max = round(t,2)
    stats.out$y1 = y1
    stats.out$y2 = y2
    stats.out$x1 = x1
    stats.out$x2 = x2
    stats.out$m = m



    if(is.null(textwet.y)){textwet.y = y-1}
    if(is.null(textobs.y)){textobs.y = t-2}
    if(is.null(textdry.y)){textdry.y = t+5}
    if(is.null(textndvi.y)){textndvi.y = min(range(LST[],na.rm=T)+0.3*c(-1,1)*diff(range(LST[],na.rm=T))) + 3}
    if(is.null(textwater.y)){textwater.y = m+t-2}
    if(is.null(arrow.y1 )){arrow.y1 = (m+t+3)}
    if(is.null(pointmin.y)){pointmin.y = t}
    if(is.null(textmin.y)){textmin.y = t+2}
    if(is.null(pointmax.y)){pointmax.x = (y-t)/m}
    if(is.null(pointmax.y)){pointmax.y = y}
    if(is.null(textmax.y)){textmax.y = y-0.8} #m+t-2}
    if(is.null(plot.ylab)){
      if(min(range(LST[],na.rm=T)+0.3*c(-1,1)*diff(range(LST[],na.rm=T)))>200){
        plot.ylab = "Surface Temperature [K]"
      }else{
        plot.ylab = "Surface Temperature [°C]"
      }
    }

    if(write == T){
      setwd(outpath)
      png(paste(outname,".png",sep=""))
      smoothScatter(class.NDVI[], LST[]
                    , transformation = function(x) x^.19# the smaller the value the more colourful are also areas with low ponit density.
                    , ylim=plot.ylim, xlim=plot.xlim
                    , xlab="", ylab="", bty="n"
                    , main=plot.main,...)

      mtext(plot.xlab,1,padj=4, cex=0.8)
      mtext(plot.ylab,2,padj=-4, cex=0.8)
      abline(v=0,col=plot.ndvi0,lwd=3)
      abline(h= y,col=plot.wetcol,lwd=3)
      abline(a=t, b=m,col=plot.obscol,lwd=3)
      abline(a=t, b=(y-t/x) ,col=plot.drycol,lwd=2,lty="dashed")
      text(textwet.x,textwet.y,labels=textwet.lab,col=textwet.col,font=textwet.font, cex=textwet.cex)
      text(textobs.x,textobs.y,labels=textobs.lab,col=textobs.col,font=textobs.font, cex=textobs.cex)
      text(textndvi.x,textndvi.y,labels=textndvi.lab,col=textndvi.col,font=textndvi.font, cex=textndvi.cex)
      text(textdry.x,textdry.y,labels=textdry.lab,col=textdry.col,font=textdry.font, cex=textdry.cex)
      text(textwater.x,textwater.y,labels=textwater.lab,col=textwater.col,font=textwater.font, cex=textwater.cex)
      arrows(arrow.x0,arrow.y0,arrow.x1,arrow.y1,code=arrow.code,length=arrow.length,lwd=arrow.lwd,col = arrow.col)
      points(pointmin.x,pointmin.y,type=pointmin.type,col=pointmin.col,lwd=pointmin.lwd, pch=pointmin.pch)
      text(textmin.x ,textmin.y,labels=textmin.lab,col=textmin.col,font=textmin.font, cex=textmin.cex)
      points(pointmax.x,pointmax.y,type=pointmax.type,col=pointmax.col,lwd=pointmax.lwd, pch=pointmax.pch)
      text(textmax.x ,textmax.y,labels=textmax.lab,col=textmax.col,font=textmax.font, cex=textmax.cex)

      dev.off()
    }else{
      smoothScatter(class.NDVI[], LST[]
                    , transformation = function(x) x^.19# the smaller the value the more colourful are also areas with low ponit density.
                    , ylim=plot.ylim, xlim=plot.xlim
                    , xlab="", ylab="", bty="n"
                    , main=plot.main, ...)

      mtext(plot.xlab,1,padj=4, cex=0.8)
      mtext(plot.ylab,2,padj=-4, cex=0.8)
      abline(v=0,col=plot.ndvi0,lwd=3)
      abline(h= y,col=plot.wetcol,lwd=3)
      abline(a=t, b=m,col=plot.obscol,lwd=3)
      abline(a=t, b=(y-t/x) ,col=plot.drycol,lwd=2,lty="dashed")
      text(textwet.x,textwet.y,labels=textwet.lab,col=textwet.col,font=textwet.font, cex=textwet.cex)
      text(textobs.x,textobs.y,labels=textobs.lab,col=textobs.col,font=textobs.font, cex=textobs.cex)
      text(textndvi.x,textndvi.y,labels=textndvi.lab,col=textndvi.col,font=textndvi.font, cex=textndvi.cex)
      text(textdry.x,textdry.y,labels=textdry.lab,col=textdry.col,font=textdry.font, cex=textdry.cex)
      text(textwater.x,textwater.y,labels=textwater.lab,col=textwater.col,font=textwater.font, cex=textwater.cex)
      arrows(arrow.x0,arrow.y0,arrow.x1,arrow.y1,code=arrow.code,length=arrow.length,lwd=arrow.lwd,col = arrow.col)
      points(pointmin.x,pointmin.y,type=pointmin.type,col=pointmin.col,lwd=pointmin.lwd, pch=pointmin.pch)
      text(textmin.x ,textmin.y,labels=textmin.lab,col=textmin.col,font=textmin.font, cex=textmin.cex)
      points(pointmax.x,pointmax.y,type=pointmax.type,col=pointmax.col,lwd=pointmax.lwd, pch=pointmax.pch)
      text(textmax.x ,textmax.y,labels=textmax.lab,col=textmax.col,font=textmax.font, cex=textmax.cex)
    }
    return(stats.out)
    invisible(y)
  }
