# R-CMplot
CMplot <- function(Pmap,col=c('red','black','green','blue','orange'),pch=19,band=1,cir.band=1,H=1,out="b",cex=c(0.5,1),r=1,outward=TRUE,line = TRUE, amplify = TRUE, cir.labels=TRUE,amplify.col=c("red","green")){
  if(out=="b") out=c("c","m")
  Pmap=Pmap[,-1]
  taxa=colnames(Pmap)[-c(1,2)]
  cir.band=cir.band/5
  if(length(cex)!=2) cex[2]=cex[1]
  R=dim(Pmap)[2]-2
  index=c(1:100)
  PmapN=Pmap[Pmap[,1] %in% index,]
  PmapXY=Pmap[!(Pmap[,1] %in% index)&Pmap[,1]!=0,]
  PmapN=matrix(as.numeric(as.matrix(PmapN)),nrow(PmapN))
  PmapXY=as.matrix(PmapXY)
  orderP=PmapN[order(PmapN[,1],PmapN[,2]),]
  orderPXY.=PmapXY[order(PmapXY[,1],PmapXY[,2]),]
  chr=c(unique(orderP[,1]),"S")
  pvalueT=as.matrix(orderP[,-c(1:2)])
  pvalueXY.=-log10(matrix(as.numeric(orderPXY.[,-c(1:2)]),nrow(orderPXY.)))
  #palette(heat.colors(1024)) #(heatmap)
  #T=floor(1024/max(pvalue))
  #plot(pvalue,pch=19,cex=0.6,col=(1024-floor(pvalue*T)))
  if(!missing(col)){
    col=col
  }else{
    col=c('red','black','green','blue','orange')
  }
  Num=as.numeric(table(PmapN[,1]))
  Nchr=length(Num)
  N=ceiling(Nchr/length(col))
  if(!missing(band)){
    band=floor(band)
    band=band*floor(dim(pvalueT)[1]/100)
  }else{
    band=floor(dim(pvalueT)[1]/100)
  }
  ticks=NULL
  NewP=NULL
  for(j in 1:R){
    pvalue=pvalueT[,j]
    for(i in 0:(Nchr-1)){
      if (i==0){
        pvalue=append(pvalue,rep(0,band),after=0)
        ticks[i+1]=band+floor((Num[i+1])/2)
      }else{
        pvalue=append(pvalue,rep(0,band),after=sum(Num[1:i])+i*band)
        ticks[i+1]=band*(i+1)+sum(Num[1:i])+floor((Num[i+1])/2)
      }
    }
    NewP[[j]]=pvalue
  }
  pvalueT=do.call(cbind,NewP)
  logpvalueT=-log10(pvalueT)
  Num=Num+band
  ticks=c(ticks,dim(pvalueT)[1]+band+dim(pvalueXY.)[1]/2)
  add=c(Num,rep(0,N*length(col)-Nchr))
  XYlim=(dim(pvalueT)[1]+band+1):((dim(pvalueT)[1]+band)+dim(pvalueXY.)[1])	
  TotalN1=dim(pvalueT)[1]
  TotalN2=dim(pvalueXY.)[1]
  TotalN=TotalN1+TotalN2+band
   if("c" %in% out){
  print("Starting Manhattan-Circle plot!",quote=F)
  jpeg("Manhattan-Plot.circle.jpg", width = 2450,height=2450,res=300,quality = 100)
  par(pty="s")
  RR=r+H*R+cir.band*R
  plot(NULL,xlim=c(-RR,RR),ylim=c(-RR,RR),axes=F,xlab="",ylab="")
  for(i in 1:R){
    pvalue=pvalueT[,i]
    logpvalue=logpvalueT[,i]
    pvalueXYi.=pvalueXY.[,i]
    Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
    Cpvalue=H*logpvalue/Max
    CpvalueXY=H*pvalueXYi./Max
	if(outward==TRUE){
    X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(1:TotalN1)/TotalN)
    Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(1:TotalN1)/TotalN)
    XX=(CpvalueXY+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*((TotalN1+band+1):TotalN)/TotalN)
    YY=(CpvalueXY+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*((TotalN1+band+1):TotalN)/TotalN)
    #plot(chrX,chrY,type="l",col="black",lwd=4)
    points(X,Y,pch=19,cex=cex[1],col=rep(rep(col,N),add))
    points(XX,YY,pch=19,cex=cex[1],col="gray")
    significantline1=H*(-log10(0.01/max(dim(Pmap))))/Max
    significantline2=H*(-log10(1/max(dim(Pmap))))/Max
	if(line==TRUE){
    s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
    s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
    s2X=(significantline2+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
    s2Y=(significantline2+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
    if(significantline1<H) lines(s1X,s1Y,type="l",col="red",lwd=1)
    if(significantline2<H) lines(s2X,s2Y,type="l",col="green",lwd=1,lty=2)
	}
	if(amplify == TRUE){
    HX1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
    HY1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
    HX2=(Cpvalue[which(Cpvalue<significantline1&Cpvalue>=significantline2)]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(which(Cpvalue<significantline1&Cpvalue>=significantline2))/TotalN)
    HY2=(Cpvalue[which(Cpvalue<significantline1&Cpvalue>=significantline2)]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(which(Cpvalue<significantline1&Cpvalue>=significantline2))/TotalN)
    points(HX1,HY1,pch=19,cex=cex[1]*1.5,col=amplify.col[1])
    points(HX2,HY2,pch=19,cex=cex[1]*1.5,col=amplify.col[2])
	}
	ticks1=(0.9*r)*sin(2*pi*ticks/TotalN)
    ticks2=(0.9*r)*cos(2*pi*ticks/TotalN)
	if(cir.labels==TRUE){
    for(i in 1:length(ticks)){
      angle=360*(1-ticks[i]/TotalN)
      text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=0.5)
    }
	}
	}
	if(outward==FALSE){
	X=(-Cpvalue+r+H*i+cir.band*(i-1))*sin(2*pi*(1:TotalN1)/TotalN)
    Y=(-Cpvalue+r+H*i+cir.band*(i-1))*cos(2*pi*(1:TotalN1)/TotalN)
    XX=(-CpvalueXY+r+H*i+cir.band*(i-1))*sin(2*pi*((TotalN1+band+1):TotalN)/TotalN)
    YY=(-CpvalueXY+r+H*i+cir.band*(i-1))*cos(2*pi*((TotalN1+band+1):TotalN)/TotalN)
	points(X,Y,pch=19,cex=cex[1],col=rep(rep(col,N),add))
    points(XX,YY,pch=19,cex=cex[1],col="gray")
    significantline1=H*(-log10(0.01/max(dim(Pmap))))/Max
    significantline2=H*(-log10(1/max(dim(Pmap))))/Max
	if(line==TRUE){
    s1X=(-significantline1+r+H*i+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
    s1Y=(-significantline1+r+H*i+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
    s2X=(-significantline2+r+H*i+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
    s2Y=(-significantline2+r+H*i+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
    if(significantline1<H) lines(s1X,s1Y,type="l",col="red",lwd=1)
    if(significantline2<H) lines(s2X,s2Y,type="l",col="green",lwd=1,lty=2)
	}
	if(amplify == TRUE){
    HX1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
    HY1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
    HX2=(-Cpvalue[which(Cpvalue<significantline1&Cpvalue>=significantline2)]+r+H*i+cir.band*(i-1))*sin(2*pi*(which(Cpvalue<significantline1&Cpvalue>=significantline2))/TotalN)
    HY2=(-Cpvalue[which(Cpvalue<significantline1&Cpvalue>=significantline2)]+r+H*i+cir.band*(i-1))*cos(2*pi*(which(Cpvalue<significantline1&Cpvalue>=significantline2))/TotalN)
    points(HX1,HY1,pch=19,cex=cex[1]*1.5,col=amplify.col[1])
    points(HX2,HY2,pch=19,cex=cex[1]*1.5,col=amplify.col[2])
	}
    ticks1=1.07*(RR-cir.band)*sin(2*pi*ticks/TotalN)
    ticks2=1.07*(RR-cir.band)*cos(2*pi*ticks/TotalN)
	if(cir.labels==TRUE){
    for(i in 1:length(ticks)){
      angle=360*(1-ticks[i]/TotalN)
      text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=1)
    }
	}	
	}
  }
  dev.off()
  print("Manhattan-Circle has been finished!",quote=F)
}  
if("m" %in% out){
  print("Starting Manhattan-Genomewise plot!",quote=F)
  for(i in 1:R){
    jpeg(paste("Manhattan-Plot.",taxa[i],".jpg",sep=""), width = 4200,height=1500,res=300,quality = 100)
	par(mar = c(3,6,5,3),xaxs="i")
    pvalue=pvalueT[,i]
    logpvalue=logpvalueT[,i]
    pvalueXYi.=pvalueXY.[,i]
    Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
    plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=c(0,Max+1),xlab="",ylab=expression(-log[10](italic(p))),
         cex.axis=1.5,cex.lab=2,font=2,axes=FALSE)
    points(pvalueXYi.~XYlim,pch=pch,cex=cex[2],col="gray")
    axis(1, at=ticks,cex.axis=1,font=2,labels=chr)
    axis(2,at=c(0:(Max+1)),cex.axis=1,font=2,labels=0:(Max+1))
	if(line==TRUE){
    abline(h=-log10(0.01/max(dim(Pmap))))
    abline(h=-log10(1/max(dim(Pmap))),lty=2)
	}
	if(amplify == TRUE){
    sgline1=-log10(0.01/max(dim(Pmap)))
    sgline2=-log10(1/max(dim(Pmap)))
    HY1=logpvalue[which(logpvalue>=sgline1)]
    HX1=which(logpvalue>=sgline1)
    HY2=logpvalue[which(logpvalue<sgline1&logpvalue>=sgline2)]
    HX2=which(logpvalue<sgline1&logpvalue>=sgline2)
    points(HX1,HY1,pch=pch,cex=cex[2]*1.5,col=amplify.col[1])
    points(HX2,HY2,pch=pch,cex=cex[2]*1.5,col=amplify.col[2])
	}
	dev.off()
	print(paste("Manhattan-",taxa[i],".plot has been finished!",sep=""),quote=F)
  }
  }
  print(paste("The plots have been stored in ","[",getwd(),"]",sep=""),quote=F)
}

