# R-CMplot
CMplot <- function(Pmap,
col=c('red','black','green','blue','orange'),
pch=19,
band=1,
cir.band=1,
H=1,
ylim=NULL,
cex.axis=1,
output="b",
multracks=FALSE,
cex=c(0.5,1),
r=1,
xlab="Chromosome",
ylab=expression(-log[10](italic(p))),
outward=TRUE,
threshold = 0.01, 
threshold.col="red",
amplify= TRUE,
signal.cex = 1.5,
signal.pch = 19,
signal.col="red",
cir.chr=TRUE,
chr.band=1,
chr.col=NULL,
cir.labels=TRUE,
plot0=FALSE,
fill.output=TRUE,
fill="jpg",
dpi=300
)
{
  if(output=="b") output=c("c","m")
  plotXY=TRUE
  Pmap=Pmap[order(Pmap[,1]),]
  Pmap=Pmap[,-1]
  taxa=colnames(Pmap)[-c(1,2)]
  chr.band=chr.band/5
  cir.band=cir.band/5
  if(length(cex)!=2) cex[2]=cex[1]
  if(length(cex)>2) stop("Less than two cexs are allowed!")
  R=dim(Pmap)[2]-2
  if(plot0==FALSE) index=c(1:100)
  if(plot0==TRUE) index=c(0:100)
  PmapN=Pmap[Pmap[,1] %in% index,]
  PmapXY=Pmap[!(Pmap[,1] %in% index)&Pmap[,1]!=0,]
  chrXY=unique(PmapXY[,1])
  if(dim(PmapXY)[1]==0) plotXY=FALSE
  PmapN=matrix(as.numeric(as.matrix(PmapN)),nrow(PmapN))
  if(plotXY==TRUE) PmapXY=as.matrix(PmapXY)
  orderP=PmapN[order(PmapN[,1],PmapN[,2]),]
  if(plotXY==TRUE) orderPXY.=PmapXY[order(PmapXY[,1],PmapXY[,2]),]
  if(plotXY==TRUE){
  chr=c(unique(orderP[,1]),"S")
  }else{
   chr=unique(orderP[,1])
  }
  pvalueT=as.matrix(orderP[,-c(1:2)])
  if(plotXY==TRUE) pvalueXY.=-log10(matrix(as.numeric(orderPXY.[,-c(1:2)]),nrow(orderPXY.)))
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
  if(plotXY==TRUE){
  ticks=c(ticks,dim(pvalueT)[1]+band+dim(pvalueXY.)[1]/2)
  add=c(Num,rep(0,N*length(col)-Nchr))
  XYlim=(dim(pvalueT)[1]+band+1):((dim(pvalueT)[1]+band)+dim(pvalueXY.)[1])	
  TotalN1=dim(pvalueT)[1]
  TotalN2=dim(pvalueXY.)[1]
  TotalN=TotalN1+TotalN2+band
  }else{
  ticks=ticks
  add=c(Num,rep(0,N*length(col)-Nchr))
  TotalN1=dim(pvalueT)[1]
  TotalN=TotalN1
  }
   if("c" %in% output){
  print("Starting Circular-Manhattan plot!",quote=F)
if(fill.output==TRUE){
  if(fill=="jpg")	jpeg("Circular-Manhattan.jpg", width = 8*dpi,height=8*dpi,res=dpi,quality = 100)
  if(fill=="pdf")	pdf("Circular-Manhattan.pdf", width = 10,height=10)
  if(fill=="tiff")	tiff("Circular-Manhattan.tiff", width = 8*dpi,height=8*dpi,res=dpi)
  }
  par(pty="s",xpd=TRUE)
  RR=r+H*R+cir.band*R
  plot(NULL,xlim=c(1.05*(-RR-4*chr.band),1.05*(RR+4*chr.band)),ylim=c(1.05*(-RR-4*chr.band),1.05*(RR+4*chr.band)),axes=F,xlab="",ylab="")
  for(i in 1:R){
    print(paste("Circular-Plotting ",taxa[i],"...",sep=""),quote=F)
    pvalue=pvalueT[,i]
    logpvalue=logpvalueT[,i]
    if(plotXY==TRUE){
	pvalueXYi.=pvalueXY.[,i]
    Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
    Cpvalue=H*logpvalue/Max
    CpvalueXY=H*pvalueXYi./Max
	}else{
    Max=ceiling(-log10(min(pvalue[pvalue!=0])))
    Cpvalue=H*logpvalue/Max
	}
	if(outward==TRUE){
	if(cir.chr==TRUE){
	XLine=(RR+chr.band)*sin(2*pi*(1:TotalN)/TotalN)
	YLine=(RR+chr.band)*cos(2*pi*(1:TotalN)/TotalN)
	lines(XLine,YLine,lwd=1.5)
	a=0
	if(plotXY==FALSE){
	for(k in 1:length(chr)){
	if(k==1){
	X1chr=(RR)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y1chr=(RR)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	X2chr=(RR+chr.band)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y2chr=(RR+chr.band)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}
	}else{
	a=a+sum(Pmap[,1]==chr[k-1])
	X1chr=(RR)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y1chr=(RR)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	X2chr=(RR+chr.band)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y2chr=(RR+chr.band)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}		
	}
	}
	}
	if(plotXY==TRUE){
	for(k in 1:(length(chr)-1)){
	if(k==1){
	X1chr=(RR)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y1chr=(RR)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	X2chr=(RR+chr.band)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y2chr=(RR+chr.band)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}	
	}else{
	a=a+sum(Pmap[,1]==chr[k-1])
	X1chr=(RR)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y1chr=(RR)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	X2chr=(RR+chr.band)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y2chr=(RR+chr.band)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}	
	}
	}
	k=length(chr)
	a=a+sum(Pmap[,1]==chr[k-1])
	X1chr=(RR)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	Y1chr=(RR)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	X2chr=(RR+chr.band)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	Y2chr=(RR+chr.band)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="gray",border="gray")	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}
	}
	}
    X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(1:TotalN1)/TotalN)
    Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(1:TotalN1)/TotalN)
    if(plotXY==TRUE){
    XX=(CpvalueXY+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*((TotalN1+band+1):TotalN)/TotalN)
    YY=(CpvalueXY+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*((TotalN1+band+1):TotalN)/TotalN)
	points(XX,YY,pch=19,cex=cex[1],col="gray")
	points(X,Y,pch=19,cex=cex[1],col=rep(rep(col,N),add))
	}else{
    #plot(chrX,chrY,type="l",col="black",lwd=4)
    points(X,Y,pch=19,cex=cex[1],col=rep(rep(col,N),add))
	}
	if(!is.null(threshold)){
	if(threshold!=0){
	significantline1=H*(-log10(threshold/max(dim(Pmap))))/Max
    s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
    s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
    if(significantline1<H){
	lines(s1X,s1Y,type="l",col=threshold.col,lwd=1)
	}else{
	warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold,"!",sep=""))
	}
	if(amplify == TRUE){
    HX1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
    HY1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
	points(HX1,HY1,pch=19,cex=cex[1],col="white")
	points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=signal.col)
	}
	}
	}
	if(cir.chr==TRUE){
	ticks1=1.07*(RR+chr.band)*sin(2*pi*ticks/TotalN)
    ticks2=1.07*(RR+chr.band)*cos(2*pi*ticks/TotalN)
	if(cir.labels==TRUE){
    for(i in 1:length(ticks)){
      angle=360*(1-ticks[i]/TotalN)
      text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=1)
    }
	}
	}else{
	ticks1=(0.9*r)*sin(2*pi*ticks/TotalN)
    ticks2=(0.9*r)*cos(2*pi*ticks/TotalN)
	if(cir.labels==TRUE){
    for(i in 1:length(ticks)){
      angle=360*(1-ticks[i]/TotalN)
      text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=0.5)
    }
	}
	}
	}
	if(outward==FALSE){
	if(cir.chr==TRUE){
	XLine=(2*cir.band+RR+chr.band)*sin(2*pi*(1:TotalN)/TotalN)
	YLine=(2*cir.band+RR+chr.band)*cos(2*pi*(1:TotalN)/TotalN)
	lines(XLine,YLine,lwd=1.5)
	a=0
	if(plotXY==FALSE){
	for(k in 1:length(chr)){
	if(k==1){
	X1chr=(2*cir.band+RR)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y1chr=(2*cir.band+RR)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	X2chr=(2*cir.band+RR+chr.band)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y2chr=(2*cir.band+RR+chr.band)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}
	}else{
	a=a+sum(Pmap[,1]==chr[k-1])
	X1chr=(2*cir.band+RR)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y1chr=(2*cir.band+RR)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	X2chr=(2*cir.band+RR+chr.band)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y2chr=(2*cir.band+RR+chr.band)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}	
	}
	}
	}
	if(plotXY==TRUE){
	for(k in 1:(length(chr)-1)){
	if(k==1){
	X1chr=(2*cir.band+RR)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y1chr=(2*cir.band+RR)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	X2chr=(2*cir.band+RR+chr.band)*sin(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	Y2chr=(2*cir.band+RR+chr.band)*cos(2*pi*((band+1):(band+sum(Pmap[,1]==chr[1])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}	
	}else{
	a=a+sum(Pmap[,1]==chr[k-1])
	X1chr=(2*cir.band+RR)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y1chr=(2*cir.band+RR)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	X2chr=(2*cir.band+RR+chr.band)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	Y2chr=(2*cir.band+RR+chr.band)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(col,ceiling(length(chr)/length(col)))[k],border=rep(col,ceiling(length(chr)/length(col)))[k])	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}	
	}
	}
	k=length(chr)
	a=a+sum(Pmap[,1]==chr[k-1])
	X1chr=(2*cir.band+RR)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	Y1chr=(2*cir.band+RR)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	X2chr=(2*cir.band+RR+chr.band)*sin(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	Y2chr=(2*cir.band+RR+chr.band)*cos(2*pi*((k*band+a+1):(k*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
	if(is.null(chr.col)){
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="gray",border="gray")	
	}else{
	polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.col,border=chr.col)
	}
	}
	}
	X=(-Cpvalue+r+H*i+cir.band*(i-1))*sin(2*pi*(1:TotalN1)/TotalN)
    Y=(-Cpvalue+r+H*i+cir.band*(i-1))*cos(2*pi*(1:TotalN1)/TotalN)
	if(plotXY==TRUE){
    XX=(-CpvalueXY+r+H*i+cir.band*(i-1))*sin(2*pi*((TotalN1+band+1):TotalN)/TotalN)
    YY=(-CpvalueXY+r+H*i+cir.band*(i-1))*cos(2*pi*((TotalN1+band+1):TotalN)/TotalN)
	points(X,Y,pch=19,cex=cex[1],col=rep(rep(col,N),add))
    points(XX,YY,pch=19,cex=cex[1],col="gray")
	}else{
	points(X,Y,pch=19,cex=cex[1],col=rep(rep(col,N),add))
	}
	if(!is.null(threshold)){
	if(threshold!=0){
	significantline1=H*(-log10(threshold/max(dim(Pmap))))/Max
    s1X=(-significantline1+r+H*i+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
    s1Y=(-significantline1+r+H*i+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
    if(significantline1<H){
	lines(s1X,s1Y,type="l",col=threshold.col,lwd=1)
	}else{
	warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold,"!",sep=""))
	}
	if(amplify == TRUE){
    HX1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
    HY1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
	points(HX1,HY1,pch=19,cex=cex[1],col="white")
	points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=signal.col)
	}
	}
	}
	if(cir.chr==TRUE){
	ticks1=1.1*(2*cir.band+RR)*sin(2*pi*ticks/TotalN)
    ticks2=1.1*(2*cir.band+RR)*cos(2*pi*ticks/TotalN)
	if(cir.labels==TRUE){
    for(i in 1:length(ticks)){
      angle=360*(1-ticks[i]/TotalN)
      text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=1)
	}
	}
	}else{
    ticks1=1.07*(RR+cir.band)*sin(2*pi*ticks/TotalN)
    ticks2=1.07*(RR+cir.band)*cos(2*pi*ticks/TotalN)
	if(cir.labels==TRUE){
    for(i in 1:length(ticks)){
      angle=360*(1-ticks[i]/TotalN)
      text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=1)
    }
	}	
	}
	}
  }
  if(fill.output==TRUE) dev.off()
  print("Circular-Manhattan has been finished!",quote=F)
}
if("m" %in% output){
if(multracks==FALSE){
  print("Starting Rectangular-Manhattan plot!",quote=F)
	for(i in 1:R){
  	print(paste("Rectangular-Plotting ",taxa[i],"...",sep=""),quote=F)
	if(fill.output==TRUE){
    if(fill=="jpg")	jpeg(paste("Rectangular-Manhattan.",taxa[i],".jpg",sep=""), width = 14*dpi,height=5*dpi,res=dpi,quality = 100)
    if(fill=="pdf")	pdf(paste("Rectangular-Manhattan.",taxa[i],".pdf",sep=""), width = 15,height=6)
    if(fill=="tiff")	tiff(paste("Rectangular-Manhattan.",taxa[i],".tiff",sep=""), width = 14*dpi,height=5*dpi,res=dpi)
	}
	if(fill.output==FALSE) {
	dev.new(width = 15, height = 6)
	}
	par(mar = c(5,6,4,3),xaxs="i")
    pvalue=pvalueT[,i]
    logpvalue=logpvalueT[,i]
    if(plotXY==TRUE){
	pvalueXYi.=pvalueXY.[,i]
	if(!is.null(threshold)){
	if(threshold!=0){
    Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.),-log10(threshold/max(dim(Pmap))))
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
	}
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
	}
	if(is.null(ylim)){
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=c(0,Max+1),ylab=ylab,
        cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
		}else{
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=ylim,ylab=ylab,
        cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))	
		}
	points(pvalueXYi.~XYlim,pch=pch,cex=cex[2],col="gray")
	}else{
	if(is.null(ylim)){
    if(!is.null(threshold)){
	if(threshold!=0){
    Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),-log10(threshold/max(dim(Pmap))))
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
	}
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
	}
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+1),ylab=ylab,
         cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab="Chromosome",main=paste("Manhattan plot of",taxa[i]))
	}else{
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+band),ylim=ylim,ylab=ylab,
         cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
	}
	}
    axis(1, at=ticks,cex.axis=cex.axis,font=2,labels=chr)
	if(is.null(ylim)){
    axis(2,at=c(0:(Max+1)),cex.axis=cex.axis,font=2,labels=0:(Max+1))
	}else{
	axis(2,at=c(0:ylim[2]),cex.axis=cex.axis,font=2,labels=0:ylim[2])
	}
	if(!is.null(threshold)){
	if(threshold!=0){
	h=-log10(threshold/max(dim(Pmap)))
	abline(h=h,col=threshold.col)
	if(amplify == TRUE){
    sgline1=-log10(threshold/max(dim(Pmap)))
    HY1=logpvalue[which(logpvalue>=sgline1)]
    HX1=which(logpvalue>=sgline1)
	points(HX1,HY1,pch=pch,cex=cex[2],col="white")
    points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=signal.col)
	}
	}
	}
	if(fill.output==TRUE)  dev.off()
  }
	print("Rectangular-Manhattan has been finished!",quote=F)
  }else{
print("Starting Rectangular-Manhattan plot!",quote=F)
print("Plotting in multiple tracks!",quote=F)
if(fill.output==TRUE){
    if(fill=="jpg")	jpeg(paste("Rectangular-Manhattan.",taxa,".jpg",sep=""), width = 14*dpi,height=5*dpi*R,res=dpi,quality = 100)
    if(fill=="pdf")	pdf(paste("Rectangular-Manhattan.",taxa,".pdf",sep=""), width = 15,height=6*R)
    if(fill=="tiff")	tiff(paste("Rectangular-Manhattan.",taxa,".tiff",sep=""), width = 14*dpi,height=5*dpi*R,res=dpi)
	}
	if(fill.output==FALSE) {
	dev.new(width = 15, height = 6)
	}
	par(mfcol=c(R,1),mar = c(5,6,4,3),xaxs="i")
  for(i in 1:R){
    print(paste("Rectangular-Plotting ",taxa[i],"...",sep=""),quote=F)
    pvalue=pvalueT[,i]
    logpvalue=logpvalueT[,i]
    if(plotXY==TRUE){
	pvalueXYi.=pvalueXY.[,i]
    if(!is.null(threshold)){
	if(threshold!=0){
    Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.),-log10(threshold/max(dim(Pmap))))
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
	}
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
	}
	if(is.null(ylim)){
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=c(0,Max+1),ylab=ylab,
        cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
		}else{
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=ylim,ylab=ylab,
        cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))	
		}
	points(pvalueXYi.~XYlim,pch=pch,cex=cex[2],col="gray")
	}else{
	if(is.null(ylim)){
	if(!is.null(threshold)){
	if(threshold!=0){
    Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),-log10(threshold/max(dim(Pmap))))
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
	}
	}else{
	Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
	}
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+1),ylab=ylab,
         cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
	}else{
	plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(col,N),add),xlim=c(0,length(logpvalue)+band),ylim=ylim,ylab=ylab,
         cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
	}
	}
    axis(1, at=ticks,cex.axis=cex.axis,font=2,labels=chr)
	if(is.null(ylim)){
    axis(2,at=c(0:(Max+1)),cex.axis=cex.axis,font=2,labels=0:(Max+1))
	}else{
	axis(2,at=c(0:ylim[2]),cex.axis=cex.axis,font=2,labels=0:ylim[2])
	}
	if(!is.null(threshold)){
	if(threshold!=0){
	h=-log10(threshold/max(dim(Pmap)))
	abline(h=h,col=threshold.col)
	if(amplify == TRUE){
    sgline1=-log10(threshold/max(dim(Pmap)))
    HY1=logpvalue[which(logpvalue>=sgline1)]
    HX1=which(logpvalue>=sgline1)
	points(HX1,HY1,pch=pch,cex=cex[2],col="white")
    points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=signal.col)
	}
	}
	}
  }
if(fill.output==TRUE) dev.off()
	print("Rectangular-Manhattan has been finished!",quote=F)
  }
  }
  print(paste("The plots have been stored in ","[",getwd(),"]",sep=""),quote=F)
}

