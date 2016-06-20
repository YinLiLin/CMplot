# R-CMplot
CMplot <-

#Adjustable parameters in CMplot
function(
Pmap,
col=c("darkgreen","darkmagenta","navy","black","orange"),
pch=19,
band=1,
cir.band=1,
H=1,
ylim=NULL,
cex.axis=1,
plot.type="b",
multracks=FALSE,
cex=c(0.5,1,1),
r=1,
xlab="Chromosome",
ylab=expression(-log[10](italic(p))),
xaxs="i",
yaxs="r",
outward=TRUE,
threshold = NULL, 
threshold.col="red",
threshold.lwd=1,
threshold.lty=2,
amplify= TRUE,
signal.cex = 1.5,
signal.pch = 19,
signal.col="red",
cir.chr=TRUE,
cir.chr.h=1,
cir.chr.col="black",
cir.chr.labels=NULL,
cir.legend=TRUE,
cir.legend.cex=0.5,
cir.legend.col="grey45",
LOG10=TRUE,
box=FALSE,
conf.int=TRUE,
conf.int.col="grey",
plot0=TRUE,
file.output=TRUE,
file="jpg",
dpi=300
)

{	
	#plot a circle with a radius of r
	circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
	{
		curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
		curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
	}
	
	#print the version of CMplot
	print(paste(paste(rep("*",5),collapse=""),"Welcome to use CMplot!",paste(rep("*",5),collapse=""),sep=""),quote=F)
	print(paste("* ","Version: ",packageVersion("CMplot"),paste(rep(" ",15),collapse=""),"*",sep=""),quote=F)
	print(paste("* "," Author: Lilin Yin",paste(rep(" ",11),collapse=""),"*",sep=""),quote=F)
	print(paste("* ","Contact: ylilin@163.com",paste(rep(" ",6),collapse=""),"*",sep=""),quote=F)
	print(paste(rep("*",32),collapse=""),quote=F)

	if(sum(plot.type %in% "b")==1) plot.type=c("c","m","q")
	plotXY=TRUE
	
	#remove the SNPs of which the pavlues equal to NA
	Pmap=na.omit(Pmap)
	
	#order Pmap by the name of SNP
	Pmap=Pmap[order(Pmap[,1]),]
	
	#delete the column of SNPs names
	Pmap=Pmap[,-1]
	
	taxa=colnames(Pmap)[-c(1,2)]
	
	#scale and adjust the parameters
	cir.chr.h=cir.chr.h/5
	cir.band=cir.band/5
	if(!is.null(threshold)){
		threshold.col=rep(threshold.col,length(threshold))
		threshold.lwd=rep(threshold.lwd,length(threshold))
		threshold.lty=rep(threshold.lty,length(threshold))
	}
	if(length(cex)!=3) cex=rep(cex,3)
	if(!is.null(ylim)){
		if(length(ylim)==1) ylim=c(0,ylim)
	}
	
	#get the number of traits
	R=dim(Pmap)[2]-2
	
	#choose whether to plot the SNPs that have no clear chromosome
	if(plot0==FALSE) index=c(1:100)
	if(plot0==TRUE) index=c(0:100)
	
	#pick the SNPs on euchromosome
	PmapN=Pmap[Pmap[,1] %in% index,]
	
	#pick the SNPs on heterosome
	PmapXY=Pmap[!(Pmap[,1] %in% index)&Pmap[,1]!=0,]
	chrXY=unique(PmapXY[,1])
	
	#check whether there exits heterosome
	if(dim(PmapXY)[1]==0) plotXY=FALSE
	
	#print(plotXY)
	PmapN=matrix(as.numeric(as.matrix(PmapN)),nrow(PmapN))
	if(plotXY==TRUE) PmapXY=as.matrix(PmapXY)
	
	#order the GWAS results by chromosome and position
	orderP=PmapN[order(PmapN[,1],PmapN[,2]),]
	if(plotXY==TRUE) orderPXY.=PmapXY[order(PmapXY[,1],PmapXY[,2]),]
	
	#get the index of chromosome
	if(plotXY==TRUE){
		chr=c(unique(orderP[,1]),"S")
	}else{
		chr=unique(orderP[,1])
	}
	
	pvalueT=as.matrix(orderP[,-c(1:2)])
	
	#change the original pvalues into LOG10 format
	if(plotXY==TRUE){
		if(LOG10 == TRUE){
			pvalueXY.=-log10(matrix(as.numeric(orderPXY.[,-c(1:2)]),nrow(orderPXY.)))
		}else{
			pvalueXY.=matrix(as.numeric(orderPXY.[,-c(1:2)]),nrow(orderPXY.))
		}
	} 
	
	#set the colors for the plot
	#palette(heat.colors(1024)) #(heatmap)
	#T=floor(1024/max(pvalue))
	#plot(pvalue,pch=19,cex=0.6,col=(1024-floor(pvalue*T)))
	if(!missing(col)){
		if(is.vector(col)){
			col=matrix(col,R,length(col),byrow=T)
		}
		if(is.matrix(col)){
			#try to transform the colors into matrix for all traits
			col=matrix(as.vector(t(col)),R,dim(col)[2],byrow=T)
		}
	}else{
		# "byrow=T" must be needed to set 
		col=matrix(c("darkgreen","darkmagenta","navy","black","orange"),R,5,byrow=T)
	}
	
	Num=as.numeric(table(PmapN[,1]))
	Nchr=length(Num)
	N=NULL
	
	#set the colors for each traits
	for(i in 1:R){
		colx=col[i,]
		colx=colx[!is.na(colx)]
		N[i]=ceiling(Nchr/length(colx))
	}
	
	#scale the space parameter between chromosomes
	if(!missing(band)){
		band=floor(band*(dim(pvalueT)[1]/100))
	}else{
		band=floor(dim(pvalueT)[1]/100)
	}
	
	#insert the space into chromosomes and return the midpoint of each chromosome
	ticks=NULL
	NewP=NULL
	for(j in 1:R){
		pvalue=pvalueT[,j]
		for(i in 0:(Nchr-1)){
			if (i==0){
				pvalue=append(pvalue,rep(Inf,band),after=0)
				ticks[i+1]=band+floor((Num[i+1])/2)
			}else{
				pvalue=append(pvalue,rep(Inf,band),after=sum(Num[1:i])+i*band)
				ticks[i+1]=band*(i+1)+sum(Num[1:i])+floor((Num[i+1])/2)
			}
		}
		NewP[[j]]=pvalue
	}
	
	#merge the pvalues of traits by column
	pvalueT=do.call(cbind,NewP)
	if(LOG10 == TRUE){
		logpvalueT=-log10(pvalueT)
	}else{
		logpvalueT=pvalueT
	}
	
	Num=Num+band
	if(plotXY==TRUE){
		ticks=c(ticks,dim(pvalueT)[1]+band+dim(pvalueXY.)[1]/2)
		
		add=list()
		for(i in 1:R){
			colx=col[i,]
			colx=colx[!is.na(colx)]
			add[[i]]=c(Num,rep(0,N[i]*length(colx)-Nchr))
		}
		
		#get the total number of points which need to plot
		XYlim=(dim(pvalueT)[1]+band+1):((dim(pvalueT)[1]+band)+dim(pvalueXY.)[1])	
		TotalN1=dim(pvalueT)[1]
		TotalN2=dim(pvalueXY.)[1]
		TotalN=TotalN1+TotalN2+band
	}else{
		ticks=ticks
		
		add=list()
		for(i in 1:R){
			colx=col[i,]
			colx=colx[!is.na(colx)]
			add[[i]]=c(Num,rep(0,N[i]*length(colx)-Nchr))
		}
		
		TotalN1=dim(pvalueT)[1]
		TotalN=TotalN1
	}
	
	#plot circle Manhattan
	if("c" %in% plot.type){
		#print("Starting Circular-Manhattan plot!",quote=F)
		if(file.output==TRUE){
			if(file=="jpg")	jpeg(paste("Circular-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = 8*dpi,height=8*dpi,res=dpi,quality = 100)
			if(file=="pdf")	pdf(paste("Circular-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = 10,height=10)
			if(file=="tiff")	tiff(paste("Circular-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = 8*dpi,height=8*dpi,res=dpi)
		}
		par(pty="s",xpd=TRUE,mar=c(1,1,1,1))
		RR=r+H*R+cir.band*R
		plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=F,xlab="",ylab="")
		
		for(i in 1:R){
		
			#get the colors for each trait
			colx=col[i,]
			colx=colx[!is.na(colx)]
			
			#debug
			#print(colx)
			
			print(paste("Circular-Plotting ",taxa[i],"...",sep=""),quote=F)
			pvalue=pvalueT[,i]
			logpvalue=logpvalueT[,i]
			if(plotXY==TRUE){
				pvalueXYi.=pvalueXY.[,i]
				if(is.null(ylim)){
					if(LOG10 == TRUE){
						Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),ceiling(max(pvalueXYi.)))
					}else{
						Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
						if(Max<=1)
						Max=max(max(pvalue[pvalue!=Inf]),max(pvalueXYi.))
					}
				}else{
					Max=ylim[2]
				}
				Cpvalue=(H*logpvalue/Max)[-c(1:round(band/2))]
				CpvalueXY=H*pvalueXYi./Max
			}else{
				if(is.null(ylim)){
					if(LOG10 == TRUE){
						Max=ceiling(-log10(min(pvalue[pvalue!=0])))
					}else{
						Max=ceiling(max(pvalue[pvalue!=Inf]))
						if(Max<=1)
						Max=max(pvalue[pvalue!=Inf])
					}
				}else{
					Max=ylim[2]
				}
				Cpvalue=(H*logpvalue/Max)[-c(1:round(band/2))]
			}
			if(outward==TRUE){
				if(cir.chr==TRUE){
					# XLine=(RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)
					circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
					circle.plot(myr=RR,lwd=1.5,add=TRUE)
					a=0
					
					#plot the boundary which represents the chromosomes
					if(plotXY==FALSE){
						for(k in 1:length(chr)){
							if(k==1){
							
								#change the axis from right angle into circle format
								X1chr=(RR)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y1chr=(RR)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								X2chr=(RR+cir.chr.h)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y2chr=(RR+cir.chr.h)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								if(is.null(cir.chr.col)){
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
								}else{
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
								}
							}else{
								a=a+sum(Pmap[,1]==chr[k-1])
								X1chr=(RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y1chr=(RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								X2chr=(RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y2chr=(RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								if(is.null(cir.chr.col)){
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])
								}else{
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
								}		
							}
						}
					}
					if(plotXY==TRUE){
						for(k in 1:(length(chr)-1)){
							if(k==1){
								X1chr=(RR)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y1chr=(RR)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								X2chr=(RR+cir.chr.h)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y2chr=(RR+cir.chr.h)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
									if(is.null(cir.chr.col)){
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
									}else{
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
									}	
							}else{
								a=a+sum(Pmap[,1]==chr[k-1])
								X1chr=(RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y1chr=(RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								X2chr=(RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y2chr=(RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								if(is.null(cir.chr.col)){
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
								}else{
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
								}	
							}
						}
						k=length(chr)
						a=a+sum(Pmap[,1]==chr[k-1])
						X1chr=(RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						Y1chr=(RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						X2chr=(RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						Y2chr=(RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						if(is.null(cir.chr.col)){
							polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="gray",border="gray")	
						}else{
							polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
						}
					}
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							significantline1=H*(-log10(threshold[thr]/max(dim(Pmap))))/Max
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
						significantline1=H*(-log10(min(threshold)/max(dim(Pmap))))/Max
					}
				}
				
				#plot the legend for each trait
				if(cir.legend==TRUE){
				
					#try to get the number after radix point
					if(Max<=1) {
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}else{
						round.n=2
					}
					segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-1)+cir.band*(i-1),H/10,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.75)+cir.band*(i-1),H/10,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.5)+cir.band*(i-1),H/10,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.25)+cir.band*(i-1),H/10,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0)+cir.band*(i-1),H/10,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				}
				X=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				Y=(Cpvalue+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				if(plotXY==TRUE){
					XX=(CpvalueXY+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(((TotalN1-round(band/2))+band+1):(TotalN-round(band/2)))/TotalN)
					YY=(CpvalueXY+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(((TotalN1-round(band/2))+band+1):(TotalN-round(band/2)))/TotalN)
						points(XX,YY,pch=19,cex=cex[1],col="gray")
						points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				}else{
						#plot(chrX,chrY,type="l",col="black",lwd=4)
						points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						if(amplify == TRUE){
							HX1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							HY1=(Cpvalue[which(Cpvalue>=significantline1)]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=19,cex=cex[1],col="white")
							if(is.null(signal.col)){
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=rep(rep(colx,N[i]),add[[i]])[which(Cpvalue>=significantline1)])
							}else{
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=signal.col)
							}
						}
					}
				}
				if(cir.chr==TRUE){
					ticks1=1.07*(RR+cir.chr.h)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.07*(RR+cir.chr.h)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(cir.chr.labels)){
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],cir.chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}else{
					ticks1=(0.9*r)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=(0.9*r)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(cir.chr.labels)){
						for(i in 1:length(ticks)){
						angle=360*(1-(ticks-round(band/2))[i]/TotalN)
						text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],cir.chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}
			}
			if(outward==FALSE){
				if(cir.chr==TRUE){
					# XLine=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
					# YLine=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
					# lines(XLine,YLine,lwd=1.5)
					circle.plot(myr=2*cir.band+RR+cir.chr.h,lwd=1.5,add=TRUE)
					circle.plot(myr=2*cir.band+RR,lwd=1.5,add=TRUE)

					a=0
					if(plotXY==FALSE){
						for(k in 1:length(chr)){
							if(k==1){
								X1chr=(2*cir.band+RR)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y1chr=(2*cir.band+RR)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
									if(is.null(cir.chr.col)){
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
									}else{
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
									}
							}else{
								a=a+sum(Pmap[,1]==chr[k-1])
								X1chr=(2*cir.band+RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y1chr=(2*cir.band+RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								if(is.null(cir.chr.col)){
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
								}else{
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
								}	
							}
						}
					}
					if(plotXY==TRUE){
						for(k in 1:(length(chr)-1)){
							if(k==1){
								X1chr=(2*cir.band+RR)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y1chr=(2*cir.band+RR)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
								Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*((round(band/2)+1):(round(band/2)+sum(Pmap[,1]==chr[1])))/TotalN)
									if(is.null(cir.chr.col)){
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
									}else{
										polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
									}	
							}else{
								a=a+sum(Pmap[,1]==chr[k-1])
								X1chr=(2*cir.band+RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y1chr=(2*cir.band+RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1]==chr[k])))/TotalN)
								if(is.null(cir.chr.col)){
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])	
								}else{
									polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
								}	
							}
						}
						k=length(chr)
						a=a+sum(Pmap[,1]==chr[k-1])
						X1chr=(2*cir.band+RR)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						Y1chr=(2*cir.band+RR)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						X2chr=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						Y2chr=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(((k-0.5)*band+a+1):((k-0.5)*band+a+sum(Pmap[,1] %in% chrXY)))/TotalN)
						if(is.null(cir.chr.col)){
							polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="gray",border="gray")	
						}else{
							polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=cir.chr.col,border=cir.chr.col)
						}
					}
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							significantline1=H*(-log10(threshold[thr]/max(dim(Pmap))))/Max
							#s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
							#s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
							if(significantline1<H){
								#lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
								circle.plot(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
							}else{
								warning(paste("No significant points for ",taxa[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
							}
						}
						significantline1=H*(-log10(min(threshold)/max(dim(Pmap))))/Max
					}
				}
				if(cir.legend==TRUE){
					
					#try to get the number after radix point
					if(Max<=1) {
						round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
					}else{
						round.n=2
					}
					segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-1)+cir.band*(i-1),H/10,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.75)+cir.band*(i-1),H/10,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.5)+cir.band*(i-1),H/10,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0.25)+cir.band*(i-1),H/10,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					segments(0,r+H*(i-0)+cir.band*(i-1),H/10,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
					text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Max*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Max*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Max*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
					text(-r/15,r+H*(i-1)+cir.band*(i-1),round(Max*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
				}
				X=(-Cpvalue+r+H*i+cir.band*(i-1))*sin(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				Y=(-Cpvalue+r+H*i+cir.band*(i-1))*cos(2*pi*(1:(TotalN1-round(band/2)))/TotalN)
				if(plotXY==TRUE){
					XX=(-CpvalueXY+r+H*i+cir.band*(i-1))*sin(2*pi*(((TotalN1-round(band/2))+band+1):(TotalN-round(band/2)))/TotalN)
					YY=(-CpvalueXY+r+H*i+cir.band*(i-1))*cos(2*pi*(((TotalN1-round(band/2))+band+1):(TotalN-round(band/2)))/TotalN)
						points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
						points(XX,YY,pch=19,cex=cex[1],col="gray")
				}else{
						points(X,Y,pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]]))
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						if(amplify == TRUE){
							HX1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*sin(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							HY1=(-Cpvalue[which(Cpvalue>=significantline1)]+r+H*i+cir.band*(i-1))*cos(2*pi*(which(Cpvalue>=significantline1))/TotalN)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=19,cex=cex[1],col="white")
							if(is.null(signal.col)){
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=rep(rep(colx,N[i]),add[[i]])[which(Cpvalue>=significantline1)])
							}else{
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[1],col=signal.col)
							}
							
						}
					}
				}
				if(cir.chr==TRUE){
					ticks1=1.1*(2*cir.band+RR)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.1*(2*cir.band+RR)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(cir.chr.labels)){
						for(i in 1:length(ticks)){
						  angle=360*(1-(ticks-round(band/2))[i]/TotalN)
						  text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],cir.chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}
				}else{
					ticks1=1.0*(RR+cir.band)*sin(2*pi*(ticks-round(band/2))/TotalN)
					ticks2=1.0*(RR+cir.band)*cos(2*pi*(ticks-round(band/2))/TotalN)
					if(is.null(cir.chr.labels)){
						for(i in 1:length(ticks)){
						
							#adjust the angle of labels of circle plot
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],chr[i],srt=angle,font=2,cex=cex.axis)
						}
					}else{
						for(i in 1:length(ticks)){
							angle=360*(1-(ticks-round(band/2))[i]/TotalN)
							text(ticks1[i],ticks2[i],cir.chr.labels[i],srt=angle,font=2,cex=cex.axis)
						}
					}	
				}
			}
		}
		if(file.output==TRUE) dev.off()
		#print("Circular-Manhattan has been finished!",quote=F)
	}

	if("m" %in% plot.type){
		if(multracks==FALSE){
			#print("Starting Rectangular-Manhattan plot!",quote=F)
			for(i in 1:R){
				colx=col[i,]
				colx=colx[!is.na(colx)]
				print(paste("Rectangular-Plotting ",taxa[i],"...",sep=""),quote=F)
					if(file.output==TRUE){
						if(file=="jpg")	jpeg(paste("Rectangular-Manhattan.",taxa[i],".jpg",sep=""), width = 14*dpi,height=5*dpi,res=dpi,quality = 100)
						if(file=="pdf")	pdf(paste("Rectangular-Manhattan.",taxa[i],".pdf",sep=""), width = 15,height=6)
						if(file=="tiff")	tiff(paste("Rectangular-Manhattan.",taxa[i],".tiff",sep=""), width = 14*dpi,height=5*dpi,res=dpi)
						par(mar = c(5,6,4,3),xaxs=xaxs,yaxs=yaxs)
					}
					if(file.output==FALSE) {
						dev.new(width = 15, height = 6)
					}
					pvalue=pvalueT[,i]
					logpvalue=logpvalueT[,i]
					if(plotXY==TRUE){
						pvalueXYi.=pvalueXY.[,i]
						if(!is.null(threshold)){
							if(sum(threshold!=0)==length(threshold)){
								if(LOG10 == TRUE){
									Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),ceiling(max(pvalueXYi.)),ceiling(-log10(min(threshold)/max(dim(Pmap)))))
								}else{
									Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)),ceiling(-log10(min(threshold)/max(dim(Pmap)))))
								}
							}else{
								if( LOG10== TRUE){
									Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),ceiling(max(pvalueXYi.)))
								}else{
									Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
									if(Max<=1)
									# {
										Max=max(max(pvalue[pvalue!=Inf]),max(pvalueXYi.))
									# }else{
										# Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
									# }
								}	
							}
						}else{
							if( LOG10== TRUE){
								Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),ceiling(max(pvalueXYi.)))
							}else{
								Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
								if(Max<=1)
								# {
									Max=max(max(pvalue[pvalue!=Inf]),max(pvalueXYi.))
								# }else{
									# Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
								# }
							}
						}
						if(is.null(ylim)){
							if(Max<=1){
								plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=c(0,Max+10^(-ceiling(-log10(Max)))),ylab=ylab,
									cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
							}else{
								plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=c(0,Max+1),ylab=ylab,
									cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
							}
						}else{
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=ylim,ylab=ylab,
								cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))	
						}
							points(pvalueXYi.~XYlim,pch=pch,cex=cex[2],col="gray")
					}else{
						if(is.null(ylim)){
							if(!is.null(threshold)){
								if(sum(threshold!=0)==length(threshold)){
									if(LOG10 == TRUE){
										Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),ceiling(-log10(min(threshold)/max(dim(Pmap)))))
									}else{
										Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(-log10(min(threshold)/max(dim(Pmap)))))
									}
								}else{
									if(LOG10 == TRUE){
										Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
									}else{
										Max=max(ceiling(max(pvalue[pvalue!=Inf])))
										if(Max<=1)
										#{
											Max=max(max(pvalue[pvalue!=Inf]))
										# }else{
											# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
										# }
									}
								}
							}else{
								if(LOG10 == TRUE){
									Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
								}else{
									Max=max(ceiling(max(pvalue[pvalue!=Inf])))
									if(Max<=1)
									#{
										Max=max(max(pvalue[pvalue!=Inf]))
									# }else{
										# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
									# }
								}
							}
							if(Max<=1){
								plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+10^(-ceiling(-log10(Max)))),ylab=ylab,
									cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
							}else{
								plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+1),ylab=ylab,
									cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
							}
						}else{
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=ylim,ylab=ylab,
								cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE,xlab=xlab,main=paste("Manhattan plot of",taxa[i]))
						}
					}
					axis(1, at=c(0,ticks),cex.axis=cex.axis,font=2,labels=c("Chr",chr))
					if(is.null(ylim)){
						if(Max>1){
							#print(seq(0,(Max+1),round((Max+1)/10)))
							axis(2,at=seq(0,(Max+1),round((Max+1)/10)),cex.axis=cex.axis,font=2,labels=seq(0,(Max+1),round((Max+1)/10)))
						}else{
							axis(2,at=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))),cex.axis=cex.axis,font=2,labels=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))))
						}
					}else{
						if(ylim[2]>1){
							axis(2,at=seq(0,(ylim[2]+1),round((ylim[2]+1)/10)),cex.axis=cex.axis,font=2,labels=seq(0,(ylim[2]+1),round((ylim[2]+1)/10)))
						}else{
							axis(2,at=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))),cex.axis=cex.axis,font=2,labels=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))))
						}
					}
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							for(thr in 1:length(threshold)){
								h=-log10(threshold[thr]/max(dim(Pmap)))
								# print(h)
								# print(threshold.col[thr])
								# print(threshold.lty[thr])
								# print(threshold.lwd[thr])
								abline(h=h,col=threshold.col[thr],lty=threshold.lty[thr],lwd=threshold.lwd[thr])
							}
							if(amplify == TRUE){
								sgline1=-log10(min(threshold)/max(dim(Pmap)))
								HY1=logpvalue[which(logpvalue>=sgline1)]
								HX1=which(logpvalue>=sgline1)
								
								#cover the points that exceed the threshold with the color "white"
								points(HX1,HY1,pch=pch,cex=cex[2],col="white")
								if(is.null(signal.col)){
									points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=rep(rep(colx,N[i]),add[[i]])[HX1])
								}else{
									points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=signal.col)
								}
							}
						}
					}
				if(box==TRUE) box()
				if(file.output==TRUE)  dev.off()
			}
			#print("Rectangular-Manhattan has been finished!",quote=F)
		}else{
			#print("Starting Rectangular-Manhattan plot!",quote=F)
			#print("Plotting in multiple tracks!",quote=F)
			if(file.output==TRUE){
				if(file=="jpg")	jpeg(paste("Rectangular-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = 14*dpi,height=5*dpi*R,res=dpi,quality = 100)
				if(file=="pdf")	pdf(paste("Rectangular-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = 15,height=6*R)
				if(file=="tiff")	tiff(paste("Rectangular-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = 14*dpi,height=5*dpi*R,res=dpi)
			}
			if(file.output==FALSE){
				dev.new(width = 15, height = 6)
			}
			par(mfcol=c(R,1),mar=c(1, 5, 1, 2),oma=c(4,0,3,0),xaxs=xaxs,yaxs=yaxs)
			for(i in 1:R){
				print(paste("Rectangular-Plotting ",taxa[i],"...",sep=""),quote=F)
				pvalue=pvalueT[,i]
				logpvalue=logpvalueT[,i]
				if(plotXY==TRUE){
					pvalueXYi.=pvalueXY.[,i]
					if(!is.null(threshold)){
						if(sum(threshold!=0)==length(threshold)){
							if(LOG10 == TRUE){
								Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.),-log10(min(threshold)/max(dim(Pmap))))
							}else{
								Max=max(ceiling(max(pvalue[pvalue!=Inf])),max(pvalueXYi.),-log10(min(threshold)/max(dim(Pmap))))
							}
						}else{
							if(LOG10==TRUE){
								Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
							}else{
								Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
								if(Max<=1)
								#{
									Max=max(max(pvalue[pvalue!=Inf]),max(pvalueXYi.))
								# }else{
									# Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
								# }
							}
						}
					}else{
						if(LOG10==TRUE){
							Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),max(pvalueXYi.))
						}else{
							Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
							if(Max<=1)
							#{
								Max=max(max(pvalue[pvalue!=Inf]),max(pvalueXYi.))
							# }else{
								# Max=max(ceiling(max(pvalue[pvalue!=Inf])),ceiling(max(pvalueXYi.)))
							# }
						}
					}
					if(is.null(ylim)){
						if(Max<=1){
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=c(0,Max+10^(-ceiling(-log10(Max)))),ylab=ylab,
								cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE)
						}else{
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=c(0,Max+1),ylab=ylab,
								cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE)
						}
					}else{
						plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+3*band+length(pvalueXYi.)),ylim=ylim,ylab=ylab,
							cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE)	
					}
					points(pvalueXYi.~XYlim,pch=pch,cex=cex[2],col="gray")
				}else{
					if(is.null(ylim)){
						if(!is.null(threshold)){
							if(sum(threshold!=0)==length(threshold)){
								if(LOG10 ==TRUE){
									Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))),-log10(min(threshold)/max(dim(Pmap))))
								}else{
									Max=max(ceiling(max(pvalue[pvalue!=Inf])),-log10(min(threshold)/max(dim(Pmap))))
								}
							}else{
								if(LOG10 == TRUE){
									Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
								}else{
									Max=max(ceiling(max(pvalue[pvalue!=Inf])))
									if(Max<=1)
									#{
										Max=max(max(pvalue[pvalue!=Inf]))
									# }else{
										# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
									# }
								}	
							}
						}else{
							if(LOG10 == TRUE){
								Max=max(ceiling(-log10(min(pvalue[pvalue!=0]))))
							}else{
							Max=max(ceiling(max(pvalue[pvalue!=Inf])))
								if(Max<=1)
								#{
									Max=max(max(pvalue[pvalue!=Inf]))
								# }else{
									# Max=max(ceiling(max(pvalue[pvalue!=Inf])))
								# }
							}
						}
						if(Max<=1){
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+10^(-ceiling(-log10(Max)))),ylab=ylab,
								cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE)
						}else{
							plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=c(0,Max+1),ylab=ylab,
								cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE)
						}
					}else{
						plot(logpvalue,pch=pch,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]]),xlim=c(0,length(logpvalue)+band),ylim=ylim,ylab=ylab,
							cex.axis=cex.axis,cex.lab=2,font=2,axes=FALSE)
					}
				}
				
				#add the names of traits on plot  
				text(ticks[1],Max,labels=taxa[i],adj=0,font=3,cex=1.5)
				axis(1, at=c(0,ticks),cex.axis=cex.axis,font=2,labels=c("Chr",chr))
				if(i==1) mtext("Manhattan plot",side=3,padj=-1,font=2,cex=1.5)
				if(is.null(ylim)){
					if(Max>1){
						axis(2,at=seq(0,(Max+1),round((Max+1)/10)),cex.axis=cex.axis,font=2,labels=seq(0,(Max+1),round((Max+1)/10)))
					}else{
						axis(2,at=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))),cex.axis=cex.axis,font=2,labels=seq(0,Max+10^(-ceiling(-log10(Max))),10^(-ceiling(-log10(Max)))))
					}
				}else{
					if(ylim[2]>1){
						axis(2,at=seq(0,(ylim[2]+1),round((ylim[2]+1)/10)),cex.axis=cex.axis,font=2,labels=seq(0,(ylim[2]+1),round((ylim[2]+1)/10)))
					}else{
						axis(2,at=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))),cex.axis=cex.axis,font=2,labels=seq(0,ylim[2]+10^(-ceiling(-log10(ylim[2]))),10^(-ceiling(-log10(ylim[2])))))
					}
				}
				if(!is.null(threshold)){
					if(sum(threshold!=0)==length(threshold)){
						for(thr in 1:length(threshold)){
							h=-log10(threshold[thr]/max(dim(Pmap)))
							abline(h=h,col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
						}
						if(amplify == TRUE){
							sgline1=-log10(min(threshold)/max(dim(Pmap)))
							HY1=logpvalue[which(logpvalue>=sgline1)]
							HX1=which(logpvalue>=sgline1)
							
							#cover the points that exceed the threshold with the color "white"
							points(HX1,HY1,pch=pch,cex=cex[2],col="white")
							if(is.null(signal.col)){
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=rep(rep(colx,N[i]),add[[i]])[HX1])
							}else{
								points(HX1,HY1,pch=signal.pch,cex=signal.cex*cex[2],col=signal.col)
							}
						}
					}
				}
				
			}
			
			#add the labels of X-axis
			mtext(xlab,side=1,padj=2.5,font=2,cex=1.5)
			if(file.output==TRUE) dev.off()
			#print("Rectangular-Manhattan has been finished!",quote=F)
		}
	}
		
	if("q" %in% plot.type){
		#print("Starting QQ-plot!",quote=F)
		for(i in 1:R){
			print(paste("QQ-Plotting ",taxa[i],"...",sep=""),quote=F)
			if(file.output==TRUE){
				if(file=="jpg")	jpeg(paste("QQplot.",taxa[i],".jpg",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi,quality = 100)
				if(file=="pdf")	pdf(paste("QQplot.",taxa[i],".pdf",sep=""), width = 5.5,height=5.5)
				if(file=="tiff")	tiff(paste("QQplot.",taxa[i],".tiff",sep=""), width = 5.5*dpi,height=5.5*dpi,res=dpi)
			}else{
				dev.new(width = 5.5, height = 5.5)
			}
			par(mar = c(5,5,4,2))
			if(plotXY==TRUE){
				P.values=c(as.numeric(PmapN[,i+2]),as.numeric(PmapXY[,i+2]))
			}else{
				P.values=as.numeric(PmapN[,i+2])
			}
			P.values=P.values[!is.na(P.values)]
			if(LOG10==TRUE){
				P.values=P.values[P.values>0]
				P.values=P.values[P.values<=1]
				N=length(P.values)
				P.values=P.values[order(P.values)]
			}else{
				N=length(P.values)
				P.values=P.values[order(P.values,decreasing=TRUE)]
			}
			p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
			log.Quantiles <- -log10(p_value_quantiles)
			if(LOG10==TRUE){
				log.P.values <- -log10(P.values)
			}else{
				log.P.values <- P.values
			}
			plot(NULL, xlim = c(0,max(log.Quantiles)), cex.axis=cex.axis, cex.lab=1.2,ylim=c(0,max(log.P.values)),xlab =expression(Expected~~-log[10](italic(p))), ylab = expression(Observed~~-log[10](italic(p))), main = paste("QQplot of",taxa[i]))
			
			#calculate the confidence interval of QQ-plot
			if(conf.int==TRUE){
				N1=length(log.Quantiles)
				c95 <- rep(NA,N1)
				c05 <- rep(NA,N1)
				for(j in 1:N1){
					xi=ceiling((10^-log.Quantiles[j])*N)
					if(xi==0)xi=1
					c95[j] <- qbeta(0.95,xi,N-xi+1)
					c05[j] <- qbeta(0.05,xi,N-xi+1)
				}
				index=length(c95):1
				
				#plot the confidence interval of QQ-plot
				polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col,border=conf.int.col)
			}
			
			if(!is.null(threshold.col))	abline(a = 0, b = 1, col = threshold.col[1],lwd=2)
			points(log.Quantiles, log.P.values, col = col[1],pch=19,cex=cex[3])
			
			if(!is.null(threshold)){
				if(sum(threshold!=0)==length(threshold)){
					thre.line=-log10(min(threshold)/N)
					if(amplify==TRUE){
						thre.index=which(log.P.values>=thre.line)
						if(length(thre.index)!=0){
						
							#cover the points that exceed the threshold with the color "white"
							points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
							if(is.null(signal.col)){
								points(log.Quantiles[thre.index],log.P.values[thre.index],col = col[1],pch=signal.pch,cex=signal.cex)
							}else{
								points(log.Quantiles[thre.index],log.P.values[thre.index],col = signal.col,pch=signal.pch,cex=signal.cex)
							}
						}
					}
				}
			}
			if(file.output==TRUE) dev.off()
		}
	}
	print(paste("The plots have been stored in ","[",getwd(),"]",sep=""),quote=F)
}
