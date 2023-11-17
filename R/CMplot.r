CMplot <- function(
    Pmap,
    col=c("#4197d8", "#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"),
    bin.size=1e6,
    bin.breaks=NULL,
    LOG10=TRUE,
    pch=19,
    type="p",
    band=1,
    H=1.5,
    ylim=NULL,
    axis.cex=1,
    axis.lwd=1.5,
    lab.cex=1.5,
    lab.font=2,
    plot.type=c("m","c","q","d"),
    multracks=FALSE,
    multracks.xaxis=FALSE,
    multraits=FALSE,
    points.alpha=100L,
    r=0.3,
    cex=c(0.5,1,1),
    outward=FALSE,
    ylab=expression(-log[10](italic(p))),
    ylab.pos=3,
    xticks.pos=1,
    mar=c(3,6,3,3),
    mar.between=0,
    threshold=NULL, 
    threshold.col="red",
    threshold.lwd=1,
    threshold.lty=2,
    amplify= TRUE,
    signal.cex=1.5,
    signal.pch=19,
    signal.col=NULL,
    signal.line=2,
    highlight=NULL,
    highlight.cex=1,
    highlight.pch=19,
    highlight.type="p",
    highlight.col="red",
    highlight.text=NULL,
    highlight.text.col="black",
    highlight.text.cex=1,
    highlight.text.font=3,
    chr.labels=NULL,
    chr.border=FALSE,
    chr.labels.angle=0,
    chr.den.col="black",
    chr.pos.max=FALSE,
    cir.band=1,
    cir.chr=TRUE,
    cir.chr.h=1.5,
    cir.axis=TRUE,
    cir.axis.col="black",
    cir.axis.grid=TRUE,
    conf.int=TRUE,
    conf.int.col=NULL,
    file.output=TRUE,
    file.name="",
    file=c("jpg","pdf","tiff","png"),
    dpi=300,
    height=NULL,
    width=NULL,
    main=NULL,
    main.cex=1.5,
    main.font=2,
    legend.ncol=NULL,
    legend.cex=1,
    legend.pos=c("left","middle","right"),
    box=FALSE,
    verbose=TRUE
)
{   

    #plot a circle with a radius of r
    circle.plot <- function(myr,type="l",x=NULL,lty=1,lwd=1,col="black",add=TRUE,n.point=1000)
    {
        curve(sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=add)
        curve(-sqrt(myr^2-x^2),xlim=c(-myr,myr),n=n.point,ylim=c(-myr,myr),type=type,lty=lty,col=col,lwd=lwd,add=TRUE)
    }

    highlight_text <- function(
        x,
        y,
        words=NULL,
        point.cex=1,
        text.cex=1,
        pch=19,
        type = "p",
        point.col = "red",
        text.col = "black",
        text.font=3,
        xlim=c(-Inf, Inf),
        ylim=c(-Inf, Inf)
    )
    {
        overlap <- function(x1, y1, sw1, sh1, boxes) {
            if (length(boxes) == 0) return(FALSE)
            for (i in c(1:length(boxes))) {
                bnds <- boxes[[i]]
                x2 <- bnds[1]
                y2 <- bnds[2]
                sw2 <- bnds[3]
                sh2 <- bnds[4]

                if (x1 < x2)
                    overlap <- x1 + sw1 > x2
                else
                    overlap <- x2 + sw2 > x1

                if (y1 < y2)
                    overlap <- overlap && (y1 + sh1 > y2)
                else
                    overlap <- overlap && (y2 + sh2 > y1)

                if (overlap) {
                    return(TRUE)
                }
            }
            return(FALSE)
        }
    
        layout <- function(x, y, words, cex=1, xlim=c(-Inf, Inf), ylim=c(-Inf, Inf), tstep = .1, rstep = .1) {
            sdx <- sd(x, na.rm=TRUE)
            sdy <- sd(y, na.rm=TRUE)
            if (sdx == 0) sdx <- 1
            if (sdy == 0) sdy <- 1
            boxes <- list()
            for(i in 1:length(words)){
                wid <- strwidth(words[i], cex=cex[i])
                ht <- strheight(words[i], cex=cex[i])
                if(i <= (length(words) / 2)){
                    boxes[[length(boxes) + 1]] <- c(x[i]-0.5*wid, y[i]-0.5*ht, wid, ht)
                }else{
                    xupdt <- xrot <- x[i]
                    yupdt <- yrot <- y[i]
                    r <- 0
                    theta <- runif(1, 0, 2 * pi)
                    ht <- 1.5 * ht
                    isOverlaped <- TRUE
                    while(isOverlaped){
                        if(
                            !overlap(xupdt-0.5*wid, yupdt-0.5*ht, wid, ht, boxes) &&
                            (xupdt-0.5*wid) > xlim[1] &&
                            (yupdt-0.5*ht) > ylim[1] &&
                            (xupdt+0.5*wid) < xlim[2] &&
                            (yupdt+0.5*ht) < ylim[2]
                        ){
                            boxes[[length(boxes) + 1]] <- c(xupdt-0.5*wid, yupdt-0.5*ht, wid, ht)
                            isOverlaped <- FALSE
                        }else{
                            theta <- theta + tstep
                            r <- r + rstep * tstep / (2 * base::pi)
                            xupdt <- xrot + sdx * r * cos(theta)
                            yupdt <- yrot + sdy * r * sin(theta)
                        }
                    }
                }
            }
            result <- do.call(rbind, boxes)
            colnames(result) <- c("x", "y", "width", "ht")
            rownames(result) <- words
            result
        }

        if(!is.null(words)){
            if(length(x) != length(words))  stop("length of highlighted labels is not equal to the highlighted SNPs.")
            indx <- order(y, decreasing=TRUE)
            x <- x[indx]
            y <- y[indx]
            words <- words[indx]
            if(length(point.cex)!=1){if(length(point.cex)==length(x)){point.cex=point.cex[indx]}else{stop("unequal length of 'cex' for highlighted points.")}}else{point.cex=rep(point.cex,length(x))}
            if(length(pch)!=1){if(length(pch)==length(x)){pch=pch[indx]}else{stop("unequal length of 'pch' for highlighted points.")}}else{pch=rep(pch,length(x))}
            if(length(point.col)!=1){if(length(point.col)==length(x)){point.col=point.col[indx]}else{stop("unequal length of 'col' for highlighted points.")}}else{point.col=rep(point.col,length(x))}
            if(length(text.col)!=1){if(length(text.col)==length(x)){text.col=text.col[indx]}else{stop("unequal length of 'col' for highlighted text.")}}else{text.col=rep(text.col,length(x))}
            if(length(text.cex)!=1){if(length(text.cex)==length(x)){text.cex=text.cex[indx]}else{stop("unequal length of 'cex' for highlighted text.")}}else{text.cex=rep(text.cex,length(x))}
            
            words_ety <- words[words == "" | is.na(words)]
            if(length(words_ety)){
                logical_idx <- words == "" | is.na(words)
                if(type=="h"){
                    points(x[logical_idx],y[logical_idx],pch=pch[logical_idx],type="h",col=point.col[logical_idx], lwd=point.cex[logical_idx]+1)
                    points(x[logical_idx],y[logical_idx],pch=pch[logical_idx],type="p",col=point.col[logical_idx], cex=point.cex[logical_idx])
                }else if(type=="l"){
                    segments(x[logical_idx], ylim[1], x[logical_idx], ylim[2], col=point.col[logical_idx], lwd=point.cex[logical_idx], lty=2)
                }else{
                    points(x[logical_idx],y[logical_idx],pch=pch[logical_idx],type="p",col=point.col[logical_idx],cex=point.cex[logical_idx])
                }
                words <- words[!logical_idx]
                x <- x[!logical_idx]
                y <- y[!logical_idx]
                point.cex <- point.cex[!logical_idx]
                pch <- pch[!logical_idx]
                point.col <- point.col[!logical_idx]
                text.col <- text.col[!logical_idx]
                text.cex <- text.cex[!logical_idx]
            }

            x1 <- x
            y1 <- y
            xadj <- rep(c(1.5, 0, -0.5), length=length(x))
            yadj <- rep(c(1.5, 0, -0.5), length=length(x))
            for(i in 1:length(x)){
                if(xadj[i] == 0){
                    if(yadj[i] == -0.5){
                        if((y[i] + 2*strheight(words[i],cex=text.cex)) > max(ylim)){
                            y[i] = y[i] - 1.5*strheight(words[i],cex=text.cex)
                        }else{
                            y[i] = y[i] + 1.5*strheight(words[i],cex=text.cex)
                        }
                    }
                    if(yadj[i] == 1.5)  y[i] = y[i] - 1.5*strheight(words[i],cex=text.cex)
                }else{
                    if(yadj[i] == -0.5){
                        if((y[i] + 1.5*strheight(words[i],cex=text.cex)) > max(ylim)){
                            y[i] = y[i] - strheight(words[i],cex=text.cex)
                        }else{
                            y[i] = y[i] + strheight(words[i],cex=text.cex)
                        }
                    }
                    if(yadj[i] == -0.5) y[i] = y[i] + strheight(words[i],cex=text.cex)
                    if(yadj[i] == 1.5)  y[i] = y[i] - strheight(words[i],cex=text.cex)
                }
                if(xadj[i] == 1.5){
                    if((x[i] - 1.2*strwidth(words[i],cex=text.cex)) < min(xlim)){
                        x[i] = x[i] + 0.6*strwidth(words[i],cex=text.cex)
                    }else{
                        x[i] = x[i] - 0.6*strwidth(words[i],cex=text.cex)
                    }
                }
                if(xadj[i] == -0.5){
                    if((x[i] + 1.2*strwidth(words[i],cex=text.cex)) > max(xlim)){
                        x[i] = x[i] - 0.6*strwidth(words[i],cex=text.cex)
                    }else{
                        x[i] = x[i] + 0.6*strwidth(words[i],cex=text.cex)
                    }
                }
            }

            x <- c(x1,x)
            y <- c(y1,y)
            words <- c(rep("OO", length(words)), as.character(words))
            lay <- layout(x=x,y=y,words=words,cex=c(rep(text.cex[1],length(x1)),text.cex),xlim=xlim,ylim=ylim)
            n <- length(x1)
            indd <- (n+1):length(x)
            for(i in indd){
                xl <- lay[i,1]
                yl <- lay[i,2]
                w <- lay[i,3]
                h <- lay[i,4]
                nx <- xl + 0.5 * w
                ny <- yl + 0.5 * h
                if((nx + 0.5 * strwidth(words[i],cex=text.cex[i-n])) < x1[i-n]){
                    nx=nx + 0.5 * strwidth(words[i],cex=text.cex[i-n])
                }else if((nx - 0.5 * strwidth(words[i],cex=text.cex[i-n])) > x1[i-n]){
                    nx=nx - 0.5 * strwidth(words[i],cex=text.cex[i-n])
                }
                if((ny + strheight(words[i],cex=text.cex[i-n])) < y1[i-n]){
                    ny=ny + 0.5 * strheight(words[i],cex=text.cex[i-n])
                }else if((ny - strheight(words[i],cex=text.cex[i-n])) > y1[i-n]){
                    ny=ny - 0.5 * strheight(words[i],cex=text.cex[i-n])
                }
                # arrows(x1[i-n], y1[i-n], nx, ny, length=.08, angle=15, code=2, col="grey", lwd=2)
                segments(x1[i-n], y1[i-n], nx, ny, col="black", lwd=text.cex[i-n])
            }
            if(type=="h"){
                points(x1,y1,pch=pch,type="h",col=point.col, lwd=point.cex+1)
                points(x1,y1,pch=pch,type="p",col=point.col, cex=point.cex)
            }else if(type=="l"){
                segments(x1, ylim[1], x1, ylim[2], col=point.col, lwd=point.cex, lty=2)
                # points(x1,y1,pch=pch,type="p",col=point.col, cex=point.cex)
            }else{
                points(x1,y1,pch=pch,type=type,col=point.col,cex=point.cex)
            }
            text(lay[indd,1]+0.5*lay[indd,3],lay[indd,2]+0.5*lay[indd,4],words[indd],xpd=TRUE,cex=text.cex,col=text.col,font=text.font)
        }else{
            if(type=="h"){
                points(x,y,pch=pch,type="h",col=point.col, lwd=point.cex+1)
                points(x,y,pch=pch,type="p",col=point.col, cex=point.cex)
            }else if(type=="l"){
                segments(x, ylim[1], x, ylim[2], col=point.col, lwd=point.cex, lty=2)
                # points(x,y,pch=pch,type="p",col=point.col, cex=point.cex)
            }else{
                points(x,y,pch=pch,type=type,col=point.col,cex=point.cex)
            }
        }
    }

    max_ylim <- function(x){
        if(x == 0) return(x)
        if(abs(x) >= 1){
            return(ceiling(x))
        }else{
            if(x < 0){
                digit <- 10^(-ceiling(-log10(abs(x))))
                return(-(floor(abs(x) / digit - 1) * digit))
            }else{
                digit <- 10^(-ceiling(-log10(x)))
                return((floor(x / digit + 1) * digit))
            }
        }
    }

    min_ylim <- function(x){
        if(x == 0) return(x)
        if(abs(x) >= 1){
            return(floor(x))
        }else{
            if(x < 0){
                digit <- 10^(-ceiling(-log10(abs(x))))
                return(-(floor(abs(x) / digit + 1) * digit))
            }else{
                digit <- 10^(-ceiling(-log10(x)))
                return((floor(x / digit - 1) * digit))
            }
        }
    }

    min_no_na <- function(x){
        return(min(x, na.rm=TRUE))
    }

    max_no_na <- function(x){
        return(max(x, na.rm=TRUE))
    }

    # created by Haohao Zhang
    filter.points <- function(x, y, w, h, dpi, scale=1) {
        x <- ceiling((x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)) * w * dpi / scale)
        y <- ceiling((y - min(y, na.rm=TRUE)) / (max(y, na.rm=TRUE) - min(y, na.rm=TRUE)) * h * dpi / scale)
        index <- !duplicated(cbind(x, y))
    }

    DensityPlot <- function(
        chr,
        pos,
        chr.orig.labels,
        col=c("darkgreen", "yellow", "red"),
        main=NULL,
        main.cex=1.2,
        main.font=2,
        chr.labels=NULL, 
        chr.pos.max=FALSE,
        bin=1e6,
        bin.breaks=NULL,
        band=3,
        width=5,
        legend.cex=1,
        legend.y.intersp=1,
        legend.x.intersp=1,
        xticks.pos=1,
        plot=TRUE,
        dpi=NULL,
        wh=NULL,
        ht=NULL
    )
    {
        legend.min <- 1
        legend.max <- NULL
        if(is.null(legend.cex)) legend.cex = 1
        if(!is.null(bin.breaks)){
            bin.breaks <- sort(bin.breaks)
            if(sum(bin.breaks < 0)) stop("breaks should not contain a negative value.")
            if(bin.breaks[1]){
                legend.min <- bin.breaks[1]
            }else{
                bin.breaks <- bin.breaks[-1]
            }
            legend.max <- bin.breaks[length(bin.breaks)]
        }
        if(is.null(col) | length(col) == 1){col=c("darkgreen", "yellow", "red")}
        max.chr <- max(chr)
        chr.num <- unique(chr)
        chorm.maxlen <- max(pos)
        bp <- ifelse(chorm.maxlen < 1e3, 1, ifelse(chorm.maxlen < 1e6, 1e3, 1e6))
        bp_label <- ifelse(bp == 1, "bp", ifelse(bp == 1e3, "Kb", "Mb"))
        if(is.null(main))   main <- paste("The number of SNPs within ", bin / bp, bp_label, " window size", sep="")
        if(plot)    plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main, cex.main=main.cex, font.main=main.font, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
        pos.x <- list()
        chr.pos.max.v <- NULL
        col.index <- list()
        maxbin.num <- NULL
        windinfo <- list()
        for(i in 1 : length(chr.num)){
            pos.x[[i]] <- pos[chr == chr.num[i]]
            maxposindx <- which.max(pos.x[[i]])
            max.pos <- pos.x[[i]][maxposindx]
            chr.pos.max.v <- c(chr.pos.max.v, max.pos)
            cut.breaks <- seq(0, max.pos, bin)
            cut.len <- length(cut.breaks)
            if(cut.breaks[length(cut.breaks)] < max.pos)  cut.breaks <- c(cut.breaks, cut.breaks[length(cut.breaks)] + bin)
            if(chr.pos.max){
                pos.x[[i]] <- pos.x[[i]][-maxposindx]
            }
            if(cut.len <= 1){
                maxbin.num <- c(maxbin.num, length(pos.x[[i]]))
                col.index[[i]] <- rep(length(pos.x[[i]]), length(pos.x[[i]]))
                names(col.index[[i]]) <- 1
            }else{
                cut.r <- cut(pos.x[[i]], cut.breaks, labels=FALSE)
                eachbin.num <- table(cut.r)
                maxbin.num <- c(maxbin.num, max(eachbin.num))
                col.index[[i]] <- rep(eachbin.num, eachbin.num)
            }
            if(plot){
                windinfo <- c(windinfo, tapply(pos.x[[i]], as.numeric(names(col.index[[i]])), function(x){
                    return(c(ifelse(!is.null(chr.labels), chr.labels[i], chr.orig.labels[i]),
                        min(x),max(x),length(x)))})
                )
            }
        }
        if(plot){
            windinfo <- as.data.frame(do.call(rbind, windinfo))
            colnames(windinfo) <- c("Chr", "Start", "End", "Num")
            rownames(windinfo) <- NULL
            for(i in 2:ncol(windinfo)){windinfo[, i]<-as.numeric(windinfo[, i])}
        }
        Maxbin.num <- max(maxbin.num)
        maxbin.num <- Maxbin.num
        if(!is.null(legend.max)){
            maxbin.num <- legend.max
        }
        if(Maxbin.num < legend.min)    stop("the maximum number of markers in windows is smaller than the lower boundary of breaks.")
        col=colorRampPalette(col)(maxbin.num - legend.min + 1)
        col.seg=NULL
        for(i in 1 : length(chr.num)){
            if(plot){
                polygon(c(0, 0, chr.pos.max.v[i], chr.pos.max.v[i]), 
                c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
                width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey95", border="grey95")
                rect(xleft=0, ybottom = -width/5 - band * (i - length(chr.num) - 1), xright=chr.pos.max.v[i], ytop=width/5 - band * (i - length(chr.num) - 1), border="grey80")
            }
            if(!is.null(legend.max)){
                if(legend.max < Maxbin.num){
                    col.index[[i]][col.index[[i]] > legend.max] <- legend.max
                }
            }
            col.index[[i]][col.index[[i]] < legend.min] <- legend.min
            if(!plot)   col.seg <- c(col.seg, col[col.index[[i]] - legend.min + 1])
            if(!is.null(ht) && !is.null(wh) && !is.null(dpi)){
                is_visable <-  filter.points(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), wh * (max(pos.x[[i]])/chorm.maxlen), ht, dpi=dpi)
                if(plot)    segments(pos.x[[i]][is_visable], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]][is_visable], width/5 - band * (i - length(chr.num) - 1), 
                                col=col[col.index[[i]][is_visable] - legend.min + 1], lwd=1)
            }else{
                if(plot)    segments(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]], width/5 - band * (i - length(chr.num) - 1), 
                                col=col[col.index[[i]] - legend.min + 1], lwd=1)
            }
        }
        
        chr.num <- rev(chr.orig.labels)
        if(plot){
            if(!is.null(chr.labels)){
                mtext(at=seq(band, length(chr.num) * band, band), text=chr.labels, side=2, las=2, font=1, cex=axis.cex*0.6, line=0.2, xpd=TRUE)
            }else{
                if(max.chr == 0)    mtext(at=seq(band, length(chr.num) * band, band), text=chr.num, side=2, las=2, font=1, cex=axis.cex*0.6, line=0.2, xpd=TRUE)
                if(max.chr != 0)    mtext(at=seq(band, length(chr.num) * band, band), text=paste("Chr", chr.num, sep=""), side=2, las=2, font=1, cex=axis.cex*0.6, line=0.2, xpd=TRUE)
            }
        }
        if(plot){
            xticks=seq(0, chorm.maxlen / bp, length=10)
            
            if(round(xticks[2]) <= 10){
                xticks=seq(0, chorm.maxlen / bp, round(xticks[2], 1))
            }else{
                xticks=seq(0, chorm.maxlen / bp, round(xticks[2]))    
            }
            
            if((chorm.maxlen/bp - max(xticks)) > 0.5*xticks[2]){
                xticks=c(xticks, round(chorm.maxlen / bp))
            }
            axis(3, mgp=c(3,xticks.pos,0), at=xticks*bp, labels=paste(xticks, bp_label, sep=""), font=1, cex.axis=axis.cex*0.8, tck=0.01, lwd=axis.lwd, padj=1.2)
            axis(3, at=c(0, chorm.maxlen), labels=c("",""), tcl=0, lwd=axis.lwd)
        }

        if(is.null(bin.breaks)){
            legend.len <- 10
            if(maxbin.num <= legend.len)    legend.len <- maxbin.num
            legend.y <- round(seq(0, maxbin.num, length=legend.len + 1))
            legend.y <- unique(legend.y)
            len <- ifelse(length(legend.y)==1, 1, legend.y[2])
            legend.y <- seq(legend.y[2], maxbin.num, len)
        }else{
            legend.y <- bin.breaks
        }
        
        if(!is.null(bin.breaks)){
            if(legend.max < Maxbin.num){
                legend.y[length(legend.y)] <- paste(">=", maxbin.num, sep="")
                legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
            }else{
                legend.y.col <- legend.y
            }
        }else{
            legend.y.col <- legend.y
        }
        if(legend.min != 1){
            legend.y[1] <- paste("<=", legend.min, sep="")
        }
        legend.y <- c("0", legend.y)
        legend.y.col <- as.numeric(legend.y.col)
        legend.col <- c("grey95", col[legend.y.col - legend.min + 1])
        if(plot){
            legend(x=(chorm.maxlen + chorm.maxlen/50), y=(-width/2.5 + band), title="", legend=legend.y, pch=15, pt.cex=legend.cex*3, col=legend.col,
            cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)
            return(windinfo)
        }else{
            return(list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y))
        }
    }

    if(!all(plot.type %in% c("c","m","q","d"))) stop("unknown 'plot.type'.")
    legend.pos <- match.arg(legend.pos)
    file <- match.arg(file)
    trait <- colnames(Pmap)[-c(1:3)]
    if(length(trait) == 0)   trait <- paste("Trait", 1:(ncol(Pmap)-3), sep="")
    taxa <- paste(trait, collapse="_")
    
    if(length(points.alpha) != 1L)   stop("invalid 'points.alpha': must be 'TRUE', 'FALSE' or an integer between 0 and 255")
    if(is.logical(points.alpha))   points.alpha <- ifelse(points.alpha, formals()$points.alpha, 255L)
    if(!is.integer(points.alpha)){
      if(is.numeric(points.alpha) && points.alpha == as.integer(points.alpha))   points.alpha <- as.integer(points.alpha)
      else   stop("invalid 'points.alpha': must an integer between")
    }
    if(!is.integer(points.alpha))    stop("invalid 'points.alpha': must an integer between")
    if(points.alpha < 0L || points.alpha > 255L)   stop("out-of range 'points.alpha': must be between 0 and 255")

    #get the number of traits
    R=ncol(Pmap)-3

    #remove illegal SNPs
    suppressWarnings(Pmap <- Pmap[Pmap[, 2] != "0", ])
    Pmap <- as.matrix(Pmap)
    Pmap <- Pmap[!is.na(Pmap[, 2]), ]
    suppressWarnings(Pmap <- Pmap[!is.na(as.numeric(Pmap[, 3])), ])

    #replace the non-euchromosome
    suppressWarnings(numeric.chr <- as.numeric(Pmap[, 2]))
    suppressWarnings(max.chr <- max(numeric.chr, na.rm=TRUE))
    if(is.infinite(max.chr))    max.chr <- 0
    suppressWarnings(map.xy.index <- which(!numeric.chr %in% c(0:max.chr)))
    if(length(map.xy.index) != 0){
        chr.xy <- unique(Pmap[map.xy.index, 2])
        for(i in 1:length(chr.xy)){
            Pmap[Pmap[, 2] == chr.xy[i], 2] <- max.chr + i
        }
    }
    SNP_id <- Pmap[,1]

    #delete the column of SNPs names
    Pmap <- Pmap[, -1]
    Pmap <- apply(Pmap, 2, as.numeric)
    order_index <- order(Pmap[, 1], Pmap[,2])

    #order the GWAS results by chromosome and position
    Pmap <- Pmap[order_index, ]
    SNP_id <- SNP_id[order_index]

    chr <- unique(Pmap[,1])
    chr.ori <- chr
    if(length(map.xy.index) != 0){
        for(i in 1:length(chr.xy)){
            chr.ori[chr.ori == max.chr + i] <- chr.xy[i]
        }
    }

    #SNP-Density plot
    wind_snp_num <- NULL
    if("d" %in% plot.type){
        if(verbose) cat(" Marker density plotting.\n")
        if(file.output){
            ht=ifelse(is.null(height), 6, height)
            wh=ifelse(is.null(width), 9, width)
            if(file=="jpg") jpeg(paste("Marker_Density.",ifelse(file.name=="",taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
            if(file=="pdf") pdf(paste("Marker_Density.",ifelse(file.name=="",taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
            if(file=="tiff")    tiff(paste("Marker_Density.",ifelse(file.name=="",taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
            if(file=="png") png(paste("Marker_Density.",ifelse(file.name=="",taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
            # par(xpd=TRUE)
            par(mar=c(mar[1]-2, mar[2]-1, mar[3]+1, mar[4]))
        }else{
            ht=ifelse(is.null(height), 6, height)
            wh=ifelse(is.null(width), 9, width)
            if(is.null(dev.list())) dev.new(width=wh,height=ht)
            # par(xpd=TRUE)
        }
        wind_snp_num <- DensityPlot(Pmap[, 1], Pmap[, 2], chr.ori, chr.pos.max=chr.pos.max, dpi=dpi, wh=wh, ht=ht, chr.labels=chr.labels, col=chr.den.col, bin=bin.size, bin.breaks=bin.breaks, main=main[1], main.cex=main.cex, main.font=main.font, legend.cex=legend.cex, xticks.pos=xticks.pos)
        if(file.output) dev.off()
    }
    
    if(length(plot.type) > 1 | (!"d" %in% plot.type)){

        #scale and adjust the parameters
        cir.chr.h <- cir.chr.h/5
        cir.band <- cir.band/5
        if(!is.null(threshold)){
            if(!is.list(threshold)){
                thresholdlist <- list()
                for(i in 1:R){
                    thresholdlist[[i]]  <- threshold
                }
                threshold <- thresholdlist
            }

            if(LOG10){
                if(sum(unlist(threshold) <= 0) != 0) stop("threshold must be greater than 0.")
            }

            threshold.col <- rep(threshold.col, max(sapply(threshold, length)))
            threshold.lwd <- rep(threshold.lwd, max(sapply(threshold, length)))
            threshold.lty <- rep(threshold.lty, max(sapply(threshold, length)))
            signal.col <- rep(signal.col, max(sapply(threshold, length)))
            signal.pch <- rep(signal.pch, max(sapply(threshold, length)))
            signal.cex <- rep(signal.cex, max(sapply(threshold, length)))
        }
        if(length(cex)!=3) cex <- rep(cex,3)

        if(!is.null(ylim)){
            if(!is.list(ylim)){
                if(R > 1)    cat(" (warning: all phenotypes will use the same ylim.)\n")
                if(length(ylim)!=2) stop("ylim for each phenotype should be assigned two values.")
                if(ylim[2] <= ylim[1])  stop("second value should be larger than the first in ylim.")
                ylimlist <- list()
                for(i in 1:R){
                    ylimlist[[i]]  <- ylim
                }
                ylim <- ylimlist
            }else{
                if(length(ylim)!=R) stop("length of list of ylim should equal to the number of phenotype.")
                for(i in 1:R){
                    if(length(ylim[[i]])!=2) stop("ylim for each phenotype should be assigned two values.") 
                    if(ylim[[i]][2] <= ylim[[i]][1])  stop("second value should be larger than the first in ylim.")
                }
            }
        }
        
        if(!is.null(conf.int.col)) conf.int.col <- rep(conf.int.col, R)
        if(!is.null(main)) main <- rep(main, R)
        if(length(mar) != 4)    stop("length of 'mar' shoud equal to 4.")
        if(chr.labels.angle > 90 | chr.labels.angle < -90)  stop("'chr.labels.angle' should be > -90 and < 90.")
        pch=rep(pch, R)
        
        if(!is.null(highlight)){
            highlight_index <- list()
            highlight_col <- list()
            if(is.list(highlight.col)){
                if(length(highlight.col) != R){stop("length of 'highlight.col' not equals to the number of traits.")}
                highlight_col=highlight.col
            }
            if(!is.list(highlight)){
                highlight <- list(highlight)
                for(i in 1:R){highlight[[i]] = highlight[[1]]}
            }else{
                if(length(highlight) != R){stop("length of 'highlight' not equals to the number of traits.")}  
            }
            length(highlight_index) <- length(highlight)
            for(i in 1:length(highlight)){
                if(sum(!is.na(highlight[[i]])) == 0 | length(highlight[[i]]) == 0){
                    highlight_index[[i]] <- NA
                    highlight_col[[i]] <- NA
                }else{
                    highlight[[i]] <- highlight[[i]][!is.na(highlight[[i]])]
                    highlight_index[[i]] <- match(as.character(as.matrix(highlight[[i]])), SNP_id)
                    if(all(is.na(highlight_index[[i]]))) stop("No shared SNPs between Pmap and highlight!")
                    highlight_index[[i]] <- na.omit(highlight_index[[i]])
                    if(!is.null(highlight.col) && !is.list(highlight.col))  highlight_col[[i]] <- highlight.col
                }
            }
        }

        if(!is.null(highlight.text)){
            if(!is.list(highlight.text)){
                highlight.text <- list(highlight.text)
                for(i in 1:R){highlight.text[[i]] = highlight.text[[1]]}
            }else{
                if(length(highlight.text) != R){stop("length of 'highlight.text' not equals to the number of traits.")}  
            }
        }

        pvalueT <- as.matrix(Pmap[,-c(1:2)])
        pvalue.pos <- Pmap[, 2]
        pvalue.pos.list <- tapply(pvalue.pos, Pmap[, 1], list)

        #scale the space parameter between chromosomes
        if(!missing(band)){
            band <- floor(band*(sum(sapply(pvalue.pos.list, max))/100))
        }else{
            band <- floor((sum(sapply(pvalue.pos.list, max))/100))
        }
        if(band==0) band=1
        
        if(LOG10){
            if(sum(pvalueT <= 0, na.rm=TRUE) != 0 || sum(pvalueT > 1, na.rm=TRUE) != 0) stop("p values should be at range of (0, 1).")
            pvalueT[pvalueT <= 0] <- NA
            pvalueT[pvalueT > 1] <- NA
        }
        Pmap[,-c(1:2)] <- pvalueT

        #set the colors for the plot
        if(is.vector(col)){
            col <- matrix(col,R,length(col),byrow=TRUE)
        }
        if(is.matrix(col)){
            #try to transform the colors into matrix for all traits
            col <- matrix(as.vector(t(col)),R,dim(col)[2],byrow=TRUE)
        }

        Num <- as.numeric(table(Pmap[,1]))
        Nchr <- length(Num)
        N <- NULL

        #set the colors for each traits
        for(i in 1:R){
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            N[i] <- ceiling(Nchr/length(colx))
        }
        
        #insert the space into chromosomes and return the midpoint of each chromosome
        ticks <- NULL
        chr.border.pos <- NULL
        pvalue.posN <- NULL
        #pvalue <- pvalueT[,j]
        if(Nchr == 1){
            bp <- ifelse((max_no_na(pvalue.pos.list[[1]]) -  min_no_na(pvalue.pos.list[[1]])) > 1000000, 1000000, 1000)
            bp_lab <- ifelse(bp == 1000000, " (Mb)", " (Kb)")
            pvalue.posN <- pvalue.pos.list[[1]] + band
            ticks <- seq(min_no_na(pvalue.pos.list[[1]]), max_no_na(pvalue.pos.list[[1]]), length=10)
            ticks <- seq(round(min_no_na(pvalue.pos.list[[1]]) / bp), round(max_no_na(pvalue.pos.list[[1]]) / bp), round((ticks[2]-ticks[1])/bp) + 0.5)
            if(!round(max_no_na(pvalue.pos.list[[1]]) / bp) %in% ticks){
                if(round(max_no_na(pvalue.pos.list[[1]]) / bp) - ticks[length(ticks)] > 0.5 * ticks[2])
                ticks <- c(ticks, round(max_no_na(pvalue.pos.list[[1]]) / bp))
            }
            ticks <- ticks[-1]
            chr.labels <- ticks
            ticks <- ticks * bp + band
            chr.border <- FALSE
        }else{
            for(i in 0:(Nchr-1)){
                if (i==0){
                    #pvalue <- append(pvalue,rep(Inf,band),after=0)
                    pvalue.posN <- pvalue.pos.list[[i+1]] + band
                    ticks[i+1] <- max_no_na(pvalue.posN)-floor(max_no_na(pvalue.pos.list[[i+1]])/2)
                    chr.border.pos[i+1] <- max_no_na(pvalue.posN) + 0.5 * band
                }else{
                    #pvalue <- append(pvalue,rep(Inf,band),after=sum(Num[1:i])+i*band)
                    pvalue.posN <- c(pvalue.posN, max_no_na(pvalue.posN) + band + pvalue.pos.list[[i+1]])
                    ticks[i+1] <- max_no_na(pvalue.posN)-floor(max_no_na(pvalue.pos.list[[i+1]])/2)
                    chr.border.pos[i+1] <- max_no_na(pvalue.posN) + 0.5 * band
                }
            }
            chr.border.pos=chr.border.pos[-length(chr.border.pos)]
        }

        if(!is.null(chr.labels) & Nchr != 1){
            chr.labels <- as.character(chr.labels)
            if(length(chr.labels) != Nchr)  stop("length of 'chr.labels' should equal to the number of chromosomes.")
            ticks.logi <- rep(TRUE, length(ticks))
            for(ti in 1:Nchr){
                if(is.na(chr.labels[ti]))    ticks.logi[ti] <- FALSE
            }
            if(!all(ticks.logi)){
                chr.labels <- chr.labels[ticks.logi]
                ticks <- ticks[ticks.logi]
            }
        }

        pvalue.posN.list <- tapply(pvalue.posN, Pmap[, 1], list)
        
        #merge the pvalues of traits by column
        if(LOG10){
            logpvalueT <- -log10(pvalueT)
        }else{
            logpvalueT <- pvalueT
        }

        add <- list()
        for(i in 1:R){
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            add[[i]] <- c(Num,rep(0,N[i]*length(colx)-Nchr))
        }

        circleMin <- (min_no_na(pvalue.posN) - band - 1)
        TotalN <- max_no_na(pvalue.posN)-circleMin
        
        if(length(chr.den.col) > 1){
            cir.density=TRUE
            den.fold <- 20
            density.list <- DensityPlot(Pmap[, 1], Pmap[, 2], chr.ori, chr.pos.max=FALSE, col=chr.den.col, plot=FALSE, bin=bin.size, bin.breaks=bin.breaks)
        }else{
            cir.density=FALSE
        }
    }

    #plot circle Manhattan
    if("c" %in% plot.type){

        signal.line.index <- NULL
        if(!is.null(threshold)){
            if(!is.null(signal.line)){
                for(l in 1:R){
                    if(!is.null(threshold[[l]])){
                        if(LOG10){
                            signal.line.index <- c(signal.line.index,which(pvalueT[,l] < min_no_na(threshold[[l]])))
                        }else{
                            signal.line.index <- c(signal.line.index,which(pvalueT[,l] > max_no_na(threshold[[l]])))
                        }
                    }
                }
                signal.line.index <- unique(signal.line.index)
            }
            signal.line.index <- pvalue.posN[signal.line.index]
        }

        if(file.output){
            ht=ifelse(is.null(height), 10, height)
            wh=ifelse(is.null(width), 10, width)
            if(file=="jpg") jpeg(paste("Cir_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
            if(file=="pdf") pdf(paste("Cir_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
            if(file=="tiff")    tiff(paste("Cir_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
            if(file=="png") png(paste("Cir_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
            par(pty="s", xpd=TRUE, mar=c(1,1,1,1))
        }
        if(!file.output){
            ht=ifelse(is.null(height), 10, height)
            wh=ifelse(is.null(width), 10, width)
            if(is.null(dev.list())) dev.new(width=wh, height=ht)
            par(pty="s", xpd=TRUE)
        }
        RR <- r+H*R+cir.band*R
        if(cir.density){
            plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.1*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
        }else{
            plot(NULL,xlim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),ylim=c(1.05*(-RR-4*cir.chr.h),1.05*(RR+4*cir.chr.h)),axes=FALSE,xlab="",ylab="")
        }
        if(!is.null(signal.line)){
            if(!is.null(signal.line.index)){
                X1chr <- (RR)*sin(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                Y1chr <- (RR)*cos(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                X2chr <- (r)*sin(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                Y2chr <- (r)*cos(2*base::pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                segments(X1chr,Y1chr,X2chr,Y2chr,lty=2,lwd=signal.line,col="grey")
            }
        }
        for(i in 1:R){
        
            #get the colors for each trait
            colx <- col[i,]
            colx <- colx[!is.na(colx)]

            if(verbose) cat(paste(" Circular Manhattan plotting ",trait[i],".\n",sep=""))
            pvalue <- pvalueT[,i]
            logpvalue <- logpvalueT[,i]
            if(is.null(ylim)){
                if(LOG10){
                    Max <- max_ylim(-log10(min_no_na(pvalue)))
                    Min <- min_ylim(-log10(max_no_na(pvalue)))
                }else{
                    Max <- max_ylim(max_no_na(pvalue))
                    #if(abs(Max)<=1)    Max <- max_no_na(pvalue)
                    Min <- min_ylim(min_no_na(pvalue))
                    #if(abs(Min)<=1)    Min <- min_no_na(pvalue)
                }
            }else{
                Max <- ylim[[i]][2]
                Min <- ylim[[i]][1]
            }
            Cpvalue <- (H*(logpvalue-Min))/(Max-Min)
            ylimIndx <- logpvalue >= Min & logpvalue <= Max
            if(outward==TRUE){
                if(cir.chr==TRUE & i == 1){
                    
                    #plot the boundary which represents the chromosomes
                    polygon.num <- 1000
                    for(k in 1:length(chr)){
                        if(k==1){
                            polygon.index <- seq(round(band/2)+1,-round(band/2)-circleMin+max_no_na(pvalue.posN.list[[1]]), length=polygon.num)
                            #change the axis from right angle into circle format
                            X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                            if(is.null(chr.den.col)){
                                polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])   
                            }else{
                                if(cir.density){
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                }else{
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                }
                            }
                        }else{
                            polygon.index <- seq(1+round(band/2)+max_no_na(pvalue.posN.list[[k-1]])-circleMin,-round(band/2)-circleMin+max_no_na(pvalue.posN.list[[k]]), length=polygon.num)
                            X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                            if(is.null(chr.den.col)){
                                polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])
                            }else{
                                if(cir.density){
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                }else{
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                }
                            }       
                        }
                    }
                    
                    if(cir.density){

                        if(file.output){
                            is_visable <- filter.points((RR+cir.chr.h)*sin(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), (RR+cir.chr.h)*cos(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), wh, ht, dpi=dpi)
                        }else{
                            is_visable <- rep(TRUE, length(pvalue.posN))
                        }
                        segments(
                            (RR)*sin(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR)*cos(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*sin(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*cos(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            col=density.list$den.col[is_visable], lwd=0.5
                        )
                        legend(
                            x=RR+4*cir.chr.h,
                            y=(RR+4*cir.chr.h)/2,
                            title="", legend=density.list$legend.y, pch=15, pt.cex=3, col=density.list$legend.col,
                            cex=legend.cex, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.3, xjust=0, xpd=TRUE
                        )
                        
                    }
                    
                    # XLine=(RR+cir.chr.h)*sin(2*base::pi*(1:TotalN)/TotalN)
                    # YLine=(RR+cir.chr.h)*cos(2*base::pi*(1:TotalN)/TotalN)
                    # lines(XLine,YLine,lwd=1.5)
                    if(cir.density){
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
                        circle.plot(myr=RR,lwd=1.5,add=TRUE,col='grey')
                    }else{
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
                        circle.plot(myr=RR,lwd=1.5,add=TRUE)
                    }

                }
                
                X=(Cpvalue[ylimIndx]+r+H*(i-1)+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                Y=(Cpvalue[ylimIndx]+r+H*(i-1)+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                if(file.output){
                    is_visable <- filter.points(X, Y, wh, ht, dpi=dpi)
                }else{
                    is_visable <- rep(TRUE, length(X))
                }

                if(cir.axis && cir.axis.grid){
                    circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                }

                points(X[is_visable],Y[is_visable],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]])[ylimIndx][is_visable])
                
                #plot the legend for each trait
                if(cir.axis==TRUE){
                    #try to get the number after radix point
                    if((Max-Min) > 1) {
                        round.n=2
                    }else{
                        if(Max == 1){
                            round.n=1
                        }else{
                            round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                        }
                    }
                    segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)

                    lab=seq(round(Min+(Max-Min)*0,round.n), round(Min+(Max-Min)*1,round.n), length=5)
                    text(-H/20,r+H*(i-0.94)+cir.band*(i-1),lab[1],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.75)+cir.band*(i-1),lab[2],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.5)+cir.band*(i-1),lab[3],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.25)+cir.band*(i-1),lab[4],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.06)+cir.band*(i-1),lab[5],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                }
                
                if(!is.null(threshold[[i]])){
                    if(sum(threshold[[i]]!=0)==length(threshold[[i]])){
                        for(thr in 1:length(threshold[[i]])){
                            significantline1=ifelse(LOG10, H*(-log10(threshold[[i]][thr])-Min)/(Max-Min), H*(threshold[[i]][thr]-Min)/(Max-Min))
                            #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*base::pi*(0:TotalN)/TotalN)
                            #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*base::pi*(0:TotalN)/TotalN)
                            if(significantline1<H){
                                #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                            }else{
                                warning(paste("No significant points for ",trait[i]," pass the threshold level using threshold=",threshold[[i]][thr],"!",sep=""))
                            }
                        }
                    }
                }
                
                if(!is.null(threshold[[i]])){
                    if(sum(threshold[[i]]!=0)==length(threshold[[i]])){
                        if(amplify==TRUE){
                            if(LOG10){
                                threshold[[i]] <- sort(threshold[[i]])
                                significantline1=H*(-log10(max_no_na(threshold[[i]]))-Min)/(Max-Min)
                            }else{
                                threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                significantline1=H*(min_no_na(threshold[[i]])-Min)/(Max-Min)
                            }
                            
                            p_amp.index <- which(Cpvalue>=significantline1)
                            HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            
                            #cover the points that exceed the threshold with the color "white"
                            points(HX1,HY1,pch=19,cex=cex[1],col="white")
                        
                            for(ll in 1:length(threshold[[i]])){
                                if(ll == 1){
                                    if(LOG10){
                                        significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                    }else{
                                        significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                    }
                                    p_amp.index <- which(Cpvalue>=significantline1)
                                    HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                    HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                }else{
                                    if(LOG10){
                                        significantline0=H*(-log10(threshold[[i]][ll-1])-Min)/(Max-Min)
                                        significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                    }else{
                                        significantline0=H*(threshold[[i]][ll-1]-Min)/(Max-Min)
                                        significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                    }
                                    p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                    HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                    HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                }
                            
                                if(is.null(signal.col)){
                                    points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
                                }else{
                                    points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=signal.col[ll])
                                }
                            }
                        }
                    }
                }

                if(!is.null(highlight)){
                    HX1=(Cpvalue[highlight_index[[i]]]+r+H*(i-1)+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[highlight_index[[i]]]-round(band/2)-circleMin)/TotalN)
                    HY1=(Cpvalue[highlight_index[[i]]]+r+H*(i-1)+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[highlight_index[[i]]]-round(band/2)-circleMin)/TotalN)
                    points(HX1,HY1[highlight_index[[i]]],pch=19,cex=cex[1],col="white")
                    if(is.null(highlight.col)){
                        points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]])
                    }else{
                        points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=highlight_col[[i]])
                    }
                }

                if(cir.chr==TRUE){
                    ticks1=(RR+1.5*cir.chr.h)*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=(RR+1.5*cir.chr.h)*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                            angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                            text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5, adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab, sep=""),srt=angle, adj=c(0.5, 0),font=lab.font,cex=lab.cex-0.5)
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5, adj=c(0.5, 0))
                            }
                        }
                    }
                }else{
                    ticks1=1.01*RR*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=1.01*RR*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    # ticks1=(0.9*r)*sin(2*base::pi*(ticks-round(band/2))/TotalN)
                    # ticks2=(0.9*r)*cos(2*base::pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                        angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                        text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab, sep=""),srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }
                    }
                }
            }
            if(outward==FALSE){
                if(cir.chr==TRUE & i == 1){
                    # XLine=(2*cir.band+RR+cir.chr.h)*sin(2*base::pi*(1:TotalN)/TotalN)
                    # YLine=(2*cir.band+RR+cir.chr.h)*cos(2*base::pi*(1:TotalN)/TotalN)
                    # lines(XLine,YLine,lwd=1.5)

                    polygon.num <- 1000
                    for(k in 1:length(chr)){
                        if(k==1){
                            polygon.index <- seq(round(band/2)+1,-round(band/2)-circleMin+max_no_na(pvalue.posN.list[[1]]), length=polygon.num)
                            X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                                if(is.null(chr.den.col)){
                                    polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])   
                                }else{
                                    if(cir.density){
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                    }else{
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                    }
                                }
                        }else{
                            polygon.index <- seq(1+round(band/2)+max_no_na(pvalue.posN.list[[k-1]]),-round(band/2)-circleMin+max_no_na(pvalue.posN.list[[k]]), length=polygon.num)
                            X1chr=(RR)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*base::pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*base::pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*base::pi*(polygon.index)/TotalN)
                            if(is.null(chr.den.col)){
                                polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=rep(colx,ceiling(length(chr)/length(colx)))[k],border=rep(colx,ceiling(length(chr)/length(colx)))[k])   
                            }else{
                                    if(cir.density){
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col="grey",border="grey")
                                    }else{
                                        polygon(c(rev(X1chr),X2chr),c(rev(Y1chr),Y2chr),col=chr.den.col,border=chr.den.col)
                                    }
                            }   
                        }
                    }
                    if(cir.density){

                        if(file.output){
                            is_visable <- filter.points((RR+cir.chr.h)*sin(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), (RR+cir.chr.h)*cos(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), wh, ht, dpi=dpi)
                        }else{
                            is_visable <- rep(TRUE, length(pvalue.posN))
                        }
                        segments(
                            (RR)*sin(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR)*cos(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*sin(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*cos(2*base::pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            col=density.list$den.col[is_visable], lwd=0.5
                        )
                        legend(
                            x=RR+4*cir.chr.h,
                            y=(RR+4*cir.chr.h)/2,
                            title="", legend=density.list$legend.y, pch=15, pt.cex=3, col=density.list$legend.col,
                            cex=legend.cex, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.3, xjust=0, xpd=TRUE
                        )
                        
                    }
                    
                    if(cir.density){
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
                        circle.plot(myr=RR,lwd=1.5,add=TRUE,col='grey')
                    }else{
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
                        circle.plot(myr=RR,lwd=1.5,add=TRUE)
                    }

                }

                X=(-Cpvalue[ylimIndx]+r+H*i+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                Y=(-Cpvalue[ylimIndx]+r+H*i+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                if(file.output){
                    is_visable <- filter.points(X, Y, wh, ht, dpi=dpi)
                }else{
                    is_visable <- rep(TRUE, length(X))
                }

                if(cir.axis && cir.axis.grid){
                    circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                }

                points(X[is_visable],Y[is_visable],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]])[ylimIndx][is_visable])
                
                if(cir.axis==TRUE){
                    
                    #try to get the number after radix point
                    if((Max-Min)<=1) {
                        if(Max == 1){
                            round.n=1
                        }else{
                            round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                        }
                    }else{
                        round.n=2
                    }
                    segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.axis.col,lwd=axis.lwd)
                    
                    lab=seq(round(Min+(Max-Min)*0,round.n), round(Min+(Max-Min)*1,round.n), length=5)
                    text(-H/20,r+H*(i-0.06)+cir.band*(i-1),lab[1],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.25)+cir.band*(i-1),lab[2],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.5)+cir.band*(i-1),lab[3],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.75)+cir.band*(i-1),lab[4],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                    text(-H/20,r+H*(i-0.94)+cir.band*(i-1),lab[5],adj=1,col=cir.axis.col,cex=axis.cex*0.5,font=lab.font)
                }
                
                if(!is.null(threshold[[i]])){
                    if(sum(threshold[[i]]!=0)==length(threshold[[i]])){
                    
                        for(thr in 1:length(threshold[[i]])){
                            significantline1=ifelse(LOG10, H*(-log10(threshold[[i]][thr])-Min)/(Max-Min), H*(threshold[[i]][thr]-Min)/(Max-Min))
                            #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
                            #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
                            if(significantline1<H){
                                #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                circle.plot(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                            }else{
                                warning(paste("No significant points for ",trait[i]," pass the threshold level using threshold=",threshold[[i]][thr],"!",sep=""))
                            }
                        }
                        if(amplify==TRUE){
                            if(LOG10){
                                threshold[[i]] <- sort(threshold[[i]])
                                significantline1=H*(-log10(max_no_na(threshold[[i]]))-Min)/(Max-Min)
                            }else{
                                threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                significantline1=H*(min_no_na(threshold[[i]])-Min)/(Max-Min)
                            }
                            p_amp.index <- which(Cpvalue>=significantline1)
                            HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            
                            #cover the points that exceed the threshold with the color "white"
                            points(HX1,HY1,pch=19,cex=cex[1],col="white")
                            
                                for(ll in 1:length(threshold[[i]])){
                                    if(ll == 1){
                                        if(LOG10){
                                            significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                        }else{
                                            significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1)
                                        HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                        HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                    }else{
                                        if(LOG10){
                                            significantline0=H*(-log10(threshold[[i]][ll-1])-Min)/(Max-Min)
                                            significantline1=H*(-log10(threshold[[i]][ll])-Min)/(Max-Min)
                                        }else{
                                            significantline0=H*(threshold[[i]][ll-1]-Min)/(Max-Min)
                                            significantline1=H*(threshold[[i]][ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                        HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                        HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                    
                                    }
                                
                                    if(is.null(signal.col)){
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rep(rep(colx,N[i]),add[[i]])[p_amp.index])
                                    }else{
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=signal.col[ll])
                                    }
                                }
                        }
                    }
                }
                
                if(!is.null(highlight)){
                    HX1=(-Cpvalue[highlight_index[[i]]]+r+H*i+cir.band*(i-1))*sin(2*base::pi*(pvalue.posN[highlight_index[[i]]]-round(band/2)-circleMin)/TotalN)
                    HY1=(-Cpvalue[highlight_index[[i]]]+r+H*i+cir.band*(i-1))*cos(2*base::pi*(pvalue.posN[highlight_index[[i]]]-round(band/2)-circleMin)/TotalN)
                    points(HX1,HY1,pch=19,cex=cex[1],col="white")
                    if(is.null(highlight.col)){
                        points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]])
                    }else{
                        points(HX1,HY1,pch=highlight.pch,cex=highlight.cex,col=highlight_col[[i]])
                    }
                }

                if(cir.chr==TRUE){
                    ticks1=(RR+1.5*cir.chr.h)*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=(RR+1.5*cir.chr.h)*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                          angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                          text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab,sep=""),srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }
                    }
                }else{
                    ticks1=1.01*RR*sin(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=1.01*RR*cos(2*base::pi*(ticks-round(band/2)-circleMin)/TotalN)
                    # ticks1=RR*sin(2*base::pi*(ticks-round(band/2))/TotalN)
                    # ticks2=RR*cos(2*base::pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                        
                            #adjust the angle of labels of circle plot
                            angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                            text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab,sep=""),srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=lab.font,cex=lab.cex-0.5,adj=c(0.5, 0))
                            }
                        }
                    }   
                }
            }
        }
        if(file.output) dev.off()
        #print("Circular-Manhattan has been finished!",quote=F)
    }

    if("m" %in% plot.type){

        is_visable <- list()
        for(i in 1:R){
            if(file.output){
                ht=ifelse(is.null(height), 6, height)
                wh=ifelse(is.null(width), 14, width)
                is_visable[[i]] <- filter.points(pvalue.posN, logpvalueT[,i], wh, ht, dpi=dpi)
            }else{
                is_visable[[i]] <- rep(TRUE, nrow(logpvalueT))
            }
        }

        if(multracks | multraits){
            if(R < 2)   stop("need more than one trait.")
            if(multracks){
                if(file.output){
                    ht=ifelse(is.null(height), 6, height)
                    wh=ifelse(is.null(width), 14, width)
                    if(file=="jpg") jpeg(paste("Multi-tracks_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi*R,res=dpi,quality=100)
                    if(file=="pdf") pdf(paste("Multi-tracks_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht*R)
                    if(file=="tiff")    tiff(paste("Multi-tracks_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi*R,res=dpi)
                    if(file=="png") png(paste("Multi-tracks_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi*R,res=dpi,bg=NA)
                    par(mfcol=c(R,1), xaxs="i")
                }
                if(!file.output){
                    ht=ifelse(is.null(height), 6, height)
                    wh=ifelse(is.null(width), 14, width)
                    if(is.null(dev.list())) dev.new(width=wh, height=ht)
                    # par(xpd=TRUE)
                }
                for(i in 1:R){
                    # Add room for x axis, if there are multiple
                    btwn_adj=if(multracks.xaxis) 2 else 0
                    if(i == 1)  par(mar=c(mar.between + btwn_adj, mar[2]+1, mar[3], 0))
                    if(i == R)  par(mar=c(mar[1]+1, mar[2]+1, 0, 0))
                    if(i != 1 & i != R) par(mar=c(mar.between + btwn_adj, mar[2]+1, 0, 0))
                    if(verbose) cat(paste(" Multi-tracks Manhattan plotting ",trait[i],".\n",sep=""))
                    colx=col[i,]
                    colx=colx[!is.na(colx)]
                    pvalue=pvalueT[,i]
                    logpvalue=logpvalueT[,i]
                    if(is.null(ylim)){
                        if(!is.null(threshold[[i]])){
                            # if(sum(threshold!=0)==length(threshold)){
                                if(LOG10){
                                    Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),-log10(min_no_na(threshold[[i]])))))
                                    Min <- min_ylim(min_no_na(c((-log10(max_no_na(pvalue))),-log10(max_no_na(threshold[[i]])))))
                                }else{
                                    Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold[[i]]))))
                                    #if(abs(Max)<=1)    Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
                                    Min<-min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold[[i]]))))
                                    #if(abs(Min)<=1)    Min=min_no_na(min_no_na(pvalue),min_no_na(threshold))
                                }
                        }else{
                            if(LOG10){
                                    Max=max_ylim((-log10(min_no_na(pvalue))))
                                    Min<-min_ylim((-log10(max_no_na(pvalue))))
                            }else{
                                    Max=max_ylim((max_no_na(pvalue)))
                                    #if(abs(Max)<=1)    Max=max_no_na(max_no_na(pvalue))
                                    Min=min_ylim((min_no_na(pvalue)))
                                    #if(abs(Min)<=1)    Min=min_no_na(min_no_na(pvalue))
                                    # }else{
                                        # Max=max_no_na(ceiling(max_no_na(pvalue)))
                                    # }
                            }
                        }
                        if((Max-Min)<=1){
                            plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=axis.cex*(R/2),font=lab.font,axes=FALSE,yaxs="r")
                        }else{
                            plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=axis.cex*(R/2),font=lab.font,axes=FALSE,yaxs="r")
                        }
                        mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex*(R/2), font=lab.font, xpd=TRUE)
                    }else{
                        Max <- max_no_na(ylim[[i]])
                        Min <- min_no_na(ylim[[i]])
                        plot(pvalue.posN[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=ylim[[i]],ann=FALSE,
                            cex.axis=axis.cex*(R/2),font=lab.font,axes=FALSE,yaxs="r")
                        mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex*(R/2), font=lab.font, xpd=TRUE)
                    }

                    if(chr.border){
                        for(b in 1:length(chr.border.pos)){
                            segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=axis.lwd, lty=2)
                        }
                    }

                    #add the names of traits on plot 
                    if(legend.pos=="left"){
                        text(min_no_na(pvalue.posN),Max,labels=trait[i],adj=c(-0.2, 1.2),font=4,cex=legend.cex*(R/2),xpd=TRUE) 
                    }else if(legend.pos=="middle"){
                        text((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))/2,Max,labels=trait[i],adj=c(0.5, 1.2),font=4,cex=legend.cex*(R/2),xpd=TRUE) 
                    }else{
                        text(max_no_na(pvalue.posN),Max,labels=trait[i],adj=c(1.2, 1.2),font=4,cex=legend.cex*(R/2),xpd=TRUE) 
                    }
                
                    if(i == R || multracks.xaxis){
                        if(chr.labels.angle == 0){
                            if(is.null(chr.labels)){
                                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd*(R/2),cex.axis=axis.cex*(R/2),font=lab.font,labels=c("Chr",chr.ori),padj=1)
                            }else{
                                if(Nchr == 1){
                                    axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd*(R/2), cex.axis=axis.cex*(R/2),font=lab.font,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                                }else{
                                    axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd*(R/2), cex.axis=axis.cex*(R/2),font=lab.font,labels=c("Chr",chr.labels))
                                }
                            }
                        }else{
                            axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd*(R/2),labels=FALSE)
                            if(is.null(chr.labels)){
                                text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex*(R/2), font=lab.font, labels=c("Chr",chr.ori), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                                # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd,cex.axis=axis.cex*(R/2),font=lab.font,labels=c("Chr",chr.ori),padj=1)
                            }else{
                                if(Nchr == 1){
                                    # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd*(R/2), cex.axis=axis.cex*(R/2),font=lab.font,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                                    text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex*(R/2), font=lab.font, labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                                }else{
                                    # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd*(R/2), cex.axis=axis.cex*(R/2),font=lab.font,labels=c("Chr",chr.labels))
                                    text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex*(R/2), font=lab.font, labels=c("Chr",chr.labels), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                                }
                            }
                        }
                        axis(1, mgp=c(3,xticks.pos,0), at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=axis.lwd*(R/2))
                    }
                    #if(i==1) mtext("Manhattan plot",side=3,padj=-1,font=lab.font,cex=xn)
                    if(is.null(ylim)){
                        if((Max-Min)>1){
                            axis(2, las=1,lwd=axis.lwd*(R/2),cex.axis=axis.cex*(R/2),font=lab.font)
                            axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=axis.lwd*(R/2))
                        }else{
                            axis(2,las=1,lwd=axis.lwd*(R/2),cex.axis=axis.cex*(R/2),font=lab.font)
                            axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=axis.lwd*(R/2))
                        }
                    }else{
                        axis(2, las=1,lwd=axis.lwd*(R/2),cex.axis=axis.cex*(R/2),font=lab.font)
                        axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=axis.lwd*(R/2))
                    }
                    if(!is.null(threshold[[i]])){
                        for(thr in 1:length(threshold[[i]])){
                            h <- ifelse(LOG10, -log10(threshold[[i]][thr]), threshold[[i]][thr])
                            segments(0, h, max_no_na(pvalue.posN), h, col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                        }
                        if(amplify==TRUE){
                            if(LOG10){
                                threshold[[i]] <- sort(threshold[[i]])
                                sgline1=-log10(max_no_na(threshold[[i]]))
                            }else{
                                threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                sgline1=min_no_na(threshold[[i]])
                            }
                            sgindex=which(logpvalue>=sgline1)
                            HY1=logpvalue[sgindex]
                            HX1=pvalue.posN[sgindex]
                            
                            #cover the points that exceed the threshold with the color "white"
                            points(HX1,HY1,pch=pch,cex=cex[2]*R,col="white")
                            
                            for(ll in 1:length(threshold[[i]])){
                                if(ll == 1){
                                    if(LOG10){
                                        sgline1=-log10(threshold[[i]][ll])
                                    }else{
                                        sgline1=threshold[[i]][ll]
                                    }
                                    sgindex=which(logpvalue>=sgline1)
                                    HY1=logpvalue[sgindex]
                                    HX1=pvalue.posN[sgindex]
                                }else{
                                    if(LOG10){
                                        sgline0=-log10(threshold[[i]][ll-1])
                                        sgline1=-log10(threshold[[i]][ll])
                                    }else{
                                        sgline0=threshold[[i]][ll-1]
                                        sgline1=threshold[[i]][ll]
                                    }
                                    sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
                                    HY1=logpvalue[sgindex]
                                    HX1=pvalue.posN[sgindex]
                                }

                                if(is.null(signal.col)){
                                    points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*R,col=rep(rep(colx,N[i]),add[[i]])[sgindex])
                                }else{
                                    points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll]*R,col=signal.col[ll])
                                }
                            }
                        }
                    }

                    if(!is.null(highlight)){
                        # points(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],pch=pch,cex=cex[2]*R,col="white")
                        if(!is.na(highlight_index[[i]][1])){
                            if(is.null(highlight.col)){
                                highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.cex*R,text.cex=highlight.text.cex*R/2, pch=highlight.pch,type=highlight.type,point.col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]],text.col=highlight.text.col,text.font=highlight.text.font)
                            }else{
                                highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.cex*R,text.cex=highlight.text.cex*R/2, pch=highlight.pch,type=highlight.type,point.col=highlight_col[[i]],text.col=highlight.text.col,text.font=highlight.text.font)
                            }
                        }
                    }
                    if(!is.null(main) & R == 1)  title(main=main[1], cex.main=main.cex, font.main= main.font)
                    if(box) box(lwd=axis.lwd)
                    #if(!is.null(threshold) & !is.null(signal.line))    abline(v=pvalue.posN[which(pvalueT[,i] < min_no_na(threshold))],col="grey",lty=2,lwd=signal.line)
                }
                if(file.output) dev.off()
            }
            if(multraits){
                if(file.output){
                    ht=ifelse(is.null(height), 6, height)
                    wh=ifelse(is.null(width), 14, width)
                    if(file=="jpg") jpeg(paste("Multi-traits_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                    if(file=="pdf") pdf(paste("Multi-traits_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
                    if(file=="tiff")    tiff(paste("Multi-traits_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                    if(file=="png") png(paste("Multi-traits_Manhtn.",ifelse(file.name=="",taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                    if(!is.null(legend.ncol) && legend.pos=="middle"){
                        mar[3] = mar[3] + ceiling(length(trait) / legend.ncol)
                    }
                    par(mar=mar,xaxs="i",yaxs="r")
                }
                if(!file.output){
                    ht=ifelse(is.null(height), 6, height)
                    wh=ifelse(is.null(width), 14, width)
                    if(is.null(dev.list())) dev.new(width=wh, height=ht)
                    # par(xpd=TRUE)
                }
                
                pvalue <- as.vector(Pmap[,3:(R+2)])
                if(is.null(ylim)){
                    if(!is.null(threshold)){
                        if(LOG10){
                            Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),-log10(min_no_na(unlist(threshold))))))
                            Min<-min_ylim(min_no_na(c((-log10(max_no_na(pvalue))),-log10(max_no_na(unlist(threshold))))))
                        }else{
                            Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(unlist(threshold)))))
                            # if(abs(Max)<=1)   Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
                            Min <- min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(unlist(threshold)))))
                            # if(abs(Min)<=1)   Min=min_no_na(c(min_no_na(pvalue),min_no_na(threshold)))
                        }
                    }else{
                        if(LOG10){
                                Max=max_ylim((-log10(min_no_na(pvalue))))
                                Min=min_ylim((-log10(max_no_na(pvalue))))
                        }else{
                                Max=max_ylim((max_no_na(pvalue)))
                                # if(abs(Max)<=1)   Max=max_no_na(max_no_na(pvalue))
                                Min<- min_ylim((min_no_na(pvalue)))
                                # if(abs(Min)<=1)   Min=min_no_na(min_no_na(pvalue))
                                # }else{
                                    # Max=max_no_na(ceiling(max_no_na(pvalue)))
                        }
                    }
                    if((Max-Min)<=1){
                        if(cir.density){
                            plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold, Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }else{
                            plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }
                    }else{
                        if(cir.density){
                            plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }else{
                            plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }
                    }
                    mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                }else{
                    Max <- max_no_na(unlist(ylim))
                    Min <- min_no_na(unlist(ylim))
                    if(cir.density){
                        plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-Max/den.fold,Max),ann=FALSE,
                            cex.axis=axis.cex,font=lab.font,axes=FALSE)
                    }else{
                        plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min, Max),ann=FALSE,
                            cex.axis=axis.cex,font=lab.font,axes=FALSE)
                    }
                    mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                }

                # Max1 <- Max
                # Min1 <- Min
                # if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
                # if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
                if(length(unique(col)) == 1 && is.null(signal.col)) stop("'signal.col' is NULL.")
                if(length(unique(col)) == 1 && amplify == FALSE)    stop("'amplify' is FALSE.")
                legend_col <- t(col)[1:R]
                if(length(unique(col)) == 1)    legend_col <- rep(signal.col, R)[1:R]
                if(legend.pos=="middle"){
                    if(is.null(legend.ncol)){
                        legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,Max,trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=TRUE,xjust=0.5,yjust=0,xpd=TRUE)
                    }else{
                        legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,Max,trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=FALSE,ncol=legend.ncol,xjust=0.5,yjust=0,xpd=TRUE)
                    }
                }else{
                    if(is.null(legend.ncol)){
                        legend(ifelse(legend.pos=="left","topleft","topright"),trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=FALSE,xpd=TRUE)
                    }else{
                        legend(ifelse(legend.pos=="left","topleft","topright"),trait,col=legend_col,pch=pch,text.font=6,cex=legend.cex,box.col=NA,horiz=FALSE,ncol=legend.ncol,xpd=TRUE)
                    }
                }

                if(chr.labels.angle == 0){
                    if(is.null(chr.labels)){
                        axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks),lwd=axis.lwd,cex.axis=axis.cex,font=lab.font,labels=c("Chr",chr.ori))
                    }else{
                        if(Nchr == 1){
                            axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=lab.font,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                        }else{
                            axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=lab.font,labels=c("Chr",chr.labels))
                        }
                    }
                }else{
                    axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd,labels=FALSE)
                    if(is.null(chr.labels)){
                        text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex, font=lab.font, labels=c("Chr",chr.ori), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                        # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks),lwd=axis.lwd,cex.axis=axis.cex,font=lab.font,labels=c("Chr",chr.ori)) 
                    }else{
                        if(Nchr == 1){
                            text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex, font=lab.font, labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                            # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=lab.font,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                        }else{
                            text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex, font=lab.font, labels=c("Chr",chr.labels), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                            # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=lab.font,labels=c("Chr",chr.labels))
                        }
                    }
                }
                axis(1, mgp=c(3,xticks.pos,0), at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=axis.lwd)
                if(is.null(ylim)){
                    if((Max-Min)>1){
                        #print(seq(0,(Max+1),ceiling((Max+1)/10)))
                        axis(2,las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                        axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                        legend.y <- Max
                    }else{
                        axis(2,las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                        axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                        legend.y <- Max
                    }
                }else{
                    axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                    axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                    legend.y <- Max
                }
                if(chr.border){
                    for(b in 1:length(chr.border.pos)){
                        segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=axis.lwd, lty=2)
                    }
                }

                if(length(unique(col)) != 1){
                    sam.index <- list()
                    trait_max_n <- 0
                    trait_max <- 0
                    for(l in 1:R){
                        sam.index[[l]] <- c(1:nrow(Pmap))[is_visable[[l]] & !is.na(logpvalueT[,l])]
                        if(length(sam.index[[l]]) >= trait_max_n){
                            trait_max_n=length(sam.index[[l]])
                            trait_max=l
                        }
                    }
                    
                    #change the sample number according to Pmap
                    #sam.num <- ceiling(nrow(Pmap)/100)
                    sam.num <- 1000
                    cat_bar <- seq(1, 100, 1)
                    trait_n <- sapply(sam.index, length)
                    trait_sams <- ceiling(trait_n / sam.num)
                    trait_max_sams <- max(trait_sams)
                    trait_1st_sam <- trait_max_sams - trait_sams + 1
                    trait_full_sams <- floor(trait_n / sam.num)
                    trait_1st_full_sam <- trait_max_sams - trait_full_sams + 1
                    for(sam in 1:trait_max_sams) {
                        for(i in 1:R){
                            if(sam < trait_1st_sam[i]){
                                # nothing
                            }else{
                                if(sam < trait_1st_full_sam[i]){
                                    plot.index <- sample(sam.index[[i]], trait_n[i] %% sam.num, replace=FALSE)
                                }else{
                                    plot.index <- sample(sam.index[[i]], sam.num, replace=FALSE)
                                }
                                sam.index[[i]] <- sam.index[[i]][-which(sam.index[[i]] %in% plot.index)]
                                logpvalue=logpvalueT[plot.index,i]
                                if(!is.null(ylim)){indexx <- logpvalue>=min_no_na(ylim[[i]])}else{indexx <- 1:length(logpvalue)}
                                points(pvalue.posN[plot.index][indexx],logpvalue[indexx],pch=pch[i],type=type,lwd=cex[2]+1,cex=cex[2],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                            }
                        }
                        if(verbose){
                            progress <- round((nrow(Pmap) - length(sam.index[[trait_max]])) * 100 / nrow(Pmap))
                            if(progress %in% cat_bar){
                                cat(" Multi-traits Rectangular plotting ... (finished ", progress, "%)\r", sep="")
                                cat_bar <- cat_bar[cat_bar != progress]
                                if(progress == 100) cat("\n")
                            }
                        }
                    }
                }else{
                    for(i in 1:R){
                        logpvalue=logpvalueT[,i]
                        if(!is.null(ylim)){indexx <- logpvalue>=min_no_na(ylim[[i]])}else{indexx <- 1:length(logpvalue)}
                        points(pvalue.posN[indexx],logpvalue[indexx],pch=pch[i],type=type,lwd=cex[2]+1,cex=cex[2],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                    }
                }
                if(!is.null(threshold)){
                    for(thr in 1:length(threshold[[i]])){
                        h <- ifelse(LOG10, -log10(threshold[[i]][thr]), threshold[[i]][thr])
                        segments(0, h, max_no_na(pvalue.posN), h, col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                    }
                    if(amplify==TRUE){
                        if(length(unique(col)) != 1){
                            for(i in 1:R){
                                logpvalue=logpvalueT[, i]
                                for(ll in 1:length(threshold[[i]])){
                                    if(ll == 1){
                                        if(LOG10){
                                            sgline1=-log10(threshold[[i]][ll])
                                        }else{
                                            sgline1=threshold[[i]][ll]
                                        }
                                        sgindex=which(logpvalue>=sgline1)
                                        HY1=logpvalue[sgindex]
                                        HX1=pvalue.posN[sgindex]
                                    }else{
                                        if(LOG10){
                                            sgline0=-log10(threshold[[i]][ll-1])
                                            sgline1=-log10(threshold[[i]][ll])
                                        }else{
                                            sgline0=threshold[[i]][ll-1]
                                            sgline1=threshold[[i]][ll]
                                        }
                                        sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
                                        HY1=logpvalue[sgindex]
                                        HX1=pvalue.posN[sgindex]
                                    }
                                    points(HX1,HY1,pch=pch[i],cex=cex[2],col="white")
                                    if(is.null(signal.col)){
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                                    }else{
                                        points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rgb(t(col2rgb(signal.col[ll])), alpha=points.alpha, maxColorValue=255))
                                    }
                                    
                                }
                            }
                        }else{
                            for(i in 1:R){
                                logpvalue=logpvalueT[, i]
                                if(LOG10){
                                    sgindex = which(logpvalue > -log10(min(unlist(threshold))))
                                }else{
                                    sgindex = which(logpvalue > max(unlist(threshold)))
                                }
                                HY1=logpvalue[sgindex]
                                HX1=pvalue.posN[sgindex]
                                points(HX1,HY1,pch=pch[i],cex=cex[2],col="white")
                                points(HX1,HY1,pch=rep(signal.pch, R)[i],cex=rep(signal.cex, R)[i],col=rgb(t(col2rgb(rep(signal.col, R)[i])), alpha=points.alpha, maxColorValue=255))
                            }
                        }
                    }
                }

                if(is.null(ylim)){ymin <- Min}else{ymin <- min_no_na(unlist(ylim))}
                if(cir.density){
                    for(yll in 1:length(pvalue.posN.list)){
                        polygon(c(min_no_na(pvalue.posN.list[[yll]]), min_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]])), 
                            c(ymin-0.5*(Max-Min)/den.fold, ymin-1.5*(Max-Min)/den.fold, 
                            ymin-1.5*(Max-Min)/den.fold, ymin-0.5*(Max-Min)/den.fold), 
                            col="grey", border="grey")
                    }
                    is_visable_den <- filter.points(pvalue.posN, ymin-0.5*(Max-Min)/den.fold, wh, ht, dpi=dpi)
                    segments(
                        pvalue.posN[is_visable_den],
                        ymin-0.5*(Max-Min)/den.fold,
                        pvalue.posN[is_visable_den],
                        ymin-1.5*(Max-Min)/den.fold,
                        col=density.list$den.col[is_visable_den], lwd=0.5
                    )
                    legend(
                        x=max_no_na(pvalue.posN)+band,
                        y=legend.y,
                        title="", legend=density.list$legend.y, pch=15, pt.cex=2.5, col=density.list$legend.col,
                        cex=legend.cex*0.8, bty="n",
                        y.intersp=1,
                        x.intersp=1,
                        yjust=0.9, xjust=0, xpd=TRUE
                    )          
                }
                if(!is.null(main))  title(main=main[1], cex.main=main.cex, font.main= main.font)
                if(box) box(lwd=axis.lwd)
                if(file.output) dev.off()
            }
        }else{
            #print("Starting Rectangular-Manhattan plot!",quote=F)
            if(file.name != "" && length(file.name) != R)   stop(paste("please provide a vector containing file names of all", R, "traits."))
            for(i in 1:R){
                colx=col[i,]
                colx=colx[!is.na(colx)]
                if(verbose) cat(paste(" Rectangular Manhattan plotting ",trait[i],".\n",sep=""))
                    if(file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(file=="jpg") jpeg(paste("Rect_Manhtn.",ifelse(file.name=="",trait[i],file.name[i]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                        if(file=="pdf") pdf(paste("Rect_Manhtn.",ifelse(file.name=="",trait[i],file.name[i]),".pdf",sep=""), width=wh,height=ht)
                        if(file=="tiff")    tiff(paste("Rect_Manhtn.",ifelse(file.name=="",trait[i],file.name[i]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                        if(file=="png") png(paste("Rect_Manhtn.",ifelse(file.name=="",trait[i],file.name[i]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                        par(mar=mar,xaxs="i",yaxs="r")
                    }
                    if(!file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(is.null(dev.list())) dev.new(width=wh, height=ht)
                        # par(xpd=TRUE)
                    }
                    
                    pvalue=pvalueT[,i]
                    logpvalue=logpvalueT[,i]
                    if(is.null(ylim)){
                        if(!is.null(threshold[[i]])){
                            if(sum(threshold[[i]]!=0)==length(threshold[[i]])){
                                if(LOG10 == TRUE){
                                    Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),(-log10(min_no_na(threshold[[i]]))))))
                                    Min <- min_ylim(min_no_na(c(-log10((max_no_na(pvalue))),-log10(max_no_na(threshold[[i]])))))
                                }else{
                                    Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold[[i]]))))
                                    #if(abs(Max)<=1)    Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
                                    Min <- min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold[[i]]))))
                                    #if(abs(Min)<=1)    Min=min_no_na(c(min_no_na(pvalue),min_no_na(threshold)))
                                }
                            }else{
                                if(LOG10){
                                    Max=max_ylim(-log10(min_no_na(pvalue)))
                                    Min<-min_ylim(-log10(max_no_na(pvalue)))
                                }else{
                                    Max=max_ylim(max_no_na(pvalue))
                                    #if(abs(Max)<=1)    Max=max_no_na(c(max_no_na(pvalue)))
                                    Min<-min_ylim(min_no_na(pvalue))
                                    #if(abs(Min)<=1)    Min=min_no_na(pvalue)
                                    # }else{
                                        # Max=max_no_na(ceiling(max_no_na(pvalue)))
                                    # }
                                }
                            }
                        }else{
                            if(LOG10){
                                    Max=max_ylim(-log10(min_no_na(pvalue)))
                                    Min<-min_ylim(-log10(max_no_na(pvalue)))
                            }else{
                                    Max=max_ylim(max_no_na(pvalue))
                                    #if(abs(Max)<=1)    Max=max_no_na(c(max_no_na(pvalue)))
                                    Min<-min_ylim(min_no_na(pvalue))
                                    #if(abs(Min)<=1)    Min=min_no_na(pvalue)
                                    # }else{
                                        # Max=max_no_na(ceiling(max_no_na(pvalue)))
                                    # }
                            }
                        }
                        if((Max-Min)<=1){
                            if(cir.density){
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold, Max),ann=FALSE,
                                    cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }else{
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }
                        }else{
                            if(cir.density){
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }else{
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=axis.cex,font=lab.font,axes=FALSE)
                            }
                        }
                        mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                    }else{
                        Max <- max_no_na(ylim[[i]])
                        Min <- min_no_na(ylim[[i]])
                        if(cir.density){
                            plot(pvalue.posN[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(min_no_na(ylim[[i]])-(Max-Min)/den.fold, max_no_na(ylim[[i]])),ann=FALSE,
                            cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }else{
                            plot(pvalue.posN[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=ylim[[i]],ann=FALSE,
                            cex.axis=axis.cex,font=lab.font,axes=FALSE)
                        }
                        mtext(side=2, text=ylab, line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                    }
                    # Max1 <- Max
                    # Min1 <- Min
                    # if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
                    # if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
                    if(chr.border){
                        for(b in 1:length(chr.border.pos)){
                            segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=axis.lwd, lty=2)
                        }
                    }

                    if(chr.labels.angle == 0){
                        if(!is.null(chr.labels)){
                            if(Nchr == 1){
                                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=lab.font,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                            }else{
                                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=lab.font,labels=c("Chr",chr.labels))
                                #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=axis.lwd)
                            }
                        }else{
                            axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=lab.font,labels=c("Chr",chr.ori))
                            #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=axis.lwd)
                        }
                    }else{
                        axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd,labels=FALSE)
                        if(!is.null(chr.labels)){
                            if(Nchr == 1){
                                text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex, font=lab.font, labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                            }else{
                                # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=2,labels=)
                                #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=axis.lwd)
                                text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex, font=lab.font, labels=c("Chr",chr.labels), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                            }
                        }else{
                            # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=axis.lwd, cex.axis=axis.cex,font=2,labels=c("Chr",chr.ori))
                            #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=axis.lwd)
                            text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=axis.cex, font=lab.font, labels=c("Chr",chr.ori), srt=chr.labels.angle, xpd=TRUE,adj=c(ifelse(chr.labels.angle < 0, 0, ifelse(chr.labels.angle == 0, 0.5, 1)), ifelse(chr.labels.angle == 0, 0.5, ifelse(abs(chr.labels.angle) > 45, 0.5, 1))))
                        }
                    }
                    axis(1, mgp=c(3,xticks.pos,0), at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=axis.lwd)
                    if(is.null(ylim)){
                        if((Max-Min)>1){
                            axis(2, las=1, lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                            axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                            legend.y <- Max
                        }else{
                            axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                            axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                            legend.y <- Max
                        }
                    }else{
                        axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex,font=lab.font)
                        axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=axis.lwd)
                        legend.y <- tail(ylim[[i]][2], 1)
                    }
                    if(!is.null(threshold[[i]])){
                        for(thr in 1:length(threshold[[i]])){
                            h <- ifelse(LOG10, -log10(threshold[[i]][thr]), threshold[[i]][thr])
                            # print(h)
                            # print(threshold.col[thr])
                            # print(threshold.lty[thr])
                            # print(threshold.lwd[thr])
                            segments(0, h, max_no_na(pvalue.posN), h,col=threshold.col[thr],lty=threshold.lty[thr],lwd=threshold.lwd[thr])
                        }
                        if(amplify == TRUE){
                            if(LOG10){
                                threshold[[i]] <- sort(threshold[[i]])
                                sgline1=-log10(max_no_na(threshold[[i]]))
                            }else{
                                threshold[[i]] <- sort(threshold[[i]], decreasing=TRUE)
                                sgline1=min_no_na(threshold[[i]])
                            }

                            sgindex=which(logpvalue>=sgline1)
                            HY1=logpvalue[sgindex]
                            HX1=pvalue.posN[sgindex]
                            
                            #cover the points that exceed the threshold with the color "white"
                            points(HX1,HY1,pch=pch,cex=cex[2],col="white")
                            
                            for(ll in 1:length(threshold[[i]])){
                                if(ll == 1){
                                    if(LOG10){
                                        sgline1=-log10(threshold[[i]][ll])
                                    }else{
                                        sgline1=threshold[[i]][ll]
                                    }
                                    sgindex=which(logpvalue>=sgline1)
                                    HY1=logpvalue[sgindex]
                                    HX1=pvalue.posN[sgindex]
                                }else{
                                    if(LOG10){
                                        sgline0=-log10(threshold[[i]][ll-1])
                                        sgline1=-log10(threshold[[i]][ll])
                                    }else{
                                        sgline0=threshold[[i]][ll-1]
                                        sgline1=threshold[[i]][ll]
                                    }
                                    sgindex=which(logpvalue>=sgline1 & logpvalue < sgline0)
                                    HY1=logpvalue[sgindex]
                                    HX1=pvalue.posN[sgindex]
                                }

                                if(is.null(signal.col)){
                                    points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=rep(rep(colx,N[i]),add[[i]])[sgindex])
                                }else{
                                    points(HX1,HY1,pch=signal.pch[ll],cex=signal.cex[ll],col=signal.col[ll])
                                }
                                
                            }
                        }

                    }

                    if(!is.null(highlight)){
                        # points(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],pch=pch,cex=cex[2],col="white")
                        if(!is.na(highlight_index[[i]][1])){
                            if(is.null(highlight.col)){
                                highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.cex,text.cex=highlight.text.cex, pch=highlight.pch,type=highlight.type,point.col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]],text.col=highlight.text.col,text.font=highlight.text.font)
                            }else{
                                highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),words=highlight.text[[i]],point.cex=highlight.cex,text.cex=highlight.text.cex, pch=highlight.pch,type=highlight.type,point.col=highlight_col[[i]],text.col=highlight.text.col,text.font=highlight.text.font)
                            }
                        }
                    }

                    #if(!is.null(threshold) & !is.null(signal.line))    abline(v=pvalue.posN[which(pvalueT[,i] < min_no_na(threshold))],col="grey",lty=2,lwd=signal.line)
            
                    if(is.null(ylim)){ymin <- Min}else{ymin <- min_no_na(ylim[[i]])}
                    if(cir.density){
                        for(yll in 1:length(pvalue.posN.list)){
                            polygon(c(min_no_na(pvalue.posN.list[[yll]]), min_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]])), 
                                c(ymin-0.5*(Max-Min)/den.fold, ymin-1.5*(Max-Min)/den.fold, 
                                ymin-1.5*(Max-Min)/den.fold, ymin-0.5*(Max-Min)/den.fold), 
                                col="grey", border="grey", xpd=TRUE)
                        }
                        is_visable_den <- filter.points(pvalue.posN, ymin-0.5*(Max-Min)/den.fold, wh, ht, dpi=dpi)
                        segments(
                            pvalue.posN[is_visable_den],
                            ymin-0.5*(Max-Min)/den.fold,
                            pvalue.posN[is_visable_den],
                            ymin-1.5*(Max-Min)/den.fold,
                            col=density.list$den.col[is_visable_den], lwd=0.5,xpd=TRUE
                        )
                        legend(
                            x=max_no_na(pvalue.posN)+band,
                            y=legend.y,
                            title="", legend=density.list$legend.y, pch=15, pt.cex=2.5, col=density.list$legend.col,
                            cex=legend.cex*0.8, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.9, xjust=0, xpd=TRUE
                        )
                        
                    }
                if(!is.null(main))  title(main=main[i], cex.main=main.cex, font.main= main.font)
                if(box) box(lwd=axis.lwd)
                if(file.output)  dev.off()
            }
        }
    }
        
    if("q" %in% plot.type){

        signal.col <- rep(signal.col,R)
        signal.pch <- rep(signal.pch,R)
        signal.cex <- rep(signal.cex*1.1,R)

        if(multracks | multraits){
            if(R < 2)   stop("need more than one trait.")
            if(multracks){
                if(file.output){
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 3.5, width)
                    if(file=="jpg") jpeg(paste("Multi-tracks_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".jpg",sep=""), width=R*wh*dpi,height=ht*dpi,res=dpi,quality=100)
                    if(file=="pdf") pdf(paste("Multi-tracks_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".pdf",sep=""), width=R*wh,height=ht)
                    if(file=="tiff")    tiff(paste("Multi-tracks_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".tiff",sep=""), width=R*wh*dpi,height=ht*dpi,res=dpi)
                    if(file=="png") png(paste("Multi-tracks_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".png",sep=""), width=R*wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                    par(mfcol=c(1,R),xpd=TRUE)
                }else{
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 3.5, width)
                    if(is.null(dev.list())) dev.new(width=wh*R, height=ht)
                    par(xpd=TRUE)
                }
                for(i in 1:R){
                    if(i == 1)  par(mar=c(mar[2], mar[2], mar[3], 0))
                    if(i == R)  par(mar=c(mar[2], 1.5, mar[3], mar[4]))
                    if(i != 1 & i != R) par(mar=c(mar[2], 1.5, mar[3], 0))
                    if(verbose) cat(paste(" Multi-tracks Q-Q plotting ",trait[i],".\n",sep=""))        
                    P.values=as.numeric(Pmap[,i+2])
                    P.values=P.values[!is.na(P.values)]
                    if(LOG10){
                        P.values=P.values[P.values>0]
                        P.values=P.values[P.values<1]
                        N=length(P.values)
                        P.values=P.values[order(P.values)]
                    }else{
                        N=length(P.values)
                        P.values=P.values[order(P.values,decreasing=TRUE)]
                    }
                    p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                    log.Quantiles <- -log10(p_value_quantiles)
                    if(LOG10){
                        log.P.values <- -log10(P.values)
                    }else{
                        log.P.values <- P.values
                    }
                    
                    #calculate the confidence interval of QQ-plot
                    if(conf.int){
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
                    }else{
                        c05 <- 1
                        c95 <- 1
                    }
                    
                    YlimMax <- max_no_na(c(floor(max_no_na(c(max_no_na(-log10(c05)), max_no_na(-log10(c95))))+1), floor(max_no_na(log.P.values)+1)))
                    if(is.null(ylim)){
                        plot(NULL, xlim=c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,YlimMax),xlab ="", ylab="")
                    }else{
                        plot(NULL, xlim=c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,max(ylim[[i]])),xlab ="", ylab="")
                    }
                    axis(1, mgp=c(3,xticks.pos,0), at=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), lwd=axis.lwd,labels=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), cex.axis=axis.cex)
                    axis(2, las=1, lwd=axis.lwd,cex.axis=axis.cex)
                    axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(ylim[[i]]))), labels=c("",""), tcl=0, lwd=axis.lwd)
                    
                    #plot the confidence interval of QQ-plot
                    if(conf.int){
                        if(is.null(conf.int.col)){
                            polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                        }else{
                            polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255))
                        }
                    }
                    if(!is.null(threshold.col)){par(xpd=FALSE); abline(a=0, b=1,lwd=threshold.lty[1], lty=threshold.lty[1], col=threshold.col[1]); par(xpd=TRUE)}
                    is_visable <- filter.points(log.Quantiles, log.P.values, wh, ht, dpi=dpi)
                    if(!is.null(threshold[[i]])){
                        # if(sum(threshold!=0)==length(threshold)){
                            thre.line=-log10(min_no_na(threshold[[i]]))
                            if(amplify==TRUE){
                                thre.index <- log.P.values<thre.line
                                if(sum(!thre.index)!=0){
                                    points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col=t(col)[i],pch=19,cex=cex[3])
                                
                                    #cover the points that exceed the threshold with the color "white"
                                    # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                    if(is.null(signal.col)){
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=t(col)[i],pch=signal.pch[i],cex=signal.cex[i])
                                    }else{
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=signal.col[i],pch=signal.pch[i],cex=signal.cex[i])
                                    }
                                }else{
                                    points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                                }
                            }else{
                                points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                            }
                        # }
                    }else{
                        points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                    }
                    mtext(side=1, text=expression(Expected~~-log[10](italic(p))), line=ylab.pos+2, cex=lab.cex, font=lab.font, xpd=TRUE)
                    if(i == 1)  mtext(side=2, text=expression(Observed~~-log[10](italic(p))), line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                    if(!is.null(main)) {
                        title(main=main[i], cex.main=main.cex, font.main= main.font)
                    }else{
                        title(main=trait[i], cex.main=main.cex, font.main= main.font) 
                    }
                    if(box) box(lwd=axis.lwd)
                }
                if(file.output) dev.off()
            }
            if(multraits){
                signal.col <- NULL
                log.Quantiles.max_no_na <- NULL
                for(i in 1:R){
                    P.values=as.numeric(Pmap[,i+2])
                    P.values=P.values[!is.na(P.values)]
                    p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                    log.Quantiles <- -log10(p_value_quantiles)
                    log.Quantiles.max_no_na <- c(log.Quantiles.max_no_na, max_no_na(log.Quantiles))
                }
                if(file.output){
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    if(file=="jpg") jpeg(paste("Multi-traits_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                    if(file=="pdf") pdf(paste("Multi-traits_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".pdf",sep=""), width=wh,height=ht)
                    if(file=="tiff")    tiff(paste("Multi-traits_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                    if(file=="png") png(paste("Multi-traits_QQplot.",ifelse(file.name=="",taxa,file.name[1]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                    par(mar=c(mar[2],mar[2],mar[3],mar[4]),xpd=TRUE)
                }else{  
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    dev.new(width=wh, height=ht)
                    par(xpd=TRUE)
                }
                p_value_quantiles=(1:nrow(Pmap))/(nrow(Pmap)+1)
                log.Quantiles <- -log10(p_value_quantiles)
                                            
                # calculate the confidence interval of QQ-plot
                if(conf.int){
                    N1=length(log.Quantiles)
                    c95 <- rep(NA,N1)
                    c05 <- rep(NA,N1)
                    for(j in 1:N1){
                        xi=ceiling((10^-log.Quantiles[j])*N1)
                        if(xi==0)xi=1
                        c95[j] <- qbeta(0.95,xi,N1-xi+1)
                        c05[j] <- qbeta(0.05,xi,N1-xi+1)
                    }
                    index=length(c95):1
                }
                
                if(!conf.int){c05 <- 1; c95 <- 1}
                
                if(is.null(ylim)){
                    Pmap.min_no_na <- Pmap[,3:(R+2)]
                    YlimMax <- max_no_na(c(floor(max_no_na(c(max_no_na(-log10(c05)), max_no_na(-log10(c95))))+1), -log10(min_no_na(Pmap.min_no_na[Pmap.min_no_na > 0]))))
                    plot(NULL, xlim=c(0,floor(max_no_na(log.Quantiles.max_no_na)+1)), axes=FALSE, xlab="", ylab="", cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0, floor(YlimMax+1)), main = "QQplot", cex.main=main.cex, font.main=main.font)
                }else{
                    plot(NULL, xlim=c(0,floor(max_no_na(log.Quantiles.max_no_na)+1)), axes=FALSE, xlab="", ylab="", cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0, max(unlist(ylim))),main = "QQplot", cex.main=main.cex, font.main=main.font)
                }
                legend("topleft",trait,col=rgb(t(col2rgb(t(col)[1:R])), alpha=points.alpha, maxColorValue=255),pch=19,cex=legend.cex,text.font=6,box.col=NA, xpd=TRUE)
                axis(1, mgp=c(3,xticks.pos,0), at=seq(0,floor(max_no_na(log.Quantiles.max_no_na)+1),ceiling((max_no_na(log.Quantiles.max_no_na)+1)/10)), lwd=axis.lwd,labels=seq(0,floor(max_no_na(log.Quantiles.max_no_na)+1),ceiling((max_no_na(log.Quantiles.max_no_na)+1)/10)), cex.axis=axis.cex)
                axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex)
                axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(unlist(ylim)))), labels=c("",""), tcl=0, lwd=axis.lwd)

                mtext(side=1, text=expression(Expected~~-log[10](italic(p))), line=ylab.pos+1, cex=lab.cex, font=lab.font, xpd=TRUE)
                mtext(side=2, text=expression(Observed~~-log[10](italic(p))), line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                
                for(i in 1:R){
                    if(verbose) cat(paste(" Multi-traits Q-Q plotting ",trait[i],".\n",sep=""))
                    P.values=as.numeric(Pmap[,i+2])
                    P.values=P.values[!is.na(P.values)]
                    if(LOG10){
                        P.values=P.values[P.values>0]
                        P.values=P.values[P.values<1]
                        N=length(P.values)
                        P.values=P.values[order(P.values)]
                    }else{
                        N=length(P.values)
                        P.values=P.values[order(P.values,decreasing=TRUE)]
                    }
                    p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                    log.Quantiles <- -log10(p_value_quantiles)
                    if(LOG10){
                        log.P.values <- -log10(P.values)
                    }else{
                        log.P.values <- P.values
                    }

                    #calculate the confidence interval of QQ-plot
                    if(conf.int){
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
                    }else{
                        c05 <- 1
                        c95 <- 1
                    }
    
                    # plot the confidence interval of QQ-plot
                    if(conf.int){
                        if(is.null(conf.int.col)){
                            polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                        }else{
                            polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255))
                        }
                    }
                       
                    if((i == R) & !is.null(threshold.col)){par(xpd=FALSE); abline(a=0, b=1,lwd=threshold.lty[1], lty=threshold.lty[1], col=threshold.col[1]); par(xpd=TRUE)}
                    # points(log.Quantiles, log.P.values, col=t(col)[i],pch=19,cex=cex[3])
                    is_visable <- filter.points(log.Quantiles, log.P.values, wh, ht, dpi=dpi)
                    if(!is.null(threshold[[i]])){
                        # if(sum(threshold!=0)==length(threshold)){
                            thre.line=-log10(min_no_na(threshold[[i]]))
                            if(amplify==TRUE){
                                thre.index <- log.P.values<thre.line
                                if(sum(!thre.index)!=0){
                                    points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                            
                                    # cover the points that exceed the threshold with the color "white"
                                    # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                    if(is.null(signal.col)){
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=signal.pch[i],cex=signal.cex[i])
                                    }else{
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=rgb(t(col2rgb(signal.col[i])), alpha=points.alpha, maxColorValue=255),pch=signal.pch[i],cex=signal.cex[i])
                                    }
                                }else{
                                    points(log.Quantiles[is_visable], log.P.values[is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                                }
                            }else{
                                points(log.Quantiles[is_visable], log.P.values[is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                            }
                        # }
                    }else{
                        points(log.Quantiles[is_visable], log.P.values[is_visable], col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),pch=19,cex=cex[3])
                    }
                }
                if(!is.null(main)) {
                    title(main=main[1], cex.main=main.cex, font.main= main.font)
                }
                if(box) box(lwd=axis.lwd)
                if(file.output) dev.off()
            }
        }else{
            if(file.name != "" && length(file.name) != R)   stop(paste("please provide a vector containing file names of all", R, "traits."))
            for(i in 1:R){
                if(verbose) cat(paste(" Q-Q plotting ",trait[i],".\n",sep=""))
                if(file.output){
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    if(file=="jpg") jpeg(paste("QQplot.",ifelse(file.name=="",trait[i],file.name[i]),".jpg",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,quality=100)
                    if(file=="pdf") pdf(paste("QQplot.",ifelse(file.name=="",trait[i],file.name[i]),".pdf",sep=""), width=wh,height=ht)
                    if(file=="tiff") tiff(paste("QQplot.",ifelse(file.name=="",trait[i],file.name[i]),".tiff",sep=""), width=wh*dpi,height=ht*dpi,res=dpi)
                    if(file=="png") png(paste("QQplot.",ifelse(file.name=="",trait[i],file.name[i]),".png",sep=""), width=wh*dpi,height=ht*dpi,res=dpi,bg=NA)
                    par(mar=c(mar[2],mar[2],mar[3],mar[4]),xpd=TRUE)
                }else{
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    if(is.null(dev.list())) dev.new(width=wh, height=ht)
                    par(xpd=TRUE)
                }
                P.values=as.numeric(Pmap[,i+2])
                P.values=P.values[!is.na(P.values)]
                if(LOG10){
                    P.values=P.values[P.values>0]
                    P.values=P.values[P.values<1]
                    N=length(P.values)
                    P.values=P.values[order(P.values)]
                }else{
                    N=length(P.values)
                    P.values=P.values[order(P.values,decreasing=TRUE)]
                }
                p_value_quantiles=(1:length(P.values))/(length(P.values)+1)
                log.Quantiles <- -log10(p_value_quantiles)
                if(LOG10){
                    log.P.values <- -log10(P.values)
                }else{
                    log.P.values <- P.values
                }
                
                #calculate the confidence interval of QQ-plot
                if(conf.int){
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
                }else{
                    c05 <- 1
                    c95 <- 1
                }
                if(is.null(ylim)){
                    YlimMax <- max_no_na(c(floor(max_no_na(c(max_no_na(-log10(c05)), max_no_na(-log10(c95))))+1), floor(max_no_na(log.P.values)+1)))
                    plot(NULL, xlim=c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,YlimMax),xlab="",ylab="")
                }else{
                    plot(NULL, xlim=c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=axis.cex, cex.lab=lab.cex,ylim=c(0,max(ylim[[i]])),xlab="",ylab="")      
                }
                axis(1, mgp=c(3,xticks.pos,0),at=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), lwd=axis.lwd,labels=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), cex.axis=axis.cex)
                axis(2, las=1,lwd=axis.lwd,cex.axis=axis.cex)
                axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(ylim[[i]]))), labels=c("",""), tcl=0, lwd=axis.lwd)

                mtext(side=1, text=expression(Expected~~-log[10](italic(p))), line=ylab.pos+1, cex=lab.cex, font=lab.font, xpd=TRUE)
                mtext(side=2, text=expression(Observed~~-log[10](italic(p))), line=ylab.pos, cex=lab.cex, font=lab.font, xpd=TRUE)
                
                #plot the confidence interval of QQ-plot
                if(conf.int){
                    if(is.null(conf.int.col)){
                        polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(t(col)[i])), alpha=points.alpha, maxColorValue=255))
                    }else{
                        polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255),border=rgb(t(col2rgb(conf.int.col[i])), alpha=points.alpha, maxColorValue=255))
                    }
                }

                if(!is.null(threshold.col)){par(xpd=FALSE); abline(a=0, b=1,lwd=threshold.lty[1], lty=threshold.lty[1], col=threshold.col[1]); par(xpd=TRUE)}
                # points(log.Quantiles, log.P.values, col=t(col)[i],pch=19,cex=cex[3])
                is_visable <- filter.points(log.Quantiles, log.P.values, wh, ht, dpi=dpi)
                if(!is.null(threshold[[i]])){
                    # if(sum(threshold!=0)==length(threshold)){
                        thre.line=-log10(min_no_na(threshold[[i]]))
                        if(amplify==TRUE){
                            thre.index <- log.P.values<thre.line
                            if(sum(!thre.index)!=0){
                                points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col=t(col)[i],pch=19,cex=cex[3])
                            
                                #cover the points that exceed the threshold with the color "white"
                                # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                # print(signal.col)
                                # print(signal.pch)
                                # print(signal.cex)
                                if(is.null(signal.col)){
                                    points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=t(col)[i],pch=signal.pch[i],cex=signal.cex[i])
                                }else{
                                    points(log.Quantiles[!thre.index],log.P.values[!thre.index],col=signal.col[i],pch=signal.pch[i],cex=signal.cex[i])
                                }
                            }else{
                                points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                            }
                        }else{
                            points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                        }
                    # }
                }else{
                    points(log.Quantiles[is_visable], log.P.values[is_visable], col=t(col)[i],pch=19,cex=cex[3])
                }
                if(!is.null(main)) {
                    title(main=main[i], cex.main=main.cex, font.main= main.font)
                }else{
                    title(main=trait[i], cex.main=main.cex, font.main= main.font) 
                }
                if(box) box(lwd=axis.lwd)
                if(file.output) dev.off()
            }
        }
    }
    if(file.output & verbose)   cat(paste(" Plots are stored in: ", getwd(), sep=""), "\n")
    if(!is.null(wind_snp_num))  return(invisible(wind_snp_num))
}
