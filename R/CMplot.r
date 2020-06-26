#Version: 3.6.2
#Data: 2020/06/26
#Author: Lilin Yin

CMplot <- function(
    Pmap,
    col=c("#4197d8", "#f8c120", "#413496", "#495226", "#d60b6f", "#e66519", "#d581b7", "#83d3ad", "#7c162c", "#26755d"),
    bin.size=1e6,
    bin.range=NULL,
    pch=19,
    type="p",
    band=1,
    H=1.5,
    ylim=NULL,
    cex.axis=1,
    lwd.axis=1.5,
    cex.lab=1.5,
    plot.type="b",
    multracks=FALSE,
    cex=c(0.5,1,1),
    r=0.3,
    outward=FALSE,
    ylab=expression(-log[10](italic(p))),
    ylab.pos=3,
    xticks.pos=1,
    mar = c(3,6,3,3),
    threshold = NULL, 
    threshold.col="red",
    threshold.lwd=1,
    threshold.lty=2,
    amplify= TRUE,
    signal.cex = 1.5,
    signal.pch = 19,
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
    highlight.text.xadj=NULL,
    highlight.text.yadj=NULL,
    highlight.text.font=3,
    chr.labels=NULL,
    chr.border=FALSE,
    chr.labels.angle=0,
    chr.den.col="black",
    cir.band=1,
    cir.chr=TRUE,
    cir.chr.h=1.5,
    cir.legend=TRUE,
    cir.legend.cex=0.6,
    cir.legend.col="black",
    LOG10=TRUE,
    box=FALSE,
    conf.int=TRUE,
    conf.int.col=NULL,
    file.output=TRUE,
    file=c("jpg","pdf","tiff"),
    dpi=300,
    height=NULL,
    width=NULL,
    memo="",
    main="",
    main.cex=1.5,
    main.font=2,
    trait.legend.ncol=NULL,
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
        xadj=NULL,
        yadj=NULL,
        words=NULL,
        point.cex = 1,
        text.cex = 1,
        pch = 19,
        type = "p",
        point.col = "red",
        text.col = "black",
        text.font = 3,
        xlim = c(-Inf, Inf),
        ylim = c(-Inf, Inf)
    )
    {
        overlap <- function(x1, y1, sw1, sh1, boxes) {
          s <- 0
          last <- 1
          if (length(boxes) == 0)
            return(FALSE)
          for (i in c(1:length(boxes))) {
            bnds <- boxes[[i]]
            x2 <- bnds[1]
            y2 <- bnds[2]
            sw2 <- bnds[3]
            sh2 <- bnds[4]
            if (x1 < x2)
              overlap <- x1 + sw1 > x2 - s
            else
              overlap <- x2 + sw2 > x1 - s
            
            if (y1 < y2)
              overlap <- overlap && (y1 + sh1 > y2 - s)
            else
              overlap <- overlap && (y2 + sh2 > y1 - s)
            if (overlap) {
              return(TRUE)
            }
          }
          FALSE
        }
    
        layout <- function(x,
                           y,
                           words,
                           cex = 1,
                           rotate90 = FALSE,
                           xlim = c(-Inf, Inf),
                           ylim = c(-Inf, Inf),
                           tstep = .1,
                           rstep = .1)
        {
          tails <- "g|j|p|q|y"
          n <- length(words)
          sdx <- sd(x, na.rm = TRUE)
          sdy <- sd(y, na.rm = TRUE)
          if (sdx == 0)
          sdx <- 1
          if (sdy == 0)
          sdy <- 1
          if (length(cex) == 1)
          cex <- rep(cex, n)
          if (length(rotate90) == 1)
          rotate90 <- rep(rotate90, n)


          boxes <- list()
          for (i in 1:length(words)) {
          rotWord <- rotate90[i]
          r <- 0
          theta <- runif(1, 0, 2 * pi)
          x1 <- xo <- x[i]
          y1 <- yo <- y[i]
          wid <- strwidth(words[i], cex = cex[i])
          ht <- 1.5 * strheight(words[i], cex = cex[i])
          #mind your ps and qs
          if (grepl(tails, words[i]))
            ht <- ht + ht * .2
          if (rotWord) {
            tmp <- ht
            ht <- wid
            wid <- tmp
          }
          isOverlaped <- TRUE
          # iter = 0
          # itermax = 20
          while (isOverlaped) {
            # iter = iter + 1
            if (!overlap(x1 - .5 * wid, y1 - .5 * ht, wid, ht, boxes) &&
                x1 - .5 * wid > xlim[1] && y1 - .5 * ht > ylim[1] &&
                x1 + .5 * wid < xlim[2] && y1 + .5 * ht < ylim[2]) {
              boxes[[length(boxes) + 1]] <- c(x1 - .5 * wid, y1 - .5 * ht, wid, ht)
              isOverlaped <- FALSE
            } else{
              theta <- theta + tstep
              r <- r + rstep * tstep / (2 * pi)
              x1 <- xo + sdx * r * cos(theta)
              y1 <- yo + sdy * r * sin(theta)
              # if(iter == itermax){
              #   boxes[[length(boxes) + 1]] <- c(x1 - .5 * wid, y1 - .5 * ht, wid, ht)
              #   isOverlaped <- FALSE
              # }
            }
          }
          }
          result <- do.call(rbind, boxes)
          colnames(result) <- c("x", "y", "width", "ht")
          rownames(result) <- words
          result
        }
        if(!is.null(words)){
            n <- length(x)
            indx <- order(y, decreasing = TRUE)
            x <- x[indx]
            y <- y[indx]
            if(length(point.cex)!=1){point.cex = rep(point.cex, n); point.cex = point.cex[indx]}
            if(length(pch)!=1){pch = rep(pch, n); pch = pch[indx]}
            if(length(point.col)!=1){point.col = rep(point.col, n); point.col = point.col[indx]}
            if(length(text.col)!=1){text.col = rep(text.col, n); text.col = text.col[indx]}
            if(length(text.cex)!=1){text.cex = rep(text.cex, n); text.cex = text.cex[indx]}   
            x1=x; y1=y
            if(length(x) != length(words))  stop("highlighted genes not equal to the highlighted SNPs.")
            words <- words[indx]
            if(is.null(xadj)){
                xadj = sample(c(1.5, 0, -0.5), n, replace=TRUE)
            }else{
                if(length(xadj) != n)   stop("length of xadj not equals to length of x")
                if(sum(!xadj %in% c(-1,0,1)) > 0)   stop("-1, 0, 1 limited for xadj")
                xadj[xadj==-1] = 1.5
                xadj[xadj==1] = -0.5
                xadj = xadj[indx]
            }
            if(is.null(yadj)){
                yadj = rep(-0.5, n)
                # yadj[y<=median(y)] = 1.5
                # yadj[y>median(y)] = -0.5
            }else{
                if(length(yadj) != n)   stop("length of yadj not equals to length of y")
                if(sum(!yadj %in% c(-1,0,1)) > 0)   stop("-1, 0, 1 limited for yadj")
                yadj[yadj==-1] = 1.5
                yadj[yadj==1] = -0.5
                yadj = yadj[indx]
            }
            for(i in 1:n){
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
                        x[i] = x[i] + 0.7*strwidth(words[i],cex=text.cex)
                    }else{
                        x[i] = x[i] - 0.7*strwidth(words[i],cex=text.cex)
                    }
                }
                if(xadj[i] == -0.5){
                    if((x[i] + 1.2*strwidth(words[i],cex=text.cex)) > max(xlim)){
                        x[i] = x[i] - 0.7*strwidth(words[i],cex=text.cex)
                    }else{
                        x[i] = x[i] + 0.7*strwidth(words[i],cex=text.cex)
                    }
                }
            }
            x <- c(x1,x)
            y <- c(y1,y)
            words <- c(rep("OO", n), as.character(words))
            lay <- layout(x,y,words,cex=text.cex,xlim=xlim,ylim=ylim)
            indd <- (n+1):length(x)
            for(i in indd){
                xl <- lay[i,1]
                yl <- lay[i,2]
                w <- lay[i,3]
                h <- lay[i,4]
                nx <- xl+.5*w
                ny <- yl+.5*h
                if((nx + 0.5 * strwidth(words[i],cex=text.cex)) < x1[i-n]){
                    nx = nx + 0.5 * strwidth(words[i],cex=text.cex)
                }else if((nx - 0.5 * strwidth(words[i],cex=text.cex)) > x1[i-n]){
                    nx = nx - 0.5 * strwidth(words[i],cex=text.cex)
                }
                if((ny + strheight(words[i],cex=text.cex)) < y1[i-n]){
                    ny = ny + 0.5 * strheight(words[i],cex=text.cex)
                }else if((ny - strheight(words[i],cex=text.cex)) > y1[i-n]){
                    ny = ny - 0.5 * strheight(words[i],cex=text.cex)
                }
                arrows(x1[i-n], y1[i-n], nx, ny, length=.08, angle=15, code=2, col="grey", lwd=2)
                # segments(x1[i], y1[i], nx, ny, col="grey", lwd=2)
            }
            if(type=="h"){
                points(x1,y1,pch = pch,type="h",col = point.col, lwd = point.cex+1)
                points(x1,y1,pch = pch,type="p",col = point.col, cex = point.cex)
            }else{
                points(x1,y1,pch = pch,type=type,col = point.col,cex = point.cex)
            }
            text(lay[indd,1]+.5*lay[indd,3],lay[indd,2]+.5*lay[indd,4],words[indd],xpd=TRUE,cex = text.cex,col=text.col,font=text.font)
        }else{
            if(type=="h"){
                points(x,y,pch = pch,type="h",col = point.col, lwd = point.cex+1)
                points(x,y,pch = pch,type="p",col = point.col, cex = point.cex)
            }else{
                points(x,y,pch = pch,type=type,col = point.col,cex = point.cex)
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
        x <- ceiling((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * w * dpi / scale)
        y <- ceiling((y - min(y, na.rm = TRUE)) / (max(y, na.rm = TRUE) - min(y, na.rm = TRUE)) * h * dpi / scale)
        index <- !duplicated(cbind(x, y))
    }

    DensityPlot <- function(
        map,
        col=c("darkgreen", "yellow", "red"),
        main="SNP Density",
        chr.labels=NULL, 
        bin=1e6,
        band=3,
        width=5,
        legend.len=10,
        legend.max=NULL,
        legend.min=1,
        legend.pt.cex=3,
        legend.cex=1,
        legend.y.intersp=1,
        legend.x.intersp=1,
        file.output=TRUE,
        plot=TRUE,
        dpi = NULL,
        wh = NULL,
        ht = NULL
    )
    {
        if(is.null(legend.min)) legend.min = 1
        if(is.null(col) | length(col) == 1){col=c("darkgreen", "yellow", "red")}
        map <- map[map[, 2] != "0", ]
        map <- as.matrix(map)
        map <- map[!is.na(map[, 2]), ]
        suppressWarnings(map <- map[!is.na(as.numeric(map[, 3])), ])
        #map <- map[map[, 3] != 0, ]
        suppressWarnings(max.chr <- max(as.numeric(map[, 2]), na.rm=TRUE))
        if(is.infinite(max.chr))    max.chr <- 0
        suppressWarnings(map.xy.index <- which(!as.numeric(map[, 2]) %in% c(0 : max.chr)))
        if(length(map.xy.index) != 0){
            chr.xy <- unique(map[map.xy.index, 2])
            for(i in 1:length(chr.xy)){
                map[map[, 2] == chr.xy[i], 2] <- max.chr + i
            }
        }
        map <- map[order(as.numeric(map[, 2]), as.numeric(map[, 3])), ]
        chr <- as.numeric(map[, 2])
        pos <- as.numeric(map[, 3])
        chr.num <- unique(chr)
        chorm.maxlen <- max(pos)
        bp <- ifelse(chorm.maxlen < 1e6, 1e3, 1e6)
        bp_label <- ifelse(bp == 1e3, "Kb", "Mb")
        main <- paste("The number of SNPs within ", bin / bp, bp_label, " window size", sep="")
        if(plot)    plot(NULL, xlim=c(0, chorm.maxlen + chorm.maxlen/10), ylim=c(0, length(chr.num) * band + band), main=main,axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
        pos.x <- list()
        col.index <- list()
        maxbin.num <- NULL
        for(i in 1 : length(chr.num)){
            pos.x[[i]] <- pos[which(chr == chr.num[i])]
            cut.len <- ceiling((max(pos.x[[i]]) - min(pos.x[[i]])) / bin)
            if(cut.len <= 1){
                maxbin.num <- c(maxbin.num,length(pos.x[[i]]))
                col.index[[i]] <- rep(length(pos.x[[i]]), length(pos.x[[i]]))
            }else{
                cut.r <- cut(pos.x[[i]], cut.len, labels=FALSE)
                eachbin.num <- table(cut.r)
                maxbin.num <- c(maxbin.num, max(eachbin.num))
                col.index[[i]] <- rep(eachbin.num, eachbin.num)
            }
        }
        Maxbin.num <- max(maxbin.num)
        maxbin.num <- Maxbin.num
        if(!is.null(legend.max)){
            maxbin.num <- legend.max
        }
        col=colorRampPalette(col)(maxbin.num - legend.min + 1)
        col.seg=NULL
        for(i in 1 : length(chr.num)){
            if(plot){
                polygon(c(0, 0, max(pos.x[[i]]), max(pos.x[[i]])), 
                c(-width/5 - band * (i - length(chr.num) - 1), width/5 - band * (i - length(chr.num) - 1), 
                width/5 - band * (i - length(chr.num) - 1), -width/5 - band * (i - length(chr.num) - 1)), col="grey95", border="grey95")
                rect(xleft = 0, ybottom = -width/5 - band * (i - length(chr.num) - 1), xright = max(pos.x[[i]]), ytop = width/5 - band * (i - length(chr.num) - 1), border="grey80")
            }
            if(!is.null(legend.max)){
                if(legend.max < Maxbin.num){
                    col.index[[i]][col.index[[i]] > legend.max] <- legend.max
                }
            }
            col.index[[i]][col.index[[i]] < legend.min] <- legend.min
            col.seg <- c(col.seg, col[round(col.index[[i]] * length(col) / maxbin.num)])
            if(file.output){
                is_visable <-  filter.points(pos.x[[i]], -width/5 - band * (i - length(chr.num) - 1), wh, ht, dpi = dpi)
            }else{
                is_visable <- rep(TRUE, length(pos.x[[i]]))
            }
            if(plot)    segments(pos.x[[i]][is_visable], -width/5 - band * (i - length(chr.num) - 1), pos.x[[i]][is_visable], width/5 - band * (i - length(chr.num) - 1), 
            col=col[round(col.index[[i]] * length(col) / (maxbin.num - legend.min + 1)) - legend.min + 1][is_visable], lwd=1)
        }
        if(length(map.xy.index) != 0){
            for(i in 1:length(chr.xy)){
                chr.num[chr.num == max.chr + i] <- chr.xy[i]
            }
        }
        chr.num <- rev(chr.num)
        if(plot){
            if(!is.null(chr.labels)){
                mtext(at=seq(band, length(chr.num) * band, band),text=chr.labels, side=2, las=2, font=1, cex=cex.axis*0.6, line=0.2, xpd=TRUE)
            }else{
                if(max.chr == 0)    mtext(at=seq(band, length(chr.num) * band, band),text=chr.num, side=2, las=2, font=1, cex=cex.axis*0.6, line=0.2, xpd=TRUE)
                if(max.chr != 0)    mtext(at=seq(band, length(chr.num) * band, band),text=paste("Chr", chr.num, sep=""), side=2, las=2, font=1, cex=cex.axis*0.6, line=0.2, xpd=TRUE)
            }
        }
        if(plot){
            xticks = seq(0, chorm.maxlen / bp, length=10)
            
            if(round(xticks[2]) <= 10){
                xticks = seq(0, chorm.maxlen / bp, round(xticks[2], 1))
            }else{
                xticks = seq(0, chorm.maxlen / bp, round(xticks[2]))    
            }
            
            if((chorm.maxlen/bp - max(xticks)) > 0.5*xticks[2]){
                xticks = c(xticks, round(chorm.maxlen / bp))
            }
            axis(3, at=xticks*bp, labels=paste(xticks, bp_label, sep=""),
            font=1, cex.axis=cex.axis*0.8, tck=0.01, lwd=lwd.axis, padj=1.2)
            axis(3, at=c(0, chorm.maxlen), labels=c("",""), tcl=0, lwd=lwd.axis)
        }
        # image(c(chorm.maxlen-chorm.maxlen * legend.width / 20 , chorm.maxlen), 
        # round(seq(band - width/5, (length(chr.num) * band + band) * legend.height / 2 , length=maxbin.num+1), 2), 
        # t(matrix(0 : maxbin.num)), col=c("white", rev(heat.colors(maxbin.num))), add=TRUE)
                
        if(maxbin.num <= legend.len)    legend.len <- maxbin.num        
        
        legend.y <- round(seq(legend.min, maxbin.num, length=legend.len))
        legend.y <- unique(legend.y)
        len <- ifelse(length(legend.y)==1, 1, legend.y[2] - legend.y[1])
        legend.y <- seq(legend.min, maxbin.num, len)
        if(!is.null(legend.max)){
            if(legend.max < Maxbin.num){
                if(!maxbin.num %in% legend.y){
                    legend.y <- c(legend.y, paste(">=", maxbin.num, sep=""))
                    legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
                }else{
                    legend.y[length(legend.y)] <- paste(">=", maxbin.num, sep="")
                    legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
                }
            }else{
                if(!maxbin.num %in% legend.y){
                    legend.y <- c(legend.y, maxbin.num)
                }
                legend.y.col <- c(legend.y)
            }
        }else{
            if(!maxbin.num %in% legend.y){
                legend.y <- c(legend.y, paste(">", max(legend.y), sep=""))
                legend.y.col <- c(legend.y[c(-length(legend.y))], maxbin.num)
            }else{
                legend.y.col <- c(legend.y)
            }
        }
        if(legend.min != 1){
            legend.y[1] <- paste("<=", legend.min, sep="")
        }
        legend.y <- c(0, legend.y)
        legend.y.col <- as.numeric(legend.y.col)
        legend.col <- c("grey95", col[round(legend.y.col * length(col) / (maxbin.num - legend.min + 1)) - legend.min + 1])
        if(plot)    legend(x=(chorm.maxlen + chorm.maxlen/100), y=( -width/2.5 - band * (length(chr.num) - length(chr.num) - 1)), title="", legend=legend.y, pch=15, pt.cex = legend.pt.cex, col=legend.col,
            cex=legend.cex, bty="n", y.intersp=legend.y.intersp, x.intersp=legend.x.intersp, yjust=0, xjust=0, xpd=TRUE)
        if(!plot)   return(list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y))
    }

    if(!all(plot.type %in% c("b","c","m","q","d"))) stop("unknown 'plot.type'.")
    if(sum(plot.type %in% "b")==1) plot.type=c("c","m","q","d")
    file=match.arg(file)
    trait=colnames(Pmap)[-c(1:3)]
    if(!is.null(memo) && memo != "")    memo <- paste("_", memo, sep="")
    if(length(trait) == 0)   trait <- paste("Col", 1:(ncol(Pmap)-3), sep="")
    taxa <- paste(trait, memo, sep="")

    #SNP-Density plot
    if("d" %in% plot.type){
        if(verbose) cat(" SNP-Density Plotting.\n")
        if(file.output){
            ht=ifelse(is.null(height), 6, height)
            wh=ifelse(is.null(width), 9, width)
            if(file=="jpg") jpeg(paste("SNP-Density.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
            if(file=="pdf") pdf(paste("SNP-Density.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
            if(file=="tiff")    tiff(paste("SNP-Density.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
            # par(xpd=TRUE)
            par(mar=c(mar[1], mar[2], mar[3]+1, mar[4]))
        }else{
            ht=ifelse(is.null(height), 6, height)
            wh=ifelse(is.null(width), 9, width)
            if(is.null(dev.list())) dev.new(width = wh,height=ht)
            # par(xpd=TRUE)
        }
        if(!is.null(bin.range)){
            if(length(bin.range) != 2)  stop("Two values (min and max) should be provided for bin.range!")
            if(bin.range[1] == 0)   stop("Min value of bin.range should be more than 1!")
        }
        DensityPlot(map=Pmap[,c(1:3)], file.output = file.output, dpi = dpi, wh = wh, ht = ht, chr.labels = chr.labels, col=chr.den.col, bin=bin.size, legend.min=bin.range[1], legend.max=bin.range[2])
        if(file.output) dev.off()
    }

    if(length(plot.type) !=1 | (!"d" %in% plot.type)){
    
        #order Pmap by the name of SNP
        #Pmap=Pmap[order(Pmap[,1]),]
        suppressWarnings(Pmap <- Pmap[Pmap[, 2] != "0", ])
        Pmap <- as.matrix(Pmap)
        Pmap <- Pmap[!is.na(Pmap[, 2]), ]
        suppressWarnings(Pmap <- Pmap[!is.na(as.numeric(Pmap[, 3])), ])

        #scale and adjust the parameters
        cir.chr.h <- cir.chr.h/5
        cir.band <- cir.band/5
        if(!is.null(threshold)){
            if(LOG10){
                if(sum(threshold <= 0) != 0) stop("threshold must be greater than 0.")
            }
            threshold.col <- rep(threshold.col,length(threshold))
            threshold.lwd <- rep(threshold.lwd,length(threshold))
            threshold.lty <- rep(threshold.lty,length(threshold))
            signal.col <- rep(signal.col,length(threshold))
            signal.pch <- rep(signal.pch,length(threshold))
            signal.cex <- rep(signal.cex,length(threshold))
        }
        if(length(cex)!=3) cex <- rep(cex,3)

        #get the number of traits
        R=ncol(Pmap)-3

        if(!is.null(ylim)){
            if(!is.list(ylim)){
                cat(" (warning: all phenotypes will use the same ylim.)\n")
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
        if(all(main != "")) main <- rep(main, R)
        if(length(mar) != 4)    stop("length of 'mar' shoud equal to 4.")
        if(chr.labels.angle > 90 | chr.labels.angle < -90)  stop("'chr.labels.angle' should be > -90 and < 90.")
        # if(is.na(conf.int.col)){
        #   conf.int=FALSE
        # }else{
        #   conf.int=TRUE
        # }

        pch=rep(pch, R)

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
        Pmap <- Pmap[,-1]

        Pmap <- matrix(as.numeric(Pmap), nrow(Pmap))
        order_index <- order(Pmap[, 1], Pmap[,2])

        #order the GWAS results by chromosome and position
        Pmap <- Pmap[order_index, ]
        SNP_id <- SNP_id[order_index]

        if(!is.null(highlight)){
            highlight_index <- list()
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
                }else{
                    highlight[[i]] <- highlight[[i]][!is.na(highlight[[i]])]
                    highlight_index[[i]] <- match(as.character(as.matrix(highlight[[i]])), SNP_id)
                    if(all(is.na(highlight_index[[i]]))) stop("No shared SNPs between Pmap and highlight!")
                    highlight_index[[i]] <- na.omit(highlight_index[[i]])
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
            for(i in 1:length(highlight.text)){
                highlight.text[[i]] <- highlight.text[[i]][!is.na(highlight.text[[i]])]
            }
        }

        if(!is.null(highlight.text.xadj)){
            if(!is.list(highlight.text.xadj)){
                highlight.text.xadj <- list(highlight.text.xadj)
                for(i in 1:R){highlight.text.xadj[[i]] = highlight.text.xadj[[1]]}
            }else{
                if(length(highlight.text.xadj) != R){stop("length of 'highlight.text.xadj' not equals to the number of traits.")}  
            }  
        }

        if(!is.null(highlight.text.yadj)){
            if(!is.list(highlight.text.yadj)){
                highlight.text.yadj <- list(highlight.text.yadj)
                for(i in 1:R){highlight.text.yadj[[i]] = highlight.text.yadj[[1]]}
            }else{
                if(length(highlight.text.yadj) != R){stop("length of 'highlight.text.yadj' not equals to the number of traits.")}  
            }  
        }

        #get the index of chromosome
        chr <- unique(Pmap[,1])
        chr.ori <- chr
        if(length(map.xy.index) != 0){
            for(i in 1:length(chr.xy)){
                chr.ori[chr.ori == max.chr + i] <- chr.xy[i]
            }
        }

        pvalueT <- as.matrix(Pmap[,-c(1:2)])
        pvalue.pos <- Pmap[, 2]
        # p0.index <- Pmap[, 1] == 0
        # if(sum(p0.index) != 0){
            # pvalue.pos[p0.index] <- 1:sum(p0.index)
        # }
        pvalue.pos.list <- tapply(pvalue.pos, Pmap[, 1], list)
        
        #scale the space parameter between chromosomes
        if(!missing(band)){
            band <- floor(band*(sum(sapply(pvalue.pos.list, max))/100))
        }else{
            band <- floor((sum(sapply(pvalue.pos.list, max))/100))
        }
        if(band==0) band=1
        
        if(LOG10){
            pvalueT[pvalueT <= 0] <- NA
            pvalueT[pvalueT > 1] <- NA
        }
        Pmap[,-c(1:2)] <- pvalueT

        #set the colors for the plot
        #palette(heat.colors(1024)) #(heatmap)
        #T=floor(1024/max(pvalue))
        #plot(pvalue,pch=19,cex=0.6,col=(1024-floor(pvalue*T)))
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
            chr.border.pos = chr.border.pos[-length(chr.border.pos)]
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
        #NewP[[j]] <- pvalue
        
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
            density.list <- DensityPlot(map=Pmap[,c(1,1,2)], file.output = FALSE, col=chr.den.col, plot=FALSE, bin=bin.size, legend.min=bin.range[1], legend.max=bin.range[2])
            #list(den.col=col.seg, legend.col=legend.col, legend.y=legend.y)
        }else{
            cir.density=FALSE
        }
        
        signal.line.index <- NULL
        if(!is.null(threshold)){
            if(!is.null(signal.line)){
                for(l in 1:R){
                    if(LOG10){
                        signal.line.index <- c(signal.line.index,which(pvalueT[,l] < min_no_na(threshold)))
                    }else{
                        signal.line.index <- c(signal.line.index,which(pvalueT[,l] > max_no_na(threshold)))
                    }
                }
                signal.line.index <- unique(signal.line.index)
            }
        }
        signal.line.index <- pvalue.posN[signal.line.index]
    }
    
    #plot circle Manhattan
    if("c" %in% plot.type){
        #print("Starting Circular-Manhattan plot!",quote=F)
        if(file.output){
            ht=ifelse(is.null(height), 10, height)
            wh=ifelse(is.null(width), 10, width)
            if(file=="jpg") jpeg(paste("Circular-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
            if(file=="pdf") pdf(paste("Circular-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
            if(file=="tiff")    tiff(paste("Circular-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
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
                X1chr <- (RR)*sin(2*pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                Y1chr <- (RR)*cos(2*pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                X2chr <- (r)*sin(2*pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                Y2chr <- (r)*cos(2*pi*(signal.line.index-round(band/2)-circleMin)/TotalN)
                segments(X1chr,Y1chr,X2chr,Y2chr,lty=2,lwd=signal.line,col="grey")
            }
        }
        for(i in 1:R){
        
            #get the colors for each trait
            colx <- col[i,]
            colx <- colx[!is.na(colx)]
            
            #debug
            #print(colx)
            
            if(verbose) cat(paste(" Circular-Manhattan Plotting ",trait[i],".\n",sep=""))
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
                            X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
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
                            X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
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
                            is_visable <- filter.points((RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), (RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), wh, ht, dpi = dpi)
                        }else{
                            is_visable <- rep(TRUE, length(pvalue.posN))
                        }
                        segments(
                            (RR)*sin(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR)*cos(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            col=density.list$den.col[is_visable], lwd=0.5
                        )
                        legend(
                            x=RR+4*cir.chr.h,
                            y=(RR+4*cir.chr.h)/2,
                            title="", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
                            cex=1, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.5, xjust=0, xpd=TRUE
                        )
                        
                    }
                    
                    # XLine=(RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
                    # YLine=(RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
                    # lines(XLine,YLine,lwd=1.5)
                    if(cir.density){
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE,col='grey')
                        circle.plot(myr=RR,lwd=1.5,add=TRUE,col='grey')
                    }else{
                        circle.plot(myr=RR+cir.chr.h,lwd=1.5,add=TRUE)
                        circle.plot(myr=RR,lwd=1.5,add=TRUE)
                    }

                }
                
                X=(Cpvalue[ylimIndx]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                Y=(Cpvalue[ylimIndx]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                if(file.output){
                    is_visable <- filter.points(X, Y, wh, ht, dpi = dpi)
                }else{
                    is_visable <- rep(TRUE, length(X))
                }

                if(cir.legend==TRUE){
                    circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                }

                points(X[is_visable],Y[is_visable],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]])[ylimIndx][is_visable])
                
                #plot the legend for each trait
                if(cir.legend==TRUE){
                    #try to get the number after radix point
                    if((Max-Min) > 1) {
                        round.n=2
                    }else{
                        round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                    }
                    segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    text(-r/15,r+H*(i-0.94)+cir.band*(i-1),round(Min+(Max-Min)*0,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Min+(Max-Min)*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Min+(Max-Min)*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Min+(Max-Min)*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.06)+cir.band*(i-1),round(Min+(Max-Min)*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                }
                
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                        for(thr in 1:length(threshold)){
                            significantline1=ifelse(LOG10, H*(-log10(threshold[thr])-Min)/(Max-Min), H*(threshold[thr]-Min)/(Max-Min))
                            #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
                            #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
                            if(significantline1<H){
                                #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                circle.plot(myr=(significantline1+r+H*(i-1)+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                            }else{
                                warning(paste("No significant points for ",trait[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
                            }
                        }
                    }
                }
                
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                        if(amplify==TRUE){
                            if(LOG10){
                                threshold <- sort(threshold)
                                significantline1=H*(-log10(max_no_na(threshold))-Min)/(Max-Min)
                            }else{
                                threshold <- sort(threshold, decreasing=TRUE)
                                significantline1=H*(min_no_na(threshold)-Min)/(Max-Min)
                            }
                            
                            p_amp.index <- which(Cpvalue>=significantline1)
                            HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            
                            #cover the points that exceed the threshold with the color "white"
                            points(HX1,HY1,pch=19,cex=cex[1],col="white")
                            
                                for(ll in 1:length(threshold)){
                                    if(ll == 1){
                                        if(LOG10){
                                            significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
                                        }else{
                                            significantline1=H*(threshold[ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1)
                                        HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                        HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                    }else{
                                        if(LOG10){
                                            significantline0=H*(-log10(threshold[ll-1])-Min)/(Max-Min)
                                            significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
                                        }else{
                                            significantline0=H*(threshold[ll-1]-Min)/(Max-Min)
                                            significantline1=H*(threshold[ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                        HX1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                        HY1=(Cpvalue[p_amp.index]+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
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
                    points(X[highlight_index[[i]]],X[highlight_index[[i]]],pch=19,cex=cex[1],col="white")
                    if(is.null(highlight.col))  highlight.col = rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]]
                    points(X[highlight_index[[i]]],Y[highlight_index[[i]]],pch=highlight.pch,cex=highlight.cex,col=highlight.col)
                }

                if(cir.chr==TRUE){
                    ticks1=(RR+1.5*cir.chr.h)*sin(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=(RR+1.5*cir.chr.h)*cos(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                            angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                            text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=2,cex=cex.axis, adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab, sep=""),srt=angle, adj=c(0.5, 0),font=2,cex=cex.axis)
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=2,cex=cex.axis, adj=c(0.5, 0))
                            }
                        }
                    }
                }else{
                    ticks1=1.01*RR*sin(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=1.01*RR*cos(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    # ticks1=(0.9*r)*sin(2*pi*(ticks-round(band/2))/TotalN)
                    # ticks2=(0.9*r)*cos(2*pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                        angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                        text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab, sep=""),srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                            }
                        }
                    }
                }
            }
            if(outward==FALSE){
                if(cir.chr==TRUE & i == 1){
                    # XLine=(2*cir.band+RR+cir.chr.h)*sin(2*pi*(1:TotalN)/TotalN)
                    # YLine=(2*cir.band+RR+cir.chr.h)*cos(2*pi*(1:TotalN)/TotalN)
                    # lines(XLine,YLine,lwd=1.5)

                    polygon.num <- 1000
                    for(k in 1:length(chr)){
                        if(k==1){
                            polygon.index <- seq(round(band/2)+1,-round(band/2)-circleMin+max_no_na(pvalue.posN.list[[1]]), length=polygon.num)
                            X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
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
                            X1chr=(RR)*sin(2*pi*(polygon.index)/TotalN)
                            Y1chr=(RR)*cos(2*pi*(polygon.index)/TotalN)
                            X2chr=(RR+cir.chr.h)*sin(2*pi*(polygon.index)/TotalN)
                            Y2chr=(RR+cir.chr.h)*cos(2*pi*(polygon.index)/TotalN)
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
                            is_visable <- filter.points((RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), (RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN), wh, ht, dpi = dpi)
                        }else{
                            is_visable <- rep(TRUE, length(pvalue.posN))
                        }
                        segments(
                            (RR)*sin(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR)*cos(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*sin(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            (RR+cir.chr.h)*cos(2*pi*(pvalue.posN-round(band/2)-circleMin)/TotalN)[is_visable],
                            col=density.list$den.col[is_visable], lwd=0.5
                        )
                        legend(
                            x=RR+4*cir.chr.h,
                            y=(RR+4*cir.chr.h)/2,
                            title="", legend=density.list$legend.y, pch=15, pt.cex = 3, col=density.list$legend.col,
                            cex=1, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=0.5, xjust=0, xpd=TRUE
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

                X=(-Cpvalue[ylimIndx]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                Y=(-Cpvalue[ylimIndx]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[ylimIndx]-round(band/2)-circleMin)/TotalN)
                if(file.output){
                    is_visable <- filter.points(X, Y, wh, ht, dpi = dpi)
                }else{
                    is_visable <- rep(TRUE, length(X))
                }

                if(cir.legend==TRUE){
                    circle.plot(myr=r+H*(i-1)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.75)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.5)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0.25)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                    circle.plot(myr=r+H*(i-0)+cir.band*(i-1),lwd=0.5,add=TRUE,col='grey')
                }

                points(X[is_visable],Y[is_visable],pch=19,cex=cex[1],col=rep(rep(colx,N[i]),add[[i]])[ylimIndx][is_visable])
                
                if(cir.legend==TRUE){
                    
                    #try to get the number after radix point
                    if((Max-Min)<=1) {
                        round.n=nchar(as.character(10^(-ceiling(-log10(Max)))))-1
                    }else{
                        round.n=2
                    }
                    segments(0,r+H*(i-1)+cir.band*(i-1),0,r+H*i+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-1)+cir.band*(i-1),H/20,r+H*(i-1)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0.75)+cir.band*(i-1),H/20,r+H*(i-0.75)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0.5)+cir.band*(i-1),H/20,r+H*(i-0.5)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0.25)+cir.band*(i-1),H/20,r+H*(i-0.25)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    segments(0,r+H*(i-0)+cir.band*(i-1),H/20,r+H*(i-0)+cir.band*(i-1),col=cir.legend.col,lwd=1.5)
                    text(-r/15,r+H*(i-0.06)+cir.band*(i-1),round(Min+(Max-Min)*0,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.25)+cir.band*(i-1),round(Min+(Max-Min)*0.25,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.5)+cir.band*(i-1),round(Min+(Max-Min)*0.5,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.75)+cir.band*(i-1),round(Min+(Max-Min)*0.75,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                    text(-r/15,r+H*(i-0.94)+cir.band*(i-1),round(Min+(Max-Min)*1,round.n),adj=1,col=cir.legend.col,cex=cir.legend.cex,font=2)
                }
                
                if(!is.null(threshold)){
                    if(sum(threshold!=0)==length(threshold)){
                    
                        for(thr in 1:length(threshold)){
                            significantline1=ifelse(LOG10, H*(-log10(threshold[thr])-Min)/(Max-Min), H*(threshold[thr]-Min)/(Max-Min))
                            #s1X=(significantline1+r+H*(i-1)+cir.band*(i-1))*sin(2*pi*(0:TotalN)/TotalN)
                            #s1Y=(significantline1+r+H*(i-1)+cir.band*(i-1))*cos(2*pi*(0:TotalN)/TotalN)
                            if(significantline1<H){
                                #lines(s1X,s1Y,type="l",col=threshold.col,lwd=threshold.col,lty=threshold.lty)
                                circle.plot(myr=(-significantline1+r+H*i+cir.band*(i-1)),col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                            }else{
                                warning(paste("No significant points for ",trait[i]," pass the threshold level using threshold=",threshold[thr],"!",sep=""))
                            }
                        }
                        if(amplify==TRUE){
                            if(LOG10){
                                threshold <- sort(threshold)
                                significantline1=H*(-log10(max_no_na(threshold))-Min)/(Max-Min)
                            }else{
                                threshold <- sort(threshold, decreasing=TRUE)
                                significantline1=H*(min_no_na(threshold)-Min)/(Max-Min)
                            }
                            p_amp.index <- which(Cpvalue>=significantline1)
                            HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                            
                            #cover the points that exceed the threshold with the color "white"
                            points(HX1,HY1,pch=19,cex=cex[1],col="white")
                            
                                for(ll in 1:length(threshold)){
                                    if(ll == 1){
                                        if(LOG10){
                                            significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
                                        }else{
                                            significantline1=H*(threshold[ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1)
                                        HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                        HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                    }else{
                                        if(LOG10){
                                            significantline0=H*(-log10(threshold[ll-1])-Min)/(Max-Min)
                                            significantline1=H*(-log10(threshold[ll])-Min)/(Max-Min)
                                        }else{
                                            significantline0=H*(threshold[ll-1]-Min)/(Max-Min)
                                            significantline1=H*(threshold[ll]-Min)/(Max-Min)
                                        }
                                        p_amp.index <- which(Cpvalue>=significantline1 & Cpvalue < significantline0)
                                        HX1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*sin(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                        HY1=(-Cpvalue[p_amp.index]+r+H*i+cir.band*(i-1))*cos(2*pi*(pvalue.posN[p_amp.index]-round(band/2)-circleMin)/TotalN)
                                    
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
                    points(X[highlight_index[[i]]],X[highlight_index[[i]]],pch=19,cex=cex[1],col="white")
                    if(is.null(highlight.col))  highlight.col = rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]]
                    points(X[highlight_index[[i]]],Y[highlight_index[[i]]],pch=highlight.pch,cex=highlight.cex,col=highlight.col)
                }

                if(cir.chr==TRUE){
                    ticks1=(RR+1.5*cir.chr.h)*sin(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=(RR+1.5*cir.chr.h)*cos(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                          angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                          text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab,sep=""),srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                            }
                        }
                    }
                }else{
                    ticks1=1.01*RR*sin(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    ticks2=1.01*RR*cos(2*pi*(ticks-round(band/2)-circleMin)/TotalN)
                    # ticks1=RR*sin(2*pi*(ticks-round(band/2))/TotalN)
                    # ticks2=RR*cos(2*pi*(ticks-round(band/2))/TotalN)
                    if(is.null(chr.labels)){
                        for(t in 1:length(ticks)){
                        
                            #adjust the angle of labels of circle plot
                            angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                            text(ticks1[t],ticks2[t],chr.ori[t],srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                        }
                    }else{
                        if(Nchr == 1){
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],paste(chr.labels[t], bp_lab,sep=""),srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
                            }
                        }else{
                            for(t in 1:length(ticks)){
                                angle=360*(1-(ticks-round(band/2)-circleMin)[t]/TotalN)
                                text(ticks1[t],ticks2[t],chr.labels[t],srt=angle,font=2,cex=cex.axis,adj=c(0.5, 0))
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
                is_visable[[i]] <- filter.points(pvalue.posN, logpvalueT[,i], wh, ht, dpi = dpi)
            }else{
                is_visable[[i]] <- rep(TRUE, nrow(logpvalueT))
            }
        }

        if(multracks==FALSE){
            #print("Starting Rectangular-Manhattan plot!",quote=F)
            for(i in 1:R){
                colx=col[i,]
                colx=colx[!is.na(colx)]
                if(verbose) cat(paste(" Rectangular-Manhattan Plotting ",trait[i],".\n",sep=""))
                    if(file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(file=="jpg") jpeg(paste("Rectangular-Manhattan.",taxa[i],".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
                        if(file=="pdf") pdf(paste("Rectangular-Manhattan.",taxa[i],".pdf",sep=""), width = wh,height=ht)
                        if(file=="tiff")    tiff(paste("Rectangular-Manhattan.",taxa[i],".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
                        par(mar = mar,xaxs="i",yaxs="r")
                    }
                    if(!file.output){
                        ht=ifelse(is.null(height), 6, height)
                        wh=ifelse(is.null(width), 14, width)
                        if(is.null(dev.list())) dev.new(width = wh, height = ht)
                        # par(xpd=TRUE)
                    }
                    
                    pvalue=pvalueT[,i]
                    logpvalue=logpvalueT[,i]
                    if(is.null(ylim)){
                        if(!is.null(threshold)){
                            if(sum(threshold!=0)==length(threshold)){
                                if(LOG10 == TRUE){
                                    Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),(-log10(min_no_na(threshold))))))
                                    Min <- min_ylim(min_no_na(c(-log10((max_no_na(pvalue))),-log10(max_no_na(threshold)))))
                                }else{
                                    Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold))))
                                    #if(abs(Max)<=1)    Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
                                    Min <- min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold))))
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
                                    cex.axis=cex.axis,font=2,axes=FALSE)
                            }else{
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=cex.axis,font=2,axes=FALSE)
                            }
                        }else{
                            if(cir.density){
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ann=FALSE,
                                cex.axis=cex.axis,font=2,axes=FALSE)
                            }else{
                                plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                                cex.axis=cex.axis,font=2,axes=FALSE)
                            }
                        }
                        mtext(side = 2, text = ylab, line = ylab.pos, cex=cex.lab, font=1, xpd=TRUE)
                    }else{
                        Max <- max_no_na(ylim[[i]])
                        Min <- min_no_na(ylim[[i]])
                        if(cir.density){
                            plot(pvalue.posN[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(min_no_na(ylim[[i]])-(Max-Min)/den.fold, max_no_na(ylim[[i]])),ann=FALSE,
                            cex.axis=cex.axis,font=2,axes=FALSE)
                        }else{
                            plot(pvalue.posN[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]+1,cex=cex[2],col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=ylim[[i]],ann=FALSE,
                            cex.axis=cex.axis,font=2,axes=FALSE)
                        }
                        mtext(side = 2, text = ylab, line = ylab.pos, cex=cex.lab, font=1, xpd=TRUE)
                    }
                    # Max1 <- Max
                    # Min1 <- Min
                    # if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
                    # if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
                    if(chr.border){
                        for(b in 1:length(chr.border.pos)){
                            segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=lwd.axis, lty=2)
                        }
                    }

                    if(chr.labels.angle == 0){
                        if(!is.null(chr.labels)){
                            if(Nchr == 1){
                                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                            }else{
                                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.labels))
                                #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
                            }
                        }else{
                            axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.ori))
                            #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
                        }
                    }else{
                        axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis,labels=FALSE)
                        if(!is.null(chr.labels)){
                            if(Nchr == 1){
                                text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis, font = 2, labels = c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %% 90 == 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                            }else{
                                # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=)
                                #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
                                text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis, font = 2, labels = c("Chr",chr.labels), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %%90== 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                            }
                        }else{
                            # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.ori))
                            #axis(1, at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
                            text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis, font = 2, labels = c("Chr",chr.ori), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %%90== 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                        }
                    }
                    axis(1, mgp=c(3,xticks.pos,0), at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
                    if(is.null(ylim)){
                        if((Max-Min)>1){
                            axis(2, las=1, lwd=lwd.axis,cex.axis=cex.axis,font=2)
                            axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
                            legend.y <- Max
                        }else{
                            axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
                            axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
                            legend.y <- Max
                        }
                    }else{
                        axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
                        axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
                        legend.y <- tail(ylim[[i]][2], 1)
                    }
                    if(!is.null(threshold)){
                        for(thr in 1:length(threshold)){
                            h <- ifelse(LOG10, -log10(threshold[thr]), threshold[thr])
                            # print(h)
                            # print(threshold.col[thr])
                            # print(threshold.lty[thr])
                            # print(threshold.lwd[thr])
                            segments(0, h, max_no_na(pvalue.posN), h,col=threshold.col[thr],lty=threshold.lty[thr],lwd=threshold.lwd[thr])
                        }
                        if(amplify == TRUE){
                            if(LOG10){
                                threshold <- sort(threshold)
                                sgline1=-log10(max_no_na(threshold))
                            }else{
                                threshold <- sort(threshold, decreasing=TRUE)
                                sgline1=min_no_na(threshold)
                            }

                            sgindex=which(logpvalue>=sgline1)
                            HY1=logpvalue[sgindex]
                            HX1=pvalue.posN[sgindex]
                            
                            #cover the points that exceed the threshold with the color "white"
                            points(HX1,HY1,pch=pch,cex=cex[2],col="white")
                            
                            for(ll in 1:length(threshold)){
                                if(ll == 1){
                                    if(LOG10){
                                        sgline1=-log10(threshold[ll])
                                    }else{
                                        sgline1=threshold[ll]
                                    }
                                    sgindex=which(logpvalue>=sgline1)
                                    HY1=logpvalue[sgindex]
                                    HX1=pvalue.posN[sgindex]
                                }else{
                                    if(LOG10){
                                        sgline0=-log10(threshold[ll-1])
                                        sgline1=-log10(threshold[ll])
                                    }else{
                                        sgline0=threshold[ll-1]
                                        sgline1=threshold[ll]
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
                        points(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],pch=pch,cex=cex[2],col="white")
                        if(!is.na(highlight_index[[i]][1])){
                            if(is.null(highlight.col)){
                                highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),xadj=highlight.text.xadj[[i]],yadj=highlight.text.yadj[[i]],words=highlight.text[[i]],point.cex=highlight.cex,text.cex=highlight.text.cex, pch=highlight.pch,type=highlight.type,point.col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]],text.col=highlight.text.col,text.font=highlight.text.font)
                            }else{
                                highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),xadj=highlight.text.xadj[[i]],yadj=highlight.text.yadj[[i]],words=highlight.text[[i]],point.cex=highlight.cex,text.cex=highlight.text.cex, pch=highlight.pch,type=highlight.type,point.col=highlight.col,text.col=highlight.text.col,text.font=highlight.text.font)
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
                        is_visable_den <- filter.points(pvalue.posN, ymin-0.5*(Max-Min)/den.fold, wh, ht, dpi = dpi)
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
                            title="", legend=density.list$legend.y, pch=15, pt.cex = 2.5, col=density.list$legend.col,
                            cex=0.8, bty="n",
                            y.intersp=1,
                            x.intersp=1,
                            yjust=1, xjust=0, xpd=TRUE
                        )
                        
                    }
                if(all(main != ""))  title(main = main[i], cex.main = main.cex, font.main= main.font)
                if(box) box(lwd=lwd.axis)
                #if(!is.null(threshold) & (length(grep("FarmCPU",taxa[i])) != 0))   abline(v=which(pvalueT[,i] < min_no_na(threshold)/max_no_na(dim(Pmap))),col="grey",lty=2,lwd=signal.line)
                if(file.output)  dev.off()
            }
            #print("Rectangular-Manhattan has been finished!",quote=F)
        }else{
            #print("Starting Rectangular-Manhattan plot!",quote=F)
            #print("Plotting in multiple tracks!",quote=F)
            if(file.output){
                ht=ifelse(is.null(height), 6, height)
                wh=ifelse(is.null(width), 14, width)
                if(file=="jpg") jpeg(paste("Multracks-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi*R,res=dpi,quality = 100)
                if(file=="pdf") pdf(paste("Multracks-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht*R)
                if(file=="tiff")    tiff(paste("Multracks-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi*R,res=dpi)
                par(mfcol=c(R,1), xaxs="i")
            }
            if(!file.output){
                ht=ifelse(is.null(height), 6, height)
                wh=ifelse(is.null(width), 14, width)
                if(is.null(dev.list())) dev.new(width = wh, height = ht)
                # par(xpd=TRUE)
            }
            for(i in 1:R){
                if(i == 1)  par(mar=c(0, mar[2]+1, mar[3], mar[4]))
                if(i == R)  par(mar=c(mar[1]+1, mar[2]+1, 0, mar[4]))
                if(i != 1 & i != R) par(mar=c(0, mar[2]+1, 0, mar[4]))
                if(verbose) cat(paste(" Multracks-Manhattan Plotting ",trait[i],".\n",sep=""))
                colx=col[i,]
                colx=colx[!is.na(colx)]
                pvalue=pvalueT[,i]
                logpvalue=logpvalueT[,i]
                if(is.null(ylim)){
                    if(!is.null(threshold)){
                        # if(sum(threshold!=0)==length(threshold)){
                            if(LOG10){
                                Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),-log10(min_no_na(threshold)))))
                                Min <- min_ylim(min_no_na(c((-log10(max_no_na(pvalue))),-log10(max_no_na(threshold)))))
                            }else{
                                Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold))))
                                #if(abs(Max)<=1)    Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
                                Min<-min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold))))
                                #if(abs(Min)<=1)    Min=min_no_na(min_no_na(pvalue),min_no_na(threshold))
                            }
                        # }else{
                        #     if(LOG10){
                        #         Max=max_ylim((-log10(min_no_na(pvalue))))
                        #         Min<-0
                        #     }else{
                        #         Max=max_ylim((max_no_na(pvalue)))
                        #         #if(abs(Max)<=1)    Max=max_no_na(max_no_na(pvalue))
                        #         Min=min_ylim((min_no_na(pvalue)))
                        #         #if(abs(Min)<=1)    Min=min_no_na(min_no_na(pvalue))
                        #         # }else{
                        #             # Max=max_no_na(ceiling(max_no_na(pvalue)))
                        #         # }
                        #     }   
                        # }
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
                            cex.axis=cex.axis*(R/2),font=2,axes=FALSE,yaxs="r")
                    }else{
                        plot(pvalue.posN[is_visable[[i]]],logpvalue[is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=c(Min,Max),ann=FALSE,
                            cex.axis=cex.axis*(R/2),font=2,axes=FALSE,yaxs="r")
                    }
                    mtext(side = 2, text = ylab, line = ylab.pos, cex=cex.lab*(R/2), font=1, xpd=TRUE)
                }else{
                    Max <- max_no_na(ylim[[i]])
                    Min <- min_no_na(ylim[[i]])
                    plot(pvalue.posN[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],logpvalue[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],pch=pch,type=type,lwd=cex[2]*(R/2)+1,cex=cex[2]*(R/2),col=rep(rep(colx,N[i]),add[[i]])[logpvalue>=min_no_na(ylim[[i]]) & is_visable[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)+band),ylim=ylim[[i]],ann=FALSE,
                        cex.axis=cex.axis*(R/2),font=1,axes=FALSE,yaxs="r")
                    mtext(side = 2, text = ylab, line = ylab.pos, cex=cex.lab*(R/2), font=1, xpd=TRUE)
                }

                if(chr.border){
                    for(b in 1:length(chr.border.pos)){
                        segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=lwd.axis, lty=2)
                    }
                }
                # Max1 <- Max
                # Min1 <- Min
                # if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
                # if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
                
                #add the names of traits on plot 
                if(!is.null(threshold)){
                    if(LOG10){
                        threshold.max <- -log10(min(threshold))
                    }else{
                        threshold.max <- max(threshold)
                    }
                    if(threshold.max == Max){
                        text(max_no_na(pvalue.posN),Max*0.98,labels=trait[i],adj=c(1, 1),font=4,cex=cex.lab*(R/2),xpd=TRUE) 
                    }else if((threshold.max + 0.5 * strheight(trait[i], cex=cex.lab*(R/2))) >= Max){
                        text(max_no_na(pvalue.posN),threshold.max*1.02,labels=trait[i],adj=c(1, 0),font=4,cex=cex.lab*(R/2),xpd=TRUE) 
                    }else{
                        text(max_no_na(pvalue.posN),Max,labels=trait[i],adj=1,font=4,cex=cex.lab*(R/2),xpd=TRUE) 
                    }
                }else{
                   text(max_no_na(pvalue.posN),Max,labels=trait[i],adj=1,font=4,cex=cex.lab*(R/2),xpd=TRUE) 
                }
                
                if(i == R){
                    if(chr.labels.angle == 0){
                        if(is.null(chr.labels)){
                            axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2,labels=c("Chr",chr.ori),padj=1)
                        }else{
                            if(Nchr == 1){
                                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2), cex.axis=cex.axis*(R/2),font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                            }else{
                                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2), cex.axis=cex.axis*(R/2),font=2,labels=c("Chr",chr.labels))
                            }
                        }
                    }else{
                        axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2),labels=FALSE)
                        if(is.null(chr.labels)){
                            text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis*(R/2), font = 2, labels = c("Chr",chr.ori), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %% 90 == 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                            # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis,cex.axis=cex.axis*(R/2),font=2,labels=c("Chr",chr.ori),padj=1)
                        }else{
                            if(Nchr == 1){
                                # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2), cex.axis=cex.axis*(R/2),font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                                text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis*(R/2), font = 2, labels = c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %% 90 == 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                            }else{
                                # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis*(R/2), cex.axis=cex.axis*(R/2),font=2,labels=c("Chr",chr.labels))
                                text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis*(R/2), font = 2, labels = c("Chr",chr.labels), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %% 90 == 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                            }
                        }
                    }
                    axis(1, mgp=c(3,xticks.pos,0), at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
                }
                #if(i==1) mtext("Manhattan plot",side=3,padj=-1,font=2,cex=xn)
                if(is.null(ylim)){
                    if((Max-Min)>1){
                        axis(2, las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
                        axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
                    }else{
                        axis(2,las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
                        axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
                    }
                }else{
                    axis(2, las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
                    axis(2, at=c((Min), Max), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))

                    # if(ylim[2]>1){
                    #     axis(2,las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
                    #     axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
                    # }else{
                    #     axis(2,las=1,lwd=lwd.axis*(R/2),cex.axis=cex.axis*(R/2),font=2)
                    #     axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis*(R/2))
                    # }
                }
                if(!is.null(threshold)){
                    for(thr in 1:length(threshold)){
                        h <- ifelse(LOG10, -log10(threshold[thr]), threshold[thr])
                        segments(0, h, max_no_na(pvalue.posN), h, col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                    }
                    if(amplify==TRUE){
                        if(LOG10){
                            threshold <- sort(threshold)
                            sgline1=-log10(max_no_na(threshold))
                        }else{
                            threshold <- sort(threshold, decreasing=TRUE)
                            sgline1=min_no_na(threshold)
                        }
                        sgindex=which(logpvalue>=sgline1)
                        HY1=logpvalue[sgindex]
                        HX1=pvalue.posN[sgindex]
                        
                        #cover the points that exceed the threshold with the color "white"
                        points(HX1,HY1,pch=pch,cex=cex[2]*R,col="white")
                        
                        for(ll in 1:length(threshold)){
                            if(ll == 1){
                                if(LOG10){
                                    sgline1=-log10(threshold[ll])
                                }else{
                                    sgline1=threshold[ll]
                                }
                                sgindex=which(logpvalue>=sgline1)
                                HY1=logpvalue[sgindex]
                                HX1=pvalue.posN[sgindex]
                            }else{
                                if(LOG10){
                                    sgline0=-log10(threshold[ll-1])
                                    sgline1=-log10(threshold[ll])
                                }else{
                                    sgline0=threshold[ll-1]
                                    sgline1=threshold[ll]
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
                    points(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],pch=pch,cex=cex[2]*R,col="white")
                    if(!is.na(highlight_index[[i]][1])){
                        if(is.null(highlight.col)){
                            highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),xadj=highlight.text.xadj[[i]],yadj=highlight.text.yadj[[i]],words=highlight.text[[i]],point.cex=highlight.cex*R,text.cex=highlight.text.cex*R/2, pch=highlight.pch,type=highlight.type,point.col=rep(rep(colx,N[i]),add[[i]])[highlight_index[[i]]],text.col=highlight.text.col,text.font=highlight.text.font)
                        }else{
                            highlight_text(x=pvalue.posN[highlight_index[[i]]],y=logpvalue[highlight_index[[i]]],xlim=c(min_no_na(pvalue.posN)-band,max_no_na(pvalue.posN)),ylim=c(Min,Max),xadj=highlight.text.xadj[[i]],yadj=highlight.text.yadj[[i]],words=highlight.text[[i]],point.cex=highlight.cex*R,text.cex=highlight.text.cex*R/2, pch=highlight.pch,type=highlight.type,point.col=highlight.col,text.col=highlight.text.col,text.font=highlight.text.font)
                        }
                    }
                }

                #if(!is.null(threshold) & !is.null(signal.line))    abline(v=pvalue.posN[which(pvalueT[,i] < min_no_na(threshold))],col="grey",lty=2,lwd=signal.line)
            }
            
            #add the labels of X-axis
            #mtext(xlab,side=1,padj=2.5,font=2,cex=R*2/3)
            if(file.output) dev.off()
            
            if(file.output){
                ht=ifelse(is.null(height), 6, height)
                wh=ifelse(is.null(width), 14, width)
                if(file=="jpg") jpeg(paste("Multraits-Manhattan.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
                if(file=="pdf") pdf(paste("Multraits-Manhattan.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
                if(file=="tiff")    tiff(paste("Multraits-Manhattan.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
                if(!is.null(trait.legend.ncol)){
                    mar[3] = mar[3] + ceiling(length(taxa) / trait.legend.ncol)
                }
                par(mar = mar,xaxs="i",yaxs="r")
            }
            if(!file.output){
                ht=ifelse(is.null(height), 6, height)
                wh=ifelse(is.null(width), 14, width)
                if(is.null(dev.list())) dev.new(width = wh, height = ht)
                # par(xpd=TRUE)
            }
            
            pvalue <- as.vector(Pmap[,3:(R+2)])
            if(is.null(ylim)){
                if(!is.null(threshold)){
                    # if(sum(threshold!=0)==length(threshold)){
                        if(LOG10){
                            Max=max_ylim(max_no_na(c((-log10(min_no_na(pvalue))),-log10(min_no_na(threshold)))))
                            Min<-min_ylim(min_no_na(c((-log10(max_no_na(pvalue))),-log10(max_no_na(threshold)))))
                        }else{
                            Max=max_ylim(max_no_na(c((max_no_na(pvalue)),max_no_na(threshold))))
                            # if(abs(Max)<=1)   Max=max_no_na(c(max_no_na(pvalue),max_no_na(threshold)))
                            Min <- min_ylim(min_no_na(c((min_no_na(pvalue)),min_no_na(threshold))))
                            # if(abs(Min)<=1)   Min=min_no_na(c(min_no_na(pvalue),min_no_na(threshold)))
                        }
                    # }else{
                    #     if(LOG10){
                    #         Max=max_ylim((-log10(min_no_na(pvalue))))
                    #         Min=0
                    #     }else{
                    #         Max=max_ylim((max_no_na(pvalue)))
                    #         # if(abs(Max)<=1)   Max=max_no_na(max_no_na(pvalue))
                    #         Min<- min_ylim((min_no_na(pvalue)))
                    #         # if(abs(Min)<=1)   Min=min_no_na(min_no_na(pvalue))
                    #         # }else{
                    #             # Max=max_no_na(ceiling(max_no_na(pvalue)))
                    #         # }
                    #     }   
                    # }
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
                            cex.axis=cex.axis,font=2,axes=FALSE)
                    }else{
                        plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                            cex.axis=cex.axis,font=2,axes=FALSE)
                    }
                }else{
                    if(cir.density){
                        plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-(Max-Min)/den.fold,Max),ann=FALSE,
                            cex.axis=cex.axis,font=2,axes=FALSE)
                    }else{
                        plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min,Max),ann=FALSE,
                            cex.axis=cex.axis,font=2,axes=FALSE)
                    }
                }
                mtext(side = 2, text = ylab, line = ylab.pos, cex=cex.lab, font=1, xpd=TRUE)
            }else{
                Max <- max_no_na(unlist(ylim))
                Min <- min_no_na(unlist(ylim))
                if(cir.density){
                    plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+1.05*max_no_na(pvalue.posN)),ylim=c(Min-Max/den.fold,Max),ann=FALSE,
                        cex.axis=cex.axis,font=2,axes=FALSE)
                }else{
                    plot(NULL,xlim=c(min_no_na(pvalue.posN)-band,band+max_no_na(pvalue.posN)),ylim=c(Min, Max),ann=FALSE,
                        cex.axis=cex.axis,font=2,axes=FALSE)
                }
                mtext(side = 2, text = ylab, line = ylab.pos, cex=cex.lab, font=1, xpd=TRUE)
            }

            # Max1 <- Max
            # Min1 <- Min
            # if(abs(Max) <= 1) Max <- round(Max, ceiling(-log10(abs(Max))))
            # if(abs(Min) <= 1) Min <- round(Min, ceiling(-log10(abs(Min))))
            if(!is.null(ylim)){
                if(is.null(trait.legend.ncol)){
                    legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,Max,trait,col=t(col)[1:R],pch=pch,text.font=6,cex=cex.lab,box.col=NA,horiz=TRUE,xjust=0.5,yjust=0,xpd=TRUE)
                }else{
                    legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,Max,trait,col=t(col)[1:R],pch=pch,text.font=6,cex=cex.lab,box.col=NA,horiz=FALSE,ncol=trait.legend.ncol,xjust=0.5,yjust=0,xpd=TRUE)
                }
            }else{
                if(is.null(trait.legend.ncol)){
                    legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,Max,trait,col=t(col)[1:R],pch=pch,text.font=6,cex=cex.lab,box.col=NA,horiz=TRUE,xjust=0.5,yjust=0,xpd=TRUE)
                }else{
                    legend((max_no_na(pvalue.posN)+min_no_na(pvalue.posN))*0.5,Max,trait,col=t(col)[1:R],pch=pch,text.font=6,cex=cex.lab,box.col=NA,horiz=FALSE,ncol=trait.legend.ncol,xjust=0.5,yjust=0,xpd=TRUE)
                }
            }
            if(chr.labels.angle == 0){
                if(is.null(chr.labels)){
                    axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks),lwd=lwd.axis,cex.axis=cex.axis,font=2,labels=c("Chr",chr.ori))
                }else{
                    if(Nchr == 1){
                        axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                    }else{
                        axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.labels))
                    }
                }
            }else{
                axis(1, mgp=c(3,xticks.pos,0), at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis,labels=FALSE)
                if(is.null(chr.labels)){
                    text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis, font = 2, labels = c("Chr",chr.ori), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %% 90== 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                    # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks),lwd=lwd.axis,cex.axis=cex.axis,font=2,labels=c("Chr",chr.ori)) 
                }else{
                    if(Nchr == 1){
                        text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis, font = 2, labels = c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %% 90== 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                        # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c(paste("Chr.", unique(Pmap[,1]), bp_lab, sep=""),chr.labels))
                    }else{
                        text(c(min_no_na(pvalue.posN)-band,ticks), par("usr")[3]*2-ifelse(cir.density, Min-(Max-Min)/den.fold, Min), cex=cex.axis, font = 2, labels = c("Chr",chr.labels), srt = chr.labels.angle, xpd = TRUE,adj=c(ifelse(chr.labels.angle %%90 == 0, 0.5, ifelse(chr.labels.angle > 0, 1, 0)), 1))
                        # axis(1, at=c(min_no_na(pvalue.posN)-band,ticks), lwd=lwd.axis, cex.axis=cex.axis,font=2,labels=c("Chr",chr.labels))
                    }
                }
            }
            axis(1, mgp=c(3,xticks.pos,0), at=c(ticks[length(ticks)], max_no_na(pvalue.posN)), labels=c("",""), tcl=0, lwd=lwd.axis)
            if(is.null(ylim)){
                if((Max-Min)>1){
                    #print(seq(0,(Max+1),ceiling((Max+1)/10)))
                    axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
                    axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
                    legend.y <- Max
                }else{
                    axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
                    axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
                    legend.y <- Max
                }
            }else{
                axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
                axis(2, at=c(Min, Max), labels=c("",""), tcl=0, lwd=lwd.axis)
                legend.y <- Max
                # if(ylim[2]>1){
                #     axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
                #     axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis)
                #     legend.y <- tail(ylim[2], 1)
                # }else{
                #     axis(2,las=1,lwd=lwd.axis,cex.axis=cex.axis,font=2)
                #     axis(2, at=c(min_no_na(ylim), ylim[2]), labels=c("",""), tcl=0, lwd=lwd.axis)
                #     legend.y <- tail(ylim[2], 1)
                # }
            }
            if(chr.border){
                for(b in 1:length(chr.border.pos)){
                    segments(chr.border.pos[b], Min, chr.border.pos[b], Max, col="grey45", lwd=lwd.axis, lty=2)
                }
            }
            do <- TRUE
            sam.index <- list()
            trait_max_n <- 0
            trait_max <- 0
            for(l in 1:R){
                sam.index[[l]] <- c(1:nrow(Pmap))[is_visable[[l]]]
                if(length(sam.index[[l]]) >= trait_max_n){
                    trait_max_n = length(sam.index[[l]])
                    trait_max = l
                }
            }
            
            #change the sample number according to Pmap
            #sam.num <- ceiling(nrow(Pmap)/100)
            sam.num <- 1000
            cat_bar <- seq(1, 100, 1)
            while(do){
                for(i in 1:R){
                    if(length(sam.index[[i]]) == 0){
                        # nothing
                    }else{
                        if(length(sam.index[[i]]) < sam.num){
                            plot.index <- sam.index[[i]]
                        }else{
                            plot.index <- sample(sam.index[[i]], sam.num, replace=FALSE)
                        }
                        sam.index[[i]] <- sam.index[[i]][-which(sam.index[[i]] %in% plot.index)]
                        logpvalue=logpvalueT[plot.index,i]
                        if(!is.null(ylim)){indexx <- logpvalue>=min_no_na(ylim[[i]])}else{indexx <- 1:length(logpvalue)}
                        points(pvalue.posN[plot.index][indexx],logpvalue[indexx],pch=pch[i],type=type,lwd=cex[2]+1,cex=cex[2],col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255))
                        #if(!is.null(threshold) & (length(grep("FarmCPU",taxa[i])) != 0))   abline(v=which(pvalueT[,i] < min_no_na(threshold)/max_no_na(dim(Pmap))),col="grey",lty=2,lwd=signal.line)
                    }
                }
                if(verbose){
                    progress <- round((nrow(Pmap) - length(sam.index[[trait_max]])) * 100 / nrow(Pmap))
                    if(progress %in% cat_bar){
                        cat(" Multraits-Rectangular Plotting...(finished ", progress, "%)\r", sep="")
                        cat_bar <- cat_bar[cat_bar != progress]
                        if(progress == 100) cat("\n")
                    }
                }
                if(length(sam.index[[trait_max]]) == 0) do <- FALSE
            }

            if(!is.null(threshold)){
                # if(sum(threshold!=0)==length(threshold)){
                    for(i in 1:R){
                        logpvalue=logpvalueT[, i]
                        if(LOG10){
                            sgindex = which(logpvalue > -log10(min(threshold)))
                        }else{
                            sgindex = which(logpvalue > max(threshold))
                        }
                        HY1=logpvalue[sgindex]
                        HX1=pvalue.posN[sgindex]
                        points(HX1,HY1,pch=pch[i],cex=cex[2],col="white")
                        if(!is.null(signal.col)){
                            points(HX1,HY1,pch=rep(signal.pch, R)[i],cex=rep(signal.cex, R)[i],col=rep(signal.col, R)[i])
                        }else{
                            points(HX1,HY1,pch=rep(signal.pch, R)[i],cex=rep(signal.cex, R)[i],col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255))
                        }
                    }
                    for(thr in 1:length(threshold)){
                        h <- ifelse(LOG10, -log10(threshold[thr]), threshold[thr])
                        segments(0, h, max_no_na(pvalue.posN), h,col=threshold.col[thr],lwd=threshold.lwd[thr],lty=threshold.lty[thr])
                    }
                # }
            }

            if(is.null(ylim)){ymin <- Min}else{ymin <- min_no_na(unlist(ylim))}
            if(cir.density){
                for(yll in 1:length(pvalue.posN.list)){
                    polygon(c(min_no_na(pvalue.posN.list[[yll]]), min_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]]), max_no_na(pvalue.posN.list[[yll]])), 
                        c(ymin-0.5*(Max-Min)/den.fold, ymin-1.5*(Max-Min)/den.fold, 
                        ymin-1.5*(Max-Min)/den.fold, ymin-0.5*(Max-Min)/den.fold), 
                        col="grey", border="grey")
                }
                is_visable_den <- filter.points(pvalue.posN, ymin-0.5*(Max-Min)/den.fold, wh, ht, dpi = dpi)
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
                    title="", legend=density.list$legend.y, pch=15, pt.cex = 2.5, col=density.list$legend.col,
                    cex=0.8, bty="n",
                    y.intersp=1,
                    x.intersp=1,
                    yjust=1, xjust=0, xpd=TRUE
                )          
            }
            if(file.output) dev.off()
            
        }
    }
        
    if("q" %in% plot.type){

        signal.col <- rep(signal.col,R)
        signal.pch <- rep(signal.pch,R)
        signal.cex <- rep(signal.cex*1.1,R)

        #print("Starting QQ-plot!",quote=F)
        if(multracks){
            if(file.output){
                ht=ifelse(is.null(height), 5.5, height)
                wh=ifelse(is.null(width), 3.5, width)
                if(file=="jpg") jpeg(paste("Multracks-QQplot.",paste(taxa,collapse="."),".jpg",sep=""), width = R*wh*dpi,height=ht*dpi,res=dpi,quality = 100)
                if(file=="pdf") pdf(paste("Multracks-QQplot.",paste(taxa,collapse="."),".pdf",sep=""), width = R*wh,height=ht)
                if(file=="tiff")    tiff(paste("Multracks-QQplot.",paste(taxa,collapse="."),".tiff",sep=""), width = R*wh*dpi,height=ht*dpi,res=dpi)
                par(mfcol=c(1,R),xpd=TRUE)
            }else{
                ht=ifelse(is.null(height), 5.5, height)
                wh=ifelse(is.null(width), 3.5, width)
                if(is.null(dev.list())) dev.new(width = wh*R, height = ht)
                par(xpd=TRUE)
            }
            log.Quantiles.max_no_na <- NULL
            for(i in 1:R){
                if(i == 1)  par(mar=c(mar[2], mar[2], mar[3], 0))
                if(i == R)  par(mar=c(mar[2], 1.5, mar[3], mar[4]))
                if(i != 1 & i != R) par(mar=c(mar[2], 1.5, mar[3], 0))
                if(verbose) cat(paste(" Multracks-QQ Plotting ",trait[i],".\n",sep=""))        
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
                log.Quantiles.max_no_na <- c(log.Quantiles.max_no_na, max_no_na(log.Quantiles))
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
                    plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,YlimMax),xlab ="", ylab="")
                }else{
                    plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,max(ylim[[i]])),xlab ="", ylab="")
                }
                axis(1, mgp=c(3,xticks.pos,0), at=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), lwd=lwd.axis,labels=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), cex.axis=cex.axis)
                axis(2, las=1, lwd=lwd.axis,cex.axis=cex.axis)
                axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(ylim[[i]]))), labels=c("",""), tcl=0, lwd=lwd.axis)
                
                #plot the confidence interval of QQ-plot
                if(conf.int){
                    if(is.null(conf.int.col)){
                        polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255),border=t(col)[i])
                    }else{
                        polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col[i],border=conf.int.col[i])
                    }
                }
                if(!is.null(threshold.col)){par(xpd=FALSE); abline(a = 0, b = 1,lwd=threshold.lty[1], lty=threshold.lty[1], col = threshold.col[1]); par(xpd=TRUE)}
                is_visable <- filter.points(log.Quantiles, log.P.values, wh, ht, dpi = dpi)
                if(!is.null(threshold)){
                    # if(sum(threshold!=0)==length(threshold)){
                        thre.line=-log10(min_no_na(threshold))
                        if(amplify==TRUE){
                            thre.index <- log.P.values<thre.line
                            if(sum(!thre.index)!=0){
                                points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col = t(col)[i],pch=19,cex=cex[3])
                            
                                #cover the points that exceed the threshold with the color "white"
                                # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                if(is.null(signal.col)){
                                    points(log.Quantiles[!thre.index],log.P.values[!thre.index],col = t(col)[i],pch=signal.pch[i],cex=signal.cex[i])
                                }else{
                                    points(log.Quantiles[!thre.index],log.P.values[!thre.index],col = signal.col[i],pch=signal.pch[i],cex=signal.cex[i])
                                }
                            }else{
                                points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                            }
                        }else{
                            points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                        }
                    # }
                }else{
                    points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                }
                if(all(main != "")) {
                    title(main = main[i], cex.main = main.cex, font.main= main.font)
                }else{
                    title(main = trait[i], cex.main = main.cex, font.main= main.font) 
                }
                mtext(side = 1, text = expression(Expected~~-log[10](italic(p))), line = ylab.pos+2, cex=cex.lab, font=1, xpd=TRUE)
                if(i == 1)  mtext(side = 2, text = expression(Observed~~-log[10](italic(p))), line = ylab.pos, cex=cex.lab, font=1, xpd=TRUE)
            }
            if(box) box(lwd=lwd.axis)
            if(file.output) dev.off()
            if(R > 1){
                signal.col <- NULL
                if(file.output){
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    if(file=="jpg") jpeg(paste("Multraits-QQplot.",paste(taxa,collapse="."),".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
                    if(file=="pdf") pdf(paste("Multraits-QQplot.",paste(taxa,collapse="."),".pdf",sep=""), width = wh,height=ht)
                    if(file=="tiff")    tiff(paste("Multraits-QQplot.",paste(taxa,collapse="."),".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
                    par(mar = c(mar[2],mar[2],mar[3],mar[4]),xpd=TRUE)
                }else{  
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    dev.new(width = wh, height = ht)
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
                    plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles.max_no_na)+1)), axes=FALSE, xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0, floor(YlimMax+1)), main = "QQplot", cex.main=main.cex, font.main=main.font)
                }else{
                    plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles.max_no_na)+1)), axes=FALSE, xlab="", ylab="", cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0, max(unlist(ylim))),main = "QQplot", cex.main=main.cex, font.main=main.font)
                }
                legend("topleft",trait,col=t(col)[1:R],pch=19,text.font=6,box.col=NA, xpd=TRUE)
                axis(1, mgp=c(3,xticks.pos,0), at=seq(0,floor(max_no_na(log.Quantiles.max_no_na)+1),ceiling((max_no_na(log.Quantiles.max_no_na)+1)/10)), lwd=lwd.axis,labels=seq(0,floor(max_no_na(log.Quantiles.max_no_na)+1),ceiling((max_no_na(log.Quantiles.max_no_na)+1)/10)), cex.axis=cex.axis)
                axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis)
                axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(unlist(ylim)))), labels=c("",""), tcl=0, lwd=lwd.axis)

                mtext(side = 1, text = expression(Expected~~-log[10](italic(p))), line = ylab.pos+1, cex=cex.lab, font=1, xpd=TRUE)
                mtext(side = 2, text = expression(Observed~~-log[10](italic(p))), line = ylab.pos, cex=cex.lab, font=1, xpd=TRUE)
                
                for(i in 1:R){
                    if(verbose) cat(paste(" Multraits-QQ Plotting ",trait[i],".\n",sep=""))
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
                            polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255),border=t(col)[i])
                        }else{
                            polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col[i],border=conf.int.col[i])
                        }
                    }
                       
                    if((i == R) & !is.null(threshold.col)){par(xpd=FALSE); abline(a = 0, b = 1,lwd=threshold.lty[1], lty=threshold.lty[1], col = threshold.col[1]); par(xpd=TRUE)}
                    # points(log.Quantiles, log.P.values, col = t(col)[i],pch=19,cex=cex[3])
                    is_visable <- filter.points(log.Quantiles, log.P.values, wh, ht, dpi = dpi)
                    if(!is.null(threshold)){
                        # if(sum(threshold!=0)==length(threshold)){
                            thre.line=-log10(min_no_na(threshold))
                            if(amplify==TRUE){
                                thre.index <- log.P.values<thre.line
                                if(sum(!thre.index)!=0){
                                    points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col = t(col)[i],pch=19,cex=cex[3])
                            
                                    # cover the points that exceed the threshold with the color "white"
                                    # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                    if(is.null(signal.col)){
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col = t(col)[i],pch=signal.pch[i],cex=signal.cex[i])
                                    }else{
                                        points(log.Quantiles[!thre.index],log.P.values[!thre.index],col = signal.col[i],pch=signal.pch[i],cex=signal.cex[i])
                                    }
                                }else{
                                    points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                                }
                            }else{
                                points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                            }
                        # }
                    }else{
                        points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                    }
                }
                if(box) box(lwd=lwd.axis)
                if(file.output) dev.off()
            }
        }else{
            for(i in 1:R){
                if(verbose) cat(paste(" QQ Plotting ",trait[i],".\n",sep=""))
                if(file.output){
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    if(file=="jpg") jpeg(paste("QQplot.",taxa[i],".jpg",sep=""), width = wh*dpi,height=ht*dpi,res=dpi,quality = 100)
                    if(file=="pdf") pdf(paste("QQplot.",taxa[i],".pdf",sep=""), width = wh,height=ht)
                    if(file=="tiff")    tiff(paste("QQplot.",taxa[i],".tiff",sep=""), width = wh*dpi,height=ht*dpi,res=dpi)
                     par(mar = c(mar[2],mar[2],mar[3],mar[4]),xpd=TRUE)
                }else{
                    ht=ifelse(is.null(height), 5.5, height)
                    wh=ifelse(is.null(width), 5.5, width)
                    if(is.null(dev.list())) dev.new(width = wh, height = ht)
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
                    plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,YlimMax),xlab="",ylab="")
                }else{
                    plot(NULL, xlim = c(0,floor(max_no_na(log.Quantiles)+1)), axes=FALSE, cex.axis=cex.axis, cex.lab=cex.lab,ylim=c(0,max(ylim[[i]])),xlab="",ylab="")      
                }
                axis(1, mgp=c(3,xticks.pos,0),at=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), lwd=lwd.axis,labels=seq(0,floor(max_no_na(log.Quantiles)+1),ceiling((max_no_na(log.Quantiles)+1)/10)), cex.axis=cex.axis)
                axis(2, las=1,lwd=lwd.axis,cex.axis=cex.axis)
                axis(2, at=c(0, ifelse(is.null(ylim), YlimMax, max(ylim[[i]]))), labels=c("",""), tcl=0, lwd=lwd.axis)

                mtext(side = 1, text = expression(Expected~~-log[10](italic(p))), line = ylab.pos+1, cex=cex.lab, font=1, xpd=TRUE)
                mtext(side = 2, text = expression(Observed~~-log[10](italic(p))), line = ylab.pos, cex=cex.lab, font=1, xpd=TRUE)
                
                #plot the confidence interval of QQ-plot
                if(conf.int){
                    if(is.null(conf.int.col)){
                        polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=rgb(col2rgb(t(col)[i])[1], col2rgb(t(col)[i])[2], col2rgb(t(col)[i])[3], 100, maxColorValue=255),border=t(col)[i])
                    }else{
                        polygon(c(log.Quantiles[index],log.Quantiles),c(-log10(c05)[index],-log10(c95)),col=conf.int.col[i],border=conf.int.col[i])
                    }
                }

                if(!is.null(threshold.col)){par(xpd=FALSE); abline(a = 0, b = 1,lwd=threshold.lty[1], lty=threshold.lty[1], col = threshold.col[1]); par(xpd=TRUE)}
                # points(log.Quantiles, log.P.values, col = t(col)[i],pch=19,cex=cex[3])
                is_visable <- filter.points(log.Quantiles, log.P.values, wh, ht, dpi = dpi)
                if(!is.null(threshold)){
                    # if(sum(threshold!=0)==length(threshold)){
                        thre.line=-log10(min_no_na(threshold))
                        if(amplify==TRUE){
                            thre.index <- log.P.values<thre.line
                            if(sum(!thre.index)!=0){
                                points(log.Quantiles[thre.index & is_visable], log.P.values[thre.index & is_visable], col = t(col)[i],pch=19,cex=cex[3])
                            
                                #cover the points that exceed the threshold with the color "white"
                                # points(log.Quantiles[thre.index],log.P.values[thre.index], col = "white",pch=19,cex=cex[3])
                                # print(signal.col)
                                # print(signal.pch)
                                # print(signal.cex)
                                if(is.null(signal.col)){
                                    points(log.Quantiles[!thre.index],log.P.values[!thre.index],col = t(col)[i],pch=signal.pch[i],cex=signal.cex[i])
                                }else{
                                    points(log.Quantiles[!thre.index],log.P.values[!thre.index],col = signal.col[i],pch=signal.pch[i],cex=signal.cex[i])
                                }
                            }else{
                                points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                            }
                        }else{
                            points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                        }
                    # }
                }else{
                    points(log.Quantiles[is_visable], log.P.values[is_visable], col = t(col)[i],pch=19,cex=cex[3])
                }
                if(all(main != "")) {
                    title(main = main[i], cex.main = main.cex, font.main= main.font)
                }else{
                    title(main = trait[i], cex.main = main.cex, font.main= main.font) 
                }
                if(box) box(lwd=lwd.axis)
                if(file.output) dev.off()
            }
        }
    }
    if(file.output & verbose)   cat(paste(" Plots are stored in: ", getwd(), sep=""), "\n")
}
