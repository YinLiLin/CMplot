CMplot v3.2.0 [![](https://img.shields.io/github/issues-raw/badges/shields/website.svg)](https://github.com/YinLiLin/R-CMplot/issues) [![](https://img.shields.io/chrome-web-store/stars/nimelepbpejjlbmoobocpfnjhihnpked.svg)](https://cran.r-project.org/web/packages/CMplot/)
=========

## A high-quality drawing tool designed for genome-wide association study

### Installation

**CMplot** is available on CRAN, so it can be installed with the following R code:

```r
install.packages("CMplot")
library("CMplot")
#(optional)if you want to use the latest version:
#source("https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/CMplot.r")
```

There are two example datasets attached in **CMplot**, users can export and view the details by following R code:

```r
data(pig60K)   #calculated p-values by MLM
data(cattle50K)   #calculated SNP effects by rrblup
head(pig60K)
head(cattle50K)
```

Total 40 parameters are available in **CMplot**, typing ```?CMplot``` can get the detail function of all parameters.

### SNP-density plot

```r
CMplot(pig60K,plot.type="d",col=c("darkgreen", "yellow", "red"),file="jpg",dpi=300)
# users can personally set the windowsize and the max of legend by:
# bin.size=1e6
# bin.max=N
# but the latest version should be sourced.
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/illumilla_60K.jpg">
<img src="Figure/illumilla_60K.jpg" height="460px" width="680px">
</a>
</p>

### Circular-Manhattan plot

#### (1) Genome-wide association study(GWAS)

```r
CMplot(pig60K,plot.type="c",chr.labels=paste("Chr",c(1:18,"X"),sep=""),threshold=c(0.05,0.01),
      cir.chr.h=1,amplify=TRUE,threshold.lty=c(2,1),threshold.col=c("blue","red"),signal.line=1,
      signal.col="red",file="jpg",dpi=300)
#Note: if signal.line=NULL, the lines that crosse circles won't be added.
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/Circular-Manhattan.jpg">
<img src="Figure/Circular-Manhattan.jpg" height="400px" width="400px">
</a>
</p>

#### (2) Genomic Selection/Prediction(GS/GP)

```r
CMplot(cattle50K,plot.type="c",LOG10=FALSE,outward=TRUE,chr.labels=paste("Chr",c(1:29),sep=""),
r=1.2,cir.chr.h=1.3,cir.legend.cex=0.5,cir.band=1,threshold=NULL,file="jpg",dpi=300)
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/Circular-Manhattan.cattle.jpg">
<img src="Figure/Circular-Manhattan.cattle.jpg" height="400px" width="400px">
</a>
</p>

### Rectangular-Manhattan plot

#### (1) Genome-wide association study(GWAS)

```r
CMplot(pig60K[,c(1:3,6)],plot.type="m",threshold=NULL,file="jpg",dpi=300)
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/Rectangular-Manhattan.trait3.jpg">
<img src="Figure/Rectangular-Manhattan.trait3.jpg" height="300px" width="900px">
</a>
</p>

#### (2) Genomic Selection/Prediction(GS/GP)

```r
CMplot(cattle50K,plot.type="m",LOG10=FALSE,ylab="SNP effect",threshold=NULL,file="jpg",dpi=300)
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/Rectangular-Manhattan.Fat percentage.jpg">
<img src="Figure/Rectangular-Manhattan.Fat percentage.jpg" height="300px" width="900px">
</a>
</p>

### Q-Q plot

```r
CMplot(pig60K[,c(1:3,6)],plot.type="q",conf.int=TRUE,conf.int.col="grey",file="jpg",dpi=300)
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/QQplot.trait3.jpg">
<img src="Figure/QQplot.trait3.jpg" height="400px" width="400px">
</a>
</p>

### Contact
Questions, suggestions, and bug reports are welcome and appreciated.
- **Author:** Lilin Yin
- **Contact:** ylilin@163.com
- **Institution:** [*Huazhong agricultural university*](http://www.hzau.edu.cn/2014/ch/)
