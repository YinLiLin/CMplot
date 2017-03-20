CMplot v3.2.0
=========

## A high-quality drawing tool designed for genome-wide association study

### Installation

**CMplot** is available on CRAN, so it can be installed with the following R code:

```r
install.packages("CMplot")
library("CMplot")
```

There are two example datasets attached in **CMplot**, users can export and view the details by following R code:

```r
data(pig60K)
data(cattle50K)
head(pig60K)
head(cattle50K)
```

Total 40 parameters are available in **CMplot**, typing ```?CMplot``` can get the detail function of all parameters.

### SNP-density plot

```r
CMplot(pig60K,plot.type="d",col=c("darkgreen", "yellow", "red"))
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/illumilla_60K.jpg">
<img src="Figure/illumilla_60K.jpg" height="460px" width="680px">
</a>
</p>

### Circular-Manhattan plot

```r
CMplot(pig60K,plot.type="c",cir.chr.h=1,chr.labels=paste("Chr",c(1:18,"X"),sep=""),threshold=c(0.05,0.01),
      amplify=T,threshold.lty=c(2,1),threshold.col=c("blue","red"))
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/Circular-Manhattan.jpg">
<img src="Figure/Circular-Manhattan.jpg" height="400px" width="400px">
</a>
</p>

### Rectangular-Manhattan plot

```r
CMplot(pig60K[,c(1:3,6)],plot.type="m")
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/Rectangular-Manhattan.trait3.jpg">
<img src="Figure/Rectangular-Manhattan.trait3.jpg" height="300px" width="900px">
</a>
</p>

### Q-Q plot

```r
CMplot(pig60K[,c(1:3,6)],plot.type="q")
```

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-CMplot/master/Figure/QQplot.trait3.jpg">
<img src="Figure/QQplot.trait3.jpg" height="400px" width="400px">
</a>
</p>

### Contact
Questions, suggestions, and bug reports are welcome and appreciated.
- **Author:** Lilin Yin
- **Contact:** yilin@163.com
- **Institution:** [*Huazhong agriculture university*](http://www.hzau.edu.cn/2014/ch/)
