##package to install for the code 
library(ramify)
library(Hmisc)
library(StratigrapheR)
#Formula of the fundamental sampling error
# S_Fe = (1/*Ml - 1/Ms) * c * f * g * l * d**3
# Mass samples Ms and the nominal size d
Ms <- c(41,6,6,1.5,1.5,0.1,0.1)
d <- c(6,6,1,1,0.02,0.02,0.075)
# Parameters of the FSE
Lm= 3.2 #density of apatite
Lg=2.7 #density of gangue
pct_ap = 0.4 #% of apatite in pure mineral 
aL=(pct_ap/100)/0.183 ##Pct P to apatite

f=0.5 # # form  factor 
g=0.25 ## granulometric distribution factor
Lgs=0.08 ## grain size of mineral in mm
Dgs =60 ## diameter of the largest fragmenet in mm   
l=(Lgs/Dgs)**0.5 ## Liberation parameter
c = ((1-aL)/aL) * ((1-aL)*Lm+aL*Lg) #The mineralogical factor, c
#Relative variation σ²[EF] by step
Sfe2 = (c) * f * g * l * d**3 / (Ms*1000)
Sfe2
#FSE σ[EF] by step
Sfe = sqrt(Sfe2)
Sfe
#Relative variances at stage of subsampling A to B
S_Fe_1 = (1/(1000*Ms[2]) - 1/(1000*Ms[1])) * (c) * f * g * l * (d[3]**3)
print(paste("1_" ,"Relative variances at stage of subsampling B to C ",round(S_Fe_1,10)))
#Relative variances at stage of subsampling C to D
S_Fe_2 = (1/(1000*Ms[4]) - 1/(1000*Ms[3])) * (c) * f * g * l * d[5]**3
print(paste("2_" ,"Relative variances at stage of subsampling D to E ",round(S_Fe_2,10)))
#Relative variances at stage of subsampling E to F
S_Fe_3 = (1/(1000*Ms[6]) - 1/(1000*Ms[5])) * (c) * f * g * l * d[6]**3
print(paste("2_" ,"Relative variances at stage of subsampling D to E ",round(S_Fe_3,10)))

#Relative variances at stage of subsampling A to F
S_FE = S_Fe_1+S_Fe_2+S_Fe_3
print(paste("3_S_Fe_1+S_Fe_2"))
print(round(S_FE,5))

#FSE at stage of subsampling A to E
print(paste("4_FSE"))
print(round(sqrt(S_FE)*100,5))
print('')


# log horizontal lines 
lseq <- function(from=10**-8, to=10**-2, length.out=9) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}
#labeling the X axis using log
ticks <- c(0.01,0.1,1,10,100,1000)
labels <- c("A", "B", "C", "D", "E", "F","G")
log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceil(lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}
## Plot the Relative variation σ²[EF] with Gy 10% safly line
plot(Ms,Sfe2, type = "b", pch = 21, cex = 3, lwd = 3, log = "xy",
     xlab = "Mass [kg]", ylab = bquote(paste('Relative variance '*'  σ'[EF]^2*'')),
     ylim=c(10**-10,10**-1),xlim=c(0.01,100),
     xaxt = 'n')
text(Ms,Sfe2, labels, col = "red", cex = 1.2)
axis(1, at=log10Tck('x','minor'), tcl= 0.1, labels=NA) # bottom
axis(side=1,at=ticks,labels=c('0.01','0.1','1','10','100','1000'))
abline(h = lseq(), v = ticks, col = "lightgray", lty = 5)
abline(h = 0.01,  col = "red", lty = 10, lw=2)
text(x=0.1,y=0.015, labels="Gy 10% safety line ", col = "red", cex=1.5)
text(x=Ms[1]/2,y=Sfe2[1]*3, labels="Sampling ", col = "red", cex=1, srt=-20)
text(x=Ms[2]*1.5,y=Sfe2[3]*8, labels="Crushing to 10mm ", col = "red", cex=0.9, srt=-90)
text(x=Ms[3]/2,y=Sfe2[3]*3, labels="Sampling ", col = "red", cex=1, srt=-20)
text(x=Ms[4]*1.2,y=Sfe2[5]*70, labels="Grinding to 0.2mm ", col = "red", cex=0.9, srt=-90)
text(x=Ms[5]/2,y=Sfe2[5]*3, labels="Sampling ", col = "red", cex=1, srt=-20)
