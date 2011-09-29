
library(rgl)  # load 3d library
library(RColorBrewer) # color library
source("/home/DS/phyjpb/Programming/R/fullrot.R")

topodata = read.table (file="filename_topologic", header=FALSE, sep='\t', quote='"\'', dec='.', fill=FALSE, comment.char="#",  na.strings = "NA", nrows = -1, skip = 5, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE)
#load data for height and phase

phasedata = read.table (file="filename_phase", header=FALSE, sep='\t', quote='"\'', dec='.', fill=FALSE, comment.char="#",  na.strings = "NA", nrows = -1, skip = 5, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE) 


P_002 = data.matrix(phasedata)	# convert tables into matricies
H_002 = data.matrix(topodata)

Hp = H_002 + abs(min(H_002))	# upshift phase data
Hpmax = max(Hp)

#x<-seq(0,20,length=512)	# define lengths of image, should match the dimensions of the scans
#y<-seq(0,20,length=512)

x<-seq(0,25,length=512)	# define lengths of image, should match the dimensions of the scans
y<-seq(0,25,length=512)

Pn_002=P_002/abs(range(P_002)[1])	#normalisation of the data

n=10					#the number of colors
#colorlut = heat.colors(n)		#get the colors from the heat palette
colorlut = brewer.pal(n,"RdYlBu")
col = colorlut[(Pn_002+1)*((n-1)*0.5)+1]  # convert the data to the nearest integer in 1 n range

Hf= 1.3  				# scaling factor

open3d(windowRect=c(0,0,1024,768),userMatrix=fullrot(pi/6,-pi/3,pi/12)) # opens a big window
surface3d(x,y,Hf*(Hp/Hpmax),color=col) 
#zlab = c(0,500,1000)		# true scale
zlab = c(0,200,400)
zat = Hf*zlab/Hpmax
yat = c(0,5,10,15,20,25)

bbox3d(color=c("white","black"), shininess=1, alpha=0.2,zat=zat,zlab=zlab,yat=yat,ylab=yat,expand = 1.1 ) # bounding box, lots of options
par3d(zoom=0.8638377)

#snapshot3d("S1_25x25_0p5.png")		# save to file

png("Real_H_004_legend.png")
phaselim = round(range(P_002)[1],1)
X2 = round(seq(phaselim,-phaselim,length.out=n),2)
#nf=layout(matrix(c(3,1,2,4,5),1,5,byrow=TRUE))
par(pin=c(0.8,6))
barplot(rep(1,n),horiz=TRUE, col=colorlut,width=0.5, axes = FALSE, names.arg = X2, las=2, ylab="Phase (°)",font=1,cex.lab=2)
dev.off()

