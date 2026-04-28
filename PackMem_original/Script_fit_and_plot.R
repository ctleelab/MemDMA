# Script fit and plot with R package
# P. Fuchs, R. Gautier 2016
# Update 2017 P. Fuchs


setwd("/scratch/local/casakurai/packmem/PackMem_original/testing-1800ext/DOPC-DOPS/analysis/seperated-leaflets")

########
# Set constants
########
# prefix name
PREFIX="membrane"
# long prefix name
LONG_PREFIX="seperated-leaflets-5fr-Total_Up"
# set precision for writing packdef cosntants (nb of decimals) in the txt 
# output file
PRECISION=2

# lowest defect area used for the fit
# (we recommand not to touch to this value)
LIMX=15
# lowest probability used for the fit
# (we recommand not to touch to this value)
LIMY=1e-4

###########
# !!! In principle no change below !!!
###########

###################
# a few functions #
###################
# this functions returns a linear fit of x & y
fit_decay = function(x,y, LIMX, LIMY){
# fit with defects above LIMX nm and proba > LIMY
RANGE=which(x >= LIMX)
x=x[RANGE]
y=y[RANGE]
RANGE=which(y >= LIMY)
x=x[RANGE]
y=y[RANGE]
FIT=lsfit(x,log(y))
return (FIT)
}

# this function divides the packdef_data in 3 blocks and returns a vector of 3 decays
block_averaging = function(vect) {
bornes = c(1, (length(vect)/3) * 1:3)
decays = c()
for (i in 1:3) {
  subvect = vect[bornes[i]:bornes[i+1]]
  H=hist(subvect,breaks=seq(0,max(vect))+0.5,plot=F)
  x=H$mids ; y=H$counts/sum(H$counts)
  FIT=fit_decay(x,y, LIMX, LIMY)
  decays = c(decays, abs(1/coef(FIT)[2]))
} # <=== end of for loop
return (decays)
}


################################
# Packdef constant calculation #
################################

# Open a pdf device
pdf(paste("Res_",".pdf",sep=""))
# Initialize a data frame to store packdef constants + errors
packdef_constants = data.frame(matrix(0,nrow=6,ncol=3),row.names=c("PackDef_cst_global","PackDef_cst_block1","PackDef_cst_block2","PackDef_cst_block3","PackDef_cst_all_blocks","error_all_blocks"))
names(packdef_constants) = c("deep","shallow","all")
# print the results also to a text file
sink(paste("Res_",".txt",sep=""))

######
# Now loop over the three default types
######
for (DEFECT in c("deep","shallow","all")) {

# we want to use 4 rainbow colors
palette(rainbow(4))

####
# prepare the plot
####
par(mar=c(5, 4.5, 4, 2) + 0.1) # mar ==>  c(bottom, left, top, right) ; default is c(5, 4, 4, 2) + 0.1
plot(1,1,type="n",xlim=c(0,100),ylim=c(-10,-3.5),ylab="Probability", xlab=expression(Defect~area~(A^2)),axes=F, main = DEFECT, cex.lab=1.6,new=T)
box(lwd=2)
box()
# x axis
axis(1,at=seq(0,200,10),labels=F, cex.lab=1.5, cex.axis=1.5,lwd.ticks=2)
axis(1,at=seq(0,200,20),labels=seq(0,200,by=20), cex.lab=1.5, cex.axis=1.5,lwd.ticks=3)
# y axis (log axis)
LABEL=c(expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),expression(10^0))
axis(2,at=log(10^(seq(-5,0,1))),labels=LABEL, cex.lab=2, cex.axis=1.5,lwd.ticks=3)
LOGAXIS=c()
# add secondary tick marks (with a log axis)
for (i in 5:1) LOGAXIS=c(LOGAXIS,log(10^-i * 1:9))
LOGAXIS=c(LOGAXIS,1)
axis(2,at=LOGAXIS,labels=FALSE, lwd.ticks=1.5)
# add legend
legend(45,-3.4,c(PREFIX),col=c(1:length(palette())),lty=1,lwd=3,cex=1.3)
# add lines
abline(h=log(LIMY),lwd=2) ; abline(v=LIMX,lty=2,lwd=2)

####
# load PackMem data
####
filename=paste(LONG_PREFIX,"_",DEFECT,".txt",sep="")
packdef_data=read.table(filename)[,2] 

####
# Compute PackDef distributions (on the whole set)
####
H = hist(packdef_data,breaks=seq(0,max(packdef_data))+0.5,plot=F)
x = H$mids ; y = H$counts/sum(H$counts)
points(x,log(y),col=1,lwd=2)
FIT = fit_decay(x,y, LIMX, LIMY) # intercept = coef(FIT)[1] ; slope = coef(FIT)[2]  
abline(FIT,col=1,lwd=3)
global_inv_decay = abs(1/coef(FIT)[2])

# print global results
print(paste("Results on","for",DEFECT,"defects"))
print(paste("Total number of defects =",length(packdef_data)))
print(paste("Using all data:", round(global_inv_decay,PRECISION),"A^2"))

# Compute PackDef distributions with block averaging method, and print results
FITS_3blocks=block_averaging(packdef_data)
print(paste("Using block 1:",round(FITS_3blocks[1],PRECISION),"A^2"))
print(paste("Using block 2:",round(FITS_3blocks[2],PRECISION),"A^2"))
print(paste("Using block 3:",round(FITS_3blocks[3],PRECISION),"A^2"))
inv_decay_block = mean(FITS_3blocks)
error_inv_decay_block = sd(FITS_3blocks)
print(paste("Mean +/- sd on 3 blocks:",round(inv_decay_block,PRECISION),"+/-",round(error_inv_decay_block,PRECISION),"A^2"))
print("")

####
# store all the results in the data frame
####
packdef_constants["PackDef_cst_global",DEFECT] = global_inv_decay
packdef_constants["PackDef_cst_block1",DEFECT] = FITS_3blocks[1]
packdef_constants["PackDef_cst_block2",DEFECT] = FITS_3blocks[2]
packdef_constants["PackDef_cst_block3",DEFECT] = FITS_3blocks[3]
packdef_constants["PackDef_cst_all_blocks",DEFECT] = inv_decay_block
packdef_constants["error_all_blocks",DEFECT] = error_inv_decay_block

} # <=== END of loop over DEFECT

####
# INTERMEDIATE plot: plot all packdef constants on a single barplot
#
# Allows to estimate the relative convergence of the simulation
#     In this plot, for each packing defects we have 5 bars: 
#     1) the 1st is the cst computed on the whole traj, 
#     2) to 4) the traj is divided in 3 blocks and the cst is calculated on each block, 
#     5) this is the mean +/- error over the three blocks
#
# PLEASE CHECK:
# The value from the 3 blocks should not be too different
# The packdef constant computed on all data (black bar) should be not too far 
# from the mean over the 3 blocks (red bar)
####
par(mar=c(7, 4.6, 4, 2) + 0.1) 
# store the five rows with packdef cst in matrix (to allow the use of barplot)
M = as.matrix(packdef_constants[1:5,])
bp = barplot(M,beside=T,ylab=expression(paste("defect size constant (",A^2,")")),cex.axis=1.4,cex.lab=1.4,col=c("black",gray(1/3),gray(2/3),gray(3/3),"red"),main= "Packing defect constants",ylim=c(0,50),legend.text=T,args.legend=list(x=10,y=25))#,xpd=FALSE)
# error bars on the last column
yvals = as.numeric(packdef_constants["PackDef_cst_all_blocks",])
perc_errors = as.numeric(packdef_constants["error_all_blocks",]) / abs(yvals)
errbar = perc_errors*yvals
segments(bp[5,],yvals+errbar,bp[5,],yvals-errbar,lwd=2)
segments(bp[5,] - 0.2, yvals - errbar, bp[5,] + 0.2, yvals - errbar, lwd=2)
segments(bp[5,] - 0.2, yvals + errbar, bp[5,] + 0.2, yvals + errbar, lwd=2)

####
# Now we plot the final decays + errors computed with block averaging 
# (for each packdef) on a barplot
####
# simplify constant name
yvals = as.numeric(packdef_constants["PackDef_cst_all_blocks",])
error_vals = as.numeric(packdef_constants["error_all_blocks",])

# Handle zero or non-finite values
yvals[!is.finite(yvals)] <- 0
error_vals[!is.finite(error_vals)] <- 0

# Calculate percentage errors safely
perc_errors = ifelse(yvals != 0, error_vals / abs(yvals), 0)

# Set ylim safely
ymin = 0
ymax = max(yvals*1.2, 1)  # 20% padding, min 1 if all zeros
par(mar=c(7, 4.6, 4, 2) + 0.1)

# Create barplot
bp = barplot(
  yvals,
  names.arg="",
  ylab=expression(paste("defect size constant (",A^2,")")),
  ylim=c(ymin, ymax),
  cex.axis=1.4,
  cex.lab=1.4,
  xaxt="n",
  col=palette(),
  main="Packing defect constants computed by block averaging",
  xpd=FALSE
)

# X-axis labels
LABELS = c("deep","shallow","all")
axis(1, at=bp, labels=FALSE, tick=FALSE)
text(bp, par("usr")[3] - 0.6, labels=LABELS, pos=1, xpd=TRUE, cex=1.3, adj=1, col=palette())

# Box around plot
box(lwd=2)

# Add text labels on bars
text_yvals = sprintf("%4.1f", round(yvals,1))
text_perc = ifelse(max(yvals) > 0,
                   sprintf("%4.1f%%", round((yvals/max(yvals))*100,1)),
                   "0.0%")

for (i in 1:length(bp)) {
  text(bp[i], yvals[i]*0.9, label=bquote(paste(.(text_yvals[i]), ~ A^2)), cex=1)
  text(bp[i], yvals[i]*0.9 - 0.8, label=text_perc[i], cex=1)
}

# Error bars (skip if error = 0)
errbar = perc_errors * yvals
segments(bp, yvals + errbar, bp, yvals - errbar, lwd=2)
segments(bp - 0.1, yvals - errbar, bp + 0.1, yvals - errbar, lwd=2)
segments(bp - 0.1, yvals + errbar, bp + 0.1, yvals + errbar, lwd=2)


# Plot average defect size
filename <- paste(LONG_PREFIX, "_", "shallow", ".txt", sep="")

# Read all lines
lines <- readLines(filename)

areas_per_frame <- list()   # store defect areas per frame
frame_areas <- c()
inside_frame <- FALSE

for (line in lines) {
  line <- trimws(line)
  
  if (startsWith(line, "## MatrixSize")) {
    # Start of new frame
    if (length(frame_areas) > 0) {
      areas_per_frame[[length(areas_per_frame) + 1]] <- frame_areas
    }
    frame_areas <- c()
    inside_frame <- TRUE
  } else if (startsWith(line, "## Total")) {
    # skip total line
    next
  } else if (startsWith(line, "##")) {
    inside_frame <- FALSE
  } else if (inside_frame && nchar(line) > 0) {
    cols <- unlist(strsplit(line, "\\s+"))
    if (length(cols) >= 2) {
      frame_areas <- c(frame_areas, as.integer(cols[2]))
    }
  }
}

# Save last frame
if (length(frame_areas) > 0) {
  areas_per_frame[[length(areas_per_frame) + 1]] <- frame_areas
}

# Compute average defect size per frame (only area >= 15)
avg_defect_size <- sapply(areas_per_frame, function(x) {
  filtered <- x[x >= 15]
  if (length(filtered) == 0) return(NA)   # no defects >= 15
  mean(filtered)
})

#compute over all frames 
# Remove NA values (frames without defects >= 15)
avg_defect_size_clean <- avg_defect_size[!is.na(avg_defect_size)]

# Overall average defect size across frames
overall_avg <- mean(avg_defect_size_clean)
message(overall_avg)

# Overall standard deviation across frames
overall_sd <- sd(avg_defect_size_clean)



# put in a data frame


# # Plot histogram of defect sizes
# hist(packdef_data, breaks=30, col="skyblue", main="Histogram of Defect Sizes",
#      xlab=expression(paste("Defect area (", A^2, ")")), ylab="Count")

# # Add a vertical line for the average
# abline(v = avg_defect_size, col="red", lwd=3)
# text(avg_defect_size, max(hist(packdef_data, plot=F)$counts) * 0.9,
#      labels=paste("Mean =", round(avg_defect_size,2)), col="red", pos=4)

# Optional: barplot for average + SD
bp <- barplot(overall_avg, beside=TRUE, ylim=c(0, max(overall_avg + overall_sd , 10)),
        names.arg="All defects", col="orange", ylab=expression(paste("Average size defect area (",A^2,")")))

# Add error bar for SD
arrows(x0 = bp, y0 = overall_avg - overall_sd, 
       x1 = bp, y1 = overall_avg + overall_sd, 
       angle = 90, code = 3, length = 0.1, col = "black")

# Plot average defect coverage 



# close the plotting device
dev.off()
# stop sinking output to the txt file
sink()