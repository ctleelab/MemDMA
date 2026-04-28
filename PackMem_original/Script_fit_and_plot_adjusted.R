################################
# Packdef constant calculation #
################################

# Open a pdf device
pdf(paste("Res_",".pdf",sep=""))

# Initialize a data frame to store packdef constants + errors
packdef_constants = data.frame(
  matrix(0, nrow=6, ncol=3),
  row.names=c("PackDef_cst_global","PackDef_cst_block1","PackDef_cst_block2","PackDef_cst_block3","PackDef_cst_all_blocks","error_all_blocks")
)
names(packdef_constants) = c("deep","shallow","all")

# print the results also to a text file
sink(paste("Res_",".txt",sep=""))

######
# Loop over defect types
######
for (DEFECT in c("deep","shallow","all")) {

  # Load PackMem data
  filename = paste(LONG_PREFIX,"_",DEFECT,".txt",sep="")
  packdef_data = read.table(filename)[,2] 
  # Remove NA/NaN and non-finite values
  packdef_data = packdef_data[is.finite(packdef_data)]

  # If no defects, create a single zero value
  if(length(packdef_data) == 0){
    packdef_data = 0
  }

  # Compute histogram
  H = hist(packdef_data, breaks=seq(0, max(packdef_data,1)) + 0.5, plot=F) 
  x = H$mids
  y = H$counts / sum(H$counts)

  # Avoid log(0)
  y[y == 0] = 1e-10

  # Compute FITs safely
  if(all(packdef_data == 0)){
    global_inv_decay = 0
    FITS_3blocks = c(0,0,0)
  } else {
    FIT = fit_decay(x, y, LIMX, LIMY)
    global_inv_decay = abs(1 / coef(FIT)[2])
    FITS_3blocks = block_averaging(packdef_data)
  }

  # Print results
  print(paste("Results for", DEFECT, "defects"))
  print(paste("Total number of defects =", length(packdef_data)))
  print(paste("Global PackDef constant:", round(global_inv_decay, PRECISION), "A^2"))
  print(paste("Block 1:", round(FITS_3blocks[1], PRECISION), "A^2"))
  print(paste("Block 2:", round(FITS_3blocks[2], PRECISION), "A^2"))
  print(paste("Block 3:", round(FITS_3blocks[3], PRECISION), "A^2"))
  inv_decay_block = mean(FITS_3blocks)
  error_inv_decay_block = sd(FITS_3blocks)
  print(paste("Mean +/- sd on 3 blocks:", round(inv_decay_block, PRECISION), "+/-", round(error_inv_decay_block, PRECISION), "A^2"))
  print("")

  # Store results
  packdef_constants["PackDef_cst_global", DEFECT] = global_inv_decay
  packdef_constants["PackDef_cst_block1", DEFECT] = FITS_3blocks[1]
  packdef_constants["PackDef_cst_block2", DEFECT] = FITS_3blocks[2]
  packdef_constants["PackDef_cst_block3", DEFECT] = FITS_3blocks[3]
  packdef_constants["PackDef_cst_all_blocks", DEFECT] = inv_decay_block
  packdef_constants["error_all_blocks", DEFECT] = error_inv_decay_block

  # -----------------------------
  # Plot histogram safely
  # -----------------------------
  # Only draw FIT if not all zeros
  par(mar=c(5,4.5,4,2)+0.1)
  plot(1,1,type="n", xlim=c(0,100), ylim=c(-10,-3.5),
       ylab="Probability", xlab=expression(Defect~area~(A^2)),
       axes=F, main=DEFECT, cex.lab=1.6, new=T)
  box(lwd=2)
  box()
  axis(1, at=seq(0,200,10), labels=F, lwd.ticks=2)
  axis(1, at=seq(0,200,20), labels=seq(0,200,by=20), lwd.ticks=3)
  LABEL = c(expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), expression(10^0))
  axis(2, at=log(10^(seq(-5,0,1))), labels=LABEL, lwd.ticks=3)
  LOGAXIS = c()
  for (i in 5:1) LOGAXIS = c(LOGAXIS, log(10^-i * 1:9))
  LOGAXIS = c(LOGAXIS, 1)
  axis(2, at=LOGAXIS, labels=FALSE, lwd.ticks=1.5)
  legend(45, -3.4, c(PREFIX), col=c(1:length(palette())), lty=1, lwd=3, cex=1.3)
  abline(h=log(LIMY), lwd=2)
  abline(v=LIMX, lty=2, lwd=2)

  # Add histogram points
  points(x, log(y), col=1, lwd=2)
  if(!all(packdef_data == 0)){
    FIT = fit_decay(x, y, LIMX, LIMY)
    abline(FIT, col=1, lwd=3)
  }

} # END DEFECT loop

# -----------------------------
# Plot barplots safely
# -----------------------------
par(mar=c(7,4.6,4,2)+0.1)
M = as.matrix(packdef_constants[1:5,])
bp = barplot(M, beside=T,
             ylab=expression(paste("defect size constant (",A^2,")")),
             cex.axis=1.4, cex.lab=1.4,
             col=c("black", gray(1/3), gray(2/3), gray(3/3), "red"),
             main="Packing defect constants",
             ylim=c(0, max(25, max(M, na.rm=TRUE))), # ensures finite ylim
             legend.text=T, args.legend=list(x=10,y=25))
# Error bars
yvals = as.numeric(packdef_constants["PackDef_cst_all_blocks",])
perc_errors = ifelse(yvals==0, 0, as.numeric(packdef_constants["error_all_blocks",]) / abs(yvals))
errbar = perc_errors * yvals
segments(bp[5,], yvals+errbar, bp[5,], yvals-errbar, lwd=2)
segments(bp[5,]-0.2, yvals-errbar, bp[5,]+0.2, yvals-errbar, lwd=2)
segments(bp[5,]-0.2, yvals+errbar, bp[5,]+0.2, yvals+errbar, lwd=2)

# Close devices
dev.off()
sink()