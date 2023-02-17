##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")}

inputDataLoc = args[1]
outputLoc = args[2]
rho = as.numeric(args[3])
shaPRSscriptLoc = args[4]
discardAmbiguousSNPs=T
if(length(args) > 4) {discardAmbiguousSNPs= args[5] =="0"}


# load exteral functions
source(shaPRSscriptLoc) # '/nfs/users/nfs_m/mk23/scripts/shaPRS.R'
library(qvalue)
# install.packages("devtools")
# library("devtools")
# install_github("jdstorey/qvalue")

# Debug

# inputDataLoc="C:/0Datasets/shaPRS/crossAncestry/EUR_JAP_asthma_SE_meta_short"  # "C:/0Datasets/_badData"
# outputLoc ="C:/0Datasets/shaPRS/crossAncestry/asthmashort"
##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load data
inputData= read.table(inputDataLoc, header = T)

# 2. lFDR estimation
results = shaPRS_adjust(inputData, rho = rho, discardAmbiguousSNPs =discardAmbiguousSNPs)
#lfdr_qvals <- results$lFDRTable$lfdr_qvals

# 3. write out a table of lFDR values for each SNP
lFDRTable <- results$lFDRTable # data.frame(inputData$SNP, results$lFDRTable$lfdr_qvals, results$lFDRTable$Q_vals)
colnames(lFDRTable) = c("SNP", "lFDR", "Qval")
filename = paste(outputLoc, "_SNP_lFDR" , sep="")
write.table(lFDRTable, filename, row.names = F, col.names = T, quote = FALSE)
print(paste("written lFDRs and Qvals for SNPs to",filename))

