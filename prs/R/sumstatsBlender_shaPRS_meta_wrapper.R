##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 6) {stop("not enough Arguments received")}

subphenoLoc = args[1]
subpheno_otherLoc = args[2]
blendFactorLoc = args[3]
rho = as.numeric(args[4]) 
outputLoc = args[5]
shaPRSscriptLoc= args[6]
discardAmbiguousSNPs=T
if(length(args) > 6) {discardAmbiguousSNPs= args[7] =="0"}

# load exteral functions
source(shaPRSscriptLoc) # '/nfs/users/nfs_m/mk23/scripts/shaPRS.R'
library(qvalue)


# Debug Vars
# subphenoLoc="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheA_sumstats"
# subpheno_otherLoc="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheB_sumstats"
# blendFactorLoc ="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheBlend_lFDR_meta_SNP_lFDR"
# outputLoc ="C:/0Datasets/shaPRS/0shaPRS/0ukbb/debug/1_pheBlend_sumstats_meta"
# rho = 0

# subphenoLoc="C:/0Datasets/shaPRS/crossAncestry/comparison/_asthmaOldPheA_5"
# subpheno_otherLoc="C:/0Datasets/shaPRS/crossAncestry/comparison/_asthmaOldPheB_5"
# blendFactorLoc ="C:/0Datasets/shaPRS/crossAncestry/comparison/_asthmaOld_lFDR_meta_SNP_lFDR_5"
# outputLoc ="C:/0Datasets/shaPRS/crossAncestry/comparison/out"
# rho = 0

##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load phenos
subpheno= read.table(subphenoLoc, header = T)

# load other subpheno to perform meta analysis
subpheno_other= read.table(subpheno_otherLoc, header = T)


# INVERSE VARIANCE FIXED EFFECT META ANALYSIS:
CombinedPheno = inverse_metaAnalaysis(subpheno, subpheno_other, rho = rho, discardAmbiguousSNPs =discardAmbiguousSNPs)

blendingFactors= read.table(blendFactorLoc, header = T)

# 4. create new data frame to store the new summary stats
#if( rho == 0) {
#  print("NO overlap shaPRS")
#  blendedSumstats = shaPRS_blend(subpheno, CombinedPheno, blendingFactors)
#} else  {
#  print("YES overlap shaPRS")
blendedSumstats = shaPRS_blend_overlap(subpheno, subpheno_other, blendingFactors,rho, discardAmbiguousSNPs =discardAmbiguousSNPs)
#}
# 5. write blended stats to disk
write.table(blendedSumstats, outputLoc, sep = "\t", row.names = F, col.names = T, quote = FALSE)
print(paste("written blended sumstats to",outputLoc))


# 6. write Combined stats to disk too
filen=paste0(outputLoc,"_combined")
write.table(CombinedPheno, filen, sep = "\t", row.names = F, col.names = T, quote = FALSE)
print(paste("written Combined sumstats to",filen))




#(SOFT)  / 

#1000G sumstats Meta-analysis:
#  / 123504

#Exome Meta-analysis
# / 78240


