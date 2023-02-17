##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")}

subphenoLoc = args[1]
subpheno_otherLoc = args[2]
outputLoc = args[3]
shaPRSscriptLoc= args[4]

# load exteral functions
source(shaPRSscriptLoc) # '/nfs/users/nfs_m/mk23/scripts/shaPRS.R'


##############################################################
#                         MAIN SCRIPT                        #
##############################################################

# 1. load phenos
subpheno= read.table(subphenoLoc, header = T)

# load other subpheno to perform meta analysis
subpheno_other= read.table(subpheno_otherLoc, header = T)


# INVERSE VARIANCE FIXED EFFECT META ANALYSIS:
CombinedPheno = inverse_metaAnalaysis(subpheno, subpheno_other)

# 6. write Combined stats to disk too
write.table(CombinedPheno, outputLoc, sep = "\t", row.names = F, col.names = T, quote = FALSE)
print(paste("written Combined sumstats to",outputLoc))





