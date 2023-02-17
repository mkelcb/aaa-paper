#############################
# Imputes N_eff from summary statistics
# using formula from "Identifying and correcting for misspecifications in GWAS summary statistics and polygenic scores" (https://www.biorxiv.org/content/10.1101/2021.03.29.437510v3)
# var_G formula is from the main LDpred2paper
#############################


# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")}


# arguments
baseLoc=args[1] # location of file with header signature:  rsID, Beta, Beta_SE, AF
outLoc=args[2] # location of my 10 col formatted sumtats

info_snp=read.table(baseLoc  ,header=T)

# impute effective sample size
var_G <- with(info_snp, 2 * AF * (1 - AF)) # varance of genotypes

n_eff = with(info_snp, ( 4/ var_G - Beta^2) /  (Beta_SE^2) )

info_snp = cbind(info_snp, n_eff) 


write.table(info_snp,file=outLoc, quote = F, row.names = F, sep = "\t",append=F)

print(paste("written num effective sample sizes", nrow(info_snp), "to", outLoc ))

