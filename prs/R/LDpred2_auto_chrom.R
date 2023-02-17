#############################
# LDpred2-auto script: takse a PLINK GWAS result and produce an LDpred2 PRS file
# make sure to run this script with ~50GB of RAM and ~30 cores for optimum speed
# according to Florian Prive, use 15 cors @ 125GBs of ram: https://github.com/privefl/bigsnpr/issues/172
#############################

# imports
library(tidyr)
library(bigsnpr)
library(dplyr)
library(R.utils)
library(bigreadr)

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 9) {stop("not enough Arguments received")}


# arguments
baseLoc=args[1] # the location of the per-chrom LDref panels and the HapMap3 map file
sumstatsLoc=args[2] # location of my 10 col formatted sumtats
chr=as.numeric(args[3])  # the chrom number
# n_eff =  as.numeric(args[4]) # number of individuals in study, UNUSED, only kept for backwards compatibility
outLoc=args[5] # output location
outName=args[6] # outputfile name
testSetLoc=args[7] # test set PLINK .bed file
shaPRSscriptLoc= args[8]
binary_outcome= args[9] =="1"
SNPQC= F

if ( is.na(args[10]) == F && args[10] == "1") { SNPQC = T}

print(paste0("SNP QC enabled: ", SNPQC))
# load exteral functions
source(shaPRSscriptLoc) # '/nfs/users/nfs_m/mk23/scripts/shaPRS.R'

#############################
# LOCAL DEBUG VARS
#baseLoc="/lustre/scratch114/projects/crohns/mk23/0Thesis/shaprs/raw/ldpred2/"
##sumstatsLoc="/lustre/scratch114/projects/crohns/mk23/missense_eqtl_epistasis/data_analysis/eqtl_missense_merged_all_sanity.assoc"
#sumstatsLoc="/lustre/scratch114/projects/crohns/mk23/0Thesis/shaprs/simPhe2/5000/0.5_1.0/A50_B50/size/4_pheA.qassoc"

#NCORES=nb_cores() # should limit this to what requested from cluster
#n_eff=17554
#outLoc= "/lustre/scratch114/projects/crohns/mk23/0Thesis/shaprs/scratch/ldpred2/"
#outName= "ldpred2test"
#testSetLoc="/lustre/scratch114/projects/crohns/mk23/0Thesis/shaprs/raw/pheA_testSet.bed"

#########################################################

minh2 = 0.001
set.seed(42)
# 1. load map of SNPs that we have available LD data for
mapLoc= paste0(baseLoc, "map.rds") # the HapMap3 map file
ldRefpanelLoc=paste0(baseLoc,"ldrefpanel") # LD-ref files for each chrom
tempDirLoc=paste0(outLoc,outName,"temp", chr,"/") # temporary dir where all dump files will go
map_ldref <- readRDS(mapLoc) # load the HapMap3


# 2. load summary stats data
sumstats <- fread2(sumstatsLoc)

# subset summary data to 
sumstats = sumstats[sumstats$chr == chr,]



# 2b. load PLINK file (this is used to obtain the a0 and a1)
bimData <- fread2(paste0(testSetLoc, ".bim")) # before loading the full file, load just the bim, to know what SNPs we will have
bimData$id = 1:nrow(bimData) # note the indices of the SNPs according to the original bimdata, that included all 22 chroms, once merged, this will be used to only load the correct chrom's variants

famData <- fread2(paste0(testSetLoc, ".fam")) # before loading the full file, find out how many people, and only load max 10K


# to make this faster, and lower memory consumption subset the individuals to a subset of a maximum of 10K
numTestIndis = min (nrow(famData), 10000 )
testIndiIndices = sample(1:nrow(famData), numTestIndis)

# find the intersection of the sumstats SNPs that we have access to on the PLINK binary, and subset the PLINK bin
bimData$SNP_chrPos = paste(bimData$V1,bimData$V4, sep=":")
sumstats$chrom_SNPs = paste(sumstats$chr,sumstats$pos, sep=":")
sumstats_bim = merge(bimData, sumstats, by.x = "SNP_chrPos", by.y ="chrom_SNPs")

sumstats_bim <- sumstats_bim[!is.na(as.numeric(as.character(sumstats_bim$b))),]
sumstats_bim$b = as.numeric(as.character(sumstats_bim$b ))

sumstats_bim <- sumstats_bim[!is.na(as.numeric(as.character(sumstats_bim$N))),]
sumstats_bim$N = as.numeric(as.character(sumstats_bim$N ))

# need to align the sumstats against the BIM, as the effect alleles may not match the A1 of the PLINK file
sumstats_bim = alignStrands(sumstats_bim, A1.x ="V5", A2.x ="V6", A1.y ="A1", A2.y ="A2")

# 2. Align PheB/B alleles
misalignedAlleleIndices = which( as.character(sumstats_bim$V5) != as.character(sumstats_bim$A1) ) # compare as character, as if we have non-SNPs with different alleles factors will break
sumstats_bim$b[misalignedAlleleIndices] = -sumstats_bim$b[misalignedAlleleIndices] # flip effects for phe B
if(length(misalignedAlleleIndices) > 0) message(paste0(length(misalignedAlleleIndices)), " misaligned allele(s) effects were reversed" )



tmp <- tempfile(tmpdir = tempDirLoc)
testSet= snp_readBed2(paste0(testSetLoc,".bed"), ind.row = testIndiIndices, ind.col = sumstats_bim$id , backingfile = tmp) # put the backing file in the temp dir so they could be deleted together with the rest
testData <- snp_attach(testSet)
testData$map$id = 1:nrow(testData$map) # save the indices that refer to the originally loaded PLINK file

# combine the genotype data with the summary statistics into the full LDpred2 required format
sumstats_processed = cbind(testData$map[-3], sumstats_bim$b,sumstats_bim$se, sumstats_bim$N ) # SE is the standard error of either the Beta (for quantitative), or the log odds (IE plink already outputs it on log scale so don't need to take logs again)
#sumstats_processed = cbind(testData$map[-3], sumstats_bim$b,sumstats_bim$se ) # SE is the standard error of either the Beta (for quantitative), or the log odds (IE plink already outputs it on log scale so don't need to take logs again)
#sumstats_processed = sumstats_processed[,c(1,3,4,5,6,7,2)]
names(sumstats_processed) <- c("chr", "rsid", "pos", "a1", "a0", "id", "beta", "beta_se", "N")
#names(sumstats_processed) <- c("chr", "rsid", "pos", "a1", "a0", "id", "beta", "beta_se")
sumstats_processed$n_eff <- sumstats_processed$N # we have per-SNP effective sample size
head(sumstats_processed)


# 3. match summary stats to the HapMap3 (will result in keeping only a subset of original variants)
info_snp <- snp_match(sumstats_processed, map_ldref)

#info_snp <- snp_match(sumstats_processed, map_ldref, match.min.prop = 0.01)

# need to keep track of the SNPs that were successfully matched to the HapMap3
SNPsKept_indices = info_snp$id  #  unlist(info_snp['_NUM_ID_']) # this doesnt work, for some reason the _NUM_ID_ refers to the wrong indices

# 4. Pre-fit QC for adjusting the betas
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB))) # standard deviation of genotypes, same as sd(Gj) in the LDpre2 paper

# apply different formula to get sd of summary statistics  for binary/continuous traits
if(binary_outcome) {
  sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2)) # standard deviation of summary stats
} else {
  sd_y_est = median( sd_ldref * info_snp$beta_se * sqrt(info_snp$n_eff))
  sd_ss = with(info_snp, sd_y_est / sqrt(n_eff * beta_se^2))
}



if(SNPQC) {
  # find SNPs that fail any of the below QC criteria
  is_bad <-  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
  print("performing SNP QC")
  df_beta <- info_snp[!is_bad, ] # exclude the bad SNPs
} else {
  print("SNP QC disabled")
  is_bad = rep(FALSE, length(sd_ss) ) # make all SNPs not bad, so that this mask can be later used
  
  df_beta <- info_snp
}
  

# 5 Create LD ref data by concatenating each chrom 
# (this step cannot be cached and has to happen after Pre-fit QC, 
# as the subset of SNPs that pass both matching against HapMap3 and Pre-fit QC will differ for each run, 
# and LDref panel has to match the surviving SNPs exactly)
tmp <- tempfile(tmpdir = tempDirLoc)


  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  
  chrLDLoc = paste0(baseLoc,"LD_chr", chr, ".rds")
  
  #if( file.exists(chrLDLoc) ) {
  
  # load LD from disk for the given subset
  corr_chr <- readRDS(chrLDLoc)[ind.chr3, ind.chr3]
  
  if( typeof(corr_chr) == "double" ) {
    print(paste0("convering dense to sparse matrix"))
    corr_chr = as(corr_chr ,"dgCMatrix")

  }
  

    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)

 # }


 

# 6. Heritability estimation via  LD score regression to be used as a starting value in LDpred2-auto
# replaced nrow(map_ldref) with  nrow(df_beta), as if there were a lot of variants excluded in the QC this resulted in abnormal h2 estimates
(ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(df_beta),chi2 = (beta / beta_se)^2,sample_size = n_eff)))
h2_est <- ldsc[["h2"]]
write.table(h2_est, paste0(outLoc,outName, "_h2"), sep = "\t", row.names = F, col.names = F, quote = FALSE)

h2_p = 2*pnorm(-abs(h2_est/ ldsc[4])) # this is the estimate / its standard error
if(h2_p > 0.05) {
  print(paste0("h2 is not significant, setting it to o"))
  h2_est = 0
}
if(h2_est < minh2) { # sometimes we get negative h2 estimates
  print(paste0("h2 estimate seems very low at ", h2_est, ", setting it to: ", minh2))
  h2_est = minh2
}
if (h2_est > 1) {
  print(paste0("h2 estimate seems too high  at ", h2_est, ", setting it to: 1"))
  h2_est = 1
}

# 7. Main LDpred2 auto fit
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,vec_p_init = seq_log(1e-4, 0.9, 30), allow_jump_sign=F, shrink_corr= 0.95, verbose = TRUE)  


beta_auto <- sapply(multi_auto, function(auto) auto$beta_est) # get the adjusted SNP estimates for all trials, this is p x NCORES 952083 x 30



# Write results for all LDpred Betas to disk
snp_and_betas = cbind(df_beta$rsid.ss,beta_auto)
gz1 <- gzfile(paste0(outLoc,outName,"all_betas.gz"), "w") # compress them to save space, as these are unlikely to be used directly
write.csv(snp_and_betas, gz1)
close(gz1)



# 8. Post-fit QC: determine best Auto beta as the average of the Betas excluding ones that produced outlier PRS
# source: https://github.com/privefl/paper-ldpred2/blob/master/code/run-ldpred2.R#L101-L108
# big_prodMat crashes with NAs, so need to fill any missing (shouldn't be many) https://github.com/privefl/bigsnpr/issues/124
imputed = snp_fastImputeSimple(testData$genotypes, method = "mean2") 

TestGenotype_hapmap3_matched = imputed[,SNPsKept_indices] # subset the genotype file to the HapMap3 matched SNPs

# subset the genotype to the QC passed SNPs
tmp <- tempfile(tmpdir = tempDirLoc)
TestGenotype_hapmap3_matched_QC = as_FBM(TestGenotype_hapmap3_matched[,!is_bad], type = "double",tmp)
pred_auto <- big_prodMat(TestGenotype_hapmap3_matched_QC, beta_auto)


sc <- apply(pred_auto, 2, sd)
final_beta_auto <- rowMeans(beta_auto[, abs(sc - median(sc)) < 3 * mad(sc)])


# 9. produce the final output in PLINK PRS format that would work for PLINK --score sum  (columns: rsid, effect_allele, beta)
PRSFILE = cbind(df_beta$rsid.ss, df_beta$a1, final_beta_auto)
write.table(PRSFILE, paste0(outLoc,outName), sep = "\t", row.names = F, col.names = F, quote = FALSE)
print(paste("written ",nrow(PRSFILE)," LDpred adjusted betas to", paste0(outLoc,outName)))


# 10. delete temp folder
unlink(tempDirLoc, recursive = TRUE)
