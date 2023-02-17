# performs SNP QC on a sumtats, given the sumstats file and an LD map file (only used to get the MAF)
# imports
library(tidyr)
library(bigsnpr)
library(ggplot2)

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 4) {stop("not enough Arguments received")}


# arguments
baseLoc=args[1] # the location of the per-chrom LDref panels and the HapMap3 map file, used to 
sumstatsLoc=args[2] # sumstats file
outLoc=args[3] # output location
binary_outcome= args[4] =="1"



#######################



# sd_ldref =1
# n_eff=100
# beta_se=2
# sd_y_est=1
# sd_ss = sd_y_est / sqrt(n_eff * beta_se^2) # standard deviation of the summary stats
# sd_ss
# bad_SNPs = sd_ss < (0.5 * sd_ldref) # this throws away ~400K SNPs
# bad_SNPs
# 

###########################

# cases= 10801
# controls=137914
# N_eff = 4 / (1 / cases + 1 / controls) 
# round(N_eff) # 40066
# 
# cases= 10801
# controls=137914
# N_eff =  cases +  controls
# round(N_eff) # 40066


mapLoc= paste0(baseLoc, "map.rds") # the HapMap3 map file
map_ldref <- readRDS(mapLoc) # load map file that has the MAF info # 1054330

sumstats <- read.table(sumstatsLoc, header = T) # 995412


# remove non numeric data
sumstats <- sumstats[!is.na(as.numeric(as.character(sumstats$se))),]
sumstats <- sumstats[!is.na(as.numeric(as.character(sumstats$N))),]
# now actually cast them to numeric
sumstats$se = as.numeric(as.character(sumstats$se ))
sumstats$N = as.numeric(as.character(sumstats$N ))


info_snp = merge(sumstats, map_ldref, by.x="SNP", by.y="rsid") # 933066

# 4. Pre-fit QC for adjusting the betas
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB))) # standard deviation of genotypes, same as sd(Gj) in the LDpre2 paper

# apply different formula to get sd of summary statistics  for binary/continuous traits
if(binary_outcome) {
  sd_ss <- with(info_snp, 2 / sqrt(N * se^2)) # standard deviation of summary stats
} else {
  sd_y_est = median( sd_ldref * info_snp$se * sqrt(info_snp$N))
  sd_ss = with(info_snp, sd_y_est / sqrt(N * se^2))
}

# find SNPs that fail any of the below QC criteria
is_bad <-  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
#is_bad <- sd_ss > (sd_ldref + 0.1)


# write out excluded SNPs
print(paste0("out of ", nrow(info_snp), " QC excluded ", sum(is_bad), ", and kept ",length(info_snp$SNP[!is_bad])," SNPs, written to: ",outLoc ))


g = qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations in the validation set",
       y = "Standard deviations derived from the summary statistics",
       color = "Removed?")

ggsave(paste0(outLoc,".png"), width = 10, height = 7)





write.table(info_snp$SNP[!is_bad], outLoc, sep = "\t", row.names = F, col.names = F, quote = FALSE)


