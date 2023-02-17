
###############################################

# AAA adjunct data and final PRS calculation
# assumes that AAA_pub.sh has been completed
# depends on variables and reusable functions loaded into the session from AAA_functions.sh

# I process adjunct datasets into the common summary data format
# I then combine these summary data
# then fit PRS and generate final model
###############################################

# find out liability scale h2 for AAA
# using 1% prevalance for AAA: https://bmccardiovascdisord.biomedcentral.com/articles/10.1186/s12872-019-1265-2


##################################
# (I) MEGASTROKE study of EUR indis
/rds/project/jmmh2/rds-jmmh2-results/public/gwas/gwas_catalog/MalikR_MEGASTROKE_NatGen_2018/MEGASTROKE.1.AS.EUR.out
# 8255861

#MarkerName Allele1 Allele2 Freq1 Effect StdErr P-value
#rs2326918 a g 0.8455 0.0002 0.0141 0.9915

#MarkerName: SNP rsID or chromosome:position if rsID not available
#Allele1: Effect allele (coded allele)
#Allele2: Non effect allele (non coded allele)
#Freq1: Effect allele frequency
#Effect: Overall estimated effect size for the effect allele
#StdErr: Overall standard error for effect size estimate
#P-value: Meta-analysis P-value using regression coefficients (beta and standard error)

# MegaStroke_AS: 40,585 cases; 406,111 controls from  MALIKR_29531354_READMEv2.docx

# cases = 40585
# controls = 406111
#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 147590

$munge_sumstats --sumstats /rds/project/jmmh2/rds-jmmh2-results/public/gwas/gwas_catalog/MalikR_MEGASTROKE_NatGen_2018/MEGASTROKE.1.AS.EUR.out --a2 Allele2 --a1 Allele1  --N 147590 --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $MegaStroke_AS_sumstatsStem

# AAA v MegaStroke_AS: 
$lsdc --rg $rawLocV4$"AAD_abdominal_hm3_sums".sumstats.gz,$MegaStroke_AS_sumstatsStem.sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $resultsLoc$'AAD_abdominal/rG_megaStroke_AS'
#Genetic Correlation: 0.3926 (0.1903)
#P: 0.0391


# MEGASTROKE AS: this is missing the chr/pos cols, so we have to match these from our hapmap3 .bim
awk 'FNR == NR { file1[ $2 ] = $1"\t"$4; next; } 
{ if (FNR == 1) {print "chr\tpos\t"$0}
else{
if ( $1 in file1) {print file1[$1]"\t"$0 }
}}   ' OFS='\t' $hapmap3_b37bim /rds/project/jmmh2/rds-jmmh2-results/public/gwas/gwas_catalog/MalikR_MEGASTROKE_NatGen_2018/MEGASTROKE.1.AS.EUR.out > $scratchLoc$'MEGASTROKE.1.AS.EUR_chr_pos'
head $scratchLoc$'MEGASTROKE.1.AS.EUR_chr_pos'
wc -l $scratchLoc$'MEGASTROKE.1.AS.EUR_chr_pos' # 1210482

#  1   2      3         4       5      6     7       8      9
#Chr Pos  MarkerName Allele1 Allele2 Freq1 Effect StdErr P-value
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
awk  -v numIndis=147590 ' FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { print $1"\t"$2"\t"$3"\t"toupper($4)"\t"toupper($5)"\t"$6"\t"$7"\t"$8"\t"$9"\t"numIndis}
} ' OFS='\t' $scratchLoc$'MEGASTROKE.1.AS.EUR_chr_pos' > $scratchLoc$'MEGASTROKE_AS_raw'
head $scratchLoc$'MEGASTROKE_AS_raw'


# perform QC
remove_ambiguous_alleles $scratchLoc$'MEGASTROKE_AS_raw'
wc -l $scratchLoc$'MEGASTROKE_AS_raw_noAmbiguousAlleles' # 1116941

# QC
arguments='/home/mk907/scripts/R/LDpred2_QC.R '$baseLDpredLoc$'/ '$scratchLoc$'MEGASTROKE_AS_raw_noAmbiguousAlleles '$rawLocV4$'MEGASTROKE_HM3_keep 1'
Rscript $arguments
wc -l $rawLocV4$'MEGASTROKE_HM3' # 1038450


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $rawLocV4$'MEGASTROKE_HM3_keep' $scratchLoc$'MEGASTROKE_AS_raw_noAmbiguousAlleles'  > $rawLocV4$'MEGASTROKE_HM3'
head $rawLocV4$'MEGASTROKE_HM3'



###############################
# (II) Nelson 2017 CAD meta analysis: $rawLocV4$'NelsonCAD_HM3'
# http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004787
# This useed the interim UKBB release

mkdir -p /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/Nelson_CAD/
cd /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/Nelson_CAD/
wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004787/harmonised/28714975-GCST004787-EFO_0001645-Build37.f.tsv.gz

wget http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST004001-GCST005000/GCST004787/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz


gunzip 28714975-GCST004787-EFO_0001645-Build37.f.tsv.gz


gunzip UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz


head 28714975-GCST004787-EFO_0001645-Build37.f.tsv

head UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt



# Nelson CAD
#    1                2                  3                    4               5                  6                    7            8               9              10           11        12  
#variant_id      chromosome      base_pair_location      effect_allele   other_allele    effect_allele_frequency    logor    standard_error     p_value       n_samples       exome   info_ukbb       ci_upper        odds_ratio      beta    ci_lower
#rs561255355          1                569406                  G               A                 0.99858            0.05191      0.27358        0.849496        154959        yes     0.41    NA      NA      NA      NA
#rs12562034           1                768448                  A               G                 0.12029            0.00096      0.0144         0.946778        269334         no      1       NA      NA      NA      NA

# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N


# Nelson CAD sometimes had NA in their CHR/BP columns so we have to replace those by RSid from hapmap3
awk 'FNR == NR { file1[ $2 ] = $1"\t"$4; next; } 
{ if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else{
if ( $1 in file1) {print file1[$1]"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }
}}   ' OFS='\t' FS="\t" $hapmap3_b37bim FS="\t" /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/Nelson_CAD/28714975-GCST004787-EFO_0001645-Build37.f.tsv > $scratchLoc$'Nelson_CAD_raw'
head $scratchLoc$'Nelson_CAD_raw'
wc -l $scratchLoc$'Nelson_CAD_raw' #1221228


# impute N_eff, as the per SNP sample sizes seem to differ from the official account 

# extract data with signature: rsID, Beta, Beta_SE, AF
awk '{if(FNR ==1) {print "SNP\tBeta\tBeta_SE\tAF\tN"}  else {print $3"\t"$7"\t"$8"\t"$6"\t"$10} }' $scratchLoc$'Nelson_CAD_raw' > $scratchLoc$'Nelson_CAD_raw_in'
head $scratchLoc$"Nelson_CAD_raw_in"

# run R script 
Rscript /home/mk907/scripts/R/ImputeNeff.R $scratchLoc$'Nelson_CAD_raw_in' $scratchLoc$'Nelson_CAD_raw_in_N_eff'
head $scratchLoc$'Nelson_CAD_raw_in_N_eff'


# produce the final summary stats file
awk 'FNR == NR { 
file1[$1] = $6; next; } 
FNR <= NR { 
gsub("\r","",$3); 
if(FNR == 1) {print $0}
else if ( $3 in file1 ) { 
$10=  sprintf("%0.f", file1[$3]) # apply rounding to get integer
print $0}  }
' OFS="\t" $scratchLoc$'Nelson_CAD_raw_in_N_eff' $scratchLoc$'Nelson_CAD_raw' > $scratchLoc$'Nelson_CAD_raw_imputedNeff'

head $scratchLoc$'Nelson_CAD_raw_imputedNeff'



# perform QC
remove_ambiguous_alleles $scratchLoc$'Nelson_CAD_raw_imputedNeff'
wc -l $scratchLoc$'Nelson_CAD_raw_imputedNeff_noAmbiguousAlleles' # 1127191


# QC
arguments='/home/mk907/scripts/R/LDpred2_QC.R '$baseLDpredLoc$'/ '$scratchLoc$'Nelson_CAD_raw_imputedNeff_noAmbiguousAlleles '$rawLocV4$'NelsonCAD_HM3_keep 1'
Rscript $arguments

# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $rawLocV4$'NelsonCAD_HM3_keep' $scratchLoc$'Nelson_CAD_raw_imputedNeff_noAmbiguousAlleles'  > $rawLocV4$'NelsonCAD_HM3'
head $rawLocV4$'NelsonCAD_HM3'
wc -l $rawLocV4$'NelsonCAD_HM3' # 1039344




# QC
# Nelson CAD: 
# average sample size: 324586
#awk '{ sum += $10; n++ } END { if (n > 0) print sum / n; }' /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/Nelson_CAD/28714975-GCST004787-EFO_0001645-Build37.f.tsv
#$munge_sumstats --sumstats /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/Nelson_CAD/28714975-GCST004787-EFO_0001645-Build37.f.tsv --N 324586 --snp variant_id --signed-sumstats logor,0 --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $NelsonCAD_sumstatsStem


# AAA v Nelson CAD: 
#$lsdc --rg $AAA_sumstatsStem.sumstats.gz,$NelsonCAD_sumstatsStem.sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $resultsLoc$'AAD_abdominal/rG_NelsonCAD'
#Genetic Correlation: 0.5805 (0.0926)
#Z-score: 6.2684
#P: 3.6473e-10



# calculate rho for Nelson CAD: assume that all controls in AAA overlap
# Main text: 
# "Association summary statistics (after λ correction) from the UKBB were combined with 
# the 1000 Genomes–imputed GWAS results2 and the Exome results3 via two separate fixed-effect 
# inverse-variance-weighted meta-analyses implemented in GWAMA24. 
# from Supplementary text, Supp figure 2, sample sizes:
#UKBB
#(HARD) 6482 / 137371
#(SOFT) 10801 / 137371
#1000G sumstats Meta-analysis:
#60801 / 123504
#Exome Meta-analysis
#42335 / 78240
# cases = 10801 + 60801 + 42335 = 113937
# controls = 137371 + 123504 + 78240 = 339115
# nkl0= 102026 # number of controls overlapping between studies
# nk1= 1068 #number of cases in study k (AAA)
# nk0= 102026 # number of controls in study k (AAA)
# nl1= 113937 # number of cases in study l
# nl0= 339115 # number of controls in study l
# shaPRS_rho(nkl0,nk1, nk0, nl1,nl0) # 0.02799685


# my AAA GWAS: 1068 cases and 102026 controls remaining after main filters.

##############################################
# III) AAD_RELATED:  $rawLocV4$'AAD_RELATED_hm3'
###############


# need to recreate a GWAS training set for AAA_RELATED, this is larger than the AAA GWAS, and also must exclude the controls from Corry's 25K 
# (because the CONTROLS were from the AAA Test set, but some of that may have been in the AAD_RELATED training set)
# 260956 

# create backup
# cp /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAD_RELATED_GWAS_TRAINING /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAD_RELATED_GWAS_TRAINING_backup
CreateBinaryGWASData_V4 $AAD_RELATED "AAD_RELATED" $rawLocV4 $AAD_ABDOMINAL_TEST_V4


# we still need remove the 25K control individuals from the AAD related training set
# this is because the Test set above is already excluded the 25K, so when we created the training set as "everyone -TEST", that TEST set did not have the 25K to be excluded
cp $rawLocV4$'AAD_RELATED_GWAS_TRAINING' $rawLocV4$'AAD_RELATED_GWAS_TRAINING_backup'

awk 'FNR == NR {  file1[ $1 ] = $1; next; }
FNR <= NR { if ($1 in file1 == 0) {print $0} }
' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry_noheader_indis' $rawLocV4$'AAD_RELATED_GWAS_TRAINING_backup' > $rawLocV4$'AAD_RELATED_GWAS_TRAINING'

wc -l $rawLocV4$'AAD_RELATED_GWAS_TRAINING_backup' # 260956
wc -l $rawLocV4$'AAD_RELATED_GWAS_TRAINING' # 235956
# so we excluded the 25K...


######

# III) GWAS:
rm -rf $scratchLocV4$'gwaslogs/'
mkdir -p $scratchLocV4$'gwaslogs/'
cd $scratchLocV4$'gwaslogs/' # so that all logs go here


for ((i=8; i<=$numChroms; i++)); do

# Run GWAS on each chrom, using the more lenient SNP QC (will filter post association)
pgen='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr'$i

pheno_name='AAD_RELATED'
mkdir -p $resultsLocV4$pheno_name$'/GWAS/'
arguments=' --pfile '$pgen$' --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --glm firth-fallback hide-covar cols=+err,+a1freq --covar-variance-standardize --mac 20 --keep '$rawLocV4$pheno_name$'_GWAS_TRAINING --pheno '$rawLocV4$pheno_name$'_all.phe --covar '$GWAS_COVARS$'_ARRAY --exclude '$scratchLoc$'ALL_FAIL_INFO_03 --extract '$hapmap3_SNPs$' --out '$resultsLocV4$pheno_name$'/GWAS/'$i$' --ci 0.95 --allow-extra-chr --snps-only --bp-space 1'
sbatch --mem ${PLINKRAM} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu -e $scratchLocV4$'gwaslogs/'$pheno_name$i$'.err' -o $scratchLocV4$'gwaslogs/'$pheno_name$i$'.out' --wrap "$plink2 $arguments"

done


# split the hapmap3 SNPs for the first 7 chroms into 3 for parallelization to prevent jobs getting killed
for ((i=1; i<=7; i++)); do

# split out the hapmap3 per chrom keeplists
awk -v i="$i" '{if($1 == i) print $2}' $hapmap3_b37bim > $hapmap3_SNPs$'_'$i

# split the per chrom into 3 files, with extension 00, 01 and 02
split -n l/3 -d $hapmap3_SNPs$'_'$i $hapmap3_SNPs$'_'$i$'_'

for ((j=0; j<=2; j++)); do

pgen='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr'$i

# AAD_RELATED
pheno_name='AAD_RELATED'
mkdir -p $resultsLocV4$pheno_name$'/GWAS/'
arguments=' --pfile '$pgen$' --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --glm firth-fallback hide-covar cols=+err,+a1freq --covar-variance-standardize --mac 20 --keep '$rawLocV4$pheno_name$'_GWAS_TRAINING --pheno '$rawLoc$pheno_name$'_all.phe --covar '$GWAS_COVARS$'_ARRAY --exclude '$scratchLoc$'ALL_FAIL_INFO_03 --extract '$hapmap3_SNPs$'_'$i$'_0'$j$' --out '$resultsLocV4$pheno_name$'/GWAS/'$i$'_'$j$' --ci 0.95 --allow-extra-chr --snps-only --bp-space 1'
sbatch --mem ${PLINKRAM} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu -e $scratchLocV4$'gwaslogs/'$pheno_name$i$'_0'$j$'.err' -o $scratchLocV4$'gwaslogs/'$pheno_name$i$'_0'$j$'.out' --wrap "$plink2 $arguments"



done # end of per split loop

done # end of chrom loop


# find out case/controls used in AAD_RELATED:
awk 'FNR == NR { file1[ $2] = $2; next; }
FNR <= NR {  
if( $2 in file1 ) {print $0} } 
' $rawLocV4$'AAD_RELATED_GWAS_TRAINING' $rawLoc$'AAD_RELATED_all.phe' > $rawLoc$'AAD_RELATED_casControl'
head $rawLoc$'AAD_RELATED_casControl'
wc -l $rawLoc$'AAD_RELATED_casControl' # 235955

awk '{count[$3]++} END {for (word in count) print word, count[word]}' $rawLoc$'AAD_RELATED_casControl'
# 1 102055
# 2 133900

# Get study overlap between AAD_RELATED and AAA
# nkl0= 102026 # number of controls overlapping between studies
# nk1= 1068 #number of cases in study k (AAA)
# nk0= 102026 # number of controls in study k (AAA)
# nl1= 133900 # number of cases in study l
# nl0= 102055 # number of controls in study l
# shaPRS_rho(nkl0,nk1, nk0, nl1,nl0) # 0.07666249

#########

# concat the GWAS results
pheno_name='AAD_RELATED'
(head -1 $resultsLocV4$pheno_name$'/GWAS/1_0'$plinkformat && tail -n +2 -q $resultsLocV4$pheno_name$'/GWAS/'*$plinkformat ) > $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat # the output must have different last name otherwise it would infinitely concat itself
head $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat 
wc -l $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat  # 1376913
# 0.920025

produce_manhattan_PLINK2 $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat $resultsLocV4$pheno_name$'/GWAS_HM3' $pheno_name$'_GWAS_HM3' 7.30103

N_eff=$(GetN_eff $rawLocV4$pheno_name$'_all.phe' $rawLocV4$pheno_name$'_GWAS_TRAINING')
ConvertPLINK2BinaryPhenoToSumtats $N_eff $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat $scratchLocV4$pheno_name$'_raw'
remove_ambiguous_alleles $scratchLocV4$pheno_name$'_raw'

# QC
arguments='/home/mk907/scripts/R/LDpred2_QC.R '$baseLDpredLoc$'/ '$scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles '$scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles_keep 1'
Rscript $arguments

# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles_keep' $scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles'  > $rawLocV4$pheno_name$'_hm3'
head $rawLocV4$pheno_name$'_hm3'
wc -l $rawLocV4$pheno_name$'_hm3' # 1048744


#####################################################################
# IV) Leicester AAA data: $rawLocV4$'Leicester_HM3'
#############################
# processing the received data 

# explore Corry's assoc files
cd '/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/corryFinal/'

# create backup
# cp UKAGS-UKB_GWAS_for_Cambridge.tar.gz UKAGS-UKB_GWAS_for_Cambridge.tar.gz_backup
tar zxvf UKAGS-UKB_GWAS_HRC_imputation.tar.gz

# check what columns we have
head $rawLocV4$'corryFinal/UKAGS-UKB_GWAS_HRC_imputation/GLM_cl1-3_CAM_PRS_HRC_3PCs.PHENO1.glm.logistic'
wc -l $rawLocV4$'corryFinal/UKAGS-UKB_GWAS_HRC_imputation/GLM_cl1-3_CAM_PRS_HRC_3PCs.PHENO1.glm.logistic' # 21,107,169



# match the MAFs to the association stats
awk 'FNR == NR { 
file1[ $2 ] = $5; next; }
FNR <= NR {  
if(FNR == 1) {print $0"\tMAF"} 
else if ($7 == "ADD" ) {
if( $3 in file1) {print $0"\t"file1[$3]} 
} }' $rawLocV4$'corryFinal/UKAGS-UKB_GWAS_HRC_imputation/GLM_cl1-3_CAM_PRS_HRC_noPCs.PHENO1.glm.logistic.MAF.afreq' $rawLocV4$'corryFinal/UKAGS-UKB_GWAS_HRC_imputation/GLM_cl1-3_CAM_PRS_HRC_3PCs.PHENO1.glm.logistic' > $rawLocV4$'corryFinal/sumstats_maf'
head $rawLocV4$'corryFinal/sumstats_maf'
wc -l $rawLocV4$'corryFinal/sumstats_maf' # 5276793

# calculate Effective sample size from cases/controls: https://privefl.github.io/bigsnpr-extdoc/polygenic-scores-pgs.html
# cases= 4777
# controls=10125
# N_eff = 4 / (1 / cases + 1 / controls) 
# round(N_eff) # 12983, 

# find all the matching ones and map
# #CHROM  POS     ID              REF     ALT     A1      TEST    OBS_CT  OR              LOG(OR)_SE      Z_STAT    P
# 1       752478  1:752478:G:A    G       A       A       ADD     14654   0.942417        0.0762627       -0.77767  0.436764
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
#  1     2     3    4     5        6         7    8   9   10
awk  -v numIndis=12983  'FNR == NR { 
file1[ $1"_"$4"_"$5"_"$6 ] = $2; file2[ $1"_"$4"_"$6"_"$5 ] = $2; next; }
FNR <= NR {  
if(FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
lookup1=$1"_"$2"_"$4"_"$5
if( lookup1 in file1 || lookup1 in file2) {

if( lookup1 in file1 ) {rsId = file1[lookup1]}
else if ( lookup1 in file2 ) {rsId = file2[lookup1]}

# find the non-effect allele, this could either be REF or ALT, as PLINK2 is shit and inconsistent
other_allele=$4 # default to REF
if(other_allele == $6) {other_allele = $5} # if REF is the same as A1, then use ALT

print $1"\t"$2"\t"rsId"\t"$6"\t"other_allele"\t"$13"\t"log($9)"\t"$10"\t"$12"\t"numIndis

} }' $hapmap3_b37bim $rawLocV4$'corryFinal/sumstats_maf' > $rawLocV4$'Corry_HM3'
head $rawLocV4$'Corry_HM3'
wc -l $rawLocV4$'Corry_HM3' # 766,970


# convert to ldsc format
$munge_sumstats --sumstats $rawLocV4$"AAD_abdominal_hm3" --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $rawLocV4$"AAD_abdominal_hm3_sums"

# convert to ldsc format
$munge_sumstats --sumstats $rawLocV4$"Corry_HM3" --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $rawLocV4$"Corry_HM3_sums"


# get rG
$lsdc --rg $rawLocV4$"AAD_abdominal_hm3_sums".sumstats.gz,$rawLocV4$"Corry_HM3_sums".sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $resultsLocV4$'/rG_UKBB_Corry'
# rg is 1.103 , IE over 1...

# perform QC
remove_ambiguous_alleles $rawLocV4$"Corry_HM3"
wc -l $rawLocV4$"Corry_HM3_noAmbiguousAlleles" # 709772

# QC
arguments='/home/mk907/scripts/R/LDpred2_QC.R '$baseLDpredLoc$'/ '$rawLocV4$'Corry_HM3_noAmbiguousAlleles '$rawLocV4$'Leicester_HM3_keep 1'
Rscript $arguments


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $rawLocV4$'Leicester_HM3_keep' $rawLocV4$'Corry_HM3_noAmbiguousAlleles'  > $rawLocV4$'Leicester_HM3'
head $rawLocV4$'Leicester_HM3'
wc -l $rawLocV4$'Leicester_HM3' 


##########################################################################################
# (V) AAAgen summary stats:$rawLocV4$'AAAGen_HM3'
##########################

cd /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAAgen/

# log into SFTP: (details redacted)


get -r PRS_weights_PRScs

get -r summary_statistics


cd /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAAgen/PRS_weights_PRScs/wo_UKBB
unzip AAAgen_wo_UKBB_PRScs_weights.zip

cd /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAAgen/PRS_weights_PRScs/wo_UKBB/AAAgen_wo_UKBB_PRScs_weights_share
cat out_pst_eff_a1_b0.5_phiauto_chr*.txt > ALL_PRSCS.txt
wc -l ALL_PRSCS.txt # 1,118,997, 1 mil SNPs

cd /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAAgen/summary_statistics/wo_UKBB
gunzip meta-AAAgen-wo_UKBB-sumstat.txt.gz



head meta-AAAgen-wo_UKBB-sumstat.txt
wc -l meta-AAAgen-wo_UKBB-sumstat.txt

#  participating studies: MVP, Geisinger, Decode
awk 'FNR == NR { 
file1[ $1"_"$4 ] = $1"_"$4; next; }
FNR <= NR {  
lookup1=$1"_"$2
if( lookup1 in file1) {print file1[lookup1]}  
}' $hapmap3_b37bim meta-AAAgen-wo_UKBB-sumstat.txt > "meta-AAAgen-wo_UKBB-sumstat.txt_hapmap3"
head "meta-AAAgen-wo_UKBB-sumstat.txt_hapmap3" # 1407351
wc -l "meta-AAAgen-wo_UKBB-sumstat.txt_hapmap3" # 1407351

# this 
#  1        2      3       4      5        6          7               8                 9           10
# chr     pos     ref     alt     N       AF      Direction       Effectsize      Effectsize_SD   pvalue
# to this
# chr	pos	SNP	A1	A2	Freq1.Hapmap	b	se	p	N
# get the SNP IDs and convert it to my format
awk 'FNR == NR { 
hap3_1[ $1"_"$4"_"$5"_"$6 ] = $2; hap3_1[ $1"_"$4"_"$6"_"$5 ] = $2; next; }
FNR <= NR {  
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { 
part1=$1"\t"$2
part2=$4"\t"$3"\t"$6"\t"$8"\t"$9"\t"$10"\t"$5
lookup1=$1"_"$2"_"$3"_"$4
lookup2=$1"_"$2"_"$4"_"$3
if( lookup1 in hap3_1) {print part1"\t"hap3_1[lookup1]"\t"part2}
else if( lookup2 in hap3_1) {print part1"\t"hap3_1[lookup2]"\t"part2}
 } }' $hapmap3_b37bim meta-AAAgen-wo_UKBB-sumstat.txt > $AAGENRAW$"AAAgenHM3"
head $AAGENRAW$"AAAgenHM3"
wc -l $AAGENRAW$"AAAgenHM3" # 1,398,937




## get rG between my sumstats and the AAAgen
# needed to redo the munge for my original UKBB association as, it for some reason did not recognise 70%+ of the SNPs in the PLINK association file...
# Removed 770966 variants that were not SNPs or were strand-ambiguous.
# Writing summary statistics for 1217311 SNPs (333605 with nonmissing beta)
myUKBBAAA="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/AAD_abdominal_hm3"
#numIndis=$(tail -n 1 '/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/AAD_abdominal_hm3' | awk '{ print $10}')
#echo $numIndis
$munge_sumstats --sumstats $myUKBBAAA --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $AAA_sumstatsStem

# convert to ldsc format
$munge_sumstats --sumstats $AAGENRAW$"AAAgenHM3" --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap --N-col N --chunksize 500000 --merge-alleles $w_hm3_snplist  --out $AAGENRAW$"AAAgenHM3_sums"
# Writing summary statistics for 1217311 SNPs (1114592 with nonmissing beta) 

# get rG between AAGen and my own AAA  sumstats
$lsdc --rg $AAA_sumstatsStem.sumstats.gz,$AAGENRAW$"AAAgenHM3_sums".sumstats.gz --ref-ld-chr $eur_w_ld_chr --w-ld-chr $eur_w_ld_chr --out $resultsLocV4$'/rG_UKBB_AAGen'
# Genetic Correlation: 0.8978 (0.0868)
# P: 4.2177e-25


# impute N_eff:
# extract data with signature: rsID, Beta, Beta_SE, AF
awk '{if(FNR ==1) {print "SNP\tBeta\tBeta_SE\tAF\tN"}  else {print $3"\t"$7"\t"$8"\t"$6"\t"$10} }' $AAGENRAW$"AAAgenHM3" > $AAGENRAW$"AAAgenHM3_in"
head $AAGENRAW$"AAAgenHM3_in"

# run R script 
Rscript /home/mk907/scripts/R/ImputeNeff.R $AAGENRAW$"AAAgenHM3_in" $AAGENRAW$"AAAgenHM3_N_eff"
head $AAGENRAW$"AAAgenHM3_N_eff"

# Find median imputed sample size
R
Neff=read.table("/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAAgen/AAAgenHM3_N_eff", header=T)
median(Neff$n_eff) # 105,871.1


# produce the final summary stats file
awk 'FNR == NR { 
file1[$1] = $6; next; } 
FNR <= NR { 
gsub("\r","",$3); 
if(FNR == 1) {print $0}
else if ( $3 in file1 ) { 
$10=  sprintf("%0.f", file1[$3]) # apply rounding to get integer
print $0}  }
' OFS="\t" $AAGENRAW$"AAAgenHM3_N_eff" $AAGENRAW$"AAAgenHM3" > $AAGENRAW$"AAAgenHM3_imputedNeff"

head $AAGENRAW$"AAAgenHM3_imputedNeff"



# perform QC
remove_ambiguous_alleles $AAGENRAW$"AAAgenHM3_imputedNeff"
wc -l $AAGENRAW$"AAAgenHM3_imputedNeff_noAmbiguousAlleles" # 1289585

# QC
arguments='/home/mk907/scripts/R/LDpred2_QC.R '$baseLDpredLoc$'/ '$AAGENRAW$'AAAgenHM3_imputedNeff_noAmbiguousAlleles '$rawLocV4$'AAAGen_HM3_keep 1'
Rscript $arguments


# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $rawLocV4$'AAAGen_HM3_keep' $AAGENRAW$'AAAgenHM3_imputedNeff_noAmbiguousAlleles'  > $rawLocV4$'AAAGen_HM3'
head $rawLocV4$'AAAGen_HM3'
wc -l $rawLocV4$'AAAGen_HM3'  # 997572





#################################################################
# (VI) Build final PRS models:
##############################
# Build AAAGen + NelsonCAD and AAAGen + MEGASTROKE via both shaPRS and MTAG

###############################
# A) Find out if adjunct datas add anything, and if they do, is shaPRS or MTAG better

#####################
# AAAGen + MEGASTROKE
phenam='AAAGen_MEGASTROKE'

# shaPRS # MTAG
shaPRS_fixed $rawLocV4$phenam $rawLocV4$'AAAGen_HM3' $rawLocV4$'MEGASTROKE_HM3' '0'
RUN_MTAG $rawLocV4$phenam $rawLocV4$'AAAGen_HM3' $rawLocV4$'MEGASTROKE_HM3' 


# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_mtag' $phenam"_mtag" $scratchLocV4
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4


# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_mtag' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4


# Sum and Evaluate
SumPRS_Agnostic $phenam$'_mtag_PRSCS' $dataLocV4
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_mtag_PRSCS_all.sscore' $resultsLocV4$phenam'_mtag' # correlation_sq 0.00536 /  AUC 0.695 // FULL: correlation_sq 0.00454 /  AUC 0.692
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' #  correlation_sq 0.00539 /  AUC 0.694 // FULL: 0.00453 /  AUC 0.69



#####################
# AAAGen + NelsonCAD
phenam='AAAGen_NelsonCAD' # NEEDS REEVAL

# shaPRS # MTAG
shaPRS_fixed $rawLocV4$phenam $rawLocV4$'AAAGen_HM3' $rawLocV4$'NelsonCAD_HM3' '0'
RUN_MTAG $rawLocV4$phenam $rawLocV4$'AAAGen_HM3' $rawLocV4$'NelsonCAD_HM3' 


# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_mtag' $phenam"_mtag" $scratchLocV4
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_mtag' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4
  
# Sum and Evaluate
SumPRS_Agnostic $phenam$'_mtag_PRSCS' $dataLocV4
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_mtag_PRSCS_all.sscore' $resultsLocV4$phenam'_mtag' # correlation_sq 0.00576 / AUC 0.703
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' # correlation_sq 0.00589 /  AUC 0.704

#################################
# AAAGen + NelsonCAD + MEGASTROKE
phenam='AAAGen_MEGASTROKE_NelsonCAD'  # NEED REEVAL

# shaPRS
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"AAAGen_NelsonCAD_shaprs" $rawLocV4$'MEGASTROKE_HM3' '0'

# MTAG via 3 way
# this assumes that the _MTAG conversions already exist, as per above
outputDir=$rawLocV4$phenam
/usr/local/software/master/miniconda/2/bin/python '/home/mk907/software/mtag/mtag/mtag.py' --sumstats $rawLocV4$'AAAGen_HM3_MTAG',$rawLocV4$'NelsonCAD_HM3_MTAG',$rawLocV4$'MEGASTROKE_HM3_MTAG'  --out $outputDir --n_min 0.0 --force

awk '{if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}
else{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12"\t"$7} }' $outputDir$"_trait_1.txt" > $outputDir$"_mtag"


# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_mtag' $phenam"_mtag" $scratchLocV4
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4


# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_mtag' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4
  
# Sum and Evaluate
SumPRS_Agnostic $phenam$'_mtag_PRSCS' $dataLocV4
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_mtag_PRSCS_all.sscore' $resultsLocV4$phenam'_mtag' # correlation_sq 0.00571 /  AUC 0.7 
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' # correlation_sq 0.00613 /  AUC 0.706 CI 

# looks like they both add to accuracy!
######################################

#######################
# B) Build PRSCS baselines for all 3 AAA sumstats to see how accurate they are

##########################
phenam="AAD_abdominal_hm3"
PRSCS_allChroms $rawLocV4$phenam $phenam $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4

# Sum and Evaluate
SumPRS_Agnostic $phenam$'_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_PRSCS_all.sscore' $resultsLocV4$phenam # correlation_sq 0.000276  AUC 0.549 // FULL: correlation_sq 0.000298 /  AUC 0.548


######################
phenam="Leicester_HM3"
PRSCS_allChroms $rawLocV4$phenam $phenam $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4

# Sum and Evaluate
SumPRS_Agnostic $phenam$'_PRSCS' $dataLocV4
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_PRSCS_all.sscore' $resultsLocV4$phenam #  correlation_sq 0.00222 AUC 0.625 // FULL:  correlation_sq 0.00185 / AUC 0.622


###################
phenam="AAAGen_HM3"
#PRSCS_allChroms $rawLocV4$phenam $phenam $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4

# Sum and Evaluate
SumPRS_Agnostic $phenam$'_PRSCS' $dataLocV4
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_PRSCS_all.sscore' $resultsLocV4$phenam # correlation_sq 0.00518  AUC 0.693 // FULL: correlation_sq 0.00439  /  AUC 0.69



#################################
# C) Get PRS for just the AAAs: UKBB+AAAgen+Leicester
# when merging with Leicester, need to add back the SNPs that were missing from Corry's work
# this is the same for both MTAG and shaPRS, so we need to merge Corry separately for MTAG too (IE cannot do it all in 4)

# merge against AAAGen as 'trait 1' as that has the most signal
# shaPRS # MTAG, just for AAAGen + Leicester
phenam='AAAGen_Corry'
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"AAAGen_HM3" $rawLocV4$'Leicester_HM3' '0'
RUN_MTAG $rawLocV4$phenam $rawLocV4$"AAAGen_HM3" $rawLocV4$'Leicester_HM3'

# shaPRS # MTAG, just for AAAGen_Leicester + UKBB
phenam='AAAGen_Leicester_UKBB'
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"AAAGen_Corry_shaprs" $rawLocV4$'AAD_abdominal_hm3' '0'
RUN_MTAG $rawLocV4$phenam $rawLocV4$"AAAGen_Corry_mtag" $rawLocV4$'AAD_abdominal_hm3'

# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_mtag' $phenam"_mtag" $scratchLocV4
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_mtag' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4



# Sum and Evaluate
SumPRS_Agnostic $phenam$'_mtag_PRSCS' $dataLocV4
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_mtag_PRSCS_all.sscore' $resultsLocV4$phenam'_mtag' # correlation_sq 0.0056 / AUC 0.698 // FULL correlation_sq: 0.00473 / AUC 0.695
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' #  correlation_sq 0.00591 /  AUC 0.705 // FULL correlation_sq 0.00499 / AUC 0.701


# B) merge against UKBB AAA as 'trait 1'
# shaPRS # MTAG, just for UKBB + Leicester
phenam='UKBB_Corry'
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"AAD_abdominal_hm3" $rawLocV4$'Leicester_HM3' '0'
RUN_MTAG $rawLocV4$phenam $rawLocV4$"AAD_abdominal_hm3" $rawLocV4$'Leicester_HM3'

# shaPRS # MTAG, just for UKBB_Leicester + AAAGen
phenam='UKBB_Leicester_AAAGen'
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"UKBB_Corry_shaprs" $rawLocV4$'AAAGen_HM3' '0'
RUN_MTAG $rawLocV4$phenam $rawLocV4$"UKBB_Corry_mtag" $rawLocV4$'AAAGen_HM3'

# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_mtag' $phenam"_mtag" $scratchLocV4
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_mtag' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4

# Sum and Evaluate
SumPRS_Agnostic $phenam$'_mtag_PRSCS' $dataLocV4 
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_mtag_PRSCS_all.sscore' $resultsLocV4$phenam'_mtag' # correlation_sq 0.0047 / AUC 0.682  // FULL: correlation_sq 0.00396 /  AUC 0.678 
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' # correlation_sq 0.00581 / AUC 0.703 // FULL: correlation_sq 0.00491 /  AUC 0.7
# - so this is worse
# It seems like there is an advantage at starting from the larger sample sized/stronger signal data




#############################
# D) FINAL PRS: Combine all AAA with the adjunct datasets
###############

# need to create _MTAG formatted sumstats for these two ,as these haven't been produced yet
awk '{if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7/$8"\t"$9"\t"$10} 
}'  $rawLocV4$'AAD_RELATED_hm3' > $rawLocV4$'AAD_RELATED_hm3_MTAG'

awk '{if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7/$8"\t"$9"\t"$10} 
}'  $rawLocV4$'AAAGen_Corry_mtag' > $rawLocV4$'AAAGen_Corry_mtag_MTAG'


###########################
# UKBB+AAAgen+Corry+RELATED

phenam='AAAGen_Leicester_UKBB_RELATED'

# work out the overlap for AAAD_RELATED, UKBB and NelsonCAD
# overlap between UKBB and AAD_RELATED_hm3 is 0.07666249, however, UKBB is now part of a much larger sumstats
# Find out exactly how big 
#Ns=$(awk '{ if(FNR > 1) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawLocV4"AAAGen_Leicester_UKBB_shaprs")
#echo $Ns # 119459
#Ns=$(awk '{ if(FNR > 1) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawLocV4"AAD_abdominal_hm3")
#echo $Ns # 4227
# 4227 / (4227 + 119459) = 0.03417525 # so our UKBB AAA cohort is only about 3% of the total so far
# 0.03417525 * 0.07666249 = 0.00261996 # which makes the total overlap 0.2%, ie ~0, so it is OK to round this to 0
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"AAAGen_Leicester_UKBB_shaprs" $rawLocV4$'AAD_RELATED_hm3' '0'



#$rawLocV4$"AAAGen_Leicester_UKBB_shaprs"

# MTAG via 4 way
# this assumes that the _MTAG conversions already exist, as per above
outputDir=$rawLocV4$phenam
/usr/local/software/master/miniconda/2/bin/python '/home/mk907/software/mtag/mtag/mtag.py' --sumstats $rawLocV4$'AAAGen_Corry_mtag_MTAG',$rawLocV4$'AAD_abdominal_hm3_MTAG',$rawLocV4$'AAD_RELATED_hm3_MTAG'  --out $outputDir --n_min 0.0 --force

awk '{if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}
else{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12"\t"$7} }' $outputDir$"_trait_1.txt" > $outputDir$"_mtag"


# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_mtag' $phenam"_mtag" $scratchLocV4
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4 


# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_mtag' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 


# Sum and Evaluate
SumPRS_Agnostic $phenam$'_mtag_PRSCS' $dataLocV4 
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_mtag_PRSCS_all.sscore' $resultsLocV4$phenam'_mtag' #  correlation_sq 0.0061 /  AUC 0.707     // FULL  correlation_sq 0.00509  / AUC 0.702
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' #  correlation_sq 0.00658 / AUC 0.714 // FULL   correlation_sq 0.00539 / AUC 0.708




################################################
# UKBB+AAAgen+Corry+Related+NelsonCAD
phenam='ALL_AAA_NelsonCAD'
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"AAAGen_Leicester_UKBB_RELATED_shaprs" $rawLocV4$'NelsonCAD_HM3' '0'

# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4



# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 # called on NON_GWAS_TEST set


# Sum and Evaluate
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' # FULL correlation_sq 0.0057 / AUC 0.712 
# so this is better than also adding in the MEGASTROKE
# make this the new champion, ie create an LDpred2 version
## Check the best PRS beats PRSCS if built in LDpred2
LDpred2_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_LDpred2 $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 # called on NON_GWAS_TEST set
# 833,480 SNPs

SumPRS_Agnostic $phenam$'_shaprs_LDpred2' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_LDpred2_all.sscore' $resultsLocV4$phenam'_shaprsEvaluate_PRS' # correlation_sq 0.00724 / AUC 0.722 / FULL test correlation_sq 0.00598 /  AUC 0.716

# check effective sample size
Ns=$(awk '{ if(FNR > 1) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawLocV4$phenam"_shaprs")
echo $Ns # 441,496

# associated SNPs
Ns=$(awk '{ if(FNR > 1 && $9 < 5 * 10e-8) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawLocV4$phenam"_shaprs")
echo $Ns # 333,371

Ns=$(awk '{ if(FNR > 1 && $9 < 5 * 10e-16) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawLocV4$phenam"_shaprs")
echo $Ns # 268,075

Ns=$(awk '{ if(FNR > 1 && $9 < 5 * 10e-16) {sum += 1; n++ } } END { printf("%0.f", sum) ; }' $rawLocV4$phenam"_shaprs")
echo $Ns # 584


###############################

# UKBB+AAAgen+Corry+Related+MEGASTROKE
phenam='ALL_AAA_MEGASTROKE'
shaPRS_fixed $rawLocV4$phenam $rawLocV4$"AAAGen_Leicester_UKBB_RELATED_shaprs" $rawLocV4$'MEGASTROKE_HM3' '0'


# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4


# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 # called on NON_GWAS_TEST set


# Sum and Evaluate
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' # FULL  correlation_sq 0.00548 / AUC 0.709



################################################
# UKBB+AAAgen+Corry+Related+NelsonCAD+MEGASTROKE

phenam='ALL_AAA_MEGASTROKE_NelsonCAD'
shaPRS_fixed $rawLocV4$'ALL_AAA_MEGASTROKE' $rawLocV4$"AAAGen_Leicester_UKBB_RELATED_shaprs" $rawLocV4$'MEGASTROKE_HM3' '0'

shaPRS_fixed $rawLocV4$phenam $rawLocV4$"ALL_AAA_MEGASTROKE_shaprs" $rawLocV4$'NelsonCAD_HM3' '0'


outputDir=$rawLocV4$phenam
/usr/local/software/master/miniconda/2/bin/python '/home/mk907/software/mtag/mtag/mtag.py' --sumstats $rawLocV4$'AAAGen_Corry_mtag_MTAG',$rawLocV4$'AAD_abdominal_hm3_MTAG',$rawLocV4$'AAD_RELATED_hm3_MTAG',$rawLocV4$'MEGASTROKE_HM3_MTAG',$rawLocV4$'NelsonCAD_HM3_MTAG'  --out $outputDir --n_min 0.0 --force

awk '{if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}
else{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12"\t"$7} }' $outputDir$"_trait_1.txt" > $outputDir$"_mtag"

# PRS-CS (wait until cluster job finishes)
PRSCS_allChroms $rawLocV4$phenam'_mtag' $phenam"_mtag" $scratchLocV4
PRSCS_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_PRSCS $phenam'_mtag' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 # called on NON_GWAS_TEST set
BuildPRS_allChroms_PRSCS $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 # called on NON_GWAS_TEST set

# Sum and Evaluate
SumPRS_Agnostic $phenam$'_mtag_PRSCS' $dataLocV4 
SumPRS_Agnostic $phenam$'_shaprs_PRSCS' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_mtag_PRSCS_all.sscore' $resultsLocV4$phenam'_mtag' # correlation_sq 0.00657 / AUC 0.713   / FULL TEST correlation_sq 0.00546 /  AUC 0.708
Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_PRSCS_all.sscore' $resultsLocV4$phenam'_shaprs' #  correlation_sq 0.00692 / AUC 0.719 / FULL correlation_sq 0.00561 /  AUC 0.711 


#####################
## Check the best PRS beats PRSCS if built in LDpred2
LDpred2_allChroms $rawLocV4$phenam'_shaprs' $phenam"_shaprs" $scratchLocV4

# Build PRS profile scores  (wait until cluster job finishes)
BuildPRS_allChroms_LDpred2 $phenam'_shaprs' $scratchLocV4 $AAD_ABDOMINAL_TEST_V4 $dataLocV4 # called on NON_GWAS_TEST set
# 833,480 SNPs

SumPRS_Agnostic $phenam$'_shaprs_LDpred2' $dataLocV4

Evaluate_PRS $dataLocV4$'PRSProfiles/'$phenam$'_shaprs_LDpred2_all.sscore' $resultsLocV4$phenam'_shaprsEvaluate_PRS' # correlation_sq 0.00724 / AUC 0.722 / FULL test correlation_sq 0.00589 /  AUC 0.715



