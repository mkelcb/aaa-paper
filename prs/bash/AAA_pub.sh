
###############################################

# AAA project main script

# depends on reusable functions and variables loaded into the session from AAA_functions.sh
###############################################

##############################################

# (I) Data pre-processing:
#######################

# 1. Get target phenotype from UKBB
# extraction scripts only work from here
# this is now moved to /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/deprecated/endpoints/output
cd /rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/
pheno_name='AAD_abdominal'
export_UKBB_pheno_nodupes2 $pheno_name $endpointLoc$'output/'$pheno_name'/'$pheno_name$'_definition.txt' $rawLoc
#0 500126
#1 2362
#0 500125
#1 2362

# also get 'AAD any, which are 
pheno_name='AAD_any'
export_UKBB_pheno_nodupes2 $pheno_name $endpointLoc$'output/'$pheno_name'/'$pheno_name$'_definition.txt' $rawLoc
#0 497559
#1 4929


# b) genetically possibly overlapping: take all the ICD10 codes from the excel sheet and export them out both fatal-non-fatal
pheno_name='AAD_related'
export_UKBB_pheno_nodupes2 $pheno_name $endpointLoc$'output/'$pheno_name'/'$pheno_name$'_definition.txt' $rawLoc
# 0 323071
# 1 179417

awk '{ if (FNR >1 && $2 == 1) {print $1} }' $rawLoc$pheno_name$'_all' > $AAD_RELATED
head $AAD_RELATED
wc -l $AAD_RELATED # 179417

#####
# 2. Get the HAPMAP3 SNPs
awk '{print $2}' $hapmap3_b37bim > $hapmap3_SNPs
head $hapmap3_SNPs
wc -l $hapmap3_SNPs # 1403851


##############################################################################
# (II) Extract Conventional predictors:
#######################################

# 1)  Conventional Predictors dataframe:  This is used only for the GWAS
# from paper: Polygenic risk scores in cardiovascular risk prediction: A cohort study and modelling analyses
# age at baseline, sex, smoking status, history of diabetes, systolic blood pressure, total cholesterol, and high-density lipoprotein (HDL) cholesterol, and also, BMI

# trying to map the above covariates to the CEU curated dataset:
# R (need to change the 'v3' manually here to current version)
R
library(readstata13)
library('data.table')
ceu <- read.dta13("/rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/20210302/STATA/analysis.dta")
setDT(ceu)
vars = attributes(ceu)$var.labels
names(vars) = names(ceu)
myData=ceu[,.(eid=idno, hypdbin)]
write.table(vars, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/scratch/CEU_covars.txt", row.names = T, col.names = T, quote = FALSE)

# I manually identified the following entries from: /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/scratch/CEU_covars.txt
# idno Study-specific subject ID
# ages Age at survey (yrs)
# sex Sex
# bmi BMI (kg/m2)
# smallbin Smoking status  // , decided not to use this, as this would make heavy ex smokers same as never smoked, instead decided to use
# smallstat   // which has Ex, Current and Never
# hxdiabbin History of diabetes
# dbp DBP (mmHg)
# sbp SBP (mmHg)
# tchol Total cholesterol (mmol/l)
# hdl HDL-C (mmol/l)
# ldl LDL-C (mmol/l)

myData=ceu[,.(eid=idno, age=ages, sex, bmi, smoking=smallstat, lipidlowering=lipdbin, hypertension=hypdbin, diabetes=hxdiabbin, dbp, sbp, total_cholesterol=tchol, HDL_cholesterol=hdl, LDL=ldl)]
write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/COVS_raw", row.names = F, col.names = T, quote = T)
quit()

head $rawLoc$'COVS_raw'


# create covariate DF for GWAS, should include CHIP now, as we split across chips
# y ~ age + sex + batch + batch 10 PCs
# eid     genotyping.array        Batch   Plate.Name      Well    Cluster.CR     dQC      Internal.Pico..ng.uL.   Submitted.Gender        Inferred.Gender   X.intensity     Y.intensity     Submitted.Plate.Name    Submitted.Well  sample.qc.missing.rate  heterozygosity  heterozygosity.pc.corrected     het.missing.outliers   putative.sex.chromosome.aneuploidy       in.kinship.table        excluded.from.kinship.inference      excess.relatives        in.white.British.ancestry.subset       used.in.pca.calculation  PC1     PC2     PC3     PC4     PC5     PC6     PC7    PC8      PC9     PC10    PC11    PC12    PC13    PC14    PC15    PC16    PC17   PC18     PC19    PC20    PC21    PC22    PC23    PC24    PC25    PC26    PC27   PC28     PC29    PC30    PC31    PC32    PC33    PC34    PC35    PC36    PC37   PC38     PC39    PC40    in.Phasing.Input.chr1_22        in.Phasing.Input.chrX  in.Phasing.Input.chrXY
#  1             2                  3        4             5           6          7              8                   9                          10               11              12               13                     14                   15                  16                  17                             18                         19                                 20                                 21                        22                           23                                    24             25     
# awk '{count[$2]++} END {for (word in count) print word, count[word]}' $QCfile

awk 'FNR == NR { 
file1[ $1 ] = $3"_"$2"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34
next; }
FNR <= NR {  
if(FNR == 1) {print "FID\tIID\tage\tsex\tbatch_chip\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10"} else {
gsub(/"/, "", $1);
gsub(/"/, "", $2);
gsub(/"/, "", $3);
if( $1 in file1) {print $1"\t"$1"\t"$2"\t"$3"\t"file1[ $1 ] } }}
' $QCfile $rawLoc$'COVS_raw' > $GWAS_COVARS$'_ARRAY'
head $GWAS_COVARS$'_ARRAY'
wc -l $GWAS_COVARS$'_ARRAY' # 487279


# 2) Generate new traditional risk factors design matrix that includes additional predictors:
# extract Stroke + family history of CHD and Alcohol  raw/ alcallbin, hxchd_p, hxstroke_p
# This is used for the final survival analysis

R
library(readstata13)
library('data.table')
ceu <- read.dta13("/rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/20210901/STATA/analysis.dta")
setDT(ceu)
vars = attributes(ceu)$var.labels
names(vars) = names(ceu)
myData=ceu[,.(eid=idno, hypdbin)]
write.table(vars, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/scratch/CEU_covars.txt", row.names = T, col.names = T, quote = FALSE)

myData=ceu[,.(eid=idno, age=ages, sex, bmi, smoking=smallstat, lipidlowering=lipdbin, hypertension=hypdbin, diabetes=hxdiabbin, dbp, sbp, total_cholesterol=tchol, HDL_cholesterol=hdl, LDL=ldl, Alcohol_consumption=alcallbin, Family_CHD=hxchd_p, Family_Stroke=hxstroke_p)]
write.table(myData, "/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/COVS_raw", row.names = F, col.names = T, quote = T, sep ="\t")
quit()

head $rawLocV4$'COVS_raw'

# PROBLEM, if we do NOT specify that field separator is tab (\t) but leave it at space, then we will get misalignment of the columns:
# (as some fields have spaces in them)

# check what kinds of values we have for CHD and Stroke
awk '{count[$15]++} END {for (word in count) print word, count[word]}' FS="\t" $rawLocV4$'COVS_raw'
#NA 52901
#"No" 250007
#"Yes" 199551
#"Family_CHD" 1


awk '{count[$15]++} END {for (word in count) print word, count[word]}' FS=" " $rawLocV4$'COVS_raw'
#NA 49335
#"No" 239085
#"Yes" 188470
#"Family_CHD" 1
#"Other" 4179
#"Current" 21390

#   1     2     3     4       5            6             7              8       9     10           11               12            13              14                    15            16
# "eid" "age" "sex" "bmi" "smoking" "lipidlowering" "hypertension" "diabetes" "dbp" "sbp" "total_cholesterol" "HDL_cholesterol" "LDL"	    "Alcohol_consumption" "Family_CHD" "Family_Stroke"
# "1000012" 58.2529983520508 "Female" 41.6100006103516 "Ex" "Other" "Other" "Other" 96 160 5.71199989318848 1.12899994850159 3.8510000705719  "Current"              "No"             "No"

# create union of CHD + Stroke # END {print "number of CVD cases "count}
awk '{
if(FNR == 1) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t\"Family_CVD\""}
else 
{
CVD=$15
if($16 == "\"Yes\"") {CVD = "\"Yes\""; count++}

print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"CVD
}
} 
' FS="\t" $rawLocV4$'COVS_raw' > $rawLocV4$'COVS_raw_CVD'
head $rawLocV4$'COVS_raw_CVD'
#I have replaced spaces with Tabs above to avoid misaligned columns


########################

# I) GWAS Sample QC:
# head $QCfile
# eid     genotyping.array        Batch   Plate.Name      Well    Cluster.CR     dQC      Internal.Pico..ng.uL.   Submitted.Gender        Inferred.Gender   X.intensity     Y.intensity     Submitted.Plate.Name    Submitted.Well  sample.qc.missing.rate  heterozygosity  heterozygosity.pc.corrected     het.missing.outliers   putative.sex.chromosome.aneuploidy       in.kinship.table        excluded.from.kinship.inference      excess.relatives        in.white.British.ancestry.subset       used.in.pca.calculation  PC1     PC2     PC3     PC4     PC5     PC6     PC7    PC8      PC9     PC10    PC11    PC12    PC13    PC14    PC15    PC16    PC17   PC18     PC19    PC20    PC21    PC22    PC23    PC24    PC25    PC26    PC27   PC28     PC29    PC30    PC31    PC32    PC33    PC34    PC35    PC36    PC37   PC38     PC39    PC40    in.Phasing.Input.chr1_22        in.Phasing.Input.chrX  in.Phasing.Input.chrXY
#  1             2                  3        4             5           6          7              8                   9                          10               11              12               13                     14                   15                  16                  17                             18                         19                                 20                                 21                        22                           23                                    24             25     


# create covariate DF for GWAS
# y ~ age + sex + batch + batch 10 PCs

awk 'FNR == NR { 
file1[ $1 ] = $3"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34
next; }
FNR <= NR {  
if(FNR == 1) {print "FID\tIID\tage\tsex\tbatch\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10"} else {
gsub(/"/, "", $1);
gsub(/"/, "", $2);
gsub(/"/, "", $3);
if( $1 in file1) {print $1"\t"$1"\t"$2"\t"$3"\t"file1[ $1 ] } }}
' $QCfile $rawLoc$'COVS_raw' > $GWAS_COVARS
head $GWAS_COVARS
wc -l $GWAS_COVARS # 487279


# II) SNP QC:

# a) exclude SNPs that failed QC in the genotype dataset

# this is the 800K chip  
# 805427

head -n 1 $GENOTYPE_QC
#rs_id affymetrix_snp_id affymetrix_probeset_id chromosome position allele1_ref allele2_alt strand array Batch_b001_qc Batch_b002_qc Batch_b003_qc Batch_b004_qc Batch_b005_qc Batch_b006_qc Batch_b007_qc Batch_b008_qc Batch_b009_qc Batch_b010_qc Batch_b011_qc Batch_b012_qc Batch_b013_qc Batch_b014_qc Batch_b015_qc Batch_b016_qc Batch_b017_qc Batch_b018_qc Batch_b019_qc Batch_b020_qc Batch_b021_qc Batch_b022_qc Batch_b023_qc Batch_b024_qc Batch_b025_qc Batch_b026_qc Batch_b027_qc Batch_b028_qc Batch_b029_qc Batch_b030_qc Batch_b031_qc Batch_b032_qc Batch_b033_qc Batch_b034_qc Batch_b035_qc Batch_b036_qc Batch_b037_qc Batch_b038_qc Batch_b039_qc Batch_b040_qc Batch_b041_qc Batch_b042_qc Batch_b043_qc Batch_b044_qc Batch_b045_qc Batch_b046_qc Batch_b047_qc Batch_b048_qc Batch_b049_qc Batch_b050_qc Batch_b051_qc Batch_b052_qc Batch_b053_qc Batch_b054_qc Batch_b055_qc Batch_b056_qc Batch_b057_qc Batch_b058_qc Batch_b059_qc Batch_b060_qc Batch_b061_qc Batch_b062_qc Batch_b063_qc Batch_b064_qc Batch_b065_qc Batch_b066_qc Batch_b067_qc Batch_b068_qc Batch_b069_qc Batch_b070_qc Batch_b071_qc Batch_b072_qc Batch_b073_qc Batch_b074_qc Batch_b075_qc Batch_b076_qc Batch_b077_qc Batch_b078_qc Batch_b079_qc Batch_b080_qc Batch_b081_qc Batch_b082_qc Batch_b083_qc Batch_b084_qc Batch_b085_qc Batch_b086_qc Batch_b087_qc Batch_b088_qc Batch_b089_qc Batch_b090_qc Batch_b091_qc Batch_b092_qc Batch_b093_qc Batch_b094_qc Batch_b095_qc UKBiLEVEAX_b1_qc UKBiLEVEAX_b2_qc UKBiLEVEAX_b3_qc UKBiLEVEAX_b4_qc UKBiLEVEAX_b5_qc UKBiLEVEAX_b6_qc UKBiLEVEAX_b7_qc UKBiLEVEAX_b8_qc UKBiLEVEAX_b9_qc UKBiLEVEAX_b10_qc UKBiLEVEAX_b11_qc in_HetMiss in_Relatedness in_PCA PC1_loading PC2_loading PC3_loading PC4_loading PC5_loading PC6_loading PC7_loading PC8_loading PC9_loading PC10_loading PC11_loading PC12_loading PC13_loading PC14_loading PC15_loading PC16_loading PC17_loading PC18_loading PC9_loading PC20_loading PC21_loading PC22_loading PC23_loading PC24_loading PC25_loading PC26_loading PC27_loading PC28_loading PC9_loading PC30_loading PC31_loading PC32_loading PC33_loading PC34_loading PC35_loading PC36_loading PC37_loading PC38_loading PC9_loading PC40_loading in_Phasing_Input
#rs28659788 Affx-13546538 AX-32115783 1 723307 C G + 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA 0

# find SNPs that failed in any of the batches
awk '{
for (i =10; i<=115 ; i++)  {
 if( $i == "0") 
 { print $1; next; }
 } }' $GENOTYPE_QC > $scratchLoc$'GENOTYPE_FAILED_QC'
 
head $scratchLoc$'GENOTYPE_FAILED_QC'
wc -l $scratchLoc$'GENOTYPE_FAILED_QC' # 51,707 , 51K SNPs failed QC in at least 1 of the batches

# UKBB SNP QC file format: $14 = MAF, $18 = INFO
#     1               2              3             4          5       6      7              8               9               10             11                 12                  13                    14                   15           16             17          18
# alternate_ids     rsid          chromosome    position   alleleA alleleB comment   HW_exact_p_value HW_lrt_p_value   alleleA_count   alleleB_count   alleleA_frequency   alleleB_frequency     minor_allele_frequency minor_allele   major_allele     info     impute_info   missing_proportion A B AA AB BB NULL total
# 21:9411239_G_A    rs559462325      21        9411239     G          A       NA            1           0.97576         974788         29.9137             0.999969              3.06865e-05          3.06865e-05          A             G            0.260815     0.260815       5.97113e-16 0 0 487379 29.8824 0.0156863 2.91038e-10 487409

rm -rf $scratchLoc$'ALL_FAIL_INFO_03'
for ((i=1; i<=$numChroms; i++)); do

chromQCFile='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/reference_files/ukb_impv3_chr'$i$'_snpstats.txt'
# b) exclude SNPs that have low imputation quality INFO < 0.3

# need to skip lines with comments, otherwise awk will default to non-numeric comparisons
awk '/^[^#]/ { print $0 }' $chromQCFile | awk '{ if ($18 < 0.3 ) {print $2} }'  > $scratchLoc$i$'_FAIL_INFO_03'

cat $scratchLoc$i$'_FAIL_INFO_03' >>  $scratchLoc$'ALL_FAIL_INFO_03'


done

# also add in the genotype failed SNPs
cat $scratchLoc$'GENOTYPE_FAILED_QC' >>  $scratchLoc$'ALL_FAIL_INFO_03'

# find out how many we are going to exclude
wc -l $scratchLoc$'ALL_FAIL_INFO_03' # 15,393,551


# create phenotype files ( need the +1 , as 1 = control, 2 = case
awk 'FNR == NR { 
phenos[$1] = $2
next; }
FNR <= NR {  
if(FNR == 1) {print "FID\tIID\tphe"} else {
if( $2 in phenos) {print $1"\t"$2"\t"phenos[$2]+1 }} }
' $rawLoc$'AAD_abdominal_all' $GWAS_COVARS > $rawLoc$'AAD_abdominal_all.phe'



# B) create lists for subsetting:
#######################


# 2) white-british: in.white.British.ancestry.subset, 0 for $23
awk '{ if ($23 == 0) {print $1} }' $QCfile > $NON_EUR
head $NON_EUR
wc -l $NON_EUR # 78,436


# 3) sex discordant or low quality genotype ($9: Submitted.Gender        $10: Inferred.Gender, or $18:het.missing.outliers)
awk '{ if ($9 != $10 || $18 == 1) {print $1} }' $QCfile > $SEX_DISC_or_LQ
head $SEX_DISC_or_LQ
wc -l $SEX_DISC_or_LQ # 374


# 4) 12K used for metaGRS
awk '{ if (FNR >1) {print $1} }' '/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/Stroke_metaGRS_training_samples.txt' > $META_GRS
head $META_GRS
wc -l $META_GRS # 11979


# 5) withdrawals
awk '{ if (FNR >1) {print $1} }' '/rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/Withdrawals/20210221_withdrawals.csv' > $WITHDRAWALS
head $WITHDRAWALS
wc -l $WITHDRAWALS # 176



# 6) Get list of cases for AAD any ( to be excluded from controls)
AAD_ANY_CASES=$rawLoc$'AAD_any_cases'
awk '{ if (FNR >1 && $2 == 1) {print $1} }' $rawLoc$'AAD_any_all' > $AAD_ANY_CASES
head $AAD_ANY_CASES
wc -l $AAD_ANY_CASES # 4866


# 7) relatedness: too closely related: 0.1875
# ID1     ID2     HetHet  IBS0    Kinship
# 1000083 5972935 0.075   0.0058  0.2353

# get all pairs whose Kinship > 0.0884 (first of second degree relatives)
awk '{ if (FNR >1 && $5 > 0.0884) {print $1"\t"$2} }' '/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/reference_files/kinship_relatedness.txt' > $RELATEDS
head $RELATEDS
wc -l $RELATEDS # 40,227

#- exclude based on relatedness (prefer to exclude controls)
awk 'FNR == NR { 
cases[ $1 ] = $1;
next; 
}
FNR <= NR {  
# if first indi is a case we exclude the other
# if second indi is a case we exclude the first
# if both are cases we just exclude the first
if( $1 in cases ) {print $2 }
else if( $2 in cases ) {print $1 }
else {print $1 }
}
' $AAD_ANY_CASES $RELATEDS > $CLOSE_REATEDS

head $CLOSE_REATEDS
wc -l $CLOSE_REATEDS # 40227

# find out how many we excluded for each category
awk 'FNR == NR { 
cases[ $1 ] = $1;
next; 
}
FNR <= NR {  
if( $1 in cases && $2 in cases ) {print "BOTH_INDI_CASES" }
else if( $1 in cases || $2 in cases ) {print "ONE_INDI_CASE" }
else {print "BOTH_INDI_CONTROL" }
}
' $AAD_ANY_CASES $RELATEDS > $scratchLoc$'howManyRelatedCases'
awk '{count[$1]++} END {for (word in count) print word, count[word]}' $scratchLoc$'howManyRelatedCases'

# ONE_INDI_CASE 753
# BOTH_INDI_CONTROL 39464
# BOTH_INDI_CASES 10
# we were able to keep an extra 761 cases in our study by preferentially excluding controls




# 9) create master exclusion lists:
# exclusion list for both cases & controls: 
cat $NON_EUR $SEX_DISC_or_LQ $META_GRS $WITHDRAWALS $CLOSE_REATEDS > $scratchLoc$'MASTER_EXCLUDE_CASECONTROL_dupes'
wc -l  $scratchLoc$'MASTER_EXCLUDE_CASECONTROL_dupes' # 131192


awk '!visited[$0]++' $scratchLoc$'MASTER_EXCLUDE_CASECONTROL_dupes' > $MASTER_EXCLUDE_CASECONTROL
wc -l $MASTER_EXCLUDE_CASECONTROL # 121879


# exclusion list for controls
#a) for GWAS we exclude the lipid/blood pressure ones
cat $MASTER_EXCLUDE_CASECONTROL $CVD $AAD_RELATED $AAD_ANY_CASES $LIPID_ANTIHYPERTENSIVE >  $scratchLoc$'MASTER_EXCLUDE_CONTROLS_GWAS_dupes'
wc -l  $scratchLoc$'MASTER_EXCLUDE_CONTROLS_GWAS_dupes' # 368108  505541

awk '!visited[$0]++' $scratchLoc$'MASTER_EXCLUDE_CONTROLS_dupes' > $MASTER_EXCLUDE_CONTROLS_GWAS
wc -l $MASTER_EXCLUDE_CONTROLS_GWAS # 287826

# b) for the incident model, we keep the hypertensives, as that variable is in the model, which would create an unbalanced scenario
cat $MASTER_EXCLUDE_CASECONTROL $CVD $AAD_RELATED $AAD_ANY_CASES  >  $scratchLoc$'MASTER_EXCLUDE_CONTROLS_dupes'
wc -l  $scratchLoc$'MASTER_EXCLUDE_CONTROLS_dupes' # 368108 

awk '!visited[$0]++' $scratchLoc$'MASTER_EXCLUDE_CONTROLS_dupes' > $MASTER_EXCLUDE_CONTROLS
wc -l $MASTER_EXCLUDE_CONTROLS # 287826


#######################


# 10) create master keep lists for Cases and Controls for both GWAS and survival analysis


# All Controls: this is the same for all 3
# Controls:   everyone - MASTER_EXCLUDE_CONTROLS

# a) GWAS
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR {  
if( $1 in file1 == 0  ) {print $1 }
}
' $MASTER_EXCLUDE_CONTROLS_GWAS $QCfile > $CONTROLS_ALL_GWAS
head $CONTROLS_ALL_GWAS
wc -l $CONTROLS_ALL_GWAS # 206617

# b) incident
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR {  
if( $1 in file1 == 0  ) {print $1 }
}
' $MASTER_EXCLUDE_CONTROLS $QCfile > $CONTROLS_ALL
head $CONTROLS_ALL
wc -l $CONTROLS_ALL # 229543



# All Cases:


# AAD_abdominal
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if (FNR >1 && $2 == 1 && $1 in file1 == 0) {print $1} }
' $MASTER_EXCLUDE_CASECONTROL $rawLoc$'AAD_abdominal_all' > $AAD_ABDOMINAL_CASES_QC
head $AAD_ABDOMINAL_CASES_QC
wc -l $AAD_ABDOMINAL_CASES_QC # , 2009 ( 700 if we excluded the anti hypertensives )



##############
# Full Cohort = Cases + Controls

# AAD_abdominal
#a) GWAS
cat $CONTROLS_ALL_GWAS $AAD_ABDOMINAL_CASES_QC > $AAD_ABDOMINAL_FULL_QC_GWAS
wc -l $AAD_ABDOMINAL_FULL_QC_GWAS # 208626

#b) incident
cat $CONTROLS_ALL $AAD_ABDOMINAL_CASES_QC > $AAD_ABDOMINAL_FULL_QC
wc -l $AAD_ABDOMINAL_FULL_QC # 231552

##########################################################

# split the UKBB, such that the Test set won't have any of the interim UKBB, so that we could use summary stats that were built using the interim UKBB
# (as the Interim included all BELIEVE people too, this also solves the problem of the test set being biased due to that study)

# create 'FID IID' list  halfUKBB=251244 # = 502488/2
awk '{if(FNR > 1) {gsub(/"/, "", $1); print  $1}}' $rawLoc$'COVS_raw' > $scratchLocV4$'all_indis'
head $scratchLocV4$'all_indis'
wc -l $scratchLocV4$'all_indis'


# take out from total the 150K interim
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 == 0) {print $0} }
' $rawLoc$'interim_indis_pubIDs' $scratchLocV4$'all_indis' > $scratchLocV4$'all_indis_no_interim'
head $scratchLocV4$'all_indis_no_interim'
wc -l $scratchLocV4$'all_indis_no_interim' # 350478

# randomly pick (total/2)-150k out again from this
arguments='/home/mk907/scripts/R/randomNums.R 350478 42 '$scratchLocV4$'randomNums'
Rscript $arguments

# 99187 = 251244-152057
awk '{if (FNR <= 99187) print $0}' $scratchLocV4$'randomNums' > $scratchLocV4$'randomNums_indis'
head $scratchLocV4$'randomNums_indis'

# pick the actual sample ids of the 99187 non Interims (these still include the non-QC people)
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $1 } }
' $scratchLocV4$'randomNums_indis' $scratchLocV4$'all_indis_no_interim' > $scratchLocV4$'random_nonInterims' 
head $scratchLocV4$'random_nonInterims' 
wc -l $scratchLocV4$'random_nonInterims' 


cat $scratchLocV4$'random_nonInterims' $rawLoc$'interim_indis_pubIDs' > $AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE
head $AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE
wc -l $AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE # 251244, we got for training, that includes the interims and a random half of the rest, (still includes the non-QC people)


# TRAINING:
#a) GWAS
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1) {print $1"\t"$1} }
' $AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE $AAD_ABDOMINAL_FULL_QC_GWAS > $AAD_ABDOMINAL_TRAINING_GWAS_V4
head $AAD_ABDOMINAL_TRAINING_GWAS_V4
wc -l $AAD_ABDOMINAL_TRAINING_GWAS_V4 # 103144 ->  189173

#  wc -l $AAD_ABDOMINAL_FULL_QC_GWAS # 208626


#b) incident: IE the traditional model
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 ) {print $1"\t"$1} }
' $AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE $AAD_ABDOMINAL_FULL_QC > $AAD_ABDOMINAL_TRAINING_V4
head $AAD_ABDOMINAL_TRAINING_V4
wc -l $AAD_ABDOMINAL_TRAINING_V4 # 114561 - > 209720



# TEST SET: all training eligible who are NOT in the GWAS set
# Note: these include the 25K controls that will be excluded in the next stage

#a) GWAS
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 == 0) {print $1"\t"$1} }
' $AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE $AAD_ABDOMINAL_FULL_QC_GWAS > $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K'
head $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K'
wc -l $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K' # 105482 -> 19453




#b) incident: IE the traditional model
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 == 0) {print $1"\t"$1} }
' $AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE $AAD_ABDOMINAL_FULL_QC > $AAD_ABDOMINAL_TEST_V4$'_extra25K'
head $AAD_ABDOMINAL_TEST_V4$'_extra25K'
wc -l $AAD_ABDOMINAL_TEST_V4$'_extra25K' # 116991 ->  21832


################################
# University of Leicester collaboration:  choose 25K random individuals from the above that are controls
# cases = 5000
# controls = 25000
#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 16667
# which will then need to be excluded from our training set (for AAD this is OK, as they come from the test set anyway, but for the AAD_RELATED, they need to be excluded from the training set too)
# $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K'

# subset to controls only
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1) {print $1} }
' $CONTROLS_ALL_GWAS $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K' > $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls'
wc -l $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls' # 104561
head $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls'


# random 25K indices
numItems=$(wc -l < "$AAD_ABDOMINAL_TEST_GWAS_V4"_extra25K_onlyControls)
mySeed=42
arguments='/home/mk907/scripts/R/randomNums.R '$numItems$' '$mySeed$' '$AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls_randomNums'
Rscript $arguments

head -n 25000 $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls_randomNums' > $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls_randomNums_25k'

# get IIDs
awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $1 } }
' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls_randomNums_25k' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControls' > $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry' 
head $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry'
wc -l $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry'

# map these indices back to the .sample file
# grab these indis in the original order that also has the sex column ( will be used later to verify match from Corry)
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1) {print $0} }
' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry' /rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/ukb_adiposity_imp_v3.sample > $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry_noheader_indis'
head $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry_noheader_indis'
wc -l $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry_noheader_indis'


# get the row indices: this is the file to be sent to Corry
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1) {print FNR} }
' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry' /rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/ukb_adiposity_imp_v3.sample > $AAD_ABDOMINAL_TEST_GWAS_V4$'_rows'
head $AAD_ABDOMINAL_TEST_GWAS_V4$'_rows'
wc -l $AAD_ABDOMINAL_TEST_GWAS_V4$'_rows'


awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if (FNR in file1) {print $0} }
' $AAD_ABDOMINAL_TEST_GWAS_V4$'_rows' /rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/ukb_adiposity_imp_v3.sample > $AAD_ABDOMINAL_TEST_GWAS_V4$'_CorryTest'
head $AAD_ABDOMINAL_TEST_GWAS_V4$'_CorryTest'
wc -l $AAD_ABDOMINAL_TEST_GWAS_V4$'_CorryTest'


# remove these individuals from our TEST sets above
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 == 0) {print $0} }
' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry_noheader_indis' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K' > $AAD_ABDOMINAL_TEST_GWAS_V4
head $AAD_ABDOMINAL_TEST_GWAS_V4
wc -l $AAD_ABDOMINAL_TEST_GWAS_V4 # 80482


# find out cases and controls in the GWAS test set
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1) {print $0} }
' $AAD_ABDOMINAL_TEST_GWAS_V4 $rawLoc$'AAD_abdominal_all.phe' > $scratchLocV4$'TEST_SET_GWAS_pheno'
head $scratchLocV4$'TEST_SET_GWAS_pheno'
wc -l $scratchLocV4$'TEST_SET_GWAS_pheno' # 80380
awk '{count[$3]++} END {for (word in count) print word, count[word]}' $scratchLocV4$'TEST_SET_GWAS_pheno'
#1 79511
#2 869


awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 == 0) {print $0} }
' $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry_noheader_indis' $AAD_ABDOMINAL_TEST_V4$'_extra25K' > $AAD_ABDOMINAL_TEST_V4
head $AAD_ABDOMINAL_TEST_V4
wc -l $AAD_ABDOMINAL_TEST_V4 # 91991

# find out cases and controls in test set
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1) {print $0} }
' $AAD_ABDOMINAL_TEST_V4 $rawLoc$'AAD_abdominal_all.phe' > $scratchLocV4$'TEST_SET_pheno'

head $scratchLocV4$'TEST_SET_pheno'
wc -l $scratchLocV4$'TEST_SET_pheno' # 91889
awk '{count[$3]++} END {for (word in count) print word, count[word]}' $scratchLocV4$'TEST_SET_pheno'

#cases = 869
#controls = 91020
#N_eff = 4 / (1 / cases + 1 / controls) 
#round(N_eff) # 3443


# ############################
# Example sctipt for CORRY
# awk 'FNR == NR { 
# file1[ $1 ] = $1;
# next; 
# }
# FNR <= NR { 
# if (FNR in file1) {print $0} }
# ' 487411_indices.txt <487411_FILE_WITH_2_HEADERLINES_v3.sample> > CONTROL_KEEPLIST.txt

# Sanity check
#awk '{if(FNR == 102131) print $0}' /rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/ukb_adiposity_imp_v3.sample


# # Check results back from Corry if order of SEX vars matches
# paste -d " " $AAD_ABDOMINAL_TEST_GWAS_V4$'_extra25K_onlyControlsCorry_noheader_indis' /rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/data/raw/CONTROL_KEEPLIST_sex_only_2.txt > $AAD_ABDOMINAL_TEST_GWAS_V4$'_FROMCORRY'
# head $AAD_ABDOMINAL_TEST_GWAS_V4$'_FROMCORRY'
# wc -l $AAD_ABDOMINAL_TEST_GWAS_V4$'_FROMCORRY'

# awk '{if ($4 != $5) {mismatch++} } END {print "number of mismatches: "mismatch}'  $AAD_ABDOMINAL_TEST_GWAS_V4$'_FROMCORRY'

# awk '{if ($4 == $5) {matched++} } END {print "number of matches: "matched}' $AAD_ABDOMINAL_TEST_GWAS_V4$'_FROMCORRY'
# # they all matched!!
# ############################


# create dataframe
pheno_name='AAD_abdominal'
head $rawLoc$pheno_name$'_survival'
# create 
awk 'FNR == NR { 
file1[ "\""$1"\"" ] = "\""$2"\"""\t""\""$3"\"""\t""\""$4"\"""\t""\""$5"\"""\t""\""$6"\"";
next; }
FNR <= NR {  
if(FNR == 1) {print "\"visit\"\t\"prev_followup\"\t\"prev_event\"\t\"inci_followup\"\t\"inci_event\"\t"$0} else {
if( $1 in file1) {print file1[ $1 ]"\t"$0 } }}
' $rawLoc$pheno_name$'_survival' $rawLocV4$'COVS_raw_CVD' > $scratchLocV4$pheno_name$'raw_COVS'
head $scratchLocV4$pheno_name$'raw_COVS'
wc -l $scratchLocV4$pheno_name$'raw_COVS' # 502489



awk 'FNR == NR { file1[ "\""$1"\"" ] = $1; next; }
FNR <= NR {  
if(FNR == 1 || $6 in file1 ) {print $0} } 
' $AAD_ABDOMINAL_TRAINING_V4 $scratchLocV4$pheno_name$'raw_COVS' > $AAD_ABDOMINAL_TRAINING_V4$'_COVS'
head $AAD_ABDOMINAL_TRAINING_V4$'_COVS'
wc -l $AAD_ABDOMINAL_TRAINING_V4$'_COVS'  # 114526


awk 'FNR == NR { file1[ "\""$1"\"" ] = $1; next; }
FNR <= NR {  
if(FNR == 1 || $6 in file1 ) {print $0} } 
' $AAD_ABDOMINAL_TEST_V4 $scratchLocV4$pheno_name$'raw_COVS' > $AAD_ABDOMINAL_TEST_V4$'_COVS'
head $AAD_ABDOMINAL_TEST_V4$'_COVS'
wc -l $AAD_ABDOMINAL_TEST_V4$'_COVS'  # 116919 -> 91934 (as we excluded the 25K)





##############################################
# III) GWAS: create GWAS of 
rm -rf $scratchLocV4$'gwaslogs/'
mkdir -p $scratchLocV4$'gwaslogs/'
cd $scratchLocV4$'gwaslogs/' # so that all logs go here


for ((i=8; i<=$numChroms; i++)); do

# Run GWAS on each chrom, using the more lenient SNP QC (will filter post association)
pgen='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr'$i

pheno_name='AAD_abdominal'
mkdir -p $resultsLocV4$pheno_name$'/GWAS/'
arguments=' --pfile '$pgen$' --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --glm firth-fallback hide-covar cols=+err,+a1freq --covar-variance-standardize --mac 20 --keep '$AAD_ABDOMINAL_TRAINING_GWAS_V4$' --pheno '$rawLoc$pheno_name$'_all.phe --covar '$GWAS_COVARS$'_ARRAY --exclude '$scratchLoc$'ALL_FAIL_INFO_03 --extract '$hapmap3_SNPs$' --out '$resultsLocV4$pheno_name$'/GWAS/'$i$' --ci 0.95 --allow-extra-chr --snps-only --bp-space 1'
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

pheno_name='AAD_abdominal'
mkdir -p $resultsLocV4$pheno_name$'/GWAS/'
arguments=' --pfile '$pgen$' --memory  '$PLINKRAM$' --threads '$SLURM_CPUS_ON_NODE$' --glm firth-fallback hide-covar cols=+err,+a1freq --covar-variance-standardize --mac 20 --keep '$AAD_ABDOMINAL_TRAINING_GWAS_V4$' --pheno '$rawLoc$pheno_name$'_all.phe --covar '$GWAS_COVARS$'_ARRAY --exclude '$scratchLoc$'ALL_FAIL_INFO_03 --extract '$hapmap3_SNPs$'_'$i$'_0'$j$' --out '$resultsLocV4$pheno_name$'/GWAS/'$i$'_'$j$' --ci 0.95 --allow-extra-chr --snps-only --bp-space 1'
sbatch --mem ${PLINKRAM} --cpus-per-task ${SLURM_CPUS_ON_NODE} --time 12:0:0 --partition cclake --account danesh-sl3-cpu -e $scratchLocV4$'gwaslogs/'$pheno_name$i$'_0'$j$'.err' -o $scratchLocV4$'gwaslogs/'$pheno_name$i$'_0'$j$'.out' --wrap "$plink2 $arguments"

# 1068 cases and 102026 controls remaining after main filters.

done # end of per split loop

done # end of chrom loop

# Concat assoc results & produce Manhattan plot
# merge partial and full chrom results : https://unix.stackexchange.com/questions/60577/concatenate-multiple-files-with-same-header
# note: this will work for partial chroms too, however, the ordering of the final files won't be per chrom, but alphabetical, eg: 2_1, 21, 2_2
# WARNING: occasionally jobs fail, which will then cause  chromosome not to have all SNPs GWAS-d, but PLINK will write partial results to file, so it won't be immediately obvious that job was killed prematurely
# this will cause problems: not all SNPs being covered, and the manhattan plot will fail, as the results will have partially completed lines, which when concated will cause errors (IE a SNP having another SNP as its P-value etc)
# to check for this, make sure that all "3.log" files have a non-zero size (as when a job gets killed this file will be 0 length)
# Manhattan
pheno_name='AAD_abdominal'
(head -1 $resultsLocV4$pheno_name$'/GWAS/1_0'$plinkformat && tail -n +2 -q $resultsLocV4$pheno_name$'/GWAS/'*$plinkformat ) > $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat # the output must have different last name otherwise it would infinitely concat itself

head $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat 
wc -l $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat  # 1351129

produce_manhattan_PLINK2 $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat $resultsLocV4$pheno_name$'/GWAS_HM3' $pheno_name$'_GWAS_HM3' 7.30103


# generate sumstats in my format
pheno_name="AAD_abdominal"
N_eff=$(GetN_eff $rawLoc$pheno_name$'_all.phe' $AAD_ABDOMINAL_TRAINING_GWAS_V4) # 4227
ConvertPLINK2BinaryPhenoToSumtats $N_eff $resultsLocV4$pheno_name$'/GWAS_HM3'$plinkformat $scratchLocV4$pheno_name$'_raw'
remove_ambiguous_alleles $scratchLocV4$pheno_name$'_raw'

# QC
arguments='/home/mk907/scripts/R/LDpred2_QC.R '$baseLDpredLoc$'/ '$scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles '$scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles_keep 1'
Rscript $arguments

# subset the SNPs to only keep the QC passed ones
awk 'FNR == NR { sums[ $1 ] = $1; next; } FNR <= NR { if( FNR == 1 || $3 in sums  ) {print $0 } }
' $scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles_keep' $scratchLocV4$pheno_name$'_raw_noAmbiguousAlleles'  > $rawLocV4$pheno_name$'_hm3'
head $rawLocV4$pheno_name$'_hm3'
wc -l $rawLocV4$pheno_name$'_hm3'  # 1049328











