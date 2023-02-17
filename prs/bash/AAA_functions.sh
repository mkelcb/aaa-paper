#############################

# screen -r -D 28731.pts-80.login-e-15
# reusable functions and variables used by other scripts
#####################################



ukbb_raw_pheno='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/output/AAD/AAD_events_and_followup.txt'
homeScratch='/rds/user/mk907/hpc-work/scratch/'
endpointLoc='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/'
baseLoc='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v3/'
dataLoc=$baseLoc$'data/'
rawLoc=$dataLoc$'raw/'
scratchLoc=$baseLoc$'scratch/'
resultsLoc=$baseLoc$'results/'
plink="/home/mk907/software/plink/plink"
plink2='/home/mk907/software/plink2/plink2'
numIndis=19318
shaPRS="_shaPRS"
NCORES_LDPred2=16
AAGENRAW="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/raw/AAAgen/"
QCfile='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/reference_files/ukb_sqc_v2.txt' # eid matched to our app
# signature:
#  1              2                 3          4            5         6            7               8                   9                        10                  11              12                 13                     14                 15                  16                    17                       18                              19                              20                                           21          22                              23                                24 
# eid     genotyping.array        Batch   Plate.Name      Well    Cluster.CR      dQC     Internal.Pico..ng.uL.   Submitted.Gender        Inferred.Gender       X.intensity     Y.intensity     Submitted.Plate.Name    Submitted.Well  sample.qc.missing.rate  heterozygosity  heterozygosity.pc.corrected  het.missing.outliers     putative.sex.chromosome.aneuploidy      in.kinship.table        excluded.from.kinship.inference excess.relatives        in.white.British.ancestry.subset      used.in.pca.calculation PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9     PC10    PC11    PC12 PC13     PC14    PC15    PC16    PC17    PC18    PC19    PC20    PC21    PC22    PC23    PC24    PC25    PC26    PC27    PC28    PC29    PC30    PC31 PC32     PC33    PC34    PC35    PC36    PC37    PC38    PC39    PC40    in.Phasing.Input.chr1_22        in.Phasing.Input.chrX   in.Phasing.Input.chrXY


NelsonCAD_sumstatsStem=$scratchLoc$'NelsonCAD_sumstatsStem'
MegaStroke_AS_sumstatsStem=$scratchLoc$'MegaStroke_AS_munge'


NON_EUR=$rawLoc$'NON_EUR' # those who are NOT European ancestry
SEX_DISC_or_LQ=$rawLoc$'SEX_DISC_or_LQ'  # Sex discordant or those with low quality genotypes ( too many missing or excess heterozygosity)
META_GRS=$rawLoc$'META_GRS' # the 12K training set for the Inouye lab's meta GRS project
WITHDRAWALS=$rawLoc$'WITHDRAWALS' # people no longer in the UKBB
CVD=$rawLoc$'CVD' # those who have some other CVD already
AAD_RELATED=$rawLoc$'AAD_RELATED' # those who have a condition that may be genetically overlapping with AAD
TEST_SET=$rawLoc$'EUR_subset' # those that are in the test set (the UKBiLEVE chip)
RELATEDS=$rawLoc$'RELATEDS' # pairs of people who are too closely related 
CLOSE_REATEDS=$rawLoc$'CLOSE_REATEDS' # the individual to exclude from the above (preferentially controls)




#############

MASTER_EXCLUDE_CASECONTROL=$rawLoc$'MASTER_EXCLUDE_CASECONTROL' # master list of everyone to be exluded from both case and control group # 159070
LIPID_ANTIHYPERTENSIVE=$rawLoc$'LIPID_ANTIHYPERTENSIVE' # list of people on hypertensive / lipid lowering drugs
MASTER_EXCLUDE_CONTROLS=$rawLoc$'MASTER_EXCLUDE_CONTROLS' # master list of everyone to be excluded from the controls # 200099

CONTROLS_ALL=$rawLoc$'CONTROLS_ALL' # controls for all 3 phenos
AAD_ABDOMINAL_CASES_QC=$rawLoc$'AAD_abdominal_cases_QC'
AAD_ABDOMINAL_FULL_QC=$rawLoc$'AAD_ABDOMINAL_FULL_QC'

AAD_ABDOMINAL_TRAINING=$rawLoc$'AAD_ABDOMINAL_TRAINING'
AAD_ABDOMINAL_TEST=$rawLoc$'AAD_ABDOMINAL_TEST'


CONTROLS_ALL_GWAS=$rawLoc$'CONTROLS_ALL_GWAS' # controls for all 3 phenos for GWAS
MASTER_EXCLUDE_CONTROLS_GWAS=$rawLoc$'MASTER_EXCLUDE_CONTROLS_GWAS' # master list of everyone to be excluded from the controls for the GWAS

AAD_ABDOMINAL_FULL_QC_GWAS=$rawLoc$'AAD_ABDOMINAL_FULL_QC_GWAS'

AAD_ABDOMINAL_TRAINING_GWAS=$rawLoc$'AAD_ABDOMINAL_TRAINING_GWAS'
AAD_ABDOMINAL_TEST_GWAS=$rawLoc$'AAD_ABDOMINAL_TEST_GWAS'

GWAS_COVARS=$rawLoc$'GWAS_COVARS'

#############
# v4
baseLocV4='/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/'
dataLocV4=$baseLocV4$'data/'
rawLocV4=$dataLocV4$'raw/'
scratchLocV4=$baseLocV4$'scratch/'
resultsLocV4=$baseLocV4$'results/'

baseLDpredLoc=$rawLocV4$'ldpred2/'
prscsrefs=$dataLocV4$'prscsrefs/'
prscsscript=$dataLocV4$'prscsscript/'

prscsxdir=$prscsscript$'PRScsx/'

UKBB_PLINK1=$dataLoc$'UKBB_PLINK1/'

AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE=$rawLoc$'AAD_ABDOMINAL_ALL_TRAINING_ELIGIBLE'
AAD_ABDOMINAL_TRAINING_GWAS_V4=$rawLoc$'AAD_ABDOMINAL_TRAINING_GWAS_V4'
AAD_ABDOMINAL_TRAINING_V4=$rawLoc$'AAD_ABDOMINAL_TRAINING_V4'
AAD_ABDOMINAL_TEST_GWAS_V4=$rawLoc$'AAD_ABDOMINAL_TEST_GWAS_V4'
AAD_ABDOMINAL_TEST_V4=$rawLoc$'AAD_ABDOMINAL_TEST_V4'







#################

hapmap3_b37bim=$dataLoc$'raw/hapmap3_r1_b37_fwd_consensus.qc.poly.recode.bim'
hapmap3_SNPs=$dataLoc$'raw/hapmap3_SNPs'
numChroms=22
plinkformat='.phe.glm.logistic.hybrid'
plinkformatlinear=".phe.glm.linear"

GENOTYPE_QC='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/COMMON/pre_qc_data/affy_ukbiobank_array/raw_data/showcase_release_19Jul/snpqc/ukb_snp_qc.txt'

shaPRSscriptLoc="/home/mk907/scripts/R/shaPRS.R"



ldscDir='/home/mk907/software/ldsc/'
w_hm3_snplist='/home/mk907/software/ldsc/w_hm3.snplist'
eur_w_ld_chr='/home/mk907/software/ldsc/eur_w_ld_chr/'
lsdc=$ldscDir$'ldsc/ldsc.py'
munge_sumstats=$ldscDir$'ldsc/munge_sumstats.py'

rhoNelsonCAD=0.02799685
rho_AAD_RELATED=0.07666249

numChroms=22

PLINKRAM=97000
SLURM_CPUS_ON_NODE=32

module load R/4.0.3
mtagscript='/home/mk907/software/mtag/'

# check for rG 
module load miniconda/2
source activate ldsc
# pip install pandas --user
# pip install scipy --user
# pip install bitarray --user

mkdir -p $rawLoc
mkdir -p $scratchLoc
mkdir -p $resultsLoc



mkdir -p $mtagscript

###############################

# install MTAG

# get the mtag repository:
cd $mtagscript

git clone https://github.com/omeed-maghzian/mtag.git
cd mtag

pip install joblib --user


# run the test
python mtag.py -h



# run tutorial: https://github.com/JonJala/mtag/wiki/Tutorial-1:-The-Basics
wget http://ssgac.org/documents/1_OA2016_hm3samp_NEUR.txt.gz

wget http://ssgac.org/documents/1_OA2016_hm3samp_SWB.txt.gz

gunzip 1_OA2016_hm3samp_NEUR.txt.gz
gunzip 1_OA2016_hm3samp_SWB.txt.gz

python $mtagscript$'mtag/mtag.py' --sumstats 1_OA2016_hm3samp_NEUR.txt,1_OA2016_hm3samp_SWB.txt --out ./tutorial_results_1.1NS --n_min 0.0 --stream_stdout &


# same thing with --no_overlap
python $mtagscript$'mtag/mtag.py' --sumstats 1_OA2016_hm3samp_NEUR.txt,1_OA2016_hm3samp_SWB.txt --no_overlap --out ./tutorial_results_1.1NS --n_min 0.0 --stream_stdout &


#####################################
outputDir=$rawLocV4$phenam
pheA_sumstats=$rawLocV4$'AAAGen_HM3'
pheB_sumstats=$rawLocV4$'NelsonCAD_HM3' 
# maps my sumstats format to the format needed for MTAG, and performs MTAG, output is  $outputDir$"_mtag" (# --use_beta_se does not work in MTAG for now...)
function RUN_MTAG { 
outputDir=$1
pheA_sumstats=$2
pheB_sumstats=$3
# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
#10      100146895       rs2296438       T       C       0.6544  0.0024  0.0082  0.7694  964057
# to MTAG
#   1             2        3              4       5        6             7           8           9             
# snpid           chr     bpos            a1      a2      freq           z          pval          n
# make sure to avoid division by zero, otherwise awk may truncate file and thus create incorrect result
awk '{
if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {
print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7/$8"\t"$9"\t"$10} 
}'  $pheA_sumstats > $pheA_sumstats$'_MTAG'


awk '{
if (FNR == 1) {print "snpid\tchr\tbpos\ta1\ta2\tfreq\tz\tpval\tn" }
is_7_numeric = $7 + 0 == $7
is_8_numeric = $8 + 0 == $8
if(is_7_numeric && is_8_numeric) {
print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7/$8"\t"$9"\t"$10} 
}'  $pheB_sumstats > $pheB_sumstats$'_MTAG'

#head $pheA_sumstats$'_MTAG'

/usr/local/software/master/miniconda/2/bin/python '/home/mk907/software/mtag/mtag/mtag.py' --sumstats $pheA_sumstats$'_MTAG',$pheB_sumstats$'_MTAG'  --out $outputDir --n_min 0.0 --force > $outputDir$'.log'


# convert the MTAG output into my own format 
# MTAG: 
#  1               2       3       4      5       6               7        8          9                     10                       11                     12
# SNP             CHR     BP      A1      A2      Z               N       FRQ     mtag_beta               mtag_se                 mtag_z                  mtag_pval
# rs4040617       1       779322  G       A       0.314132        10874   0.1134  0.0059077854624592615   0.030391777030177466    0.19438762848888813     0.8458723762003009

# map my sumstats format:
# 1       2        3      4      5            6          7        8      9       10
#chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N

awk '{if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN"}
else{print $2"\t"$3"\t"$1"\t"$4"\t"$5"\t"$8"\t"$9"\t"$10"\t"$12"\t"$7} }' $outputDir$"_trait_1.txt" > $outputDir$"_mtagMerged"

# in case there were SNPs that were not present in trait 2 these would get removed by the above, we use trait 1's values
awk 'FNR == NR { file1[ $3 ] = $0; next;  }
FNR <= NR { if ($3 in file1) {print file1[$3] } else {print $0} }
' $outputDir$'_mtagMerged' $pheASumstats > $outputDir$'_mtag'


}
export -f RUN_MTAG # this makes local functions executable when bsubbed



# exports single pheno from multiple codes, also filters out duplicates
function export_UKBB_pheno_nodupes2 { 
pheno_name=$1
definition=$2
rawLoc=$3

# load the right version of R that has all the packages
module load R/4.0.3

# extraction script only works from this location
#endpointLoc='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/endpoints/'

endpointLoc='/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/deprecated/endpoints/'

cd $endpointLoc

mkdir -p $endpointLoc$'output/'$pheno_name

# run script to extract pheno
Rscript src/curate_endpoints.R $definition

ukbb_raw_pheno=$endpointLoc$'output/'$pheno_name'/'$pheno_name$'_events_and_followup.txt'
# signature:
#  1        2          3                  4                  5               6
# eid     visit   prev_followup        prev_event      inci_followup     inci_event
# 5401266 0       -15.9069130732375       0             11.523613963039   0

# The above no longer works as Scott changes the script, will need to update the below to accomodate the
# NEW signature:

#  1        2             3              4             5                              6                                   7                                 8                       9                10            11               12              13
# eid     visit   date_assessment age_assessment  assessment_centre       assessment_centre_nation        hospital_records_nation_at_max_followup prev_followupprev_date        prev_event      inci_followup   inci_date       inci_event      inci_fatal      all_cause_mortality     lost_to_followup        linkage_withdrawn


# filter duplicates, we want the baseline entry for each individual only IE, where visit == 0, (otherwise the inci_followup will be from their current visit which is not the total)
awk '
FNR <= NR {  
if( FNR == 1 || $2 == 0 ) {print $0 } 
}
' $ukbb_raw_pheno > $rawLoc$pheno_name$'_survival'

original=$(wc -l $ukbb_raw_pheno)
noDupes=$(wc -l $rawLoc$pheno_name$'_survival')
echo "before/after filtering for dupes: "$original$" / "$noDupes

outExtract=$endpointLoc$'output/'$pheno_name'/'$pheno_name$'_filtered'
awk '{ if(FNR == 1) {print "eid\tpheno"} else { if ($4 == "1" || $6 == "1" ) { print $1"\t1" } else {print $1"\t0"} } }' $rawLoc$pheno_name$'_survival' > $rawLoc$pheno_name$'_all'
awk '{count[$2]++} END {for (word in count) print word, count[word]}' $rawLoc$pheno_name$'_all'

# export out the incident/prevalent phenos separately too
outExtract=$endpointLoc$'output/'$pheno_name'/'$pheno_name$'_filtered_prevalent'
awk '{ if(FNR == 1) {print "eid\tpheno"} else { if ($4 == "1") { print $1"\t1" } else {print $1"\t0"} } }' $rawLoc$pheno_name$'_survival' > $rawLoc$pheno_name$'_prevalent'


outExtract=$endpointLoc$'output/'$pheno_name'/'$pheno_name$'_filtered_incident'
awk '{ if(FNR == 1) {print "eid\tpheno"} else { if ( $6 == "1" ) { print $1"\t1" } else {print $1"\t0"} } }' $rawLoc$pheno_name$'_survival' > $rawLoc$pheno_name$'_incident'


}
export -f export_UKBB_pheno_nodupes2 # this makes local functions executable when bsubbed



# 3) exclude ambiguous alleles:  A/T or C/G
function remove_ambiguous_alleles { 
sumstatsLoc=$1
#  1     2     3    4     5
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
awk ' 
FNR <= NR {
if (FNR == 1) {print $0 }
else { 
if( toupper($4) == "A" && toupper($5) == "T" || toupper($5) == "A" && toupper($4) == "T" || toupper($4) == "G" && toupper($5) == "C" || toupper($5) == "G" && toupper($4) == "C") {}
else {print $0}
}
}
' $sumstatsLoc > $sumstatsLoc$'_noAmbiguousAlleles'

head $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc$'_noAmbiguousAlleles'
wc -l $sumstatsLoc
}
export -f remove_ambiguous_alleles # this makes local functions executable when bsubbed




function Add_PRS_TO_COVS { 
prsLoc=$1
covsLoc=$2
prsName=$3


awk -v prsName=$prsName 'FNR == NR { file1[ "\""$1"\""] = $3; next; }
FNR <= NR {  
if(FNR == 1) {print $0"\t\""prsName"\""}
else if( $6 in file1 ) {print $0"\t"file1[$6]} } 
' $prsLoc $covsLoc > $covsLoc$'_'$prsName


head $covsLoc$'_'$prsName
wc -l $covsLoc$'_'$prsName
}
export -f Add_PRS_TO_COVS # this makes local functions executable when bsubbed





# Binary traits: converts the results of a PLINK2 to the common sumstats format
function ConvertPLINK2BinaryPhenoToSumtats { 
numIndis=$1 # this is the standard plink .phe.glm.logistic.hybrid file
plink2file=$2
outlocation=$3

# AAA PLINK2
#   1 		2     3        4       5       6         7       8     9        10    11          12           13      14       15    16     17
# #CHROM  POS     ID      REF     ALT     A1      A1_FREQ FIRTH?  TEST    OBS_CT  OR      LOG(OR)_SE      L95     U95     Z_STAT  P    ERRCODE
# TO:
# chr   pos   SNP   A1   A2   Freq1.Hapmap   b   se   p   N
#  1     2     3    4     5        6         7    8   9   10

awk  -v numIndis=$numIndis ' FNR <= NR {
if (FNR == 1) {print "chr\tpos\tSNP\tA1\tA2\tFreq1.Hapmap\tb\tse\tp\tN" }
else { 
# find the non-effect allele, this could either be REF or ALT, as PLINK2 is shit and inconsistent
other_allele=$4 # default to REF
if(other_allele == $6) {other_allele = $5} # if REF is the same as A1, then use ALT
print $1"\t"$2"\t"$3"\t"$6"\t"other_allele"\t"$7"\t"log($11)"\t"$12"\t"$16"\t"numIndis
}
} ' OFS='\t' $plink2file > $outlocation
}
export -f ConvertPLINK2BinaryPhenoToSumtats # this makes local functions executable when bsubbed


# generates manhattan plot from a PLINK2 GWAS association file
function produce_manhattan_PLINK2 { 
plinkAssocFile=$1 # this is the standard plink .phe.glm.logistic.hybrid file
output_loc=$2
plottitle=$3
sigTreshold=$4
linear=$5


# map the PLINK2 output format depending on if results file is binary or linear
if [[ "$linear" == '1' ]]; then
echo 'linear'
#   1 		2     3        4       5       6           7       8       9        10    11       12       13      14       15    16  
#CHROM	   POS	  ID	  REF	   ALT	   A1   	A1_FREQ	 TEST	 OBS_CT	  BETA	  SE	  L95	    U95	 T_STAT	  P	    ERRCODE
# TO
# extract summary stats into format into the following signature: # SNP CHR BP         P    zscore
awk '{ 
if ( FNR == 1) {print "SNP\tCHR\tBP\tP\tzscore"}
else  { print $3"\t"$1"\t"$2"\t"$15"\t"$14 }
 }' $plinkAssocFile > $output_loc$'_gwasResults'
else
echo 'binary'
# 1      2         3         4   5  6      7         8        9      10       11            12        13          14          15         16             17
#CHROM	POS	      ID	    REF	ALT	A1	A1_FREQ	   FIRTH?	TEST	OBS_CT	OR	        LOG(OR)_SE	L95       	U95	        Z_STAT	     P	        ERRCODE 
#9	  22115286	rs944797	C	T	T	0.493031	N	     ADD	190724	0.875327	0.0254157	0.832792	0.920035	-5.2392 	1.61274e-07	. 
# TO
# extract summary stats into format into the following signature: # SNP CHR BP         P    zscore
awk '{ 
if ( FNR == 1) {print "SNP\tCHR\tBP\tP\tzscore"}
else  { print $3"\t"$1"\t"$2"\t"$16"\t"$15 }
 }' $plinkAssocFile > $output_loc$'_gwasResults'
fi

# call Rscript to produce plot
arguments='/home/mk907/scripts/R/make_manhattan.R '$output_loc$'_gwasResults '$output_loc$' '$plottitle$' '$sigTreshold
Rscript $arguments # has to be run on farm5, as dgx-server has the "X11 is not available" error

#rm -rf $output_loc$'_gwasResults'
}
export -f produce_manhattan_PLINK2 # this makes local functions executable when bsubbed




# converts my sumstats format into the PRS-CS format, postfixing input with '_PRSCS' (specify if binary)
function convertToPRSCS { 
sumstatsLoc=$1
isBinary=$2

# convert my _sumstats format:
# chr     pos     SNP     A1      A2      Freq1.Hapmap    b       se      p       N
# to PRS-CS:
# SNP          A1   A2   BETA/OR      P

awk -v isBinary="$isBinary" ' 
FNR <= NR {
if (FNR == 1) {
coefName="BETA"
if(isBinary == "1") {coefName="OR"}
print "SNP\tA1\tA2\t"coefName"\tP" 
}
else { 
COEF=$7
if(isBinary == "1") {COEF=exp($7)}
print $3"\t"$4"\t"$5"\t"COEF"\t"$9 
}
}
' OFS='\t' $sumstatsLoc > $sumstatsLoc$'_PRSCS'

#head $sumstatsLoc$'_PRSCS'

}
export -f convertToPRSCS # this makes local functions executable when bsubbed


# Performs PRS-CS for 1 population (EUR)
function performPRSCS { 
eursums=$1
Neur=$2
validBim=$3
pheno=$4
chrom=$5
ncores=$6


export MKL_NUM_THREADS=$ncores
export NUMEXPR_NUM_THREADS=$ncores
export OMP_NUM_THREADS=$ncores

echo "Made directory at "$scratchLocV4$pheno$"/"
prscsxdir="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/prscsscript/PRScsx/"
prscsrefs="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/data/prscsrefs/"
scratchLocV4="/rds/project/jmmh2/rds-jmmh2-projects/aneurysm/v4/scratch/"
mkdir -p $scratchLocV4$pheno$"/"
/usr/local/software/master/miniconda/2/bin/python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$validBim --sst_file=$eursums --n_gwas=$Neur --pop=EUR --out_dir=$scratchLocV4$pheno/ --out_name=$pheno --chrom=$chrom  --seed=42

}
export -f performPRSCS # this makes local functions executable when bsubbed




# performs PRS-CS for each of the 22 chroms
function PRSCS_allChroms { 
rawPRSLoc=$1
pheno_name=$2
scratchLoca=$3

convertToPRSCS $rawPRSLoc '1'

# get AVERAGE sample size, this was better than minimum
Neur=$(awk '{ if(FNR > 1) {sum += $10; n++ } } END { printf("%0.f", sum / n) ; }' $rawPRSLoc)
echo $Neur


validBim=$rawPRSLoc$'_validBim'
# create a 'fake' .bim file for PRSx from my sumstats file that has all the info
awk '{if (FNR !=1) {print $1"\t"$3"\t0\t"$2"\t"$4"\t"$5} }' $rawPRSLoc > $validBim$'.bim'


for ((i=1; i<=22; i++)); do # $numChroms
mkdir -p $scratchLoca$pheno_name$'/'
# 1) PRS-CS
# check if output doesn't already exist
outfile=$scratchLoca$pheno_name$'/'$pheno_name$'_EUR_pst_eff_a1_b0.5_phiauto_chr'$i$'.txt'
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
b=1
echo "PRS-CS SUBMITTING: "$pheno_name$' '$i

sbatch --mem 43000 --cpus-per-task 16 --time 4:0:0 --partition cclake --account danesh-sl3-cpu -e $scratchLoca$'logs/'$pheno_name$'_'$i$'.err' -o $scratchLoca$'logs/'$pheno_name$'_'$i$'.out' --wrap "performPRSCS $rawPRSLoc$'_PRSCS' $Neur $validBim $pheno_name $i 16"
#performPRSCS $rawPRSLoc$'_PRSCS' $Neur $validBim $pheno_name $i 16

fi

done # end of chrom loop

} 
export -f PRSCS_allChroms # this makes local functions executable when bsubbed









# performs LDpred2 for each of the 22 chromsa
function LDpred2_allChroms { 
rawPRSLoc=$1
pheno_name=$2
scratchLoca=$3
shaPRSscriptLoc="/home/mk907/scripts/R/shaPRS.R"


for ((i=1; i<=22; i++)); do # $numChroms
mkdir -p $scratchLoca$pheno_name$'/'

# 2) LDpred2
# LDPred2 Auto does not really use the test set, the 'test set' I used was for CD (it is only used to eliminate extreme PRS, without looking at the phenotypes)

LdpredOutDir=$scratchLoca$pheno_name$'_ldpred2/'
mkdir -p $LdpredOutDir
# III) Build PRS via LDpred2
outfile=$LdpredOutDir$pheno_name$'_'$i
#rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
echo LDPRED ${pheno_name}_chrom${i}

arguments='/home/mk907/scripts/R/LDpred2_auto_chrom.R '$baseLDpredLoc$' '$rawLocV4$pheno_name$' '$i$' -1 '$LdpredOutDir$' '$pheno_name$'_'$i$' '$UKBB_PLINK1$'_EUR10K '$shaPRSscriptLoc$' 1'

#Rscript $arguments
sbatch --mem 20000 --cpus-per-task 4 --time 12:0:0 --partition cclake --account danesh-sl3-cpu -e $scratchLoca$'logs/'$pheno_name$'_ldpred2'$i$'.err' -o $scratchLoca$'logs/'$pheno_name$'_ldpred2'$i$'.out' --wrap "Rscript $arguments"

fi


done # end of chrom loop

} 
export -f LDpred2_allChroms # this makes local functions executable when bsubbed






















# evaluates the r^2 and AUC of a PRS on a test set
function Evaluate_PRS { 
allScoreLoc=$1
outResultLoc=$2


#rm -rf $outResultLoc$'_AUC'
#arguments='/home/mk907/scripts/R/AUC_genr.R '$allScoreLoc$' '$outResultLoc$'_AUC'
#Rscript $arguments
#rm -rf $outResultLoc$'_r2'
arguments='/home/mk907/scripts/R/correlator.R '$allScoreLoc$' '$outResultLoc$'_r2 1 1 0 1'
Rscript $arguments
}
export -f Evaluate_PRS # this makes local functions executable when bsubbed



# Performs the shaPRS step that produces a _sumstats formatted file for both the blended and the combined sumstats: output $outPutLocation'_shaprs'
function shaPRS_fixed { 
outPutLocation=$1
pheASumstats=$2
pheBSumstats=$3
rho=$4

# Create input file for adjusting script:
#       SNP	CHR	BP	Beta_A	SE_A	Beta_B	SE_B
# rs4040617   1  779322 -0.0017630 0.008608 -0.010990 0.008592
awk 'FNR == NR { file1[ $3 ] = $7"\t"$8"\t"$4"\t"$5; next; } 
FNR <= NR { { 
if (FNR == 1) {print "SNP\tCHR\tBP\tBeta_A\tSE_A\tA1.x\tA2.x\tBeta_B\tSE_B\tA1.y\tA2.y" } else {
if ( $3 in file1) { print $3"\t"$1"\t"$2"\t"file1[$3]"\t"$7"\t"$8"\t"$4"\t"$5} }   }
}
' $pheBSumstats $pheASumstats > $outPutLocation$'_SE_meta'


# Export out the lFDR values
arguments='/home/mk907/scripts/R/shaPRS_adjust_wrapper.R '$outPutLocation$'_SE_meta '$outPutLocation$'_lFDR_meta /home/mk907/scripts/R/shaPRS.R'
Rscript $arguments



# R script to blend between the two association stats
arguments='/home/mk907/scripts/R/sumstatsBlender_shaPRS_meta_wrapper.R '$pheASumstats$' '$pheBSumstats$' '$outPutLocation$'_lFDR_meta_SNP_lFDR '$rho$' '$outPutLocation$'_shaprsMerged /home/mk907/scripts/R/shaPRS.R'
Rscript $arguments

# in case there were SNPs that were not present in trait 2 these would get removed by the above, we use trait 1's values
awk 'FNR == NR { file1[ $3 ] = $0; next;  }
FNR <= NR { if ($3 in file1) {print file1[$3] } else {print $0} }
' $outPutLocation$'_shaprsMerged' $pheASumstats > $outPutLocation$'_shaprs'


# find out the max Q-val of the run
awk 'BEGIN {max = -1; maxSNP="";} FNR <= NR { { if ( FNR > 1 && max < $3) { max = $3; maxSNP=$1 } } } END {print "SNP with max Qval: "maxSNP"\t"max;}' FS=" " $outPutLocation$'_lFDR_meta_SNP_lFDR'

}
export -f shaPRS_fixed # this makes local functions executable when bsubbed



# calculates N_eff based on a .phe file and a list of individuals used for assoc
function GetN_eff {
pheFile=$1 # stores the cases/controls counts
totalsFile=$2 # stores the actualy people used for training (this may be more than the phefile)

# 3) filter for QC failed in UKBB
awk 'FNR == NR { file1[$2] = $2; next; }
FNR <= NR { 

if($2 in file1) {

if($3 == 1) controls++
else if ($3 == 2) cases++
 }
 } END { 
# print "numCases "cases" / numcontrols: "controls
 # calculate N_eff
N_eff = 4 / (1 / cases + 1 / controls) 
print int(N_eff) 
} ' $totalsFile $pheFile
}


PRS_rawLoc=$1
phe_name=$2
testPheno_name=$3
PRSCS=$4
keeplist=$AAD_ABDOMINAL_TEST_GWAS_V4
datLo=$dataLocV4
# builds a PRS from PRS-CS or LDPred2 formatted score file
function BuildPRS_cluster_PRSCS_V5 {
PRS_rawLoc=$1
phe_name=$2
testPheno_name=$3
PRSCS=$4
keeplist=$5
datLo=$6


i=21

# generate PRS per each chrom
mkdir -p $datLo$'PRSProfiles/'
mkdir -p $datLo$'PRSProfiles/'$phe_name$'/'
for ((i=1; i<=22; i++)); do

# only submit if it doesn't exist yet
outfile=$datLo$'PRSProfiles/'$phe_name$'/'$i$'.sscore'
rm -rf $outfile
if [ -s "$outfile" ] ; then
#echo "ALREADY EXISTS: "$outfile
b=1
else
b=1
pgen='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr'$i

#PRS-CS
if [ $PRSCS == "1" ] ; then
echo "PRS-CS PRS"
arguments=' --pfile '$pgen$' --memory  43000 --threads 16 --score '$PRS_rawLoc$' 2 4 6 ignore-dup-ids cols=-scoreavgs,+scoresums --keep '$keeplist$' --pheno '$rawLoc$testPheno_name$'_all.phe --out '$datLo$'PRSProfiles/'$phe_name$'/'$i

else
# LDpred2/Rapido
echo "LDpred2 PRS"
arguments=' --pfile '$pgen$' --memory  43000 --threads 16 --score '$PRS_rawLoc$' ignore-dup-ids cols=-scoreavgs,+scoresums --keep '$keeplist$' --pheno '$rawLoc$testPheno_name$'_all.phe --out '$datLo$'PRSProfiles/'$phe_name$'/'$i
fi 

sbatch --mem 43000 --cpus-per-task 16 --time 4:0:0 --partition cclake --account danesh-sl3-cpu --wrap "/home/mk907/software/plink2/plink2 $arguments"


#$plink2 $arguments

fi # if score exists

done

}
export -f BuildPRS_cluster_PRSCS_V5 # this makes local functions executable when bsubbed


# Builds the individual profile scores for all chromosomes for PRS-CS
function BuildPRS_allChroms_PRSCS { 
penoNam=$1
scratchLoca=$2
keeplist=$3
datLo=$4

#  concat all chrom PRS into single files
cat $scratchLoca$penoNam$'/'$penoNam$'_EUR_pst_eff_a1_b0.5_phiauto_chr'*'.txt' > $scratchLoca$penoNam$'/'$penoNam$'_PRSCS'
head $scratchLoca$penoNam$'/'$penoNam$'_PRSCS'
wc -l $scratchLoca$penoNam$'/'$penoNam$'_PRSCS' 

# remove dupes:
awk ' { a[$2]++; if( a[$2] <= 1){ print $0} }' $scratchLoca$penoNam$'/'$penoNam$'_PRSCS' > $scratchLoca$penoNam$'/'$penoNam$'_PRSCS_no_dupes'

BuildPRS_cluster_PRSCS_V5 $scratchLoca$penoNam$'/'$penoNam$'_PRSCS_no_dupes' $penoNam$'_PRSCS' 'AAD_abdominal' '1' $keeplist $datLo

}
export -f BuildPRS_allChroms_PRSCS # this makes local functions executable when bsubbed


# Builds the individual profile scores for all chromosomes for LDpred2
function BuildPRS_allChroms_LDpred2 { 
penoNam=$1
scratchLoca=$2
keeplist=$3
datLo=$4

# Concat and build PRS
LdpredOutDir=$scratchLoca$penoNam$'_ldpred2/'
#  concat all chrom PRS into single files
cat $LdpredOutDir$penoNam$'_'[!a-z] $LdpredOutDir$penoNam$'_'[!a-z][!a-z] > $LdpredOutDir$penoNam$'_ldpred2'

head $LdpredOutDir$penoNam$'_ldpred2'
wc -l $LdpredOutDir$penoNam$'_ldpred2'

# remove dupes:
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $LdpredOutDir$penoNam$'_ldpred2' > $LdpredOutDir$penoNam$'_ldpred2_no_dupes'

BuildPRS_cluster_PRSCS_V5 $LdpredOutDir$penoNam$'_ldpred2_no_dupes' $penoNam$'_LDpred2' 'AAD_abdominal' '0' $keeplist $datLo
}
export -f BuildPRS_allChroms_LDpred2 # this makes local functions executable when bsubbed





function CreateBinaryGWASData_V4 { 
casesRawFile=$1
phen_nam=$2
rawlocation=$3
testSetLocation=$4

# exclude the QC failed indis
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 == 0) {print $0} }
' $MASTER_EXCLUDE_CASECONTROL $casesRawFile > $scratchLoc$phen_nam$'_temp'

cat $CONTROLS_ALL_GWAS $scratchLoc$phen_nam$'_temp' > $rawlocation$phen_nam$'_GWAS_temp'
rm -rf $scratchLoc$phen_nam$'_temp'

# possible that the conttols incldued some of the cases, so we need to filter duplicates
awk ' { a[$1]++; if( a[$1] <= 1){ print $0} }' $rawlocation$phen_nam$'_GWAS_temp' > $rawlocation$phen_nam$'_GWAS'
wc -l $rawlocation$phen_nam$'_GWAS_temp'
wc -l $rawlocation$phen_nam$'_GWAS' # 

rm -rf $rawlocation$phen_nam$'_GWAS_temp'




# All Training = All - TEST_SET
#a) GWAS
awk 'FNR == NR { 
file1[ $1 ] = $1;
next; 
}
FNR <= NR { 
if ($1 in file1 == 0) {print $1"\t"$1} }
' $testSetLocation $rawlocation$phen_nam$'_GWAS' > $rawlocation$phen_nam$'_GWAS_TRAINING'
head $rawlocation$phen_nam$'_GWAS_TRAINING'
wc -l $rawlocation$phen_nam$'_GWAS_TRAINING' # 190960


# create phenotype files ( need the +1 , as 1 = control, 2 = case
awk 'FNR == NR { 
phenos[$1] = $1
next; }
FNR <= NR {  
if(FNR == 1) {print "FID\tIID\tphe"} else {
if( $1 in phenos) { print $1"\t"$1"\t2" }
else { print $1"\t"$1"\t1" }} }
' $casesRawFile $rawlocation$phen_nam$'_GWAS' > $rawlocation$phen_nam$'_all.phe'
head $rawlocation$phen_nam$'_all.phe'

}
export -f CreateBinaryGWASData_V4 # this makes local functions executable when bsubbed



# sums each 22 chroms PRS into a single profile score, is not tied to V4
function SumPRS_Agnostic {
phe_name=$1
datLoc=$2


i=21

# sum across each chrom
for ((i=1; i<=22; i++)); do

# if it is first iteration, we just copy the score file for 1st chrom, as 
if [ $i == 1 ] ; then
echo "first chrom"
# 3rd col is the true pheno, 6th col is the scoresum, we only want a IID,pheno,PRS file

awk '{print $2"\t"$3"\t"$6 }' $datLoc$'PRSProfiles/'$phe_name$'/'$i$'.sscore' > $datLoc$'PRSProfiles/'$phe_name$'_all.sscore'
cp $datLoc$'PRSProfiles/'$phe_name$'_all.sscore' $datLoc$'PRSProfiles/'$phe_name$'_all.temp'

else

awk 'FNR == NR { sums[ $1 ] = $3; next; }
FNR <= NR { if( $2 in sums ) {print $2"\t"$3"\t"$6+sums[$2] } }
' $datLoc$'PRSProfiles/'$phe_name$'_all.temp' $datLoc$'PRSProfiles/'$phe_name$'/'$i$'.sscore' > $datLoc$'PRSProfiles/'$phe_name$'_all.sscore'
cp $datLoc$'PRSProfiles/'$phe_name$'_all.sscore' $datLoc$'PRSProfiles/'$phe_name$'_all.temp'

fi 

done

rm -rf $datLoc$'PRSProfiles/'$phe_name$'_all.temp'


}
export -f SumPRS_Agnostic # this makes local functions executable when bsubbed


GWAS3_IBD_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[2]}PLINK.qassoc
GWAS3_CD_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[3]}PLINK.qassoc
GWAS3_UC_assoc=${ibdResultsLoc}GWAS/${ibdarray_target_extended[4]}PLINK.qassoc

$GWAS3_IBD_assoc
$GWAS3_CD_assoc
$GWAS3_UC_assoc 


outPutLocation=$scratchLocV4$'rhotest'
pheASumstats=$rawLocV4$"AAAGen_Leicester_UKBB_shaprs"
pheBSumstats=$rawLocV4$'AAD_RELATED_hm3'
pheASumstats=$rawLocV4$'AAD_RELATED_hm3'
# calculates empriical overlap between studies: load the rho into session via: source $outPutLocation$'rho'
function GetEmpiricalRho {
outPutLocation=$1
pheASumstats=$2
pheBSumstats=$3

# perform meta analysis
arguments='/home/mk907/scripts/R/metaAnalysis.R '$pheASumstats$' '$pheBSumstats$' '$outPutLocation$'_meta /home/mk907/scripts/R/shaPRS.R'
Rscript $arguments


# in GWAS3 IBD naive GWAS, find SNPs with p > 0.01 # this is to find the empirical correlation, we find the 'mildly to non associated' SNPS
# otherwise we would get 
awk '{ if ( $9 > 0.01 && NR >1) {print $3} }' $outPutLocation$'_meta' > $outPutLocation$'_unassociated'
wc -l $outPutLocation$'_unassociated'
# 6,933,546

# Get same SNPs in the two subphenos
# create a file with signature RSiD, Beta_CD, SE_CD, Beta_UC, SE_UC
awk 'FNR == NR { file1[ $1 ] = $1; next; } 
FNR <= NR { { if ( $3 in file1) { print $3"\t"$7"\t"$8 } }  }
'  $outPutLocation$'_unassociated' $pheASumstats > $outPutLocation$'A_unassociated'

awk 'FNR == NR { file1[ $1 ] = $0; next; } 
FNR <= NR { { 
if (FNR == 1) {print "RSid\tBeta_CD\tSE_CD\tBeta_UC\tSE_UC" }
if ( $3 in file1) { print file1[$3]"\t"$7"\t"$8 } }  
}
'  $outPutLocation$'A_unassociated' $pheBSumstats > $outPutLocation$'A_B_unassociated'


# run R script to get rho (correlation between studies,  this is a scalar)
# quantifies the linear association between estimated coefficients, for NON-disease SNPs
# rho = cor(Beta_CD/SE_CD, Beta_UC/SE_UC)
arguments='/home/mk907/scripts/R/rho.R '$outPutLocation$'A_B_unassociated '$outPutLocation$'rho'
Rscript $arguments

# 
#source $outPutLocation$'rho'
}
export -f GetEmpiricalRho # this makes local functions executable when bsubbed


#########################################
# install PRS-CS(x)

######################
mkdir -p $prscsrefs
mkdir -p $prscsscript




# download ref panels for EUR and EAS (need to add 'dl=1' to the file page to make wget work: https://superuser.com/questions/470664/how-to-download-dropbox-files-using-wget-command
cd $prscsrefs
wget https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=1
wget https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=1
wget https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=1


# decompress them
mv snpinfo_mult_ukbb_hm3?dl=1 snpinfo_mult_ukbb_hm3
mv ldblk_ukbb_eur.tar.gz?dl=1 ldblk_ukbb_eur.tar.gz
mv ldblk_ukbb_eas.tar.gz?dl=1 ldblk_ukbb_eas.tar.gz

tar -xvf ldblk_ukbb_eas.tar.gz
tar -xvf ldblk_ukbb_eur.tar.gz

# get the prsCSx repository:
cd $prscsscript

git clone https://github.com/getian107/PRScsx.git



# run test to confirm that installation works
pip install h5py
python $prscsxdir/PRScsx.py --help --user
python $prscsxdir/PRScsx.py --ref_dir=$prscsrefs --bim_prefix=$prscsxdir/test_data/test --sst_file=$prscsxdir/test_data/EUR_sumstats.txt,$prscsxdir/test_data/EAS_sumstats.txt --n_gwas=200000,100000 --pop=EUR,EAS --chrom=22 --phi=1e-2 --out_dir=$scratchLocV4 --out_name=test



#########################################
# install LDSC
cd $ldscDir

git clone https://github.com/bulik/ldsc.git
cd ldsc
# need to activate the python2 miniconda to install it
# https://docs.hpc.cam.ac.uk/hpc/software-tools/python.html#using-anaconda-python
module load miniconda/2
conda env create --file environment.yml

git clone https://github.com/bulik/ldsc.git
cd ldsc
# need to activate the python2 miniconda to install it
# https://docs.hpc.cam.ac.uk/hpc/software-tools/python.html#using-anaconda-python
module load miniconda/2
conda env create --file environment.yml
#  --a1 A1               Name of A1 column , this is described as the 
# "A1: Allele 1, interpreted as ref allele for signed sumstat.", the 'ref' here means that LDSC thinks its the effect allele 

###################################

module load miniconda/2
source activate ldsc

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2

tar -jxvf eur_w_ld_chr.tar.bz2
bunzip2 w_hm3.snplist.bz2



##################################
#########################################

# install LDpred

# Download the LDpred ref data:
cd $baseLDpredLoc

# https://figshare.com/articles/dataset/European_LD_reference/13034123
wget https://figshare.com/ndownloader/articles/13034123/versions/3
unzip 3


# install LDpred2
R
install.packages("R.utils")
install.packages("bigsparser")
remotes::install_github("privefl/bigsnpr")



# Need to have a dummy validation set for LDpred2 in PLINK1 format..
# convert UKBB each Chrom to PLINK1 for hapmap3 vars
# generate PRS per each chrom
mkdir -p $UKBB_PLINK1

for ((i=1; i<=22; i++)); do

# only submit if it doesn't exist yet
outfile=$UKBB_PLINK1$i$'.bim'
if [ -s "$outfile" ] ; then
echo "ALREADY EXISTS: "$outfile
else
pgen='/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr'$i
arguments=' --pfile '$pgen$' --memory  43000 --maf 0.001 --hard-call-threshold 0.1 --geno 0.02  --extract '$hapmap3_SNPs$' --make-bed --allow-extra-chr --allow-no-sex --snps-only --bp-space 1 --out '$UKBB_PLINK1$i
#sbatch --mem 43000 --cpus-per-task 16 --time 4:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$plink2 $arguments"
$plink2 $arguments

fi


done



# merge chroms into one
plinkFileList=$scratchLocV4'/mergelist.txt'
rm -f $plinkFileList

for ((j=1; j<=22; j++)); do
echo $UKBB_PLINK1$j >> ${plinkFileList}

done

arguments=' --memory  43000 --merge-list '$plinkFileList$' --make-bed --out '$UKBB_PLINK1$'ALL --allow-extra-chr --allow-no-sex'
#$plink $arguments

sbatch --mem 43000 --cpus-per-task 16 --time 12:0:0 --partition cclake --account danesh-sl3-cpu --wrap "$plink $arguments"


# subset it into a random 10K people from the EUR subset
numItems=$(wc -l < "$CONTROLS_ALL_GWAS")
mySeed=42
arguments='/home/mk907/scripts/R/randomNums.R '$numItems$' '$mySeed$' '$scratchLocV4$'randomNums'
Rscript $arguments

head -n 10000 $scratchLocV4$'randomNums' > $scratchLocV4$'randomNums_10k'

awk 'FNR == NR { test[ $0 ] = $0; next; } FNR <= NR { if( FNR in test ) {print $1"\t"$1 } }
' $scratchLocV4$'randomNums_10k' $CONTROLS_ALL_GWAS > $scratchLocV4$'randomNums_10k_indis' 
head $scratchLocV4$'randomNums_10k_indis' 
wc -l $scratchLocV4$'randomNums_10k_indis' 


arguments=' --bfile '$UKBB_PLINK1$'ALL --keep '$scratchLocV4$'randomNums_10k_indis --make-bed --out '$UKBB_PLINK1$'_EUR10K --allow-extra-chr --allow-no-sex'
$plink $arguments

########################

# I was getting this error for the Nelson CAD, complaining that the median Zscore is too large
# but lookoing at  the manhattan plots: https://www.nature.com/articles/ng.3913/figures/2
# this very large GWAS seems plausible that we have so many significant SNPs
# since there is no way to disable this error: https://github.com/JonJala/mtag/issues/86
# I had to manually edit the "/home/mk907/software/mtag/mtag/mtag_munge.py"
# commenting this out
#       if not args.a1_inc:
#            logging.info(
#                check_median(dat.SIGNED_SUMSTAT, signed_sumstat_null, args.median_z_cutoff, sign_cname))


# Traceback (most recent call last):
  # File "/home/mk907/software/mtag/mtag/mtag.py", line 1577, in <module>
    # mtag(args)
  # File "/home/mk907/software/mtag/mtag/mtag.py", line 1343, in mtag
    # DATA_U, DATA, args = load_and_merge_data(args)
  # File "/home/mk907/software/mtag/mtag/mtag.py", line 273, in load_and_merge_data
    # GWAS_d[p], sumstats_format[p] = _perform_munge(args, GWAS_d[p], gwas_dat_gen, p)
  # File "/home/mk907/software/mtag/mtag/mtag.py", line 166, in _perform_munge
    # munged_results = munge_sumstats.munge_sumstats(argnames, write_out=False, new_log=False)
  # File "/home/mk907/software/mtag/mtag/mtag_munge.py", line 881, in munge_sumstats
    # check_median(dat.SIGNED_SUMSTAT, signed_sumstat_null, args.median_z_cutoff, sign_cname))
  # File "/home/mk907/software/mtag/mtag/mtag_munge.py", line 525, in check_median
    # raise ValueError(msg.format(F=name, M=expected_median, V=round(m, 2)))
# ValueError: WARNING: median value of SIGNED_SUMSTAT is 0.69 (should be close to 0.0). This column may be mislabeled.









