
# Code repository for "Integrated genetic and conventional risk factor model for predicting abdominal aortic aneurysm"

DOI: XXX

This respository represents the last snapshot of the bash and R scripts used to generate our results and is provided as-is. As the analysis involved a lot of input/output operations on very large files, these were generated asynchronously on a cluster. Thus these scripts are meant to be executed on the command line manually, block-by-block, waiting for the remote jobs to finish and verifying the integrity of the resulting files at each step. To reduce code duplication, certaint functions that were reused multiple times are defined only once across all files, however, they may be called from different scripts.

Code for the AAA PRS development in the paper can be found under /prs/: 

1. AAA_functions.sh: variables and reusable functions
2. AAA_pub.sh: data pre-processing and UKB analyses
3. AAA_adjunct.sh: shaPRS processing and generation of final PRS

Auxilliary R scripts used by the above scripts can be found under prs/R/.

Code for the Time-To-Event analyses in the paper can be found under /time_to_event/: 

