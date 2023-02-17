# gets the rho, empirical correlation between studies from unassociated SNPs

# grab command line arguments
args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 2) {stop("not enough Arguments received")} 
 
inputDataLoc = args[1]
outputLoc = args[2]



# Debug wars
# inputDataLoc = "C:/softwares/Cluster/GIANT/miniPRS/PRS/CD_UC_unassociated"
# outputLoc = "C:/softwares/Cluster/GIANT/miniPRS/PRS/rho"





# 1. load phenos
inputData= read.table(inputDataLoc, header = T) 


rho = cor(inputData$Beta_CD/inputData$SE_CD, inputData$Beta_UC/inputData$SE_UC)
print(paste("rho is:",round(rho,4) ))
if ( rho < 0 ) rho = 0

outdata=paste("rho=",rho,sep="")
write.table(outdata, outputLoc, sep="\t", row.names = FALSE, col.names = F, quote = FALSE)
print(paste("saved rho to: ", outputLoc ))




# find out how/why the inputs are negatively associated

inputData_orig = inputData

random = sample(1:nrow(inputData_orig))
random_10K = random[1:10000]

inputData = inputData_orig[random_10K,]
Beta_CD_D_SE_CD = inputData$Beta_CD/inputData$SE_CD
Beta_UC_D_SE_UC = inputData$Beta_UC/inputData$SE_UC


filename = paste(outputLoc, "_scaled.pdf" , sep="")
pdf(filename, width=6.4 , height=6.4);
plot(Beta_CD_D_SE_CD,Beta_UC_D_SE_UC,main = paste("Relationship of Beta/SE for  A vs B (rho:",round(rho,3),")")  )
dev.off()


filename = paste(outputLoc, "_Betas.pdf" , sep="")
pdf(filename, width=6.4 , height=6.4);
plot(inputData$Beta_CD,inputData$Beta_UC,main =  paste("Relationship of Betas for  A vs B (rho:",round(rho,3),")"))
dev.off()






#################################


# approx correlation:

# nkl0 = 9492 # number of controls overlapping
# 
# nk1 = 3810 # number of cases in k
# nk0 = 9492 # number of controls in k
# nk=nk1+nk0 # total num indis in k
# 
# nl1= 3765 # number of cases in l
# nl0= 9492 # number of controls in l
# nl= nl1+nl0 # number of indis in l

# formula:
# approx_cor= (nkl0 * sqrt(nk1*nl1 / (nk0*nl0) )   ) / sqrt(nk*nl)
# approx_cor


# (-0.5258)^2

