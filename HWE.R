# setting working directory :

#fhwe <- read.csv("Filtered HWE estimates.csv", na.strings = '') # assigning NAs to missing values - IMPORTANT
#head(fhwe, 3)
#str(fhwe)
library(HardyWeinberg) #loading HWE package
fhwe=t(combined_poly[,c("PATIENT","RS8192678_CC1_CT2_TT3")])

# use column one (loci ID) as rownames
rownames(fhwe) <- fhwe[, 1] # use names in col 1 to set row names
fhwe <-fhwe[, -1] # remove column one
fhwe[2,] = as.numeric(fhwe[2,])
# Count homo- and heterozygotes ...
AA <- rowSums(x = fhwe == '1', na.rm = T)
AB <- rowSums(x = fhwe == '2', na.rm = T)
BB <- rowSums(x = fhwe == '3', na.rm = T)


# Create condensed version of the fhwe data frame
dim(fhwe) # query dimensions of the original data frame
#fhwe1 <- fhwe[, 14:16] # select only the last three cols holding the values of interest
fhwe1 = c(as.data.frame(AA)[2,], as.data.frame(AB)[2,], as.data.frame(BB)[2,])
names(fhwe1)=c("AA","AB","BB")
HWExact(fhwe1)
head(fhwe1)
str(fhwe1)
fhwe2 <- data.matrix(fhwe1) # it is required to convert the data frame into a data matrix to make the HWExact function work

HWExact(fhwe2)$pval

