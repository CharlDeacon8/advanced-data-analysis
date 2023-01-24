#############
#Calculate functional trait diversity metrics
#############

library(FD)
library(vcd)

###Overall Odonata

#Read in dataset

odonataMatrix <- read.delim("odonataMatrix.txt", row.names = 1)   #reads the matrix table
odonataTraits <- read.delim("odonataTraitsRaw.txt", row.names = 1)   #reads the matrix table

#Convert categorical values to dummy/binary values
OdoTraitsDummy <- odonataTraits

#Check if imported correctly
print(OdoTraitsDummy)


##USE THE DUMMY DATA FOR ALL FD CALCULATIONS
##Calculate the FD metrics from the combined species and trait data - initially, leave out the 'corr="cailliez"' term. It will tell you
## if you need to add a correction. There are a few different correction options, see Laliberte et al. for more details if needed.
OverallFD<-dbFD(OdoTraitsDummy, odonataMatrix, corr="cailliez", stand.FRic = TRUE, scale.RaoQ = TRUE)
print(OverallFD)

OverallFD$FRic
OverallFD$FEve
OverallFD$FDiv
OverallFD$FDis
OverallFD$RaoQ

write.table(OverallFD$FRic, file="OverallFDRich.txt", row.names=T, sep="\t")
write.table(OverallFD$FEve, file="OverallFDEve.txt", row.names=T, sep="\t")
write.table(OverallFD$FDiv, file="OverallFDDiv.txt", row.names=T, sep="\t")
write.table(OverallFD$FDis, file="OverallFDDis.txt", row.names=T, sep="\t")
write.table(OverallFD$RaoQ, file="OverallFDRaoQ.txt", row.names=T, sep="\t")


#############
#Test distribution of newly created trait variables (can also test other diversity metric this way)
#############

library(car)
library(MASS)

odoDummy<- read.delim("odoDummy.txt")
dmslDummy<- read.delim("dmslDummy.txt")


gamma <- fitdistr(odoDummy$odoFRic, "gamma")
qqp(odoDummy$odoFRic, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

gamma <- fitdistr(odoDummy$odoRaoQ, "gamma")
qqp(odoDummy$odoRaoQ, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


##Variables then to be used in linear modelling
