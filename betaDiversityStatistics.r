##################
###Beta diversity
##################

library(betapart)
library(vegan)
library(tidyverse)
library(ape)

bflMatrix <- read.delim("bflOverallMatrix.txt") #reads the matrix table

rownames(bflMatrix) <- bflMatrix[,1] #Assigning row names from 1st column
bflMatrix[,1] <- NULL #Removing the first column
bflMatrix #Check if the data imported correctly



BFLcomp <- read.delim("bflSiteDummyData.txt") #reads the site info table

rownames(BFLcomp) <- BFLcomp[,1] #Assigning row names from 1st column
BFLcomp[,1] <- NULL #Removing the first column
BFLcomp #Check if the data imported correctly

ENV_sel <- select(BFLcomp, featureSize_area:connToPAT)    #include only selected environmental variables (chem_d)
ENV_sel                 #Check if selection was done correctly






multi <- beta.multi.abund(bflMatrix, index.family="bray") #calculate overall/turnover/nestedness components
multi


dist<-beta.pair.abund(bflMatrix, index.family="bray")    #calculate pairwise total beta diversity and components

dist		### to see all parts (this is in triangular distance matrix format)
dist$beta.bray.bal		###to see the turnover component
dist$beta.bray.gra		###to see the nestedness component
dist$beta.bray		###to see the total beta diversity



#### 3.1 PCoA for total beta diversity
TOT <- pcoa(dist$beta.bray, correction = "lingoes")		###lingoes correction to account for negative eigenvalues
TOT

biplot(TOT, plot.axes = c(1,2))


#### 3.2 PCoA for turnover component of beta diversity
TUR <- pcoa(dist$beta.bray.bal, correction = "lingoes")			###lingoes correction to account for negative eigenvalues
TUR

biplot(TUR, plot.axes = c(1,2))


#### 3.3 PCoA for nestedness component of beta diversity
NES <- pcoa(dist$beta.bray.gra, correction = "lingoes")		###lingoes correction to account for negative eigenvalues
NES

biplot(NES, plot.axes = c(1,2))


### 4.1.1 Forward selection for environmental drivers of total beta diversity

tot.int <- rda(TOT$vectors.cor ~ 1, ENV_sel) 			###intercept-only model

tot.glob <- rda(TOT$vectors.cor ~ #featureSize_area +
#distToNearest +
numCorrConn +
#avgSlope +
elevation +
#corrWidth +
#windSpeed +
#airTemp +
vegHeight +
percGrasses +
percWoodyHerb +
#percAlien +
#percRocky +
percBareGround +
avgGrassSpecies +
avgHerbSpecies +
avgFlowerHeadCount +
foodPlants +
featureTypeDong +
featureTypeOutcr +
featureTypeWetl +
simplAspectE +
simplAspectN +
simplAspectNE +
simplAspectNW +
simplAspectSE +
simplAspectSW +
corrOrientNS +
connToNearestT +
connToPAT,
ENV_sel) 	 		###global model with all variables ( . after ~ means ‘include all variables’)



fw.sel <- ordiR2step(tot.int, scope = formula(tot.glob), direction = 'forward', permutations = 9999)


###foodPlants + elevation + featureTypeWetl + simplAspectSE + simplAspectSW + avgHerbSpecies are important variables



### 4.1.2 Forward selection for environmental drivers of turnover component

tur.int <- rda(TUR$vectors.cor ~ 1, ENV_sel) 			###intercept-only model

tur.glob <- rda(TUR$vectors.cor ~ #featureSize_area +
#distToNearest +
numCorrConn +
#avgSlope +
elevation +
#corrWidth +
#windSpeed +
#airTemp +
vegHeight +
percGrasses +
percWoodyHerb +
#percAlien +
#percRocky +
percBareGround +
avgGrassSpecies +
avgHerbSpecies +
avgFlowerHeadCount +
foodPlants +
featureTypeDong +
featureTypeOutcr +
featureTypeWetl +
simplAspectE +
simplAspectN +
simplAspectNE +
simplAspectNW +
simplAspectSE +
simplAspectSW +
corrOrientNS +
connToNearestT +
connToPAT,
ENV_sel) 	 		###global model with all variables ( . after ~ means ‘include all variables’)



fw.selTUR <- ordiR2step(tur.int, scope = formula(tur.glob), direction = 'forward', permutations = 9999)


## foodPlants + elevation + featureTypeWetl are important variables for the turnover component


### 4.1.3 Forward selection for environmental drivers of nestedness component

nes.int <- rda(NES$vectors.cor ~ 1, ENV_sel) 			###intercept-only model

nes.glob <- rda(NES$vectors.cor ~ #featureSize_area +
#distToNearest +
numCorrConn +
#avgSlope +
elevation +
#corrWidth +
#windSpeed +
#airTemp +
vegHeight +
percGrasses +
percWoodyHerb +
#percAlien +
#percRocky +
percBareGround +
avgGrassSpecies +
avgHerbSpecies +
avgFlowerHeadCount +
foodPlants +
featureTypeDong +
featureTypeOutcr +
featureTypeWetl +
simplAspectE +
simplAspectN +
simplAspectNE +
simplAspectNW +
simplAspectSE +
simplAspectSW +
corrOrientNS +
connToNearestT +
connToPAT,
ENV_sel) 	 		###global model with all variables ( . after ~ means ‘include all variables’)



fw.selNES <- ordiR2step(nes.int, scope = formula(nes.glob), direction = 'forward', permutations = 9999)


### no important variables were selected important variable



##  Testing all the selected environmental variables for significance

## TOTAL

###foodPlants + elevation + featureTypeWetl + simplAspectSE + simplAspectSW + avgHerbSpecies

TotRDA <- rda(TOT$vectors.cor ~
#featureSize_area +
#distToNearest +
#numCorrConn +
#avgSlope +
elevation +
#corrWidth +
#windSpeed +
#airTemp +
#vegHeight +
#percGrasses +
#percWoodyHerb +
#percAlien +
#percRocky +
#percBareGround +
avgGrassSpecies +
#avgHerbSpecies +
#avgFlowerHeadCount +
foodPlants +
#featureTypeDong +
#featureTypeOutcr +
featureTypeWetl +
#simplAspectE +
#simplAspectN +
#simplAspectNE +
#simplAspectNW +
simplAspectSE +
simplAspectSW,
#corrOrientNS +
#connToNearestT +
#connToPAT,
ENV_sel)

plot(TotRDA) # use base plot, might be done with ggplot2
anova(TotRDA) # is the model significant?

anova(TotRDA, by="axis", perm.max=999) ## test axes for significance
anova(TotRDA, by="terms", permu=999) ## test for sig. environ. Variables

overallTotal <- anova(TotRDA, by="terms", permu=999) ## test for sig. environ. Variables
write.table(overallTotal, file="overallTotal.txt", row.names=T, sep="\t")


summary(TotRDA)$cont %>%
as.data.frame %>%
round(2) %>%
pander::pander()    #get eigenvalues, all >1 are important to variation. The fewer axes it takes for eigenvalues to be < 1, the better the ordination fits the distance matrix.



screeplot(TotRDA, type="lines")
abline(a=1, b=0, lty=2)       #visualise results

#variance partitioning

###foodPlants + elevation + featureTypeWetl + simplAspectSW

### constrained = variance explained

### unconstrained = unexplained variance

### conditioned = Variance explained by left out variables (covariates)

partElevation <- rda(TOT$vectors.cor ~ elevation + Condition(foodPlants + featureTypeWetl + simplAspectSW), data = ENV_sel)
summary(partElevation)

R2partElevation <- RsquareAdj(rda(TOT$vectors.cor ~ elevation + Condition(foodPlants + featureTypeWetl + simplAspectSW), data = ENV_sel))
R2partElevation

anova.cca(partElevation, step = 1000)


partFoodPlants <- rda(TOT$vectors.cor ~ foodPlants + Condition(elevation + featureTypeWetl + simplAspectSW), data = ENV_sel)
summary(partFoodPlants)

R2partFoodPlants <- RsquareAdj(rda(TOT$vectors.cor ~ foodPlants + Condition(elevation + featureTypeWetl + simplAspectSW), data = ENV_sel))
R2partFoodPlants

anova.cca(partFoodPlants, step = 1000)


partFtype <- rda(TOT$vectors.cor ~ featureTypeWetl + Condition(foodPlants + elevation + simplAspectSW), data = ENV_sel)
summary(partFtype)

R2partFtype <- RsquareAdj(rda(TOT$vectors.cor ~ featureTypeWetl + Condition(foodPlants + elevation + simplAspectSW), data = ENV_sel))
R2partFtype

anova.cca(partFtype, step = 1000)


partAspect <- rda(TOT$vectors.cor ~ simplAspectSW + Condition(foodPlants + featureTypeWetl + elevation), data = ENV_sel)
summary(partAspect)

R2partAspect <- RsquareAdj(rda(TOT$vectors.cor ~ simplAspectSW + Condition(foodPlants + featureTypeWetl + elevation), data = ENV_sel))
R2partAspect

anova.cca(partAspect, step = 1000)



partElevationFP <- rda(TOT$vectors.cor ~ elevation + foodPlants + Condition(featureTypeWetl + simplAspectSW), data = ENV_sel)
summary(partElevationFP)

R2partElevationFP <- RsquareAdj(rda(TOT$vectors.cor ~ elevation + foodPlants + Condition(featureTypeWetl + simplAspectSW), data = ENV_sel))
R2partElevationFP

anova.cca(partElevationFP, step = 1000)


partElevationFT <- rda(TOT$vectors.cor ~ elevation + featureTypeWetl + Condition(foodPlants + simplAspectSW), data = ENV_sel)
summary(partElevationFT)

R2partElevationFT <- RsquareAdj(rda(TOT$vectors.cor ~ elevation + featureTypeWetl + Condition(foodPlants + simplAspectSW), data = ENV_sel))
R2partElevationFT

anova.cca(partElevationFT, step = 1000)




partElevationAspect <- rda(TOT$vectors.cor ~ elevation + simplAspectSW + Condition(foodPlants + featureTypeWetl), data = ENV_sel)
summary(partElevationAspect)

R2partElevationAspect <- RsquareAdj(rda(TOT$vectors.cor ~ elevation + simplAspectSW + Condition(foodPlants + featureTypeWetl), data = ENV_sel))
R2partElevationAspect

anova.cca(partElevationAspect, step = 1000)





partFoodPlantsFT <- rda(TOT$vectors.cor ~ foodPlants + featureTypeWetl + Condition(elevation + simplAspectSW), data = ENV_sel)
summary(partFoodPlantsFT)

R2partFoodPlantsFT <- RsquareAdj(rda(TOT$vectors.cor ~ foodPlants + featureTypeWetl + Condition(elevation + simplAspectSW), data = ENV_sel))
R2partFoodPlantsFT

anova.cca(partFoodPlantsFT, step = 1000)




partFoodPlantsAspect <- rda(TOT$vectors.cor ~ foodPlants + simplAspectSW + Condition(elevation + featureTypeWetl), data = ENV_sel)
summary(partFoodPlantsAspect)

R2partFoodPlantsAspect <- RsquareAdj(rda(TOT$vectors.cor ~ foodPlants + simplAspectSW + Condition(elevation + featureTypeWetl), data = ENV_sel))
R2partFoodPlantsAspect

anova.cca(partFoodPlantsAspect, step = 1000)




partFtypeAspect <- rda(TOT$vectors.cor ~ featureTypeWetl + simplAspectSW + Condition(foodPlants + elevation), data = ENV_sel)
summary(partFtypeAspect)

R2partFtypeAspect <- RsquareAdj(rda(TOT$vectors.cor ~ featureTypeWetl + simplAspectSW + Condition(foodPlants + elevation), data = ENV_sel))
R2partFtypeAspect

anova.cca(partFtypeAspect, step = 1000)




partElevationFTAspect <- rda(TOT$vectors.cor ~ elevation + featureTypeWetl + simplAspectSW + Condition(foodPlants), data = ENV_sel)
summary(partElevationFTAspect)

R2partElevationFTAspect <- RsquareAdj(rda(TOT$vectors.cor ~ elevation + featureTypeWetl + simplAspectSW + Condition(foodPlants), data = ENV_sel))
R2partElevationFTAspect

anova.cca(partElevationFTAspect, step = 1000)




partElevationFPAspect <- rda(TOT$vectors.cor ~ elevation + foodPlants + simplAspectSW + Condition(featureTypeWetl), data = ENV_sel)
summary(partElevationFPAspect)

R2partElevationFPAspect <- RsquareAdj(rda(TOT$vectors.cor ~ elevation + foodPlants + simplAspectSW + Condition(featureTypeWetl), data = ENV_sel))
R2partElevationFPAspect

anova.cca(partElevationFPAspect, step = 1000)




partElevationFPFT <- rda(TOT$vectors.cor ~ elevation + foodPlants + featureTypeWetl + Condition(simplAspectSW), data = ENV_sel)
summary(partElevationFPFT)

R2partElevationFPFT <- RsquareAdj(rda(TOT$vectors.cor ~ elevation + foodPlants + featureTypeWetl + Condition(simplAspectSW), data = ENV_sel))
R2partElevationFPFT

anova.cca(partElevationFPFT, step = 1000)





partFPFTAspect <- rda(TOT$vectors.cor ~ foodPlants + featureTypeWetl + simplAspectSW + Condition(elevation), data = ENV_sel)
summary(partFPFTAspect)

R2partFPFTAspect <- RsquareAdj(rda(TOT$vectors.cor ~ foodPlants + featureTypeWetl + simplAspectSW + Condition(elevation), data = ENV_sel))
R2partFPFTAspect

anova.cca(partFPFTAspect, step = 1000)





partFPFTAspectElev <- rda(TOT$vectors.cor ~ foodPlants + featureTypeWetl + simplAspectSW + elevation, data = ENV_sel)
summary(partFPFTAspectElev)

R2partFPFTAspectElev <- RsquareAdj(rda(TOT$vectors.cor ~ foodPlants + featureTypeWetl + simplAspectSW + elevation, data = ENV_sel))
R2partFPFTAspectElev

anova.cca(partFPFTAspectElev, step = 1000)



library(eulerr)

overallBeta <-  c(
A = 1.6,
B = 1.2,
C = 1.1,
D = 0.6,
"A&B" = 2.8,
"B&C" = 2.2,
"B&D" = 1.7,
"A&C" = 2.5,
"A&D" = 2.2,
"C&D" = 1.6,
"B&C&D" = 2.7,
"A&B&D" = 3.3,
"A&B&C" = 3.7,
"A&C&D" = 3,
"A&B&C&D" = 30.2)

fit1 <- euler(overallBeta, shape = "ellipse")

plot(fit1,
quantities = TRUE,
labels = list(font = 4))

###


#TURNOVER COMPONENT

## foodPlants + elevation + featureTypeWetl are important variables for the turnover component


TurRDA <- rda(TUR$vectors.cor ~
#featureSize_area +
#distToNearest +
#numCorrConn +
#avgSlope +
elevation +
#corrWidth +
#windSpeed +
#airTemp +
#vegHeight +
#percGrasses +
#percWoodyHerb +
#percAlien +
#percRocky +
#percBareGround +
#avgGrassSpecies +
#avgHerbSpecies +
#avgFlowerHeadCount +
foodPlants +
#featureTypeDong +
#featureTypeOutcr +
featureTypeWetl,
#simplAspectE +
#simplAspectN +
#simplAspectNE +
#simplAspectNW +
#simplAspectSE +
#simplAspectSW +
#corrOrientNS +
#connToNearestT +
#connToPAT,
ENV_sel)

plot(TurRDA) # use base plot, might be done with ggplot2
anova(TurRDA) # is the model significant?

anova(TurRDA, by="axis", perm.max=999) ## test axes for significance
anova(TurRDA, by="terms", permu=999) ## test for sig. environ. Variables

overallTurnover <- anova(TurRDA, by="terms", permu=999) ## test for sig. environ. Variables
write.table(overallTurnover, file="overallTurnover.txt", row.names=T, sep="\t")

summary(TurRDA)$cont %>%
as.data.frame %>%
round(2) %>%
pander::pander()    #get eigenvalues, all >1 are important to variation. The fewer axes it takes for eigenvalues to be < 1, the better the ordination fits the distance matrix.



screeplot(TurRDA, type="lines")
abline(a=1, b=0, lty=2)       #visualise results

#variance partitioning

###foodPlants + elevation + featureTypeWetl


partElevation <- rda(TUR$vectors.cor ~ elevation + Condition(foodPlants + featureTypeWetl), data = ENV_sel)
summary(partElevation)

R2partElevation <- RsquareAdj(rda(TUR$vectors.cor ~ elevation + Condition(foodPlants + featureTypeWetl), data = ENV_sel))
R2partElevation

anova.cca(partElevation, step = 1000)


partFP <- rda(TUR$vectors.cor ~ foodPlants + Condition(elevation + featureTypeWetl), data = ENV_sel)
summary(partFP)

R2partFP <- RsquareAdj(rda(TUR$vectors.cor ~ foodPlants + Condition(elevation + featureTypeWetl), data = ENV_sel))
R2partFP

anova.cca(partFP, step = 1000)


partFT <- rda(TUR$vectors.cor ~ featureTypeWetl + Condition(foodPlants + elevation), data = ENV_sel)
summary(partFT)

R2partFT <- RsquareAdj(rda(TUR$vectors.cor ~ featureTypeWetl + Condition(foodPlants + elevation), data = ENV_sel))
R2partFT

anova.cca(partFT, step = 1000)


partElevationFP <- rda(TUR$vectors.cor ~ elevation + foodPlants + Condition(featureTypeWetl), data = ENV_sel)
summary(partElevationFP)

R2partElevationFP <- RsquareAdj(rda(TUR$vectors.cor ~ elevation + foodPlants + Condition(featureTypeWetl), data = ENV_sel))
R2partElevationFP

anova.cca(partElevationFP, step = 1000)



partElevationFT <- rda(TUR$vectors.cor ~ elevation + featureTypeWetl + Condition(foodPlants), data = ENV_sel)
summary(partElevationFT)

R2partElevationFT <- RsquareAdj(rda(TUR$vectors.cor ~ elevation + featureTypeWetl + Condition(foodPlants), data = ENV_sel))
R2partElevationFT

anova.cca(partElevationFT, step = 1000)



partFPFT <- rda(TUR$vectors.cor ~ foodPlants + featureTypeWetl + Condition(elevation), data = ENV_sel)
summary(partFPFT)

R2partFPFT <- RsquareAdj(rda(TUR$vectors.cor ~ foodPlants + featureTypeWetl + Condition(elevation), data = ENV_sel))
R2partFPFT

anova.cca(partFPFT, step = 1000)



partFPFTElev <- rda(TUR$vectors.cor ~ foodPlants + featureTypeWetl + elevation, data = ENV_sel)
summary(partFPFTElev)

R2partFPFTElev <- RsquareAdj(rda(TUR$vectors.cor ~ foodPlants + featureTypeWetl + elevation, data = ENV_sel))
R2partFPFTElev

anova.cca(partFPFT, step = 1000)





library(eulerr)

overallTurn <-  c(
A = 2,
B = 1.8,
C = 1.6,
"A&B" = 3.8,
"B&C" = 3.3,
"A&C" = 3.5,
"A&B&C" = 16)
fit2 <- euler(overallTurn, shape = "ellipse")

plot(fit2,
quantities = TRUE,
labels = list(font = 4))




#NESTEDNESS COMPONENT

## No important variables (skip this step)


NesRDA <- rda(NES$vectors.cor ~
#featureSize_area +
#distToNearest +
#numCorrConn +
#avgSlope +
#elevation +
#corrWidth +
#windSpeed +
#airTemp +
#vegHeight +
#percGrasses +
#percWoodyHerb +
#percAlien +
#percRocky +
#percBareGround +
#avgGrassSpecies +
#avgHerbSpecies +
#avgFlowerHeadCount +
foodPlants,
#featureTypeDong +
#featureTypeOutcr +
#featureTypeWetl,
#simplAspectE +
#simplAspectN +
#simplAspectNE +
#simplAspectNW +
#simplAspectSE +
#simplAspectSW +
#corrOrientNS +
#connToNearestT +
#connToPAT,
ENV_sel)

plot(NesRDA) # use base plot, might be done with ggplot2
anova(NesRDA) # is the model significant?

anova(NesRDA, by="axis", perm.max=999) ## test axes for significance
anova(NesRDA, by="terms", permu=999) ## test for sig. environ. Variables

overallNestedness <- anova(NesRDA, by="terms", permu=999) ## test for sig. environ. Variables
write.table(overallNestedness, file="overallNestedness.txt", row.names=T, sep="\t")

summary(NesRDA)$cont %>%
as.data.frame %>%
round(2) %>%
pander::pander()    #get eigenvalues, all >1 are important to variation. The fewer axes it takes for eigenvalues to be < 1, the better the ordination fits the distance matrix.



screeplot(NesRDA, type="lines")
abline(a=1, b=0, lty=2)


#variance partitioning

## In this case, none selected
