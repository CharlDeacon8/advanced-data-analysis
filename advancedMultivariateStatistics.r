#############
#Checking similarity/dissimilarity of environmental variables between sites
#############

library(vegan)
library(tidyverse)
library(pander)
library(mvabund)


BFLcomp <- read.delim("bflSiteDummyData.txt") #reads the matrix table

rownames(BFLcomp) <- BFLcomp[,1] #Assigning row names from 1st column
BFLcomp[,1] <- NULL #Removing the first column
BFLcomp #Check if the data imported correctly

ENV_sel <- select(BFLcomp, featureSize_area:connToPAT)    #include only selected environmental variables (chem_d)
ENV_sel                 #Check if selection was done correctly
as_tibble(ENV_sel)                  #store dataframe as a as_tibble

ENV_m <- vegdist(ENV_sel, method="euclidean")
ENV_clust <- hclust(ENV_m, method="average")


ENV_selScaled <- scale(ENV_sel)   #scale all the variables, becuase they are all measured on different scales (chem_s)

ENV_pca <- capscale(ENV_selScaled ~ 1, method = "euc")      #Runs a PCA/distLM with euclidian distance for all environmental variables

biplot(ENV_pca, col=c("blue", "red"))

scores(ENV_pca, display = "sites", choices = c(1:10)) %>%
  head %>% round(1) %>% pander::pander()              #See scores for study sites


  scores(ENV_pca, display = "species", choices = c(1:10)) %>%
    head %>% round(1) %>% pander::pander()              #See scores for variables


    summary(ENV_pca)$cont %>%
      as.data.frame %>%
        round(2) %>%
        pander::pander()            #get eigenvalues, all >1 are important to variation. The fewer axes it takes for eigenvalues to be < 1, the better the ordination fits the distance matrix.
                                    #MDS1-10 are important for variation, explaining 81% of the total proportion, cumulatively

screeplot(ENV_pca, type="lines")


par(mgp=c(4, 1, 0),
    mar=c(6, 6, 1, 1),
     mfrow=c(1,2),
    las=1, cex.lab=1, cex.axis=1)
plot(ENV_clust, las=1,
      labels=ENV_sel$name,
      xlab="Sample",
     ylab="Euclidean distance")
plot(ENV_pca, display = "sites", las=1)
ordicluster(ENV_pca, ENV_clust)


#############
#Environmental variables as drivers of compositional change
#############


bflMatrix <- read.delim("bflOverallMatrix.txt") #reads the matrix table

rownames(bflMatrix) <- bflMatrix[,1] #Assigning row names from 1st column
bflMatrix[,1] <- NULL #Removing the first column
bflMatrix #Check if the data imported correctly


BFLcomp <- read.delim("bflSiteDummyData.txt") #reads the matrix table

rownames(BFLcomp) <- BFLcomp[,1] #Assigning row names from 1st column
BFLcomp[,1] <- NULL #Removing the first column
BFLcomp #Check if the data imported correctly

ENV_sel <- select(BFLcomp, featureSize_area:connToPAT)    #include only selected environmental variables (chem_d)
ENV_sel


rankindex(ENV_sel, bflMatrix, indices = c("euc", "man", "gow", "bra", "kul"),stepacross= FALSE, method = "spearman")    #Check which distribution has the best fit


mod0 <- capscale(bflMatrix ~ 1, ENV_sel, dist="bray")       #Model with intercept only
#mod1 <- capscale(bflMatrix ~ ., ENV_sel, dist="bray")       #Model with all variables
mod1 <- capscale(bflMatrix ~
  featureSize_area +
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
  percRocky +
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
  ENV_sel,
  dist="bray")       #Model with all selected variables (corresponding to univariate - VIF checked)


step.res <- ordiR2step(mod0, scope = formula(mod1), perm.max = 9999, direction = "forward")

step.res$anova  # Summary table for selection

#Final model with only important variables

dbRDA = capscale(bflMatrix ~
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
  featureTypeWetl +
  #simplAspectE,
  #simplAspectN +
  #simplAspectNE +
  #simplAspectNW +
  simplAspectSE,
  #simplAspectSW +
  #corrOrientNS +
  #connToNearestT +
  #connToPAT,
  ENV_sel,
  dist="bray")

plot(dbRDA) # use base plot, might be done with ggplot2
anova(dbRDA) # is the model significant?

anova(dbRDA, by="axis", perm.max=9999) ## test axes for significance
anova(dbRDA, by="terms", permu=9999) ## test for sig. environ. Variables

dbRDAresults <- anova(dbRDA, by="terms", permu=9999)
write.table(dbRDAresults, file="dbRDAresults.txt", row.names=T, sep="\t")


summary(dbRDA)$cont %>%
  as.data.frame %>%
    round(2) %>%
    pander::pander()    #get eigenvalues, all >1 are important to variation. The fewer axes it takes for eigenvalues to be < 1, the better the ordination fits the distance matrix.
                        #Cap1 = 1.52, explaining 20% of variation


    screeplot(dbRDA, type="lines")
      abline(a=1, b=0, lty=2)       #visualise results


#############
#Attractive plots of composition
#############

bflMatrix <- read.delim("bflOverallMatrix.txt") #reads the matrix table

rownames(bflMatrix) <- bflMatrix[,1] #Assigning row names from 1st column
bflMatrix[,1] <- NULL #Removing the first column
bflMatrix #Check if the data imported correctly


MESpp<-mvabund(bflMatrix)

BFLcomp <- read.delim("bflSiteDummyData.txt") #reads the matrix table

rownames(BFLcomp) <- BFLcomp[,1] #Assigning row names from 1st column
BFLcomp[,1] <- NULL #Removing the first column
BFLcomp #Check if the data imported correctly


vare.cap <- capscale(bflMatrix ~ featureType, BFLcomp,
                     dist="bray")
vare.cap
#plot(vare.cap)
#anova(vare.cap)

### PLOTS WITH ELLIPSES SHOWING THE STANDARD ERRORS (RUN THIS ENTIRE CHUNK AT THE SAME TIME)

plot(vare.cap, type="points", display = "sites")
with(BFLcomp, text(vare.cap, display="sites", labels = as.character(featureType),
                    col=as.numeric(featureType)))

pl <- with(BFLcomp, ordiellipse(vare.cap, featureType, kind="se", conf=0.95, lwd=2,
                                 draw = "polygon", col=1:4, border=1:4,
                                 alpha=50))

###NOT RUN, just a variation of above
### PLOTS WITH POLYGONS SHOWING ALL ITEMS IN A CLASS/TREATMENT (RUN THIS ENTIRE CHUNK AT THE SAME TIME)

plot(vare.cap, type="points", display = "sites")
with(BFLcomp, text(vare.cap, display="sites", labels = as.character(featureType),
                      col=as.numeric(featureType)))

pl <- with(BFLcomp, ordihull(vare.cap, featureType, lwd=2,
                                  draw = "polygon", col=1:4, border=1:4,
                                  alpha=50))


#############
#Variance partitioning
#############







#############
#Simper analysis
#############

library(vegan)

simpBio <- read.delim("simpBio.txt") #reads the matrix table

rownames(simpBio) <- simpBio[,1] #Assigning row names from 1st column
simpBio[,1] <- NULL #Removing the first column
simpBio #Check if the data imported correctly

simpEnv <- read.delim("simpEnv.txt")

(sim <- with(simpEnv, simper(simpBio, biotope)))

summary(sim)
