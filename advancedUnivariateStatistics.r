#############
#Testing for covariation
#############

CorDF <- cbind(UniD$CorrW, UniD$DistP) #add continuous variables to new dataframe
CorDF <- cbind(CorDF, with(UniD, model.matrix(~ Graze + 0)), with(UniD, model.matrix(~ Dung + 0))) #add the categorical variables as dummy variables

colnames(CorDF)[1] <- "CorrW"
colnames(CorDF)[2] <- "DistP"

print(CorDF)
cor(CorDF)


#############
#Testing for spatial autocorrelation
#############

library(ape)
library(ade4)

station.dist <- dist(cbind(UniD$Long, UniD$Lat))
sppr.dist <- dist (UniD$SppR)

mantel.rtest(station.dist, sppr.dist, nrepet = 9999)


#############
#Linear modelling
#############

ModLM <-lm(SppR ~ Dung + Graze + DistP + CorrW, data=UniD)
summary(ModLM)


#############
#Linear mixed modelling
#############

library(lme4)
modLMM<-lmer(SppR~Dung+Graze+DistP+CorrW+(1|Type), data=UniD)
summary(modLMM) # does not give us p-values, so will have to use the anova option

modDung<-lmer(SppR~Graze+DistP+CorrW+(1|Type), data=UniD)
modGraze<-lmer(SppR~Dung+DistP+CorrW+(1|Type), data=UniD)
modDistP<-lmer(SppR~Dung+Graze+CorrW+(1|Type), data=UniD)
modCorrW<-lmer(SppR~Dung+Graze+DistP+(1|Type), data=UniD)

anova (modLMM, modDung)
anova (modLMM, modGraze)
anova (modLMM, modDistP)
anova (modLMM, modCorrW)

#Rescale variables as per warning (only continuous variables)
UniD$DistPS<-scale (UniD$DistP)
UniD$CorrWS<-scale (UniD$CorrW)

#rerun with scaled variables:
modLMM<-lmer(SppR~Dung+Graze+DistPS+CorrWS+(1|Type), data=UniD)
modDung<-lmer(SppR~Graze+DistPS+CorrWS+(1|Type), data=UniD)
modGraze<-lmer(SppR~Dung+DistPS+CorrWS+(1|Type), data=UniD)
modDistP<-lmer(SppR~Dung+Graze+CorrWS+(1|Type), data=UniD)
modCorrW<-lmer(SppR~Dung+Graze+DistPS+(1|Type), data=UniD)

anova (modLMM, modDung)
anova (modLMM, modGraze)
anova (modLMM, modDistP)
anova (modLMM, modCorrW)


#############
#Generalized linear modelling
#############

library(lme4)
modGLM<-glm(SppR~Dung+Graze+DistP+CorrW, family = poisson, data=UniD)
summary (modGLM)


#############
#Generalized linear mixed modelling
#############

library(lme4)

modGLMM<-glmer(SppR~Dung+Graze+DistP+CorrW+(1|Type), family = poisson, data=UniD)
modDung<-glmer(SppR~Graze+DistP+CorrW+(1|Type), family = poisson, data=UniD)
modGraze<-glmer(SppR~Dung+DistP+CorrW+(1|Type), family = poisson, data=UniD)
modDistP<-glmer(SppR~Dung+Graze+CorrW+(1|Type), family = poisson, data=UniD)
modCorrW<-glmer(SppR~Dung+Graze+DistP+(1|Type), family = poisson, data=UniD)

anova (modGLMM, modDung)
anova (modGLMM, modGraze)
anova (modGLMM, modDistP)
anova (modGLMM, modCorrW)

#Rescale issues as above so will rerun the model with the rescaled variables
modGLMM<-glmer(SppR~Dung+Graze+DistPS+CorrWS+(1|Type), family = poisson, data=UniD)
modDung<-glmer(SppR~Graze+DistPS+CorrWS+(1|Type), family = poisson, data=UniD)
modGraze<-glmer(SppR~Dung+DistPS+CorrWS+(1|Type), family = poisson, data=UniD)
modDistP<-glmer(SppR~Dung+Graze+CorrWS+(1|Type), family = poisson, data=UniD)
modCorrW<-glmer(SppR~Dung+Graze+DistPS+(1|Type), family = poisson, data=UniD)

anova (modGLMM, modDung)
anova (modGLMM, modGraze)
anova (modGLMM, modDistP)
anova (modGLMM, modCorrW)


#############
#Posthoc testing for categorical variables if significant through modelling
#############

library (multcomp)
summary(glht(modGLMM, linfct=mcp(Dung="Tukey")))
summary(glht(modGLMM, linfct=mcp(Graze="Tukey")))

windows()
par(mfrow = c(2,1))
attach(UniD)
boxplot(SppR~Dung, cex.lab=1.5, cex.axis=1.5, xlab= 'Dung type', ylab= 'Number of species')
boxplot(SppR~Graze, cex.lab=1.5, cex.axis=1.5, xlab= 'Grazing type', ylab= 'Number of species')


#############
#Variance inflation
#############

library(car)
library(dplyr)
library(caret)
library(lme4)

training.samples <- BFLdummy$lepiRareNew %>%
  createDataPartition(p = 0.8, list = FALSE)

  train.data  <- BFLdummy[training.samples, ]
  test.data <- BFLdummy[-training.samples, ]

  # Build the model
    all.parmsBFLdummy <- glm(lepiRareNew ~
      #featureSize_area +
      #distToNearest +
      #numCorrConn +
      #avgSlope +
      elevation +
      #corrWidth +
      #windSpeed +
      #airTemp +
      vegHeight +
      percGrasses +
      #percWoodyHerb +
      #percAlien +
      #percRocky +
      #percBareGround +
      #avgGrassSpecies +
      #avgHerbSpecies +
      #avgFlowerHeadCount +
      foodPlants +
      featureTypeDong +
      featureTypeOutcr +
      featureTypeWetl +
      #simplAspect +
      #corrOrientNS +
      connToNearestT +
      connToPAT,
      family = poisson,
      data = train.data)


  # Make predictions
  predictions <-   all.parmsBFLdummy %>% predict(test.data)
  # Model performance
  data.frame(
    RMSE = RMSE(predictions, test.data$lepiRareNew),
    R2 = R2(predictions, test.data$lepiRareNew)
  )

vif(all.parmsBFLdummy)


#############
#Model averaging (Only variables with VIF < 3)
#############


options(na.action = "na.fail")

  SpROdonata <- lm(overallSpR ~
    elev +
    pH +
    waterTemp +
    conductivity +
    vegHeight +
    grassy +
    woody +
    reedy +
    bareGround +
    ppCatSPSC +
    ppCatSPLC +
    #ppCatLPSC +
    ppCatLPLC,
    #connTrue,
    data = odoDummy)


  resultsminOdonata <- dredge(SpROdonata, trace = 2)

  # grab best supported models
  subset(resultsminOdonata, delta <4)
  minOdonataBest <- subset(resultsminOdonata, delta <4)

  MA.ests<-model.avg(resultsminOdonata, subset= delta < 4, revised.var = TRUE)

  confint(MA.ests, level = 0.95)
  coefTable(MA.ests)
  sw(MA.ests)

  minOdonataResults <-cbind(coefTable(MA.ests), confint(MA.ests, level = 0.95))
  write.table(minOdonataResults, file="SpROdonataResults.txt", row.names=T, sep="\t")

  weights <- sw(MA.ests)

  write.table(weights, file="WeightsSpROdonataResults.txt", row.names=T, sep="\t")



  ######################
  ### Relationships between environmental variables
  ######################

  pairs.panels(BFL[, 6:28], method = "spearman", density = FALSE, ellipses = FALSE, hist.col = "grey", lm = TRUE, stars = TRUE)


  ###OR

  library(GGally)
  library(multcomp)

  BFL<- read.delim("bflSiteData.txt")

  ggpairs(BFL,
          columns = 6:24,
          aes(color = featureType,
              alpha = 0.5),
            upper = list(continuous = wrap("cor", size = 1.5)))



  BFL$featureType <- factor(BFL$featureType)

  modALL <- lm(featureSize_area ~ featureType, data = BFL) #
  summary(glht(modALL, linfct=mcp(featureType = "Tukey")))


#for categorical variables
  library(chisq.posthoc.test)

  Chi<- read.delim("Chi.txt")
  Chi


  chisq.test(Chi$PA, Chi$FT)

  tablePA <- table(Chi$PA, Chi$FT)
  tablePA

  chisq.posthoc.test(tablePA,
                     method = "bonferroni")


  ####################END####################
