library(stringr)
library(dplyr)
library(biomod2)
setwd("H:/HOTR_models")
pth <- "response/"
BRTpth <- paste0(pth, "response.GBMtrial1.models.out")
RFpth <- paste0(pth, "response.RFtrial1.models.out")

### BRT ###
load(BRTpth)
BRTmod <- response.GBMtrial1.models.out
BRTmod
BRTeval <- get_evaluations(BRTmod)
BRTeval[, "Testing.data", "GBM",,]
BRTeval[,c("Cutoff", "Sensitivity", "Specificity"), "GBM",,]

eval <- get_evaluations(BRTmod,as.data.frame=TRUE)
eval$Model <- "BRT"

load(paste0(pth, "models/GBMtrial1/response_AllData_Full_GBM"))
get_formal_model(response_AllData_Full_GBM)
response_AllData_Full_GBM@model_options
summary(get_formal_model(response_AllData_Full_GBM))

### RF ###
load(RFpth)
RFmod <- response.RFtrial1.models.out
RFmod
RFeval <- get_evaluations(RFmod)
RFeval[, "Testing.data", "RF",,]
RFeval[,c("Cutoff", "Sensitivity", "Specificity"), "RF",,]

eRF <- get_evaluations(RFmod,as.data.frame=TRUE)
eRF$Model <- "RF"

load(paste0(pth, "models/RFtrial1/response_AllData_Full_RF"))
get_formal_model(response_AllData_Full_RF)
response_AllData_Full_RF@model_options
summary(get_formal_model(response_AllData_Full_RF))

### Graphs ###
eval <- bind_rows(eval, eRF)
eval <- eval %>%
  mutate(
    strat = case_when(
      str_detect(Model.name, "Full") ~ "Full",
      str_detect(Model.name, "RUN1") ~ "RUN1",
      str_detect(Model.name, "RUN2") ~ "RUN2",
      str_detect(Model.name, "RUN3") ~ "RUN3"
    )
  )

boxplot(eval$Testing.data[which(eval$Eval.metric=="ROC")] ~ eval$Model[which(eval$Eval.metric=="ROC")],
        ylab="ROC AUC", xlab="Model", ylim=c(0.9, 1), xpd=FALSE)

RFvarimp <- response_AllData_Full_RF@model$importance
barplot(RFvarimp ~ row.names(RFvarimp), horiz=TRUE)

gp1 <- models_scores_graph(BRTmod, by="models", metrics = c("ROC", "TSS", "KAPPA"))
