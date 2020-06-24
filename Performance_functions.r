##############################################################################
# Functions for performance metrics at two different thresholds, either
# sensitivity = specificity threshold or sensitivity = 0.95 threshold.
#
# Michelle M. Fink, michelle.fink@colostate.edu
# Colorado Natural Heritage Program, Colorado State University
# Code Last Modified 06/22/2020
#
# Code licensed under the GNU General Public License version 3.
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/
#############################################################################

library(ROCR)
library(PresenceAbsence)
library(ggplot2)

opt.cut <- function(perf, pred, type="sens95"){
  # Returns [Sensitivity, Specificity, Threshold]
  # adapted from:
  #https://www.r-bloggers.com/a-small-introduction-to-the-rocr-package/
  #type="senspec" calculates sensitivity = specificity threshold
  #type="sens95" caclulates sensitivity = 0.95 threshold

  mapply(FUN=function(x, y, p){
    if(type == "senspec"){
      d = (x - 0)^2 + (y-1)^2
      ind = which(d == min(d))
      c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
        cutoff = p[[ind]])
    } else if(type == "sens95"){
      d = which(y >= 0.95)
      ind = which(y == min(y[d]))
      if(length(ind > 1)){ # if ties, take first
        idx <- ind[1]
      } else {
        idx <- ind
      }
      c(sensitivity = y[[idx]], specificity = 1-x[[idx]],
        cutoff = p[[idx]])
    }
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

testEval <- function(evobj, run_name, cutype="sens95"){
  #evobj == data.table(prob, resp)
  #cutype -> see 'type' in functon opt.cut above
  pred.mod <- prediction(evobj$prob, evobj$resp)
  perf.roc <- performance(pred.mod, "tpr", "fpr")
  perf.err <- performance(pred.mod, "err")
  AUC <- performance(pred.mod, measure = "auc")@y.values[[1]]
  sst <- opt.cut(perf.roc, pred.mod, type = cutype)
  TSS <- sst[1] + sst[2] - 1
  err_rate <- perf.err@y.values[[1]][which(perf.err@x.values[[1]] == sst[3])]

  #Presence.Absence insists on its own data structure
  padat <- data.frame(ID=seq(1:length(pred.mod@predictions[[1]])),
                      OBS=as.numeric(levels(pred.mod@labels[[1]]))[pred.mod@labels[[1]]],
                      PRED=pred.mod@predictions[[1]])

  pacmx <- cmx(padat, threshold = sst[3])
  kappa <- Kappa(pacmx, st.dev = FALSE)
  PCC <- pcc(pacmx, st.dev = FALSE)

  return(list(run_name, AUC, TSS, err_rate, kappa, PCC, sst[1], sst[2], sst[3]))
}

get.pred <- function(modobj, modtype="GLM"){
  if(modtype=="RF"){
    p_rf.mod <- predict(modobj, type = "prob")
    pred.mod <- prediction(p_rf.mod[,2], modobj$y)
  } else if(modtype=="BRT") {
    p_brt.mod <- predict(modobj, type = "response")
    pred.mod <- prediction(p_brt.mod, modobj$data$y)
  } else {
    pred.mod <- prediction(modobj@resp$mu, modobj@resp$y)
  }
  return(pred.mod)
}