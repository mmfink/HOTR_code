# Parameter reduction functions
require(gam)
require(dplyr)

comb <- function(x, ...) {
  #https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

devianceExp <- function(dat, y, vars){
  wgt <- rep(1, times = length(y))
  outlists <- foreach(i=vars, .combine="comb", .multicombine=TRUE,
              .packages=c("gam", "stats", "dplyr"),
              .init=list(list(), list())) %dopar% {
    x <- dat %>% pull(var = i)
    g <- try(gam(y ~ x, family = "binomial", weights = wgt), silent = TRUE)
    dev.broke <- try((1-g$dev/g$null.deviance)<0, silent = TRUE)
    if(class(dev.broke) == "try-error") {dev.broke = TRUE}
    if("try-error"%in%class(g) | dev.broke){
      #gam.failed=TRUE
      g <- glm(y~x+x^2, weights=wgt, family = "binomial")
      y.fit <- predict(g, type = "response")
    } else {
      y.fit <- predict.Gam(g, type = "response")
      #gam.failed=FALSE
    }
    list(i, 100*(1-g$dev/g$null.deviance))
  }
  names(outlists) <- c("invar", "dev_exp")
  vartbl <- as.data.frame(mapply(rbind, outlists))
  return(vartbl)
}

numCorr <- function(dat, max_cor, dev_exp){
  # dat has to be just the input vars, nothing else
  cmat <- cor(dat, use = "pairwise.complete.obs")
  smat <- cor(dat, method = "spearman", use = "pairwise.complete.obs")

  if (dim(dat)[1] < 2000) {
    kmat <- cor(dat, method = "kendall", use = "pairwise.complete.obs")
  } else {
    s <- sample(seq(1:dim(dat)[1]), size = 2000, replace = FALSE)
    kmat <- cor(dat[s, ], method = "kendall", use = "pairwise.complete.obs")
  }

  cmat = pmax(abs(cmat), abs(smat), abs(kmat), na.rm = TRUE)
  cmat[is.na(cmat)] <- 0
  High_cor <- apply(abs(cmat) > max_cor, 2, sum) - 1
  corIssues <- tibble(invar = attr(High_cor, "names"), num_cor = High_cor)
  corIssues$dev_exp <- dev_exp
  bmat <- apply(abs(cmat) > max_cor, 2, identity)
  diag(bmat) <- FALSE
  corIssues$keep <- 0
  vartbl <- corIssues %>% arrange(desc(dev_exp))
  vartbl$keep[1] <- 1
  if(vartbl$num_cor[1] != (nrow(vartbl) - 1)){
    for(i in 1:nrow(vartbl)){
      if((vartbl$num_cor[i] > 0) & (vartbl$keep[i] != 99)){
        chk.var <- as.character(vartbl$invar[i])
        out <- which((bmat)[, chk.var])
        vartbl$keep <- if_else(vartbl$invar %in% names(out), 99, vartbl$keep)
      }
      if((vartbl$keep[i] == 0) & (vartbl$dev_exp > 0.00000001)){vartbl$keep[i] <- 1}
    }
  }
  return(vartbl)
}
