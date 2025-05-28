rm(list=ls())
devtools::load_all()
load("C:/Users/go73jec/Documents/R/dynConfiR/Spielplatz/Forster_EEGWEV/first_fit_res_sbj22.RData")
beta <- res_22
model_matrix <- NULL
fixed <- NULL
n  <- 10000
method="simulation_matrix"
pred_matrix <- NULL
maxt0=NULL
restr_tau=NULL
simult_conf = FALSE

predict_dynaViTE_formula <- function(beta, model_matrix=NULL, fixed=NULL, #method="simulation_matrix",
                                     n=NULL, pred_matrix=NULL, maxt0=NULL, restr_tau=NULL,
                                     simult_conf = FALSE) {
  if (is.list(beta) && all(c("beta", "model_matrix", "fixed", "restr_tau", "maxt0") %in% names(beta))) {
    model_matrix <- beta$model_matrix
    fixed <- beta$fixed
    maxt0 <- beta$maxt0
    restr_tau <- beta$restr_tau
    beta <- beta$beta
  }

  # turn 'fixed' from character to a list
  eval(parse(text=paste0("fixed <- list(",fixed,")")))
  #betatrans <- transform_const_pars_to_scale(beta, fixed, maxt0=res_22$maxt0)

  manipulations <- list()
  man_pars <- grep(pattern="_", names(beta), value=TRUE)
  man_pars <- strsplit(man_pars, "_")
  for (i in 1:length(man_pars)) {
    manipulations[[man_pars[[i]][1]]] <- c(manipulations[[man_pars[[i]][1]]],
                                           man_pars[[i]][2])
  }


  # if (method=="simulation_matrix") {     ## So far the only method
    if (is.null(pred_matrix)) {
      if (is.null(n)) stop("For the given method 'simulation_matrix', please specify the n argument")
      pred_matrix <- model_matrix[sample(1:nrow(model_matrix), replace=TRUE, size = n), ]
    }

    #parnames <- c("v", "z", "a", "sz", "t0", "st0", "d", "sv", "tau", "w", "svis", "sigvis", "lambda", "s", "th1", "th2")
    parnames <- c("a", "v", "t0", "d","z",  "sz", "sv", "st0","tau", "th1", "th2", "lambda", "w", "muvis", "sigvis", "svis", "s")
    fixed01 <- c("z", "sz", "w", "d", "st0")
    fixedpos <- c("v", "a", "t0", "sv", "svis", "sigvis", "lambda", "s")
    parammatrix <-  matrix(NA, nrow=nrow(pred_matrix), ncol=length(parnames))
    colnames(parammatrix) <- parnames
    for (i in 1:length(parnames)) {
      if (parnames[i] %in% names(fixed)) {   ### If parameter was fixed, use fixed value for column (no transformation, then!)
        parammatrix[,parnames[i]] <- fixed[[parnames[i]]]
      } else {
        # If parameter was not fixed, then either it was manipulated:
        if (parnames[i] %in% names(manipulations)) {
          parammatrix[,parnames[i]] <-
            pred_matrix[,manipulations[[parnames[i]]], drop=FALSE] %*% # take columns from prediction model matrix and multiply
            beta[grepl(names(beta), pattern=paste0(parnames[i], "_"))] # with respective beta-values
        } else { # or it was fit as a constant, then take parameter as fitted:
          parammatrix[,parnames[i]] <- beta[parnames[i]]
        }
        ## Afterwards, for all parameters involved in the fitting procedure,
        # we have to transform them back to their model scale:
        if (parnames[i] %in% fixed01) parammatrix[,parnames[i]] <- pnorm(parammatrix[,parnames[i]])
        if (parnames[i] %in% fixedpos) parammatrix[,parnames[i]] <- exp(parammatrix[,parnames[i]])
        if (parnames[i] == "sz") parammatrix[,"sz"] <- (pmin(parammatrix[,"z"], (1-parammatrix[,"z"]))*2)*parammatrix[,"sz"]
        if (parnames[i] == "t0") parammatrix[,"t0"] <- maxt0*parammatrix[,"t0"]
        if (parnames[i] == "d") parammatrix[,"d"] <- parammatrix[,"t0"]*parammatrix[,"d"]
        if (parnames[i] == "tau") {
          if (restr_tau == Inf) {
            parammatrix[,parnames[i]] <- exp(parammatrix[,parnames[i]])
          } else if (simult_conf) {
            parammatrix[,parnames[i]] <- pnorm(parammatrix[,parnames[i]])*(maxt0-parammatrix[,"t0"])
          } else {
            parammatrix[,parnames[i]] <- restr_tau * pnorm(parammatrix[,parnames[i]])
          }
        }
      }

    }
    # Fill single missing muvis (index: 14) values with absolute value of decision drift v
    parammatrix[is.na(parammatrix[,14]), 14] <-
      abs(parammatrix[is.na(parammatrix[,14]), 2])

    # Scale the parameters scaled by s by s
    # In the actual function:
    # cbind (a/s, v/s, t0, d, sz, sv/s, st0, z,
    #        tau, th1/s, th2/s, lambda, w, muvis/s, sigvis/s, svis/s, numeric_bounds)
    parammatrix[, c(1, 2, 7, 10, 11, 14, 15, 16)] <-
      parammatrix[, c(1, 2,7, 10, 11, 14, 15, 16)]/parammatrix[,17]

    # bound between trial variability st0 to 2 (seconds)
    parammatrix[,8] <-  parammatrix[,8]*2


    # t0 <- t0+st0/2
    parammatrix[,3] <- parammatrix[,3] + parammatrix[,8]/2

    ## Simulation:
    res <- r_WEV_matrix(parammatrix[,c(1:4, 6:8, 5, 9, 12:16)], delta=0.01, maxT = 15, stop_on_error = TRUE)[,c(1:3)]
    colnames(res) <- c("rt", "response", "conf")

    # Compute confidence rating:
    thetas <- beta[grep(x=names(beta), pattern="theta")]
    thetasUpper <- thetas[grep(x=names(thetas), pattern="Upper")]
    thetasLower <- thetas[grep(x=names(thetas), pattern="Lower")]
    if (length(thetasLower)==0 && length(thetasUpper)==0) {
      thetasLower <- thetas
      thetasUpper <- thetas
    }
    nRatings <- length(thetasUpper)+1
    thetasUpper <- c(-Inf,thetasUpper, Inf)
    thetasLower <- c(-Inf,thetasLower, Inf)


    levels_lower <- cumsum(as.numeric(table(thetasLower)))
    levels_lower <- levels_lower[-length(levels_lower)]
    levels_upper <- cumsum(as.numeric(table(thetasUpper)))
    levels_upper <- levels_upper[-length(levels_upper)]
    thetasLower <- unique(thetasLower)
    thetasUpper <- unique(thetasUpper)

    res <- cbind(res,1)
    res[res[,2]==1, 4] <- as.numeric(as.character(cut(res[res[,2]==1, 3],
                                                                   breaks=thetasUpper, labels = levels_upper)))
    res[res[,2]==-1,4] <- as.numeric(as.character(cut(res[res[,2]==-1, 3],
                                                                    breaks=thetasLower, labels = levels_lower)))

    if (simult_conf) {
      res[,1] <- res[,1] + parammatrix[,'tau']
    }
    colnames(res)[4] <- "rating"
    res <- cbind(pred_matrix, res)
    return(res)
}

transform_const_pars_to_scale <- function(beta, fixed=list(), maxt0=NULL, restr_tau=Inf, simult_conf=FALSE) {
  parnames <- c("a", "v", "t0", "d","z",  "sz", "sv", "st0","tau", "th1", "th2", "lambda", "w", "muvis", "sigvis", "svis", "s")
  fixed01 <- c("z", "sz", "w", "d", "st0")
  fixedpos <- c("v", "a", "t0", "sv", "svis", "sigvis", "lambda", "s")
  for (i in 1:length(parnames)) {
    if (parnames[i] %in% names(fixed)) {
      beta[parnames[i]] <- fixed[[parnames[i]]]
    } else if (parnames[i] %in% names(beta)) {
      if (parnames[i] %in% fixed01) beta[parnames[i]] <- pnorm(beta[parnames[i]])
      if (parnames[i] %in% fixedpos) beta[parnames[i]] <- exp(beta[parnames[i]])
      if (parnames[i] == "sz") beta["sz"] <- (pmin(beta["z"], (1-beta["z"]))*2)*beta["sz"]
      if (parnames[i] == "t0") beta["t0"] <- maxt0*beta["t0"]
      if (parnames[i] == "d") beta["d"] <- beta["t0"]*beta["d"]
      if (parnames[i] == "tau") {
        if (restr_tau == Inf) {
          beta[parnames[i]] <- exp(beta[parnames[i]])
        } else if (simult_conf) {
          beta[parnames[i]] <- pnorm(beta[parnames[i]])*(maxt0-beta["t0"])
        } else {
          beta[parnames[i]] <- restr_tau * pnorm(beta[parnames[i]])
        }
      }
    }

  }
  return(beta)
}

#   beta[parnames[i]] <- beta[parnames[i]]
# if (parnames[i] %in% names(fit_pars_columns)) parammatrix[,parnames[i]] <-
#     model_matrix[,fit_pars_columns[[parnames[i]]], drop=FALSE] %*%
#     beta[grepl(names(beta), pattern=paste0(parnames[i], "_"))]
