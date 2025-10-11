# "2DSDT", "DDMConf" ??
DYNWEV_MODELS <- c("2DSD", "dynWEV", "dynaViTE")
RM_MODELS <- c("IRM", "PCRM", "IRMt", "PCRMt")
MODELS <- c(DYNWEV_MODELS, RM_MODELS)
DYNWEV_PARNAMES <- c("a", "v", "t0", "d", "sz", "sv", "st0", "z", "tau", "lambda", "w", "muvis", "sigvis", "svis", "s")
RM_PARNAMES <- c("a", "b", "mu1", "mu2", "s1", "s2", "rho", "wx", "wrt", "wint", "t0", "st0")
OPTIM_METHODS <- c("Nelder-Mead", "bobyqa") # "L-BFGS-B"
DEFAULT_OPTS <- list(
  n_attempts = 5,
  n_restarts = 5,
  maxfun = 8000,
  maxit = 2000,
  reltol = 1e-6,
  factr = 1e-10
)

utils::globalVariables(c(
  "DYNWEV_MODELS",
  "RM_MODELS",
  "MODELS",
  "DYNWEV_PARNAMES",
  "RM_PARNAMES",
  "OPTIM_METHODS",
  "DEFAULT_OPTS",
  "PARALLEL_MODES"
))
