#include "validate_params.h"

#define CHECK_PARAM(cond, paramName, value, message)                           \
  if (!(cond)) {                                                               \
    valid = false;                                                             \
    logError(paramName, value, message);                                       \
  }

void logError(const std::string &paramName, double value,
              const std::string &message) {
  Rcpp::Rcout << "[param validation error]: parameter '" << paramName
              << "' with value " << value << " " << message << std::endl;
}

bool validate_params(Rcpp::NumericVector params) {
  bool valid = true;
  // a: must be > 0
  CHECK_PARAM(params[0] > 0, "a", params[0], "must be > 0")
  // sz: must be >= 0
  CHECK_PARAM(params[4] >= 0, "sz", params[4], "must be >= 0")
  // sv: must be >= 0
  CHECK_PARAM(params[5] >= 0, "sv", params[5], "must be >= 0")
  // st0: must be >= 0
  CHECK_PARAM(params[6] >= 0, "st0", params[6], "must be >= 0")
  // z: must be in [0, 1]
  CHECK_PARAM(params[7] >= 0 && params[7] <= 1, "z", params[7],
              "must be in [0, 1]")

  // t0 < 0.5*(d+st0)
  if (params[2] - fabs(0.5 * params[3]) - 0.5 * params[6] < 0) {
    valid = false;
    Rcpp::Rcout
        << "[param validation error]: invalid parameter combination t0 = "
        << params[2] << ", d = " << params[3] << ", st0 =" << params[6]
        << ", t0 < 0.5*(d+st0) must hold" << std::endl;
  }
  // z - sz/2 >= 0 AND z + sz/2 <= 1
  if (params[7] - 0.5 * params[4] < 0 || params[7] + 0.5 * params[4] > 1) {
    valid = false;
    Rcpp::Rcout << "[param validation error]: z = " << params[7]
                << ", sz = " << params[4]
                << " must satisfy: z - sz/2 >= 0 AND z + sz/2 <= 1"
                << std::endl;
  }

  if (params.size() <= 10) { // DDMConf model
    // th2 < th1
    if (params[9] < params[8]) {
      valid = false;
      Rcpp::Rcout
          << "[param validation error]: invalid parameter combination th1 = "
          << params[8] << ", th2 = " << params[9] << ", th2 must be < th1"
          << std::endl;
    }
  } else { // 2DSDT and dynaViTE
    // tau: must be >= 0
    CHECK_PARAM(params[8] >= 0, "tau", params[8], "must be >= 0")

    // th2 < th1
    if (params[10] < params[9]) {
      valid = false;
      Rcpp::Rcout
          << "[param validation error]: invalid parameter combination th1 = "
          << params[9] << ", th2 = " << params[10] << ", th2 must be < th1"
          << std::endl;
    }

    // lambda: must be >= 0
    CHECK_PARAM(params[11] >= 0, "lambda", params[11], "must be >= 0")

    if (params.size() > 12) { // only dynaViTE
      // w: must be in [0, 1]
      CHECK_PARAM(params[12] >= 0 && params[12] <= 1, "w", params[12],
                  "must be in [0, 1]")
      // sigvis: must be >= 0
      CHECK_PARAM(params[14] >= 0, "sigvis", params[14], "must be >= 0")
      // svis: must be > 0
      CHECK_PARAM(params[15] > 0, "svis", params[15], "must be > 0")
    }
  }

  return valid;
}
