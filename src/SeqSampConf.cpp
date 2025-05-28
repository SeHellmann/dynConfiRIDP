/* SeqSampConf.cpp - Main source file for the RCpp implementation of the
 * sequential sampling models contained in the dynConfiR package
 *
 * Copyright (C) 2022 Sebastian Hellmann.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 */


#include <Rcpp.h>
//#include <R_ext/Rdynload.h>
#include <iostream>
#include <sstream>

#include "SeqSampConf.h"

using namespace Rcpp;



// [[Rcpp::export]]
NumericVector d_WEVmu_fit (NumericVector rts, LogicalVector boundary,
                           NumericMatrix parammatrix, double precision=3)
{
  int length = rts.length();
  if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }
  NumericVector out(length, -10.0);  // Should default to 0s when creating NumericVector, but just in case..


  NumericVector params(20, 0.0);

  // Add tuning values for numerical integrations at the end of parameters
  // ToDo: Optimize and check precision values
  params[16] = 0.0089045 * exp(-1.037580*precision); // TUNE_INT_T0
  params[17] = 0.0508061 * exp(-1.022373*precision); // TUNE_INT_Z
  //     These have been added to optimise code paths by treating very small variances as 0
  //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
  params[18] = pow (10, -(precision+2.0)); // TUNE_SZ_EPSILON
  params[19] = pow (10, -(precision+2.0)); // TUNE_ST0_EPSILON

      // for each entry in rts and each row in parammatrix, define params as the row of parammatrix and do the rest
    for (int i = 0; i < length; i++)
    {
      params[Rcpp::Range(0, 15)] = parammatrix(i, _);
      if (!ValidateParams(params, true)) Rcpp::stop("Error validating parameters.\n");

      if (boundary[i]) {
        params[7] = 1- params[7]; // z -> 1 - z
        params[1] = - params[1]; // v  -> - v
        params[3] = - params[3]; // d  -> - d
      }
      out[i] = fabs(g_minus_WEVmu(rts[i], params));
    }
    return out;
}



// [[Rcpp::export]]
NumericMatrix r_WEV_matrix (NumericMatrix params,
                     double delta=0.01, double maxT=9,
                     bool stop_on_error=true)
{
  if (params.ncol()<13) { Rcpp::stop("Not enough parameters supplied.\n"); }
  NumericMatrix out(params.nrow(), 6);
  out = RNG_WEV_matrix(params, delta, maxT, stop_on_error);

  return out;
}











//// The following are just functions used for the other version, i.e.
//// To check, how I used the other models in the package


// R-callable PDF for 2DSD - pass boundary to retrieve (1 = lower, 2 = upper)
// [[Rcpp::export]]
NumericVector d_2DSD (NumericVector rts, NumericVector params, double precision=1e-5, int boundary=2, bool stop_on_error=true, int stop_on_zero=false)
{
  int length = rts.length();
  if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }

  if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

  if (!ValidateParams(params, true))
  {
    if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
    else { return out; }
  }

  // Add tuning values for numerical integrations at the end of parameters
  // ToDo: Optimize and check precision values
  if (precision >= 1) {
    params.push_back(0.0089045 * exp(-1.037580*precision)); // TUNE_INT_T0
    params.push_back(0.0508061 * exp(-1.022373*precision)); // TUNE_INT_Z
    //     These have been added to optimise code paths by treating very small variances as 0
    //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
    params.push_back(pow (10, -(precision+2.0))); // TUNE_SZ_EPSILON
    params.push_back(pow (10, -(precision+2.0))); // TUNE_ST0_EPSILON
  } else {
    params.push_back(precision); // TUNE_INT_T0
    params.push_back(precision); // TUNE_INT_Z
    params.push_back(0); // TUNE_SZ_EPSILON
    params.push_back(0); // TUNE_ST0_EPSILON
  }

  out = density_2DSD (rts, params, boundary-1, stop_on_zero);

  //delete g_Params;
  return out;
}


// [[Rcpp::export]]
NumericVector d_WEVmu (NumericVector rts, NumericMatrix parammatrix, double precision=1e-5, int boundary=2,
                       bool stop_on_error=true, int stop_on_zero=false)
{
  int length = rts.length();
  if (length > MAX_INPUT_VALUES) { Rcpp::stop("Number of RT values passed in exceeds maximum of %d.\n", MAX_INPUT_VALUES); }

  if ((boundary < 1) || (boundary > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..
  NumericVector params(20, 0.0);

  // Add tuning values for numerical integrations at the end of parameters
  // ToDo: Optimize and check precision values
  params[16] = 0.0089045 * exp(-1.037580*precision); // TUNE_INT_T0
  params[17] = 0.0508061 * exp(-1.022373*precision); // TUNE_INT_Z
  //     These have been added to optimise code paths by treating very small variances as 0
  //     e.g. with precision = 3, sv or sz values < 10^-5 are considered 0
  params[18] = pow (10, -(precision+2.0)); // TUNE_SZ_EPSILON
  params[19] = pow (10, -(precision+2.0)); // TUNE_ST0_EPSILON

  if (parammatrix.nrow()==1) {
    params[Rcpp::Range(0, 15)] = parammatrix(0, _);
    if (!ValidateParams(params, true))
    {
      if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
      else { return out; }
    }
    out = density_WEVmu (rts, params, boundary-1, stop_on_zero);
  } else {
    // for each entry in rts and each row in parammatrix, define params as the row of parammatrix and do the rest
    for (int i = 0; i < length; i++)
    {
      params[Rcpp::Range(0, 15)] = parammatrix(i, _);
      if (!ValidateParams(params, true))
      {
        if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
        else { return out; }
      }
      out[i] = density_WEVmu2 (rts[i], params, boundary-1, false);
    }
  }
  return out;
}


// [[Rcpp::export]]
NumericVector d_IRM (NumericVector rts, NumericVector params, int win=1,  double step_width = 0.0001)
{
  int length = rts.length();

  if (params.length()!= 12) {Rcpp::stop ("Wrong number of parameters given. (Must be 12)\n"); }
  if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

  out = density_IRM (rts, params, win, step_width);

  return out;
}



// [[Rcpp::export]]
NumericVector d_PCRM (NumericVector rts, NumericVector params, int win=1, double step_width = 0.0001)
{
  int length = rts.length();

  if (params.length()!= 12) {Rcpp::stop ("Wrong number of parameters given. (Must be 12)\n"); }
  if ((win < 1) || (win > 2)) { Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n"); }

  NumericVector out(length, 0.0);  // Should default to 0s when creating NumericVector, but just in case..

  out = density_PCRM (rts, params, win, step_width);

  return out;
}



// End CPP
