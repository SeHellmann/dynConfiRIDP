/* RNG_WEV.h - Functions for generating random trials in the WEV model
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

#pragma once
#include <Rcpp.h>

Rcpp::NumericMatrix RNG_WEV(int n, Rcpp::NumericVector params,
                            double delta = 0.01, double maxT = 9,
                            bool stop_on_error = true);

Rcpp::NumericMatrix RNG_WEV_matrix(Rcpp::NumericMatrix params,
                                   double delta = 0.01, double maxT = 9,
                                   bool stop_on_error = true);
