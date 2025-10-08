/* Density_IRM.h - Functions for PDF calculation in the independent race model
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

Rcpp::NumericVector density_IRM(Rcpp::NumericVector rts,
                                Rcpp::NumericVector params, int win = 1,
                                double step_width = 0.0001);
