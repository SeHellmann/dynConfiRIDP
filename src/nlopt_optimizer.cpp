// [[Rcpp::depends(RcppArmadillo, nloptr)]]
#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <Rcpp.h>
#include <vector>
#include "densities/density_WEVmu.h"
#include "model_context.hpp"

double neg_loglikelihood_formula(
    unsigned n,
    const double* x,
    double* grad,
    void* func_data
) {
    Rcpp::checkUserInterrupt();

    EstimationContext* estimation_context = static_cast<EstimationContext*>(func_data);
    
    arma::vec beta(const_cast<double*>(x), n, false, true);

    int n_trials = estimation_context->dependent_vars.n_rows;
    double total_logl = 0.0;

    for (int i = 0; i < n_trials; ++i) {
        ModelParameters params = estimation_context->create_trial_params(i, beta);
        params.apply_transformations(*estimation_context);

        if (estimation_context->dependent_vars.n_cols > 3) {
            params.v *= estimation_context->dependent_vars(i, 3);
        }

        bool boundary = estimation_context->dependent_vars(i, 1);
        double rt = estimation_context->dependent_vars(i, 0);

        Rcpp::NumericVector fit_vector = params.to_density_vector(boundary, estimation_context->precision);

        double prob = std::abs(g_minus_WEVmu(rt, fit_vector));

        if (prob == 0) {
            total_logl += std::log(std::numeric_limits<double>::min());
        } else if (R_IsNaN(prob)) {
            return 1e12;
        } else {
            total_logl += std::log(prob);
        }
    }

    if (!std::isfinite(total_logl)) {
        return 1e12;
    }

    return -total_logl;
}

// [[Rcpp::export]]
Rcpp::List nlopt_optimizer(const Rcpp::List& context, Rcpp::NumericVector start_params) {
    EstimationContext estimation_context(context);
    
    unsigned n_params = start_params.size();
    nlopt_opt opt;

    std::string optim_method = Rcpp::as<std::string>(context["optim_method"]);
    if (optim_method == "Nelder-Mead") {
        opt = nlopt_create(NLOPT_LN_NELDERMEAD, n_params);
    } else if (optim_method == "bobyqa") {
        opt = nlopt_create(NLOPT_LN_BOBYQA, n_params);
    } else {
        // should not happen, already validated in utils_validate_args.R
        Rcpp::stop("Unsupported optimization method provided.");
    }

    nlopt_set_min_objective(opt, neg_loglikelihood_formula, &estimation_context);
    
    Rcpp::List opts = Rcpp::as<Rcpp::List>(context["opts"]);
    nlopt_set_xtol_rel(opt, Rcpp::as<double>(opts["reltol"]));
    nlopt_set_maxeval(opt, Rcpp::as<int>(opts["maxfun"]));

    std::vector<double> x = Rcpp::as<std::vector<double>>(start_params);
    double minf;

    nlopt_result result = nlopt_optimize(opt, x.data(), &minf);

    nlopt_destroy(opt);

    if (result < 0) {
        Rcpp::warning("NLopt failed with error code: %d", result);
        return Rcpp::List::create(
            Rcpp::Named("value") = NA_REAL,
            Rcpp::Named("par") = Rcpp::NumericVector(n_params, NA_REAL)
        );
    }

    Rcpp::NumericVector final_params = Rcpp::wrap(x);
    final_params.names() = start_params.names();

    return Rcpp::List::create(
        Rcpp::Named("value") = minf,
        Rcpp::Named("par") = final_params
    );
}

// [[Rcpp::export]]
double grid_search_worker(const Rcpp::List& context, Rcpp::NumericVector params) {
    EstimationContext estimation_context(context);
    return neg_loglikelihood_formula(params.size(), params.begin(), NULL, &estimation_context);
}