// [[Rcpp::depends(RcppArmadillo, nloptr)]]
#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <sstream>
#include <vector>
#include "densities/density_WEVmu.h"
#include "utils/optimization_context.hpp"
#include "utils/logger.hpp"

double calculate_neg_loglikelihood(OptimizationContext& optimization_context, const arma::vec& beta) {
    int n_trials = optimization_context.dependent_vars.n_rows;
    double total_logl = 0.0;

    for (int i = 0; i < n_trials; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        ModelParameters params = optimization_context.get_trial_params(i, beta);
        params.apply_transformations(optimization_context);

        if (optimization_context.dependent_vars.n_cols > 3) {
            params.v *= optimization_context.dependent_vars(i, 3);
        }

        bool boundary = optimization_context.dependent_vars(i, 2);
        double rt = optimization_context.dependent_vars(i, 0);

        Rcpp::NumericVector fit_vector = params.to_density_vector(boundary, optimization_context.precision);
        double prob = std::abs(g_minus_WEVmu(rt, fit_vector));

        if (prob == 0) {
            total_logl += std::log(std::numeric_limits<double>::min());
        } else if (R_IsNaN(prob)) {
            Logger::error("NaN probability at trial " + std::to_string(i));
            return 1e12;
        } else {
            total_logl += std::log(prob);
        }
    }

    if (!std::isfinite(total_logl)) {
        Logger::error("Non-finite log-likelihood encountered");
        return 1e12;
    }

    return -total_logl;
}

double neg_loglikelihood(
    unsigned n,
    const double* x,
    double* grad,
    void* func_data
) {
    OptimizationContext* optimization_context = static_cast<OptimizationContext*>(func_data);
    arma::vec beta(const_cast<double*>(x), n, false, true);

    double negLogLik = calculate_neg_loglikelihood(*optimization_context, beta);

    if (++(optimization_context->eval_count) % 20 == 0) {
        std::stringstream ss;
        ss << "Evaluation " << optimization_context->eval_count 
           << " - negLogLik: " << negLogLik
           << " (previous best: " << optimization_context->last_best_logl << ")\n";
           
        Logger::info(ss.str());
        
        if (negLogLik < optimization_context->last_best_logl) {
            optimization_context->last_best_logl = negLogLik;
        }
    }

    return negLogLik;
}

// [[Rcpp::export]]
Rcpp::List nlopt_optimizer(const Rcpp::List& optimization_context_dto, Rcpp::NumericVector start_params) {
    try {
        OptimizationContext optimization_context(optimization_context_dto);
        Logger::init(optimization_context.logging);

        unsigned n_params = start_params.size();

        std::stringstream ss;
        ss << "Starting NLopt optimization with method " << optimization_context.optim_method 
           << " and " << n_params << " parameters";
        Logger::info(ss.str());
        
        nlopt_opt opt;
        if (optimization_context.optim_method == "Nelder-Mead") {
            opt = nlopt_create(NLOPT_LN_NELDERMEAD, n_params);

            std::vector<double> step(n_params, 0.02);
            nlopt_set_initial_step(opt, step.data());
        } else if (optimization_context.optim_method == "bobyqa") {
            opt = nlopt_create(NLOPT_LN_BOBYQA, n_params);

            std::vector<double> step(n_params, 2);
            nlopt_set_initial_step(opt, step.data());
        } else {
            Logger::error("Unsupported optimization method: " + optimization_context.optim_method);
            Rcpp::stop("Unsupported optimization method provided.");
        }

        nlopt_set_min_objective(opt, neg_loglikelihood, &optimization_context);

        nlopt_set_ftol_rel(opt, optimization_context.opts.reltol);
        nlopt_set_xtol_rel(opt, optimization_context.opts.reltol);
        nlopt_set_ftol_abs(opt, 1e-10);
        nlopt_set_xtol_abs1(opt, 1e-10);
        nlopt_set_maxeval(opt, optimization_context.opts.maxfun);

        std::vector<double> x = Rcpp::as<std::vector<double>>(start_params);
        double minf;
        
        ss.str("");
        ss << "Initial parameters: ";
        for (size_t i = 0; i < x.size(); ++i) {
            ss << x[i] << " ";
        }
        Logger::info(ss.str());

        nlopt_result result = nlopt_optimize(opt, x.data(), &minf);
        nlopt_destroy(opt);

        if (result < 0) {
            ss.str("");
            ss << "NLopt failed with error code: " << result;
            Logger::error(ss.str());
            
            return Rcpp::List::create(
                Rcpp::Named("value") = NA_REAL,
                Rcpp::Named("params") = NULL
            );
        }

        ss.str("");
        ss << "Optimization complete after " << optimization_context.eval_count << " evaluations"
           << " - Final negLogLik: " << minf;
        Logger::info(ss.str());

        Rcpp::NumericVector final_params = Rcpp::wrap(x);
        final_params.names() = start_params.names();
        
        return Rcpp::List::create(
            Rcpp::Named("value") = minf,
            Rcpp::Named("params") = final_params
        );

    } catch (const std::exception& e) {
        Logger::error(std::string("Error in nlopt_optimizer: ") + e.what());
        return Rcpp::List::create(
            Rcpp::Named("value") = NA_REAL,
            Rcpp::Named("params") = NULL
        );
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector grid_search_worker(const Rcpp::List& optimization_context_dto, Rcpp::NumericMatrix inits_matrix) {
    OptimizationContext optimization_context(optimization_context_dto);
    Logger::init(optimization_context.logging);

    int n_sets = inits_matrix.nrow();
    Rcpp::NumericVector results(n_sets);

    arma::mat params_matrix(inits_matrix.begin(), inits_matrix.nrow(), inits_matrix.ncol(), false);

    for (int i = 0; i < n_sets; ++i) {
        arma::vec beta = arma::trans(params_matrix.row(i));
        results[i] = calculate_neg_loglikelihood(optimization_context, beta);

        if ((i + 1) % 20 == 0 || (i + 1) == n_sets) {
            std::stringstream ss;
            ss << "Grid Search Batch Progress: [" 
               << (i + 1) << " / " << n_sets 
               << "] - Last negLogLik: " << results[i];
            
            Logger::info(ss.str());
        }
    }

    return results;
}