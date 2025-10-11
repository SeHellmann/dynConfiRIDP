// [[Rcpp::depends(RcppArmadillo, nloptr)]]
#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include "densities/density_WEVmu.h"
#include "utils/optimization_context.hpp"
#include "utils/logger.hpp"

static int eval_count = 0;
static double last_best_logl = std::numeric_limits<double>::infinity();

double neg_loglikelihood_formula(
    unsigned n,
    const double* x,
    double* grad,
    void* func_data
) {
    OptimizationContext* optimization_context = static_cast<OptimizationContext*>(func_data);
    arma::vec beta(const_cast<double*>(x), n, false, true);

    int n_trials = optimization_context->dependent_vars.n_rows;
    double total_logl = 0.0;

    for (int i = 0; i < n_trials; ++i) {
        if (i % 100 == 0) Rcpp::checkUserInterrupt();

        ModelParameters params = optimization_context->get_trial_params(i, beta);
        params.apply_transformations(*optimization_context);

        if (optimization_context->dependent_vars.n_cols > 3) {
            params.v *= optimization_context->dependent_vars(i, 3);
        }

        bool boundary = optimization_context->dependent_vars(i, 2);
        double rt = optimization_context->dependent_vars(i, 0);

        Rcpp::NumericVector fit_vector = params.to_density_vector(boundary);
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

    double negLogLik = -total_logl;
    if (++eval_count % 20 == 0) {
        std::stringstream ss;
        ss << "Evaluation " << eval_count 
           << " - negLogLik: " << negLogLik
           << " (previous best: " << last_best_logl << ")\n";
           
        Logger::info(ss.str());
        
        if (negLogLik < last_best_logl) {
            last_best_logl = negLogLik;
        }
    }

    return negLogLik;
}

// [[Rcpp::export]]
Rcpp::List nlopt_optimizer(const Rcpp::List& optimization_context_dto, Rcpp::NumericVector start_params) {
    try {
        OptimizationContext optimization_context(optimization_context_dto);
    
        Logger::init(Rcpp::as<bool>(optimization_context_dto["logging"]));

        std::string optim_method = Rcpp::as<std::string>(optimization_context_dto["optim_method"]);
        unsigned n_params = start_params.size();

        std::stringstream ss;
        ss << "Starting NLopt optimization with method " << optim_method 
           << " and " << n_params << " parameters";
        Logger::info(ss.str());
        
        nlopt_opt opt;
        if (optim_method == "Nelder-Mead") {
            opt = nlopt_create(NLOPT_LN_NELDERMEAD, n_params);

            std::vector<double> step(n_params, 0.02);
            nlopt_set_initial_step(opt, step.data());
        } else if (optim_method == "bobyqa") {
            opt = nlopt_create(NLOPT_LN_BOBYQA, n_params);

            std::vector<double> step(n_params, 2);
            nlopt_set_initial_step(opt, step.data());
        } else {
            Logger::error("Unsupported optimization method: " + optim_method);
            Rcpp::stop("Unsupported optimization method provided.");
        }

        nlopt_set_min_objective(opt, neg_loglikelihood_formula, &optimization_context);
        
        Rcpp::List opts = Rcpp::as<Rcpp::List>(optimization_context_dto["opts"]);
        double reltol = Rcpp::as<double>(opts["reltol"]);

        nlopt_set_ftol_rel(opt, reltol);
        nlopt_set_xtol_rel(opt, reltol);
        nlopt_set_ftol_abs(opt, 1e-10);
        nlopt_set_xtol_abs1(opt, 1e-10);
        nlopt_set_maxeval(opt, Rcpp::as<int>(opts["maxfun"]));
        nlopt_set_maxtime(opt, 3600);

        std::vector<double> x = Rcpp::as<std::vector<double>>(start_params);
        double minf;

        eval_count = 0;
        last_best_logl = std::numeric_limits<double>::infinity();
        
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
                Rcpp::Named("par") = Rcpp::NumericVector(n_params, NA_REAL)
            );
        }

        ss.str("");
        ss << "Optimization complete after " << eval_count << " evaluations"
           << " - Final negLogLik: " << minf;
        Logger::info(ss.str());

        Rcpp::NumericVector final_params = Rcpp::wrap(x);
        final_params.names() = start_params.names();
        
        return Rcpp::List::create(
            Rcpp::Named("value") = minf,
            Rcpp::Named("par") = final_params
        );

    } catch (const std::exception& e) {
        Logger::error(std::string("Error in nlopt_optimizer: ") + e.what());
        return Rcpp::List::create(
            Rcpp::Named("value") = NA_REAL,
            Rcpp::Named("par") = Rcpp::NumericVector(start_params.size(), NA_REAL)
        );
    }
}

// [[Rcpp::export]]
double grid_search_worker(const Rcpp::List& optimization_context_dto, Rcpp::NumericVector params) {
    OptimizationContext optimization_context(optimization_context_dto);
    return neg_loglikelihood_formula(params.size(), params.begin(), NULL, &optimization_context);
}