// [[Rcpp::depends(RcppArmadillo, nloptr)]]
#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <sstream>
#include <vector>
#include "densities/density_WEVmu.h"
// #include "densities/legacy_density_WEVmu.h"
#include "utils/optimization_context.hpp"
#include "utils/logger.hpp"

// these are used by the grid search to quickly reject bad parameter sets
// without calculating the likelihood for all trials
// can and should be fine-tuned for better overall performance
const int max_bad_trials = 15;           // absolute limit: 15 bad trials
const double max_bad_ratio = 0.15;       // relative limit: 15% bad trials
const int min_trials_for_ratio = 30;     // only check ratio after 30 trials
const double bailout_threshold = -1e6;   // logLikelihood threshold
const int bailout_check_after = 50;      // check likelihood after 50 trials

double calculate_neg_loglikelihood(
    OptimizationContext& optimization_context, 
    const arma::vec& beta,
    bool use_early_rejection = false
) {
    int n_trials = optimization_context.dependent_vars.n_rows;
    double total_logl = 0.0;
    
    int bad_trials = 0;

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

        if (prob == 0 || R_IsNaN(prob)) {
            bad_trials++;
            
            if (use_early_rejection) {
                if (bad_trials > max_bad_trials) {
                    if (optimization_context.logging) {
                        std::stringstream ss;
                        ss << "Early rejection: " << bad_trials 
                           << " bad trials (limit: " << max_bad_trials << ")";
                        Logger::warn(ss.str());
                    }
                    return 1e12;
                }
                
                if (i >= min_trials_for_ratio) {
                    double bad_ratio = (double)bad_trials / (i + 1);
                    if (bad_ratio > max_bad_ratio) {
                        if (optimization_context.logging) {
                            std::stringstream ss;
                            ss << "Early rejection: bad ratio " 
                               << std::fixed << std::setprecision(1)
                               << (bad_ratio * 100) << "% (limit: " 
                               << (max_bad_ratio * 100) << "%) at trial " << (i + 1);
                            Logger::warn(ss.str());
                        }
                        return 1e12;
                    }
                }
            }

            total_logl += std::log(std::numeric_limits<double>::min());
        } else {
            total_logl += std::log(prob);
        }
        
        if (use_early_rejection && i >= bailout_check_after) {
            if (total_logl < bailout_threshold) {
                if (optimization_context.logging) {
                    std::stringstream ss;
                    ss << "Early rejection: logLik " << total_logl 
                       << " too low (threshold: " << bailout_threshold 
                       << ") at trial " << (i + 1);
                    Logger::warn(ss.str());
                }
                return 1e12;
            }
        }
    }

    if (!std::isfinite(total_logl)) {
        Logger::error("Non-finite log-likelihood encountered");
        return 1e12;
    }

    return -total_logl;
}

// the objective function required by NLopt's C API
double neg_loglikelihood(
    unsigned n,
    const double* x,
    double* grad,
    void* func_data
) {
    // cast void* data back to our C++ context object
    OptimizationContext* optimization_context = static_cast<OptimizationContext*>(func_data);
    arma::vec beta(const_cast<double*>(x), n, false, true);

    double negLogLik = calculate_neg_loglikelihood(*optimization_context, beta);

    if (++(optimization_context->eval_count) % 50 == 0) {
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

// runs the main optimization for a single set of start_params
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
        
        // setup nlopt_opt (e.g., "Nelder-Mead" or "bobyqa")
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

        // set optimizer controls
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

// runs the grid search over a matrix of initial parameter sets
// [[Rcpp::export]]
Rcpp::NumericVector grid_search_worker(const Rcpp::List& optimization_context_dto, Rcpp::NumericMatrix inits_matrix) {
    OptimizationContext optimization_context(optimization_context_dto);
    Logger::init(optimization_context.logging);

    int n_sets = inits_matrix.nrow();
    Rcpp::NumericVector results(n_sets);

    arma::mat params_matrix(inits_matrix.begin(), inits_matrix.nrow(), inits_matrix.ncol(), false);

    for (int i = 0; i < n_sets; ++i) {
        arma::vec beta = arma::trans(params_matrix.row(i));
        results[i] = calculate_neg_loglikelihood(optimization_context, beta, true);
    }

    return results;
}