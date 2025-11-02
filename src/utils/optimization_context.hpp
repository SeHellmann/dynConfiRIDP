#pragma once

#include <RcppArmadillo.h>
#include <vector>
#include <regex>

const double MIN_THETA = -1e32;
const double MAX_THETA = 1e32;
// regex to find 'theta' params but not 'thetaLower' or 'thetaUpper'
const std::regex THETA_REGEX("^(d)?theta(?!Lower|Upper)");

class OptimizationContext;

// enum to map parameter names (strings) to a safe type
enum class ParamType {
    a, v, t0, d, sz, sv, st0, z, tau, lambda, w, muvis, sigvis, svis, s
};

/**
 * @struct ModelParameters
 * @brief Holds all model parameters for a SINGLE TRIAL.
 * This struct is populated by OptimizationContext::get_trial_params()
 * and contains the final, transformed parameter values ready for
 * the density function.
 */
struct ModelParameters {
    double a = 0.0, v = 0.0, t0 = 0.0, d = 0.0, sz = 0.0;
    double sv = 0.0, st0 = 0.0, z = 0.0, tau = 0.0;
    double lambda = 0.0, w = 0.0, muvis = 0.0, sigvis = 0.0, svis = 0.0, s = 1.0;
    
    double th1 = 0.0, th2 = 0.0;

    void apply_transformations(const OptimizationContext& optimization_context);

    // converts the struct into the raw Rcpp::NumericVector
    // required by the g_minus_WEVmu density function
    Rcpp::NumericVector to_density_vector(bool boundary, double precision) const;
};

/**
 * @class OptimizationContext
 * @brief The main data transfer object (DTO) from R to C++.
 *
 * This class holds all data for the ENTIRE optimization run.
 * It is initialized once from the R list and then passed to the
 * optimizer. Its main job is to efficiently calculate the
 * trial-specific ModelParameters for any given trial index.
 */
class OptimizationContext {
public:
    struct OptimizationOptions { int maxfun = 8000; double reltol = 1e-6; };

    const arma::mat dependent_vars;
    const std::string optim_method;
    const double maxt0;
    const double restr_tau;
    const double precision;
    const bool simult_conf;
    const bool logging;
    const OptimizationOptions opts;

    int eval_count = 0;
    double last_best_logl = std::numeric_limits<double>::infinity();

    explicit OptimizationContext(const Rcpp::List& optimization_context);

    ModelParameters get_trial_params(int trial_idx, const arma::vec& beta) const;
private:
    /**
     * These structs define the *parameter mapping* strategy.
     * They store the indices needed to build the trial-specific
     * 'ModelParameters' struct from the main 'beta'
     * optimization vector and the 'model_matrix'.
     */
    struct ThetaIndices {
        std::vector<int> sym;
        std::vector<int> lower;
        std::vector<int> upper;
    };

    struct FixedParam { ParamType type; double value; };
    struct EstimatedParam { ParamType type; int beta_index; };
    struct FormulaParam {
        ParamType type;
        arma::uvec beta_indices;
        arma::uvec model_matrix_indices;
    };

    const arma::mat model_matrix;
    const int n_ratings;
    const bool sym_thetas;
    const ThetaIndices theta_indices;

    const std::vector<FixedParam> fixed_params;
    const std::vector<EstimatedParam> estimated_params;
    const std::vector<FormulaParam> formula_params;

    // private parsers (called by constructor)
    static OptimizationOptions parse_opts(const Rcpp::List& optimization_context_dto);
    static ThetaIndices parse_theta_indices(const Rcpp::List& optimization_context_dto);
    static std::vector<FixedParam> parse_fixed_params(const Rcpp::List& optimization_context_dto);
    static std::vector<EstimatedParam> parse_estimated_params(const Rcpp::List& optimization_context_dto);
    static std::vector<FormulaParam> parse_formula_params(const Rcpp::List& optimization_context_dto);

    arma::vec calculate_thetas(const arma::vec& beta) const;
    arma::vec calculate_sym_thetas(const arma::vec& beta) const;
    arma::vec calculate_asym_thetas(const arma::vec& beta) const;

    // helper to get a direct pointer to a member of ModelParameters
    double* get_param_pointer(ModelParameters& params, ParamType type) const;
};