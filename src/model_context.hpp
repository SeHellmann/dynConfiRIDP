#pragma once

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <string>

class EstimationContext;

struct ModelParameters {
    double a = 0.0, v = 0.0, t0 = 0.0, d = 0.0, z = 0.0;
    double sz = 0.0, sv = 0.0, st0 = 0.0, tau = 0.0;
    double lambda = 0.0, w = 0.0, muvis = 0.0, sigvis = 0.0, svis = 0.0, s = 1.0;
    
    double th1 = 0.0, th2 = 0.0;

    void apply_transformations(const EstimationContext &context);
    Rcpp::NumericVector to_density_vector(bool boundary, double precision) const;
};

class EstimationContext {
public:
    explicit EstimationContext(const Rcpp::List& context);

    ModelParameters create_trial_params(int trial_idx, const arma::vec& beta) const;
    arma::vec calculate_thetas(const arma::vec& beta) const;

    const arma::mat dependent_vars;
    const arma::mat model_matrix;
    const Rcpp::List beta_map;
    const Rcpp::List fixed;
    const Rcpp::CharacterVector beta_names;
    const double maxt0;
    const double restr_tau;
    const double precision;
    const int n_ratings;
    const bool simult_conf;
    const bool sym_thetas;

private:
    arma::vec calculate_sym_thetas(const arma::vec& beta) const;
    arma::vec calculate_asym_thetas(const arma::vec& beta) const;
    double* get_param_pointer(ModelParameters& params, const std::string& param_name) const;

    std::vector<int> sym_theta_indices;
    std::vector<int> lower_theta_indices;
    std::vector<int> upper_theta_indices;
};