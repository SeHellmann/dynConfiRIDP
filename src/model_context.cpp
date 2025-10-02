#include "model_context.hpp"
#include "Rcpp/vector/instantiation.h"
#include "current/armadillo"
#include <cstddef>
#include <map>
#include <vector>

const double MIN_THETA = -1e32;
const double MAX_THETA = 1e32;

enum class ParamType { 
    a, v, t0, d, z, sz, sv, st0, tau, lambda, w, muvis, sigvis, svis, s, unknown 
};

static const std::map<std::string, ParamType> param_map = {
    {"a", ParamType::a}, {"v", ParamType::v}, {"t0", ParamType::t0}, 
    {"d", ParamType::d}, {"z", ParamType::z}, {"sz", ParamType::sz}, 
    {"sv", ParamType::sv}, {"st0", ParamType::st0}, {"tau", ParamType::tau}, 
    {"lambda", ParamType::lambda}, {"w", ParamType::w}, {"muvis", ParamType::muvis}, 
    {"sigvis", ParamType::sigvis}, {"svis", ParamType::svis}, {"s", ParamType::s}
};

static const std::vector<std::string> param_names = {
    "a", "v", "t0", "d", "z",
    "sz", "sv", "st0", "tau", "lambda", 
    "w", "muvis", "sigvis", "svis", "s"
};

EstimationContext::EstimationContext(const Rcpp::List& context) :
    dependent_vars(Rcpp::as<arma::mat>(context["dependent_vars"])),
    model_matrix(Rcpp::as<arma::mat>(context["model_matrix"])),
    beta_map(Rcpp::as<Rcpp::List>(context["beta_map"])),
    fixed(Rcpp::as<Rcpp::List>(context["fixed"])),
    beta_names(Rcpp::as<Rcpp::CharacterVector>(context["beta_names"])),
    maxt0(Rcpp::as<double>(context["maxt0"])),
    restr_tau(Rcpp::as<double>(context["restr_tau"])),
    precision(Rcpp::as<double>(context["precision"])),
    n_ratings(Rcpp::as<int>(context["n_ratings"])),
    simult_conf(Rcpp::as<bool>(context["simult_conf"])),
    sym_thetas(Rcpp::as<bool>(context["sym_thetas"]))
{
    for (int i = 0; i < beta_names.size(); ++i) {
        std::string name(beta_names[i]);
        if (name.find("thetaLower") != std::string::npos) {
            lower_theta_indices.push_back(i);
        } else if (name.find("thetaUpper") != std::string::npos) {
            upper_theta_indices.push_back(i);
        } else if (name.find("theta") != std::string::npos) {
            sym_theta_indices.push_back(i);
        }
    }
}

ModelParameters EstimationContext::create_trial_params(int trial_idx, const arma::vec& beta) const {
    ModelParameters params = {};

    arma::vec thetas = calculate_thetas(beta);

    std::map<std::string, double> beta_lookup;
    for (unsigned int i = 0; i < beta.n_elem; ++i) {
        beta_lookup[std::string(beta_names[i])] = beta[i];
    }

    for (const auto& param_name : param_names) {
        double* target_param = get_param_pointer(params, param_name);
        if (!target_param) continue;
        
        const char* param_name_str = param_name.c_str();

        if (fixed.containsElementNamed(param_name_str)) {
            *target_param = Rcpp::as<double>(fixed[param_name]);
            continue;
        }

        if (beta_map.containsElementNamed(param_name_str)) {
            Rcpp::IntegerVector beta_indices = beta_map[param_name];
            arma::uvec u_beta_indices = Rcpp::as<arma::uvec>(beta_indices) - 1;
            arma::rowvec trial_row = model_matrix.row(trial_idx);
            *target_param = arma::dot(trial_row.elem(u_beta_indices), beta.elem(u_beta_indices));
            continue;
        }

        auto it = beta_lookup.find(param_name);
        if (it != beta_lookup.end()) {
            *target_param = it->second;
        } else {
            *target_param = arma::datum::nan;
        }
    }

    if (!thetas.is_empty()) {
        int response = static_cast<int>(dependent_vars(trial_idx, 1));
        int rating = static_cast<int>(dependent_vars(trial_idx, 2));
        int idx = response * (n_ratings + 2) + rating;

        if (idx < static_cast<int>(thetas.n_elem) && (idx + 1) < static_cast<int>(thetas.n_elem)) {
            params.th1 = thetas(idx);
            params.th2 = thetas(idx + 1);
        }
    }

    return params;
}

arma::vec EstimationContext::calculate_thetas(const arma::vec& beta) const {
    return sym_thetas ? calculate_sym_thetas(beta) : calculate_asym_thetas(beta);
}

arma::vec EstimationContext::calculate_sym_thetas(const arma::vec& beta) const {
    if (sym_theta_indices.empty()) return arma::vec();

    std::vector<double> relevant_betas;
    for (int idx : sym_theta_indices) {
        relevant_betas.push_back(beta[idx]);
    }

    if (relevant_betas.empty()) return arma::vec();

    arma::vec thetas(n_ratings + 2);
    thetas[0] = MIN_THETA;
    thetas[1] = relevant_betas[0];
    
    size_t i = 1;
    for (; i < relevant_betas.size() && (i + 1) < static_cast<size_t>(n_ratings + 2); ++i) {
        thetas[i + 1] = thetas[i] + std::exp(relevant_betas[i]);
    }
    
    thetas[n_ratings + 1] = MAX_THETA;
    
    return arma::join_cols(thetas, thetas);
}

arma::vec EstimationContext::calculate_asym_thetas(const arma::vec& beta) const {
    auto process_theta_vector = [&](const std::vector<int>& indices) {
        arma::vec thetas(n_ratings + 2);
        thetas[0] = MIN_THETA;
        thetas[n_ratings + 1] = MAX_THETA;
        
        if (!indices.empty()) {
            thetas[1] = beta[indices[0]];
            for (size_t i = 1; i < indices.size() && (i + 1) < static_cast<size_t>(n_ratings + 2); ++i) {
                thetas[i + 1] = thetas[i] + std::exp(beta[indices[i]]);
            }
        }
        
        return thetas;
    };
    
    arma::vec lower_thetas = process_theta_vector(lower_theta_indices);
    arma::vec upper_thetas = process_theta_vector(upper_theta_indices);
    
    return arma::join_cols(lower_thetas, upper_thetas);
}

void ModelParameters::apply_transformations(const EstimationContext& context) {
    // fixed01
    z = R::pnorm(z, 0, 1, 1, 0);
    sz = R::pnorm(sz, 0, 1, 1, 0);
    w = R::pnorm(w, 0, 1, 1, 0);
    d = R::pnorm(d, 0, 1, 1, 0);
    t0 = R::pnorm(t0, 0, 1, 1, 0);
    st0 = R::pnorm(st0, 0, 1, 1, 0);
    
    // fixedpos
    a = std::exp(a);
    v = std::exp(v);
    sv = std::exp(sv);
    svis = std::exp(svis);
    sigvis = std::exp(sigvis);
    lambda = std::exp(lambda);
    s = std::exp(s);
    
    sz *= (2.0 * std::min(z, 1.0 - z));
    t0 *= context.maxt0;
    d *= t0;
    
    if (std::isinf(context.restr_tau)) {
        tau = std::exp(tau);
    } else if (context.simult_conf) {
        tau *= (context.maxt0 - t0);
    } else {
        tau *= context.restr_tau;
    }
    
    if (R_IsNaN(muvis)) {
        muvis = std::abs(v);
    }
    
    a /= s;
    v /= s;
    sv /= s;
    th1 /= s;
    th2 /= s;
    muvis /= s;
    sigvis /= s;
    svis /= s;
    
    st0 *= 2.0;
    t0 += st0 / 2.0;
}

Rcpp::NumericVector ModelParameters::to_density_vector(bool boundary, double precision) const {
    Rcpp::NumericVector p(20);
    
    p[0] = a;
    p[1] = v;
    p[2] = t0;
    p[3] = d;
    p[4] = sz;
    p[5] = sv;
    p[6] = st0;
    p[7] = z;
    p[8] = tau;
    p[9] = th1;
    p[10] = th2;
    p[11] = lambda;
    p[12] = w;
    p[13] = muvis;
    p[14] = sigvis;
    p[15] = svis;
    
    if (boundary) {
        p[7] = 1.0 - p[7]; // z -> 1 - z
        p[1] = -p[1];      // v -> -v
        p[3] = -p[3];      // d -> -d
    }
    
    p[16] = 0.0089045 * std::exp(-1.037580 * precision); // TUNE_INT_T0
    p[17] = 0.0508061 * std::exp(-1.022373 * precision); // TUNE_INT_Z
    p[18] = std::pow(10, -(precision + 2.0)); // TUNE_SZ_EPSILON
    p[19] = std::pow(10, -(precision + 2.0)); // TUNE_ST0_EPSILON
    
    return p;
}

double* EstimationContext::get_param_pointer(ModelParameters& params, const std::string& param_name) const {
    auto it = param_map.find(param_name);
    if (it == param_map.end()) return nullptr;
    
    switch (it->second) {
        case ParamType::a: return &params.a;
        case ParamType::v: return &params.v;
        case ParamType::t0: return &params.t0;
        case ParamType::d: return &params.d;
        case ParamType::z: return &params.z;
        case ParamType::sz: return &params.sz;
        case ParamType::sv: return &params.sv;
        case ParamType::st0: return &params.st0;
        case ParamType::tau: return &params.tau;
        case ParamType::lambda: return &params.lambda;
        case ParamType::w: return &params.w;
        case ParamType::muvis: return &params.muvis;
        case ParamType::sigvis: return &params.sigvis;
        case ParamType::svis: return &params.svis;
        case ParamType::s: return &params.s;
        default: return nullptr;
    }
}
