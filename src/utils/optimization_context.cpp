#include "optimization_context.hpp"
#include "validate_params.h"
#include <map>
#include <string>
#include <regex>

static const std::map<std::string, ParamType> param_map = {
    {"a", ParamType::a}, {"v", ParamType::v}, {"t0", ParamType::t0},
    {"d", ParamType::d}, {"sz", ParamType::sz}, {"sv", ParamType::sv},
    {"st0", ParamType::st0}, {"z", ParamType::z}, {"tau", ParamType::tau}, 
    {"lambda", ParamType::lambda}, {"w", ParamType::w}, {"muvis", ParamType::muvis}, 
    {"sigvis", ParamType::sigvis}, {"svis", ParamType::svis}, {"s", ParamType::s}
};

OptimizationContext::OptimizationContext(const Rcpp::List& optimization_context) :
    dependent_vars(Rcpp::as<arma::mat>(optimization_context["dependent_vars"])),
    maxt0(Rcpp::as<double>(optimization_context["maxt0"])),
    restr_tau(Rcpp::as<double>(optimization_context["restr_tau"])),
    precision(Rcpp::as<double>(optimization_context["precision"])),
    simult_conf(Rcpp::as<bool>(optimization_context["simult_conf"])),
    model_matrix(Rcpp::as<arma::mat>(optimization_context["model_matrix"])),
    n_ratings(Rcpp::as<int>(optimization_context["n_ratings"])),
    sym_thetas(Rcpp::as<bool>(optimization_context["sym_thetas"]))
{
    Rcpp::List fixed_list = optimization_context["fixed_params"];
    Rcpp::CharacterVector fixed_names = fixed_list.names();
    for (int i = 0; i < fixed_list.size(); ++i) {
        auto it = param_map.find(std::string(fixed_names[i]));
        if (it != param_map.end()) {
            fixed_params.push_back({ it->second, Rcpp::as<double>(fixed_list[i]) });
        }
    }

    Rcpp::List estimated_list = optimization_context["estimated_params"];
    Rcpp::CharacterVector estimated_names = estimated_list.names();
    for (int i = 0; i < estimated_list.size(); ++i) {
        auto it = param_map.find(std::string(estimated_names[i]));
        if (it != param_map.end()) {
            estimated_params.push_back({ it->second, Rcpp::as<int>(estimated_list[i]) - 1 });
        }
    }

    Rcpp::List formula_list = optimization_context["formula_params"];
    Rcpp::CharacterVector formula_names = formula_list.names();
    for (int i = 0; i < formula_list.size(); ++i) {
        auto it = param_map.find(std::string(formula_names[i]));
        if (it != param_map.end()) {
            FormulaParam formula_param = { it->second };

            Rcpp::List components = Rcpp::as<Rcpp::List>(formula_list[i]);
            if (components.containsElementNamed("beta_indices")) {
                formula_param.beta_indices = Rcpp::as<arma::uvec>(components["beta_indices"]) - 1;
            }

            if (components.containsElementNamed("model_matrix_indices")) {
                formula_param.model_matrix_indices = Rcpp::as<arma::uvec>(components["model_matrix_indices"]) - 1;
            }

            formula_params.push_back(formula_param);
        }
    }

    const Rcpp::CharacterVector beta_names = optimization_context["beta_names"];
    const std::regex sym_theta_regex("^(d)?theta(?!Lower|Upper)");
    for (int i = 0; i < beta_names.size(); ++i) {
        std::string name(beta_names[i]);
        if (name.find("thetaLower") != std::string::npos) {
            lower_theta_indices.push_back(i);
        } else if (name.find("thetaUpper") != std::string::npos) {
            upper_theta_indices.push_back(i);
        } else if (std::regex_search(name, sym_theta_regex)) {
            sym_theta_indices.push_back(i);
        }
    }
}

ModelParameters OptimizationContext::get_trial_params(int trial_idx, const arma::vec& beta) const {
    ModelParameters params = {};

    for (const auto& fixed_param : fixed_params) {
        *get_param_pointer(params, fixed_param.type) = fixed_param.value;
    }

    for (const auto& estimated_param : estimated_params) {
        *get_param_pointer(params, estimated_param.type) = beta[estimated_param.beta_index];
    }

    if (!formula_params.empty()) {
        arma::rowvec trial_row = model_matrix.row(trial_idx);
        for (const auto& formula_param : formula_params) {
            double accumulated_value = arma::dot(
                trial_row.elem(formula_param.model_matrix_indices),
                beta.elem(formula_param.beta_indices)
            );
            
            *get_param_pointer(params, formula_param.type) = accumulated_value;
        }
    }

    arma::vec thetas = calculate_thetas(beta);
    if (!thetas.is_empty()) {
        int rating = static_cast<int>(dependent_vars(trial_idx, 1));
        int response = static_cast<int>(dependent_vars(trial_idx, 2));
        
        int block_size = n_ratings + 1;
        int offset = response * block_size;

        int th1_idx = offset + rating - 1;
        int th2_idx = offset + rating;

        if (th2_idx < static_cast<int>(thetas.n_elem)) {
            params.th1 = thetas(th1_idx);
            params.th2 = thetas(th2_idx);
        }
    }
    
    if (params.muvis == 0.0) {
        bool was_set = false;
        for (const auto& instr : fixed_params) if (instr.type == ParamType::muvis) was_set = true;
        if (!was_set) {
          for (const auto& instr : estimated_params) if (instr.type == ParamType::muvis) was_set = true;
        }
        if(!was_set) params.muvis = arma::datum::nan;
    }
    
    return params;
}

arma::vec OptimizationContext::calculate_thetas(const arma::vec& beta) const {
    return sym_thetas ? calculate_sym_thetas(beta) : calculate_asym_thetas(beta);
}

arma::vec OptimizationContext::calculate_sym_thetas(const arma::vec& beta) const {
    if (sym_theta_indices.empty()) return arma::vec();

    std::vector<double> relevant_betas;
    for (int idx : sym_theta_indices) {
        relevant_betas.push_back(beta[idx]);
    }

    if (relevant_betas.empty()) return arma::vec();

    arma::vec thetas(n_ratings + 1);
    thetas[0] = MIN_THETA;
    thetas[n_ratings] = MAX_THETA;

    if(relevant_betas.size() > 0) thetas[1] = relevant_betas[0];

    for (size_t i = 1; i < relevant_betas.size(); ++i) {
        thetas[i + 1] = thetas[i] + std::exp(relevant_betas[i]);
    }

    return arma::join_cols(thetas, thetas);
}

arma::vec OptimizationContext::calculate_asym_thetas(const arma::vec& beta) const {
    auto process_theta_vector = [&](const std::vector<int>& indices) {
        arma::vec thetas(n_ratings + 1);
        thetas[0] = MIN_THETA;
        thetas[n_ratings] = MAX_THETA;

        if (!indices.empty()) {
            std::vector<double> relevant_betas;
            for (int idx : indices) {
                relevant_betas.push_back(beta[idx]);
            }

            if(relevant_betas.size() > 0) thetas[1] = relevant_betas[0];

            for (size_t i = 1; i < relevant_betas.size(); ++i) {
                thetas[i + 1] = thetas[i] + std::exp(relevant_betas[i]);
            }
        }

        return thetas;
    };

    arma::vec lower_thetas = process_theta_vector(lower_theta_indices);
    arma::vec upper_thetas = process_theta_vector(upper_theta_indices);
    
    return arma::join_cols(lower_thetas, upper_thetas);
}

void ModelParameters::apply_transformations(const OptimizationContext& optimization_context) {
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
    t0 *= optimization_context.maxt0;
    d *= t0;

    if (std::isinf(optimization_context.restr_tau)) {
        tau = std::exp(tau);
    } else if (optimization_context.simult_conf) {
        tau = R::pnorm(tau, 0, 1, 1, 0) * (optimization_context.maxt0 - t0);
    } else {
        tau = optimization_context.restr_tau * R::pnorm(tau, 0, 1, 1, 0);
    }

    if (R_IsNaN(muvis)) muvis = std::abs(v);

    a /= s;
    v /= s;
    sv /= s;
    th1 /= s;
    th2 /= s;
    muvis /= s;
    sigvis /= s; 
    svis /= s;

    t0 += st0;
    st0 *= 2.0;
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
    
    p[16] = 0.0089045 * std::exp(-1.037580 * precision); // TUNE_INT_T0
    p[17] = 0.0508061 * std::exp(-1.022373 * precision); // TUNE_INT_Z
    p[18] = std::pow(10, -(precision + 2.0));       // TUNE_SZ_EPSILON
    p[19] = std::pow(10, -(precision + 2.0));       // TUNE_ST0_EPSILON

    if(!validate_params(p)) Rcpp::stop("Error validating params.\n");
    
    if (boundary) {
        p[7] = 1.0 - p[7]; // z -> 1 - z
        p[1] = -p[1];      // v -> -v
        p[3] = -p[3];      // d -> -d
    }

    return p;
}

double* OptimizationContext::get_param_pointer(ModelParameters& params, ParamType type) const {
    switch (type) {
        case ParamType::a: return &params.a;
        case ParamType::v: return &params.v;
        case ParamType::t0: return &params.t0;
        case ParamType::d: return &params.d;
        case ParamType::sz: return &params.sz;
        case ParamType::sv: return &params.sv;
        case ParamType::st0: return &params.st0;
        case ParamType::z: return &params.z;
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