#include "optimization_context.hpp"
#include "validate_params.h"
#include <map>

static const std::map<std::string, ParamType> param_map = {
    {"a", ParamType::a}, {"v", ParamType::v}, {"t0", ParamType::t0},
    {"d", ParamType::d}, {"sz", ParamType::sz}, {"sv", ParamType::sv},
    {"st0", ParamType::st0}, {"z", ParamType::z}, {"tau", ParamType::tau}, 
    {"lambda", ParamType::lambda}, {"w", ParamType::w}, {"muvis", ParamType::muvis}, 
    {"sigvis", ParamType::sigvis}, {"svis", ParamType::svis}, {"s", ParamType::s}
};

OptimizationContext::OptimizationOptions OptimizationContext::parse_opts(const Rcpp::List& optimization_context_dto) {
    OptimizationOptions opts;
    if (optimization_context_dto.containsElementNamed("opts")) {
        Rcpp::List opts_list = Rcpp::as<Rcpp::List>(optimization_context_dto["opts"]);
        if (opts_list.containsElementNamed("maxfun")) {
            opts.maxfun = Rcpp::as<int>(opts_list["maxfun"]);
        }
        if (opts_list.containsElementNamed("reltol")) {
            opts.reltol = Rcpp::as<double>(opts_list["reltol"]);
        }
    }
    return opts;
}

OptimizationContext::ThetaIndices OptimizationContext::parse_theta_indices(const Rcpp::List& optimization_context_dto) {
    ThetaIndices theta_indices;
    if (optimization_context_dto.containsElementNamed("beta_names")) {
        const Rcpp::CharacterVector beta_names_list = optimization_context_dto["beta_names"];
        for (int i = 0; i < beta_names_list.size(); ++i) {
            std::string name(beta_names_list[i]);
            if (name.find("thetaLower") != std::string::npos) {
                theta_indices.lower.push_back(i);
            } else if (name.find("thetaUpper") != std::string::npos) {
                theta_indices.upper.push_back(i);
            } else if (std::regex_search(name, THETA_REGEX)) {
                theta_indices.sym.push_back(i);
            }
        }
    }
    return theta_indices;
}

std::vector<OptimizationContext::FixedParam> OptimizationContext::parse_fixed_params(const Rcpp::List& optimization_context_dto) {
    std::vector<OptimizationContext::FixedParam> fixed_params;
    if (optimization_context_dto.containsElementNamed("fixed_params")) {
        Rcpp::List fixed_params_list = optimization_context_dto["fixed_params"];
        if (fixed_params_list.size() > 0) {
            Rcpp::CharacterVector fixed_params_names = fixed_params_list.names();
            for (int i = 0; i < fixed_params_list.size(); ++i) {
                auto it = param_map.find(std::string(fixed_params_names[i]));
                if (it != param_map.end()) {
                    fixed_params.push_back({ it->second, Rcpp::as<double>(fixed_params_list[i]) });
                }
            }
        }
    }
    return fixed_params;
}

std::vector<OptimizationContext::EstimatedParam> OptimizationContext::parse_estimated_params(const Rcpp::List& optimization_context_dto) {
    std::vector<OptimizationContext::EstimatedParam> estimated_params;
    if (optimization_context_dto.containsElementNamed("estimated_params")) {
        Rcpp::List estimated_params_list = optimization_context_dto["estimated_params"];
        if (estimated_params_list.size() > 0) {
            Rcpp::CharacterVector estimated_params_names = estimated_params_list.names();
            for (int i = 0; i < estimated_params_list.size(); ++i) {
                auto it = param_map.find(std::string(estimated_params_names[i]));
                if (it != param_map.end()) {
                    estimated_params.push_back({
                        it->second, 
                        Rcpp::as<int>(estimated_params_list[i]) - 1
                    });
                }
            }
        }
    }    
    return estimated_params;
}

std::vector<OptimizationContext::FormulaParam> OptimizationContext::parse_formula_params(const Rcpp::List& optimization_context_dto) {
    std::vector<OptimizationContext::FormulaParam> formula_params;
    if (optimization_context_dto.containsElementNamed("formula_params")) {
        Rcpp::List formula_params_list = optimization_context_dto["formula_params"];
        if (formula_params_list.size() > 0) {
            Rcpp::CharacterVector formula_names = formula_params_list.names();
            for (int i = 0; i < formula_params_list.size(); ++i) {
                auto it = param_map.find(std::string(formula_names[i]));
                if (it != param_map.end()) {
                    FormulaParam formula_param = { it->second };

                    Rcpp::List components = Rcpp::as<Rcpp::List>(formula_params_list[i]);
                    if (components.containsElementNamed("beta_indices")) {
                        formula_param.beta_indices = Rcpp::as<arma::uvec>(components["beta_indices"]) - 1;
                    }

                    if (components.containsElementNamed("model_matrix_indices")) {
                        formula_param.model_matrix_indices = Rcpp::as<arma::uvec>(components["model_matrix_indices"]) - 1;
                    }

                    formula_params.push_back(formula_param);
                }
            }
        }
    }
    return formula_params;
}


OptimizationContext::OptimizationContext(const Rcpp::List& optimization_context_dto) :
    dependent_vars(Rcpp::as<arma::mat>(optimization_context_dto["dependent_vars"])),
    optim_method(Rcpp::as<std::string>(optimization_context_dto["optim_method"])),
    maxt0(Rcpp::as<double>(optimization_context_dto["maxt0"])),
    restr_tau(Rcpp::as<double>(optimization_context_dto["restr_tau"])),
    precision(Rcpp::as<double>(optimization_context_dto["precision"])),
    simult_conf(Rcpp::as<bool>(optimization_context_dto["simult_conf"])),
    logging(Rcpp::as<bool>(optimization_context_dto["logging"])),
    opts(parse_opts(optimization_context_dto)),
    model_matrix(Rcpp::as<arma::mat>(optimization_context_dto["model_matrix"])),
    n_ratings(Rcpp::as<int>(optimization_context_dto["n_ratings"])),
    sym_thetas(Rcpp::as<bool>(optimization_context_dto["sym_thetas"])),
    theta_indices(parse_theta_indices(optimization_context_dto)),
    fixed_params(parse_fixed_params(optimization_context_dto)),
    estimated_params(parse_estimated_params(optimization_context_dto)),
    formula_params(parse_formula_params(optimization_context_dto))
{}

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
    if (theta_indices.sym.empty()) return arma::vec();

    std::vector<double> relevant_betas;
    for (int idx : theta_indices.sym) {
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

    arma::vec lower_thetas = process_theta_vector(theta_indices.lower);
    arma::vec upper_thetas = process_theta_vector(theta_indices.upper);
    
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