// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <Rcpp.h>
#include <vector>
#include "densities/density_WEVmu.h"

struct ObjectiveData {
    arma::mat dependent_vars;
    arma::mat model_matrix;
    Rcpp::List beta_map;
    Rcpp::List fixed;
    Rcpp::CharacterVector beta_names;
    double maxt0;
    double restr_tau;
    double precision;
    int n_ratings;
    bool simult_conf;
    bool sym_thetas;
};

struct ModelParameters {
    inline static const std::vector<std::string> parnames = {
        "a", "v", "t0", "d", "z", 
        "sz", "sv", "st0", "tau", "th1", 
        "th2", "lambda", "w", "muvis", "sigvis", 
        "svis", "s"
    };

    double a, v, t0, d, z, sz, sv, st0, tau, th1, th2, lambda, w, muvis, sigvis, svis, s;

    Rcpp::NumericVector to_vector(bool boundary, double precision) const {
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

        p[16] = 0.0089045 * exp(-1.037580 * precision); // TUNE_INT_T0
        p[17] = 0.0508061 * exp(-1.022373 * precision); // TUNE_INT_Z
        p[18] = pow(10, -(precision + 2.0));       // TUNE_SZ_EPSILON
        p[19] = pow(10, -(precision + 2.0));       // TUNE_ST0_EPSILON
        
        return p;
    }
};

enum Param { a, v, t0, d, z, sz, sv, st0, tau, lambda, w, muvis, sigvis, svis, s, unknown };
std::map<std::string, Param> param_map = {
    {"a", a}, {"v", v}, {"t0", t0}, {"d", d}, {"z", z},
    {"sz", sz}, {"sv", sv}, {"st0", st0}, {"tau", tau}, {"lambda", lambda},
    {"w", w}, {"muvis", muvis}, {"sigvis", sigvis}, {"svis", svis}, {"s", s}
};

arma::vec calculate_thetas(const arma::vec& beta, const ObjectiveData* objective_data) {
    arma::vec thetas;
    std::vector<double> relevant_betas;
    Rcpp::CharacterVector beta_names = objective_data->beta_names;

    if (objective_data->sym_thetas) {
        for (int i = 0; i < beta_names.size(); ++i) {
            std::string name(beta_names[i]);
            if (name.find("theta") != std::string::npos) {
                relevant_betas.push_back(beta[i]);
            }
        }
        if (relevant_betas.empty()) return thetas;

        arma::vec cumsum_vec(relevant_betas.size() - 1);
        for (size_t i = 0; i < cumsum_vec.n_elem; ++i) {
            cumsum_vec[i] = std::exp(relevant_betas[i + 1]);
        }

        thetas.set_size(objective_data->n_ratings + 2);
        thetas[0] = -1e32;
        thetas[1] = relevant_betas[0];
        for (size_t i = 0; i < cumsum_vec.n_elem; ++i) {
            thetas[i + 2] = thetas[i + 1] + cumsum_vec[i];
        }
        thetas[objective_data->n_ratings + 1] = 1e32;

        arma::vec full_thetas = arma::join_cols(thetas, thetas);
        return full_thetas;
    } else {
        std::vector<double> lower_betas, upper_betas;
        for (int i = 0; i < beta_names.size(); ++i) {
            std::string name(beta_names[i]);
            if (name.find("thetaLower") != std::string::npos) lower_betas.push_back(beta[i]);
            if (name.find("thetaUpper") != std::string::npos) upper_betas.push_back(beta[i]);
        }

        auto process_asym_thetas = [&](const std::vector<double>& betas_in) {
            arma::vec thetas_out(objective_data->n_ratings + 2);
            thetas_out[0] = -1e21;
            if (!betas_in.empty()) {
                thetas_out[1] = betas_in[0];
                for (size_t i = 1; i < betas_in.size(); ++i) {
                    thetas_out[i + 1] = thetas_out[i] + std::exp(betas_in[i]);
                }
            }
            thetas_out[objective_data->n_ratings + 1] = 1e21;
            return thetas_out;
        };

        arma::vec thetas_lower = process_asym_thetas(lower_betas);
        arma::vec thetas_upper = process_asym_thetas(upper_betas);
        arma::vec full_thetas = arma::join_cols(thetas_lower, thetas_upper);
        return full_thetas;
    }
}

ModelParameters calculate_trial_params(
    int trial_idx,
    const ObjectiveData* objective_data,
    const arma::vec& beta,
    const std::map<std::string, double>& beta_map,
    const arma::vec& Thetas
) {
    ModelParameters p;

    for (const auto& par_name : ModelParameters::parnames) {
        double* target_param = nullptr;

        auto it = param_map.find(par_name);
        Param p_enum = (it != param_map.end()) ? it->second : unknown;

        switch (p_enum) {
            case a:      target_param = &p.a;      break;
            case v:      target_param = &p.v;      break;
            case t0:     target_param = &p.t0;     break;
            case d:      target_param = &p.d;      break;
            case z:      target_param = &p.z;      break;
            case sz:     target_param = &p.sz;     break;
            case sv:     target_param = &p.sv;     break;
            case st0:    target_param = &p.st0;    break;
            case tau:    target_param = &p.tau;    break;
            case lambda: target_param = &p.lambda; break;
            case w:      target_param = &p.w;      break;
            case muvis:  target_param = &p.muvis;  break;
            case sigvis: target_param = &p.sigvis; break;
            case svis:   target_param = &p.svis;   break;
            case s:      target_param = &p.s;      break;
            case unknown: continue;
        }
            
        if (objective_data->fixed.containsElementNamed(par_name.c_str())) {
            *target_param = Rcpp::as<double>(objective_data->fixed[par_name]);
        } else {
            double value = arma::datum::nan;
            if (objective_data->beta_map.containsElementNamed(par_name.c_str())) {
                Rcpp::IntegerVector beta_indices = objective_data->beta_map[par_name];
                arma::uvec u_beta_indices = Rcpp::as<arma::uvec>(beta_indices) - 1;
                arma::rowvec trial_row = objective_data->model_matrix.row(trial_idx);
                value = arma::dot(trial_row.elem(u_beta_indices), beta.elem(u_beta_indices));
            } else if (beta_map.count(par_name)) {
                value = beta_map.at(par_name);
            }
            *target_param = value;
        }
    }

    // fixed01
    p.z = R::pnorm(p.z, 0, 1, 1, 0);
    p.sz = R::pnorm(p.sz, 0, 1, 1, 0);
    p.w = R::pnorm(p.w, 0, 1, 1, 0);
    p.d = R::pnorm(p.d, 0, 1, 1, 0);
    p.t0 = R::pnorm(p.t0, 0, 1, 1, 0);
    p.st0 = R::pnorm(p.st0, 0, 1, 1, 0);

    // fixedpos
    p.a = exp(p.a);
    p.v = exp(p.v);
    p.sv = exp(p.sv);
    p.svis = exp(p.svis);
    p.sigvis = exp(p.sigvis);
    p.lambda = exp(p.lambda);
    p.s = exp(p.s);

    p.sz *= (2.0 * std::min(p.z, 1.0 - p.z));
    p.t0 *= objective_data->maxt0;
    p.d *= p.t0;

    if (std::isinf(objective_data->restr_tau)) {
        p.tau = exp(p.tau);
    } else if (objective_data->simult_conf) {
        p.tau *= (objective_data->maxt0 - p.t0);
    } else {
        p.tau *= objective_data->restr_tau;
    }

    if (objective_data->dependent_vars.n_cols > 3) {
        p.v *= objective_data->dependent_vars(trial_idx, 3);
    }

    int response = objective_data->dependent_vars(trial_idx, 1);
    int rating = objective_data->dependent_vars(trial_idx, 2);
    int idx = response * (objective_data->n_ratings + 2) + rating;
    p.th1 = Thetas(idx);
    p.th2 = Thetas(idx + 1);

    if (R_IsNaN(p.muvis)) p.muvis = std::abs(p.v);

    p.a /= p.s;
    p.v /= p.s;
    p.sv /= p.s;
    p.th1 /= p.s;
    p.th2 /= p.s;
    p.muvis /= p.s;
    p.sigvis /= p.s;
    p.svis /= p.s;
    p.st0 *= 2.0;
    p.t0 += p.st0 / 2.0;

    return p;
}

double neg_loglikelihood_formula(
    unsigned n,
    const double* x,
    double* grad,
    void* func_data
) {
    Rcpp::checkUserInterrupt();

    ObjectiveData* objective_data = static_cast<ObjectiveData*>(func_data);
    arma::vec beta(const_cast<double*>(x), n, false, true);

    int n_trials = objective_data->dependent_vars.n_rows;

    std::map<std::string, double> beta_map;
    for (unsigned int i = 0; i < n; ++i) {
        beta_map[std::string(objective_data->beta_names[i])] = beta[i];
    }

    arma::vec Thetas = calculate_thetas(beta, objective_data);

    double total_logl = 0.0;

    for (int i = 0; i < n_trials; ++i) {
        ModelParameters params = calculate_trial_params(i, objective_data, beta, beta_map, Thetas);

        bool boundary = objective_data->dependent_vars(i, 1);
        Rcpp::NumericVector fit_vector = params.to_vector(boundary, objective_data->precision);

        double rt = objective_data->dependent_vars(i, 0);
        double prob = fabs(g_minus_WEVmu(rt, fit_vector));

        if (prob == 0) {
            total_logl += std::log(std::numeric_limits<double>::min());
        } else if (R_IsNaN(prob) || prob < 0) {
            return 1e12;
        } else {
            total_logl += std::log(prob);
        }
    }

    if (!std::isfinite(total_logl)) return 1e12;

    return -total_logl;
}

ObjectiveData pack_objective_data(const Rcpp::List& context) {
    ObjectiveData objective_data;
    objective_data.dependent_vars = Rcpp::as<arma::mat>(context["dependent_vars"]);
    objective_data.model_matrix = Rcpp::as<arma::mat>(context["model_matrix"]);
    objective_data.beta_map = Rcpp::as<Rcpp::List>(context["beta_map"]);
    objective_data.fixed = Rcpp::as<Rcpp::List>(context["fixed"]);
    objective_data.beta_names = Rcpp::as<Rcpp::CharacterVector>(context["beta_names"]);
    objective_data.maxt0 = Rcpp::as<double>(context["maxt0"]);
    objective_data.restr_tau = Rcpp::as<double>(context["restr_tau"]);
    objective_data.precision = Rcpp::as<double>(context["precision"]);
    objective_data.n_ratings = Rcpp::as<int>(context["n_ratings"]);
    objective_data.simult_conf = Rcpp::as<bool>(context["simult_conf"]);
    objective_data.sym_thetas = Rcpp::as<bool>(context["sym_thetas"]);
    return objective_data;
}

// [[Rcpp::export]]
Rcpp::List nlopt_optimizer(const Rcpp::List& context, Rcpp::NumericVector start_params) {
    ObjectiveData objective_data = pack_objective_data(context);

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

    nlopt_set_min_objective(opt, neg_loglikelihood_formula, &objective_data);
    
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
    ObjectiveData objective_data = pack_objective_data(context);
    return neg_loglikelihood_formula(params.size(), params.begin(), NULL, &objective_data);
}