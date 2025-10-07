#pragma once

#include <Rcpp.h>
#include <string>

namespace Logger {
    static bool logging = false;

    inline void init(bool enabled) {
        logging = enabled;
    }

    inline void log(const std::string& level, const std::string& msg) {
        if (!logging) return;

        try {
            Rcpp::Function log_fn = Rcpp::Environment::namespace_env("logger")[std::string("log_") + level];
            log_fn(msg);
        } catch (...) {
            // Fallback to Rcout if logger fails
            Rcpp::Rcout << "[" << level << "] " << msg << std::endl;
        }
    }

    inline void info(const std::string& msg) { log("info", msg); }
    inline void warn(const std::string& msg) { log("warn", msg); }
    inline void error(const std::string& msg) { log("error", msg); }
}