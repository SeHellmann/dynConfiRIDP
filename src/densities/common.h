#pragma once
#include <cmath>

namespace constants {
constexpr double EPSILON = 1e-6;
constexpr double fac_errf = 2.170803763674803; // sqrt(M_PI*3/2); used in densPCRM
constexpr float inv_sqrt_2pi = 0.3989422804014327; // used in density_IRM
} // namespace constants

namespace utils {
inline double Phi(double x) { return 0.5 * (1 + erf(x / M_SQRT2)); }
} // namespace utils