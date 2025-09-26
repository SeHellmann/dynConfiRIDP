# Find the include paths for the required R packages
rcpp_path <- system.file("include", package = "Rcpp")
rcpp_armadillo_path <- system.file("include", package = "RcppArmadillo")
nloptr_path <- system.file("include", package = "nloptr")

r_path <- R.home("include")

# Combine all paths for the compiler flags
include_paths <- c(
  paste0("-I", rcpp_path),
  paste0("-I", rcpp_armadillo_path),
  paste0("-I", nloptr_path),
  paste0("-I", r_path),
  "-std=c++17"
)

add_flags <- paste0("\t\t\"", include_paths, "\"", collapse = ",\n")

clangd_content <- paste0(
  "CompileFlags:\n",
  "\tAdd: [\n",
  add_flags,
  "\n\t]"
)

writeLines(clangd_content, here::here(".clangd"))

message("Successfully created .clangd configuration file in project root.")
