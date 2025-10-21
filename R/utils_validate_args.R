#' @keywords internal
validate_rtconf_models_args <- function(
  data,
  models,
  optim_method,
  fixed,
  manipulations,
  n_ratings,
  restr_tau,
  sym_thetas,
  precision,
  opts,
  grid_search,
  logging,
  parallel,
  n_cores
) {
  # models
  assert_that(
    is.character(models),
    all(vapply(models, is.string, logical(1))),
    msg = sprintf(
      "`models` must be a character vector\nGot: %s",
      describe(models)
    )
  )
  assert_that(
    all(models %in% MODELS),
    msg = sprintf(
      "`models` must all be one of the supported models:\n%s\nGot: %s",
      paste(MODELS, collapse = ", "),
      paste(models, collapse = ", ")
    )
  )
  # manipulations
  assert_that(
    is.list(manipulations),
    msg = sprintf(
      "`manipulations` must be a list\nGot: %s",
      describe(manipulations)
    )
  )
  # either a single list of formulas (each element is a formula), or
  # a list of lists where each inner list is a list of formulas matching models length
  if (length(manipulations) > 0) {
    if (inherits(manipulations[[1]], "formula")) {
      # single list of formulas, check each element
      assert_that(
        all(vapply(manipulations, inherits, what = "formula", FUN.VALUE = logical(1))),
        all(vapply(manipulations, function(x) length(all.vars(x[[3]])) > 0, FUN.VALUE = logical(1))),
        msg = sprintf(
          "`manipulations` must be a list of formulas of the form LHS ~ RHS\nGot: %s",
          describe(manipulations)
        )
      )
    } else {
      # manipulations is a list of lists
      assert_that(
        length(manipulations) == length(models),
        msg = sprintf(
          "`manipulations` must be a list of lists matching length of models\nGot length: %d, Expected length: %d",
          length(manipulations), length(models)
        )
      )
      for (m in seq_along(manipulations)) {
        inner_manip <- manipulations[[m]]
        assert_that(
          is.list(inner_manip),
          msg = sprintf(
            "Each element of manipulations must be a list, but manipulations[[%d]] is not\nGot: %s",
            m,
            describe(inner_manip)
          )
        )
        if (length(inner_manip) > 0) {
          assert_that(
            all(vapply(inner_manip, inherits, "formula", logical(1))),
            all(vapply(inner_manip, function(x) length(all.vars(x[[3]])) > 0, logical(1))),
            msg = sprintf(
              "All manipulations in manipulations[[%d]] must be formulas of the form LHS ~ RHS\nGot: %s",
              m,
              describe(inner_manip)
            )
          )
        }
      }
    }
  }

  validate_rtconf_common_args(
    data,
    optim_method,
    fixed,
    n_ratings,
    restr_tau,
    sym_thetas,
    precision,
    opts,
    grid_search,
    logging,
    parallel,
    n_cores
  )
}

#' @keywords internal
validate_rtconf_args <- function(
  data,
  model,
  optim_method,
  fixed,
  manipulations,
  n_ratings,
  restr_tau,
  sym_thetas,
  precision,
  opts,
  grid_search,
  logging,
  parallel,
  n_cores
) {
  # model
  assert_that(
    is.string(model),
    msg = sprintf(
      "`model` must be a string\nGot: %s",
      describe(model)
    )
  )
  assert_that(
    model %in% MODELS,
    msg = sprintf(
      "`model` must be one of the supported models:\n%s\nGot: %s",
      paste(MODELS, collapse = ", "),
      model
    )
  )
  # manipulations
  assert_that(
    is.list(manipulations),
    all(vapply(manipulations, inherits, what = "formula", FUN.VALUE = logical(1))),
    all(vapply(manipulations, function(x) length(all.vars(x[[3]])) > 0, FUN.VALUE = logical(1))),
    msg = sprintf(
      "`manipulations` must be a list of formulas of the form LHS ~ RHS\nGot: %s",
      describe(manipulations)
    )
  )

  validate_rtconf_common_args(
    data,
    optim_method,
    fixed,
    n_ratings,
    restr_tau,
    sym_thetas,
    precision,
    opts,
    grid_search,
    logging,
    parallel,
    n_cores
  )
}

#' @keywords internal
validate_rtconf_common_args <- function(
  data,
  optim_method,
  fixed,
  n_ratings,
  restr_tau,
  sym_thetas,
  precision,
  opts,
  grid_search,
  logging,
  parallel,
  n_cores
) {
  # data
  assert_that(
    is.data.frame(data),
    msg = sprintf(
      "`data` must be a data frame\nGot: %s",
      describe(data)
    )
  )
  # optim_method
  assert_that(
    is.string(optim_method),
    msg = sprintf(
      "`optim_method` must be a string\nGot: %s",
      describe(optim_method)
    )
  )
  assert_that(
    optim_method %in% OPTIM_METHODS,
    msg = sprintf(
      "`optim_method` must be one of the supported optim_methods:\n%s\nGot: %s",
      paste(OPTIM_METHODS, collapse = ", "),
      optim_method
    )
  )
  # fixed
  assert_that(
    is.list(fixed),
    !is.null(names(fixed)),
    all(names(fixed) != ""),
    all(vapply(fixed, function(x) is.numeric(x) && length(x) == 1 && !is.na(x), logical(1))),
    msg = sprintf(
      "`fixed` must be a named list\nGot: %s",
      describe(fixed)
    )
  )
  # n_ratings
  assert_that(
    is.null(n_ratings) || (
      is.numeric(n_ratings) &&
      length(n_ratings) == 1 &&
      n_ratings %% 1 == 0 &&
      n_ratings >= 2
    ),
    msg = sprintf(
      "`n_ratings` must be NULL or an integer >= 2\nGot: %s",
      describe(n_ratings)
    )
  )
  # restr_tau
  assert_that(
    length(restr_tau) == 1,
    (is.numeric(restr_tau) && restr_tau > 0) || (is.string(restr_tau) && restr_tau == "simult_conf"),
    msg = sprintf(
      "`restr_tau` must be a positive numeric scalar (or Inf), or the string 'simult_conf'\nGot: %s",
      describe(restr_tau)
    )
  )
  # sym_thetas
  assert_that(
    is.logical(sym_thetas),
    length(sym_thetas) == 1,
    msg = sprintf(
      "`sym_thetas` must be a single logical value\nGot: %s",
      describe(sym_thetas)
    )
  )
  # precision
  assert_that(
    is.numeric(precision),
    length(precision) == 1,
    precision %% 1 == 0,
    precision >= 1,
    msg = sprintf(
      "`precision` must be an integer >= 1\nGot: %s",
      describe(precision)
    )
  )
  # opts
  assert_that(
    is.list(opts),
    length(opts) == 0 || (!is.null(names(opts)) && all(names(opts) != "")),
    msg = sprintf(
      "`opts` must be an empty list or a non-empty named list\nGot: %s",
      describe(opts)
    )
  )
  if (length(opts) > 0) {
    invalid_opts <- setdiff(names(opts), names(DEFAULT_OPTS))
    assert_that(
      length(invalid_opts) == 0,
      msg = sprintf(
        "`opts` must be a list of supported options settings:\n%s\nGot: %s",
        paste(names(DEFAULT_OPTS), collapse = ", "),
        describe(opts)
      )
    )
    if (!is.null(opts$n_attempts)) {
      assert_that(
        is.numeric(opts$n_attempts),
        length(opts$n_attempts) == 1,
        opts$n_attempts %% 1 == 0,
        opts$n_attempts >= 1,
        msg = sprintf(
          "`opts$n_attempts` must be an integer >= 1\nGot: %s",
          describe(opts$n_attempts)
        )
      )
    }
    if (!is.null(opts$n_restarts)) {
      assert_that(
        is.numeric(opts$n_restarts),
        length(opts$n_restarts) == 1,
        opts$n_restarts %% 1 == 0,
        opts$n_restarts >= 1,
        msg = sprintf(
          "`opts$n_restarts` must be an integer >= 1\nGot: %s",
          describe(opts$n_restarts)
        )
      )
    }
    if (!is.null(opts$maxfun)) {
      assert_that(
        is.numeric(opts$maxfun),
        length(opts$maxfun) == 1,
        opts$maxfun %% 1 == 0,
        opts$maxfun >= 1,
        msg = sprintf(
          "`opts$maxfun` must be an integer >= 1\nGot: %s",
          describe(opts$maxfun)
        )
      )
    }
    if (!is.null(opts$reltol)) {
      assert_that(
        is.numeric(opts$reltol),
        length(opts$reltol) == 1,
        opts$reltol > 0,
        msg = sprintf(
          "`opts$reltol` must be a single positive numeric value\nGot: %s",
          describe(opts$reltol)
        )
      )
    }
  }
  # grid_search
  assert_that(
    is.logical(grid_search),
    length(grid_search) == 1,
    msg = sprintf(
      "`grid_search` must be a single logical value\nGot: %s",
      describe(grid_search)
    )
  )
  # logging
  assert_that(
    is.logical(logging),
    length(logging) == 1,
    msg = sprintf(
      "`logging` must be a single logical value\nGot: %s",
      describe(logging)
    )
  )
  # parallel
  assert_that(
    is.logical(parallel),
    length(parallel) == 1,
    msg = sprintf(
      "`parallel` must be a single logical value\nGot: %s",
      describe(parallel)
    )
  )
  # n_cores
  assert_that(
    is.null(n_cores) || (
      is.numeric(n_cores) &&
      length(n_cores) == 1 &&
      n_cores %% 1 == 0 &&
      n_cores >= 1
    ),
    msg = sprintf(
      "`n_cores` must be NULL or an integer >= 1\nGot: %s",
      describe(n_cores)
    )
  )
}

#' @keywords internal
validate_input_data <- function(
  cols,
  rt,
  rating,
  response,
  stimulus,
  correct
) {
  # cols
  assert_that(
    "rt" %in% cols,
    msg = sprintf(
      "`data` must contain a `rt` column\nGot: %s",
      describe(cols)
    )
  )
  assert_that(
    "rating" %in% cols,
    msg = sprintf(
      "`data` must contain a `rating column\nGot: %s",
      describe(cols)
    )
  )
  assert_that(
    ("response" %in% cols) || ("correct" %in% cols),
    msg = sprintf(
      "`data` must contain at least `response` or `correct`\nGot: %s",
      describe(cols)
    )
  )
  # rt
  assert_that(
    is.numeric(rt),
    all(rt > 0),
    msg = sprintf(
      "`rt` must be a positive numeric vector\nGot: %s",
      describe(rt)
    )
  )
  # rating
  assert_that(
    is.numeric(rating),
    all(rating >= 0),
    all(rating %% 1 == 0),
    length(unique(rating)) >= 2,
    msg = sprintf(
      "`rating` must be a numeric integer vector >= 0 with at least 2 unique levels\nGot: %s",
      describe(rating)
    )
  )
  # response
  assert_that(
    is.null(response) || length(unique(response)) == 2,
    msg = sprintf(
      "`response` must have exactly 2 unique values\nGot: %s",
      describe(response)
    )
  )
  # stimulus
  assert_that(
    is.null(stimulus) || length(unique(stimulus)) == 2,
    msg = sprintf(
      "`stimulus` must have exactly 2 unique values\nGot: %s",
      describe(stimulus)
    )
  )
  # correct
  assert_that(
    is.null(correct) || all(correct %in% c(0, 1)),
    msg = sprintf(
      "`correct` must contain only 0 and 1 values\nGot: %s",
      describe(correct)
    )
  )
}

#' @keywords internal
describe <- function(x, max_lines = 5) {
  out <- capture.output(str(x))
  if (length(out) > max_lines) {
    out <- c(out[1:max_lines], sprintf("... [%d more lines]", length(out) - max_lines))
  }
  paste(out, collapse = "\n")
}
