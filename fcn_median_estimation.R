###################################
# Deductive estimation of median
# Source code of functions
#      Tianchen Qian
#        2014.9.18
###################################

deduct_est_median <- function(x_df, r_vec, y_vec,
                            e_vec, py_vec_fcn,
                            beta_init, free_idx = 1, 
                            eps = 0.00001, p_vec = NULL, jack_se = FALSE){
  ##########
  # Return deductive mean estimator and its jackknife standard error.
  # The mean is for everyone being assigned R = 1.
  # Estimated by freeing-up the intercept in beta_init, and solve for
  # delta that makes sum of EIF zero.
  #
  # Parameters:
  # -----------
  # x_df: matrix or data frame of nrow = n (n = total number of subjects, don't need to specify it)
  #     Covariate matrix, each row being an observation
  # r_vec: binary vector of length n
  #     Vector of assignment, 1 for treatment, 0 for control
  # y_vec: vector of length n
  #     Vector of observed outcome
  #     (It doesn't matter what you put for r_vec = 0; put 0 or NA if you like)
  # p_vec: vector of length n, or default = NULL
  #     Marginal distribution of X of the working model
  #     Put discrete mass on each row of x_df
  #     If NULL, put 1/n on each point
  # e_vec: vector of length n
  #     Propensity score of the working model
  #     Each entry corresponds to the row in x_df
  # py_vec_fcn: function that returns vector value
  #     P{Y<=t | X, beta} of the working model
  #     Takes in arguments as (t, x_df, beta_init) and returns a vector value of P{Y<=t | x_df, beta_init}
  # beta_init: vector
  #     Coefficients in expected outcome model, with first entry being intercept
  # free_idx: number, default = 1; or a character indicating the name of the variable in x_df whose coefficient is to be free-up
  #     Which coefficient in beta to free up? Default is intercept
  # eps: number, default = 0.00001
  #     The small epsilon used in calculating numerical Gateaux derivative,
  # jack.se: logical, detault = FALSE
  #     The function returns also jackknife standard error if jack.se is TRUE
  #     This may be very time consuming
  ##########
  
  # assign default p_vec = 1/n
  if (is.null(p_vec)){
    n <- nrow(x_df)
    p_vec <- rep(1/n, n)
  }
  
  # check agreement of length
  len <- c(nrow(x_df), length(r_vec), length(y_vec), length(p_vec), length(e_vec))
  if ( !all(len[1] == len) ){
    stop("The lengths of r_vec, y_vec, p_vec, e_vec and nrow(x_df) don't agree!")
  }
  
  # check no NA's in beta_init (otherwise delete corresponding columns in x_df)
  # only do the deletion when dealing with linear term x*beta, with beta including an intercept
  if (length(beta_init) - ncol(x_df) == 1){
    if ( any(is.na(beta_init)) ){
      warning("There are NA's in beta_init. Deleted them with corresponding columns in x_df.\n  Assuming linear term x*beta in y_vec_fcn.\n",
              immediate. = TRUE)
      x_df <- x_df[, -(which(is.na(beta_init))-1)]
      beta_init <- beta_init[-which(is.na(beta_init))]
    }
  } else if ( any(is.na(beta_init)) ){
    stop("There are NA's in beta_init.")
  }
  
  if (class(x_df) == "data.frame") {
    x_df <- data.matrix(x_df)
  }
  
  if (class(free_idx) == "numeric"){
    message(paste0("The coefficient to be free-up is ", names(beta_init)[free_idx]))
  } else if (class(free_idx) == "character"){
    free_idx <- which(names(beta_init) == free_idx)
    message(paste0("The coefficient to be free-up is ", names(beta_init)[free_idx]))
  } else {
    stop("free_idx should be either numeric or character!")
  }

  y_vec[r_vec == 0] <- 0
  
  ### Compute deductive mean estimate ###
  F_beta <- list(p_vec = p_vec, e_vec = e_vec, y_vec_fcn = py_vec_fcn)
  
  # Searching for delta to make sum(EIF) zero
  message("\nSearching for delta to make sum(EIF) zero...")

  delta <- tryCatch({
    uniroot(function(delta){
      beta <- beta_init
      beta[free_idx] <- beta[free_idx] + delta
      return(sum_EIF(beta, F_beta, median_, median_Di_eps, x_df, r_vec, y_vec, eps))
    }, interval = c(-10,10), extendInt = "yes", maxiter = 200)
  },
  error = function(cond){
    message("\nError in solving delta.")
    message("Here's the original error message:")
    message(cond)
    return("error_delta")
  }
  )
  
  if (identical(delta, "error_delta")){ # if an error occurred in solving for delta, return NA
    if (jack_se){
      output <- list(estimate = NA, std_err = NA)
    } else {
      output <- NA
    }    
    return(output)
  }
  
  message(sprintf("    delta: %f; sum(EIF): %f; iter: %d",
                  delta$root, delta$f.root, delta$iter))  
  beta <- beta_init
  beta[free_idx] <- beta[free_idx] + delta$root
  estimate <- median_(beta, F_beta, x_df)
  message(sprintf("    Deductive median estimate: %f\n", estimate))
  
  if (jack_se){
    ### Compute jackknife standard error of the deductive mean estimate ###  
    fcn_jack <- function(x_df, r_vec, y_vec,
                         p_vec, e_vec, y_vec_fcn){
      F_beta <- list(p_vec = p_vec, e_vec = e_vec, y_vec_fcn = y_vec_fcn)    
      delta <- uniroot(function(delta){
        beta <- beta_init
        beta[free_idx] <- beta[free_idx] + delta
        return(sum_EIF(beta, F_beta, median_, median_Di_eps, x_df, r_vec, y_vec, eps))
      }, interval = c(-10,10), extendInt = "yes")
      beta <- beta_init
      beta[free_idx] <- beta[free_idx] + delta$root
      estimate <- median_(beta, F_beta, x_df)
    }    
    message("Calculating Jackknife standard error of the deductive estimator...\n    This may take a while...")
    std_err <- jackknife(1:nrow(x_df),
                             function(rows, x_df, r_vec, y_vec, p_vec, e_vec, y_vec_fcn){
                               fcn_jack(x_df[rows,], r_vec[rows], y_vec[rows], p_vec[rows], e_vec[rows], y_vec_fcn)
                             },
                             x_df = x_df, r_vec = r_vec, y_vec = y_vec,
                             p_vec = p_vec, e_vec = e_vec, y_vec_fcn = py_vec_fcn)
    message(sprintf("    Jackknife standard error: %f\n", std_err$jack.se))
  }
  
  # Construct output list
  if (jack_se){
    output <- list(estimate = estimate, std_err = std_err$jack.se)
  } else {
    output <- estimate
  }
  
  invisible(output)
}


median_ <- function(beta, F_beta, x_df){
  ##########
  # Return a vector p_Di_eps_vec and a function y_Di_eps_fcn(t, x, beta)
  #
  # Keyword Arguments:
  # beta   -- vector with length p+1, first element being intercept
  # F_beta -- list of three elements, as working distribution
  # x_df   -- data frame, each row being covariates for an observation
  ##########
  
  # unpack parameters
  p_vec <- F_beta$p_vec
  e_vec <- F_beta$e_vec
  y_vec_fcn <- F_beta$y_vec_fcn
  
  # uniroot solving for median
  objective_f <- function(t){
    sum(y_vec_fcn(t, x_df, beta)*p_vec) - 0.5
  }
  
  answer <- uniroot(objective_f, interval = c(-10, 10), extendInt = "yes")
  
  return(answer$root)
}

median_Di_eps <- function(beta, F_beta, x_df, D_i, i, eps){
  ##########
  # Return a vector p_Di_eps_vec and a function y_Di_eps_fcn(t, x, beta)
  #
  # Keyword Arguments:
  # beta   -- vector with length p+1, first element being intercept
  # F_beta -- list of three elements, as working distribution
  # x_df   -- data frame, each row being covariates for an observation
  # D_i    -- list consisting (X_i, R_i, Y_i = Y_i*R_i) as a perturbing data point
  # i      -- index of D_i
  # eps    -- a small number
  ##########
  
  # unpack parameters
  p_vec <- F_beta$p_vec
  e_vec <- F_beta$e_vec
  y_vec_fcn <- F_beta$y_vec_fcn
  X_i <- D_i$X_i
  R_i <- D_i$R_i
  Y_i <- D_i$Y_i
  
  # compute p_Di_eps_vec
  p_Di_eps_vec <- (1-eps)*p_vec
  p_Di_eps_vec[i] <- p_Di_eps_vec[i] + eps
  
  # compute y_Di_eps_vec_fcn
  y_Di_eps_vec_fcn <- function(t){
    
    # j != i
    tmp <- (1-eps)*p_vec*e_vec
    return_vec <- ( tmp*y_vec_fcn(t, x_df, beta) ) / tmp
    if (length(return_vec) != nrow(x_df)) {
      stop("Error in y_Di_eps_vec_fcn()! length(return_vec) != nrow(x_df)")
    }
    
    # j == i
    tmp <- (1-eps)*p_vec[i]*e_vec[i]
    return_vec[i] <- (eps*R_i*(Y_i <= t) + tmp*y_vec_fcn(t, X_i, beta)) / (eps*R_i + tmp)
    
    return(return_vec)
  }
  
  # try using uniroot
  objective_f <- function(t){
    sum(y_Di_eps_vec_fcn(t)*p_Di_eps_vec) - 0.5
  }
  
  answer <- uniroot(objective_f, interval = c(-10, 10), extendInt = "yes")
  
  return(answer$root)
}

sum_EIF <- function(beta, F_beta, tau_fcn, tau_fcn_eps, x_df, r_vec, y_vec, eps){
  ##########
  # Return the sum of EIFs (sum over i = 1:n)
  #
  # Keyword Arguments:
  # beta        -- vector with length p+1, first element being intercept
  # tau_fcn     -- function to compute estimator tau under working model
  # tau_fcn_eps -- function to compute estimator tau under perturbed working model
  # F_beta      -- list of three elements, as working distribution
  # x_df        -- data frame, each row being covariates for an observation
  # r_vec       -- vector of R_i
  # y_vec  -- vector of continuous Y_i
  # eps         -- small number
  ##########
  
  phi <- rep(NA, nrow(x_df))
  for(i in 1:nrow(x_df)){
    D_i <- list(X_i = x_df[i, ], R_i = r_vec[i], Y_i = y_vec[i])
    tmp1 <- tau_fcn_eps(beta, F_beta, x_df, D_i, i, eps)
    tmp2 <- tau_fcn(beta, F_beta, x_df)
    phi[i] <- (tmp1 - tmp2) / eps
  }
  
  return(sum(phi))
}
