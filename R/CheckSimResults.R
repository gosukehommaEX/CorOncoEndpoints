#' Check Simulation Results Against Theoretical Values
#'
#' @description
#' Compares empirical simulation results from \code{\link{rOncoEndpoints}} with
#' their corresponding theoretical values. This function calculates standard
#' performance metrics (Bias, Relative Bias, SE, MSE, RMSE) for validating
#' simulation implementations and assessing the accuracy of the proposed
#' random number generation method.
#'
#' @param dataset A data frame. The output from \code{\link{rOncoEndpoints}}
#'   containing simulated endpoint data with columns simID, Group, and endpoint
#'   variables (OS, PFS, and/or Response).
#' @param p Named numeric vector or NULL. The true probabilities of the binary
#'   response endpoint for each group. Names must match group names in dataset.
#'   Required if Response endpoint is present. Example: c(Treatment = 0.4, Control = 0.3)
#' @param hazard_OS Named numeric vector or NULL. The hazard rates for OS for
#'   each group. Names must match group names in dataset. Required if OS
#'   endpoint is present. Example: c(Treatment = 0.05, Control = 0.07)
#' @param hazard_PFS Named numeric vector or NULL. The hazard rates for PFS for
#'   each group. Names must match group names in dataset. Required if PFS
#'   endpoint is present. Example: c(Treatment = 0.08, Control = 0.10)
#' @param rho_tte_resp Named numeric vector or NULL. The specified correlations
#'   between TTE and Response for each group. Names must match group names in
#'   dataset. Required if both TTE and Response endpoints are present.
#'   Example: c(Treatment = 0.3, Control = 0.2)
#' @param copula Character or NULL. The copula family used for modeling dependence.
#'   Options are "Clayton" or "Frank". Required if both TTE and Response
#'   endpoints are present.
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{Group}{Treatment group name}
#'   \item{Endpoint}{Name of the endpoint or statistic}
#'   \item{Empirical}{Mean of estimates across simulations}
#'   \item{Theoretical}{Expected theoretical value}
#'   \item{Bias}{Bias = Empirical - Theoretical (signed difference showing
#'     direction of systematic error)}
#'   \item{Relative_Bias}{Relative bias as percentage: 100 × Bias / Theoretical.
#'     Positive values indicate overestimation, negative values indicate
#'     underestimation}
#'   \item{SE}{Empirical standard error (SD of estimates across simulations)}
#'   \item{MSE}{Mean squared error = Bias² + SE²}
#'   \item{RMSE}{Root mean squared error = √MSE (overall accuracy in original scale)}
#'   \item{Assessment}{Quick interpretation: "Excellent" (|Relative_Bias| < 5%),
#'     "Acceptable" (5% ≤ |Relative_Bias| < 10%), "Review" (|Relative_Bias| ≥ 10%)}
#' }
#'
#' @details
#' This function calculates both empirical and theoretical values for:
#'
#' **Time-to-event endpoints (OS/PFS):**
#' \itemize{
#'   \item Mean: Empirical vs 1/hazard
#'   \item Median: Empirical vs log(2)/hazard
#' }
#'
#' **Binary endpoint (Response):**
#' \itemize{
#'   \item Proportion: Empirical vs p
#' }
#'
#' **Correlations:**
#' \itemize{
#'   \item Corr(OS, Response): Empirical vs rho_tte_resp (specified)
#'   \item Corr(PFS, Response): Empirical vs value calculated by \code{\link{CorResponsePFS}}
#'   \item Corr(OS, PFS): Empirical vs hazard_OS / hazard_PFS (Fleischer model)
#' }
#'
#' **Performance metrics:**
#' \itemize{
#'   \item **Bias**: Measures systematic error. Should be close to 0 for unbiased
#'     methods. Sign indicates direction: positive = overestimation,
#'     negative = underestimation
#'   \item **Relative Bias (%)**: Bias relative to true value. Recommended
#'     interpretation: <5% excellent, <10% acceptable
#'   \item **SE (Standard Error)**: Measures variability of estimates. Decreases
#'     with √nsim. Smaller is better
#'   \item **MSE**: Combines bias and variance into single metric. Smaller is better
#'   \item **RMSE**: MSE in original scale. Directly comparable to SE
#'   \item **Assessment**: Automatic interpretation based on relative bias
#' }
#'
#' The function is particularly useful for:
#' \itemize{
#'   \item Validating that \code{\link{rOncoEndpoints}} correctly implements the models
#'   \item Demonstrating unbiasedness of the proposed method for publication
#'   \item Checking if sufficient simulations have been run (SE should be small)
#'   \item Understanding the relationship between parameters and resulting correlations
#'   \item Quality control in simulation studies
#' }
#'
#' @note
#' **Interpretation guidelines:**
#' \itemize{
#'   \item **Bias close to 0**: Method is unbiased (desirable for publication)
#'   \item **Relative Bias < 5%**: Excellent performance
#'   \item **Relative Bias < 10%**: Acceptable performance
#'   \item **Small SE**: Stable estimates (increase nsim to reduce SE)
#'   \item **RMSE ≈ SE when Bias ≈ 0**: Indicates unbiased estimator
#'   \item **Correlations typically show higher variability** (larger SE/RMSE)
#'     than means/medians, especially with smaller sample sizes
#' }
#'
#' **For publication:**
#' \itemize{
#'   \item Report Bias and Relative Bias to demonstrate unbiasedness
#'   \item Report SE to show precision
#'   \item Report MSE or RMSE for overall accuracy
#'   \item Typical table format: Group | Endpoint | Theoretical | Empirical |
#'     Bias | Relative Bias (%) | SE | RMSE | Assessment
#' }
#'
#' @references
#' Fleischer, F., Gaschler-Markefski, B., & Bluhmki, E. (2009). A statistical
#' model for the dependence between progression-free survival and overall
#' survival. Statistics in Medicine, 28(21), 2669-2686.
#'
#' @importFrom stats cor median sd
#'
#' @examples
#' # Example 1: OS and Response with Clayton copula
#' set.seed(123)
#' sim_data1 <- rOncoEndpoints(
#'   nsim = 1000,
#'   group = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   p = c(0.4, 0.3),
#'   hazard_OS = c(0.05, 0.07),
#'   rho_tte_resp = c(0.3, 0.2),
#'   copula = "Clayton"
#' )
#'
#' check1 <- CheckSimResults(
#'   dataset = sim_data1,
#'   p = c(Treatment = 0.4, Control = 0.3),
#'   hazard_OS = c(Treatment = 0.05, Control = 0.07),
#'   rho_tte_resp = c(Treatment = 0.3, Control = 0.2),
#'   copula = "Clayton"
#' )
#' print(check1, n = Inf)
#'
#' # Interpretation:
#' # - Bias close to 0: method is unbiased
#' # - Relative_Bias < 5%: excellent performance
#' # - Small SE: precise estimates
#' # - RMSE ≈ SE when Bias ≈ 0: confirms unbiasedness
#'
#' # Example 2: All three endpoints (OS, PFS, Response) with Frank copula
#' set.seed(456)
#' sim_data2 <- rOncoEndpoints(
#'   nsim = 1000,
#'   group = c("Experimental", "Standard"),
#'   n = c(150, 150),
#'   p = c(0.5, 0.35),
#'   hazard_OS = c(0.04, 0.06),
#'   hazard_PFS = c(0.08, 0.10),
#'   rho_tte_resp = c(0.4, 0.25),
#'   copula = "Frank"
#' )
#'
#' check2 <- CheckSimResults(
#'   dataset = sim_data2,
#'   p = c(Experimental = 0.5, Standard = 0.35),
#'   hazard_OS = c(Experimental = 0.04, Standard = 0.06),
#'   hazard_PFS = c(Experimental = 0.08, Standard = 0.10),
#'   rho_tte_resp = c(Experimental = 0.4, Standard = 0.25),
#'   copula = "Frank"
#' )
#' print(check2, n = Inf)
#'
#' # Note: PFS-Response correlation theoretical values are calculated
#' # using CorResponsePFS function, demonstrating the consistency
#' # of the three-endpoint model
#'
#' @seealso
#' \code{\link{rOncoEndpoints}} for generating correlated oncology endpoints,
#' \code{\link{CorResponsePFS}} for calculating PFS-Response correlation,
#' \code{\link{CorBoundResponseTTE}} for correlation bounds,
#' \code{\link{CopulaParamResponseTTE}} for copula parameters
#'
#' @export
CheckSimResults <- function(dataset, p = NULL, hazard_OS = NULL, 
                            hazard_PFS = NULL, rho_tte_resp = NULL, 
                            copula = NULL) {
  
  # Input validation
  if (!is.data.frame(dataset)) {
    stop("'dataset' must be a data frame")
  }
  
  required_cols <- c("simID", "Group")
  if (!all(required_cols %in% names(dataset))) {
    stop("'dataset' must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Determine which endpoints are present
  has_os <- "OS" %in% names(dataset)
  has_pfs <- "PFS" %in% names(dataset)
  has_response <- "Response" %in% names(dataset)
  
  if (!any(has_os, has_pfs, has_response)) {
    stop("'dataset' must contain at least one endpoint column (OS, PFS, or Response)")
  }
  
  # Validate parameters match endpoints
  if (has_os && is.null(hazard_OS)) {
    stop("'hazard_OS' must be provided when OS endpoint is present")
  }
  if (has_pfs && is.null(hazard_PFS)) {
    stop("'hazard_PFS' must be provided when PFS endpoint is present")
  }
  if (has_response && is.null(p)) {
    stop("'p' must be provided when Response endpoint is present")
  }
  
  # Check if TTE and Response are both present
  has_tte <- has_os || has_pfs
  if (has_response && has_tte) {
    if (is.null(rho_tte_resp)) {
      stop("'rho_tte_resp' must be provided when both TTE and Response endpoints are present")
    }
    if (is.null(copula)) {
      stop("'copula' must be provided when both TTE and Response endpoints are present")
    }
  }
  
  # Get unique groups and preserve order
  group_order <- unique(dataset$Group)
  
  # Validate that parameter vectors have names matching groups
  validate_named_vector <- function(vec, name) {
    if (!is.null(vec)) {
      if (is.null(names(vec))) {
        stop("'", name, "' must be a named vector with names matching group names.\n",
             "Example: c(", paste(group_order, "= value", collapse = ", "), ")")
      }
      missing_groups <- setdiff(group_order, names(vec))
      if (length(missing_groups) > 0) {
        stop("'", name, "' is missing values for groups: ", 
             paste(missing_groups, collapse = ", "), "\n",
             "Required groups: ", paste(group_order, collapse = ", "))
      }
      extra_groups <- setdiff(names(vec), group_order)
      if (length(extra_groups) > 0) {
        warning("'", name, "' contains values for non-existent groups: ",
                paste(extra_groups, collapse = ", "), ". These will be ignored.")
      }
    }
  }
  
  validate_named_vector(p, "p")
  validate_named_vector(hazard_OS, "hazard_OS")
  validate_named_vector(hazard_PFS, "hazard_PFS")
  validate_named_vector(rho_tte_resp, "rho_tte_resp")
  
  # Process each group
  result_list <- lapply(group_order, function(grp) {
    grp_data <- dataset[dataset$Group == grp, ]
    grp_data_by_sim <- split(grp_data, grp_data$simID)
    nsim <- length(grp_data_by_sim)
    
    # Initialize results
    empirical_values <- list()
    theoretical_values <- list()
    endpoint_names <- character(0)
    
    # Calculate endpoint statistics
    if (has_os) {
      # OS Mean
      os_mean_sims <- sapply(grp_data_by_sim, function(x) mean(x$OS))
      os_mean_emp <- mean(os_mean_sims)
      os_mean_se <- sd(os_mean_sims)
      os_mean_theo <- 1 / hazard_OS[grp]
      empirical_values <- c(empirical_values, list(OS_Mean = c(value = os_mean_emp, se = os_mean_se)))
      theoretical_values <- c(theoretical_values, OS_Mean = os_mean_theo)
      endpoint_names <- c(endpoint_names, "OS_Mean")
      
      # OS Median
      os_median_sims <- sapply(grp_data_by_sim, function(x) median(x$OS))
      os_median_emp <- mean(os_median_sims)
      os_median_se <- sd(os_median_sims)
      os_median_theo <- log(2) / hazard_OS[grp]
      empirical_values <- c(empirical_values, list(OS_Median = c(value = os_median_emp, se = os_median_se)))
      theoretical_values <- c(theoretical_values, OS_Median = os_median_theo)
      endpoint_names <- c(endpoint_names, "OS_Median")
    }
    
    if (has_pfs) {
      # PFS Mean
      pfs_mean_sims <- sapply(grp_data_by_sim, function(x) mean(x$PFS))
      pfs_mean_emp <- mean(pfs_mean_sims)
      pfs_mean_se <- sd(pfs_mean_sims)
      pfs_mean_theo <- 1 / hazard_PFS[grp]
      empirical_values <- c(empirical_values, list(PFS_Mean = c(value = pfs_mean_emp, se = pfs_mean_se)))
      theoretical_values <- c(theoretical_values, PFS_Mean = pfs_mean_theo)
      endpoint_names <- c(endpoint_names, "PFS_Mean")
      
      # PFS Median
      pfs_median_sims <- sapply(grp_data_by_sim, function(x) median(x$PFS))
      pfs_median_emp <- mean(pfs_median_sims)
      pfs_median_se <- sd(pfs_median_sims)
      pfs_median_theo <- log(2) / hazard_PFS[grp]
      empirical_values <- c(empirical_values, list(PFS_Median = c(value = pfs_median_emp, se = pfs_median_se)))
      theoretical_values <- c(theoretical_values, PFS_Median = pfs_median_theo)
      endpoint_names <- c(endpoint_names, "PFS_Median")
    }
    
    if (has_response) {
      # Response Proportion
      response_sims <- sapply(grp_data_by_sim, function(x) mean(x$Response))
      response_emp <- mean(response_sims)
      response_se <- sd(response_sims)
      response_theo <- p[grp]
      empirical_values <- c(empirical_values, list(Response = c(value = response_emp, se = response_se)))
      theoretical_values <- c(theoretical_values, Response = response_theo)
      endpoint_names <- c(endpoint_names, "Response")
    }
    
    # Calculate correlations
    if (has_os && has_pfs) {
      # Corr(OS, PFS)
      cor_os_pfs_sims <- sapply(grp_data_by_sim, function(x) cor(x$OS, x$PFS))
      cor_os_pfs_emp <- mean(cor_os_pfs_sims)
      cor_os_pfs_se <- sd(cor_os_pfs_sims)
      cor_os_pfs_theo <- hazard_OS[grp] / hazard_PFS[grp]  # Fleischer model
      empirical_values <- c(empirical_values, list(Cor_OS_PFS = c(value = cor_os_pfs_emp, se = cor_os_pfs_se)))
      theoretical_values <- c(theoretical_values, Cor_OS_PFS = cor_os_pfs_theo)
      endpoint_names <- c(endpoint_names, "Cor_OS_PFS")
    }
    
    if (has_os && has_response) {
      # Corr(OS, Response)
      cor_os_resp_sims <- sapply(grp_data_by_sim, function(x) cor(x$OS, x$Response))
      cor_os_resp_emp <- mean(cor_os_resp_sims)
      cor_os_resp_se <- sd(cor_os_resp_sims)
      cor_os_resp_theo <- rho_tte_resp[grp]  # User-specified
      empirical_values <- c(empirical_values, list(Cor_OS_Response = c(value = cor_os_resp_emp, se = cor_os_resp_se)))
      theoretical_values <- c(theoretical_values, Cor_OS_Response = cor_os_resp_theo)
      endpoint_names <- c(endpoint_names, "Cor_OS_Response")
    }
    
    if (has_pfs && has_response) {
      # Corr(PFS, Response)
      cor_pfs_resp_sims <- sapply(grp_data_by_sim, function(x) cor(x$PFS, x$Response))
      cor_pfs_resp_emp <- mean(cor_pfs_resp_sims)
      cor_pfs_resp_se <- sd(cor_pfs_resp_sims)
      
      # Calculate theoretical value
      if (has_os && has_pfs && has_response) {
        # OS + PFS + Response: Use CorResponsePFS
        cor_pfs_resp_theo <- CorResponsePFS(
          p = p[grp],
          hazard_OS = hazard_OS[grp],
          hazard_PFS = hazard_PFS[grp],
          rho_OS_Response = rho_tte_resp[grp],
          copula = copula
        )
      } else {
        # PFS + Response only: rho_tte_resp directly specifies PFS-Response correlation
        cor_pfs_resp_theo <- rho_tte_resp[grp]
      }
      
      empirical_values <- c(empirical_values, list(Cor_PFS_Response = c(value = cor_pfs_resp_emp, se = cor_pfs_resp_se)))
      theoretical_values <- c(theoretical_values, Cor_PFS_Response = cor_pfs_resp_theo)
      endpoint_names <- c(endpoint_names, "Cor_PFS_Response")
    }
    
    # Create data frame for this group
    result_df <- data.frame(
      Group = grp,
      Endpoint = endpoint_names,
      Empirical = sapply(empirical_values, function(x) x["value"]),
      Theoretical = unlist(theoretical_values),
      SE = sapply(empirical_values, function(x) x["se"]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
    # Calculate performance metrics
    # Bias = Empirical - Theoretical (signed difference)
    result_df$Bias <- result_df$Empirical - result_df$Theoretical
    
    # Relative Bias (%) = 100 * Bias / Theoretical
    result_df$Relative_Bias <- 100 * result_df$Bias / result_df$Theoretical
    
    # MSE = Bias² + SE² (variance of the estimator)
    result_df$MSE <- result_df$Bias^2 + result_df$SE^2
    
    # RMSE = √MSE (root mean squared error in original scale)
    result_df$RMSE <- sqrt(result_df$MSE)
    
    # Add assessment column for quick interpretation
    result_df$Assessment <- ifelse(
      abs(result_df$Relative_Bias) < 5, "Excellent",
      ifelse(abs(result_df$Relative_Bias) < 10, "Acceptable", "Review"))
    
    return(result_df)
  })
  
  # Combine all groups
  result <- do.call(rbind, result_list)
  rownames(result) <- NULL
  
  # Reorder columns for publication-ready format
  result <- result[, c("Group", "Endpoint", "Empirical", "Theoretical", 
                       "Bias", "Relative_Bias", "SE", "MSE", "RMSE", "Assessment")]
  
  # Convert to tibble
  result <- tibble::as_tibble(result)
  
  return(result)
}