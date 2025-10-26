#' Generate Correlated Oncology Endpoints for Multiple Groups
#'
#' @description
#' A unified function that generates random samples of oncology endpoints with
#' flexible configurations. This function supports seven different patterns:
#' \enumerate{
#'   \item OS only
#'   \item PFS only
#'   \item Response only
#'   \item OS + Response (correlated)
#'   \item PFS + Response (correlated)
#'   \item OS + PFS (correlated, with PFS <= OS constraint)
#'   \item OS + PFS + Response (all correlated, with PFS <= OS constraint)
#' }
#'
#' This function is designed specifically for oncology clinical trial simulations,
#' particularly for generating objective response (Response), progression-free
#' survival (PFS), and overall survival (OS) data with appropriate correlation
#' structure.
#'
#' @param nsim Integer. The number of simulations to perform. Must be a
#'   positive integer.
#' @param group Character vector. Names of the treatment groups. If NULL
#'   (default), groups will be named "Group1", "Group2", etc. Length must
#'   match the lengths of other vector parameters.
#' @param n Integer vector. Sample sizes for each group. All elements must be
#'   positive integers. Length must match the number of groups.
#' @param p Numeric vector or NULL. Marginal probabilities of objective
#'   response for each group. If NULL, no Response endpoint is generated. When
#'   provided, all elements must be between 0 and 1 (exclusive). Length must
#'   match the number of groups. Default is NULL.
#' @param hazard_OS Numeric vector or NULL. Hazard rates for OS for each group.
#'   If NULL, no OS endpoint is generated. When provided, all elements must be
#'   positive. Length must match the number of groups. Default is NULL.
#' @param hazard_PFS Numeric vector or NULL. Hazard rates for PFS for each group.
#'   If NULL, no PFS endpoint is generated. When provided, all elements must
#'   be positive. If both hazard_OS and hazard_PFS are provided, hazard_PFS
#'   must be strictly greater than hazard_OS to ensure PFS <= OS. Length must
#'   match the number of groups. Default is NULL.
#' @param rho_tte_resp Numeric vector or NULL. Desired correlation coefficients
#'   between a time-to-event endpoint and Response. Required when both a TTE
#'   endpoint (OS or PFS) and p are provided. All elements must be
#'   within the feasible range determined by Fréchet-Hoeffding bounds (which
#'   depend only on p, not on hazard rates). Length must match the number of
#'   groups. Default is NULL.
#' @param copula Character or NULL. The copula family to use for modeling the
#'   dependence between TTE and Response. Options are "Clayton" or "Frank".
#'   Required when p and at least one TTE endpoint are provided (defaults to
#'   "Clayton" if not specified). Default is NULL.
#'
#' @return A data frame with nsim * sum(n) rows and columns depending on the
#'   configuration:
#'   \itemize{
#'     \item OS only: simID, Group, OS
#'     \item PFS only: simID, Group, PFS
#'     \item Response only: simID, Group, Response
#'     \item OS + Response: simID, Group, OS, Response
#'     \item PFS + Response: simID, Group, PFS, Response
#'     \item OS + PFS: simID, Group, OS, PFS
#'     \item OS + PFS + Response: simID, Group, OS, PFS, Response
#'   }
#'
#' @details
#' This function provides seven modes of operation based on which parameters are
#' provided:
#'
#' **Single endpoint modes:**
#' \itemize{
#'   \item OS only (hazard_OS provided)
#'   \item PFS only (hazard_PFS provided)
#'   \item Response only (p provided)
#' }
#'
#' **Two endpoint modes:**
#' \itemize{
#'   \item OS + Response (hazard_OS + p): Correlated via copula specified
#'     by rho_tte_resp. The copula parameter is determined solely by p and
#'     rho_tte_resp, independent of hazard_OS.
#'   \item PFS + Response (hazard_PFS + p): Correlated via copula specified
#'     by rho_tte_resp. The copula parameter is determined solely by p and
#'     rho_tte_resp, independent of hazard_PFS.
#'   \item OS + PFS (hazard_OS + hazard_PFS): Correlated using the Fleischer model,
#'     ensuring PFS <= OS. The correlation is directly determined by the hazard
#'     ratio: Corr(OS, PFS) = hazard_OS / hazard_PFS.
#' }
#'
#' **Three endpoint mode:**
#' \itemize{
#'   \item OS + PFS + Response (all three provided): The PFS <= OS constraint is
#'     maintained using the Fleischer model. The Response endpoint is correlated
#'     with OS (not PFS) via the copula specified by rho_tte_resp. The copula
#'     parameter depends only on p and rho_tte_resp.
#' }
#'
#' @note
#' **CRITICAL: Correlation structure in three-endpoint mode (OS + PFS + Response):**
#' \itemize{
#'   \item Response is correlated with OS via rho_tte_resp (user-specified)
#'   \item Response is NOT directly correlated with PFS
#'   \item The correlation between PFS and Response emerges automatically from
#'     the model structure and is NOT user-specified
#'   \item To predict the resulting PFS-Response correlation, use
#'     \code{\link{CorResponsePFS}} with your specified parameters
#'   \item This design reflects the biological assumption that response and
#'     survival are linked through the underlying disease process
#' }
#'
#' **Important implementation details:**
#' \itemize{
#'   \item The correlation between a binary Response endpoint and an exponentially
#'     distributed TTE endpoint depends only on the response probability (p) and
#'     the copula parameter, not on the TTE hazard rate. This is a key theoretical
#'     property that simplifies the correlation structure.
#'   \item The feasible range of correlations between Response and TTE endpoints
#'     is determined by Fréchet-Hoeffding bounds, which depend only on p.
#'   \item For the OS-PFS relationship (Fleischer model), the correlation is
#'     directly determined by the ratio of hazard rates.
#' }
#'
#' **Copula selection:**
#' The function uses copula-based methods to model dependencies:
#' \itemize{
#'   \item Clayton copula: Suitable for lower tail dependence (strong correlation
#'     at low values). Cannot model negative dependence.
#'   \item Frank copula: Flexible for both positive and negative dependence with
#'     symmetric tail behavior.
#' }
#'
#' @references
#' Fleischer, F., Gaschler-Markefski, B., & Bluhmki, E. (2009). A statistical
#' model for the dependence between progression-free survival and overall
#' survival. Statistics in Medicine, 28(21), 2669-2686.
#'
#' Trivedi, P. K., & Zimmer, D. M. (2005). Copula modeling: an introduction
#' for practitioners. Foundations and Trends in Econometrics, 1(1), 1-111.
#'
#' @importFrom stats runif rbinom
#'
#' @examples
#' # Example 1: OS only
#' set.seed(123)
#' data1 <- rOncoEndpoints(
#'   nsim = 100,
#'   group = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   hazard_OS = c(0.05, 0.07)
#' )
#' head(data1)
#'
#' # Example 2: Response only
#' set.seed(456)
#' data2 <- rOncoEndpoints(
#'   nsim = 100,
#'   group = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   p = c(0.4, 0.3)
#' )
#' head(data2)
#'
#' # Example 3: OS and Response with Clayton copula
#' set.seed(789)
#' data3 <- rOncoEndpoints(
#'   nsim = 1000,
#'   group = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   p = c(0.4, 0.3),
#'   hazard_OS = c(0.05, 0.07),
#'   rho_tte_resp = c(0.3, 0.2),
#'   copula = "Clayton"
#' )
#' head(data3)
#'
#' # Example 4: PFS and Response with Frank copula
#' set.seed(101)
#' data4 <- rOncoEndpoints(
#'   nsim = 1000,
#'   group = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   p = c(0.5, 0.4),
#'   hazard_PFS = c(0.08, 0.10),
#'   rho_tte_resp = c(0.25, 0.15),
#'   copula = "Frank"
#' )
#' head(data4)
#'
#' # Example 5: OS and PFS (Fleischer model)
#' # Correlation between OS and PFS = hazard_OS / hazard_PFS
#' set.seed(112)
#' data5 <- rOncoEndpoints(
#'   nsim = 1000,
#'   group = c("Treatment", "Control"),
#'   n = c(100, 100),
#'   hazard_OS = c(log(2) / 20, log(2) / 16),
#'   hazard_PFS = c(log(2) / 14, log(2) / 10)
#' )
#' head(data5)
#' # Verify correlation
#' with(subset(data5, Group == "Treatment"), cor(OS, PFS))
#'
#' # Example 6: All three endpoints with Frank copula
#' # IMPORTANT: Response is correlated with OS, not PFS!
#' set.seed(131)
#' data6 <- rOncoEndpoints(
#'   nsim = 1000,
#'   group = "Experimental",
#'   n = 150,
#'   p = 0.5,
#'   hazard_OS = 0.04,
#'   hazard_PFS = 0.06,
#'   rho_tte_resp = 0.3,  # This is the correlation between OS and Response
#'   copula = "Frank"
#' )
#' head(data6)
#' 
#' # Verify correlations
#' cat("OS-Response correlation:", cor(data6$OS, data6$Response), "\n")
#' cat("PFS-Response correlation:", cor(data6$PFS, data6$Response), "\n")
#' 
#' # Predict PFS-Response correlation using CorResponsePFS
#' predicted_cor <- CorResponsePFS(
#'   p = 0.5,
#'   hazard_OS = 0.04,
#'   hazard_PFS = 0.06,
#'   rho_OS_Response = 0.3,
#'   copula = "Frank"
#' )
#' cat("Predicted PFS-Response correlation:", predicted_cor, "\n")
#'
#' @seealso \code{\link{CopulaParamResponseTTE}}, \code{\link{CorBoundResponseTTE}},
#'   \code{\link{CorResponsePFS}}, \code{\link{CheckSimResults}}
#'
#' @export
rOncoEndpoints <- function(nsim, group = NULL, n, p = NULL,
                           hazard_OS = NULL, hazard_PFS = NULL,
                           rho_tte_resp = NULL, copula = NULL) {
  
  # Input validation for nsim
  if (!is.numeric(nsim) || length(nsim) != 1 || nsim <= 0 ||
      nsim != as.integer(nsim)) {
    stop("'nsim' must be a single positive integer")
  }
  nsim <- as.integer(nsim)
  
  # Validate that at least one endpoint is provided
  if (is.null(hazard_OS) && is.null(hazard_PFS) && is.null(p)) {
    stop("At least one of 'hazard_OS', 'hazard_PFS', or 'p' must be provided")
  }
  
  # Determine number of groups
  ngroups <- length(n)
  
  # Input validation for n
  if (!is.numeric(n) || ngroups < 1 || any(n <= 0) ||
      any(n != as.integer(n))) {
    stop("'n' must be a vector of positive integers with at least one element")
  }
  n <- as.integer(n)
  
  # Validate group names
  if (is.null(group)) {
    group <- paste0("Group", seq_len(ngroups))
  } else {
    if (!is.character(group) || length(group) != ngroups) {
      stop(paste0("'group' must be a character vector of length ", ngroups))
    }
  }
  
  # Determine which endpoints to generate
  has_response <- !is.null(p)
  has_os <- !is.null(hazard_OS)
  has_pfs <- !is.null(hazard_PFS)
  
  # Validate hazard_OS
  if (has_os) {
    if (!is.numeric(hazard_OS) || length(hazard_OS) != ngroups || any(hazard_OS <= 0)) {
      stop(paste0("'hazard_OS' must be a numeric vector of length ", ngroups,
                  " with all positive elements"))
    }
  }
  
  # Validate hazard_PFS
  if (has_pfs) {
    if (!is.numeric(hazard_PFS) || length(hazard_PFS) != ngroups || any(hazard_PFS <= 0)) {
      stop(paste0("'hazard_PFS' must be a numeric vector of length ", ngroups,
                  " with all positive elements"))
    }
    
    # Validate hazard_PFS > hazard_OS if both are provided
    if (has_os && any(hazard_PFS <= hazard_OS)) {
      invalid_groups <- which(hazard_PFS <= hazard_OS)
      stop(paste0("'hazard_PFS' must be strictly greater than 'hazard_OS' for all groups.\n",
                  "Violation found in group(s): ",
                  paste(group[invalid_groups], collapse = ", "), "\n",
                  "Current values:\n",
                  paste(sprintf("  %s: hazard_PFS = %.4f, hazard_OS = %.4f", 
                                group[invalid_groups], 
                                hazard_PFS[invalid_groups], 
                                hazard_OS[invalid_groups]), 
                        collapse = "\n"), "\n",
                  "This constraint ensures PFS <= OS in the Fleischer model."))
    }
  }
  
  # Validate p
  if (has_response) {
    if (!is.numeric(p) || length(p) != ngroups ||
        any(p <= 0) || any(p >= 1)) {
      stop(paste0("'p' must be a numeric vector of length ", ngroups,
                  " with all elements in (0, 1)"))
    }
  }
  
  # Validate rho_tte_resp and copula when Response and TTE are both present
  has_tte <- has_os || has_pfs
  if (has_response && has_tte) {
    # Validate rho_tte_resp
    if (is.null(rho_tte_resp)) {
      stop("'rho_tte_resp' must be provided when both Response and a TTE endpoint are provided")
    }
    if (!is.numeric(rho_tte_resp) || length(rho_tte_resp) != ngroups) {
      stop(paste0("'rho_tte_resp' must be a numeric vector of length ", ngroups))
    }
    
    # Validate and set copula
    if (is.null(copula)) {
      copula <- "Clayton"
      message("'copula' was NULL; using default 'Clayton'")
    }
    if (!is.character(copula) || length(copula) != 1 ||
        !copula %in% c("Clayton", "Frank")) {
      stop("'copula' must be either 'Clayton' or 'Frank'")
    }
    
    # Additional check for Clayton copula with negative correlations
    if (copula == "Clayton" && any(rho_tte_resp < 0)) {
      negative_groups <- which(rho_tte_resp < 0)
      stop("Clayton copula cannot model negative dependence.\n",
           "Groups with negative rho_tte_resp: ",
           paste(sprintf("%s (rho = %.3f)", group[negative_groups], 
                         rho_tte_resp[negative_groups]), collapse = ", "), "\n",
           "Use Frank copula instead for negative correlations.")
    }
    
  } else {
    # When Response and TTE are not both present, rho_tte_resp and copula must be NULL
    if (!is.null(rho_tte_resp)) {
      stop("'rho_tte_resp' must be NULL when Response and a TTE endpoint are not both provided")
    }
    if (!is.null(copula)) {
      stop("'copula' must be NULL when Response and a TTE endpoint are not both provided")
    }
  }
  
  # Define conditional copula functions if Response endpoint is included with TTE
  # These are conditional transforms C(u2|u1) derived from the copula
  # (see Trivedi and Zimmer (2005) Table A.1)
  if (has_response && has_tte) {
    if (copula == "Clayton") {
      # Clayton copula conditional transform: C(u2|u1) = ∂C(u1,u2)/∂u1
      fun_cond_copula <- function(v1, v2, theta) {
        (v1 ^ (-theta) * (v2 ^ (-theta / (theta + 1)) - 1) + 1) ^ (-1 / theta)
      }
    } else if (copula == "Frank") {
      # Frank copula conditional transform: C(u2|u1) = ∂C(u1,u2)/∂u1
      fun_cond_copula <- function(v1, v2, theta) {
        -(1 / theta) * log(1 + (v2 * (1 - exp(-theta))) /
                             (v2 * (exp(-theta * v1) - 1) - exp(-theta * v1)))
      }
    }
    
    # Pre-calculate association parameters for all groups using CopulaParamResponseTTE
    # Note: The copula parameter depends only on p and rho, not on the TTE hazard rate
    # This allows early error detection if any rho is infeasible
    theta_vec <- vapply(seq_len(ngroups), function(j) {
      CopulaParamResponseTTE(
        p = p[j],
        rho = rho_tte_resp[j],
        copula = copula
      )
    }, numeric(1))
  }
  
  # Generate data for all groups (vectorized for efficiency)
  group_data_list <- lapply(seq_len(ngroups), function(j) {
    # Total number of observations for this group
    n_total <- n[j] * nsim
    
    # Create simID for this group
    simID <- rep(seq_len(nsim), each = n[j])
    
    # Determine the generation pattern based on which endpoints are requested
    if (has_os && has_pfs && has_response) {
      # Pattern 7: OS + PFS + Response
      # Generate OS, PFS (with PFS <= OS constraint), and Response (correlated with OS)
      
      u1 <- runif(n_total)
      u2 <- runif(n_total)
      
      # Generate OS ~ Exp(hazard_OS)
      OS <- -log(1 - u1) / hazard_OS[j]
      
      # Generate TTP (time to progression) ~ Exp(hazard3)
      # where hazard3 = hazard_PFS - hazard_OS
      # This ensures PFS = min(OS, TTP) ~ Exp(hazard_PFS) (Fleischer model)
      hazard3 <- hazard_PFS[j] - hazard_OS[j]
      u3 <- runif(n_total)
      TTP <- -log(1 - u3) / hazard3
      
      # PFS = min(OS, TTP) ~ Exp(hazard_PFS)
      PFS <- pmin(OS, TTP)
      
      # Generate Response endpoint (correlated with OS via copula)
      # Response = 1 if conditional copula transform > 1 - p, else 0
      Response <- as.double(fun_cond_copula(u1, u2, theta_vec[j]) > 1 - p[j])
      
      data.frame(
        simID = simID,
        Group = group[j],
        OS = OS,
        PFS = PFS,
        Response = Response,
        stringsAsFactors = FALSE
      )
      
    } else if (has_os && has_pfs) {
      # Pattern 6: OS + PFS
      # Generate OS and PFS with PFS <= OS constraint (Fleischer model)
      
      u1 <- runif(n_total)
      
      # Generate OS ~ Exp(hazard_OS)
      OS <- -log(1 - u1) / hazard_OS[j]
      
      # Generate TTP (time to progression) ~ Exp(hazard3)
      # where hazard3 = hazard_PFS - hazard_OS
      hazard3 <- hazard_PFS[j] - hazard_OS[j]
      u3 <- runif(n_total)
      TTP <- -log(1 - u3) / hazard3
      
      # PFS = min(OS, TTP) ~ Exp(hazard_PFS)
      # Correlation: Corr(OS, PFS) = hazard_OS / hazard_PFS
      PFS <- pmin(OS, TTP)
      
      data.frame(
        simID = simID,
        Group = group[j],
        OS = OS,
        PFS = PFS,
        stringsAsFactors = FALSE
      )
      
    } else if (has_os && has_response) {
      # Pattern 4: OS + Response
      # Generate OS and Response with specified correlation via copula
      
      u1 <- runif(n_total)
      u2 <- runif(n_total)
      
      # Generate OS ~ Exp(hazard_OS)
      OS <- -log(1 - u1) / hazard_OS[j]
      
      # Generate Response endpoint (correlated with OS via copula)
      Response <- as.double(fun_cond_copula(u1, u2, theta_vec[j]) > 1 - p[j])
      
      data.frame(
        simID = simID,
        Group = group[j],
        OS = OS,
        Response = Response,
        stringsAsFactors = FALSE
      )
      
    } else if (has_pfs && has_response) {
      # Pattern 5: PFS + Response
      # Generate PFS and Response with specified correlation via copula
      
      u1 <- runif(n_total)
      u2 <- runif(n_total)
      
      # Generate PFS ~ Exp(hazard_PFS)
      PFS <- -log(1 - u1) / hazard_PFS[j]
      
      # Generate Response endpoint (correlated with PFS via copula)
      Response <- as.double(fun_cond_copula(u1, u2, theta_vec[j]) > 1 - p[j])
      
      data.frame(
        simID = simID,
        Group = group[j],
        PFS = PFS,
        Response = Response,
        stringsAsFactors = FALSE
      )
      
    } else if (has_os) {
      # Pattern 1: OS only
      # Generate independent OS observations
      
      u1 <- runif(n_total)
      OS <- -log(1 - u1) / hazard_OS[j]
      
      data.frame(
        simID = simID,
        Group = group[j],
        OS = OS,
        stringsAsFactors = FALSE
      )
      
    } else if (has_pfs) {
      # Pattern 2: PFS only
      # Generate independent PFS observations
      
      u1 <- runif(n_total)
      PFS <- -log(1 - u1) / hazard_PFS[j]
      
      data.frame(
        simID = simID,
        Group = group[j],
        PFS = PFS,
        stringsAsFactors = FALSE
      )
      
    } else {
      # Pattern 3: Response only
      # Generate independent binary Response observations
      
      Response <- rbinom(n_total, 1, p[j])
      
      data.frame(
        simID = simID,
        Group = group[j],
        Response = Response,
        stringsAsFactors = FALSE
      )
    }
  })
  
  # Combine all groups into a single data frame
  result <- do.call(rbind, group_data_list)
  
  # Preserve group order by converting Group to factor with specified levels
  result$Group <- factor(result$Group, levels = group)
  
  # Sort by simulation ID and then by Group (preserving original group order)
  result <- result[order(result$simID, result$Group), ]
  
  # Convert Group back to character
  result$Group <- as.character(result$Group)
  
  # Reset row names
  rownames(result) <- NULL
  
  return(result)
}