#' Calculate Correlation Between Response and PFS in OS-PFS-Response Framework
#'
#' @description
#' Calculates the correlation coefficient between a binary response endpoint
#' and progression-free survival (PFS) in the specific context where all three
#' endpoints (OS, PFS, Response) are modeled together using the Fleischer model
#' and copula-based dependence. Given a user-specified correlation between OS
#' and Response, this function predicts the resulting correlation between PFS
#' and Response.
#'
#' @param p Numeric. The true probability of the binary response endpoint.
#'   Must be between 0 and 1 (exclusive).
#' @param hazard_OS Numeric. The hazard rate for overall survival (OS).
#'   Must be positive.
#' @param hazard_PFS Numeric. The hazard rate for progression-free survival (PFS).
#'   Must be positive and strictly greater than hazard_OS to ensure PFS <= OS.
#' @param rho_OS_Response Numeric. The desired correlation coefficient between
#'   OS and Response. Must be within the feasible range determined by
#'   Fréchet-Hoeffding bounds for the given p.
#' @param copula Character. The copula family to use for modeling the dependence
#'   between OS and Response. Options are "Clayton" or "Frank".
#'
#' @return A numeric value representing the correlation coefficient between
#'   PFS and Response.
#'
#' @details
#' This function is designed for the specific setting where:
#' \itemize{
#'   \item OS ~ Exp(hazard_OS)
#'   \item TTP (time to progression) ~ Exp(hazard_TTP) where hazard_TTP = hazard_PFS - hazard_OS
#'   \item PFS = min(OS, TTP) ~ Exp(hazard_PFS) (Fleischer model)
#'   \item Response is correlated with OS via the specified copula with correlation rho_OS_Response
#' }
#'
#' The calculation proceeds in three steps:
#' \enumerate{
#'   \item Calculate the copula parameter theta from the specified rho_OS_Response
#'     using \code{\link{CopulaParamResponseTTE}}
#'   \item Define the first partial derivative of the copula function cy(u, v; theta)
#'   \item Calculate Corr(PFS, Response) using the formula:
#'     \deqn{\text{Corr}(PFS, Response) = \sqrt{\frac{1-p}{p}} \frac{1}{\lambda_{TTP}} \left[\frac{\lambda_{PFS}}{1-p} \int_0^1 c_y(u, 1-p; \theta) (1-u)^{\lambda_{TTP}/\lambda_{OS}} du - \lambda_{OS}\right]}
#' }
#'
#' where \eqn{\lambda_{TTP} = \lambda_{PFS} - \lambda_{OS}} and
#' \eqn{c_y(u, v; \theta) = \frac{\partial C(u, v; \theta)}{\partial u}}
#' is the first partial derivative of the copula function.
#'
#' **Key insight**: When generating OS + PFS + Response simultaneously, users
#' only specify the correlation between OS and Response. The correlation between
#' PFS and Response is automatically determined by the model structure. This
#' function calculates what that resulting correlation will be.
#'
#' @note
#' **Practical usage**:
#' \itemize{
#'   \item In \code{\link{rOncoEndpoints}}, when all three endpoints are generated,
#'     users specify rho_tte_resp which represents the correlation between OS
#'     and Response
#'   \item The correlation between PFS and Response is not specified but emerges
#'     from the model structure
#'   \item This function predicts what the PFS-Response correlation will be, which
#'     is useful for:
#'     - Trial design and planning
#'     - Understanding the correlation structure among all three endpoints
#'     - Sensitivity analyses
#'     - Power calculations involving PFS
#' }
#'
#' **Copula families**:
#' \itemize{
#'   \item Clayton copula: Exhibits lower tail dependence, cannot model negative
#'     dependence
#'   \item Frank copula: Flexible for both positive and negative dependence with
#'     symmetric tail behavior
#' }
#'
#' **Numerical considerations**:
#' \itemize{
#'   \item The function uses numerical integration which may be less accurate for
#'     extreme parameter combinations (e.g., very large hazard ratios)
#'   \item Warnings are issued when hazard ratios are extreme or when integration
#'     encounters difficulties
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
#' @importFrom stats integrate
#'
#' @examples
#' # Example 1: Calculate PFS-Response correlation given OS-Response correlation
#' # Median OS = 30 months, Median PFS = 18 months
#' # Response rate = 40%, OS-Response correlation = 0.3
#' rho_pfs_resp <- CorResponsePFS(
#'   p = 0.4,
#'   hazard_OS = log(2) / 30,
#'   hazard_PFS = log(2) / 18,
#'   rho_OS_Response = 0.3,
#'   copula = "Clayton"
#' )
#' cat("PFS-Response correlation:", rho_pfs_resp, "\n")
#'
#' # Example 2: Compare correlations for different copulas
#' p <- 0.5
#' hazard_OS <- 0.04
#' hazard_PFS <- 0.08
#' rho_OS <- 0.25
#'
#' rho_pfs_clayton <- CorResponsePFS(p, hazard_OS, hazard_PFS, rho_OS, "Clayton")
#' rho_pfs_frank <- CorResponsePFS(p, hazard_OS, hazard_PFS, rho_OS, "Frank")
#'
#' cat("Clayton copula - PFS-Response correlation:", rho_pfs_clayton, "\n")
#' cat("Frank copula - PFS-Response correlation:", rho_pfs_frank, "\n")
#'
#' # Example 3: Sensitivity analysis - how does PFS-Response correlation
#' # change with different OS-Response correlations?
#' p <- 0.4
#' hazard_OS <- log(2) / 30
#' hazard_PFS <- log(2) / 18
#'
#' # Calculate feasible range for OS-Response correlation
#' rho_bounds <- CorBoundResponseTTE(p = p)
#'
#' # For Clayton copula: use only positive correlations
#' rho_OS_seq <- seq(0.05, rho_bounds[2] - 0.05, length.out = 10)
#'
#' # Calculate corresponding PFS-Response correlations
#' rho_PFS_seq <- sapply(rho_OS_seq, function(rho) {
#'   CorResponsePFS(p, hazard_OS, hazard_PFS, rho, copula = "Clayton")
#' })
#'
#' # Plot the relationship
#' plot(rho_OS_seq, rho_PFS_seq, type = "b",
#'      xlab = "Correlation: OS-Response",
#'      ylab = "Correlation: PFS-Response",
#'      main = "Relationship between OS-Response and PFS-Response Correlations",
#'      ylim = range(c(rho_OS_seq, rho_PFS_seq)))
#' abline(a = 0, b = 1, lty = 2, col = "gray")  # Identity line
#' legend("topleft", 
#'        legend = c("PFS-Response", "Identity line (for reference)"),
#'        lty = c(1, 2), col = c("black", "gray"), pch = c(1, NA))
#'
#' # Example 4: Effect of hazard ratio on PFS-Response correlation
#' # Fix OS-Response correlation, vary PFS hazard
#' p <- 0.5
#' hazard_OS <- 0.04
#' rho_OS <- 0.3
#' hazard_PFS_vec <- seq(0.06, 0.15, by = 0.01)
#'
#' rho_PFS_vec <- sapply(hazard_PFS_vec, function(h_pfs) {
#'   CorResponsePFS(p, hazard_OS, h_pfs, rho_OS, copula = "Clayton")
#' })
#'
#' # Plot
#' hazard_ratio <- hazard_PFS_vec / hazard_OS
#' plot(hazard_ratio, rho_PFS_vec, type = "b",
#'      xlab = "Hazard Ratio (PFS/OS)",
#'      ylab = "Correlation",
#'      main = "Effect of Hazard Ratio on PFS-Response Correlation",
#'      ylim = range(c(rho_PFS_vec, rho_OS)))
#' abline(h = rho_OS, lty = 2, col = "red", lwd = 2)
#' legend("topright",
#'        legend = c("PFS-Response correlation", 
#'                   paste0("OS-Response correlation (", rho_OS, ")")),
#'        lty = c(1, 2), col = c("black", "red"), pch = c(1, NA), lwd = c(1, 2))
#'
#' # Example 5: Verify calculated correlation matches simulation
#' \dontrun{
#' set.seed(123)
#' p <- 0.4
#' hazard_OS <- 0.05
#' hazard_PFS <- 0.08
#' rho_OS_Response <- 0.3
#'
#' # Predict PFS-Response correlation
#' predicted_rho <- CorResponsePFS(
#'   p, hazard_OS, hazard_PFS, rho_OS_Response, copula = "Clayton"
#' )
#'
#' # Generate data and calculate empirical correlation
#' data <- rOncoEndpoints(
#'   nsim = 1000,
#'   n = 500,
#'   p = p,
#'   hazard_OS = hazard_OS,
#'   hazard_PFS = hazard_PFS,
#'   rho_tte_resp = rho_OS_Response,
#'   copula = "Clayton"
#' )
#'
#' empirical_rho <- cor(data$PFS, data$Response)
#'
#' cat("Predicted PFS-Response correlation:", predicted_rho, "\n")
#' cat("Empirical PFS-Response correlation:", empirical_rho, "\n")
#' cat("Difference:", abs(predicted_rho - empirical_rho), "\n")
#' }
#'
#' @seealso
#' \code{\link{CopulaParamResponseTTE}} for computing copula parameters,
#' \code{\link{CorBoundResponseTTE}} for general TTE-Response correlation bounds,
#' \code{\link{CorBoundResponsePFS}} for theoretical bounds of PFS-Response correlation,
#' \code{\link{rOncoEndpoints}} for generating correlated oncology endpoints
#'
#' @export
CorResponsePFS <- function(p, hazard_OS, hazard_PFS, rho_OS_Response, copula) {
  
  # Input validation
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
    stop("'p' must be a single numeric value between 0 and 1 (exclusive)")
  }
  
  if (!is.numeric(hazard_OS) || length(hazard_OS) != 1 || hazard_OS <= 0) {
    stop("'hazard_OS' must be a single positive numeric value")
  }
  
  if (!is.numeric(hazard_PFS) || length(hazard_PFS) != 1 || hazard_PFS <= 0) {
    stop("'hazard_PFS' must be a single positive numeric value")
  }
  
  # Validate hazard_PFS > hazard_OS to ensure PFS <= OS
  if (hazard_PFS <= hazard_OS) {
    stop("'hazard_PFS' (", hazard_PFS, ") must be strictly greater than 'hazard_OS' (", 
         hazard_OS, ") to ensure PFS <= OS.\n",
         "In the Fleischer model, PFS = min(OS, TTP), which requires hazard_PFS > hazard_OS.")
  }
  
  if (!is.numeric(rho_OS_Response) || length(rho_OS_Response) != 1 || 
      abs(rho_OS_Response) >= 1) {
    stop("'rho_OS_Response' must be a single numeric value between -1 and 1 (exclusive)")
  }
  
  if (!is.character(copula) || length(copula) != 1 || 
      !copula %in% c("Clayton", "Frank")) {
    stop("'copula' must be either 'Clayton' or 'Frank'")
  }
  
  # Check if rho_OS_Response is within the feasible range
  rho_bounds <- CorBoundResponseTTE(p = p)
  if (rho_OS_Response <= rho_bounds[1] || rho_OS_Response >= rho_bounds[2]) {
    stop("'rho_OS_Response' (", round(rho_OS_Response, 3), 
         ") must be within the feasible range (",
         round(rho_bounds[1], 3), ", ", round(rho_bounds[2], 3), 
         ") for the given 'p' (", p, ").\n",
         "Use CorBoundResponseTTE(p = ", p, ") to check feasible bounds.")
  }
  
  # Additional check for Clayton copula (cannot handle negative dependence)
  if (copula == "Clayton" && rho_OS_Response < 0) {
    stop("Clayton copula cannot model negative dependence (rho_OS_Response = ", 
         rho_OS_Response, "). Use Frank copula instead.")
  }
  
  # Validate hazard ratio is in reasonable range
  hazard_ratio <- hazard_PFS / hazard_OS
  if (hazard_ratio > 10) {
    warning("hazard_PFS/hazard_OS ratio (", round(hazard_ratio, 2), 
            ") is very large. This indicates most patients progress before death. ",
            "Results may be numerically unstable. Consider checking your hazard rate specifications.")
  } else if (hazard_ratio < 1.01) {
    warning("hazard_PFS/hazard_OS ratio (", round(hazard_ratio, 2), 
            ") is very close to 1. This indicates very few progressions before death. ",
            "Consider if this reflects your intended clinical scenario.")
  }
  
  # Step 1: Calculate copula parameter θ_0 from ρ_OS_Response
  # This uses the relationship between OS and Response (Theorem 1)
  theta_0 <- CopulaParamResponseTTE(
    p = p,
    rho = rho_OS_Response,
    copula = copula
  )
  
  # Step 2: Define first partial derivative of copula function cy(u, v; θ_0)
  # This is ∂C(u, v; θ)/∂u evaluated at θ = θ_0
  if (copula == "Clayton") {
    # Clayton copula: first derivative with respect to first argument
    cy_copula <- function(u1, u2) {
      (u1 ^ (-theta_0) + u2 ^ (-theta_0) - 1) ^ (-1 / theta_0 - 1) * 
        u1 ^ (-theta_0 - 1)
    }
  } else if (copula == "Frank") {
    # Frank copula: first derivative with respect to first argument
    cy_copula <- function(u1, u2) {
      exp(-theta_0 * u1) * (exp(-theta_0 * u2) - 1) / 
        ((exp(-theta_0) - 1) + (exp(-theta_0 * u1) - 1) * (exp(-theta_0 * u2) - 1))
    }
  }
  
  # Calculate hazard rate for TTP (time to progression)
  # In Fleischer model: PFS = min(OS, TTP)
  # hazard_PFS = hazard_OS + hazard_TTP
  hazard_TTP <- hazard_PFS - hazard_OS
  
  # Calculate power exponent for the integrand
  power_exponent <- hazard_TTP / hazard_OS
  
  # Step 3: Calculate Corr(PFS, Response) using Theorem 4
  # The formula involves an integral of cy(u, 1-p) weighted by (1-u)^(hazard_TTP/hazard_OS)
  integral_result <- integrate(
    function(u) {
      # Add small safeguards for numerical stability
      u_safe <- pmin(pmax(u, 1e-10), 1 - 1e-10)
      cy_copula(u_safe, 1 - p) * (1 - u_safe) ^ power_exponent
    },
    lower = 0,
    upper = 1,
    stop.on.error = FALSE,
    subdivisions = 1000L  # Increase subdivisions for better accuracy
  )
  
  # Check integration status and provide warnings if needed
  if (integral_result$message != "OK") {
    warning("Numerical integration encountered an issue: ", integral_result$message, 
            ". Results may have reduced accuracy. ",
            "Absolute error estimate: ", round(integral_result$abs.error, 8))
  }
  
  # Check if integral value is reasonable
  if (!is.finite(integral_result$value)) {
    stop("Numerical integration failed to produce a finite result. ",
         "This may occur with extreme parameter combinations. ",
         "Please check your input parameters.")
  }
  
  integral_value <- integral_result$value
  
  # Calculate the correlation coefficient between PFS and Response
  # Formula: sqrt((1-p)/p) * (1/hazard_TTP) * [(hazard_PFS/(1-p)) * integral - hazard_OS]
  rho_PFS_Response <- sqrt((1 - p) / p) * (1 / hazard_TTP) * 
    ((hazard_PFS / (1 - p)) * integral_value - hazard_OS)
  
  # Validate result is within theoretical bounds
  if (!is.finite(rho_PFS_Response)) {
    stop("Calculated correlation is not finite. This indicates numerical issues. ",
         "Please check your input parameters.")
  }
  
  # Check against theoretical bounds
  theoretical_bounds <- CorBoundResponsePFS(p = p, hazard_OS = hazard_OS, 
                                            hazard_PFS = hazard_PFS)
  
  if (rho_PFS_Response < theoretical_bounds[1] - 0.01 || 
      rho_PFS_Response > theoretical_bounds[2] + 0.01) {
    warning("Calculated PFS-Response correlation (", round(rho_PFS_Response, 4), 
            ") is outside theoretical bounds [", 
            round(theoretical_bounds[1], 4), ", ", 
            round(theoretical_bounds[2], 4), "]. ",
            "This may indicate numerical issues or extreme parameter values.")
  }
  
  return(rho_PFS_Response)
}