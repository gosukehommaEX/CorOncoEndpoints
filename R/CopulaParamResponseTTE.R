#' Compute Association Parameter for Copulas Between Response and Time-to-Event Endpoints
#'
#' @description
#' Computes the copula association parameter (theta) that corresponds to a
#' specified correlation coefficient between a binary response endpoint and a
#' time-to-event (TTE) endpoint. This function solves for theta using the
#' relationship between the correlation coefficient and the copula function
#' via numerical integration and the secant method.
#'
#' @param p Numeric. The true probability of the binary response endpoint.
#'   Must be between 0 and 1 (exclusive).
#' @param rho Numeric. The desired correlation coefficient between the response
#'   and TTE endpoints. Must be within the feasible range determined by
#'   Fréchet-Hoeffding bounds for the given p.
#' @param copula Character. The copula family to use. Options are "Clayton"
#'   or "Frank" (default is "Clayton").
#'
#' @return A numeric value representing the association parameter theta that
#'   achieves the specified correlation rho.
#'
#' @details
#' This function calculates the copula association parameter theta by solving
#' the equation that relates the correlation coefficient to the copula function:
#' \deqn{\rho = \sqrt{\frac{1-p}{p}} \left[1 + \frac{1}{1-p} \int_0^1 \log(1-u) \cdot c_y(u, 1-p; \theta) \, du\right]}
#'
#' where \eqn{c_y(u, v; \theta) = \frac{\partial C(u, v; \theta)}{\partial u}}
#' is the first partial derivative of the copula function C with respect to
#' its first argument.
#'
#' **Clayton copula**:
#' \deqn{C(u_1, u_2; \theta) = (u_1^{-\theta} + u_2^{-\theta} - 1)^{-1/\theta}}
#' \itemize{
#'   \item Parameter range: θ ∈ (0, ∞)
#'   \item Cannot model negative dependence
#'   \item Exhibits strong lower tail dependence
#'   \item As θ → 0: independence; as θ → ∞: Fréchet upper bound
#' }
#'
#' **Frank copula**:
#' \deqn{C(u_1, u_2; \theta) = -\frac{1}{\theta} \log\left[1 + \frac{(e^{-\theta u_1} - 1)(e^{-\theta u_2} - 1)}{e^{-\theta} - 1}\right]}
#' \itemize{
#'   \item Parameter range: θ ∈ (-∞, ∞)
#'   \item Can model both positive and negative dependence
#'   \item Exhibits symmetric tail dependence
#'   \item As θ → -∞: Fréchet lower bound; θ = 0: independence; θ → ∞: Fréchet upper bound
#' }
#'
#' The function first checks whether the specified rho is within the feasible
#' range using \code{\link{CorBoundResponseTTE}}. It then uses the secant
#' method to iteratively solve for theta. The secant method uses linear
#' interpolation between two points to find where the objective function
#' crosses zero.
#'
#' @note
#' The correlation between a binary response endpoint and an exponentially
#' distributed TTE endpoint depends only on the response probability p and
#' the copula parameter theta, not on the hazard rate of the TTE endpoint.
#' This is because the hazard rate only affects the scale of the TTE variable,
#' while the correlation structure is determined by the copula and the response
#' probability.
#'
#' The secant method algorithm:
#' \enumerate{
#'   \item Start with two initial values theta_0 and theta_1 (adaptively chosen based on rho)
#'   \item Evaluate the objective function at both points: f_0 = f(theta_0), f_1 = f(theta_1)
#'   \item Find the next approximation using linear interpolation:
#'         theta_new = theta_1 - f_1 * (theta_1 - theta_0) / (f_1 - f_0)
#'   \item Check for convergence (either |f_new| or |theta_new - theta_1| is small)
#'   \item If not converged, update: theta_0 <- theta_1, theta_1 <- theta_new, and repeat
#' }
#'
#' @references
#' Trivedi, P. K., & Zimmer, D. M. (2005). Copula modeling: an introduction
#' for practitioners. Foundations and Trends in Econometrics, 1(1), 1-111.
#'
#' Hofert, M., Kojadinovic, I., Maechler, M., & Yan, J. (2018). Elements of
#' copula modeling with R. Springer.
#'
#' Burden, R. L., & Faires, J. D. (2010). Numerical Analysis (9th ed.).
#' Brooks/Cole.
#'
#' @importFrom stats integrate
#'
#' @examples
#' # Clayton copula with positive correlation
#' CopulaParamResponseTTE(p = 0.4, rho = 0.3, copula = "Clayton")
#'
#' # Frank copula with positive correlation
#' CopulaParamResponseTTE(p = 0.4, rho = 0.3, copula = "Frank")
#'
#' # Frank copula with negative correlation
#' CopulaParamResponseTTE(p = 0.4, rho = -0.2, copula = "Frank")
#'
#' # Compare theta values across different copulas for the same correlation
#' rho_target <- 0.3
#' theta_clayton <- CopulaParamResponseTTE(p = 0.4, rho = rho_target, 
#'                                          copula = "Clayton")
#' theta_frank <- CopulaParamResponseTTE(p = 0.4, rho = rho_target, 
#'                                        copula = "Frank")
#' cat("Clayton theta:", theta_clayton, "\nFrank theta:", theta_frank, "\n")
#'
#' # Explore relationship between correlation and theta for different p values
#' rho_seq <- seq(0.1, 0.5, by = 0.1)
#' theta_p40 <- sapply(rho_seq, function(r) {
#'   CopulaParamResponseTTE(p = 0.4, rho = r, copula = "Clayton")
#' })
#' theta_p60 <- sapply(rho_seq, function(r) {
#'   CopulaParamResponseTTE(p = 0.6, rho = r, copula = "Clayton")
#' })
#' 
#' plot(rho_seq, theta_p40, type = "b", col = "blue",
#'      xlab = "Correlation", ylab = "Theta",
#'      main = "Clayton Copula: Correlation vs Theta",
#'      ylim = range(c(theta_p40, theta_p60)))
#' lines(rho_seq, theta_p60, type = "b", col = "red")
#' legend("topleft", legend = c("p = 0.4", "p = 0.6"), 
#'        col = c("blue", "red"), lty = 1, pch = 1)
#'
#' @seealso \code{\link{CorBoundResponseTTE}}
#'
#' @export
CopulaParamResponseTTE <- function(p, rho, copula = "Clayton") {
  # Input validation
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1) {
    stop("'p' must be a single numeric value between 0 and 1 (exclusive)")
  }
  if (!is.numeric(rho) || length(rho) != 1 || abs(rho) >= 1) {
    stop("'rho' must be a single numeric value between -1 and 1 (exclusive)")
  }
  if (!copula %in% c("Clayton", "Frank")) {
    stop("'copula' must be either 'Clayton' or 'Frank'")
  }
  
  # Calculate feasible range for rho using Fréchet-Hoeffding bounds
  rho_bounds <- CorBoundResponseTTE(p = p)
  
  # Check whether rho is within the feasible range
  if (rho <= rho_bounds[1] || rho >= rho_bounds[2]) {
    stop(paste0("The value of rho must be in (",
                round(rho_bounds[1], 3), ", ", round(rho_bounds[2], 3), "). ",
                "Use CorBoundResponseTTE(p = ", p, ") to check feasible bounds."))
  }
  
  # Additional check for Clayton copula (cannot handle negative dependence)
  if (copula == "Clayton" && rho < 0) {
    stop("Clayton copula cannot model negative dependence (rho = ", rho, "). ",
         "Use Frank copula instead for negative correlations.")
  }
  
  # Define first partial derivative of copula function: cy(u, v) = ∂C(u, v)/∂u
  if (copula == "Clayton") {
    ## Clayton copula: first derivative with respect to first argument
    cy_copula <- function(u1, u2, theta) {
      (u1 ^ (-theta) + u2 ^ (-theta) - 1) ^ (-1 / theta - 1) * u1 ^ (-theta - 1)
    }
  } else if (copula == "Frank") {
    ## Frank copula: first derivative with respect to first argument
    cy_copula <- function(u1, u2, theta) {
      exp(-theta * u1) * (exp(-theta * u2) - 1) / 
        ((exp(-theta) - 1) + (exp(-theta * u1) - 1) * (exp(-theta * u2) - 1))
    }
  }
  
  # Define the objective function f(theta) = calculated_rho - target_rho
  # We seek theta such that f(theta) = 0
  objective_function <- function(theta) {
    integral_result <- integrate(
      function(u, theta) {
        log(1 - u) * cy_copula(u, 1 - p, theta)
      },
      lower = 0,
      upper = 1,
      stop.on.error = FALSE,
      theta = theta
    )
    
    calculated_rho <- sqrt((1 - p) / p) * (1 + 1 / (1 - p) * integral_result$value)
    return(calculated_rho - rho)
  }
  
  # Solve for theta using the secant method
  # Step 0: Set algorithm parameters
  max_iter <- 100
  tolerance <- 1e-6
  
  # Adaptive initial value selection based on correlation magnitude
  if (copula == "Clayton") {
    # For Clayton copula: theta must be positive
    # Larger rho requires larger theta
    theta_0 <- 0.1
    theta_1 <- ifelse(abs(rho) < 0.3, 0.5, 
                      ifelse(abs(rho) < 0.6, 2.0, 5.0))
  } else if (copula == "Frank") {
    # For Frank copula: theta can be any real number
    # Choose initial values based on sign and magnitude of target correlation
    if (rho < 0) {
      theta_0 <- ifelse(abs(rho) < 0.3, -0.5, 
                        ifelse(abs(rho) < 0.6, -2.0, -5.0))
      theta_1 <- -0.1
    } else {
      theta_0 <- 0.1
      theta_1 <- ifelse(rho < 0.3, 0.5, 
                        ifelse(rho < 0.6, 2.0, 5.0))
    }
  }
  
  # Step 1: Evaluate objective function at initial values
  f_0 <- objective_function(theta_0)
  f_1 <- objective_function(theta_1)
  
  # Iterate using the secant method until convergence
  for (iter in 1:max_iter) {
    # Check if denominator is too small (avoid division by zero)
    if (abs(f_1 - f_0) < 1e-10) {
      # Expand the search range if function values are too close
      if (copula == "Clayton" || (copula == "Frank" && rho > 0)) {
        theta_1 <- theta_1 * 2
      } else {
        theta_1 <- theta_1 * 2  # More negative for Frank copula with rho < 0
      }
      f_1 <- objective_function(theta_1)
      next
    }
    
    # Step 2: Linear interpolation to find the next approximation
    # Find where the line through (theta_0, f_0) and (theta_1, f_1) crosses zero
    theta_new <- theta_1 - f_1 * (theta_1 - theta_0) / (f_1 - f_0)
    
    # For Clayton copula, ensure theta remains positive
    if (copula == "Clayton" && theta_new <= 0) {
      theta_new <- theta_1 * 0.5
    }
    
    # Step 3: Evaluate objective function at new theta
    f_new <- objective_function(theta_new)
    
    # Step 4: Check for convergence
    if (abs(f_new) < tolerance || abs(theta_new - theta_1) < tolerance) {
      return(theta_new)
    }
    
    # Update values for next iteration (discard oldest point, shift values)
    theta_0 <- theta_1
    f_0 <- f_1
    theta_1 <- theta_new
    f_1 <- f_new
  }
  
  # Step 5: If maximum iterations reached without convergence
  # Return best approximation with warning instead of stopping
  warning("Secant method did not converge after ", max_iter, " iterations ",
          "(final |f(theta)| = ", round(abs(f_1), 6), "). ",
          "Returning last iteration value. Consider: ",
          "(1) checking if rho is within feasible bounds using CorBoundResponseTTE(p = ", p, "), ",
          "(2) trying a different copula, or (3) adjusting tolerance.")
  return(theta_1)
}