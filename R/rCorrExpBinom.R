#' Generate Correlated Time-to-Event and Binary Outcomes
#'
#' @description
#' This function generates pseudorandom numbers of correlated time-to-event and binary outcomes
#' using copula-based approaches. It allows researchers to simulate datasets with specified
#' correlation structures between survival times and binary response variables.
#'
#' @param hazard Hazard rate (= log(2) / median survival time).
#' @param p Probability of patients being a responder.
#' @param rho Correlation between time-to-event and binary outcomes.
#' @param n Sample size.
#' @param copula.fun Copula function to generate pseudorandom numbers (options: 'Clayton' or 'Frank').
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{hazard}{True hazard rate (= log(2) / median survival time).}
#'   \item{MST}{True median survival time.}
#'   \item{p}{True probability of responder.}
#'   \item{rho}{True correlation between two outcomes.}
#'   \item{UL.rho}{Upper bound of the correlation.}
#'   \item{n}{Sample size.}
#'   \item{copula.fun}{Specified copula function.}
#'   \item{theta}{Estimated association parameter.}
#'   \item{simdata}{Simulated dataset including correlated time-to-event and binary outcomes.}
#'   \item{hat.MST}{Estimated median survival time for simulated dataset.}
#'   \item{hat.p}{Estimated probability of patients being a responder for simulated dataset.}
#'   \item{hat.rho}{Estimated correlation between two outcomes for simulated dataset.}
#' }
#'
#' @details
#' Since correlation of paired time-to-event and binary outcomes is not free to range over [-1, 1],
#' the boundary of correlation can be calculated by Frechet-Hoeffding bounds. 
#' This function currently only covers positive correlation cases.
#'
#' @note
#' The current implementation supports only Clayton and Frank copulas.
#'
#' @references
#' Trivedi PK and Zimmer DM. (2007) Copula modeling: an introduction for practitioners. 
#' Found Trends Econ, 1: 1-111.
#'
#' Hofert M, Kojadinovic I, Machler M, et al. (2018) Elements of Copula Modeling with R. 
#' Use R!, 1: 1-267.
#'
#' Bender R, Augustin T, and Blettner M. (2005) Generating survival times to simulate 
#' Cox proportional hazards models. Statistics in Medicine, 24: 1713-1723.
#'
#' @examples
#' \dontrun{
#' # Generate 1000 observations with median survival time of 12 months,
#' # response probability of 0.2, and correlation of 0.5 using Clayton copula
#' sim_data <- rCorrExpBinom(hazard = log(2) / 12, p = 0.2, rho = 0.5, 
#'                           n = 1000, copula.fun = 'Clayton')
#' 
#' # Access the simulated dataset
#' head(sim_data$simdata[[1]])
#' 
#' # View summary statistics
#' sim_data[c("hat.MST", "hat.p", "hat.rho")]
#' }
#'
#' @export
#'
#' @importFrom dplyr %>% tibble select reframe
#' @importFrom purrr pmap
rCorrExpBinom = function(hazard, p, rho, n, copula.fun) {
  # Upper bound of the correlation (UL.rho) calculated by Frechet-Hoeffding bounds
  # (see Section 2.1.2 in Trivedi and Zimmer (2005))
  UL.cov = '+'(
    '*'(
      pbinom(0, 1, p),
      integrate(
        function(y) 1 - pexp(y, hazard), 
        lower = qexp(pbinom(0, 1, p), hazard),  
        upper = Inf,
        stop.on.error = FALSE
      )[['value']]
    ),
    '*'(
      1 - pbinom(0, 1, p),
      integrate(
        function(y) pexp(y, hazard), 
        lower = 0, 
        upper = qexp(pbinom(0, 1, p), hazard),
        stop.on.error = FALSE
      )[['value']]
    ))
  UL.rho = UL.cov / (sqrt(1 / (hazard ^ 2)) * sqrt(p * (1 - p)))
  # Check whether rho is in (0, UL.rho)
  if(rho <= 0 | rho >= UL.rho) {
    stop(paste0('The value of rho must be rho ', 'in ', '(', 0, ',', round(UL.rho, 3), ').'))
  }
  # Generate pseudorandom numbers
  tibble(
    hazard = hazard,
    MST = log(2) / hazard,
    p = p,
    rho = rho,
    UL.rho = UL.rho,
    n = n,
    copula.fun = copula.fun,
    # Copula functions 
    # (see Table 2.1 in Trivedi and Zimmer (2005))
    f.copula = pmap(
      list(copula.fun), ~ ifelse(
        # Clayton copula
        copula.fun == 'Clayton',
        function(u1, u2, theta) { 
          (u1 ^ (-theta) + u2 ^ (-theta) - 1) ^ (-1 / theta) 
        },
        # Frank copula
        function(u1, u2, theta) { 
          -(1 / theta) * log(1 + (exp(-theta * u1) - 1) * (exp(-theta * u2) - 1) / (exp(-theta) - 1)) 
        }
      )
    ),
    # Solve an equation of 'rho = correlation' by Hoeffding's formula 
    # (see formula (2.19) in Hofert et al. (2018)) for theta
    solve.theta = pmap(
      list(copula.fun), ~ function(theta) {
        '-'(
          '/'(
            integrate(function(x2, theta) {
              '+'(
                f.copula[[1]](1 - p, pexp(x2, hazard), theta),
                f.copula[[1]](1, pexp(x2, hazard), theta) - (2 - p) * pexp(x2, hazard)
              )
            }, 
            lower = 0, 
            upper = Inf,
            stop.on.error = FALSE,
            theta = theta)[['value']],
            (sqrt(1 / (hazard ^ 2)) * sqrt(p * (1 - p)))
          ),
          rho
        )
      }
    ),
    # Return association parameter (i.e., theta), which is corresponding to assumed correlation value via a copula function
    # For Clayton copula 
    # (Note: it is described in Section 2.3.5 in Trivedi and Zimmer (2005) that 0 < theta < Inf,
    # and "The Clayton copula cannot account for negative dependence.")
    theta = if(copula.fun == 'Clayton') {
      theta0 = NA
      U.theta = 1
      while(is.na(theta0)) {
        tryCatch(
          { theta0 = uniroot(solve.theta[[1]], interval = c(1e-5, U.theta))[['root']] }
          , error = function(e) { theta0 = NA }
        )
        U.theta = U.theta + 1
      }
      theta0
      # For Frank copula
      # (Note: it is described in Section 2.3.6 in Trivedi and Zimmer (2005) that -Inf < theta < Inf,
      # and "[...] unlike some other copulas, it permits negative dependence between the marginals.")
    } else if(copula.fun == 'Frank') {
      theta1 = theta2 = NA
      L.theta = -1
      U.theta = 1
      while(is.na(theta1) & is.na(theta2)) {
        tryCatch(
          { theta1 = uniroot(solve.theta[[1]], interval = c(L.theta, -1e-5))[['root']] }
          , error = function(e) { theta1 = NA }
        )
        L.theta = L.theta - 1
        tryCatch(
          { theta2 = uniroot(solve.theta[[1]], interval = c(1e-5, U.theta))[['root']] }
          , error = function(e) { theta2 = NA }
        )
        U.theta = U.theta + 1
      }
      c(theta1, theta2)[!is.na(c(theta1, theta2))]
    },
    # Conditional copula 
    # (see Table A.1 in Trivedi and Zimmer (2005))
    cond.copula = pmap(
      list(copula.fun), ~ ifelse(
        # Clayton copula
        copula.fun == 'Clayton',
        function(v1, v2, theta) {
          (v1 ^ (-theta) * (v2 ^ (-theta / (theta + 1)) - 1) + 1) ^ (-1 / theta)
        },
        # Frank copula
        function(v1, v2, theta) {
          -(1 / theta) * log(1 + (v2 * (1 - exp(-theta))) / (v2 * (exp(-theta * v1) - 1) - exp(-theta * v1)))
        }
      )
    ),
    # Simulated dataset of correlated time-to-event and binary outcomes for each patient
    simdata = list(
      tibble(
        u1 = runif(n, 0, 1),
        u2 = runif(n, 0, 1),
        # Generated exponential distributed random number by uniform distribution 
        # (see Bender et al. (2005)) 
        survtime = -log(1 - u1) / hazard,
        binary = as.double(cond.copula[[1]](u1, u2, theta) > 1 - p)
      ) %>% 
        select(survtime, binary)
    ),
    # Summary statistics (median survival, proportion fo responders and correlation) for simulated dataset
    simdata[[1]] %>% 
      reframe(
        hat.MST = median(survtime), 
        hat.p = mean(binary), 
        hat.rho = cor(survtime, binary)
      )
  ) %>% 
    select(-f.copula, -solve.theta, -cond.copula)
}