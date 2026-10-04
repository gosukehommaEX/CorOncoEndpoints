#' CorOncoEndpoints: Correlated Progression-Free Survival, Overall Survival and
#' Objective Response
#'
#' @description
#' The package generates progression-free survival (PFS), overall survival (OS)
#' and binary objective response for oncology trial simulations so that the
#' three endpoints of a patient are dependent in a clinically interpretable way.
#'
#' Each treatment group is described by an object created with
#' \code{\link{OncoArm}}. PFS follows an exponential distribution exactly and the
#' response probability equals the specified objective response rate exactly.
#' PFS and response are linked through a Gaussian copula. OS is obtained from an
#' illness-death structure: a PFS event is a death with a specified probability
#' and otherwise a progression, after which the patient survives for a
#' post-progression time. Two models for OS are available.
#' \describe{
#'   \item{\code{os.model = "idm"}}{Time-constant pre-progression hazards and an
#'     exponential post-progression survival whose hazard may differ from the
#'     pre-progression death hazard and may depend on response. OS is then a
#'     mixture of exponential distributions with a closed-form or one-dimensional
#'     integral survival function.}
#'   \item{\code{os.model = "expexp"}}{Both PFS and OS are exactly exponential.
#'     The pre-progression death hazard decreases over time and the
#'     post-progression hazard is determined by the exponential OS distribution,
#'     so it cannot be specified separately and does not depend on response.}
#' }
#'
#' Main functions:
#' \describe{
#'   \item{\code{\link{OncoArm}}}{Define one treatment group from design inputs
#'     (medians, response rate, association) and calibrate the generator.}
#'   \item{\code{\link{rOncoEndpoints}}}{Generate patient-level data for one or
#'     more groups and many simulated trials, with accrual and dropout.}
#'   \item{\code{\link{CorEndpoints}}}{Pearson correlations among PFS, OS and
#'     response implied by a group.}
#'   \item{\code{\link{SurvEndpoint}}, \code{\link{QuantileEndpoint}}}{Survival,
#'     density, hazard and quantiles of PFS and OS, overall or by response.}
#'   \item{\code{\link{CorBoundPFSResponse}}}{Attainable range of the PFS and
#'     response correlation.}
#'   \item{\code{\link{ExpectedEvents}}, \code{\link{AverageHR}},
#'     \code{\link{RequiredEvents}}}{Expected events, average hazard ratio and
#'     the number of events required by the log-rank test.}
#'   \item{\code{\link{CutoffData}}, \code{\link{EventTime}}}{Data at an analysis
#'     cutoff and event-driven cutoff times.}
#' }
#'
#' @references
#' Fleischer, F., Gaschler-Markefski, B. and Bluhmki, E. (2009). A statistical
#' model for the dependence between progression-free survival and overall
#' survival. \emph{Statistics in Medicine}, 28, 2669--2686.
#' \doi{10.1002/sim.3637}
#'
#' Meller, M., Beyersmann, J. and Rufibach, K. (2019). Joint modeling of
#' progression-free and overall survival and computation of correlation
#' measures. \emph{Statistics in Medicine}, 38, 4270--4289.
#' \doi{10.1002/sim.8295}
#'
#' @keywords internal
#' @useDynLib CorOncoEndpoints, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
