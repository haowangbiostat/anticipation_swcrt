#' Compute Statistical Power for Stepped Wedge Cluster Randomized Trials
#' with/without Anticipation
#'
#' `Trt.Ant.Power()` implements the closed‑form power calculation derived from
#' Theorem 4.1 for the four working models in the main article. It supports
#' four models, including HH, HH-ANT, ETI, and ETI-ANT.
#'
#' The analytical variance expressions are independent of the true treatment
#' effect and depend only on trial design parameters and the intra‑cluster
#' correlation coefficient (ICC).
#'
#' @section Referenced helpers:
#' The function relies on the following lower‑level helpers, which must be
#' available in `utils_power.R`:
#' * `build_Z_HH()`, `build_A_HH()`, `build_Z_ETI()`, `build_A_ETI()`
#' * `general_HH()`, `general_HH_ANT()`, `general_ETI()`, `general_ETI_ANT()`
#' * `power_calculation()`, `get_tau()`
#'
#' @param model String. One of `"HH"`, `"HH-ANT"`, `"ETI"`, or `"ETI-ANT"`.
#' @param trt Numeric.  Effect size: constant treatment effect for HH/HH‑ANT;
#'   time-average treatment effect for ETI/ETI‑ANT.
#' @param I Integer.  Number of clusters.
#' @param J Integer.  Number of periods.
#' @param K Integer.  Number of individuals per cluster‑period (assumed equal).
#' @param rho Numeric.  Outcome ICC.
#' @param sigma_sq Numeric. Individual‑level error variance.
#' @param alpha Numeric. Two‑sided type‑I error rate. Default 0.05.
#'
#' @return A scalar of class `"TrtAntPower"`, which contains the estimated power.
#'   The object also carries attributes with the call parameters and the
#'   computed variance, accessible through `attr(x, "details")`.
#' @export
#'
#' @examples
#' Trt.Ant.Power(model    = "ETI-ANT",
#'               trt      = 0.299,
#'               I        = 18,
#'               J        = 7,
#'               K        = 50,
#'               rho      = 0.10,
#'               sigma_sq = 1)
#' Stepped Wedge CRT Power Analysis
#'   Model            : ETI-ANT 
#'   Treatment Effect : 0.299 
#'   Clusters (I)     : 18 
#'   Periods  (J)     : 7 
#'   Individuals (K)  : 50 
#'   ICC (rho)        : 0.1 
#'   Sigma^2          : 1 
#'   Alpha            : 0.05 
#'   
#' Estimated power    : 0.8012236 


source("../code/utils_power.R")
Trt.Ant.Power <- function(model, trt, I, J, K, rho, sigma_sq, alpha = 0.05) {
  
  model <- match.arg(model, c("HH","HH-ANT","ETI","ETI-ANT"))
  
  scalar <- function(x) is.numeric(x) && length(x) == 1L
  stopifnot(
    scalar(trt), scalar(I) && I > 0,
    scalar(J) && J > 0, scalar(K) && K > 0,
    scalar(rho) && rho >= 0 && rho < 1,
    scalar(sigma_sq) && sigma_sq > 0,
    scalar(alpha) && alpha > 0 && alpha < 1
  )
  
  tau <- get_tau(rho, sigma_sq)
  
  builders <- list(
    "HH"      = list(Z = build_Z_HH,  A = NULL,        f = general_HH),
    "HH-ANT"  = list(Z = build_Z_HH,  A = build_A_HH,  f = general_HH_ANT),
    "ETI"     = list(Z = build_Z_ETI, A = NULL,        f = general_ETI),
    "ETI-ANT" = list(Z = build_Z_ETI, A = build_A_ETI, f = general_ETI_ANT)
  )
  
  b   <- builders[[model]]
  Z   <- b$Z(I, J)
  var <- if (is.null(b$A))
    b$f(tau, sqrt(sigma_sq), I, J, K, Z)
  else
    b$f(tau, sqrt(sigma_sq), I, J, K, Z, b$A(I, J))
  
  structure(
    power_calculation(trt, var, alpha),
    class   = c("TrtAntPower","numeric"),
    details = list(model = model, trt = trt, I = I, J = J,
                   K = K, rho = rho, sigma_sq = sigma_sq,
                   alpha = alpha, variance = var)
  )
}


#' @export
print.TrtAntPower <- function(x, digits = getOption("digits"), ...) {
  cat("\nStepped Wedge CRT Power Analysis\n")
  det <- attr(x, "details")
  cat("  Model            :", det$model, "\n")
  cat("  Treatment Effect :", det$trt,   "\n")
  cat("  Clusters (I)     :", det$I,      "\n")
  cat("  Periods  (J)     :", det$J,      "\n")
  cat("  Individuals (K)  :", det$K,      "\n")
  cat("  ICC (rho)        :", det$rho,    "\n")
  cat("  Sigma^2          :", det$sigma_sq, "\n")
  cat("  Alpha            :", det$alpha,  "\n\n")
  cat("Estimated power    :",
      formatC(unclass(x), digits = digits, format = "f"), "\n")
  invisible(x)
}