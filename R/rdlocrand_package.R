###################################################################
# rdlocrand: Local Randomization Methods for RD Designs
# !version 1.1 22-May-2025
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################


#' rdlocrand: Local Randomization Methods for RD Designs
#'
#' The regression discontinuity (RD) design is a popular quasi-experimental design
#' for causal inference and policy evaluation. Under the local randomization approach,
#' RD designs can be interpreted as randomized experiments inside a window around the
#' cutoff. The \code{rdlocrand} package provides tools to analyze RD designs under local
#' randomization: \code{\link{rdrandinf}} to perform hypothesis
#' testing using randomization inference, \code{\link{rdwinselect}} to select a window
#' around the cutoff in which randomization is likely to hold, \code{\link{rdsensitivity}}
#' to assess the sensitivity of the results to different window lengths and null hypotheses
#' and \code{\link{rdrbounds}} to construct Rosenbaum bounds for sensitivity to
#' unobserved confounders. For more details, and related \code{Stata} and \code{R} packages
#'  useful for analysis of RD designs, visit \url{https://rdpackages.github.io/}.
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}
#'
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}
#'
#' Gonzalo Vazquez-Bare, UC Santa Barbara. \email{gvazquez@econ.ucsb.edu}
#'
#' @references
#' Cattaneo, M.D., B. Frandsen and R. Titiunik. (2015).  \href{https://rdpackages.github.io/references/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to Party Advantages in the U.S. Senate}. \emph{Journal of Causal Inference} 3(1): 1-24.
#'
#' Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2016). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2016_Stata.pdf}{Inference in Regression Discontinuity Designs under Local Randomization}. \emph{Stata Journal} 16(2): 331-367.
#'
#' Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2017). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2017_JPAM.pdf}{Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality}. \emph{Journal of Policy Analysis and Management} 36(3): 643-681.
#'
#' Rosenbaum, P. (2002). Observational Studies. Springer.
#'
#' @importFrom grDevices gray.colors
#' @importFrom graphics filled.contour
#' @importFrom graphics title
#' @importFrom stats binom.test
#' @importFrom stats complete.cases
#' @importFrom stats cov
#' @importFrom stats ks.test
#' @importFrom stats lm
#' @importFrom stats pf
#' @importFrom stats pnorm
#' @importFrom stats poly
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats var
#' @importFrom stats vcov
#' @importFrom stats wilcox.test
#'
#' @aliases rdlocrand_package
"_PACKAGE"
