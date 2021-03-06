% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdwinselect.R
\name{rdwinselect}
\alias{rdwinselect}
\title{Window selection for RD designs under local randomization}
\usage{
rdwinselect(
  R,
  X,
  cutoff = 0,
  obsmin = NULL,
  wmin = NULL,
  wobs = NULL,
  wstep = NULL,
  wasymmetric = FALSE,
  wmasspoints = FALSE,
  dropmissing = FALSE,
  nwindows = 10,
  statistic = "diffmeans",
  p = 0,
  evalat = "cutoff",
  kernel = "uniform",
  approx = FALSE,
  level = 0.15,
  reps = 1000,
  seed = 666,
  plot = FALSE,
  quietly = FALSE,
  obsstep = NULL
)
}
\arguments{
\item{R}{a vector containing the values of the running variable.}

\item{X}{the matrix of covariates to be used in the balancing tests. The matrix is optional but the recommended window is only provided when at least one covariate is specified. This should be a matrix of size n x k where n is the total sample size and $k$ is the number of covariates.}

\item{cutoff}{the RD cutoff (default is 0).}

\item{obsmin}{the minimum number of observations above and below the cutoff in the smallest window. Default is 10.}

\item{wmin}{the smallest window to be used.}

\item{wobs}{the number of observations to be added at each side of the cutoff at each step. Default is 5.}

\item{wstep}{the increment in window length.}

\item{wasymmetric}{allows for asymmetric windows around the cutoff when (\code{wobs} is specified).}

\item{wmasspoints}{specifies that the running variable is discrete and each masspoint should be used as a window.}

\item{dropmissing}{drop rows with missing values in covariates when calculating windows.}

\item{nwindows}{the number of windows to be used. Default is 10.}

\item{statistic}{the statistic to be used in the balance tests. Allowed options are \code{diffmeans} (difference in means statistic), \code{ksmirnov} (Kolmogorov-Smirnov statistic), \code{ranksum} (Wilcoxon-Mann-Whitney standardized statistic) and \code{hotelling} (Hotelling's T-squared statistic). Default option is \code{diffmeans}. The statistic \code{ttest} is equivalent to \code{diffmeans} and included for backward compatibility.}

\item{p}{the order of the polynomial for outcome adjustment model (for covariates). Default is 0.}

\item{evalat}{specifies the point at which the adjusted variable is evaluated. Allowed options are \code{cutoff} and \code{means}. Default is \code{cutoff}.}

\item{kernel}{specifies the type of kernel to use as weighting scheme. Allowed kernel types are \code{uniform} (uniform kernel), \code{triangular} (triangular kernel) and \code{epan} (Epanechnikov kernel). Default is \code{uniform}.}

\item{approx}{forces the command to conduct the covariate balance tests using a large-sample approximation instead of finite-sample exact randomization inference methods.}

\item{level}{the minimum accepted value of the p-value from the covariate balance tests. Default is .15.}

\item{reps}{number of replications. Default is 1000.}

\item{seed}{the seed to be used for the randomization tests.}

\item{plot}{draws a scatter plot of the minimum p-value from the covariate balance test against window length.}

\item{quietly}{suppress output}

\item{obsstep}{the minimum number of observations to be added on each side of the cutoff for the sequence of fixed-increment nested windows. This option is deprecated and only included for backward compatibility.}
}
\value{
\item{window}{recommended window (NA is covariates are not specified)}
\item{wlist}{list of window lengths}
\item{results}{table including window lengths, minimum p-value in each window, corresponding number of the variable with minimum p-value (i.e. column of covariate matrix), Binomial test p-value and sample sizes to the left and right of the cutoff in each window.}
\item{summary}{summary statistics.}
}
\description{
\code{rdwinselect} implements the window-selection procedure
 based on balance tests for RD designs under local randomization.
 Specifically, it constructs a sequence of nested windows around the RD cutoff
 and reports binomial tests for the running variable runvar and covariate balance
 tests for covariates covariates (if specified). The recommended window is the
 largest window around the cutoff such that the minimum p-value of the balance test
 is larger than a prespecified level for all nested (smaller) windows. By default,
 the p-values are calculated using randomization inference methods.
}
\examples{
# Toy dataset
X <- array(rnorm(200),dim=c(100,2))
R <- X[1,] + X[2,] + rnorm(100)
# Window selection adding 5 observations at each step
# Note: low number of replications to speed up process.
tmp <- rdwinselect(R,X,obsmin=10,wobs=5,reps=500)
# Window selection setting initial window and step
# The user should increase the number of replications.
tmp <- rdwinselect(R,X,wmin=.5,wstep=.125,reps=500)
# Window selection with approximate (large sample) inference and p-value plot
tmp <- rdwinselect(R,X,wmin=.5,wstep=.125,approx=TRUE,nwin=80,quietly=TRUE,plot=TRUE)


}
\references{
Cattaneo, M.D., R. Titiunik and G. Vazquez-Bare. (2016). \href{https://rdpackages.github.io/references/Cattaneo-Titiunik-VazquezBare_2016_Stata.pdf}{Inference in Regression Discontinuity Designs under Local Randomization}. \emph{Stata Journal} 16(2): 331-367.
}
\author{
Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}

Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}

Gonzalo Vazquez-Bare, UC Santa Barbara. \email{gvazquez@econ.ucsb.edu}
}
