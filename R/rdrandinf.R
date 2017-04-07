
###################################################################
# rdrandinf: randomization inference in RD window
# !version 0.1 4-Apr-2017
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################

#' Randomization Inference for RD Designs under Local Randomization
#'
#' \code{rdrandinf} implements randomization inference and related methods for RD designs,
#' using observations in a specified or data-driven selected window around the cutoff where
#' local randomization is assumed to hold.
#'
#'
#' @author
#' Matias Cattaneo, University of Michigan. \email{cattaneo@umich.edu}
#'
#' Rocio Titiunik, University of Michigan. \email{titiunik@umich.edu}
#'
#' Gonzalo Vazquez-Bare, University of Michigan. \email{gvazquez@umich.edu}
#'
#' @references
#' M.D. Cattaneo, B. Frandsen and R. Titiunik. (2015).  \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to Party Advantages in the U.S. Senate.} \emph{Journal of Causal Inference} 3(1): 1-24.
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2016). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Titiunik-VazquezBare_2016_Stata.pdf}{Inference in Regression Discontinuity Designs under Local Randomization.} \emph{Stata Journal} 16(2): 331-367.
#'
#' M.D. Cattaneo, R. Titiunik and G. Vazquez-Bare. (2017). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Titiunik-VazquezBare_2017_JPAM.pdf}{Comparing Inference Approaches for RD Designs: A Reexamination of the Effect of Head Start on Child Mortality.} \emph{Journal of Policy Analysis and Management}, forthcoming.
#'
#'
#'
#'
#' @param Y a vector containing the values of the outcome variable.
#' @param R a vector containing the values of the running variable.
#' @param cutoff the RD cutoff (default is 0).
#' @param wl the left limit of the window. The default takes the minimum of the running variable.
#' @param wr the right limit of the window. The default takes the maximum of the running variable.
#' @param reps the number of replications (default is 1000).
#' @param statistic the statistic to be used in randomization inference. Options are \code{ttest} (difference in means), \code{ksmirnov} (Kolmogorov-Smirnov statistic), \code{ranksum} (Wilcoxon-Mann-Whitney standardized statistic) or \code{all}, which gives all three statistics. Default option is \code{ttest}.
#' @param p the order of the polynomial for outcome transformation model (default is 0).
#' @param nulltau the value of the treatment effect under the null hypothesis (default is 0).
#' @param evall the point at the left of the cutoff at which to evaluate the transformed outcome is evaluated. Default is the cutoff value.
#' @param evalr specifies the point at the right of the cutoff at which the transformed outcome is evaluated. Default is the cutoff value.
#' @param kernel specifies the type of kernel to use as weighting scheme. Allowed kernel types are \code{uniform} (uniform kernel), \code{triangular} (triangular kernel) and \code{epan} (Epanechnikov kernel). Default is \code{uniform}.
#' @param ci calculates a confidence interval for the treatment effect by test inversion. \code{ci} can be specified as a scalar or a vector, where the first element indicates the level of the confidence interval and the remaining elements, if specified, indicate the grid of treatment effects to be evaluated. This option uses \code{rdsensitivity} to calculate the confidence interval. See corresponding help for details.
#' @param interfci the level for Rosenbaum's confidence interval under arbitrary interference between units.
#' @param seed the seed to be used for the randomization test.
#' @param fuzzy indicates that the RD design is fuzzy. \code{fuzzy} can be specified as a vector containing the values of the endogenous treatment variable, or as a list where the first element is the vector of endogenous treatment values and the second element is a string containing the name of the statistic to be used. Allowed statistics are \code{ar} (Anderson-Rubin statistic) and \code{tsls} (2SLS statistic). Default statistic is \code{ar}. The \code{tsls} statistic relies on large-sample approximation.
#' @param d the effect size for asymptotic power calculation. Default is 0.5 * standard deviation of outcome variable for the control group.
#' @param dscale the fraction of the standard deviation of the outcome variable for the control group used as alternative hypothesis for asymptotic power calculation. Default is 0.5.
#' @param bernoulli the probabilities of treatment for each unit when assignment mechanism is a Bernoulli trial. This option should be specified as a vector of length equal to the length of the outcome and running variables.
#' @param quietly suppresses the output table.
#' @param covariates the covariates used by \code{rdwinselect} to choose the window when \code{wl} and \code{wr} are not specified. This should be a matrix of size n x k where n is the total sample size and k is the number of covariates.
#' @param obsmin the minimum number of observations above and below the cutoff in the smallest window employed by the companion command \code{rdwinselect}. Default is 10.
#' @param obsstep the minimum number of observations to be added on each side of the cutoff for the sequence of nested windows constructed by the companion command \code{rdwinselect}. Default is 2.
#' @param wmin the smallest window to be used (if \code{minobs} is not specified) by the companion command \code{rdwinselect}. Specifying both \code{wmin} and \code{obsmin} returns an error.
#' @param wstep the increment in window length (if \code{obsstep} is not specified) by the companion command \code{rdwinselect}.  Specifying both \code{obsstep} and \code{wstep} returns an error.
#' @param nwindows the number of windows to be used by the companion command \code{rdwinselect}. Default is 10.
#' @param rdwstat the statistic to be used by the companion command \code{rdwinselect} (see corresponding help for options). Default option is \code{ttest}.
#' @param approx forces the companion command \code{rdwinselect} to conduct the covariate balance tests using a large-sample approximation instead of finite-sample exact randomization inference methods.
#' @param rdwreps the number of replications to be used by the companion command \code{rdwinselect}. Default is 1000.
#' @param level the minimum accepted value of the p-value from the covariate balance tests to be used by the companion command \code{rdwinselect}. Default is .15.
#' @param plot draws a scatter plot of the minimum p-value from the covariate balance test against window length implemented by the companion command \code{rdwinselect}.
#'
#' @return
#' \item{sumstats}{summary statistics}
#' \item{obs.stat}{observed statistic(s)}
#' \item{p.value}{randomization p-value(s)}
#' \item{asy.pvalue}{asymptotic p-value(s)}
#' \item{window}{chosen window}
#' \item{ci}{confidence interval (only if \code{ci} option is specified)}
#' \item{interf.ci}{confidence interval under interferecen (only if \code{interfci} is specified)}
#'
#' @examples
#' # Toy dataset
#' X <- array(rnorm(200),dim=c(100,2))
#' R <- X[1,] + X[2,] + rnorm(100)
#' Y <- 1 + R -.5*R^2 + .3*R^3 + (R>=0) + rnorm(100)
#' # Randomization inference in window (-.75,.75)
#' tmp <- rdrandinf(Y,R,wl=-.75,wr=.75)
#' # Randomization inference in window (-.75,.75), all statistics
#' tmp <- rdrandinf(Y,R,wl=-.75,wr=.75,statistic='all')
#' # Randomization inference with window selection
#' # Note: low number of replications to speed up process.
#' # The user should increase the number of replications.
#' tmp <- rdrandinf(Y,R,statistic='all',covariates=X,wmin=.5,wstep=.125,rdwreps=500)
#'
#'
#'
#' @export


rdrandinf = function(Y,R,
                     cutoff = 0,
                     wl = '',
                     wr = '',
                     reps = 1000,
                     statistic = 'ttest',
                     p = 0,
                     nulltau = 0,
                     evall = '',
                     evalr = '',
                     kernel = 'uniform',
                     ci,
                     interfci = '',
                     seed = '',
                     fuzzy = '',
                     d = '',
                     dscale = '',
                     bernoulli,
                     quietly = FALSE,

                     covariates,
                     obsmin = '',
                     obsstep = '',
                     wmin = '',
                     wstep = '',
                     nwindows = 10,
                     rdwstat = 'ttest',
                     approx = FALSE,
                     rdwreps = 1000,
                     level = .15,
                     plot = FALSE){


  #################################################################
  # Parameters and error checking
  #################################################################

  randmech = 'fixed margins'

  Rc.long = R - cutoff

  if(missing(bernoulli)){
    data = cbind(Y,R)
    data = data[complete.cases(data),]
    Y = data[,1]
    R = data[,2]
  } else{
    data = cbind(Y,R,bernoulli)
    data = data[complete.cases(data),]
    Y = data[,1]
    R = data[,2]
    bernoulli = data[,3]
  }

  if (cutoff<=min(R,na.rm=TRUE) | cutoff>=max(R,na.rm=TRUE)){stop('Cutoff must be within the range of the running variable')}
  if (p<0){stop('p must be a positive integer')}
  if (statistic!='ttest' & statistic!='ksmirnov' & statistic!='ranksum' & statistic!='all'){stop(paste(statistic,'not a valid statistic'))}
  if (kernel!='uniform' & kernel!='triangular' & kernel!='epan'){stop(paste(kernel,'not a valid kernel'))}
  if (kernel!='uniform' & evall!='' & evalr!=''){
    if (evalr!=cutoff | evalr !=cutoff) {stop('kernel only allowed when evall=evalr=cutoff')}
  }
  if (kernel!='uniform' & statistic!='ttest'){stop('kernel only allowed for ttest')}
  if (!missing(ci)){if (ci[1]>1 | ci[1]<0){stop('ci must be in [0,1]')}}
  if (interfci!=''){
    if (interfci>1 | interfci<0){stop('interfci must be in [0,1]')}
    if (statistic!='ttest' & statistic!='ksmirnov' & statistic!='ranksum'){stop('interfci only allowed with ttest, ksmirnov or ranksum')}
  }
  if (!missing(bernoulli)){
    randmech = 'Bernoulli'
    if (max(bernoulli,na.rm=TRUE)>1 | min(bernoulli,na.rm=TRUE)<0){stop('bernoulli probabilities must be in [0,1]')}
    if (length(bernoulli)!=length(R)){stop('bernoulli should have the same length as the running variable')}
  }
  if (wl!='' & wr!=''){
    wselect = 'set by user'
    if (wl>=wr){stop('wl has to be smaller than wr')}
    if (wl>=cutoff | wr<=cutoff){stop('window does not include cutoff')}
  }
  if (wl=='' & wr!=''){stop('wl not specified')}
  if (wl!='' & wr==''){stop('wr not specified')}
  if (evall!='' & evalr==''){stop('evalr not specified')}
  if (evall=='' & evalr!=''){stop('evall not specified')}
  if (d!='' & dscale!=''){stop('cannot specify both d and dscale')}

  Rc = R - cutoff
  D = as.numeric(Rc >= 0)

  n = length(D)
  n1 = sum(D)
  n0 = n - n1

  if (seed!=''){set.seed(seed)}

  if (fuzzy!=''){
    if (class(fuzzy)!='list'){
      fuzzy.stat = 'ar'
      fuzzy.tr = fuzzy
    } else {
      fuzzy.tr = fuzzy[[1]]
      if (fuzzy[[2]]!='ar' & fuzzy[[2]]!='tsls'){
        stop('fuzzy statistic not valid')
      } else {
        fuzzy.stat = 'wald'
      }
    }
  } else {
    fuzzy.stat = ''
  }


  #################################################################
  # Window selection
  #################################################################

  if (wl=='' & wr==''){
    if (missing(covariates)){
      wl = min(R,na.rm=TRUE)
      wr = max(R,na.rm=TRUE)
      wselect = 'run. var. range'
    } else {
      wselect = 'rdwinselect'
      if (quietly==FALSE){cat('\nRunning rdwinselect...\n')}
      rdwlength = rdwinselect(Rc.long,covariates,obsmin=obsmin,obsstep=obsstep,wmin=wmin,wstep=wstep,
                  nwindows=nwindows,statistic=rdwstat,approx=approx,
                  reps=rdwreps,plot=plot,level=level,quietly=TRUE)
      wl = rdwlength$window[1]
      wr = rdwlength$window[2]
      if (quietly==FALSE){cat('\nrdwinselect complete.\n')}
    }
  }
  if (quietly==FALSE){cat(paste0('\nSelected window = [',wl,';',wr,'] \n'))}


  if (evall!=''&evalr!=''){if (evall<wl | evalr>wr){stop('evall and evalr need to be inside window')}}

  ww = (Rc >= wl) & (Rc <= wr)

  Yw = Y[ww]
  Rw = Rc[ww]
  Dw = D[ww]

  if (fuzzy!=''){
      Tw = fuzzy.tr[ww]
  }

  if (missing(bernoulli)){
    data = cbind(Yw,Rw,Dw)
    data = data[complete.cases(data),]
    Yw = data[,1]
    Rw = data[,2]
    Dw = data[,3]
  } else {
    Bew = bernoulli[ww]
    data = cbind(Yw,Rw,Dw,Bew)
    data = data[complete.cases(data),]
    Yw = data[,1]
    Rw = data[,2]
    Dw = data[,3]
    Bew = data[,4]
  }


  n.w = length(Dw)
  n1.w = sum(Dw)
  n0.w = n.w - n1.w


  #################################################################
  # Summary statistics
  #################################################################

  sumstats = array(NA,dim=c(5,2))
  sumstats[1,] = c(n0,n1)
  sumstats[2,] = c(n0.w,n1.w)
  mean0 = mean(Yw[Dw==0],na.rm=TRUE)
  mean1 = mean(Yw[Dw==1],na.rm=TRUE)
  sd0 = sd(Yw[Dw==0],na.rm=TRUE)
  sd1 = sd(Yw[Dw==1],na.rm=TRUE)
  sumstats[3,] = c(mean0,mean1)
  sumstats[4,] = c(sd0,sd1)
  sumstats[5,] = c(wl,wr)

  if (d=='' & dscale==''){
    delta = .5*sd0
  }
  if (d!='' & dscale==''){
    delta = d
  }
  if (d=='' & dscale!=''){
    delta = dscale*sd0
  }

  #################################################################
  # Weights
  #################################################################

  kweights = rep(1,n.w)

  if (kernel=='triangular'){
    bwt = wr - cutoff
    bwc = wl - cutoff
    kweights[Dw==1] = (1-abs(Rw[Dw==1]/bwt))*(abs(Rw[Dw==1]/bwt)<1)
    kweights[Dw==0] = (1-abs(Rw[Dw==0]/bwc))*(abs(Rw[Dw==0]/bwc)<1)
  }
  if (kernel=='epan'){
    bwt = wr - cutoff
    bwc = wl - cutoff
    kweights[Dw==1] = .75*(1-(Rw[Dw==1]/bwt)^2)*(abs(Rw[Dw==1]/bwt)<1)
    kweights[Dw==0] = .75*(1-(Rw[Dw==0]/bwc)^2)*(abs(Rw[Dw==0]/bwc)<1)
  }


  #################################################################
  # Outcome adjustment: model and null hypothesis
  #################################################################

  Y.adj = Yw

  if (p>0){
    if (evall=='' & evalr==''){
      evall = cutoff
      evalr = cutoff
    }
    R.adj = Rw + cutoff - Dw*evalr - (1-Dw)*evall
    Rpoly = poly(R.adj,order=p,raw=TRUE)
    lfit.t = lm(Yw[Dw==1] ~ Rpoly[Dw==1,],weights=kweights[Dw==1])
    Y.adj[Dw==1] = lfit.t$residuals + lfit.t$coefficients[1]
    lfit.c = lm(Yw[Dw==0] ~ Rpoly[Dw==0,],weights=kweights[Dw==0])
    Y.adj[Dw==0] = lfit.c$residuals + lfit.c$coefficients[1]
  }

  Y.adj.null = Y.adj - nulltau*Dw


  #################################################################
  # Observed statistics and asymptotic p-values
  #################################################################


  if (fuzzy==''){
    results = rdrandinf.model(Y.adj.null,Dw,statistic=statistic,pvalue=TRUE,kweights=kweights,delta=delta)
  } else {
    results = rdrandinf.model(Y.adj.null,Dw,statistic=fuzzy.stat,endogtr=Tw,pvalue=TRUE,kweights=kweights,delta=delta)
  }

  obs.stat = as.numeric(results$statistic)

  if (p==0){
    if (fuzzy.stat=='wald'){
      aux = AER::ivreg(Yw ~ Tw,~ Dw,weights=kweights)
      obs.stat = aux$coefficients["Tw"]
      se = sqrt(diag(sandwich::vcovHC(aux,type='HC1'))['Tw'])
      tstat = aux$coefficients['Tw']/se
      asy.pval = as.numeric(2*pnorm(-abs(tstat)))
      asy.power = as.numeric(1-pnorm(1.96-delta/se)+pnorm(-1.96-delta/se))
    } else {
      asy.pval = as.numeric(results$p.value)
      asy.power = as.numeric(results$asy.power)
    }

  } else {
    if (statistic=='ttest'|statistic=='all'){
      lfit = lm(Yw ~ Dw + Rpoly + Dw*Rpoly,weights=kweights)
      se = sqrt(diag(sandwich::vcovHC(lfit,type='HC2'))['Dw'])
      tstat = lfit$coefficients['Dw']/se
      asy.pval = as.numeric(2*pnorm(-abs(tstat)))
      asy.power = as.numeric(1-pnorm(1.96-delta/se)+pnorm(-1.96-delta/se))
    }
    if (statistic=='ksmirnov'|statistic=='ranksum'){
      asy.pval = NA
      asy.power = NA
    }
    if (statistic=='all'){
      asy.pval = c(as.numeric(asy.pval),NA,NA)
      asy.power = c(as.numeric(asy.power),NA,NA)
    }

    if (fuzzy.stat=='wald'){
      inter = Rpoly*Dw
      aux = AER::ivreg(Yw ~ Rpoly + inter + Tw,~ Rpoly + inter + Dw,weights=kweights)
      obs.stat = aux$coefficients["Tw"]
      se = sqrt(diag(sandwich::vcovHC(aux,type='HC1'))['Tw'])
      tstat = aux$coefficients['Tw']/se
      asy.pval = as.numeric(2*pnorm(-abs(tstat)))
      asy.power = as.numeric(1-pnorm(1.96-delta/se)+pnorm(-1.96-delta/se))
    }
  }


  #################################################################
  # Randomization-based inference
  #################################################################


  if (statistic == 'all'){
    stats.distr = array(NA,dim=c(reps,3))
  } else{
    stats.distr = array(NA,dim=c(reps,1))
  }

  if (quietly==FALSE){cat('\nRunning randomization-based test...\n')}

  if (fuzzy.stat!='wald'){
    if (missing(bernoulli)){

      max.reps = choose(n.w,n1.w)
      reps = min(reps,max.reps)
      if (max.reps<reps){
        warning(paste0('Chosen no. of reps > total no. of permutations.\n reps set to ',reps,'.'))
      }

      for (i in 1:reps) {
        D.sample = sample(Dw,replace=FALSE)
        if (missing(fuzzy) |  fuzzy==''){
          obs.stat.sample = as.numeric(rdrandinf.model(Y.adj.null,D.sample,statistic,kweights=kweights,delta=delta)$statistic)
        } else {
          obs.stat.sample = as.numeric(rdrandinf.model(Y.adj.null,D.sample,statistic=fuzzy.stat,endogtr=Tw,kweights=kweights,delta=delta)$statistic)
        }
        stats.distr[i,] = obs.stat.sample
      }

    } else {

      for (i in 1:reps) {
        D.sample = as.numeric(runif(n.w)<=Bew)
        if (mean(D.sample)==1 | mean(D.sample)==0){
          stats.distr[i,] = NA # ignore cases where bernoulli assignment mechanism gives no treated or no controls
        } else {
          obs.stat.sample = as.numeric(rdrandinf.model(Y.adj.null,D.sample,statistic,kweights=kweights,delta=delta)$statistic)
          stats.distr[i,] = obs.stat.sample
        }
      }

    }

    if(quietly==FALSE){cat('\nRandomization-based test complete. \n')}

    if (statistic == 'all'){
      p.value1 = mean(abs(stats.distr[,1]) >= abs(obs.stat[1]),na.rm=TRUE)
      p.value2 = mean(abs(stats.distr[,2]) >= abs(obs.stat[2]),na.rm=TRUE)
      p.value3 = mean(abs(stats.distr[,3]) >= abs(obs.stat[3]),na.rm=TRUE)
      p.value = c(p.value1,p.value2,p.value3)
    } else{
      p.value = mean(abs(stats.distr) >= abs(obs.stat),na.rm=TRUE)
    }
  } else {
    p.value = NA
  }

  #################################################################
  # Confidence interval
  #################################################################

  if (!missing(ci)){
    ci.level = ci[1]
    if (length(ci)>1){
      tlist = ci[-1]
      aux = rdsensitivity(Y,Rc,wlist=wr,tlist=tlist,ci=c(wr,ci.level),
                          reps=reps,nodraw=TRUE)

    } else {
      aux = rdsensitivity(Y,Rc,wlist=wr,ci=c(wr,ci.level),
                          reps=reps,nodraw=TRUE)

    }
    conf.int = aux$ci
  }


  #################################################################
  # Confidence interval under interference
  #################################################################

  if (interfci!=''){
    p.low = interfci/2
    p.high = 1-interfci/2
    qq = quantile(stats.distr,probs=c(p.low,p.high))
    interf.ci = c(obs.stat-as.numeric(qq[2]),obs.stat-as.numeric(qq[1]))
  }


  #################################################################
  # Output and display results
  #################################################################


  if (missing(ci) & interfci==''){
    output = list(sumstats = sumstats,
                  obs.stat = obs.stat,
                  p.value = p.value,
                  asy.pvalue = asy.pval,
                  window = c(wl,wr))
  }
  if (!missing(ci) & interfci==''){
    output = list(sumstats = sumstats,
                  obs.stat = obs.stat,
                  p.value = p.value,
                  asy.pvalue = asy.pval,
                  window = c(wl,wr),
                  ci = conf.int)
  }
  if (missing(ci) & interfci!=''){
    output = list(sumstats = sumstats,
                  obs.stat = obs.stat,
                  p.value = p.value,
                  asy.pvalue = asy.pval,
                  window = c(wl,wr),
                  interf.ci = interf.ci)
  }
  if (!missing(ci) & interfci!=''){
    output = list(sumstats = sumstats,
                  obs.stat = obs.stat,
                  p.value = p.value,
                  asy.pvalue = asy.pval,
                  window = c(wl,wr),
                  ci = conf.int,
                  interf.ci = interf.ci)
  }


  if (quietly==FALSE){
    if (statistic=='ttest'){statdisp = 'Diff. in means'}
    if (statistic=='ksmirnov'){statdisp = 'Kolmogorov-Smirnov'}
    if (statistic=='ranksum'){statdisp = 'Rank sum z-stat'}
    if (fuzzy.stat=='ar'){
      statdisp = 'Anderson-Rubin'
      obs.stat = NA
      output[[2]] = NA
    }
    if (fuzzy.stat=='wald'){statdisp = 'TSLS'}

    cat('\n\n')
    cat(paste0(format('Number of obs =', width=22), toString(n))); cat("\n")
    cat(paste0(format('Order of poly =', width=22), p)); cat("\n")
    cat(paste0(format('Kernel type   =', width=22), kernel)); cat("\n")
    cat(paste0(format('Reps          =', width=22), reps)); cat("\n")
    cat(paste0(format('Window        =', width=22), wselect)); cat("\n")
    cat(paste0(format('H0:       tau =', width=22), nulltau)); cat("\n")
    cat(paste0(format('Randomization =', width=22), randmech))
    cat('\n\n')

    cat(paste0(format(paste0("Cutoff c = ", toString(round(cutoff, 3))), width=22), format("Left of c", width=16), format("Right of c", width=16))); cat("\n")
    cat(paste0(format("Number of obs",       width=22), format(toString(n0),             width=16), format(toString(n1),             width=16))); cat("\n")
    cat(paste0(format("Eff. number of obs",  width=22), format(toString(n0.w),           width=16), format(toString(n1.w),           width=16))); cat("\n")
    cat(paste0(format("Mean of outcome",     width=22), format(toString(round(mean0,3)), width=16), format(toString(round(mean1,3)), width=16))); cat("\n")
    cat(paste0(format("S.d. of outcome",     width=22), format(toString(round(sd0,3)),   width=16), format(toString(round(sd1,3)),   width=16))); cat("\n")
    cat(paste0(format("Window",              width=22), format(toString(round(wl,3)),    width=16), format(toString(round(wr,3)),    width=16)))
    cat("\n\n")

    cat(paste0(format('',width=34),format('Finite sample',width=21),format('Large sample')));cat('\n\n')

    cat(paste0(format('Statistic',    width=21),
               format('T',            width=15),
               format('P>|T|',        width=15),
               format('P>|T|',        width=10),
               format('Power vs d = ',width=3),
               format(toString(round(delta,3))))); cat('\n\n')

    if (statistic!='all'){

      cat(paste0(format(statdisp,width=21),
                 format(toString(round(obs.stat,3)),width=15),
                 format(toString(round(p.value,3)), width=15),
                 format(toString(round(asy.pval,3)),width=10),
                 format(toString(round(asy.power,3))))); cat('\n\n')

    }

    if (statistic=='all'){

      cat(paste0(format('Diff. in means',width=21),
                 format(toString(round(obs.stat[1],3)),width=15),
                 format(toString(round(p.value[1],3)), width=15),
                 format(toString(round(asy.pval[1],3)),width=10),
                 format(toString(round(asy.power[1],3))))); cat('\n')

      cat(paste0(format('Kolmogorov-Smirnov',width=21),
                 format(toString(round(obs.stat[2],3)),width=15),
                 format(toString(round(p.value[2],3)), width=15),
                 format(toString(round(asy.pval[2],3)),width=10),
                 format(toString(round(asy.power[2],3))))); cat('\n')

      cat(paste0(format('Rank sum z-stat',width=21),
                 format(toString(round(obs.stat[3],3)),width=15),
                 format(toString(round(p.value[3],3)), width=15),
                 format(toString(round(asy.pval[3],3)),width=10),
                 format(toString(round(asy.power[3],3))))); cat('\n\n')

    }

    if (!missing(ci)){
      cat('\n')
      cat(paste0((1-ci.level)*100,'% confidence interval: [',conf.int[1],',',conf.int[2],']\n'))
    }

    if (interfci!=''){
      cat('\n')
      cat(paste0((1-interfci)*100,'% confidence interval under interference: [',round(interf.ci[1],3),';',round(interf.ci[2],3),']')); cat('\n')
    }
  }
  return(output)
}