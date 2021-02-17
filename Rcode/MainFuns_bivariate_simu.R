library(splines)
library(copula)
library(cubature) ###multivariate integration

###############################################################
################  Estimate the Marginal CDF distirbution
###############################################################
MarginCDF <- function(Y, ID, Age, X, numKnot, deg, Normal)
{
   ## estimate the marginal CDF
   Age_lower <- min(Age)
   Age_upper <- max(Age)
   KN0 <- seq(from=Age_lower, to=Age_upper, length=numKnot[1])[-c(1,numKnot[1])]
   KN1 <- seq(from=Age_lower, to=Age_upper, length=numKnot[2])[-c(1,numKnot[2])]
   KN2 <- seq(from=Age_lower, to=Age_upper, length=numKnot[3])[-c(1,numKnot[3])]
   T0_bs <- bs(x=Age, knots=KN0, degree=deg, intercept=T)
   T1_bs <- bs(x=Age, knots=KN1, degree=deg, intercept=T)
   T2_bs <- bs(x=Age, knots=KN2, degree=deg, intercept=T)
   X0_bs <- T0_bs  ## the intercept term
   X1_bs <- matrix(X[,1], nrow=length(X[,1]), ncol=ncol(T1_bs))*T1_bs
   X2_bs <- matrix(X[,2], nrow=length(X[,2]), ncol=ncol(T2_bs))*T2_bs
  ######### Fit the WLS estimate to the the coefficient for the BS, 
  ###### then with r_ls* BS(ls), we can obtain  point estimate of the coefficient
   Ni <- table(ID)
   Nsub <- length(unique(ID))
   wi <- 1/(Nsub*rep(Ni, Ni))
   ### outcome=sbp
   fit_WLS <- lm(Y~0+X0_bs+X1_bs+X2_bs, weights=wi)
   p0 <- ncol(X0_bs); p1 <- ncol(X1_bs); p2 <- ncol(X2_bs)
   gamm0 <- fit_WLS$coef[1:p0]
   gamm1 <- fit_WLS$coef[(p0+1):(p0+p1)]
   gamm2 <- fit_WLS$coef[(p0+p1+1):(p0+p1+p2)]
   gamm_est <- list(gamm0, gamm1, gamm2)
   #### estimators of the coefficient
   bet0 <- T0_bs%*%matrix(gamm0, ncol=1)
   bet1 <- T1_bs%*%matrix(gamm1, ncol=1)
   bet2 <- T2_bs%*%matrix(gamm2, ncol=1)
   bet_est <- cbind(bet0, bet1, bet2)
   ############## estimate the variance parameter (based on MLE estimate)
   resError <- Y-bet0-X[,1]*bet1-X[,2]*bet2
   sigma2Est <- sum(wi*(resError^2))/sum(wi)
   ###### The obtained density and CDF estimator
   temp_mu <- bet0+X[,1]*bet1+X[,2]*bet2
   var_est <- sigma2Est
   #### for the marginal normal/lognormal distribution (if normal=F, Y should be log(Y))
   density_est <- dnorm(Y-temp_mu, mean=0, sd=sqrt(var_est))
   cdf_est <- pnorm(Y-temp_mu, mean=0, sd=sqrt(var_est))
   list(gamm_est=gamm_est, bet_est=bet_est, var_est=var_est, density_est=density_est, cdf_est=cdf_est)
}

#################################################################################################
################  Estimate the bivariate copula parameter for one outcome at two time points 
#################################################################################################

################################## For the same outcome at two time points
#####################################################################################################################
### first define the loglikehood of the bivariate coupla function
############################# for the spline function with one inner knot and quadric spline (number of parameter=10) 
BiCop_EqualOut_LogLik <- function(lambda1, lambda2, lambda3, basisFun, cdfMargin, ID)
{
    tempPara <- c(lambda1, lambda2, lambda3)
    tempMat <- matrix(0, nrow=ncol(basisFun), ncol=ncol(basisFun))
    tempMat[upper.tri(tempMat, diag=T)] <- tempPara
    tempMat <- tempMat + t(tempMat) - diag(diag(tempMat))
    lamPara <- c(tempMat)
    uniID <- unique(ID)
    Nsub <- length(uniID)
    Ni <- table(ID)
    dd <- cumsum(c(0,Ni))
    SamLoc <- cbind(dd[1:(length(dd)-1)], dd[2:length(dd)])
    Samlist <- list()
    for(kk in 1:Nsub)  Samlist[[kk]] <- (SamLoc[kk,1]+1):SamLoc[kk,2]
    logCopFun <- function(paraCop)
    {
      tempCop <- normalCopula(param=paraCop[3], dim=2, dispstr="un") 
      return(log(dCopula(u=c(paraCop[1], paraCop[2]), copula=tempCop)))
    }
    loglikSam <- function(tempLoc, basisFun, cdfMargin, lamPara)
    {
      cdfSam <- cdfMargin[tempLoc] 
      tempBasis <- basisFun[tempLoc, ]
      NumRepeat <- length(tempLoc)
      AllCom0 <- expand.grid(1:NumRepeat,1:NumRepeat)
      selectLoc <- which(AllCom0[,1]<AllCom0[,2])
      AllCom <- AllCom0[selectLoc,]
      NumCom <- nrow(AllCom)
      cdf1 <- cdfSam[AllCom[,1]]
      cdf2 <- cdfSam[AllCom[,2]]
      temp.delta <- apply(AllCom, MARGIN=1, function(u) sum(kronecker(tempBasis[u[1],], tempBasis[u[2],])*lamPara))
      if(sum(abs(temp.delta)>=1)>0)  wLogLik <- -1e200
      else
      {
        wLogLik <- 1/(NumCom)*sum(apply(cbind(cdf1, cdf2, temp.delta), MARGIN=1, logCopFun))
      }
      return(wLogLik)
    }
    Res_logLik <- do.call(sum, lapply(Samlist, loglikSam, basisFun=basisFun, cdfMargin=cdfMargin, lamPara=lamPara))
    return(Res_logLik)
}

#######################################################################################
BiCop_EqualOut_paraEst <- function(Y, ID, Age, X, numKnot, deg, Normal)
{
  tempEst <- MarginCDF(Y=Y, ID=ID, Age=Age, X=X, numKnot=numKnot, deg=deg, Normal=Normal)
  gammEst <- tempEst$gamm_est 
  gamm0 <- gammEst[[1]]
  gamm1 <- gammEst[[2]]
  gamm2 <- gammEst[[3]]
  densityEst <- tempEst$density_est
  cdfEst <- tempEst$cdf_est
  varEst <- tempEst$var_est
  #######
  Age_lower <- min(Age)
  Age_upper <- max(Age)
  KN0 <- seq(from=Age_lower, to=Age_upper, length=numKnot[1])[-c(1,numKnot[1])]
  #### use the 2-degree spline to estimate
  bs_Fun <- bs(x=Age, knots=KN0, degree=deg, intercept=T)
  para_dim <- ncol(bs_Fun)*(ncol(bs_Fun)+1)/2
  
  #### the iteration process (the initial value is 10-dimensional)
  likDiff <- 999
  iter <- 0
  lam1_old <- rep(0.5,3);lam2_old <- rep(0.5,3); lam3_old <- rep(0.5,4)
  while(likDiff>1e-3)
  {
    loglik_old <- BiCop_EqualOut_LogLik(lambda1=lam1_old, lambda2=lam2_old, lambda3=lam3_old, 
                                        basisFun=bs_Fun, cdfMargin=cdfEst, ID=ID)
    OptimEst1 <- optim(par=lam1_old, fn=BiCop_EqualOut_LogLik, gr=NULL, control=list(fnscale=-1, reltol=0.001), lambda2=lam2_old,
                       lambda3=lam3_old, basisFun=bs_Fun, cdfMargin=cdfEst, ID=ID)
    lam1_hat <- OptimEst1$par
    OptimEst2 <- optim(par=lam2_old, fn=BiCop_EqualOut_LogLik, gr=NULL, control=list(fnscale=-1, reltol=0.001), lambda1=lam1_hat,
                       lambda3=lam3_old, basisFun=bs_Fun, cdfMargin=cdfEst, ID=ID)
    lam2_hat <- OptimEst2$par
    OptimEst3 <- optim(par=lam3_old, fn=BiCop_EqualOut_LogLik, gr=NULL, control=list(fnscale=-1, reltol=0.001), lambda1=lam1_hat,
                       lambda2=lam2_hat, basisFun=bs_Fun, cdfMargin=cdfEst, ID=ID)
    lam3_hat <- OptimEst3$par
    ##### calculate the log-likelihood at each iteration
    loglik_hat <- BiCop_EqualOut_LogLik(lambda1=lam1_hat, lambda2=lam2_hat, lambda3=lam3_hat,
                                        basisFun=bs_Fun, cdfMargin=cdfEst, ID=ID)
    likDiff <- abs(loglik_hat-loglik_old)
    lam1_old <- lam1_hat; lam2_old <- lam2_hat; lam3_old <- lam3_hat
    iter <- iter + 1
  }
  temp_lamEst <- c(lam1_hat, lam2_hat, lam3_hat)
  paraMat <- matrix(0, nrow=ncol(bs_Fun), ncol=ncol(bs_Fun))
  paraMat[upper.tri(paraMat, diag=T)] <- temp_lamEst
  paraMat <- paraMat + t(paraMat) - diag(diag(paraMat))
  lamEst <- c(paraMat)
  return(matrix(lamEst, nrow=1))
}

#######################################################################################
######## calculate the estimated bivariate conditional density of Y1 and Y2 given tt, xx
BiCop_EqualOut_disEst <- function(yy1, yy2, gammList, varEst, copulaPara, tt1, tt2, xx, AgeGrid, numKnot, deg, Normal)
{
    ### the parameter beta
    temp_gamma0 <- gammList[[1]]
    temp_gamma1 <- gammList[[2]]
    temp_gamma2 <- gammList[[3]]
    #####
    Age_lower <- min(AgeGrid)
    Age_upper <- max(AgeGrid)
    KN0 <- seq(from=Age_lower, to=Age_upper, length=numKnot[1])[-c(1,numKnot[1])]
    KN1 <- seq(from=Age_lower, to=Age_upper, length=numKnot[2])[-c(1,numKnot[2])]
    KN2 <- seq(from=Age_lower, to=Age_upper, length=numKnot[3])[-c(1,numKnot[3])]
    Bs0 <- bs(AgeGrid, knots=KN0, degree=deg, intercept=T)
    Bs1 <- bs(AgeGrid, knots=KN1, degree=deg, intercept=T)
    Bs2 <- bs(AgeGrid, knots=KN2, degree=deg, intercept=T)
    #############################################################################
    ##### estimate the marginal density and cdf for tt1
    loc_Spe1 <- which(AgeGrid==tt1)
    Bs0_Spe1 <- Bs0[loc_Spe1,]
    Bs1_Spe1 <- Bs1[loc_Spe1,]
    Bs2_Spe1 <- Bs2[loc_Spe1,]
    bet0_time1 <- sum(Bs0_Spe1*temp_gamma0)  ## the intercept
    bet1_time1 <- sum(Bs1_Spe1*temp_gamma1)  ## for the first covariate
    bet2_time1 <- sum(Bs2_Spe1*temp_gamma2)  ## for the second covariate
    temp_mu_time1 <- bet0_time1+xx[1]*bet1_time1+xx[2]*bet2_time1
    cdf_time1 <- pnorm(yy1-temp_mu_time1, mean=0, sd=sqrt(varEst))
    #############################################################################
    ##### estimate the marginal density and cdf for tt2
    loc_Spe2 <- which(AgeGrid==tt2)
    Bs0_Spe2 <- Bs0[loc_Spe2,]
    Bs1_Spe2 <- Bs1[loc_Spe2,]
    Bs2_Spe2 <- Bs2[loc_Spe2,]
    bet0_time2 <- sum(Bs0_Spe2*temp_gamma0)  ## the intercept
    bet1_time2 <- sum(Bs1_Spe2*temp_gamma1)  ## for the first covariate
    bet2_time2 <- sum(Bs2_Spe2*temp_gamma2)  ## for the second covariate
    temp_mu_time2 <- bet0_time2+xx[1]*bet1_time2+xx[2]*bet2_time2
    cdf_time2 <- pnorm(yy2-temp_mu_time2, mean=0, sd=sqrt(varEst))

    ######## estimate the bivariate cdf function   ############################
    lamPara <- c(copulaPara)
    temp_delta <- sum(kronecker(Bs0_Spe1, Bs0_Spe2)*lamPara)
    tempCop <- normalCopula(param=temp_delta, dim=2, dispstr="un") 
    CopCDF <- pCopula(u=c(cdf_time1, cdf_time2), copula=tempCop)
    return(CopCDF)
}

################################## For two outcome at the same time points
#####################################################################################################################
##########calculate the log-likelihood function of copula (Gaussian) for two outcomes at the same time ########
BiCop_EqualTime_LogLik <- function(lambda, basisFun, cdfMargin1, cdfMargin2, ID)
{
    uniID <- unique(ID)
    Nsub <- length(uniID)
    Ni <- table(ID)
    logCopFun <- function(paraCop)
    {
        tempCop <- normalCopula(param=paraCop[3], dim=2, dispstr="un") 
        return(log(dCopula(u=c(paraCop[1], paraCop[2]), copula=tempCop)))
    }
    delta_All <- basisFun%*%matrix(lambda, ncol=1)
    LogLik_sig <- apply(cbind(cdfMargin1, cdfMargin2, delta_All), MARGIN=1, logCopFun)
    w.vec <- 1/rep(Ni, Ni)
    Res_logLik <- sum(w.vec*LogLik_sig)
    return(Res_logLik)
}

#########################################################################################
####### estimate the copula parameter by optimation
BiCop_EqualTime_paraEst <- function(Y1, Y2, ID, Age, X, numKnot, deg, Normal1, Normal2)
{
  ### for the first outcome
  tempEst1 <- MarginCDF(Y=Y1, ID=ID, Age=Age, X=X, numKnot=numKnot, deg=deg, Normal=Normal1) 
  cdfEst1 <- tempEst1$cdf_est
  ### for the second outcome
  tempEst2 <- MarginCDF(Y=Y2, ID=ID, Age=Age, X=X, numKnot=numKnot, deg=deg, Normal=Normal2) 
  cdfEst2 <- tempEst2$cdf_est
  ##################
  Age_lower <- min(Age)
  Age_upper <- max(Age)
  KN_bi <- seq(from=Age_lower, to=Age_upper, length=numKnot[1])[-c(1,numKnot)]
  bs_Fun <- bs(x=Age, knots=KN_bi, degree=deg, intercept=T)
  para.old <- rep(0.5, ncol(bs_Fun))
  ResEst <- optim(par=para.old, fn=BiCop_EqualTime_LogLik, gr=NULL, control=list(fnscale=-1, reltol=0.001), basisFun=bs_Fun,
                  cdfMargin1=cdfEst1, cdfMargin2=cdfEst2, ID=ID)
  lamEst <- ResEst$par
  return(lamEst)
}

#########################################################################################
###########################
BiCop_EqualTime_disEst <- function(yy1, yy2, gammList1, gammList2, lamEst, varEst1, varEst2,
                                   tt, xx, AgeGrid, numKnot, deg, Normal1, Normal2)
{
    #############################################################################
    #####
    Age_lower <- min(AgeGrid)
    Age_upper <- max(AgeGrid)
    KN0 <- seq(from=Age_lower, to=Age_upper, length=numKnot[1])[-c(1,numKnot[1])]
    KN1 <- seq(from=Age_lower, to=Age_upper, length=numKnot[2])[-c(1,numKnot[2])]
    KN2 <- seq(from=Age_lower, to=Age_upper, length=numKnot[3])[-c(1,numKnot[3])]
    Bs0 <- bs(AgeGrid, knots=KN0, degree=deg, intercept=T)
    Bs1 <- bs(AgeGrid, knots=KN1, degree=deg, intercept=T)
    Bs2 <- bs(AgeGrid, knots=KN2, degree=deg, intercept=T)
    loc_Spe <- which(AgeGrid==tt)
    Bs0_Spe <- Bs0[loc_Spe,]
    Bs1_Spe <- Bs1[loc_Spe,]
    Bs2_Spe <- Bs2[loc_Spe,]
    
    #############################################################################
    ### for the first outcome
    temp1_gamma0 <- gammList1[[1]]
    temp1_gamma1 <- gammList1[[2]]
    temp1_gamma2 <- gammList1[[3]]
    ##### estimate the marginal density and cdf for tt
    bet0_outcome1 <- sum(Bs0_Spe*temp1_gamma0)  ## the intercept
    bet1_outcome1 <- sum(Bs1_Spe*temp1_gamma1)  ## for the first covariate
    bet2_outcome1 <- sum(Bs2_Spe*temp1_gamma2)  ## for the second covariate
    temp_mu_outcome1 <- bet0_outcome1+bet1_outcome1*xx[1]+bet2_outcome1*xx[2]
    cdf_outcome1 <- pnorm(yy1-temp_mu_outcome1, mean=0, sd=sqrt(varEst1))

    ### for the second outcome
    temp2_gamma0 <- gammList2[[1]]
    temp2_gamma1 <- gammList2[[2]]
    temp2_gamma2 <- gammList2[[3]]
    ##### estimate the marginal density and cdf for tt
    bet0_outcome2 <- sum(Bs0_Spe*temp2_gamma0)  ## the intercept
    bet1_outcome2 <- sum(Bs1_Spe*temp2_gamma1)  ## for the first covariate
    bet2_outcome2 <- sum(Bs2_Spe*temp2_gamma2)  ## for the second covariate
    temp_mu_outcome2 <- bet0_outcome2+xx[1]*bet1_outcome2+xx[2]*bet2_outcome2
    cdf_outcome2 <- pnorm(yy2-temp_mu_outcome2, mean=0, sd=sqrt(varEst2))

    ### the basis function at the interested time point
    tempdelta <- sum(Bs0_Spe*lamEst)
    tempCop <- normalCopula(param=tempdelta, dim=2, dispstr="un") 
    CopCDF <- pCopula(u=c(cdf_outcome1, cdf_outcome2), copula=tempCop)
    return(CopCDF)
}


