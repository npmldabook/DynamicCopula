Get_epsilon <- function(Y,  ID, Age, X, numKnot, deg)
{
  ## estimate the residual
  Age_lower <- min(Age)
  Age_upper <- max(Age)
  #Get knots and b-spline basis for 3 covariates
  KN0 <- seq(from=Age_lower, to=Age_upper, length=numKnot[1])[-c(1,numKnot[1])]
  KN1 <- seq(from=Age_lower, to=Age_upper, length=numKnot[2])[-c(1,numKnot[2])]
  KN2 <- seq(from=Age_lower, to=Age_upper, length=numKnot[3])[-c(1,numKnot[3])]
  KN3 <- seq(from=Age_lower, to=Age_upper, length=numKnot[4])[-c(1,numKnot[4])]
  T0_bs <- bs(x=Age, knots=KN0, degree=deg, intercept=T)
  T1_bs <- bs(x=Age, knots=KN1, degree=deg, intercept=T)
  T2_bs <- bs(x=Age, knots=KN2, degree=deg, intercept=T)
  T3_bs <- bs(x=Age, knots=KN3, degree=deg, intercept=T)
  X0_bs <- T0_bs  ## the intercept term
  X1_bs <- matrix(X[,1], nrow=length(X[,1]), ncol=ncol(T1_bs))*T1_bs
  X2_bs <- matrix(X[,2], nrow=length(X[,2]), ncol=ncol(T2_bs))*T2_bs
  X3_bs <- matrix(X[,3], nrow=length(X[,3]), ncol=ncol(T3_bs))*T3_bs
  ######### Fit the WLS estimate to the the coefficient for the BS, 
  ######obtain the residual of WLS
  Ni <- table(ID)
  Nsub <- length(unique(ID))
  wi <- 1/(Nsub*rep(Ni, Ni))
  ###model
  fit_WLS <- lm(Y~0+X0_bs+X1_bs+X2_bs+X3_bs, weights=wi)
  resError <- fit_WLS$residuals
  return(resError)
}


ObjectiveFun_LSE=function(lambda, basisFun,resError){
  tempPara <- lambda
  #Get the residuals
  uniID <- unique(ID)
  Nsub <- length(uniID)
  Ni <- table(ID)
  dd <- cumsum(c(0,Ni))
  SamLoc <- cbind(dd[1:(length(dd)-1)], dd[2:length(dd)])
  Samlist <- list()
  for(kk in 1:Nsub)  Samlist[[kk]] <- (SamLoc[kk,1]+1):SamLoc[kk,2]
  sum_of_squre<- rep(0,Nsub)
  NumCom<- rep(0,Nsub)
  #For every subject
  for(i in 1:Nsub){
    ind=Samlist[[i]]
    residuals_i<- resError[ind]
    #Get delta
    tempMat <- matrix(0, nrow=ncol(basisFun), ncol=ncol(basisFun))
    tempMat[upper.tri(tempMat, diag=T)] <- tempPara
    tempMat <- tempMat + t(tempMat) - diag(diag(tempMat))
    lamPara <- c(tempMat)
    tempBasis <- basisFun[ind, ]
    NumRepeat <- length(ind)
    AllCom0 <- expand.grid(1:NumRepeat,1:NumRepeat)
    selectLoc <- which(AllCom0[,1]<AllCom0[,2])
    AllCom <- AllCom0[selectLoc,]
    NumCom[i] <- nrow(AllCom)
    #Calculate delta and product of residuals
    temp.delta <- apply(AllCom, MARGIN=1, function(u) sum(kronecker(tempBasis[u[1],], tempBasis[u[2],])*lamPara))
    prod_res<- apply(AllCom,MARGIN=1, function(u) residuals_i[u[1]]*residuals_i[u[2]])
    sum_of_squre[i]<- sum((prod_res-temp.delta)^2)
  }
  LSCV=sum(sum_of_squre)/sum(NumCom)
  return(LSCV)
}


