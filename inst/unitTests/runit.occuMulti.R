test.occuMulti.fit.simple.1 <- function() {

  y <- list(matrix(rep(1,10),5,2),
            matrix(rep(1,10),5,2))
  umf <- unmarkedFrameOccuMulti(y = y)
  fm <- occuMulti(detformulas=rep("~1",2),
                  stateformulas=rep("~1",3), data = umf, se=FALSE)

  #Probably should not be calling predict here b/c unit test
  #but complicated to get actual occupancy prob otherwise
  occ <- predict(fm,'state')$Predicted[1,1]
  checkEqualsNumeric(occ,1, tolerance = 1e-4)

  detlist <- predict(fm,'det')
  det <- sapply(detlist,function(x) x[1,1])
  checkEqualsNumeric(det, rep(1,length(detlist)), tolerance= 1e-4)

  #Check fitList
  fl <- fitList(fm, fm)
  checkEquals(class(fl)[1],"unmarkedFitList")
  checkEqualsNumeric(length(fl@fits), 2)

  #Check error when random effect in formula
  checkException(occuMulti(detformulas=rep("~1",2),
                           stateformulas=c("~(1|group)",rep("~1",2)), umf))
}

test.occuMulti.fit.simple.0 <- function() {

  y <- list(matrix(rep(0,10),5,2),
            matrix(rep(0,10),5,2))
  umf <- unmarkedFrameOccuMulti(y = y)
  fm <- occuMulti(detformulas=rep("~1",2),
                  stateformulas=rep("~1",3), data = umf, se=FALSE)

  occ <- predict(fm,'state')$Predicted[1,1]
  checkEqualsNumeric(occ,0, tolerance = 1e-4)

  detlist <- predict(fm,'det')
  det <- sapply(detlist,function(x) x[1,1])
  checkEqualsNumeric(det, rep(0,length(detlist)), tolerance= 1e-4)


}

test.occuMulti.fit.covs <- function() {

  y <- list(matrix(rep(0:1,10),5,2),
            matrix(rep(0:1,10),5,2))

  set.seed(123)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)
  stateformulas <- c('~occ_cov1','~occ_cov2','~occ_cov3')
  detformulas <- c('~det_cov1','~det_cov2')

  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)

  occ <- fm['state']
  det <- fm['det']

  checkEqualsNumeric(coef(occ), c(5.36630,0.79876,5.45492,-0.868451,9.21242,1.14561),
                     tolerance = 1e-4)
  checkEqualsNumeric(coef(det), c(-0.27586,-0.81837,-0.09537,0.42334), tolerance = 1e-4)

  fit <- fitted(fm)
  checkEqualsNumeric(length(fit),2)
  checkEqualsNumeric(sapply(fit,function(x) x[1,1]),c(0.14954,0.30801), tol = 1e-4)

  res <- residuals(fm)
  checkEqualsNumeric(length(res),2)
  checkEqualsNumeric(sapply(res,function(x) x[1,1]),c(-0.14954,-0.30801), tol= 1e-4)

  gp <- getP(fm)
  checkEqualsNumeric(length(gp), 2)
  checkEqualsNumeric(dim(gp[[1]]), c(N,J))

  #Check site cov can be used in detection formula
  detformulas <- c('~occ_cov1','~det_cov2')
  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)
  checkEqualsNumeric(coef(fm,'det')[2],3.355328e-05, tol=1e-4)
}

test.occuMulti.fit.NA <- function() {

  y <- list(matrix(rep(0:1,10),5,2),
            matrix(rep(0:1,10),5,2))

  set.seed(456)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~occ_cov1','~occ_cov2','~occ_cov3')
  detformulas <- c('~det_cov1','~det_cov2')

  #Check error thrown when missing site covariates
  occ_covsNA <- occ_covs
  occ_covsNA[1,1] <- NA
  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covsNA, obsCovs = det_covs)
  checkException(occuMulti(detformulas, stateformulas, data=umf, se=FALSE))

  #Check for warning when missing detection
  yna <- y
  yna[[1]][1,1] <- NA
  umf <- unmarkedFrameOccuMulti(y = yna, siteCovs = occ_covs, obsCovs = det_covs)

  options(warn=2)
  checkException(occuMulti(detformulas, stateformulas, data=umf, se=FALSE))

  options(warn=1)

  #Check correct answer given when missing detection
  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)
  checkEqualsNumeric(coef(fm)[c(1,7)], c(6.63207,0.35323), tol= 1e-4)

  fit <- fitted(fm)
  checkTrue(is.na(fit[[1]][1,1]))

  res <- residuals(fm)
  checkTrue(is.na(res[[1]][1,1]))

  gp <- getP(fm)
  checkTrue(is.na(gp[[1]][1,1]))

  #Check error thrown when all detections are missing
  yna[[1]][1,] <- NA
  umf <- unmarkedFrameOccuMulti(y = yna, siteCovs = occ_covs, obsCovs = det_covs)
  checkException(occuMulti(detformulas, stateformulas, data=umf, se=FALSE))

  #Check warning when missing covariate value on detection
  det_covsNA <- det_covs
  det_covsNA[1,1] <- NA
  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covsNA)
  options(warn=2)
  checkException(occuMulti(detformulas,stateformulas,data=umf, se=FALSE))
  options(warn=1)
}

test.occuMulti.fit.fixed0 <- function(){

  y <- list(matrix(rep(0:1,10),5,2),
            matrix(rep(0:1,10),5,2))

  set.seed(123)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~occ_cov1','~occ_cov2','0')
  detformulas <- c('~det_cov1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)

  occ <- fm['state']
  checkEqualsNumeric(length(coef(occ)),4)
  checkEqualsNumeric(coef(occ),c(12.26043,0.61183,12.41110,0.18764),tol=1e-4)


  stateformulas <- c('~occ_cov1','~occ_cov2')
  fm2 <- occuMulti(detformulas, stateformulas, data = umf, maxOrder=1,se=FALSE)

  occ <- fm2['state']
  checkEqualsNumeric(length(coef(occ)),4)
  checkEqualsNumeric(coef(occ),c(12.26043,0.61183,12.41110,0.18764),tol=1e-4)

}

test.occuMulti.predict <- function(){

  set.seed(123)
  y <- list(matrix(rbinom(40,1,0.2),20,2),
            matrix(rbinom(40,1,0.3),20,2))

  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~occ_cov1','~occ_cov2','0')
  detformulas <- c('~det_cov1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf)

  prState <- predict(fm, type='state')
  checkEqualsNumeric(sapply(prState,function(x) x[1,1]),
                     c(0.30807707,0.20007250,0.04234835,0.73106618),tol=1e-4)
  prDet <- predict(fm, type='det')
  checkEqualsNumeric(as.numeric(prDet$sp2[1,]),
                     c(0.190485,0.12201,0.0475270,0.525988), tol=1e-4)

  #Check with newdata
  nd <- siteCovs(umf)[1:2,]
  pr_nd <- predict(fm, type='state', newdata=nd)$Predicted
  checkEqualsNumeric(pr_nd[,1],c(0.3080771,0.3196486), tol=1e-4)
  nd <- siteCovs(umf)[1:2,]
  pr_nd <- predict(fm, type='state', newdata=nd, species=1, cond=2)$Predicted
  checkEqualsNumeric(pr_nd,c(0.3858233,0.5402935), tol=1e-4)
  #Make sure it works with newdata having only one row
  nd <- siteCovs(umf)[1,]
  pr_nd <- predict(fm, type='state', newdata=nd)$Predicted
  checkEqualsNumeric(pr_nd[,1],c(0.3080771), tol=1e-4)
  pr_nd <- predict(fm, type='state', newdata=nd, species=1, cond=2)$Predicted
  checkEqualsNumeric(pr_nd,c(0.3858233), tol=1e-4)

  stateformulas <- c('~1','~1','0')
  detformulas <- c('~1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf)

  prState <- predict(fm, type='state')
  checkEqualsNumeric(sapply(prState,function(x) x[1,1]),
                     c(0.475928,0.2548407,0.01496681,0.86713789),tol=1e-4)
  prDet <- predict(fm, type='det')
  checkEqualsNumeric(as.numeric(prDet$sp2[1,]),
                     c(0.20494,0.11865,0.0582563,0.517888), tol=1e-4)

  #Check predicting co-occurrence
  nd <- siteCovs(umf)[1:2,]
  pr_all <- predict(fm, type='state', se=F)$Predicted[1:2,1]
  pr_nd <- predict(fm, type='state', newdata=nd, species=c(1,2))$Predicted
  checkEqualsNumeric(pr_nd,pr_all, tol=1e-4)

  #Check with site cov in detection formula
  stateformulas <- c('~occ_cov2','~1','0')
  detformulas <- c('~occ_cov1','~det_cov2')
  fm <- occuMulti(detformulas, stateformulas, data = umf)
  pr_state_actual <- predict(fm, "state")
  checkEqualsNumeric(length(pr_state_actual), 4)
  checkEqualsNumeric(pr_state_actual$Predicted[1,1], 0.729927907, tol=1e-5)
  checkEqualsNumeric(nrow(pr_state_actual$Predicted), 20)

  pr_det_actual <- predict(fm, "det")
  checkEqualsNumeric(length(pr_det_actual), 2)
  checkEqualsNumeric(pr_det_actual$sp1$Predicted[1], 0.1448311, tol=1e-5)
  checkEqualsNumeric(nrow(pr_det_actual$sp1), 20*2)

  #with newdata
  pr_state_nd <- predict(fm, "state", newdata=data.frame(occ_cov2=0))
  checkEqualsNumeric(length(pr_state_nd), 4)
  checkEqualsNumeric(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  checkEqualsNumeric(nrow(pr_state_nd$Predicted), 1)

  pr_det_nd <- predict(fm, "det", newdata=data.frame(occ_cov1=0, det_cov2=0))
  checkEqualsNumeric(length(pr_det_nd), 2)
  checkEqualsNumeric(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  checkEqualsNumeric(nrow(pr_state_nd$Predicted), 1)

  #With maxOrder set
  stateformulas <- c('~occ_cov2','~1')
  detformulas <- c('~occ_cov1','~det_cov2')

  fm <- occuMulti(detformulas, stateformulas, data = umf, maxOrder=1)

  pr_state <- predict(fm, "state")
  checkEqualsNumeric(length(pr_state), 4)
  checkEqualsNumeric(pr_state$Predicted[1,1], 0.729927907, tol=1e-5)
  checkEqualsNumeric(nrow(pr_state$Predicted), 20)

  pr_state_nd <- predict(fm, "state", newdata=data.frame(occ_cov2=0))
  checkEqualsNumeric(length(pr_state_nd), 4)
  checkEqualsNumeric(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  checkEqualsNumeric(nrow(pr_state_nd$Predicted), 1)

  pr_det <- predict(fm, "det")
  checkEqualsNumeric(length(pr_det), 2)
  checkEqualsNumeric(pr_det$sp1$Predicted[1], 0.1448311, tol=1e-5)
  checkEqualsNumeric(nrow(pr_det$sp1), 20*2)

  pr_det_nd <- predict(fm, "det", newdata=data.frame(occ_cov1=0, det_cov2=0))
  checkEqualsNumeric(length(pr_det_nd), 2)
  checkEqualsNumeric(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  checkEqualsNumeric(nrow(pr_state_nd$Predicted), 1)

  #getP with maxOrder set
  gp <- getP(fm)
  checkEquals(length(gp), 2)
  checkEquals(dim(gp[[1]]), c(20,2))

  #simulate with maxOrder set
  s <- simulate(fm, 2)
  checkTrue(inherits(s, "list"))
  checkEquals(length(s), 2)
  checkEquals(dim(s[[1]][[1]]), c(N, J))

  #fitList with maxOrder set
  fm2 <- occuMulti(c("~1","~1"), c("~1","~1"), umf, maxOrder=1)
  fl2 <- fitList(fm, fm2)
  checkTrue(inherits(fl2, "unmarkedFitList"))
  ms <- modSel(fl2)
  checkTrue(inherits(ms, "unmarkedModSel"))

  #fitted with maxOrder set
  ft <- fitted(fm)
  checkEquals(length(ft), 2)

  #parboot with maxOrder set
  pb <- parboot(fm, nsim=2)
  checkTrue(inherits(pb, "parboot"))
}

test.occuMulti.predict.NA <- function(){

  set.seed(123)
  y <- list(matrix(rbinom(40,1,0.2),20,2),
            matrix(rbinom(40,1,0.3),20,2))

  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')
  det_covs[1,1] <- NA

  stateformulas <- c('~occ_cov1','~occ_cov2','0')
  detformulas <- c('~det_cov1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf)

  prDet <- predict(fm, type='det')
  checkTrue(all(is.na(prDet$sp1[1,])))
  checkEqualsNumeric(as.numeric(prDet$sp1[2,]),
                     c(0.49781,0.19621,0.175514,0.8219401), tol=1e-4)

  #Check that you can predict with NAs in siteCovs
  newdata <- siteCovs(umf)
  newdata[1,1] <- NA
  prOcc <- predict(fm, type='state', newdata=newdata)
  checkTrue(all(is.na(prOcc$Predicted[1,])))
  checkTrue(all(!is.na(sapply(prOcc,`[`,2,1))))
  prOcc_sp <- predict(fm, type='state', species=1, newdata=newdata)
  checkTrue(all(is.na(prOcc_sp[1,])))
  checkTrue(all(!is.na(prOcc_sp[2,])))
  checkEqualsNumeric(prOcc_sp$Predicted[2],0.4731427, tol=1e-4)
  prOcc_cond <- predict(fm, type='state', species=1, cond=2, newdata=newdata)
  checkTrue(all(is.na(prOcc_cond[1,])))
  checkTrue(all(!is.na(prOcc_cond[2,])))
  checkEqualsNumeric(prOcc_sp$Predicted[2],0.4731427, tol=1e-4)
}


test.occuMulti.predict.complexFormulas <- function(){

  #Check scale(), etc
  set.seed(123)
  y <- list(matrix(rbinom(40,1,0.2),20,2),
            matrix(rbinom(40,1,0.3),20,2))

  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3, mean=2),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2, mean=3),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~scale(occ_cov1)','~1','0')
  detformulas <- c('~scale(det_cov1)','~1')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf)

  #Check with newdata; contents of newdata should not
  #effect resulting predictions (scale should be based on
  #original data)
  nd <- siteCovs(umf)[1:5,]
  pr_nd <- predict(fm, type='state', newdata=nd, se=F)$Predicted
  nd <- siteCovs(umf)[1:2,]
  pr_nd2 <- predict(fm, type='state', newdata=nd, se=F)$Predicted
  nd <- siteCovs(umf)[c(1,1),]
  pr_nd3 <- predict(fm, type='state', newdata=nd, se=F)$Predicted

  checkEqualsNumeric(pr_nd[1:2,], pr_nd2)
  checkEqualsNumeric(pr_nd[c(1,1),], pr_nd3)

  #Check for factor level handling
  occ_covs$occ_fac <- factor(sample(c('a','b','c'),N,replace=T))

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)
  stateformulas <- c('~occ_fac','~1','~1')
  fm <- occuMulti(detformulas, stateformulas, data = umf)

  nd <- siteCovs(umf)[1:2,]
  pr_nd <- predict(fm, type='state', newdata=nd, se=F)$Predicted

  nd2 <- data.frame(occ_fac=factor(c('a','b'),levels=c('a','b','c')))
  pr_nd2 <- predict(fm, type='state', newdata=nd2, se=F)$Predicted

  checkEqualsNumeric(pr_nd, pr_nd2[c(2,1),])

  nd3 <- data.frame(occ_fac=c('a','b'))
  pr_nd3 <- predict(fm, type='state', newdata=nd3, se=F)$Predicted

  checkEqualsNumeric(pr_nd, pr_nd3[c(2,1),])

  nd4 <- data.frame(occ_fac=factor(c('a','d'),levels=c('a','d')))
  checkException(predict(fm, type='state', newdata=nd4, se=F))

  #Check that predicting detection also works
  nd5 <- data.frame(det_cov1 = rnorm(5))
  pr_nd5 <- predict(fm, type='det', newdata=nd5)
  checkEqualsNumeric(sapply(pr_nd5, nrow), c(5,5))
  checkEqualsNumeric(pr_nd5$sp1$Predicted[1], 0.1680881)
}

test.occuMulti.fix.NA.mismatch <- function(){
  y1 <- matrix(rbinom(10,1,0.5), nrow=5)
  y1[1,1] <- NA

  y2 <- matrix(rbinom(10,1,0.5), nrow=5)
  y2[1,1] <- NA

  y3 <- matrix(rbinom(10,1,0.5), nrow=5)
  y3[1,1] <- NA
  y3[5,1] <- NA

  ylist <- list(y1=y1,y2=y2,y3=y3)

  umf <- unmarkedFrameOccuMulti(y=ylist)

  options(warn=2)
  checkException(unmarkedFrameOccuMulti(y=ylist))
  options(warn=0)

  pre_na <- sapply(ylist, function(x) sum(is.na(x)))
  post_na <- sapply(umf@ylist, function(x) sum(is.na(x)))

  checkTrue(any(pre_na[1] != pre_na))
  checkTrue(!any(post_na[1] != post_na))
}
