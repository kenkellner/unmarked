test.fitModels.occu <- function() {

  set.seed(123)
  data(frogs)
  pferUMF <- unmarkedFrameOccu(pfer.bin)
  siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
  #Add some correlated covariates
  obsvar2 = rnorm(numSites(pferUMF) * obsNum(pferUMF))
  obsvar3 = rnorm(numSites(pferUMF) * obsNum(pferUMF),mean=obsvar2,sd=0.5)
  obsCovs(pferUMF) <- data.frame(
                        obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)),
                        obsvar2=obsvar2,obsvar3=obsvar3)
     
  fit <- occu(~ obsvar1+obsvar2 ~sitevar1, pferUMF)
  
  #Check model name
  checkEquals(unmarked:::modName(fit), "psi(sitevar1)p(obsvar1+obsvar2)")
  
  #Check fitting 2 models
  fm <- fitModels(fit, formula=list(~1~1, ~1~sitevar1), quiet=T)
  checkEquals(class(fm)[1], "unmarkedFitList")
  checkEquals(length(fm@fits), 2)
  checkEquals(lapply(fm@fits, function(x) x@formula),
              list("psi(.)p(.)"=~1~1, "psi(sitevar1)p(.)"=~1~sitevar1))

  #Check dredge
  dr <- dredge(fit, quiet=T)
  checkEquals(class(fm)[1], "unmarkedFitList")
  checkEquals(length(dr@fits), 8)
  m <- modSel(dr)
  top <- m@Full$model[1]
  checkEquals(top, "psi(sitevar1)p(.)")
}

test.fitModels.pcountOpen <- function() {

set.seed(3)
M <- 50
T <- 5
lambda <- 4
gamma <- 1.5
omega <- 0.8
p <- 0.7
y <- N <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
    S[,t] <- rbinom(M, N[,t], omega)
    G[,t] <- rpois(M, gamma)
    N[,t+1] <- S[,t] + G[,t]
}
y[] <- rbinom(M*T, N, p)
     
# Prepare data
sc <- as.data.frame(matrix(rnorm(50*2),ncol=2))
names(sc) <- c('sc1','sc2')
umf <- unmarkedFramePCO(y = y, siteCovs=sc, numPrimary=T)

# Fit model and backtransform
m1 <- pcountOpen(~sc1+sc2, ~1, ~1, ~1, umf, K=20)

dr <- dredge(m1, quiet=T)
checkEquals(class(dr)[1], "unmarkedFitList")
checkEquals(length(dr@fits), 4)

m <- modSel(dr)
top <- m@Full$model[1]
checkEquals(top, "lam(sc2)gamConst(.)omega(.)p(.)")
}
