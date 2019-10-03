\name{dredge}
\alias{dredge}
\alias{dredge-methods}
\alias{dredge,unmarkedFit-method}
\alias{dredge,unmarkedFitGDS-method}
\alias{dredge,unmarkedFitGMM-method}
\alias{dredge,unmarkedFitGPC-method}
\alias{dredge,unmarkedFitOccuFP-method}
\alias{dredge,unmarkedFitColExt-method}
\alias{dredge,unmarkedFitPCO-method}
\alias{dredge,unmarkedFitOccuTTD-method}
\alias{dredge,unmarkedFitOccuMulti-method}
\alias{dredge,unmarkedFitOccuMS-method}

\title{Generate models for all possible subsets of covariates}

\description{Using a fitted global model \code{globalModel}, finds
  all possible combinations of covariates across all submodels (state,
  detection, etc.) and fits a model for each combination. Model fits can 
  be done in parallel.
}

\usage{
\S4method{dredge}{unmarkedFit}(globalModel, quiet=FALSE, parallel=FALSE, ...)
}

\arguments{
  \item{globalModel}{A fitted global model inheriting class \code{unmarkedFit}, 
    with all covariates of interest included}
  \item{quiet}{Suppress progress bar}
  \item{parallel}{If \code{TRUE}, fit models in parallel}
  \item{...}{Other arguments passed to the method}
} 

\value{\code{fitList} object containing model fits for each combination of 
  covariates. Use \code{modSel} on the output to rank models by AIC.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{fitModels}},
  \code{\link{fitList}},
  \code{\link{unmarkedFit}}
}

\examples{

\dontrun{
#Get data
data(ovendata)
#Add some fake obs covs
fakeOC <- data.frame(matrix(rnorm(70*4*2),ncol=2))
names(fakeOC) <- c('fakeA','fakeB')

#Build frame
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
    obsCovs=fakeOC, type = "removal")

#Fit global model
fm1 <- multinomPois(~ fakeA + fakeB ~ ufc + trba, ovenFrame)

#Dredge global model
dr <- dredge(fm1)

#Get AIC table
modSel(dr)
}
}
