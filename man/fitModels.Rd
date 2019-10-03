\name{fitModels}
\alias{fitModels}
\title{Fit Multiple \code{unmarked} Models}

\description{
  Fit multiple variations on an \code{unmarked} model with different
  arguments, e.g. different formulas, with a single function call. Start
  by fitting a single "reference" model (\code{refModel}) with an unmarked
  fitting function. Then in the call to \code{fitModels}, provide this fitted 
  model along with additional named arguments that match the names of arguments 
  for the fitting function. Each of these provided arguments should be lists
  of equal length. \code{fitModels} will iterate through each combination
  of provided arguments and fit a new model. Any arguments to the fitting 
  function that are not provided to \code{fitModels} default to the values 
  in the reference model. You cannot provide different datasets; all models 
  use the original dataset in the reference model. Models may optionally be fit 
  in parallel.
}

\usage{
fitModels(refModel, quiet=FALSE, parallel=FALSE, ...)
}

\arguments{
  \item{refModel}{A fitted reference model inheriting class \code{unmarkedFit}}
  \item{quiet}{Suppress progress bar}
  \item{parallel}{If \code{TRUE}, fit models in parallel}
  \item{...}{Arguments passed to the fitting function. Must be lists of 
  equal length}
}

\value{\code{fitList} object containing model fits for each combination of 
  covariates. Use \code{modSel} on the output to rank models by AIC.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\examples{

\dontrun{
#Get data
data(ovendata)
#Add some fake obs covs

#Build frame
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")

#Fit reference model
fm1 <- multinomPois(~ 1 ~ 1, ovenFrame)

#Fit two new models with different formulas
fits <- fitModels(fm1, formula=list(~1~ufc, ~1~ufc+trba))
}
}
