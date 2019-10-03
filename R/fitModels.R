
#Method to generate short model names in form psi(.)p(.)-----------------------

format_forms <- function(sh, forms){
  forms <- sapply(forms, function(x) paste(x, collapse=""))
  forms <- gsub("~","", forms)
  forms <- gsub(" ","", forms)
  forms <- ifelse(forms == "1", ".", forms)
  out <- sapply(1:length(forms), 
                function(i) paste0(sh[i],'(',forms[i],')'))
  paste(out, collapse="")
}

setGeneric("modName", function(mod, ...) standardGeneric("modName"))

setMethod("modName", "unmarkedFit", function(mod){
  sh <- sapply(mod@estimates@estimates, function(x) x@short.name)
  forms <- split_form(mod@formula)
  format_forms(sh, forms[c(2,1)])
}) 
          
.modName_gmods <- function(mod, ...){
  sh <- sapply(mod@estimates@estimates, function(x) x@short.name)
  format_forms(sh, mod@formlist)
}
          
setMethod("modName", "unmarkedFitGDS", .modName_gmods)
setMethod("modName", "unmarkedFitGMM", .modName_gmods)
setMethod("modName", "unmarkedFitGPC", .modName_gmods)

setMethod("modName", "unmarkedFitOccuFP", function(mod, ...){
  sh <- c("psi", "p", "FP", "B")
  formlist <- list(mod@stateformula, mod@detformula, mod@FPformula,
                mod@Bformula)
  format_forms(sh, formlist)
})

setMethod("modName", "unmarkedFitColExt", function(mod, ...){
  sh <- sapply(mod@estimates@estimates, function(x) x@short.name)
  formlist <- list(mod@psiformula, mod@gamformula, mod@epsformula,
                   mod@detformula)
  format_forms(sh, formlist)
})

setMethod("modName", "unmarkedFitPCO", function(mod, ...){
  sh <- sapply(mod@estimates@estimates, function(x) x@short.name)
  sh <- sh[!sh%in%c("psi", "alpha")]
  formlist <- mod@formlist
  if(!"iota"%in%sh) formlist <- formlist[1:(length(formlist)-1)]
  format_forms(sh, formlist)
})

setMethod("modName", "unmarkedFitOccuTTD", function(mod, ...){ 
  sh <- sapply(mod@estimates@estimates, function(x) x@short.name)
  formlist <- list(mod@psiformula, mod@gamformula, mod@epsformula,
                   mod@detformula)
  if(!"col"%in%sh){
    formlist <- formlist[c(1,4)]
  }
  format_forms(sh, formlist)
})

#Generate formulas for all variable subsets in a global model------------------

#Extract all terms from a formula
get_terms <- function(formula){ 
  cf <- as.character(formula)
  cf <- cf[length(cf)]
  if(grepl("*", cf, fixed=TRUE)){
    stop("Interactions not supported", call.=FALSE)
  }
  strsplit(gsub(" ", "", cf), '+',fixed=TRUE)[[1]]
}

#Get all combinations of terms
get_term_combs <- function(formula){
  terms <- get_terms(formula)
  maxn <- length(terms)
  out <- lapply(1:maxn, function(i) utils::combn(terms, m=i, simplify=FALSE))
  lapply(rapply(out, enquote, how="unlist"), eval)
}

#Generate all possible subset formulas from a global formula
#Output is a list of strings of formulas
all_forms <- function(formula){
  gt <- get_term_combs(formula)
  dotmod <- NULL
  if(gt[[1]]!="1") dotmod <- list("1")
  forms_raw <- c(dotmod, get_term_combs(formula))
  lapply(forms_raw, function(x) paste("~",paste(x, collapse="+")))
}

#Count number of formulas in a "squashed" formula object like with occu()
nforms <- function(formula){
  sum(grepl("~", as.character(formula)))
}

#Split apart a "squashed" formula
split_form <- function(formula){
  if(nforms(formula)!=2) stop("Requires 2 formulas")
  form1 <- as.formula(formula[[2]])
  form2 <- as.formula(paste("~", formula[3], sep=""))
  list(form1,form2)
}

#Take one or more formulas
#Output list(s) of all combinations of variables in those formulas
#If merge=TRUE, squash pairs of formulas together to be used eg by occu()
formulaCombs <- function(..., merge=FALSE){

  fl <- list(...)
  if(length(fl)==1){
    if(nforms(fl[[1]])==2){
      fl <- split_form(fl[[1]])
      merge <- TRUE
    } else{
      return(lapply(all_forms(fl[[1]]), as.formula))
    }
  }
 
  eg <- expand.grid(lapply(fl, all_forms))
  if(merge){
    out <- lapply(1:nrow(eg), 
                  function(i) paste(unlist(eg[i,]), collapse=" "))
    out <- lapply(out, as.formula)
  } else{
    out <- lapply(eg, as.list)
    out <- lapply(out, function(x) lapply(x, as.formula))
  }
  out
}

#Run a bunch of models based on lists of parameter inputs----------------------

#Convert a series of arguments in ... to a list
args_to_list <- function(...){
  lens <- sapply(list(...), length)
  if(any(lens!=lens[1])){
    stop("All argument lists must be equal length")
  }

  temp <- cbind(...)
  lapply(1:nrow(temp), FUN=function(x) temp[x,])
}

#Fit series of models
#Arguments passed to ... should be lists of equal length
#corresponding to arguments of original "reference" model
fitModels <- function(refModel, quiet=FALSE, parallel=FALSE, ...){
  
  if(!inherits(refModel, "unmarkedFit")){
    stop("refModel must be a fitted unmarked model", call.=FALSE)
  }

  if(quiet){
    old_opt <- pbapply::pboptions()
    pbapply::pboptions(type="none")
    on.exit(pbapply::pboptions(old_opt))
  }

  plist <- args_to_list(...)  

  new_fit <- function(plist, fit){
    plist$object <- fit
    plist$data <- getData(fit)
    utils::capture.output(out <- do.call(unmarked::update, plist))
    out@call$data <- fit@call$data
    out
  }

  cl <- NULL
  if(parallel){
    if(!quiet) cat("Setting up cluster\n")
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(unmarked))
    if(methods::.hasSlot(getData(refModel), "piFun")){
      if(getData(refModel)@piFun %in% ls(envir=.GlobalEnv)){
        parallel::clusterExport(cl, getData(refModel)@piFun) 
      }
    }
  }
  
  if(!quiet) cat("Fitting", length(plist), "models\n")
  
  mods <- pbapply::pblapply(plist, new_fit, refModel, cl=cl)
  names(mods) <- sapply(mods, modName)
  unmarked::fitList(fits = mods)
}

#"dredge" method for fitting all subsets of a global model---------------------

setGeneric("dredge", function(globalModel, quiet=FALSE, parallel=FALSE, ...)
           standardGeneric("dredge"))

setMethod("dredge", "unmarkedFit", 
          function(globalModel, quiet=FALSE, parallel=FALSE, ...){
  fc <- formulaCombs(globalModel@formula)
  fitModels(globalModel, formula=fc, quiet=quiet, parallel=parallel)
})

.dredge_gmods <- function(globalModel, quiet=FALSE, parallel=FALSE, ...){
  f <- globalModel@formlist
  fc <- formulaCombs(lf=f$lambdaformula, phf=f$phiformula, pf=f$pformula)
  fitModels(globalModel, lambdaformula=fc$lf, phiformula=fc$phf,
            pformula=fc$pf, quiet=quiet, parallel=parallel)
}

setMethod("dredge", "unmarkedFitGDS", .dredge_gmods)
setMethod("dredge", "unmarkedFitGMM", .dredge_gmods)
setMethod("dredge", "unmarkedFitGPC", .dredge_gmods)

#Currently broken because fitList doesn't work with occuFP
setMethod("dredge", "unmarkedFitOccuFP", 
          function(globalModel, quiet=FALSE, parallel=FALSE, ...){
  fc <- formulaCombs(sf=globalModel@stateformula, df=globalModel@detformula,
                     fpf=globalModel@FPformula, bf=globalModel@Bformula)
  fitModels(globalModel, detformula=fc$df, FPformula=fc$fpf, Bformula=fc$bf,
            stateformula=fc$sf, quiet=quiet, parallel=parallel)
})

setMethod("dredge", "unmarkedFitColExt",
          function(globalModel, quiet=FALSE, parallel=FALSE, ...){
  fc <- formulaCombs(pf=globalModel@psiformula, gf=globalModel@gamformula,
                     ef=globalModel@epsformula, df=globalModel@detformula)
  fitModels(globalModel, pformula=fc$df, gammaformula=fc$gf, epsilonformula=fc$ef,
            psiformula=fc$pf, quiet=quiet, parallel=parallel)
})

setMethod("dredge", "unmarkedFitPCO",
          function(globalModel, quiet=FALSE, parallel=FALSE, ...){
  f <- globalModel@formlist
  fc <- formulaCombs(lf=f$lambdaformula, gf=f$gammaformula,
                     of=f$omegaformula, pf=f$pformula, iof=f$iotaformula)
  
  if(is.null(globalModel@call[["iotaformula"]])){
    fitModels(globalModel, lambdaformula.=fc$lf, gammaformula.=fc$gf,
              omegaformula.=fc$of, pformula.=fc$pf,
              quiet=quiet, parallel=parallel)
  } else {
    fitModels(globalModel, lambdaformula.=fc$lf, gammaformula.=fc$gf,
            omegaformula.=fc$of, pformula.=fc$pf, iotaformula.=fc$iof,
            quiet=quiet, parallel=parallel)
  }
})

setMethod("dredge", "unmarkedFitOccuTTD",
          function(globalModel, quiet=FALSE, parallel=FALSE, ...){
  fc <- formulaCombs(pf=globalModel@psiformula, gf=globalModel@gamformula,
                     ef=globalModel@epsformula, df=globalModel@detformula)
  fitModels(globalModel, psiformula=fc$pf, gammaformula=fc$gf,
            epsilonformula=fc$ef, detformula=fc$df, 
            quiet=quiet, parallel=parallel)
})

#Unsupported fit types---------------------------------------------------------

.dredge_nosupport <- function(globalModel, quiet=FALSE, parallel=FALSE, ...){
  stop(paste("dredge does not yet support", class(globalModel)[1]))
}

#May never be practical due to complex formula structure
setMethod("dredge", "unmarkedFitOccuMulti", .dredge_nosupport)
setMethod("dredge", "unmarkedFitOccuMS", .dredge_nosupport)
