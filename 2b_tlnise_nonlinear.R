###################
### To organize poisson regression results with nonlinear temporal trend
### and to apply TLNISE to aggregate to state level
### written by Chen Chen
### Note: search for "please update before run" to make edits where necessary;
###################

## Installation of tlnise package
## Here is a link to the package: https://github.com/rdpeng/tlnise
## You will also need a Fortran compiler installed to compile the Fortran code
# install.packages("remotes")
# remotes::install_github("rdpeng/tlnise")
library(tlnise)
library(MASS)
library(splines)

setwd("") ## please update before run
years <- 2006:2019

## create a directory for storing organized estimates from first-stage model
outdir  <- "results/estimates"    
if (!file.exists(outdir)) dir.create(outdir)
## create a directory for storing estimates from second-stage model
outdir2  <- "results/estimates/TLNISE"    
if (!file.exists(outdir2)) dir.create(outdir2)
## create a directory for storing organized estimates from second-stage model
outdir3  <- "results/summary"    
if (!file.exists(outdir3)) dir.create(outdir3)

## read in list of zctas
siteList <- read.csv(file.path("data", "Countylist_50yearp5w_192_clean.csv"),
                     colClasses=c(zcta5="character"))

# Default computes computes effect for moving average exposures of days 0-2, no regional analysis
##################
cmdOpts <- c(
  'lagtest',  ## Sensititivty: 0-3 single lags, m01 and m02
  'region' ## whether to include analysis for regions
)

## Process any command line options
opts <- cmdOpts %in% commandArgs()
names(opts) <- cmdOpts

stopifnot(length(cmdOpts) == length(opts))

## Create some global variables w.r.t. the command line arguments
for(i in seq(along = cmdOpts)) {
  assign(cmdOpts[i], opts[i], globalenv())
}
##################

## Set the stage
## need to specify whether to run multiple lags or regional analysis, 
## and which subgroup to analyze
##################
## please update before run
lagtest <- T ## calcualte single lags, m01 and m02
region <- T ## include analysis for regions (ie, with and without intense wildfire)

## please update before run
## can only analyze one subgroup in a big loop due to data structure
## please only leave one subgroup here
nm <- c ( 
  # "age65up"# look at age 65 and older
  # "age15to64"
  # "age0to14"
  # "female"
  # "male"
  # "white"
  # "black"
  # "hispanic"
  # "asian"
  # "native"
  # "other"
  "" # look at total group
) 

## End of updates before run


outcomes <- c(  # Outcomes to consider
  "circulatory", 
  "respiratory"
)

yr <- c(  # decide which years to analyze
  "06-19"
)

if (lagtest) {
  n.lag <- c("m02mi", "m01mi", 0:2)
} else {
  n.lag <- "m02mi"
}
nltrend.df <- 4 ## since I manually calculated the knots, df=# of knots=terms-2=ns df -1. In theory number of internal knots = df - 1, terms = df + 1
nltrend.knots <- lapply(nltrend.df, function(a) {
  ceiling(quantile(0:13, probs=(1:a)/(a+1)))
})

bar <- lapply(nltrend.df, function(a) {
  lapply(n.lag, function(b) {
    paste0("nonlinear", "-lag", b, "-df", a)
  })
})

tt <- rep("nltrend", 2)

tt[2] <- paste0(tt[2], "df13")
dfVec <- 13

regions <- if (region) {
  c("State",
    "SomeWFabove35", "NoWFabove35"
  )
} else {
  c("State")
}
##################

#Orgnize first level results -- non linear trend
##################
ncounty <- numeric()
est.title <- numeric()

for (i in 1:length(nltrend.df)) {
  lags <- unlist(bar[i])
  out<-as.numeric()
  outfal<-numeric()
  for(nm_ in nm){
    for(outcome in outcomes) {
      for(lag_ in lags){
        for(yr_ in yr){
          for(df_ in dfVec){  #}}}}}
            if (nm_=="") {
              results <- readRDS(file.path("results", "raw-results", 
                                           paste0("results.", tt[1], ".", outcome,
                                                  ".", yr_, ".rds")))
            } else {
              results <- readRDS(file.path("results", "raw-results", 
                                           paste0("results.", tt[1], ".", outcome, "_",
                                                  nm_, ".", yr_, ".rds")))
            }
            
            lag1 <- lapply(results[[lag_]], "[[", paste("df",df_,sep=""))
            lagfal<-lag1[sapply(lag1, inherits, what = "condition")]
            lag1 <- lag1[!sapply(lag1, inherits, what = "condition")]
            
            ## Make a table
            loc.suc <- which(sapply(lag1, length)==3)
            loc.sing <- which(sapply(lag1, length)==4) # counties with singularity
            tab <- local({
              d <- as.data.frame(do.call("rbind", lapply(lag1[loc.suc], unlist)))
              nbeta <- grep("beta", names(d))
              names(d)[nbeta] <- paste("Beta", nbeta-min(nbeta), sep="")
              d$zcta5 <- row.names(d)
              d$outcome <- outcome
              d$dfTime <- df_
              d$lag <- lag_
              d$subgroup <- nm_
              d$period <- yr_
              d$singularity <- NA
              d
            })
            if (length(loc.sing) != 0) {
              tab.sing <- local({
                d <- as.data.frame(do.call("rbind", lapply(lag1[loc.sing], function(a) {
                  unlist(a[1:3])
                })))
                nbeta <- grep("beta", names(d))
                names(d)[nbeta] <- paste("Beta", nbeta-min(nbeta), sep="")
                d$zcta5 <- row.names(d)
                d$outcome <- outcome
                d$dfTime <- df_
                d$lag <- lag_
                d$subgroup <- nm_
                d$period <- yr_
                d$singularity <- sapply(lag1[loc.sing], function(a) {
                  unlist(a[c(4)])
                })
                d
              }) 
            } else {
              tab.sing <- numeric(0)
            }
            out<-rbind(out, tab, tab.sing)
            if (length(lagfal)!=0){
              tabfal <- local({
                d <- as.data.frame(do.call("rbind", lapply(lagfal, unlist)))
                d$zcta5 <- row.names(d)
                d$outcome <- outcome
                d$dfTime<-df_
                d$lag<-lag_
                d$subgroup<-nm_
                d$period<-yr_
                d
              })
              outfal<-rbind(outfal,tabfal)
            }
          }
        }
      }
    }
  }
  nc <- NULL
  write.csv(out, row.names = F,
            file = paste(outdir, "/",  nm, "_", "df", nltrend.df[i], "_", "estimates_", tt[2], "_full.csv", sep=""))
  if (length(outfal)!=0){
    write.csv(outfal, row.names = F,
              file = paste(outdir, "/",  nm, "_", "df", nltrend.df[i], "_", "estimates_", tt[2], "_fails.csv", sep=""))
    outall <- out[which(!out$zcta5 %in% outfal$zcta5), ]
    nc <- length(unique(c(out$zcta5, outfal$zcta5)))
  } else {
    outall <- out 
    nc <- length(unique(out$zcta5))
  }
  write.csv(outall, row.names = F,
            file = paste(outdir, "/",  nm, "_", "df", nltrend.df[i], "_", "estimates_", tt[2], "_completeset.csv", sep=""))
  results <- merge(outall, siteList, by="zcta5", all.x=TRUE)
  ncounty <- c(ncounty, length(unique(results$zcta5)))
  est.title <- c(est.title, paste0(outdir, "/",  nm, "_", "df", nltrend.df[i], "_", 
                                   "estimates_", tt[2], "_50yearp5w_", nc, "_", ncounty[i], ".csv"))
  write.csv(results, est.title[i], row.names=FALSE)
}
ncounty
est.title
##################

##
#RunTLNISE and orgnize results -- nltrend
##################
## assign zctas to different regions (eg, with and without intense wildfire)
meta <- siteList
meta$region <- "State"
meta <- rbind(meta, siteList)

for (j in 1:length(nltrend.df)) {
  lags <- unlist(bar[j])
  ## load estimates of first-stage model parameters
  results <- read.csv(file = est.title[j],
                      stringsAsFactors = F, 
                      colClasses = c(zcta5="character"))
  
  ### estimation of second-stage model parameters
  for(nm_ in nm){
    for(outcome_ in outcomes){   
      for(lag_ in lags){
        for(df_ in dfVec){
          for(yr_ in yr){
            for (region_ in regions){ #}}}}}}
              filename <- paste(ncounty[j], nm_,outcome_,yr_,region_,lag_,df_,"rds",sep=".")
              outfile <- file.path(outdir2, filename)
              
              zcta5.include <- meta[which(meta$region==region_),"zcta5"]
              all_results <- data.frame(results, include = 0)
              all_results$include[all_results$zcta5%in%zcta5.include] <-1
              loc <- grep("Beta|var|zcta5", names(all_results))
              
              ## run tlnise and if it does not produce an estimate,
              ## remove the county with the largest standard error,
              ## then repeat call to tlnise with reduced data set
              
              params <- all_results[all_results$outcome == outcome_ &
                                      all_results$lag == lag_ &
                                      all_results$dfTime == df_ &
                                      all_results$period == yr_ &
                                      # all_results$subgroup == nm_ &
                                      all_results$include == 1, loc]
              
              params <- na.omit(params)                      ## rm.na for tlnise
              loc.beta <- grep("Beta", names(params))
              loc.var <- grep("var", names(params))
              params <- params[order(params[, max(loc.var)],decreasing=T),]       ## sort by variance
              
              estimate <- NULL
              
              estimate <- tryCatch(
                { f <- tlnise(Y=params[, loc.beta], 
                              V=array(unlist(lapply(1:nrow(params), function(a) {
                                params[a, loc.var]})), 
                                c(sqrt(length(loc.var)), sqrt(length(loc.var)), dim(params)[1])),
                              w=rep(1, nrow(params)),
                              seed=1234, maxiter=5000, prnt = FALSE
                )
                f
                },
                condition = function(cond) {
                  cat("\t", paste(ncounty[j], nm_,outcome_,yr_,region_,lag_,df_, sep="."))
                  cat("\t", as.character(cond))
                  cond
                }
              ) ## end tryCatch
              
              if (!inherits(estimate, what = "condition")) {
                gs <- tlnise:::sampleGamma(estimate,                                          
                                           V=array(unlist(lapply(1:nrow(params), function(a) {
                                             params[a, loc.var]})),
                                             c(sqrt(length(loc.var)), sqrt(length(loc.var)), dim(params)[1])),
                                           as.matrix(params[, loc.beta])) ## generate samples of gammastar, Dstar from equation (14) of Everson & Morris
                sampsGamma <- tryCatch({
                  set.seed(1234)
                  tlnise:::drawSample0(gs, n=10000) 
                },
                condition = function(cond) {
                  cat("\t", paste(ncounty[j], nm_,outcome_,yr_,region_,lag_,df_, sep="."))
                  cat("\t", as.character(cond))
                  cond$call <- NULL
                  cond
                }) ## end tryCatch
                
                estimate <- list( gamma=estimate$gamma, theta=estimate$theta,
                                  SDtheta=estimate$SDtheta, A=estimate$A,
                                  rtA=estimate$rtA, Dgamma=estimate$Dgamma,
                                  Vtheta=estimate$Vtheta, B0=estimate$B0,
                                  lr=estimate$lr, zcta5.used= as.character(params$zcta5), 
                                  beta=params[, loc.beta], var=params[, loc.var],
                                  sampsGamma=sampsGamma
                )
              }
              
              ## save estimates in external file
              saveRDS(estimate,file=outfile)
            } ## end for region_
          } ## end for yr_
        } ## end for df_
      } ## end for lag_
    } ## end for outcome_
  } ## end for nm_
} ## end for different df of ns, different length for beta
##################


## create a directory for storing estimates from second-stage model
##################
## retrieve second-stage estimates from previous step
for (j in 1:length(nltrend.df)) {
  lags <- unlist(bar[j])
  outcome <- lag <- df <- years <- NULL
  pooled.est <- pooled.var <- NULL
  n <- NULL
  region <- subs <- NULL
  change.est <- change.lb <- change.ub <- dif.year <- NULL
  samps <- NULL
  out <- outfal <- NULL
  est_st <- est_st.lb <- est_st.ub <- NULL
  est_end <- est_end.lb <- est_end.ub <- NULL
  est_min <- est_min.lb <- est_min.ub <- min.year <- NULL
  est_max <- est_max.lb <- est_max.ub <- max.year <- NULL
  
  for(nm_ in nm){
    for(outcome_ in outcomes){   
      for(lag_ in lags){
        for(df_ in dfVec){
          for(yr_ in yr){
            for (region_ in regions){  #}}}}}}
              filename <- paste(ncounty[j], nm_,outcome_,yr_,region_,lag_,df_,"rds",sep=".")
              infile <- file.path(outdir2, filename)
              estimate <- readRDS(infile)
              
              if (inherits(estimate, what = "condition")) {
                d <- data.frame(Message=as.character(estimate), 
                                subgroup=nm_, outcome=outcome_, years=yr_, 
                                regions=region_, lag=lag_, df=df_)
                outfal <- rbind(outfal, d)
              } else {
                if (inherits(estimate$sampsGamma, what = "condition")) {
                  d <- data.frame(Message=paste("Samps error", as.character(estimate$sampsGamma)), 
                                  subgroup=nm_, outcome=outcome_, years=yr_, 
                                  regions=region_, lag=lag_, df=df_)
                  outfal <- rbind(outfal, d)
                } else {
                  pooled.est <- rbind(pooled.est, estimate$gamma[, 1])
                  pooled.var <- rbind(pooled.var, as.vector(estimate$Dgamma))
                  n <- c(n, length(estimate$zcta5.used))
                  outcome <- c(outcome, outcome_)
                  lag <- c(lag, lag_)
                  df <- c(df, df_)
                  years <- c(years, yr_)
                  region <- c(region, region_) 
                  subs <- c(subs, nm_)
                  yrstart <- ifelse(substring(yr_, 1, 2)=="99",
                                    as.numeric(paste("19", substring(yr_, 1, 2), sep="")),
                                    as.numeric(paste("20", substring(yr_, 1, 2), sep="")))
                  yrend <- ifelse(substring(yr_, 4, 5)=="99",
                                  as.numeric(paste("19", substring(yr_, 4, 5), sep="")),
                                  as.numeric(paste("20", substring(yr_, 4, 5), sep="")))
                  newx <- seq(0, yrend-yrstart)
                  XX <- cbind(1, ns(newx, knots=nltrend.knots[[j]]))
                  samps <- estimate$sampsGamma %*% t(XX)
                  
                  min.year <- c(min.year, yrstart-1+which.min(apply(samps, 2, mean)))
                  max.year <- c(max.year, yrstart-1+which.max(apply(samps, 2, mean)))
                  dif.year <- c(dif.year,paste(yrstart-1+which.min(apply(samps, 2, mean)),
                                               yrstart-1+which.max(apply(samps, 2, mean)), sep="-"))
                  
                  ## PM2.5
                  diffrisk.samps <- -100 * (exp(10*samps[,which.min(apply(samps, 2, mean))])-1) +
                    100 * (exp(10*samps[, which.max(apply(samps, 2, mean))])-1)

                  change.est <- c(change.est, mean(diffrisk.samps))
                  change.lb <- c(change.lb, quantile(diffrisk.samps, probs = 0.025))
                  change.ub <- c(change.ub, quantile(diffrisk.samps, probs = 0.975))

                  ## PM2.5
                  samps.st <- 100 * (exp(10*samps[, 1])-1)
                  samps.end <- 100 * (exp(10*samps[, yrend-yrstart+1])-1)
                  samps.min <- 100 * (exp(10*samps[, which.min(apply(samps, 2, mean))])-1)
                  samps.max <- 100 * (exp(10*samps[, which.max(apply(samps, 2, mean))])-1)
                  
                  est_st <- c(est_st, mean(samps.st))
                  est_st.lb <- c(est_st.lb, quantile(samps.st, probs = 0.025))
                  est_st.ub <- c(est_st.ub, quantile(samps.st, probs = 0.975))
                  est_end <- c(est_end, mean(samps.end))
                  est_end.lb <- c(est_end.lb, quantile(samps.end, probs = 0.025))
                  est_end.ub <- c(est_end.ub, quantile(samps.end, probs = 0.975))
                 est_min <- c(est_min, mean(samps.min))
                  est_min.lb <- c(est_min.lb, quantile(samps.min, probs = 0.025))
                  est_min.ub <- c(est_min.ub, quantile(samps.min, probs = 0.975))
                  est_max <- c(est_max, mean(samps.max))
                  est_max.lb <- c(est_max.lb, quantile(samps.max, probs = 0.025))
                  est_max.ub <- c(est_max.ub, quantile(samps.max, probs = 0.975))
                  
                } ## end for samps fail or success
              } ## end for estimate fail or success
            } ## end for region_
          } ## end for yr_
        } ## end for df_
      } ## end for lag_
    } ## end for outcome_
  } ## end for nm_
  
  out <- data.frame(pooled.est, pooled.var, 
                    subgroup=subs, outcome=outcome, years=years, 
                    regions=region, lag=lag, df=df, zcta5.used=n, 
                    min.year=min.year, max.year=max.year,
                    dif.year=dif.year, change.est=change.est, change.lb=change.lb, change.ub=change.ub,
                    est_st = est_st, est_st.lb = est_st.lb, est_st.ub = est_st.ub,
                    est_end = est_end, est_end.lb = est_end.lb, est_end.ub = est_end.ub,
                    est_min = est_min, est_min.lb = est_min.lb, est_min.ub = est_min.ub,
                    est_max = est_max, est_max.lb = est_max.lb, est_max.ub = est_max.ub)
  nbeta <- ncol(pooled.est)
  nvar <- ncol(pooled.var)
  names(out)[1:(nbeta+nvar)] <- c(paste("beta", 1:nbeta-1, sep=""), 
                                  paste0("cov", rep(1:sqrt(nvar)-1, times=sqrt(nvar)), 
                                         rep(1:sqrt(nvar)-1, each=sqrt(nvar))))
  write.csv(out, paste(outdir3, "/", ncounty[j], nm, "_", "estimates", nltrend.df[j], "df.", tt[2], ".state.csv", sep=""), row.names=FALSE)
  if (!is.null(outfal)) {
    write.csv(outfal, paste(outdir3, "/", "Fail", ncounty[j], nm, "_", "estimates", nltrend.df[j], "df.", tt[2], ".state.csv", sep=""), row.names=FALSE)
  }
}
####################