###################
### To organize poisson regression results without considering temporal trend
### and to apply TLNISE to aggregate to state level
### written by Chen Chen 
###################
## Installation of tlnise--might require installation of rtools as well
# install.packages("remotes")
# remotes::install_github("rdpeng/tlnise")
library(tlnise)

setwd("")
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

siteList <- read.csv(file.path("data", "Countylist_50yearp5w_192_clean.csv"),
                     colClasses=c(zcta5="character"))

# Default computes df=8 no subgroup analysis with 2 regions, no comparison
##################
cmdOpts <- c(
  'lagtest', ## whether to calcualte single lags and m01
  'comp2p', ## decide whether to compare two periods
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

# Set the stage
##################
lagtest <- T
comp2p <- c("06-12", "13-19")
region <- T #I have not incorporated this feature into analysis 102720

nm <- c ( # can only analyze one subgroup at a time due to data structure
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

yr <- c(  # decide which years to analyze
  "06-19",
  "06-12",
  "13-19"
)

outcomes <- c(  # Outcomes to consider
  "circulatory", 
  "respiratory"
)

if (lagtest) {
  lags <- paste0("lag", c("m02mi", "m01mi", 0:2))
} else {
  lags <- "lagm02mi"
}

tt <- rep("", 2)

if (lagtest) tt[2] <- paste0("lagtest.", tt[2])

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

#Orgnize first level results -- basic model
################################
out <- numeric()
outfal <- numeric()
for(nm_ in nm){
  for(outcome in outcomes) {
    for(lag_ in lags){
      for(yr_ in yr){
        for(df_ in dfVec){
          if (nm_=="") {
            results <- readRDS(file.path("results", "raw-results", 
                                         paste0("results.", tt[1], outcome,
                                                ".", yr_, ".rds")))
          } else {
            results <- readRDS(file.path("results", "raw-results", 
                                         paste0("results.", tt[1], outcome, "_", 
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
          file =  paste0(outdir, "/", tt[2], nm, "_", "estimates_full.csv"))
if (length(outfal)!=0){
  write.csv(outfal, row.names = F,
            file =  paste0(outdir, "/", tt[2], nm, "_", "estimates_fails.csv"))
  outall <- out[which(!out$zcta5 %in% outfal$zcta5), ]
  nc <- length(unique(c(out$zcta5, outfal$zcta5)))
} else {
  outall <- out 
  nc <- length(unique(out$zcta5))
}
write.csv(outall, row.names = F,
          file =  paste0(outdir, "/", tt[2], nm, "_", "estimates_completeset.csv"))
results <- merge(outall, siteList, by="zcta5", all.x=TRUE)
ncounty <- length(unique(results$zcta5))
est.title <- paste0(outdir, "/", tt[2], nm, "_", "estimates_50yearp5w_", nc, "_",
                    ncounty,".csv")
write.csv(results, est.title, row.names=FALSE)
ncounty
est.title
##################

##
#RunTLNISE and orgnize results
##################
results <- read.csv(est.title, colClasses = c(zcta5="character"))

## load zcta5s to include and regions of each county--need to update if region is considered 102720
meta <- siteList
meta$region <- "State"
meta <- rbind(meta, siteList)

### estimation of second-stage model parameters
for(nm_ in nm){
  for(outcome_ in outcomes){   
    for(lag_ in lags){
      for(df_ in dfVec){
        for(yr_ in yr){
          for (region_ in regions){
            filename <- paste(ncounty, nm_,outcome_,yr_,region_,lag_,df_,"rds",sep=".")
            outfile <- file.path(outdir2, filename)
            
            zcta5.include<-meta[which(meta$region==region_),"zcta5"]
            all_results <- data.frame(results, include = 0)
            all_results$include[all_results$zcta5%in%zcta5.include] <-1
            loc <- grep("beta|var|zcta5", names(all_results))
            
            params <- all_results[all_results$outcome == outcome_ &
                                    all_results$lag == lag_ &
                                    all_results$dfTime == df_ &
                                    all_results$period == yr_ &
                                    all_results$include == 1,
                                  loc]
            
            params <- na.omit(params)                      ## rm.na for tlnise
            
            loc.beta <- grep("beta", names(params))
            loc.var <- grep("var", names(params))
            params <- params[order(params[, max(loc.var)],decreasing=T),]       ## sort by variance
            
            estimate <- NULL
            
            estimate <- tryCatch(
              { f <- tlnise(Y=params[, loc.beta], 
                            V=array(unlist(lapply(1:nrow(params), function(a) {
                              params[a, loc.var]})), 
                              c(sqrt(length(loc.var)), sqrt(length(loc.var)), dim(params)[1])),
                            w=rep(1, nrow(params)),
                            seed=1234, prnt = FALSE
              )
              f
              },
              condition = function(cond) {
                cat("\t", paste(ncounty, nm_,outcome_,yr_,region_,lag_,df_, sep="."))
                cat("\t", as.character(cond), "\n")
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
                                beta=params$beta, var=params$var,
                                sampsGamma=sampsGamma)  
            }
            ## save estimates in external file
            saveRDS(estimate, file=outfile)
          } ## end for region_
        } ## end for yr_
      } ## end for df_
    } ## end for lag_
  } ## end for outcome_
} ## end for nm_
##################

## create a directory for storing estimates from second-stage model
## I calculated intervals percentage change in risk of HA with 10ug/m^3 increase 
## in PM2.5 using random samples (PI)
##################
## retrieve second-stage estimates from previous step
outcome <- region <- NULL
lag <- df <- years <- n <- NULL
pooled.est <- pooled.se <- NULL
dif2p <- NULL
est <- est.lb <- est.ub <- NULL

for(nm_ in nm){
  for(lag_ in lags){
    for(df_ in dfVec){
      for(outcome_ in outcomes){   
        for (region_ in regions){
          comp <- difsamps <- samps <- NULL
          for(yr_ in yr){
            filename <- paste(ncounty, nm_,outcome_,yr_,region_,lag_,df_,"rds",sep=".")
            infile <- file.path(outdir2, filename)
            estimate<-readRDS(infile)
            samps <- estimate$sampsGamma
            if (!is.null(comp2p)){
              comp <- c(comp, structure(list(samps), names=yr_))
            }
            samps <- 100*(exp(10*samps)-1)
            est <- c(est, mean(samps))
            est.lb <- c(est.lb, quantile(samps, probs = 0.025))
            est.ub <- c(est.ub, quantile(samps, probs = 0.975))
            
            pooled.est <- c(pooled.est, estimate$gamma[1])
            pooled.se <- c(pooled.se, estimate$gamma[2])
            n <- c(n, length(estimate$zcta5.used))
            outcome <- c(outcome, outcome_)
            lag <- c(lag, lag_)
            df <- c(df, df_)
            years <- c(years, yr_)
            region <- c(region, region_)
          } ## end for yr_
          if (!is.null(comp2p)) {
            difsamps <- 100*(exp(10*(comp[[comp2p[2]]]))-1) - 
              100*(exp(10*(comp[[comp2p[1]]]))-1)
            temp <- data.frame(period=names(comp), setNames(replicate(3, NA, simplify = F), 
                                                            c("dif2p", "dif2p.lb", "dif2p.ub")))
            temp$dif2p[temp$period %in% comp2p] <- mean(difsamps)
            temp$dif2p.lb[temp$period %in% comp2p] <- quantile(difsamps, probs = 0.025)
            temp$dif2p.ub[temp$period %in% comp2p] <- quantile(difsamps, probs = 0.975)
            dif2p <- rbind(dif2p, temp[, c("dif2p", "dif2p.lb", "dif2p.ub")])
          }
        } ## end for region_
      } ## end for df_
    } ## end for lag_
  } ## end for outcome_
} ## end for nm_

out <- data.frame(pooled.est, pooled.se, subgroup=nm, outcome=outcome, years=years, 
                  regions=region, lag=lag, df=df, zcta5.used=n, 
                  dif2p, est, est.lb, est.ub)

write.csv(out, paste0(outdir3, "/", ncounty, tt[2], nm, "_", "estimates.state.csv"), 
          row.names=FALSE)
#####################


