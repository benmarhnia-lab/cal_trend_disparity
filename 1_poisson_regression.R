##################
# Time-series trend first stage model for PM2.5 in California
# written by Chen Chen, based on Roger Peng's codes for NMMAPS
##################

library(zoo)
library(tsModel)
library(splines)
library(data.table)

## read in site list
##################
setwd("")
dataset <- "data/mock_data"

years <- 2006:2019

siteList <- read.csv(file.path("data", "Countylist_50yearp5w_192_clean.csv"),
                     colClasses=c(zcta5="character"))
##################

# Default computes national basic estimate for 2006-2019
##################
cmdOpts <- c(
  'nl.trend', ## Sensitivity: Non-linear PM2.5-HA trend with internal knot=2 at quantile locations
  'lagtest' ## estimates with 8 df and 0-2 single lags, and moving average of lag0-1, lag0-2
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

# Set the stage--please update before run
##################
lagtest <- T # whether to run m01, single lag 0-3
nl.trend <- T  # estimate non-linear temporal trend with knots at quantile

nms <- c ( # can only analyze one subgroup at a time due to data structure
  "age65up", # look at age 65 and older
  "age15to64",
  "age0to14",
  "female",
  "male",
  "white",
  "black",
  "hispanic",
  "asian",
  "native",
  "other",
  "" # look at total group
) 

yr <- if (nl.trend) {
  "06-19"
} else {
  c(  # decide which years to analyze
    "06-19",
    "06-12",
    "13-19"
  )
}

outcomes <- c(  # Outcomes to consider
  "circulatory", 
  "respiratory"
)
##################

# Define pollutant variable
##################
if (lagtest) {
  n.lag <- c("m02mi", "m01mi", 0:3)  # test 0 to 2 lags, moving averages of 01 and 02
} else {
  n.lag <- "m02mi"  # use moving average of lag0 to 2
}

dfVec <- if (nl.trend) {
  c(3, 8, 13, 18)
} else {
  13
}

if (nl.trend) {
  nltrend.df <- 1:5 ## since I manually calculated the knots, df=# of knots=terms-2=ns df -1. In theory number of internal knots = df - 1, terms = df + 1
  nltrend.knots <- sapply(nltrend.df, function(a) {
    paste(ceiling(quantile(0:13, probs=(1:a)/(a+1))), collapse = ",")
  })
}

poll.basic <- if (lagtest) {
  structure(c("mean_pm25m02mi", "mean_pm25m01mi", 
              "mean_pm25", paste0("Lag(mean_pm25, ", 1:3, ")")),
            names=paste0("lag", n.lag))
} else {
  structure("mean_pm25m02mi", names=paste0("lag", n.lag))
}

polls <- if (nl.trend) {
  structure(unlist(lapply(poll.basic, function(a) {
    lapply(nltrend.knots, function(b) {
      c(a, paste(a, ":ns(year, knots=c(", b,"))", sep=""))
    })
  }), recursive = FALSE), 
  names = unlist(lapply(names(poll.basic), function(a) {
    lapply(nltrend.df, function(b) {
      paste0("nonlinear", "-", a, "-df", b)
    })
  })))
} else {
  structure(as.list(poll.basic),
            names =names(poll.basic))
}
##################

# Prepare folders
##################
## Then we specify the directory for storing the outcome specific
## results.  We store a separate output file for each outcome---the file
## will have the results for all the pollutant Lags, degrees of freedom,
## and counties.  If the directory doesn't exist it will be created.

## Create output directory if it doesn't already exist
outdir <- "results/raw-results"

if (!file.exists(outdir)) dir.create(outdir)

## Cache directory
if (!file.exists(file.path(outdir, "cache"))) dir.create(file.path(outdir, "cache"))

makeCacheDir <- function(site) {
  file.path(outdir, "cache", site)
}

makeCacheFilename <- function(siteName, outcome, poll, df.Time, yr_) {
  if(length(poll) > 1) {
    poll <- paste(poll, collapse = "&")
  }
  poll <- gsub(":", "", poll)
  n <- paste("cache", outcome, yr_, gsub(" +", "", poll), siteName,
             formatC(df.Time, flag = "0", width = 2), "rds", sep = ".")
  file.path(makeCacheDir(siteName), n)
}

for(cachedir in makeCacheDir(siteList$zcta5)) {
  if (!file.exists(cachedir)) dir.create(cachedir, showWarnings = FALSE)
}
##################

## Modeling functions
##################
setupFormula <- function(
  outcome,
  pollutant, 
  df.Time, 
  df.Temp, 
  df.Dew
  ) {
  reformulate(c(
    "offset(log(denominator))", 
    pollutant,
    "dow",
    paste("ns(tmpt,", df.Temp, ")"),
    paste("ns(rmtmpt,", df.Temp, ")"),
    paste("ns(dptp,", df.Dew, ")"),
    paste("ns(rmdptp,", df.Dew, ")"),
    paste("ns(date,", df.Time, ")")),
    response = outcome
  )
}

fitSingleSite <- function(
  data,  ## A data frame for a site
  outcome,  ## character, name of outcome
  pollutant,  ## character, name of pollutant variable
  df.Time,  ## df for smooth function of time
  df.Temp = 3,  ## df for temp smooth function
  df.Dew = 3,
  fitModel = TRUE, ...) {
  
  stopifnot(is.character(outcome), is.character(pollutant))
  form <- setupFormula(outcome, pollutant, df.Time, df.Temp, df.Dew)
  
  rval <- if(fitModel)
    glm(form, data = data, family = quasipoisson,
        control = glm.control(maxit = 1000))
  else
    form
  rval
}
##################

# Utilities
##################
preProcess <- function(
  dataraw, 
  siteName,
  denominator, 
  outcome, 
  yrstart, 
  yrend, 
  poll) {
  varList <- c("zcta5", "date", outcome, denominator,
               "tmpt", "rmtmpt", "dptp","rmdptp", 
               "mean_pm25")
  data <- dataraw[dataraw$zcta5==siteName, varList]
  names(data)[4] <- "denominator"
  data$date <- as.Date(data$date)
  data$year <- format(data$date, format="%Y")
  data$dow <- as.factor(format(data$date, "%w")) #reference category is Sunday
  
  if (yrstart!=2006 | yrend!=2019) {
    data <- subset(data, year <= yrend & year >= yrstart)
  }
  
  temp <- unique(data[, c("date", "mean_pm25")])
  temp$mean_pm25m02 <- rollapply(temp$mean_pm25, FUN=mean, fill=NA, align="right", width=3)
  temp$mean_pm25m01 <- rollapply(temp$mean_pm25, FUN=mean, fill=NA, align="right", width=2)
  ## the calculation below ignored missing values in exposure--it will be NA only
  ## when data were missing for 3 or 2 consecutive days
  temp$mean_pm25m02mi <- rollapply(temp$mean_pm25, FUN=mean, fill=NA, align="right", width=3, na.rm=TRUE)
  temp$mean_pm25m01mi <- rollapply(temp$mean_pm25, FUN=mean, fill=NA, align="right", width=2, na.rm=TRUE)
  temp$mean_pm25 <- NULL
  data <- merge(data, temp, by="date")
  
  data$date <- as.numeric(difftime(data$date, as.Date(paste(yrstart, "-01-01", sep="")), units="days"))
  data$year <- as.numeric(data$year)-yrstart
  
  data$tmpt <- data$tmpt-mean(data$tmpt, na.rm = TRUE)
  data$rmtmpt <- data$rmtmpt-mean(data$rmtmpt, na.rm = TRUE)
  data$dptp <- data$dptp-mean(data$dptp, na.rm = TRUE)
  data$rmdptp <- data$rmdptp-mean(data$rmdptp, na.rm = TRUE)
  
  data$dow <- as.factor(data$dow)

  data[order(data$date), ]
}

postProcess <- function(  # Post process glm object; extract various entities
  glmObject, 
  poll) {
  ## Extract coefficients
  cc <- summary(glmObject)$coefficients
  
  ## Extract covariance matrix
  V <- vcov(glmObject)
  
  ## Only keep coefficients corresponding to pollutant variable
  if (seasonal) {
    i <- grep("mean_pm25", rownames(cc), fixed = TRUE)
    j <- grep("mean_pm25", rownames(V), fixed = TRUE)
  } else {
    i <- lapply(poll, function(x) grep(x, rownames(cc), fixed = TRUE))
    i <- unique(unlist(i))
    j <- lapply(poll, function(x) grep(x, rownames(V), fixed = TRUE))
    j <- unique(unlist(j))
  }
  
  if(length(i) == 0) ## No match
    stop("pollutant coefficients not found in model fit")
  
  if(length(i)!=length(j)) ## NA exists in coefficients
    stop("pollutant coefficients with NA")
  
  if (sum(i!=j)==0) {
    rval <- list(beta = cc[i, 1], var = V[i, i], 
                 dispersion=summary(glmObject)$dispersion)
  } else { ## this is to record "singularity"
    rval <- list(beta = cc[i, 1], var = V[j, j], 
                 dispersion=summary(glmObject)$dispersion,
                 condit = paste0("different dimension for coef(", 
                                 paste0(i, collapse = ","), ") and cov(", 
                                 paste0(j, collapse = ","), ")"))
  } 
  
  if(is.null(rval) || length(rval) == 0)
    stop("problem with return value from 'fitSingleSite'")
  
  else
    rval
}
##################

# Running the models:  The big loop
##################
for (nm in nms) {
  for(outcome in outcomes) {  # Cycle over outcomes
    denominator <- ifelse(nm=="", "pop", paste("pop", nm, sep="_"))
    outcome <- ifelse(nm=="", outcome, paste(outcome, nm, sep="_"))
    for (yr_ in yr) {  # Cycle over years
      yrstart <- ifelse(substring(yr_, 1, 2)=="99", 
                        as.numeric(paste("19", substring(yr_, 1, 2), sep="")),
                        as.numeric(paste("20", substring(yr_, 1, 2), sep="")))
      yrend <- ifelse(substring(yr_, 4, 5)=="99", 
                      as.numeric(paste("19", substring(yr_, 4, 5), sep="")),
                      as.numeric(paste("20", substring(yr_, 4, 5), sep="")))
      nyears <- yrend - yrstart + 1
      results <- lapply(polls, function(poll) {  # Cycle over pollutant lags
        out <- lapply(1:nrow(siteList), function(nn) {  # Cycle over counties
          siteName <- siteList$zcta5[nn]
          # Some data frame preprocessing
          data.raw <- readRDS(file.path(dataset, paste0(siteName, ".rds")))
          data <- preProcess(data.raw, siteName, denominator, outcome,
                             yrstart, yrend, poll)
          
          dfVec.use <- dfVec * nyears  # Cycle over degrees of freedom
          
          # Run each model in a tryCatch() --- if anything goes wrong
          # (i.e. a warning or error is thrown), then return a
          # "condition" object and discard the model fit
          out <- lapply(dfVec.use, function(d) {
            # -- Crash Recovery --
            # Cache results in a file.  If the file already
            # exists, just read the cached results.
            cachefile <- makeCacheFilename(siteName, outcome, poll, d/nyears, yr_)
            if(file.exists(cachefile)) {
              message("Using ", cachefile, "\n")
              rval <- readRDS(cachefile)
            } else {
              rval <- tryCatch({
                f <- fitSingleSite(data,
                                   outcome = outcome,
                                   pollutant = poll,
                                   df.Time = d
                )
                postProcess(f, poll)
              }, condition = function(cond) {
                cat("\t", trunc(d / nyears), as.character(cond))
                cond$call <- NULL
                cond
              })
              saveRDS(rval, file = cachefile)
            }
            rval
          })  # End lapply over DFs
          if(length(out) == length(dfVec)) 
            structure(out, names = paste("df", dfVec, sep = ""))
          else
            out
        })  # End cycle over counties
        structure(out, names = siteList$zcta5)
      })  # End cycle over pollutant Lags
      names(results) <- names(polls)
      
      outfile <- if (nl.trend) {
        file.path(outdir, paste("results", "nltrend", outcome, yr_,"rds", sep = "."))
      } else if (dftest) {
        file.path(outdir, paste("results", "dftest", outcome, yr_,"rds", sep = "."))
      } else {
        file.path(outdir, paste("results", outcome, yr_,"rds", sep = "."))
      }
      message("Saving outcome results....\n")
      print(outfile)
      saveRDS(results, file = outfile, compress = TRUE)
    }
  }  ## End cycle over outcomes
} ## End cycle over nms

##################