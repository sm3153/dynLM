#########################################################################
#### BUILDING SUPER DATASET ####
#########################################################################

# -----------------------------------------------------------------------
# cutLMsuper: Build stacked super dataset from original dataset
#              (which can be either in wide or long format - see details in cutLM)
# -----------------------------------------------------------------------
# Input:
# - data             : Data frame from which to construct landmark super dataset
# - outcome          : List with items time and status, containing character strings 
#                      identifying the names of time and status variables, respectively, 
#                      of the survival outcome
# - LMs              : vector, the value of the landmark time points
# - w          : Scalar, the value of the window TODO describe better
# - covs             : List with items fixed and varying, containing character strings 
#                      specifying column names in the data containing time-fixed and 
#                      time-varying covariates, respectively
# - format	         : Character string specifying whether the original data are in wide (default) or in long format
# - id	             : Character string specifying the column name in data containing the subject id; only needed if format="long"
# - rtime	           : Character string specifying the column name in data containing the (running) time variable associated; only needed if format="long"
# - right	           : Boolean (default=TRUE), indicating if the intervals for the time-varying covariates are closed on the 
#                      right (and open on the left) or vice versa, see cut
# -----------------------------------------------------------------------
# Output: 
# LMdata             : An object of class "LM.data.frame". This has various attributes
# - $LMdata: containing the stacked data set, i.e., the outcome and the values of time-fixed and time-varying covariates taken at the landmark time points. The value of the landmark time point is stored in column LM.
# - $outcome: same as input
# - $w: same as input
# -----------------------------------------------------------------------
cutLMsuper <- function(data, outcome, LMs, w, covs, format = c("wide", "long"), id, rtime, right=T){
  if (format == "wide"){
    LMdata <- cutLM(data=data,
                    outcome=outcome,
                    LM=LMs[1],
                    horizon=LMs[1]+w,
                    covs=covs, 
                    format="wide",
                    right=right)
    if (length(LMs) > 1){
      for (i in 2:length(LMs))
        LMdata <- rbind(LMdata,cutLM(data=data,
                                     outcome=outcome,
                                     LM=LMs[i],
                                     horizon=LMs[i]+w,
                                     covs=covs, 
                                     format="wide",
                                     right=right))
    }
    
  } else if (format == "long"){
    LMdata <- cutLM(data=data,
                    outcome=outcome,
                    LM=LMs[1],
                    horizon=LMs[1]+w,
                    covs=covs, 
                    format, id, rtime, right)
    if (length(LMs) > 1){
      for (i in 2:length(LMs))
        LMdata <- rbind(LMdata,cutLM(data=data,
                                     outcome=outcome,
                                     LM=LMs[i],
                                     horizon=LMs[i]+w,
                                     covs=covs, 
                                     format, id, rtime, right))
    }
  }
  out=list(LMdata=LMdata, outcome=outcome, w=w)
  class(out)="LM.data.frame"
  return(out)
}


# -----------------------------------------------------------------------
# addLMtime : Add LM-time interations to a dataset that contains a column LM
# TODO: add _1, _2 explanation and attributes
# -----------------------------------------------------------------------
# Input:
# - LMdata           : An object of class "LM.data.frame", this can be created by running cutLMsuper, or creating a stacked data set and storing it in a list with attributes outcome and w
# - LMcovars         : List of ovariates that are to have a LM interaction
# - func_covars      : A list of functions to use for interactions between LMs and covariates. If fitting a coxph model, the list has length 1, e.g. list(c(f1,f2,f3)). If fitting a CSC model, the list can have length 1 or length=number of causes, if different interactions are desired for different causes. 
# - func_LMs         : A list of functions to use for transformations of the LMs. Its form is analogous to func_covars. 
# -----------------------------------------------------------------------
# Output:
# LMdata             :  An object of class "LM.data.frame" which has the following components: 
# - LMdata: The LM super dataset which now also contains LM time-interactions.
# - w, outcome: (already)
# - func_covars: as the input
# - func_LMs: as the input
# - LMcovars: as the input
# - allLMcovars: a list of covariates that include LM-time interactions, i.e. if cov was in LMcovars, allLMcovars will contain cov_1, cov_2, ..., cov_i if there are i func_covars interactions
# -----------------------------------------------------------------------
addLMtime <- function(LMdata, LMcovars, func_covars=NULL, func_LMs=NULL){
  data <- LMdata$LMdata
  if (is.null(func_covars)){
    # f gives covariate-time interactions
    f1 <- function(t) 1
    f2 <- function(t) t
    f3 <- function(t) t^2
    func_covars <- list(c(f1,f2,f3))
  }
  if (is.null(func_LMs)){
    # g lets the hazard depend on time
    g1 <- function(t) f2(t)
    g2 <- function(t) f3(t)
    func_LMs <- list(c(g1,g2))
  }
  
  allLMcovars <- c()
  # Add func_covarss: covariate LM interactions
  for(i in 1:length(pred_covars)){
    for (j in 1:length(func_covars[[1]])){
      f <- func_covars[[1]][[j]]
      name <- paste(pred_covars[i],"_",j,sep="")
      data[[name]]  <- data[[pred_covars[i]]]*f(data$LM)
      allLMcovars <- c(allLMcovars, name)
    }
  }
  # Add func_LMs: LM interactions
  for (k in 1:length(func_LMs[[1]])){
    g <- func_LMs[[1]][[k]]
    name <- paste("LM_",k,sep="")
    data[[name]]  <- g(data$LM)
    allLMcovars <- c(allLMcovars, name)
  }
  
  data <- data %>% select(-all_of(pred_covars))
  LMdata$LMdata <- data
  
  LMdata$func_covars <- func_covars
  LMdata$func_LMs <- func_LMs
  LMdata$LMcovars <- LMcovars
  LMdata$allLMcovars <- allLMcovars
  
  return(LMdata)
}


head.LM.data.frame <- function(LMdata){
  print(head(LMdata$LMdata))
}

#########################################################################
#### FIT LM-ING MODEL  ####
#########################################################################

# -----------------------------------------------------------------------
# fitLM : fit a coxph or CSC model to a LM super dataset 
# -----------------------------------------------------------------------
# Input:
# - type             : "coxph" or "CSC"/"CauseSpecificCox"
# - formula          : formula to be used. 
# - LMdata           : An object of class "LM.data.frame", this can be created by running cutLMsuper and addLMtime
# -----------------------------------------------------------------------
# Output:
# An object of class "LMcoxph" or "LMCSC" with components
# - superfm: fitted model
# - type: as input
# - w, func_covars, func_LMs, LMcovars, allLMcovars, outcome: as in LMdata
# - LHS: the LHS of the input formula, stored for later use in calibration/discrimination
# -----------------------------------------------------------------------
# TODO: add FGR?
# TODO: handle formulas that use ~ ., etc (use all with LM-interaction)
# TODO: store LMdata somehow for NULL prediction later
fitLM <- function(type, formula, LMdata, ...){
  LHS = Reduce(paste, deparse(formula[[2]]))

  if(class(LMdata)!="LM.data.frame"){
    stop("data must be of type LM.data.frame with attributes w, func_covars, func_LMs, etc") #TODO: improve
  }
  data=LMdata$LMdata
  
  if(type=="coxph"){
    superfm <- coxph(formula, data, ...)
    cl <- "LMcoxph"
    
  } else if (type=="CauseSpecificCox" | type=="CSC"){
    superfm <- CSC(formula, data, ...)
    cl <- "LMCSC"
  }
  
  out=list(superfm=superfm, type=type, 
           w = LMdata$w,
           func_covars=LMdata$func_covars, func_LMs=LMdata$func_LMs, 
           LMcovars=LMdata$LMcovars, allLMcovars=LMdata$allLMcovars,
           outcome=LMdata$outcome, LHS=LHS)
  class(out)=cl
  
  return(out)
}


#########################################################################
#### VISUALING HAZARD RATIOS  ####
#########################################################################

# ----------------------------------------------------------
# find_se: Helper function to calculate SE for time varying log HR
# of the form coef[1] + t*coef[2] + t^2*coef[2] + ..
# ----------------------------------------------------------
find_se <- function(t, coefs, covar, func_covars){
  form <- "x1"
  if (length(coefs) > 1){
    for (i in 2:length(coefs)){
      form <- paste(form, sprintf("%s * %f", paste("x",i,sep=""), func_covars[[i]](t)),sep=" + ")
    }
  }
  form <- paste("~", form) 
  se <- deltamethod(as.formula(form), coefs, covar)
  return(se)
}


# -----------------------------------------------------------------------
# plot_dynamic_HR: Plots the dynamic hazard ratio
# -----------------------------------------------------------------------
# Input:
# - superfm          : An object of class "LMcoxph" or "LMCSC", i.e. a fitted supermodel
# - covars           : variables to plot the HR of  (note these must be given without time interaction label, for e.g., as in LMcovars)
# - end_time         : final time point to plot HR
# - CI               : include confidence intervals or not
# - cause            : cause of interest if class(fm) is CSC
# -----------------------------------------------------------------------
# TODO: for non-binary variables allow for choice of value (instead of assumed=1)
# TODO: change to base R plots
plot_dynamic_HR <- function(superfm, covars=NULL, end_time=3, CI=FALSE, cause=NULL){
  fm = superfm$superfm
  
  if (is.null(covars)){ covars <- superfm$LMcovars }
  
  if (superfm$type == "coxph") {
    if (!is.null(cause)) {stop("no cause should be input for coxph supermodels.")}
    num_causes <- 1
    bet <- fm$coefficients
    func_covars <- superfm$func_covars[[1]]
    if(CI){ sig <- vcov(fm) }
    
  } else if (superfm$type == "CauseSpecificCox" | superfm$type == "CSC") {
    if (is.null(cause)) { cause <- as.numeric(fm$theCause) }
    
    num_causes <- length(fm$models) 
    bet <- fm$models[[cause]]$coefficients
    
    
    models <- fm$models
    num_causes <- length(models)
    nb_func_covars <- length(superfm$func_covars)
    nb_func_LM <- length(superfm$func_LM)

    if (nb_func_covars == 1){
      superfm$func_covars <- rep(superfm$func_covars, num_causes)
    } else if (nb_func_covars != num_causes){
      stop("Error in func_covar argument. Wrong number of sets of functions.
           \nThere should as many sets of functions as causes or one set that is used for all causes.")
    }
    if (nb_func_LM == 1){
      superfm$func_LM <- rep(superfm$func_LM, num_causes)
    } else if (nb_func_LM != num_causes){
      stop("Error in func_covar argument. Wrong number of sets of functions.
           \nThere should as many sets of functions as causes or one set that is used for all causes.")
    }
    func_covars <- superfm$func_covars[[cause]]
    
    if(CI){ sig <- vcov(fm$models[[cause]]) }
  }
  

  t <- seq(0, end_time, by=0.1)
  
  plots <- lapply(1:length(covars), function(i){
    idx <- startsWith(names(bet), covars[i])
    bet_var <- bet[idx]

    HR <- sapply(t, function(x){ # eval HR over times x in t
      sum(sapply(1:length(bet_var), function(i){
        bet_var[i] * func_covars[[i]](x) 
      })) # bet0 + bet1*x + bet2*x^2 + ... 
    }) 
    
    # plot HR vs baseline (=1)
    g <- ggplot(data.frame(t=t, HR=HR), aes(t, HR))+
        geom_line(color="darkblue")+
        geom_line(aes(t, rep(0,end_time/0.1+1)), color="dodgerblue", alpha=0.4) +
        labs(x="Prediction time", y="log HR",title=covars[i])
    
    if(CI){
      se <- sapply(t, find_se, bet_var, sig[idx,idx], func_covars)
      lower <- HR - 1.96*se
      upper <- HR + 1.96*se
      g <- g +
        geom_line(aes(t,lower), color="darkblue", lty=2) +
        geom_line(aes(t,upper), color="darkblue", lty=2)
    }
    return(g)
  })

  return(plot_grid(plotlist=plots))
}


#########################################################################
#### MAKING PREDICTIONS USING THE SUPERMODEL ####
#########################################################################

# -----------------------------------------------------------------------------
# riskScore: Helper function for predLMsurv & predLMrisk
# Calcutes dynamic risk score at a time for an individual
# -----------------------------------------------------------------------------
# INPUT:
# - fm      :  fitted model (coxph)
# - tt      :  time point at which to calculate risk score
# - data    :  dataframe (single row) of individual. Must contain the original covariates.
# - func_covars, func_LM : Time-interaction functions for the specific model 
# -----------------------------------------------------------------------------
riskScore <- function(fm, tt, data, func_covars, func_LM, cause=NULL)
{
  coefs <- fm$coefficients
  pred_covars <- names(coefs)
  idx_LM_covars <- grep("LM",pred_covars, fixed=TRUE)
  LM_covars <- pred_covars[idx_LM_covars]
  bet_covars <- pred_covars[-idx_LM_covars]

  # coef_LM1*g1(t) + coef_LM2*g2(t) + ...
  risk <- sum(
    sapply(LM_covars, function(coef_name){
      # Get associated function
      n <- nchar(coef_name)
      idx <- as.numeric(substr(coef_name,n,n))
      g <- func_LM[[idx]]
      # Multiply by coef
      return(g(tt) * coefs[coef_name])
    })
    # X1*coef + X1*t*coef + X1*t^2*coef + ..
  ) + sum(
    sapply(bet_covars, function(coef_name){
      # Get associated function
      n <- nchar(coef_name)
      idx <- as.numeric(substr(coef_name,n,n))
      f <- func_covars[[idx]]
      # Get associated covariate info (remove _i from the name)
      covar <- substr(coef_name,1,n-2)
      # Multiply both by coef
      return(f(tt) * coefs[coef_name] * data[,covar])
    })
  )
  return(risk)
}


# -----------------------------------------------------------------------------
# predLMsurv: Custom-made function to calculate w-year survival from a LM point
# -----------------------------------------------------------------------------
# INPUT:
# - fm      :  fitted LM supermodel 
# - newdata :  dataframe of individuals to make predictions for. 
#              Must contain the original covariates.
# - tt      :  time points at which to predict risk of w more years. Note tt must be one value for the whole dataframe newdata or must have the same length as the number of rows of newdata (each datapoint is associated with one LM/prediction time point).
# -----------------------------------------------------------------------------
# OUTPUT: An object of class "LMpred" with components:
# - preds: a dataframe with columns time and risk, each entry corresponds to one individual and prediction time point
# - w, type, LHS: as in the fitted super model
# - data: the newdata, stored for calibration/discrimination 
# -----------------------------------------------------------------------------
# TODO: handle null newdata
# TODO: run and test again
# TODO: allow for different prediction window than w (with a warning)
predLMsurv <- function(superfm, newdata, tt)
{
  func_covars <- superfm$func_covars
  func_LM <- superfm$func_LM
  w <- superfm$w
  fm <- superfm$superfm
  type <- superfm$type
  
  num_preds <- nrow(newdata)
  if(!(length(tt)==num_preds)){
    if (length(tt) == 1){
      tt<-rep(tt,num_preds)
    } else {
      stop("Error in newdata or tt. Must have length(tt) == nrow(newdata) or tt be one landmarking point.")
    }
  }
  if (type == "coxph") {
    models <- list(fm)
    num_causes <- 1
    if (length(func_covars) != 1){
      stop("Error in func_covars argument. There should be only one set of functions. 
           \nDid you forget list() around your sets of functions?")
    }
    if (length(func_LM) != 1){
      stop("Error in func_covars argument. There should be only one set of functions. 
           \nDid you forget list() around your sets of functions?")
    }
  } else if (type == "CauseSpecificCox" | type == "CSC") {
    models <- fm$models
    num_causes <- length(models) 
    nb_func_covars <- length(func_covars)
    nb_func_LM <- length(func_LM)
    
    if (nb_func_covars == 1){
      func_covars <- rep(func_covars, num_causes)
    } else if (nb_func_covars != num_causes){
      stop("Error in func_covar argument. Wrong number of sets of functions. 
           \nThere should as many sets of functions as causes or one set that is used for all causes.")
    }
    if (nb_func_LM == 1){
      func_LM <- rep(func_LM, num_causes)
    } else if (nb_func_LM != num_causes){
      stop("Error in func_covar argument. Wrong number of sets of functions. 
           \nThere should as many sets of functions as causes or one set that is used for all causes.")
    }
  } else {
    stop("Error in fm argument. Input model is of the wrong type. 
         \nOnly supported classes are coxph and CauseSpecificCox")
  }
  
  bet <- lapply(models,function(model) model$coefficient)
  num_covarsLM <- sapply(bet, length) 
  
  # Create a baseline individual for each cause
  base_data <- lapply(1:num_causes, function(i){
    data.frame(matrix(rep(0,num_covarsLM[i]), ncol = num_covarsLM[i], dimnames=list(c(""),names(bet[[i]]))))
  })
  # Baseline hazards
  sf <- lapply(1:num_causes,function(i) survfit(models[[i]], newdata=base_data[[i]]))
  sf <- lapply(sf, function(s) data.frame(time=s$time,surv=s$surv,Haz=-log(s$surv)))
  Fw <- rep(NA,length(tt))
  
  for (i in 1:num_preds) {
    tti <- tt[i]
    newdatai <- newdata[i,]
    risks <-  sapply(1:num_causes, function(c) riskScore(models[[c]], tti, newdatai, func_covars[[c]], func_LM[[c]]))
    s<-0
    for (c in 1:num_causes){
      sfc <- sf[[c]]
      sfc$Haz <- sfc$Haz * exp(risks[c])
      f<- stepfun(sfc$time, c(0,sfc$Haz), right=FALSE)
      s<-s+(f(tti+wj)-f(tti))
    }
    Fw[i]<-s
  }
  preds <- data.frame(time=tt,Fw=exp(-Fw))
  
  out = list(preds=preds, w=w, type=type, data=newdata, LHS=superfm$LHS)
  return(out)}


# -----------------------------------------------------------------------------
# predLMrisk: Custom-made function to calculate (cause-specific) w-year risk from a LM point
# -----------------------------------------------------------------------------
# INPUT:
# - fm      :  fitted LM supermodel 
# - newdata :  dataframe of individuals to make predictions for. 
#              Must contain the original covariates.
# - tt      :  time points at which to predict risk of w more years. Note tt must be one value for the whole dataframe newdata or must have the same length as the number of rows of newdata (each datapoint is associated with one LM/prediction time point).
# - cause   :  cause of interest
# -----------------------------------------------------------------------------
# OUTPUT: An object of class "LMpred" with components:
# - preds: a dataframe with columns time and risk, each entry corresponds to one individual and prediction time point
# - w, type, LHS: as in the fitted super model
# - data: the newdata, stored for calibration/discrimination 
# -----------------------------------------------------------------------------
# TODO: handle null newdata
# TODO: allow for different prediction window than w (with a warning)
predLMrisk <- function(superfm, newdata, tt, cause=NULL)
{
  # TODO: check superfm/args are correct 
  
  func_covars <- superfm$func_covars
  func_LM <- superfm$func_LM
  w <- superfm$w
  fm <- superfm$superfm
  type <- superfm$type
  
  num_preds <- nrow(newdata)
  if(!(length(tt)==num_preds)){
    if (length(tt) == 1){
      tt<-rep(tt,num_preds)
    } else {
      stop("Error in newdata or tt. Must have length(tt) == nrow(newdata) or tt be one landmarking point.")
    }
  }
  if (type == "coxph") {
    models <- list(fm)
    num_causes <- 1
    if (length(func_covars) != 1){
      stop("Error in func_covars argument. There should be only one set of functions. 
           \nDid you forget list() around your sets of functions?")
    }
    if (length(func_LM) != 1){
      stop("Error in func_covars argument. There should be only one set of functions. 
           \nDid you forget list() around your sets of functions?")
    }
    if (!is.null(cause)){stop("No cause should be specified for a coxph model.")}
  } else if (type == "CauseSpecificCox" | type =="CSC") {
    models <- fm$models
    num_causes <- length(models) 
    nb_func_covars <- length(func_covars)
    nb_func_LM <- length(func_LM)
    
    if (is.null(cause)) cause <- 1
    if (!(cause %in% fm$causes)) stop("Error in cause. Not one of the causes in the model.")
    if (nb_func_covars == 1){
      func_covars <- rep(func_covars, num_causes)
    } else if (nb_func_covars != num_causes){
      stop("Error in func_covar argument. Wrong number of sets of functions. 
           \nThere should as many sets of functions as causes or one set that is used for all causes.")
    }
    if (nb_func_LM == 1){
      func_LM <- rep(func_LM, num_causes)
    } else if (nb_func_LM != num_causes){
      stop("Error in func_covar argument. Wrong number of sets of functions. 
           \nThere should as many sets of functions as causes or one set that is used for all causes.")
    }
  } else {
    stop("Error in fm argument. Input model is of the wrong type. 
         \nOnly supported classes are CauseSpecificCox and coxph")
  }
  
  # Create a baseline individual for this cause
  bet <- lapply(models,function(model) model$coefficient)
  num_covarsLM <- sapply(bet, length) 
  # Create a baseline individual for each cause
  base_data <- lapply(1:num_causes, function(i){
    data.frame(matrix(rep(0,num_covarsLM[i]), 
                      ncol = num_covarsLM[i], 
                      dimnames=list(c(""),names(bet[[i]]))))
  })
  
  # Baseline hazards
  sf <- lapply(1:num_causes,function(i) survfit(models[[i]], newdata=base_data[[i]]))
  sf <- lapply(sf, function(s) data.frame(time=s$time,surv=s$surv,Haz=-log(s$surv)))
  Fw <- rep(0, length(tt))
  
  sf1 <-sf[[cause]]
  Fw <- rep(NA, num_preds)
  for( i in 1:num_preds) {
    tti <- tt[i] # TODO: don't repeat for the same time!?
    newdatai <- newdata[i,]
    risk <- sapply(1:num_causes, function(c) riskScore(models[[c]], tti, newdatai, func_covars[[c]], func_LM[[c]]))
    
    pred_window <- (tti <= sf1$time & sf1$time <= tti+w)
    n_times <- sum(pred_window)
    
    haz <- sf1$Haz[pred_window] * exp(risk[cause])
    instHaz <- haz[2:n_times]-haz[1:n_times-1]
    idx <- (instHaz != 0)
    instHaz <- instHaz[idx]
    
    times <- sf1$time[pred_window][2:n_times]
    times <- times[idx]
    
    
    w_adj <- times-tti
    surv <- c()
    for (j in 1:length(w_adj)){
      s <- 0 
      wj<-w_adj[j]
      for (c in 1:num_causes){
        sfc <- sf[[c]]
        sfc$Haz <- sfc$Haz * exp(risk[c])
        f <- stepfun(sfc$time, c(0,sfc$Haz), right=FALSE)
        s <- s + (f(tti+wj)-f(tti))
      }
      surv <- c(surv,exp(-s))
    }
    Fw[i] <- sum(instHaz*surv)
  }
  
  preds = data.frame(time=tt,risk=Fw)
  
  out = list(preds=preds, w=w, type=type, LHS=superfm$LHS, data=newdata)
  class(out) = "LMpred"
  return(out)
}

#########################################################################
#### CALIBRATION ####
#########################################################################


# -----------------------------------------------------------------------------
# LMcalPlot: Calibration plots for dynamic risk prediction models
# -----------------------------------------------------------------------------
# INPUT:
# - preds: A named list of prediction models, where allowed entries are outputs from predLMsurv/predlMrisk
# - cause: Cause of interest
# - plot : If FALSE, do not plot the results, just return a plottable object.
# - sub  : If TRUE, add a subheading with the number of individuals at risk, and the number that under the event of interest
# - ...  : Additional arguments to pass to calPlot
# -----------------------------------------------------------------------------
# OUTPUT:
# - List of plots of w-year risk, one entry per LM/prediction time point
# -----------------------------------------------------------------------------
# TODO: Add option to show AUCt, Brier on plot
LMcalPlot <- function(preds,cause=1,..., plot=T,main=NULL,sub=T){
  # TODO: check data+formula+times+w for all preds is the same
  
  times = unique(preds[[1]]$preds$time)
  w = preds[[1]]$w
  LHS = preds[[1]]$LHS
  
  add_title=T
  if(!is.null(main)){add_title=F}
  
  outlist = list()
  for(t in 1:length(times)){
    tLM = times[t]
    
    idx = preds[[1]]$preds$time == tLM
    data_to_test = preds[[1]]$data[idx, ]
    risks_to_test = lapply(preds, function(p) p$preds$risk[idx])
    
    if (nrow(data_to_test)!=length(risks_to_test[[1]])){
      stop("nrow(data_to_test)!=length(risks_to_test)")
    }

    x = NULL
    x <- calPlot(
      risks_to_test, 
      time=tLM+w, 
      formula=as.formula(paste0(LHS,"~1")),
      data=data_to_test,
      cause=cause,
      plot=plot,
      ...
    )
    
    if (plot){
      if(add_title){ title(main = paste0("Calibration of ",w, "-yr risk \n measured at LM time ",tLM)) }
      else { title(main=main) }
      
      if(sub){ 
        num_patients = length( risks_to_test[[1]] )
        num_events = sum( data_to_test$event==cause)
        subtitle = paste0("#at risk=",num_patients,", #that undergo event=",num_events)
        title(sub = substitute(paste(italic(subtitle) )))
      }
    }
    
    
    outlist[[t]] <- x 	
  }
  return(outlist)
}

#########################################################################
#### DISCRIMINATION ####
#########################################################################

# -----------------------------------------------------------------------------
# print.LMScore: custom print function for objects of class "LMScore"
# -----------------------------------------------------------------------------
print.LMScore <- function(x,digits=3){
  
  if(nrow(x$auct)>0){
    cat(paste0("\nMetric: Time-dependent AUC for ",x$w,"-year risk prediction\n"))
    cat("\nResults by model:\n")
    
    AUC=se=lower=upper=delta.AUC=NULL
    fmt <- paste0("%1.",digits[[1]],"f")
    X <- copy(x$auct)
    X[,AUC:=sprintf(fmt=fmt,100*AUC)]
    if (match("se",colnames(X),nomatch=0)) X[,se:=NULL]
    if (match("lower",colnames(X),nomatch=0)) X[,lower:=sprintf(fmt=fmt,100*lower)]
    if (match("upper",colnames(X),nomatch=0)) X[,upper:=sprintf(fmt=fmt,100*upper)]
    
    print(X,digits=digits)
    
    cat("\nResults by comparison: \nTODO\n")
    
    message("\nNOTE: Values are multiplied by 100 and given in %.")
    message("NOTE: The higher AUC the better.")
    message(paste0("NOTE: Predictions are made at time tLM for ",x$w,"-year risk"))
  }
  if(nrow(x$briert)>0){
    cat(paste0("\nMetric: Brier Score for ",x$w,"-year risk prediction\n"))
    cat("\nResults by model:\n")
    
    Brier=se.conservative=se=lower=upper=delta.Brier=NULL
    fmt <- paste0("%1.",digits[[1]],"f")
    X <- copy(x$briert)
    X[,Brier:=sprintf(fmt=fmt,100*Brier)]
    if (match("se",colnames(X),nomatch=0)) X[,se:=NULL]
    if (match("se.conservative",colnames(X),nomatch=0)) X[,se.conservative:=NULL]
    if (match("lower",colnames(X),nomatch=0)) X[,lower:=sprintf(fmt=fmt,100*lower)]
    if (match("upper",colnames(X),nomatch=0)) X[,upper:=sprintf(fmt=fmt,100*upper)]
    print(X)
    
    cat("\nResults by comparison: \nTODO\n")
    
    message("\nNOTE: Values are multiplied by 100 and given in %.") 
    message("NOTE: The lower Brier the better.")
    message(paste0("NOTE: Predictions are made at time tLM for ",x$w,"-year risk"))
  }
  
}

# -----------------------------------------------------------------------------
# LMScore: Methods (AUCt, Brier) to score the predictive performance of dynamic risk markers from LM super models
# -----------------------------------------------------------------------------
# INPUT:
# - preds: A named list of prediction models, where allowed entries are outputs from predLMsurv/predlMrisk
# - cause: Cause of interest
# - metrics : Character vector specifying which metrics to apply. Choices are "auc" and "brier". Case matters.
# - ...  : Additional arguments to pass to Score
# -----------------------------------------------------------------------------
# OUTPUT: An object of class "LMScore", which has components:
# - dataframe auct if "auc" was a metric
# - dataframe briert if "brier" was a metric
# - w: from preds, prediction window of interest
# -----------------------------------------------------------------------------
LMScore <- function(preds,cause=1,metrics=c("auc","brier"), ...){
  # TODO: check data+formula+times+w for all preds is the same

  times = unique(preds[[1]]$preds$time)
  w = preds[[1]]$w
  LHS = preds[[1]]$LHS
  
  auct <- data.table()
  briert <- data.table()
  for(t in 1:length(times)){
    tLM = times[t]
    
    idx = preds[[1]]$preds$time == tLM
    data_to_test = preds[[1]]$data[idx, ]
    risks_to_test = lapply(preds, function(p) p$preds$risk[idx])
    
    if (nrow(data_to_test)!=length(risks_to_test[[1]])){
      stop("nrow(data_to_test)!=length(risks_to_test)")
    }
    
    score_t <- Score(risks_to_test, 
                     formula=as.formula(paste0(LHS,"~1")),
                     data=data_to_test, 
                     metrics=metrics, 
                     cause=cause, 
                     times=c(tLM+w),
                     ...)
    
    if ("auc" %in% metrics) auct <- rbind(auct, cbind(tLM,score_t$AUC$score)) 
    if ("brier" %in% metrics) briert <- rbind(briert, cbind(tLM,score_t$Brier$score)) 
    
  }
  outlist = list(auct=auct,briert=briert,w=w)
  class(outlist) = "LMScore"
  return(outlist)
}



#########################################################################
#### VISUALING RISK SCORE TRAJECTORIES  ####
#########################################################################

# -----------------------------------------------------------------------------
# plotRisk: Plots the absolute risk of individuals for different LM points for
#           an event of interest within a window 
# -----------------------------------------------------------------------------
# INPUT:
# - fm      :  Fitted LM super model 
# - data    :  Data frame of individuals from which to plot risk 
# - w       :  Prediction window
# - format  :  Character string specifying whether the data are in wide (default) or in long format
# - LM_col  :  Character string specifying the column name in data containing the (running) time variable 
#              associated with the time-varying covariate(s); only needed if format="long"
# - id_col  :  Character string specifying the column name in data containing the subject id; only needed if format="long"
# - cause   :  The cause we are looking at, only needed if class(fm) == "LMCSC"
# - varing  :  character string specifying column name in the data containing time-varying covariates; only needed if format="wide"
# -----------------------------------------------------------------------------
plotRisk <- function(fm, data, format, LM_col=NULL, id_col=NULL, cause=1, varying,
                     pch,lty,lwd,col,main,xlab,ylab,...){
  if(format=="long"){
    if(missing(id_col)) stop("argument 'id_col' should be specified for long format data")
    if(missing(LM_col)) stop("argument 'LM_col' should be specified for long format data")
    if(! id_col %in% colnames(data)) stop("arg 'id_col' is not a column in data")
    if(! LM_col %in% colnames(data)) stop("arg 'LM_col' is not a column in data")
    unique_ids <- unique(data[[id_col]])
    NF <- length(unique_ids)
  } else if (format=="wide"){
    if(missing(varing)) message("NOTE: wide format but no column with varying covariates is given.")
    NF <- nrow(data)
  } else {stop("format must be wide or long.")}
  
  ## Set up plot params
  if (missing(lwd)) 
    lwd <- rep(2, NF)
  if (missing(col)) {col <- 1:NF}
  if (missing(lty)) 
    lty <- rep(1, NF)
  if (missing(pch)) 
    pch <- rep(NA_integer_, NF)
  if (length(lwd) < NF) 
    lwd <- rep(lwd, NF)
  if (length(lty) < NF) 
    lty <- rep(lty, NF)
  if (length(col) < NF) 
    col <- rep(col, NF)
  if (length(pch) < NF) 
    pch <- rep(pch, NF)
  if(missing(main))
    main <- paste0(fm$w,"-year dynamic risk prediction")
  if(missing(xlab))
    xlab <- "LM prediction time"
  if(missing(ylab))
    ylab <- "Risk"

  ## Create plot
  if (format == "long") {
    for (i in 1:NF){
      id = unique_ids[i]
      data_ind <- data[data[[id_col]] == id,]
      x <- data_ind[[LM_col]]
      y <- predLMrisk(fm, data_ind, x, cause)$preds$risk
      if(i==1) plot(stepfun(x,c(y[1],y)), xlab=xlab, ylab=ylab,main=main,pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
      else lines(stepfun(x,c(y[1],y)),pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
    }
    legend(x = "topright",
           legend = unique_ids,
           lty = lty,
           col = col,
           lwd = lwd)
    
  } else if (format == "wide") { # TODO: test this with data
    if (is.null(end_time)) end_time = max(data[[outcome$time]])
    time_col=fm$outcome$time
    
    for (row in 1:NF){
      ind = data[row,] 
      if(ind[[covs$varying]] != ind[[time_col]]){
        x <- c(0, ind[[varying]])
        ind = cutLMsuper(ind, outcome, x, end_time, covs, format="wide") 
      } else {
        x <- c(0)
        ind = cutLM(ind, outcome, 0, end_time, covs, format="wide")
      }
      
      y <- sapply(1:nrow(ind), function(row) {
        predLMrisk(fm,  ind[row,], x[row])
      })
      if (x[1]!=end_time){
        x <- c(x, end_time)
        y <- c(y, y[length(y)])
      }
      if(i==1) plot(stepfun(x,c(y[1],y)), xlab=xlab, ylab=ylab,main=main,pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
      else lines(stepfun(x,c(y[1],y)),pch=pch[i],lty=lty[i],lwd=lwd[i],col=col[i],...)
    }
  
  } 
}

