#########################################################################
#### BUILDING SUPER DATASET ####
#########################################################################

# -----------------------------------------------------------------------
# packyears2: helper function to calculate packyears, only used if fup
#             occurs between IPLC and prediction time point
# -----------------------------------------------------------------------

pack_years <- function(IPLC_py, smkstatus, smkstatus_fup, cigday, time_IPLC_to_fup, pred_time){
  time_fup_to_pred = pred_time - time_IPLC_to_fup
  
  ## if fup info is missing use BL info
  cigday_fup = cigday
  smkstatus_fup = if_else(is.na(smkstatus_fup), smkstatus, smkstatus_fup)
  
  ## smkstatus=2,smkstatus_fup=2 => py = IPLC_py
  ## smkstatus=2,smkstatus_fup=3 => py = IPLC_py + (time_fup to pred point)*cigday_fup/20 
  ## smkstatus=3,smkstatus_fup=2 => py = IPLC_py + (time IPLC to fup)*cigday/20
  ## smkstatus=3,smkstatus_fup=3 => py = IPLC_py + (time IPLC to fup)*cigday/20 + (time_fup to pred point)*cigday_fup/20 
  
  py = if_else(smkstatus==2,
                if_else(smkstatus_fup==2,
                        IPLC_py,
                        IPLC_py + time_fup_to_pred*cigday_fup/20),
                #smkstatus==3
                if_else(smkstatus_fup==2,
                        IPLC_py + time_IPLC_to_fup*cigday/20,
                        IPLC_py + time_IPLC_to_fup*cigday/20 + time_fup_to_pred*cigday_fup/20)
                )
  return(py)
}


# -----------------------------------------------------------------------
# update: updates covariates in LM data frame to be LM dependent
# ----------------------------------------------------------------------- 
update <- function(LMDAT){
  agevar="age_ix"
  pyvar ="packyears2"
  quitvar = "quityears2"
  statvar="smkstatus2"
  startage = 55
  stopage = 80
  py.thred = 30
  
  # TODO: lots of NA handling to think about
  
  if (class(LMDAT)=="LM.data.frame") { DAT <- LMDAT$LMdata }
  else { DAT <- LMDAT }

  DAT <- DAT %>%
    mutate(
      age_ix = age_ix + LM,
      
        # diff_fup_IPLC = time fup->IPLC
        # so -diff_fup_IPLC gives time IPLC(Time)->fup (>0 is relevant)
        # and -diff_fup_IPLC-LM gives time from pred time point -> fup
      fup_af_IPLC_bef_LM = if_else(is.na(diff_fup_IPLC),          # no fup
                                   F, if_else(-diff_fup_IPLC>0 & -diff_fup_IPLC-LM<0,  # fup after IPLC fup before LM
                                              T,F)),
       
      smkstatus_changed = if_else(fup_af_IPLC_bef_LM, 
                                  if_else(is.na(smkstatus_fup),F,
                                          if_else(smkstatus2!=smkstatus_fup, T, F)), F),
      smkstatus2 = if_else(fup_af_IPLC_bef_LM,
                           if_else(is.na(smkstatus_fup),
                                   smkstatus2,
                                   smkstatus_fup),
                           smkstatus2
      ),
      
      # stick to cigday2 - do not update
      # cigday2 =    if_else(fup_af_IPLC_bef_LM,
      #                      if_else(is.na(cigday_fup),
      #                              cigday2,
      #                              cigday_fup),
      #                      cigday2
      # ),
      
      ## smkstatus2 has been updated
      
      quityears2 = if_else(fup_af_IPLC_bef_LM,
                          
                          ## Need to update
                          if_else(is.na(smkstatus2),
                                  NA_real_, ## assume current smoker
                                  if_else(smkstatus2 == 2, ## former smoker
                                          ## now need to check if they become a former smoker or have always been
                                          if_else(smkstatus_changed,
                                                  diff_fup_IPLC+LM, ## quityears are from fup time to pred time
                                                  if_else(is.na(quityears2), LM, quityears2+LM)), ## otherwise prev amount + time to pred
                                          NA_real_)), ## current smoker
                          
                          ## only using info from baseline
                          if_else(is.na(smkstatus2),
                                  NA_real_, # current smoker
                                  if_else(smkstatus2 == 2, # former smoker
                                          if_else(is.na(quityears2),# no BL quit info
                                                  LM,               # assume at least LM
                                                  quityears2+LM),   # BL info+LM
                                          NA_real_)) # current smoker
      ),
      
      packyears2 = if_else(fup_af_IPLC_bef_LM,
                           pack_years(packyears2, smkstatus2, smkstatus_fup, cigday2, -diff_fup_IPLC, LM),
                           if_else(smkstatus2==2, # TODO: add an if_else for NA
                                   packyears2,
                                   packyears2 + LM * (cigday2/20))
      )
      
    ) %>% select(-fup_af_IPLC_bef_LM,-smkstatus_changed)
  
  DAT[["USPSTF2"]] = myEligibility2(agevar, pyvar, quitvar, statvar, DAT, startage, stopage, py.thred)
  DAT[["USPSTF.stage"]] = DAT$USPSTF2 * DAT$stage2.ix

  TESTDAT = FALSE
  if (TESTDAT) { 
    if (missing(add_vars)){
      stop("If wanting test data, add_vars that need to be removed must be specified.")
    }
    DAT <- DAT %>% select(-all_of(add_vars),-LM) 
  }
  
  if (class(LMDAT)=="LM.data.frame") { LMDAT$LMdata <- DAT }
  else { LMDAT <-  DAT }
  
  
  return(LMDAT)
}

#############################################################################
## CREATING NEW DATA 
## This is the same at cutLM but without censoring time at prediction window
#############################################################################

cutLMnewdata <- function (data, outcome, LM, horizon, covs, format = c("wide", 
                                                       "long"), id, rtime, right = TRUE) 
{
  format <- match.arg(format)
  if (format == "wide") {
    LMdat <- data
    if (!is.null(covs$varying)) 
      LMdat[[covs$varying]] <- 1 - as.numeric(LMdat[[covs$varying]] > 
                                                 LM)
  }
  else {
    if (missing(id)) 
      stop("argument 'id' should be specified for long format data")
    if (missing(rtime)) 
      stop("argument 'rtime' should be specified for long format data")
    ord <- order(data[[id]], data[[rtime]])
    data <- data[ord, ]
    ids <- unique(data[[id]])
    n <- length(ids)
    LMdat <- data[which(!duplicated(data[[id]])), ]
    for (i in 1:n) {
      wh <- which(data[[id]] == ids[i])
      di <- data[wh, ]
      idx <- cut(LM, c(data[[rtime]][wh], Inf), right = right, 
                 labels = FALSE)
      if (!is.na(idx)) 
        LMdat[i, ] <- di[idx, ]
      else {
        LMdat[i, ] <- di[1, ]
        LMdat[[covs$varying]][i] <- NA
        LMdat[[rtime]][i] <- NA
      }
    }
  }
  LMdat <- LMdat[LMdat[[outcome$time]] > LM, ]
  if (format == "long") 
    LMdat <- LMdat[!is.na(LMdat[[id]]), ]
  LMdat[outcome$status] <- LMdat[[outcome$status]] * as.numeric(LMdat[[outcome$time]] <= 
                                                                    horizon)
  # LMdat[outcome$time] <- pmin(as.vector(LMdat[[outcome$time]]), 
  #                              horizon)
  LMdat$LM <- LM
  if (format == "long") 
    cols <- match(c(id, outcome$time, outcome$status, covs$fixed, 
                    covs$varying, rtime, "LM"), names(LMdat))
  else cols <- match(c(outcome$time, outcome$status, covs$fixed, 
                       covs$varying, "LM"), names(LMdat))
  return(LMdat[, cols])
}