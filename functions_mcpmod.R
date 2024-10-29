# R functions


################################################################################################################




logit <- function(p) log(p / (1 - p))
inv_logit <- function(y) 1 / (1 + exp(-y))


################################################################################################################


# binary EP

analy_MCPMod2 <-
  function(sim.obj = NULL,
           dose = NULL,
           resp = NULL,
           data = NULL,
           Mods.obj.cand,
           S = NULL,
           type = "binomial",
           Delta,
           selModel = "aveAIC",
           direction = "decreasing",
           alpha = 0.025,
           ...) {
    if (is.null(sim.obj) == F &
        is.null(dose) == T & is.null(resp) == T & is.null(data) == T) {
      data <-  sim.obj
      dose <- data$dose
    } else if (is.null(sim.obj) == T &
               is.null(dose) == F & is.null(resp) == F & is.null(data) == T) {
      if (length(dose) != length(resp)) {
        stop(
          "Dose and response are required to be vectors of equal length specifying dose and response values, if 'data' is not specified."
        )
      }
    }
    
    dls <- unique(dose)
    
    if (type == "normal") {
      fmod <-
        DoseFinding::MCPMod(
          dose,
          resp,
          data,
          models = Mods.obj.cand,
          S = S,
          type = type,
          Delta = Delta,
          selModel = selModel,
          ...
        )
      poc <- ifelse(min(attr(fmod$MCTtest$tStat, "pVal")) < alpha, 1, 0)
      dseq <- seq(dls[1], dls[length(dls)], by = 0.005)
      if (length(fmod$selMod) > 0) {
        mod.pred <-
          do.call("cbind",
                  predict(
                    fmod,
                    doseSeq = dseq,
                    predType = "effect-curve",
                    se.fit = F
                  ))
        wtd.pred <- mod.pred %*% fmod$selMod
        med <- dseq[which(wtd.pred >= Delta)[1]]
      } else {
        med <- NA
      }
    } else if (type == "binomial") {
      fmod <-
        MCPModGeneral::MCPModGen(
          family = "binomial",
          link = "logit",
          returnS = TRUE,
          dose = "dose",
          resp = "resp",
          data = data,
          models = mods5,
          selModel = "aveAIC",
          doseType = "TD",
          pVal = TRUE,
          Delta = abs(Delta)
        )
      
      poc <-
        ifelse(min(attr(fmod$MCPMod$MCTtest$tStat, "pVal")) < alpha, 1, 0)
      
      dseq <- seq(dls[1], dls[length(dls)], by = 0.005)
      
      if (length(fmod$MCPMod$selMod) > 0) {
        mod.pred <-
          do.call(
            "cbind",
            predict(
              fmod$MCPMod,
              doseSeq = dseq,
              predType = "effect-curve",
              se.fit = F
            )
          )
        wtd.pred <- mod.pred %*% fmod$MCPMod$selMod
        
        if (direction == "increasing") {
          med <- dseq[which(wtd.pred >= Delta)[1]]
        }
        else if (direction == "decreasing") {
          med <- dseq[which(wtd.pred <= Delta)[1]]
        }
      } else {
        med <- NA
      }
    }
    
    res <- data.frame(poc, med)
    return(res)
  }



#######
##
##


doses <- c(0, 1, 2, 3)/3 # normalized equidistant doses
p.ocr <- 0.43
p.magli <- p.ocr * 0.5

# Candidate models
mods5 <- Mods(emax = c(0.25, 0.5), 
              sigEmax = rbind(c(0.25, 3), c(0.5, 4)), 
              betaMod = c(1.5, 0.1),
              placEff = logit(p.ocr), 
              maxEff = logit(p.magli) - logit(p.ocr), 
              doses = doses,
              direction = "decreasing")

mod.true <- Mods(emax = c(0.5), 
                 placEff = logit(p.ocr), 
                 maxEff = logit(p.magli) - logit(p.ocr), 
                 doses = doses,
                 direction = "decreasing")

