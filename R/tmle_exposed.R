#'@title Targeted minimum-loss based estimator for stochastic and deterministic interventions among the exposed
#'@description The \code{tmle_exposed} is a Targeted Minumum-loss based
#'estimator (TMLE) that can be used to estimate effects of interventions
#'targeting an intermediate factor among the exposed. The intermediate factor
#'(Z) is an effect of the exposure (A) and a cause of the outcome (Y), A -> Z ->
#'Y. According to the intervention of interest, different parameters and
#'estimators can be chosen using \code{intervention}. If
#'\code{intervention='unexposed'}, which is the default, then the estimated
#'parameter is the expected change in outcome risk among the exposed (A=1) if
#'hypothetically the exposed had the same probability of the intermediate (Z)
#'variable as the unexposed (A=0). Other parameteres can be estimated by setting
#'\code{intervention} to a number between 0 and 1. If \code{intervention=1},
#'then the estimated parameter is the expected outcome risk among the exposed
#'(A=1) if all exposed received the intermediate variable (Z).
#'\code{intervention} can be fixed to any probability of interest, for example
#'if \code{intervention=0.6}, then the estimated parameter would be the expected
#'outcome risk among the exposed (A=1) if all exposed had 60% chance of the
#'intermediate variable (Z). Regardless of which \code{intervention} is chosen,
#'the expected outcome risk under intervention is compared to the outcome risk
#'among the exposed when the distribution of the intermediate variable was set
#'to the observed level among the exposed. Importantly, for this estimator the
#'exposure, intermediate variable, and outcome must all be binary. The
#'underlying model for the exposure, intermediate, and outcome, which are needed
#'to estimate any of the parameters, can be modelled using Super learning. Super
#'learning can be used to produce a weighted combination of candidate algorithms
#'that optimize the cross-validated loss-function or to select the single best
#'perfoming algorithm among the candidate algorithms, also known as the discrete
#'Super learner. One must define a library of candidate algortihms which should
#'be considered by the Super learner. If the Super learner library contains only
#'one algorithm, results will be estimated based on this algorithm alone, and
#'thus, not using Super Learning.
#'
#'@name tmle_exposed
#'
#'@author Amalie Lykkemark Moller \email{amalielykkemark@@live.dk}
#'
#'@usage tmle_exposed(data=data, intervention='unexposed', discrete.SL=TRUE,
#'  exposure.A=NA, intermediate.Z=NA, outcome.Y=NA, cov.A, cov.Z, cov.Y,
#'  SL.lib.A=FALSE, SL.lib.Z=FALSE, SL.lib.Y=FALSE, iterations=10)
#'
#'@param data  A data frame/data table with a binary exposure, a binary
#'  intermediate variable, a binary outcome, and covariates.
#'@param intervention  Either the character string \code{unexposed}, in which
#'  case ...  or a value, X, which is between 0 and 1, indicating the average
#'  probability of having the intermediate = 1. See details.... bla If
#'  \code{unexposed} then the intervention will be set to the distribution of
#'  the intermediate among the unexposed (A=0).
#'@param discrete.SL  If \code{FALSE} TMLE will use the weighted combination of
#'  the algorithms in the super learner library that minimizes the loss
#'  function. Defaults to \code{TRUE}, in which case the discrete super learner
#'  i used, i.e. the best performing algorithm.
#'@param exposure.A  Name of the binary exposure.
#'@param intermediate.Z  Name of the binary intermediate variable, which is the
#'  target of the hypothetical intervention.
#'@param outcome.Y  Name of the binary outcome.
#'@param cov.A  A vector containing names of possible confounders which should
#'  be included in models of the exposure.
#'@param cov.Z  A vector of confounders which should be included in models of
#'  the intermediate. Do not include the exposure as the function does this.
#'@param cov.Y  A vector of confounders which should be included in models of
#'  the outcome. Do not include the exposure and the intermediate as the
#'  function does this.
#'@param SL.lib.A  A vector of algorithms that should be considered by the super
#'  learner when modelling the exposure. All algorithms must be specified as
#'  Super Learner objects.
#'@param SL.lib.Z  A vector of algorithms for modelling the intermediate. All
#'  algorithms must be specified as Super Learner objects.
#'@param SL.lib.Y  A vector of algorithms for modelling the outcome. All
#'  algorithms must be specified as Super Learner objects.
#'@param iterations  Number of iterations for the updating step in TMLE.
#'  Defaults to 10.
#'
#'@details The structure of the data should be as follows: \item For the binary
#'exposure (\code{exposure.A) 1 = exposed and 0 = unexposed. \item For the
#'binary intermediate variable (\code{intermediate.Z}) 1 = treatment and 0 = no
#'treatment. \item For the binary outcome (\code{outcome.Y}) 1 = event and 0 =
#'no event.
#'
#'@return The function outputs the absolute outcome risk among the exposed under
#'  a hypothetical intervention ('unexposed' or fixed chance of intermediate)
#'  (psi0 or psi0star), the absolute outcome risk among the exposed under no
#'  intervention, where the probability of the intermediate variable is as
#'  observed (psi1), and the absolute risk difference between the two (psi or
#'  psistar) along with the standard errors for psi0/psi0star, psi1, and
#'  psi/psistar.
#'
#' @examples
#'library(data.table)
#'library(SuperLearner)
#'require(tmleExposed)
#'n=5000
#'set.seed(1)
#'sex <- rbinom(n,1,0.4)
#'age <- rnorm(n,65,sd=5)
#'disease <- rbinom(n,1,0.6)
#'
#'A <- rbinom(n, 1, plogis(-3+0.05*age+1*sex))
#'Z <- rbinom(n, 1, plogis(5-0.08*age+1*sex-1.2*disease-0.8*A+0.01*A*disease))
#'Y <- rbinom(n, 1, plogis(-9+0.09*age+0.5*sex+0.8*disease-1.2*Z+0.7*A))
#'
#'d <- data.table(id=1:n, exposure=as.integer(A), intermediate=as.integer(Z), outcome=as.integer(Y), age, sex, disease)
#'
#'##### Define algorithms for the Super Learner library #####
#'lib = c('SL.glm','SL.step.interaction')
#'
#' #intervention: changing probability of the intermediate (Z=1) among the exposed (A=1)
#' #to what it would have been had they been unexposed (A=0).
#' #target parameter: the change in outcome among the exposed (A=1) had their chance of
#' #the intermediate (Z=1) been as among the unexposed (A=0).
#'
#'res<-tmle_exposed(data=d,
#'                  intervention = 'unexposed',
#'                  exposure.A='exposure',
#'                  intermediate.Z='intermediate',
#'                  outcome.Y='outcome',
#'                  cov.A=c('sex','age'),
#'                  cov.Z =c('sex','age','disease'),
#'                  cov.Y=c('sex','age','disease'),
#'                  SL.lib.A = lib,
#'                  SL.lib.Z = lib,
#'                  SL.lib.Y = lib,
#'                  discrete.SL = FALSE)
#'
#'summary(res)
#'
#'#intervention: all exposed (A=1) receives the intermediate (Z=1).
#'#target parameter: the change in outcome among the exposed had all
#'#exposed received the intermediate (Z=1).
#'tmle_exposed(data=d,
#'             intervention = 1,
#'             exposure.A='exposure',
#'             intermediate.Z='intermediate',
#'             outcome.Y='outcome',
#'             cov.A=c('sex','age'),
#'             cov.Z =c('sex','age','disease'),
#'             cov.Y=c('sex','age','disease'),
#'             SL.lib.A = lib,
#'             SL.lib.Z = lib,
#'             SL.lib.Y = lib,
#'             discrete.SL = FALSE)
#'
#'intervention: all exposed (A=1) had 60% chance of reciving
#'#the intermediate (Z=1).
#'target parameter: the change in outcome among the exposed had
#'#60% exposed received the intermediate (Z=1).
#'tmle_exposed(data=d,
#'             intervention = 0.6,
#'             exposure.A='exposure',
#'             intermediate.Z='intermediate',
#'             outcome.Y='outcome',
#'             cov.A=c('sex','age'),
#'             cov.Z =c('sex','age','disease'),
#'             cov.Y=c('sex','age','disease'),
#'             SL.lib.A = lib,
#'             SL.lib.Z = lib,
#'             SL.lib.Y = lib,
#'             discrete.SL = FALSE)
#'
#'@export
tmle_exposed<-function(data=data,
                       intervention='unexposed',
                       discrete.SL=TRUE,
                       exposure.A=NA,
                       intermediate.Z=NA,
                       outcome.Y=NA,
                       cov.A,
                       cov.Z,
                       cov.Y,
                       SL.lib.A=FALSE,
                       SL.lib.Z=FALSE,
                       SL.lib.Y=FALSE,
                       iterations=10){

  require(data.table)
  require(SuperLearner)
  require(riskRegression)

  #variable and input check
  if(is.na(exposure.A)|is.na(intermediate.Z)|is.na(outcome.Y)) {
    stop(paste('Please speficy names for the exposure, intermediate, and outcome variables'))
  }
  dt<-copy(data)
  data.table::setDT(dt)

  if(exists('id',dt)==FALSE){
    dt[,id:=1:.N]
  }

  #convert variable names
  if (exposure.A!='A'){
    setnames(dt,exposure.A,'A')
  }

  if (intermediate.Z!='Z'){
    setnames(dt,intermediate.Z,'Z')
  }

  if (outcome.Y!='Y'){
    setnames(dt,outcome.Y,'Y')
  }

  # #stop if bigger than 2 instead of this
  # if (length(unique(dt[,A]))!=2 | length(unique(dt[,Z]))!=2 | length(unique(dt[,Y]))!=2) {
  #   stop('Exposure, interediate, and outcome must be binary.')
  # }

  if (!is.integer(dt[,A])){
    dt[,A:=as.integer(A)]
    warning(paste('The exposure variable has been converted to integer'))
  }

  if (!is.integer(dt[,Z])) {
    dt[,Z:=as.integer(Z)]
    warning(paste('The intermediate variable has been converted to integer'))
  }

  if (!is.integer(dt[,Y])) {
    dt[,Y:=as.integer(Y)]
    warning(paste('The outcome variable has been converted to integer'))
  }


  pibar <- dt[,mean(A==1)]
  pifit <- SuperLearner::SuperLearner(Y=dt[,as.numeric(A==1)],
                                      X=dt[,.SD,.SDcols=cov.A],
                                      family = binomial(),
                                      SL.library = SL.lib.A)

  gammafit <-SuperLearner::SuperLearner(Y=dt[,as.numeric(Z==1)],
                                        X=dt[,.SD,.SDcols=c(cov.Z,'A')],
                                        family = binomial(),
                                        SL.library = SL.lib.Z)

  Qfit <-SuperLearner::SuperLearner(Y=dt[,as.numeric(Y==1)],
                                    X=dt[,.SD,.SDcols=c(cov.Y,'A','Z')],
                                    family = binomial(),
                                    SL.library = SL.lib.Y)

  if (discrete.SL==FALSE){
    dt[, pihat:=predict(pifit, newdata=dt[,.SD,.SDcols=cov.A], onlySL = T)$pred]

    if (intervention!='unexposed'){
      dt[, gamma.star:=as.numeric(intervention)]
    }

    if (intervention=='unexposed'){
      dt[, gammahat.a0:=predict(gammafit, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[, A:=0], onlySL = T)$pred]
    }

    dt[, gammahat.a1:=predict(gammafit, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[,A:=1], onlySL = T)$pred]

    dt.full <- data.table(rbind(copy(dt)[, Z:=1],
                                copy(dt)[, Z:=0]),
                          Z.obs=c(dt[, Z], dt[, Z]))

    dt.full[, Qhat:=predict(Qfit, newdata=dt.full[,.SD,.SDcols=c(cov.Y,'A','Z')], onlySL = T)$pred]
    dt.full[, Qhat.a1:=predict(Qfit, newdata=copy(dt.full[,.SD,.SDcols=c(cov.Y,'Z')])[,A:=1], onlySL = T)$pred]
    dt.full[, Qhat.a1.z0:=predict(Qfit, newdata=copy(dt.full[,.SD,.SDcols=c(cov.Y)])[, `:=`(A=1, Z=0)], onlySL = T)$pred]
    dt.full[, Qhat.a1.z1:=predict(Qfit, newdata=copy(dt.full[,.SD,.SDcols=c(cov.Y)])[, `:=`(A=1, Z=1)], onlySL = T)$pred]
  }

  if (discrete.SL==TRUE){
    p<-pifit$libraryNames[which.max(pifit$coef)]
    pifit.discrete<-pifit$fitLibrary[[p]]
    g<-gammafit$libraryNames[which.max(gammafit$coef)]
    gammafit.discrete<-gammafit$fitLibrary[[g]]
    Q<-Qfit$libraryNames[which.max(Qfit$coef)]
    Qfit.discrete<-Qfit$fitLibrary[[Q]]
    dt[, pihat:=predict(pifit.discrete, newdata=copy(dt[,.SD,.SDcols=c(cov.A)]), type="response")]

    if (intervention!='unexposed'){
      dt[, gamma.star:=as.numeric(intervention)]
    }

    if (intervention=='unexposed'){
      dt[, gammahat.a0:=predict(gammafit.discrete, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[, A:=0], type="response")]
    }
    dt[, gammahat.a1:=predict(gammafit.discrete, newdata=copy(dt[,.SD,.SDcols=c(cov.Z)])[, A:=1], type="response")]

    dt.full <- data.table(rbind(copy(dt)[, Z:=1],
                                copy(dt)[, Z:=0]),
                          Z.obs=c(dt[, Z], dt[, Z]))

    dt.full[, Qhat:=predict(Qfit.discrete, newdata=dt.full, type="response")]
    dt.full[, Qhat.a1:=predict(Qfit.discrete, newdata=copy(dt.full)[, A:=1], type="response")]
    dt.full[, Qhat.a1.z0:=predict(Qfit.discrete, newdata=copy(dt.full)[, `:=`(A=1, Z=0)], type="response")]
    dt.full[, Qhat.a1.z1:=predict(Qfit.discrete, newdata=copy(dt.full)[, `:=`(A=1, Z=1)], type="response")]
  }

  # no intervention
  dt.full[, psi.1:=sum(Qhat.a1*(gammahat.a1*Z+(1-gammahat.a1)*(1-Z))), by="id"]
  # as among the unexposed
  if (intervention=='unexposed'){
    dt.full[, psi.0:=sum(Qhat.a1*(gammahat.a0*Z+(1-gammahat.a0)*(1-Z))), by="id"]
  }

  # intervention
  if(intervention!='unexposed'){
    dt.full[, psi.0.star:=sum(Qhat.a1*(gamma.star*Z+(1-gamma.star)*(1-Z))), by="id"]
  }

  psihat.1 <- dt.full[Z==Z.obs, mean((A==1)/pibar*(Y - Qhat.a1) +
                                       (A==1)/pibar*(Qhat.a1 - psi.1) +
                                       (A==1)/pibar*(psi.1))]

  dt.full[Z==Z.obs,eic1:=((A==1)/pibar*(Y - Qhat.a1) +
                            (A==1)/pibar*(Qhat.a1 - psi.1) +
                            (A==1)/pibar*(psi.1 - psihat.1))]

  # intervention (as among the exposed)
  # intervention as for chest pain

  if (intervention=='unexposed'){
    psihat.0 <-dt.full[Z==Z.obs, mean((A==1)/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                                      (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                                        (Y - Qhat.a1) +
                                        (A==0)/(1-pihat)*(pihat)/pibar*(Qhat.a1 - psi.0) +
                                        (A==1)/pibar*(psi.0))]

    dt.full[Z==Z.obs,eic0:=((A==1)/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                            (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                              (Y - Qhat.a1) +
                              (A==0)/(1-pihat)*(pihat)/pibar*(Qhat.a1 - psi.0) +
                              (A==1)/pibar*(psi.0 - psihat.0))]

    se0 <- sqrt(mean(dt.full[Z==Z.obs,eic0]^2)/nrow(dt))
    se1 <- sqrt(mean(dt.full[Z==Z.obs,eic1]^2)/nrow(dt))
    dt.full[Z==Z.obs,eic:=eic0-eic1]
    se.diff <-  sqrt(mean((dt.full[Z==Z.obs,eic])^2)/nrow(dt))
  }

  # intervention
  if (intervention!='unexposed'){
    psihat.0.star <- dt.full[Z==Z.obs,mean((A==1)/pibar*((Z*gamma.star+(1-Z)*(1-gamma.star))/
                                                           (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                                             (Y - Qhat.a1) +
                                             (A==1)/pibar*(psi.0.star))]

    dt.full[Z==Z.obs,eic0.star:= ((A==1)/pibar*((Z*gamma.star+(1-Z)*(1-gamma.star))/
                                                  (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                                    (Y - Qhat.a1) +
                                    (A==1)/pibar*(psi.0.star - psihat.0.star))]

    se0.star <- sqrt(mean(dt.full[Z==Z.obs,eic0.star]^2)/nrow(dt))
    se1 <- sqrt(mean(dt.full[Z==Z.obs,eic1]^2)/nrow(dt))
    dt.full[Z==Z.obs,eic:=eic0.star-eic1]
    se.diff <-  sqrt(mean((dt.full[Z==Z.obs,eic])^2)/nrow(dt))
  }

  dt.full.copy <- copy(dt.full)

  #-------- tmle
  #-- psi0;
  if (intervention=='unexposed'){

    for (iter in 1:iterations) {  #-- iterative updating;
      #-- update Q;
      dt.full[, H.Y:=1/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))]
      eps.Y <- coef(glm(Y ~ offset(qlogis(Qhat.a1))+H.Y-1, data=dt.full[Z==Z.obs & A==1], family=binomial()))

      dt.full[, Qhat.a1:=plogis(qlogis(Qhat.a1) + eps.Y/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                                                 (Z*gammahat.a1+(1-Z)*(1-gammahat.a1))))]
      dt.full[, Qhat.a1.z0:=plogis(qlogis(Qhat.a1.z0) + eps.Y/pibar*(((1-gammahat.a0))/
                                                                       ((1-gammahat.a1))))]
      dt.full[, Qhat.a1.z1:=plogis(qlogis(Qhat.a1.z1) + eps.Y/pibar*((gammahat.a0)/
                                                                       (gammahat.a1)))]
      #-- update Z;
      dt.full[, H.Z:=1/(1-pihat)*(pihat)/pibar*(Qhat.a1.z1 - Qhat.a1.z0)]
      eps.Z <- coef(glm(Z ~ offset(qlogis(gammahat.a0))+H.Z-1, data=dt.full[Z==Z.obs & A==0], family=binomial()))

      dt.full[, gammahat.a0:=plogis(qlogis(gammahat.a0) + eps.Z*H.Z)]
      dt.full[, psi.0:=sum(Qhat.a1*(gammahat.a0*Z+(1-gammahat.a0)*(1-Z))), by="id"]

      tmle.est0 <- dt.full[Z==Z.obs, (1/pibar)*mean((pihat)*psi.0)]

      #-- update A;
      dt.full[, H.A:=(psi.0-tmle.est0)/pibar]

      eps.A <- coef(glm(A==1 ~ offset(qlogis(pihat))+H.A-1, data=dt.full[Z==Z.obs], family=binomial()))
      dt.full[, pihat:=plogis(qlogis(pihat) + eps.A*H.A)]

      #updated estimate
      tmle.est0 <- dt.full[Z==Z.obs, (1/pibar)*mean((A)*psi.0)]


      solve.eic0 <- abs(dt.full[Z==Z.obs, mean((A==1)/pibar*((Z*gammahat.a0+(1-Z)*(1-gammahat.a0))/
                                                               (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                                                 (Y - Qhat.a1) +
                                                 (A==0)/(1-pihat)*(pihat)/pibar*(Qhat.a1.z1 - Qhat.a1.z0)*(Z-gammahat.a0) +
                                                 (psi.0 - tmle.est0)/pibar*((A==1) - pihat) +
                                                 pihat/pibar*(psi.0-tmle.est0))])

      if (solve.eic0<=se0/(log(nrow(dt))*sqrt(nrow(dt)))) break

      if(iter==iterations) {
        warning(paste('Efficient influence function for psi0 was not solved in',iterations,'iterations'))
      }
    }
  }
  #-- psi0star;
  if (intervention!='unexposed'){

    dt.full <- copy(dt.full.copy)

    for (iter in 1:iterations) {  #-- iterative updating;

      #-- update Q;
      dt.full[, H.Y:=1/pibar*((Z*gamma.star+(1-Z)*(1-gamma.star))/
                                (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))]
      eps.Y <- coef(glm(Y ~ offset(qlogis(Qhat.a1))+H.Y-1, data=dt.full[Z==Z.obs & A==1], family=binomial()))

      dt.full[, Qhat.a1:=plogis(qlogis(Qhat.a1) + eps.Y/pibar*((Z*gamma.star+(1-Z)*(1-gamma.star))/
                                                                 (Z*gammahat.a1+(1-Z)*(1-gammahat.a1))))]
      #-- update psi.0.star with the updated Q;
      dt.full[, psi.0.star:=sum(Qhat.a1*(gamma.star*Z+(1-gamma.star)*(1-Z))), by="id"]

      tmle.est0.star <- dt.full[Z==Z.obs, (1/pibar)*mean((pihat)*psi.0.star)]

      #-- update A;
      dt.full[, H.A:=(psi.0.star-tmle.est0.star)/pibar]
      eps.A <- coef(glm(A==1 ~ offset(qlogis(pihat))+H.A-1, data=dt.full[Z==Z.obs], family=binomial()))
      dt.full[, pihat:=plogis(qlogis(pihat) + eps.A*H.A)]

      #-- updated estimate
      tmle.est0.star <- dt.full[Z==Z.obs, (1/pibar)*mean((A)*psi.0.star)]

      solve.eic0.star <- abs(dt.full[Z==Z.obs, mean((A==1)/pibar*((Z*gamma.star+(1-Z)*(1-gamma.star))/
                                                                    (Z*gammahat.a1+(1-Z)*(1-gammahat.a1)))*
                                                      (Y - Qhat.a1) +
                                                      (A==1)/pibar*(psi.0.star - tmle.est0.star))])
      if (solve.eic0.star<=se0.star/(log(nrow(dt))*sqrt(nrow(dt)))) break
      if(iter==iterations) {
        warning(paste('Efficient influence function for psi0star was not solved in',iterations,'iterations'))
      }
    }
  }


  #-- psi1;
  dt.full <- copy(dt.full.copy)

  for (iter in 1:iterations) {  #-- iterative updating;

    #-- update Q;
    dt.full[, H.Y:=1/pibar]
    eps.Y <- coef(glm(Y ~ offset(qlogis(Qhat.a1))+H.Y-1, data=dt.full[Z==Z.obs & A==1], family=binomial()))

    dt.full[, Qhat.a1:=plogis(qlogis(Qhat.a1) + eps.Y/pibar)]
    dt.full[, Qhat.a1.z0:=plogis(qlogis(Qhat.a1.z0) + eps.Y/pibar)]
    dt.full[, Qhat.a1.z1:=plogis(qlogis(Qhat.a1.z1) + eps.Y/pibar)]

    dt.full[, psi.1:=sum(Qhat.a1*(gammahat.a1*Z+(1-gammahat.a1)*(1-Z))), by="id"]

    tmle.est1 <- dt.full[Z==Z.obs, (1/pibar)*mean((pihat)*psi.1)]

    #-- update A;
    dt.full[, H.A:=(psi.1-tmle.est1)/pibar]
    eps.A <- coef(glm(A==1 ~ offset(qlogis(pihat))+H.A-1, data=dt.full[Z==Z.obs], family=binomial()))
    dt.full[, pihat:=plogis(qlogis(pihat) + eps.A*H.A)]

    #-- updated estimate
    tmle.est1 <- dt.full[Z==Z.obs, (1/pibar)*mean((A)*psi.1)]

    solve.eic1 <- abs(dt.full[Z==Z.obs, mean((A==1)/pibar*(Y - Qhat.a1) +
                                               (psi.1 - tmle.est1)/pibar*((A==1) - (pihat)) +
                                               (pihat)/pibar*(psi.1-tmle.est1))])

    if (solve.eic1<=se1/(log(nrow(dt))*sqrt(nrow(dt)))) break
    if(iter==iterations) {
      warning(paste('Efficient influence function for psi1 was not solved in',iterations,'iterations'))
    }
  }

  # prepare output
  out<-list()
  if(intervention=='unexposed'){
    psi.diff.tmle<-tmle.est0-tmle.est1

    out$estimate=list(psi0=tmle.est0,psi1=tmle.est1,psi=psi.diff.tmle)
    out$se=list(se0=se0,se1=se1,se.diff=se.diff)

    #CV risk for exposure, intermediate, and outcome
    out$superlearner.CVrisk$A.exposure<-pifit$cvRisk
    out$superlearner.CVrisk$Z.intermediate<-gammafit$cvRisk
    out$superlearner.CVrisk$Y.outcome<-Qfit$cvRisk
    #weights assigned to each algorithm in the super learner library
    if (discrete.SL==TRUE){
      out$superlearner.discrete$A.exposure<-p
      out$superlearner.discrete$Z.intermediate<-g
      out$superlearner.discrete$Y.outcome<-Q
    }
    if (discrete.SL==FALSE){
      out$superlearner.weight$A.exposure<-pifit$coef
      out$superlearner.weight$Z.intermediate<-gammafit$coef
      out$superlearner.weight$Y.outcome<-Qfit$coef
    }
    out$distributions=rbind(distribution.A1=dt.full[Z==Z.obs & A==1, summary(pihat)],
                            distribution.Z.a1=dt.full[Z==Z.obs & A==1, summary(gammahat.a1)],
                            distribution.Z.a0=dt.full[Z==Z.obs & A==1, summary(gammahat.a0)],
                            distribution.Y=dt.full[Z==Z.obs & A==1, summary(Qhat)],
                            distribution.Y.a1=dt.full[Z==Z.obs & A==1, summary(Qhat.a1)])
  }

  if(intervention!='unexposed'){
    psi.diff.tmle<-tmle.est0.star-tmle.est1
    out$estimate=list(psi0star=tmle.est0.star,psi1=tmle.est1,psistar=psi.diff.tmle)
    out$se=list(se0.star=se0.star,se1=se1,se.diff.star=se.diff)

    #CV risk for exposure, intermediate, and outcome
    out$superlearner.CVrisk$A.exposure<-pifit$cvRisk
    out$superlearner.CVrisk$Z.intermediate<-gammafit$cvRisk
    out$superlearner.CVrisk$Y.outcome<-Qfit$cvRisk
    #weights assigned to each algorithm in the super learner library
    if (discrete.SL==TRUE){
      out$superlearner.discrete$A.exposure<-p
      out$superlearner.discrete$Z.intermediate<-g
      out$superlearner.discrete$Y.outcome<-Q
    }
    if (discrete.SL==FALSE){
      out$superlearner.weight$A.exposure<-pifit$coef
      out$superlearner.weight$Z.intermediate<-gammafit$coef
      out$superlearner.weight$Y.outcome<-Qfit$coef
    }
    out$distributions=rbind(distribution.A1=dt.full[Z==Z.obs & A==1, summary(pihat)],
                            distribution.Z.gamma.star=dt.full[Z==Z.obs & A==1, summary(gamma.star)],
                            distribution.Z.a1=dt.full[Z==Z.obs & A==1, summary(gammahat.a1)],
                            distribution.Y=dt.full[Z==Z.obs & A==1, summary(Qhat)],
                            distribution.Y.a1=dt.full[Z==Z.obs & A==1, summary(Qhat.a1)])

  }

  out$output.dataset<-dt.full

  class(out)<-'tmle_exposed'
  return(out)
}


