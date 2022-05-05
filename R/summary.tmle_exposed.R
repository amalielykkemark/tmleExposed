#'@rdname summary.tmle_exposed
#'@method summary tmle_exposed
#'@export
summary.tmle_exposed <- function(x,...) {
  if(identical(class(x), "tmle_exposed")){

    if(!is.null(x$estimate$psi)){
      cat("\nRisk difference among the exposed under stochastic intervention")
      cat("\n     Parameter estimate: ", signif(x$estimate$psi,5))
      cat("\n     Estimated variance: ", signif(x$se$se.diff,5))
      cat("\n     95% confidence interval:",paste("(", signif(x$estimate$psi-1.96*x$se$se.diff,5),", ", signif(x$estimate$psi+1.96*x$se$se.diff,5), ")", sep=""),"\n")
      cat('\n')
      cat("\nRisk among the exposed under stochastic intervention")
      cat("\n     Parameter estimate: ", signif(x$estimate$psi0,5))
      cat("\n     Estimated variance: ", signif(x$se$se0,5))
      cat("\n     95% confidence interval:",paste("(", signif(x$estimate$psi0-1.96*x$se$se0,5),", ", signif(x$estimate$psi0+1.96*x$se$se0,5), ")", sep=""),"\n")
      cat("\n     Distribution of chance of the intermediate under intervention among the exposed:\n")
      print(round(x$distributions['distribution.Z.a0',],4))
      cat('\n')

      cat("\nRisk among the exposed without intervention")
      cat("\n     Parameter estimate: ", signif(x$estimate$psi1,5))
      cat("\n     Estimated variance: ", signif(x$se$se1,5))
      cat("\n     95% confidence interval:",paste("(", signif(x$estimate$psi1-1.96*x$se$se1,5),", ", signif(x$estimate$psi1+1.96*x$se$se1,5), ")", sep=""),"\n")
      cat("\n     Distribution of chance of the intermediate without intervention among the exposed:\n")
      print(round(x$distributions['distribution.Z.a1',],4))
      cat('\n')
    }

    if(!is.null(x$estimate$psistar)){
      cat(paste0("\nRisk difference among the exposed for deterministic intervention (the probability of the intermediate was ",round(x$distributions['distribution.Z.gamma.star','Mean']*100),'%)'))
      cat("\n     Parameter estimate: ", signif(x$estimate$psistar,5))
      cat("\n     Estimated variance: ", signif(x$se$se.diff,5))
      cat("\n     95% confidence interval:",paste("(", signif(x$estimate$psistar-1.96*x$se$se.diff,5),", ", signif(x$estimate$psistar+1.96*x$se$se.diff,5), ")", sep=""),"\n")
      cat('\n')

      cat("\nRisk  among the exposed under deterministic intervention")
      cat("\n     Parameter estimate: ", signif(x$estimate$psi0star,5))
      cat("\n     Estimated variance: ", signif(x$se$se0.star,5))
      cat("\n     95% confidence interval:",paste("(", signif(x$estimate$psi0-1.96*x$se$se0.star,5),", ", signif(x$estimate$psi0+1.96*x$se$se0.star,5), ")", sep=""))
      cat("\n     Distribution of the intermediate under intervention among the exposed:\n")
      print(round(x$distributions['distribution.Z.gamma.star',],4))
      cat('\n')

      cat("\nRisk  among the exposed without intervention")
      cat("\n     Parameter estimate: ", signif(x$estimate$psi1,5))
      cat("\n     Estimated variance: ", signif(x$se$se1,5))
      cat("\n     95% confidence interval:",paste("(", signif(x$estimate$psi1-1.96*x$se$se1,5),", ", signif(x$estimate$psi1+1.96*x$se$se1,5), ")", sep=""))
      cat("\n     Distribution of chance of the intermediate without intervention among the exposed:\n")
      print(round(x$distributions['distribution.Z.a1',],4))
      cat('\n')
    }

    if (!is.null(x$superlearner.discrete)){
      cat("\nDiscrete super learner")
      cat("\n     Algorithm chosen for modelling the exposure: ", x$superlearner.discrete$A.exposure)
      cat("\n     Algorithm chosen for modelling the intermediate: ", x$superlearner.discrete$Z.intermediate)
      cat("\n     Algorithm chosen for modelling the outcome:", x$superlearner.discrete$Y.outcome,"\n")
    }

    if (!is.null(x$superlearner.weight)){
      cat("\nSuper learner weights")
      cat("\n     Weights for algorithms for the exposure model:\n")
      print(x$superlearner.weight$A.exposure)
      cat("\n     Weights for algorithms for the intermediate model:\n")
      print(x$superlearner.weight$Z.intermediate)
      cat("\n     Weights for algorithms for the outcome model:\n")
      print(x$superlearner.weight$Y.outcome)

    }

  }

}
