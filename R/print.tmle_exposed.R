#'@title Print function for \code{tmle_exposed}
#'@description Print function for \code{tmle_exposed} a targeted minimum-loss based estimator for stochastic interventions among the exposed.
#'@name print.tmle_exposed
#'@author Amalie Lykkemark Moller \email{amalielykkemark@@live.dk}
#'@rdname print.tmle_exposed
#'@method print tmle_exposed
#'@export
print.tmle_exposed <- function(x,...) {
  if(identical(class(x), "tmle_exposed")){

    table<-data.table(data.frame(Estimate=unlist(x$estimate),
                                 SE=unlist(x$se)))
    table[,CI95lower:=Estimate-1.96*SE]
    table[,CI95upper:=Estimate+1.96*SE]

    if(!is.null(x$estimate$psi)){
      table[,Parameter:=c('Risk under stochastic intervention','Risk under no intervetion','Risk difference')]
      table<-table[,.(Parameter,Estimate,SE,CI95lower,CI95upper)]

    }

    if(!is.null(x$estimate$psistar)){
      table[,Parameter:=c('Risk under deterministic intervention','Risk under no intervetion','Risk difference')]
      table<-table[,.(Parameter,Estimate,SE,CI95lower,CI95upper)]
    }

    print(table)
  }}
