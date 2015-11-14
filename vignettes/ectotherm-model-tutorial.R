## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ------------------------------------------------------------------------
library(NicheMapR)

## ---- echo=FALSE---------------------------------------------------------
source('../R/ectotherm.R')

## ------------------------------------------------------------------------
micro<-micro_global(loc="Kuranda, Queensland")
ecto<-ectotherm()

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(ecto$environ[,1:13], 12), digits = 2)
knitr::kable(head(ecto$environ[,14:22], 12), digits = 2)

