% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COND_ecto.R
\encoding{UTF-8}
\name{COND_ecto}
\alias{COND_ecto}
\title{COND_ecto}
\usage{
COND_ecto(AV = 0.001325006, DEP = 0.025, TSKIN = 25, TSUBST = 10, SUBTK = 0.1)
}
\arguments{
\item{AV}{ventral surface area (m2)}

\item{DEP}{depth of substrate across which conduction is occurring (m)}

\item{TSKIN}{solar absorptivity of ventral fur (fractional, 0-1)}

\item{TSUBST}{substrate temperature (°C)}

\item{SUBTK}{skin temperature (°C)}
}
\description{
R version of Fortran COND.f (ectotherm model) for calculating conductive heat exchange.
}
