% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RADIN_ecto.R
\encoding{UTF-8}
\name{RADIN_ecto}
\alias{RADIN_ecto}
\title{RADIN_ecto}
\usage{
RADIN_ecto(
  ATOT = 0.01325006,
  AV = 0.001325006,
  AT = 0,
  FATOSK = 0.4,
  FATOSB = 0.4,
  FATOBJ = 0,
  EMISAN = 0.95,
  EMISSB = 0.95,
  EMISSK = 0.8,
  TSKY = 10,
  TGRD = 30
)
}
\arguments{
\item{ATOT}{total body surface area (m2)}

\item{AV}{ventral surface area (m2)}

\item{AT}{body surface area contacting another organism of same temperature (m2)}

\item{FATOSK}{configuration factor to sky (-)}

\item{FATOSB}{configuration factor to substrate (-)}

\item{FATOBJ}{configuration factor to nearby object (-) (not functional at the moment)}

\item{EMISAN}{emissivity of animal (fractional, 0-1)}

\item{EMISSB}{emissivity of substrate (fractional, 0-1)}

\item{EMISSK}{emissivity of sky (fractional, 0-1)}

\item{TSKY}{sky temperature (°C)}

\item{TGRD}{ground temperature (°C)}
}
\description{
R version of Fortran RADIN.f (ectotherm model) for calculating absorbed long-wave radiation.
}
