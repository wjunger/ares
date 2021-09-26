\name{admrio}
\alias{admrio}
\docType{data}
\title{Hospital Admissions in Rio de Janeiro}
\description{
This is a typical dataset used in environment air pollution epidemiology. It contains outcome, confounders and exposure variables, i.e., event counts, meteorological condition and air pollutants concentrations.
}
\usage{data(admrio)}
\format{
  A data frame with 730 observations on the following 20 variables.
  \describe{
    \item{\code{date}}{dates string formatted as \emph{dd/mm/yyyy}}
    \item{\code{resp65}}{counts respiratory hospital admissions over 65 years old}
    \item{\code{circ65}}{counts of circulatory hospital admissions over 65 years old}
    \item{\code{copd65}}{counts of COPD hospital admissions over 65 years old}
    \item{\code{pneum65}}{counts of pneumonia hospital admissions over 65 years old}
    \item{\code{myinf65}}{counts of myocardial infarction hospital admissions over 65 years old}
    \item{\code{resp5}}{counts of respiratory hospital admissions under 5 years old}
    \item{\code{pneum5}}{counts of pneumonia hospital admissions under 5 years old}
    \item{\code{warmseason}}{an indicator variable for months with higher temperature}
    \item{\code{tmpmed}}{mean of measured temperature}
    \item{\code{tmpmin}}{minimum of measured temperature}
    \item{\code{tmpmax}}{maximum of measured temperature}
    \item{\code{humid}}{mean of measured relative humidity}
    \item{\code{rain}}{daily rainfall in \emph{mm}}
    \item{\code{rainy}}{an indicator for days with rainfall over 50 milimetres}
    \item{\code{pm10}}{daily averaged concentrations of particulate matter up to 10 microns in micrograms per metre cubic}
    \item{\code{so2}}{daily averaged concentrations of sulphur dioxide in micrograms per metre cubic}
    \item{\code{co}}{daily averaged concentrations of carbon monoxide in micrograms per metre cubic}
    }
}
\details{
This dataset is a subset of a larger dataset used in an actual study in Rio de Janeiro.
}
% \source{
%
%}
% \references{
% Air quality and effects of the Rio de Janeiro municipality population - Final research report (2005). It can be downloaded from \url{http://www.ims.uerj.br/ares-rio/download/relatorio_ares-rio_mma_jan2005.pdf} - Portuguese version only.
%}
\examples{
data(admrio)
names(admrio)
}
\keyword{datasets}
