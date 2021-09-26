\name{pdf_report}
\alias{pdf_report}
\title{Output a report to PDF}
\description{
Output the analysis results and diagnostics to a PDF file
}
\usage{
pdf_report(model, file, pollutants, method = "pdlm", labels = toupper(pollutants),
	unit=10, outcome.label = NULL, city = NULL, df = 0, ...)
}
\arguments{
  \item{model}{a model fitted by \code{\link{fit_core}}}
  \item{file}{a string with a file name to save the output, i.e., \code{"output.pdf"}}
  \item{pollutants}{a vector with the names of the variables to estimate the effects}
  \item{method}{estimation method to be used to estimate the pollutant effect. Default is \code{pdlm} for distributed lag models. Use \code{both} for both single lag and PDLM estimates}
  \item{labels}{a vector of quoted strings with alternate labels for the pollutants. Default is the names of the variables in \code{pollutant}}
  \item{unit}{a vector indicating the units for relative risk computation. Default is 10 for all pollutants. See Details}
  \item{outcome.label}{an alternate label for the outcome variable}
  \item{city}{a string indicating the city. It is for the header}
  \item{df}{the number  of degrees of freedom for the outcome variable smoothing. Default is 0 for no smoothing}
  \item{\dots}{further options for \code{\link{estimate_risks}}}
}
\details{
This function outputs several diagnostics plots and statistics of the core model. It also outputs the effects estimates. The results are piped to a PDF device. See \code{\link{pdf}} for details.

}
\value{
This function does not return a value.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

Schwartz, J. (2000) The distributed lag between air pollution and daily deaths. \emph{Epidemiology} \bold{11(3)}, 320--326.

McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{fit_core}},\code{\link{print_risk}},\code{\link{plot_risk}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)

## pdlm effect estimation
pdf_report(m,"test.pdf",pollutants=c("so2","pm10"),labels=c("SO2","PM10"),unit=10)
}

\keyword{models}
\keyword{regression}
