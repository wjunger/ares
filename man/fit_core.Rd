\name{fit_core}
\alias{fit_core}
\title{Fit Core Model}
\description{
Fit the core model
}
\usage{
fit_core(formula, class = "gam", ...)
}
\arguments{
  \item{formula}{the formula for the model. See \code{\link{formula}}}
  \item{class}{a class for the model to be used in the fitting process. Default is \code{gam}. See Details}
\item{\dots}{further arguments for the model fitting function}
}
\details{
This function is primarily a wrapper for the real fitting function. For now, \code{gam} and \code{glm} are the allowed classes. However, other fitting methods can be added in the future. Some information is added to the resulting object and therefore the class \code{ares} is added to its inheritance. For details on GAM options and output see \code{\link{gam}} for GLM options and output see \code{\link{glm}}.
}
\value{
An object of the same class as the one specified by \code{class}. See each class help for a list of objects returned.
}
\references{
Schwartz, J., Spix, C., Touloumi, G. et al. (1996) Methodological issues in studies of air pollution and daily counts of deaths or hospital admissions. \emph{J Epidemiol. Community Health} \bold{50 (suppl 1)}, S12--S18.

McGullagh, P., Nelder, J. A. (1989) \emph{Generalized linear models}. Chapman and Hall.

Hastie, T., Tibshirani, R. (1990) \emph{Generalized additive models}. 2 ed. Chapman and Hall.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{plot}},\code{\link{smooth.spline}},\code{\link{setup}}}
\examples{
data(admrio)
setup(admrio,"date")
f <- resp5~s(time)+weekdays+s(tmpmax)+s(humid)
m <- fit_core(f)
print_summary(m)
}
\keyword{regression}
