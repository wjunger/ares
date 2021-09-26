\name{bubbleplot}
\alias{bubbleplot}
\title{Bubble Plot}
\description{
Produce a bubble type plot
}
\usage{
bubbleplot(x, y, z, bs=0.1, ...)
}
\arguments{
  \item{x}{a vector of x-axis values}
  \item{y}{a vector of y-axis values}
  \item{z}{a vector of values representing the radius of the bubbles}
  \item{bs}{an expansion factor for the bubbles size}
  \item{\dots}{further options for \code{\link{plot}}}
}
\details{
This function implements a bubble type plot. One example of its use is plotting some characteristics of the population and pass the weights of each observation as the radius of the bubbles.
}
\value{
This function does not return a value.
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}

\seealso{\code{\link{plot}}}

\keyword{hplot}
