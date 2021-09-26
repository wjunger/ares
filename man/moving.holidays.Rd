\name{moving.holidays}
\alias{moving.holidays}
\title{Moving Holidays}
\description{
Compute moving holidays based on Easter
}
\usage{
moving.holidays(year)
}
\arguments{
  \item{year}{a numeric or character vector of 4-digits years. It will be coerced to numeric}
}
\details{
Ecclesiastical holidays like Easter, Carnaval, Passion and Corpus Christi are efficiently computed.  
}
\value{
This function returns a list of vectors with the same length as number of years for each date described in Details.
}
\references{
Astronomical Society of South Australia at \url{http://www.assa.org.au/edm.html}
}
\author{Washington Junger \email{wjunger@ims.uerj.br} and Antonio Ponce de Leon \email{ponce@ims.uerj.br}

}
\seealso{\code{\link{attach}},\code{\link{search}},\code{\link{as.Date}}}
\keyword{data}
\keyword{datagen}
