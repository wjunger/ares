lspline <- function(var,knots=NULL,nknots=NULL,percentiles=NULL,marginal=FALSE,names=NULL)
# generate the basis for a piecewise linear spline
# var : variable
# knots : vector positions for the knots
# nknots : (optional) number of knots
# percentiles : (optional) vector of percentiles for knots positions
# marginal : (optional) marginal effect
# names : (optional) names for the new created variables
{
# initialization
n <- length(var)
if (!is.null(nknots))
    if ((nknots > n/2) || (length(knots) > n/2) || (length(percentiles) > n/2))
        stop("Less than two points between two knots")
if (is.null(nknots) && is.null(knots) && is.null(percentiles))
    stop("Either knots position, number of knots, or percentiles for the knots must be supplied")
if (sum(!is.null(nknots),!is.null(knots),!is.null(percentiles))>1)
    stop("Knots position, number of knots and percentiles are mutually excludent")
if (!is.null(knots))
	{
	if (any(knots<min(var,na.rm=TRUE)) || any(knots>max(var,na.rm=TRUE)))
		stop("At least one knot is out of bounds")
	w <- "knots"
    k <- length(knots)+1
    }
else if (!is.null(percentiles))
	{
	w <- "percentiles"
	k <- length(percentiles)+1
	}
else
	{
	w <- "nknots"
    k <- nknots+1
    }
# if knots position are not supplied then create 
if (w=="nknots")
	knots <- quantile(var,probs=seq(0,1,by=1/k)[2:k],na.rm=TRUE)
else if (w=="percentiles")
	knots <- quantile(var,probs=percentiles,na.rm=TRUE)
# else knots <- knots # just for reminder    

# creating new variables
knots <- c(knots,max(var,na.rm=TRUE))
if (marginal)
    {
    new.vars <- var
    for (i in 2:k)
        {
        tmp.var <- double(n)
        tmp.var <- pmax(0,(var-knots[i-1]))
        new.vars <- cbind(new.vars,tmp.var)
        }
    }
else
    {
    new.vars <- pmin(var,knots[1],na.rm=FALSE)
    for (i in 2:k)
        {
        tmp.var <- double(n)
        tmp.var <- pmax(pmin(var,knots[i],na.rm=FALSE),knots[i-1],na.rm=FALSE)-knots[i-1]
        new.vars <- cbind(new.vars,tmp.var)
        }
    }
# naming new variables
if (length(names) == k)
    dimnames(new.vars)[[2]] <- names
else
    dimnames(new.vars)[[2]] <- paste(deparse(substitute(var)),".",1:k,sep="")
# returning
retval <- new.vars
class(retval) <- "lspline"
return(retval)
}


sincos <- function(period,n=NULL,largest.period=365)
# generate a sinusoidal curve
{
if (period<2)
	stop("Period less than 2")
if (period>largest.period)
	stop("Period is larger than the largest period")
if (is.null(n))
	n <- length(.aresEnv$ares.selection)
k <- largest.period/period
wk <- 2*pi*(k/largest.period)
t<- seq(1,n)
curve <- cbind(sin(wk*t),cos(wk*t))
colnames(curve) <- c(paste("sin",period,sep=""),paste("cos",period,sep=""))
return(curve)
}
