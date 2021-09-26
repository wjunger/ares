# R package for air pollution and health effects analysis
# inspired on S+ Steve toolbar
# AreS-Rio Program
# IMS - UERJ - Brasil
# written by Washington Junger <wjunger@ims.uerj.br>
# started on 30/04/2003
#---------------------------------
# last changes list: (date first)
# 30/05/2003	functions for daily mean of pollutants 
# 30/05/2003	embedding of EM algorithm routines
# 05/02/2004	minor changes in envelope routine
# 17/06/2004    update for R 1.9.0
# 04/08/2004	dailymean runs for meteorological means too
# 22/10/2004	EM imputation routines removed
# 12/08/2004	Complete redesign
# 18/05/2006	ported to use gam package, several bugs removed
# --------------------------------------------------------------
# 14/09/2006	Reescrita em portugues para o VIGIAR
# 10/05/2007	Back to English for the ESCALA Project
# 13/06/2007	Finished documentation, some improves, and catched some bugs. It's working! version 0.5.0
# ...
# 26/11/2007	No more dependence on GhostScript
# ...
# 30/04/2008	Added AIC in the temperature/humidity exposure/response plot and save model info after calling estimate.risks()
# ...
# 03/05/2009	exposure-response curves, overdispersion moved to estimate.risks so AIC can be computed in explore.*, estimate.risks is basically a generic function (with methods)
# 24/08/2010	house keeping and make-up for the workshop at ISEE 2010. Improvement in the docs. Namespace. Some minor changes. Version: 0.7.0
# after some time off
# 01/04/2014	renaming, fix new R CRAN restrictions, minor bugs, etc version 0.8.0
#
#---------------------------------
# To do list:
# likelihood ratio test
# residuals bubble plot
# heteroscedasticity test
# enhance plot_evelope performance (it is full of for blocks)
# add autoregressive Poisson
#
# 
# 
# 

# main functions
#---------------------------------

plot_event <- function(x,df=4,gaps=FALSE,type="p",title=NULL,date.format="%d/%m/%Y",new=TRUE,...)
# plot events
{
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
if (missing(x))
	stop("Missing event counts variable")
if(!is.character(x))
	stop("Argument x must be a literal")
# getting what is needed
label <- x
x <- eval(parse(text=x),envir=.aresEnv$ares.active.dataset)
x <- x[.aresEnv$ares.selection==1]
doe <- eval(parse(text='doe'),envir=.aresEnv$ares.active.dataset)
time <- eval(parse(text='time'),envir=.aresEnv$ares.active.dataset)

if(gaps)
	time <- time[.aresEnv$ares.selection==1]
else
	time <- seq(1,length(time[.aresEnv$ares.selection==1]))
doe <- doe[.aresEnv$ares.selection==1]
maxt <- length(time)
set_graph_window(new)
plot(time,x,type=type,xlab="",ylab=label,bg="white",xaxt="n",...)

if (df > 0)
	lines(smooth.spline(na.omit(cbind(time,x)),df=df),col="red",lwd=2,xaxt="n")
axis(1,labels=format(doe,format=date.format)[seq(1,maxt,by=maxt*.05)],at=time[seq(1,maxt,by=maxt*.05)])

if (is.null(title))
	{
	if (df > 0)
		title(main=paste("Daily counts of",label),sub=paste("Spline with",df,"degrees of freedom"),...)
	else
		title(main=paste("Daily counts of",label),...)
	}
else
	title(main=title)
}


plot_pollutant <- function(x,df=4,gaps=FALSE,type="l",title=NULL,date.format="%d/%m/%Y",new=TRUE,...)
# plot pollutants
{
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
if (missing(x))
	stop("Missing pollutant concentration variable")
if(!is.character(x))
	stop("Argument x must be a literal")
# getting what is needed
label <- x
x <- eval(parse(text=x),envir=.aresEnv$ares.active.dataset)
x <- x[.aresEnv$ares.selection==1]
doe <- eval(parse(text='doe'),envir=.aresEnv$ares.active.dataset)
time <- eval(parse(text='time'),envir=.aresEnv$ares.active.dataset)

if(gaps)
	time <- time[.aresEnv$ares.selection==1]
else
	time <- seq(1,length(time[.aresEnv$ares.selection==1]))
doe <- doe[.aresEnv$ares.selection==1]
maxt <- length(time)
set_graph_window(new)
plot(time,x,type=type,xlab="",ylab=label,bg="white",xaxt="n",...)
	
if (df > 0)
	lines(smooth.spline(na.omit(cbind(time,x)),df=df),col="red",lwd=2,xaxt="n")
axis(1,labels=format(doe,format=date.format)[seq(1,maxt,by=maxt*.05)], at=time[seq(1,maxt,by=maxt*.05)])

if (is.null(title))
	{
	if (df > 0)
		title(main=paste("Daily concentrations of",label),sub=paste("Spline with",df,"degrees of freedom"),...)
	else
		title(main=paste("Daily concentrations of",label),...)
	}
else
	title(main=title)
}

plot_pacf <- function(x,lags=25,acf.too=FALSE,type="deviance",new=TRUE,...)
# plot residuals pacf
{
if (missing(x))
	stop("Model is missing")
set_graph_window(new)
resid <- na.exclude(get_residuals(x,type))
if (acf.too)
	{
	par(mfrow=c(1,2))
	res.acf <- acf(resid,lag.max=lags,main=paste("ACF of the residuals of",formula(x)[2]),...)
	}
else
	res.acf <- NULL
res.pacf <- pacf(resid,lag.max=lags,main=paste("PACF of the residuals of",formula(x)[2]),...)

retval <- list(acf=res.acf,pacf=res.pacf)
invisible(retval)
}


plot_qq <- function(x,type="deviance",new=TRUE,...)
# plot residuals qqplot
{
if (missing(x))
	stop("Model is missing")
set_graph_window(new)
resid <- na.exclude(get_residuals(x,type))
qqnorm(resid,main=paste("Normality plot of residuals of",formula(x)[2]),xlab="Standard Normal Quantiles",ylab=attr(resid,"type"),bg="white",...)
qqline(resid,col="red",lty=1,lwd=2)
}


print_summary <- function(x,digits=getOption("digits"),...)
# model summary information
{
scale <- dispersion(x)
pgps <- pgps(x)
res.df <- resdf(x)
deviance <- x$deviance
mod.sum.glm <- summary.glm(x,...)
print(mod.sum.glm)
if (inherits(x,"gam"))
	{
	mod.sum.gam <- summary.gam(x,...)
	mod.sum.pdlm <- NULL
	print(mod.sum.gam)
	}
else if (inherits(x,"pdlm"))
	{
	mod.sum.pdlm <- summary.pdlm(x,...)
	mod.sum.gam <- NULL
	print(mod.sum.pdlm)
	}
else
	{
	mod.sum.pdlm <- NULL
	mod.sum.gam <- NULL
	}
cat("\nResidual deviance",round(deviance,digits),"with",round(res.df,digits),"degrees of freedom.\n",sep=" ")
cat("\nDispersion parameter",round(scale,digits),"based on Pearson\'s statistics (",round(pgps,digits),").\n",sep=" ")
retval <- list(summary.glm=mod.sum.glm,summary.gam=mod.sum.gam,summary.pdlm=mod.sum.pdlm,dipersion=scale,pearson=pgps,residuals.df=res.df,deviance=deviance)
class(retval) <- "summary.ares"
invisible(retval)
}


periodogram <- function(x,type="deviance",print=TRUE,rows=20,test=TRUE,new=TRUE,digits=getOption("digits"))
# plot periodogram of residuals
{
	pgram.iomega <- function(x,n,res)
	# spectral decomposition of residulas
	{
	t <- seq(1:n)
	sp <- ((sum(res*cos(x*t),na.rm=TRUE))^2+(sum(res*sin(x*t),na.rm=TRUE))^2)/n
	return(sp)
	}

if(typeof(x)=="list")
	resid <- get_residuals(x,type)
else
	resid <- x
n <- length(resid)
i <- seq(1:trunc(n/2-1))
t <- seq(1:n)
omega <- (2*pi*i)/n

Iomega <- sapply(i,function(x){pgram.iomega(omega[x],n=n,res=resid)})
period <- (2*pi)/omega
period.max <- round(max(period),2)
period.min <- round(min(period),2)

# plot the periodogram
set_graph_window(new)
plot(omega, Iomega, xlab="Angular frequency (rad) / [Period on the top axis]",ylab="Intensity",bg="white")
axis(3,at=c(min(omega),0.5,1.0,1.5,2.0,2.5,3.0,max(omega)),
	labels=c(period.max,12.57,6.28,4.19,3.14,2.51,2.09,period.min))
title(main=paste("Periodogram of",attr(resid,"type"),sep=" "))
lines(omega,Iomega,type="h")

# display periodogram
periodogram <- cbind.data.frame("period"=period,"frequency"=omega,"intensity"=Iomega)
periodogram <- periodogram[order(periodogram$intensity,decreasing=TRUE),]
if(print)
	print(round(periodogram[1:rows,],digits))
class(periodogram) <- c(class(periodogram),"periodogram")
if(test)
	periodogram_test(periodogram,plot=FALSE)
invisible(periodogram)
}


periodogram_test <- function(object,plot=TRUE)
# test uniformity of a periodogram
{
if(!inherits(object,"periodogram"))
	stop("This object is not a periodogram")
Iomega <- object$intensity    
k <- length(Iomega)
print(ks.test(Iomega,y="punif",alternative="two.sided"))
if(plot)
	{
	plot(ecdf(Iomega),do.points=FALSE,verticals=TRUE,xlab="Intensity",main="Empirical distribution of the periodogram",ylab="F(Intensity)")
	x <- seq(min(Iomega)-1000,max(Iomega)+1000,by=.1)
	lines(x,punif(x,min(Iomega),max(Iomega)),lty=3,col="red")
	}
}            


explore_temp <- function(model,var,df=4,type="deviance",new=TRUE,...)
# plot smoothed residuals against temperature
{
# kept for back compatibility only
explore_meteorology(model,var,df=4,type="deviance",new=TRUE,...)
}


explore_humid <- function(model,var,df=4,type="deviance",new=TRUE,...)
# plot smoothed residuals against humidity
{
# kept for back compatibility only
explore_meteorology(model,var,df=4,type="deviance",new=TRUE,...)
}


explore_meteorology <- function(model,var,df=4,type="deviance",new=TRUE,...)
# plot smoothed residuals against temperature
{
# getting what is needed
resid <- get_residuals(model,type)

if (missing(var))
	stop("Missing meteorology variable")
if(!is.character(var))
	stop("Argument x must be a literal")
	
label <- var
if(exists("ares.active.dataset",envir=.aresEnv) & (var %in% names(.aresEnv$ares.active.dataset)))
var <- eval(parse(text=var),envir=.aresEnv$ares.active.dataset)
else
	stop("Setup() has not been used yet")

set_graph_window(new)
par(mfrow=c(2,3))
new.model <- update(model,paste(".~. + ns(l(",label,",0),",df,")"),sep="")
plot(l(var,0,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," (deg)",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,0,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(l(",label,",1),",df,")"),sep="")
plot(l(var,1,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," lag1",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,1,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(l(",label,",2),",df,")"),sep="")
plot(l(var,2,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," lag2",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,2,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(ma(",label,",1,0),",df,")"),sep="")
plot(ma(var,1,0,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," ma01",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(ma(var,1,0,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(ma(",label,",2,0),",df,")"),sep="")
plot(ma(var,2,0,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," ma02",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(ma(var,2,0,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(ma(",label,",2,1),",df,")"),sep="")
plot(ma(var,2,1,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," ma12",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(ma(var,2,1,TRUE),resid)), df = df),col="red",lwd=2)
par(mfrow=c(1,1))
title(main=paste("Residuals of series",formula(model)[2],"against meteoroly -",df,"d.f."),ylab=attr(resid,"type"))
}


exposure_response <- function(model,var,df=4,type="deviance",new=TRUE,...)
# plot smoothed residuals against exposure
{
# getting what is needed
resid <- get_residuals(model,type)

if (missing(var))
	stop("Missing exposure variable")
if(!is.character(var))
	stop("Argument x must be a literal")
	
label <- var
if(exists("ares.active.dataset",envir=.aresEnv) & (var %in% names(.aresEnv$ares.active.dataset)))
var <- eval(parse(text=var),envir=.aresEnv$ares.active.dataset)
else
	stop("Setup() has not been used yet")

set_graph_window(new)
par(mfrow=c(2,3))
new.model <- update(model,paste(".~. + ns(l(",label,",0),",df,")"),sep="")
plot(l(var,0,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," (lag 0)",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,0,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(l(",label,",1),",df,")"),sep="")
plot(l(var,1,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," lag 1",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,1,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(l(",label,",2),",df,")"),sep="")
plot(l(var,2,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," lag 2",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,2,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(l(",label,",3),",df,")"),sep="")
plot(l(var,3,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," lag 3",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,3,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(l(",label,",4),",df,")"),sep="")
plot(l(var,4,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," lag 4",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,4,TRUE),resid)), df = df),col="red",lwd=2)
new.model <- update(model,paste(".~. + ns(l(",label,",5),",df,")"),sep="")
plot(l(var,5,TRUE),resid,type="p",main="",sub=paste("AIC:",round(AIC(new.model))),xlab=paste(label," lag 5",sep=""),ylab="",bg="white",...)
rug(var,quiet=TRUE)
lines(smooth.spline(na.omit(cbind(l(var,5,TRUE),resid)), df = df),col="red",lwd=2)
par(mfrow=c(1,1))
title(main=paste("Exposure-response between",label,"and",formula(model)[2],"-",df,"d.f."),ylab=attr(resid,"type"))
}


plot_envelope <- function(x,rep=39,type="deviance",new=TRUE,...)
# simulate and plot an envelope of residuals
{
if (rep < 39)
	stop("Number of replications must be at least 39")

if (family(x)$family=="quasi")
	warning("This model is not Poisson. It is a quasi-poisson model")
else if (family(x)$family!="poisson")
	stop("This model is neither Poisson nor Quasi")

formula.core <- as.character(x$formula)
formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3]))
model <- gam(formula,family=x$family,data=x$data,control=x$control)
formula <- as.formula(paste("sim.response",x$formula[1],x$formula[3]))
fitted <- na.exclude(fitted(x))
resid <- na.exclude(get_residuals(x))
n <- length(resid)
e <- matrix(NA,n,rep)

for (i in 1:rep)
	{
	sim.response <- napredict(attr(fitted,"na.action"),rpois(length(fitted),fitted))
	sim.resid <- get_residuals(gam(formula,data=x$data,family=x$family,control=x$control),type)
	e[,i] <- sort(sim.resid)
	}
e1 <- numeric(n)
e2 <- numeric(n)
if (rep == 39)
	for (i in 1:n)
	{
	eo <- sort(e[i,])
	e1[i] <- min(eo)
	e2[i] <- max(eo)
	}
else
	for (i in 1:n)
	{
	eo <- sort(e[i,])
	e1[i] <- eo[round(0.025*rep,0)]
	e2[i] <- eo[round(0.975*rep,0)]
	}
residmean <- apply(e,1,mean)
band <- range(resid,e1,e2)

# plotting the envelope
set_graph_window(new)
qqnorm(e1,axes=FALSE,main="",xlab="",ylab="",type="l",ylim=band,lty=1,lwd=1,col="red",bg="white")
par(new=TRUE)
qqnorm(e2,axes=FALSE,main="",xlab="",ylab="",type="l",ylim=band,lty=1,lwd=1,col="red",,bg="transparent")
par(new=TRUE)
qqnorm(residmean,axes=FALSE,main="",xlab="",ylab="",type="l",ylim=band,lty=2,col="blue",bg="transparent")
par(new=TRUE)
qqnorm(resid,main=paste("Simulated envelope of series",formula(x)[2]), xlab="Standard Normal Quantiles",ylab=attr(resid,"type"),ylim=band,bg="transparent",...)
}


plot_fitted <- function(x,gaps=FALSE,date.format="%d/%m/%Y",new=TRUE,...)
# plot observed and fitted values
{
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
if (missing(x))
	stop("Fitted model is missing")
# getting what is needed
doe <- eval(parse(text='doe'),envir=.aresEnv$ares.active.dataset)
time <- eval(parse(text='time'),envir=.aresEnv$ares.active.dataset)

if(gaps)
	time <- time[.aresEnv$ares.selection==1]
else
	time <- seq(1,length(time[.aresEnv$ares.selection==1]))
doe <- doe[.aresEnv$ares.selection==1]
y <- eval(parse(text=as.character(x$formula[2])),envir=.aresEnv$ares.active.dataset)[.aresEnv$ares.selection==1]
maxt <- length(time)

set_graph_window(new)
plot(time,y,type="p",xlab="",ylab="Observed and fitted values",bg="white",xaxt="n")
lines(time,fitted(x),type="l",col="red",lwd=1,xaxt="n")
axis(1,labels=format(doe,format=date.format)[seq(1,maxt,by=maxt*.05)], at=time[seq(1,maxt,by=maxt*.05)])
title(main=paste("Observed and predicted daily counts of",formula(x)[2]),sub="(fitted values in red)")
}


plot_residuals <- function(x,gaps=FALSE,type="deviance",band=c(-3,3),date.format="%d/%m/%Y",new=TRUE,...)
# plot residuals
{
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
if (missing(x))
	stop("Fitted model is missing")
# getting what is needed
doe <- eval(parse(text='doe'),envir=.aresEnv$ares.active.dataset)
time <- eval(parse(text='time'),envir=.aresEnv$ares.active.dataset)

if(gaps)
	time <- time[.aresEnv$ares.selection==1]
else
	time <- seq(1,length(time[.aresEnv$ares.selection==1]))
doe <- doe[.aresEnv$ares.selection==1]
maxt <- length(time)

set_graph_window(new)
if(inherits(x,"residuals"))
	{
	resid <- x
	varname <- attr(x,"varname")
	}
else
	{
	resid <- get_residuals(x,type)
	varname <- formula(x)[2]
	}
ymin <- min(min(resid,na.rm=TRUE),-4)*1.2
ymax <- max(max(resid,na.rm=TRUE),4)*1.2

plot(time,resid,type="p",ylim=c(ymin,ymax),xlab="",ylab=attr(resid,"type"),,bg="white",xaxt="n",...)
if (!is.null(band))
	{
	lines(time,rep(band[1],length(time)),col="red",lty=1,lwd=2,xaxt="n")
	lines(time,rep(band[2],length(time)),col="red",lty=1,lwd=2,xaxt="n")
	}
axis(1,labels=format(doe,format=date.format)[seq(1,maxt,by=maxt*.05)], at=time[seq(1,maxt,by=maxt*.05)])
title(main=paste("Residuals of series",varname))
invisible(resid)
}


plot_constinfo <- function(x,type="deviance",new=TRUE,...)
# plot residuals against 2 times the squared root of fitted values (measure of constant information)
{
resid <- get_residuals(x,type)

set_graph_window(new)
plot(2*sqrt(fitted(x)),resid,type="p",xlab="2*sqrt(fitted)",ylab=attr(resid,"type"),bg="white",...)
title(main=paste("Constant information on the scale of series",formula(x)[2]))
}


fit_core <- function(formula,class="gam",...)
# fit the core model and save the results
{
if (missing(formula))
	stop("Formula is missing")
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")

if(class=="gam")
	model <- gam(as.formula(formula),family=poisson(link=log),data=subset(.aresEnv$ares.active.dataset,.aresEnv$ares.selection==1),weights=.aresEnv$ares.weights,control=gam.control(trace=FALSE),...)
else if(class=="glm")
	model <- glm(as.formula(formula),family=poisson(link=log),data=subset(.aresEnv$ares.active.dataset,.aresEnv$ares.selection==1),weights=.aresEnv$ares.weights,control=glm.control(trace=FALSE),...)
else
	stop("Model class not allowed yet")

class(model) <- c(class(model),"ares")
return(model)
}


estimate_risks <- function(model,pollutant,unit=10,confidence.level=.95,method="singlelag",perc.rr=TRUE,interaction=NULL,lag.struc=list(l=0:5,ma=NULL),pdlm.struc=list(l=5,d=2,overall=TRUE),overdispersion=FALSE,labels=NULL,print=TRUE,digits=getOption("digits"),plot=TRUE,new=TRUE,graph.scale=FALSE,verbose=TRUE,...)
# estimate the relative risks
{
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
if (missing(model))
	stop("Model is missing")
if (missing(pollutant))
	stop("Pollutant is missing")
else
	if (class(pollutant)=="formula")
		pollutant <- all.vars(pollutant)
if (is.null(labels))
	labels <- toupper(pollutant)
if(length(pollutant)!=length(labels))
	stop("Length of pollutant and labels differs")
if(is.null(unit))
	unit <- perc.range(pollutant,min.perc=10,max.perc=90)
else
	{
	if (length(unit)==1)
		unit <- rep(unit,length(pollutant))
	else
		if (length(unit)!=length(pollutant))
			stop("Length of pollutant and unit differs")
	}
if(!is.null(interaction))
	{
	interaction <- as.factor(eval(parse(text=interaction),envir=.aresEnv$ares.active.dataset))
	if(nlevels(interaction)>2)
		stop("Only 2-level interaction supported")
	}

# calling specific functions
if ((tolower(method)=="singlelag") & (is.null(interaction)))
	# estimate effects using the single lag framework
	estimates <- estimate.risks.singlelag(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
else if (tolower(method)=="singlelag.auto")
	# estimate effects using the single lag framework and controlling for autocorrelation
	stop("Not implemented yet")
	#estimates <- estimate.risks.singlelag.auto(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
else if ((tolower(method)=="singlelag") & (!is.null(interaction)))
	# estimate effects using the single lag framework and interaction term
	estimates <- estimate.risks.singlelag.interaction(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
else if (tolower(method)=="pdlm")
	# estimate effects using the pdlm framework
	estimates <- estimate.risks.pdlm(model,pollutant,pdlm.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
else if (tolower(method)=="dual")
	# estimate effects using dual pollutant models
	estimates <- estimate.risks.dual(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
#else if (tolower(method)=="pdlm.dual")
#	# estimate effects using the pdlm framework
#	estimates <- estimate.risks.pdlm.dual(model,pollutant,pdlm.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
else
	stop("Method not implemented yet")

if(print)
	print_risk(estimates,digits)

# plotting RR
if (plot)
	plot_risk(estimates,labels=attr(estimates,"plot.labels"),new=new,graph.scale=graph.scale,...)

invisible(estimates)
}


estimate.risks.singlelag <- function(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
# estimate effects using the single lag framework
{
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
for (p in pollutant)
	assign(p,.aresEnv$ares.active.dataset[p])
rr2plot <- NULL
rr10902plot <- NULL
npoll <- length(pollutant)
formula.core <- as.character(model$formula)
#	if(missing(lag.struc))
#		stop("Lag structure is missing")
	if(any(duplicated(pollutant)))
		warning("It seems that at least one pollutant is duplicated")
	lags <- lag.struc$l # lags
	mavs <- lag.struc$ma # moving averages
	# need to redefine labels in order to dimensions to agree in plot_risk()
	poll.labels <- labels
	if(!is.null(lag.struc$ma.base))
		ma.base <- lag.struc$ma.base
	else
		ma.base <- 0
	if(any(ma.base>mavs))
		stop("Moving average start lag is greater than the end lag")
	if(!is.null(lag.struc$labels))
		labels <- lag.struc$labels
	else
		{
		labels <- NULL
		if(!is.null(lags))
			labels <- c(labels,paste("Lag",lags))
		if(!is.null(mavs))
			labels <- c(labels,paste("MA",ma.base,mavs))
		}
	nrows <- length(lags)+length(mavs)
	estimates <- array(NA,dim=c(nrows,4,npoll),dimnames=list(labels,c("RR","LBRR","UBRR","p.value"),poll.labels))
	beta.estimates <- array(NA,dim=c(nrows,2,npoll),dimnames=list(labels,c("beta","se"),poll.labels))
	cat("Working...",npoll,"pollutants to estimate risks...\n")
	for (i in 1:npoll)
		{
		if (verbose)
			{
			lag.formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",pollutant[i]))
			cat(i,": ",sep="")
			print(lag.formula,showEnv=FALSE)
			}
		if(!is.null(lags))
			lags.exp <- paste("l(",pollutant[i],", ",lags,")",sep="")
		else
			lags.exp <- NULL
		if(!is.null(mavs))
			mavs.exp <- paste("ma(",pollutant[i],", ",ma.base,", ",mavs,")",sep="")
		else
			mavs.exp <- NULL
		exposure <- c(lags.exp,mavs.exp)
		estimate.temp <- NULL
		beta.estimate.temp <- NULL
		for(j in 1:nrows)
			{	
			formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",exposure[j]))
			if (inherits(model,"gam"))
				run.model <- gam(formula,family=model.family(overdispersion),data=model$data,control=model$control)
			else
				run.model <- glm(formula,family=model.family(overdispersion),data=model$data,control=model$control)
			# last.coeff <- length(coef(run.model))
			beta <- coef(run.model)[exposure[j]]
			se <- sqrt(diag(vcov(run.model)))[exposure[j]]
			poll <- eval(parse(text=pollutant[i]))
			rr <- rr.eval(beta,se,unit[i],confidence.level)
			if(perc.rr)
				rr <- (rr-1)*100
			pvalue <- 2*(1-pnorm(abs(beta/se)))
			estimate.vector <- c(rr,pvalue)
			estimate.temp <- rbind(estimate.temp,t(estimate.vector))
			beta.estimate.temp <- rbind(beta.estimate.temp,c(beta,se))
			}
		estimates[,,i] <- estimate.temp
		beta.estimates[,,i] <- beta.estimate.temp
		}
if (verbose)
	cat("\n")
attr(estimates,"unit") <- unit
attr(estimates,"labels") <- poll.labels
attr(estimates,"plot.labels") <- labels
attr(estimates,"core.model") <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3]))
attr(estimates,"perc.rr") <- perc.rr
attr(estimates,"beta") <- beta.estimates
class(estimates) <- c(class(estimates),"risk","lag.risk")

return(estimates)
}


# estimate.risks.singlelag.auto <- function(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
# # estimate effects using the single lag framework and controlling for autocorrelation
# {	
# if(exists("ares.active.dataset",envir=.aresEnv))
#	ares.active.dataset <- get("ares.active.dataset",envir=.aresEnv)
# else
	# stop("Setup() has not been used yet")
# for (p in pollutant)
	# assign(p,ares.active.dataset[p])
# rr2plot <- NULL
# rr10902plot <- NULL
# npoll <- length(pollutant)
# formula.core <- as.character(model$formula)
# #	if(missing(lag.struc))
# #		stop("Lag structure is missing")
# 	if(any(duplicated(pollutant))) # I don't know how long it should be
# 		warning("It seems that at least one pollutant is duplicated")
# 	lags <- lag.struc$l # lags
# 	mavs <- lag.struc$ma # moving averages
# 	# need to redefine labels in order to dimensions to agree in plot_risk()
# 	poll.labels <- labels
# 	if(!is.null(lag.struc$ma.base))
# 		ma.base <- lag.struc$ma.base
# 	else
# 		ma.base <- 0
# 	if(any(ma.base>mavs))
# 		stop("Moving average start lag is greater than the end lag")
# 	if(!is.null(lag.struc$labels))
# 		labels <- lag.struc$labels
# 	else
# 		{
# 		labels <- paste("Lag",lags)
# 		if(!is.null(mavs))
# 			labels <- c(labels,paste("MA",ma.base,mavs))
# 		}
# 	nrows <- length(lags)+length(mavs)
# 	estimates <- array(NA,dim=c(nrows,4,npoll),dimnames=list(labels,c("RR","LBRR","UBRR","p.value"),poll.labels))
# 	beta.estimates <- array(NA,dim=c(nrows,2,npoll),dimnames=list(labels,c("beta","se"),poll.labels))
# 	cat("Working...",npoll,"pollutants to estimate risks...\n")
# 	for (i in 1:npoll)
# 		{
# 		if (verbose)
# 			{
# 			lag.formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",pollutant[i]),env=.aresEnv)
# 			cat(i,": ",sep="")
# 			print(lag.formula,showEnv=FALSE)
# 			}
# 		if(!is.null(lags))
# 			lags.exp <- paste("l(",pollutant[i],", ",lags,")",sep="")
# 		else
# 			lags.exp <- NULL
# 		if(!is.null(mavs))
# 			mavs.exp <- paste("ma(",pollutant[i],", ",ma.base,", ",mavs,")",sep="")
# 		else
# 			mavs.exp <- NULL
# 		exposure <- c(lags.exp,mavs.exp)
# 		estimate.temp <- NULL
# 		beta.estimate.temp <- NULL
# 		for(j in 1:nrows)
# 			{	
# 			formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",exposure[j]),env=.aresEnv)
# 			if (inherits(model,"gam"))
# 				run.model <- gam(formula,family=model.family(overdispersion),data=model$data,control=model$control)
# 			else
# 				run.model <- glm(formula,family=model.family(overdispersion),data=model$data,control=model$control)
# 			beta <- coef(run.model)[exposure[j]]
# 			se <- sqrt(diag(vcov(run.model)))[exposure[j]]
# 			poll <- eval(parse(text=pollutant[i]))
# 			rr <- rr.eval(beta,se,unit[i],confidence.level)
# 			if(perc.rr)
# 				rr <- (rr-1)*100
# 			pvalue <- 2*(1-pnorm(abs(beta/se)))
# 			estimate.vector <- c(rr,pvalue)
# 			estimate.temp <- rbind(estimate.temp,t(estimate.vector))
# 			beta.estimate.temp <- rbind(beta.estimate.temp,c(beta,se))
# 			}
# 		estimates[,,i] <- estimate.temp
# 		beta.estimates[,,i] <- beta.estimate.temp
# 		}
# if (verbose)
# 	cat("\n")
# attr(estimates,"unit") <- unit
# attr(estimates,"labels") <- poll.labels
# attr(estimates,"plot.labels") <- labels
# attr(estimates,"core.model") <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3]))
# attr(estimates,"perc.rr") <- perc.rr
# attr(estimates,"beta") <- beta.estimates
# class(estimates) <- c(class(estimates),"risk","lag.risk.auto")
# 
# return(estimates)
# }


estimate.risks.singlelag.interaction <- function(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
# estimate effects using the single lag framework and interaction term
{	
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
for (p in pollutant)
	assign(p,.aresEnv$ares.active.dataset[p])
rr2plot <- NULL
rr10902plot <- NULL
npoll <- length(pollutant)
formula.core <- as.character(model$formula)
#	if(missing(lag.struc))
#		stop("Lag structure is missing")
	if(any(duplicated(pollutant)))
		warning("It seems that at least one pollutant is duplicated")
	lags <- lag.struc$l # lags
	mavs <- lag.struc$ma # moving averages
	# need to redefine labels in order to dimensions to agree in plot_risk()
	poll.labels <- labels
	if(!is.null(lag.struc$ma.base))
		ma.base <- lag.struc$ma.base
	else
		ma.base <- 0
	if(any(ma.base>mavs))
		stop("Moving average start lag is greater than the end lag")
	if(!is.null(lag.struc$labels))
		labels <- lag.struc$labels
	else
		{
		labels <- paste("Lag",lags)
		if(!is.null(mavs))
			labels <- c(labels,paste("MA",ma.base,mavs))
		}
	nrows <- length(lags)+length(mavs)
	if(is.null(interaction))
		{
		stop("Interaction term must be specified")
		}
	else
		{
		interaction.levels <- levels(interaction)
		estimates <- array(NA,dim=c(nrows,8,npoll),dimnames=list(labels,c(paste(c("RR","LBRR","UBRR","p.value"),interaction.levels[1],sep="_i="),paste(c("RR","LBRR","UBRR","p.value"),interaction.levels[2],sep="_i=")),poll.labels))
		cat("Working...",npoll,"pollutants to estimate risks...\n")
		for (i in 1:npoll)
			{
			if (verbose)
				{
				lag.formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",pollutant[i],"*",interaction))
				cat(i,": ",sep="")
				print(lag.formula,showEnv=FALSE)
				}
			if(!is.null(lags))
				lags.exp <- paste("l(",pollutant[i],", ",lags,")",sep="")
			else
				lags.exp <- NULL
			if(!is.null(mavs))
				mavs.exp <- paste("ma(",pollutant[i],", ",ma.base,", ",mavs,")",sep="")
			else
				mavs.exp <- NULL
			exposure <- c(lags.exp,mavs.exp)
			estimate.temp <- NULL
			for(j in 1:nrows)
				{	
				formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",exposure[j],"*",interaction))
				if (inherits(model,"gam"))
					run.model <- gam(formula,family=model.family(overdispersion),data=model$data,control=model$control)
				else
					run.model <- glm(formula,family=model.family(overdispersion),data=model$data,control=model$control)
				
				beta <- coef(run.model)[exposure[j]]  # beta of the exposure term
				var <- diag(vcov(run.model))[exposure[j]]
				se <- sqrt(var)
				beta2 <- coef(run.model)[paste(exposure[j],":",interaction,sep="")]  # beta of the interaction term
				var2 <- diag(vcov(run.model))[paste(exposure[j],":",interaction,sep="")]
				covar <- vcov(run.model)[exposure[j],paste(exposure[j],":",interaction,sep="")]
				beta.i <- beta+beta2  # compute total effect
				se.i <- sqrt(var+var2+2*covar)
				poll <- eval(parse(text=pollutant[i]))
				rr <- rr.eval(beta,se,unit[i],confidence.level)
				rr.i <- rr.eval(beta.i,se.i,unit[i],confidence.level)
				if(perc.rr)
					{
					rr <- (rr-1)*100
					rr.i <- (rr.i-1)*100
					}
				pvalue <- 2*(1-pnorm(abs(beta/se)))
				pvalue.i <- 2*(1-pnorm(abs(beta.i/se.i)))
				estimate.vector <- c(rr,pvalue,rr.i,pvalue.i)
				estimate.temp <- rbind(estimate.temp,t(estimate.vector))
				}
			estimates[,,i] <- estimate.temp
			}
		}
if (verbose)
	cat("\n")
attr(estimates,"unit") <- unit
attr(estimates,"labels") <- poll.labels
attr(estimates,"plot.labels") <- labels
attr(estimates,"core.model") <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3]))
attr(estimates,"perc.rr") <- perc.rr
attr(estimates,"interaction") <- interaction
attr(estimates,"interaction.levels") <- interaction.levels
class(estimates) <- c(class(estimates),"risk","lag.risk.interaction")

return(estimates)
}


estimate.risks.pdlm <- function(model,pollutant,pdlm.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
# estimate effects using the pdlm framework
{
#if(missing(pdlm.struc))
#	stop("Polynomial structure is missing")
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
for (p in pollutant)
	assign(p,.aresEnv$ares.active.dataset[p])
npoll <- length(pollutant)
if(any(duplicated(pollutant)))
	warning("It seems that at least one pollutant is duplicated")
lags <- pdlm.struc$l
if (length(pdlm.struc$d)==1)
	degrees <- rep(pdlm.struc$d,npoll)
else
	degrees <- pdlm.struc$d
# need to redefine labels in order to dimensions agree in plot_risk()
poll.labels <- labels
if(!is.null(pdlm.struc$labels))
	labels <- pdlm.struc$labels
else
	labels <- paste("Lag",seq(0,lags))
overall <- ifelse(is.null(pdlm.struc$overall),TRUE,pdlm.struc$overall)
if(overall)
	{
	overall <- TRUE
	labels <- c(labels,"Overall")
	ovr <- 1
	}
else
	{
	overall <- FALSE
	ovr <- 0
	}
formula.core <- as.character(model$formula)
estimates <- array(NA,dim=c(1+lags+ovr,4,npoll),dimnames=list(labels,c("RR","LBRR","UBRR","p.value"),poll.labels))
beta.estimates <- array(NA,dim=c(1+lags+ovr,2,npoll),dimnames=list(labels,c("beta","se"),poll.labels))
cat("Working...",npoll,"pollutants to estimate risks...\n")
for (i in 1:npoll)
	{
	if (verbose)
		{
		cat(i,": ",sep="")
		print(as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+pdl(",pollutant[i],",",lags,",",degrees[i],")")),showEnv=FALSE)
		}
	run.model <- pdlm(model,pollutant[i],lags,degrees[i],family=model.family(overdispersion))
	if(overall)
		{
		beta <- c(run.model$beta$beta,run.model$beta$overall.beta)
		se <- c(run.model$beta$se,run.model$beta$overall.se)
		}
	else
		{
		beta <- run.model$beta$beta
		se <- run.model$beta$se
		}
	poll <- eval(parse(text=pollutant[i]))
	rr <- rr.eval(beta,se,unit[i],confidence.level)
	if(perc.rr)
		rr <- (rr-1)*100
	pvalue <- 2*(1-pnorm(abs(beta/se)))
	estimates[,,i] <- cbind(rr,pvalue)
	beta.estimates[,,i] <- cbind(beta,se)
	}
attr(estimates,"unit") <- unit
attr(estimates,"labels") <- poll.labels
attr(estimates,"plot.labels") <- labels
attr(estimates,"core.model") <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3]))
attr(estimates,"perc.rr") <- perc.rr
attr(estimates,"beta") <- beta.estimates
class(estimates) <- c(class(estimates),"risk","pdlm.risk")

return(estimates)
}


estimate.risks.dual <- function(model,pollutant,lag.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
# estimate effects using dual pollutant models
{	
	get.pairs <- function(x)
	# generate labels given a list of pairs
	{
	l <- NULL
	for (i in 1:length(x))
		l <- c(l,paste(x[[i]][1],x[[i]][2],sep="+"))
	return(l)
	}
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
for (p in pollutant)
	assign(p,.aresEnv$ares.active.dataset[p])
estimates <- NULL
rr2plot <- NULL
rr10902plot <- NULL
formula.core <- as.character(model$formula)
#if(missing(lag.struc))
#	stop("Lag structure is missing")
if(any(duplicated(pollutant))) 
	warning("It seems that at least one of the pollutants is duplicated")
lags <- lag.struc$l # lags
mavs <- lag.struc$ma # moving averages
# need to redefine labels in order to dimensions to agree in plot_risk()
single.labels <- labels # labels for each pollutant
if(!is.null(lag.struc$ma.base))
	ma.base <- lag.struc$ma.base
else
	ma.base <- 0
if(any(ma.base>mavs))
	stop("Moving average start lag is greater than the end lag")
if(!is.null(lag.struc$labels))
	labels <- lag.struc$labels
else
	labels <- c(paste("Lag",lags),paste("MA",ma.base,mavs))
nrows <- length(lags)+length(mavs)

pollutant <- combn(pollutant,2,simplify=FALSE)
#pollutant <- get.pairs(poll.pairs)
lab.pairs <- combn(single.labels,2,simplify=FALSE)
poll.labels <- get.pairs(lab.pairs)
npoll <- length(pollutant)
estimates <- array(NA,dim=c(nrows,4,2*npoll),dimnames=list(labels,c("RR","LBRR","UBRR","p.value"),c(paste(poll.labels,1),paste(poll.labels,2))))
cat("Working...",npoll,"pollutants to estimate risks...\n")
for (i in 1:npoll)
	{
	if (verbose)
		{
		lag.formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",pollutant[[i]][1],"+",pollutant[[i]][2]))
		cat(i,": ",sep="")
		print(lag.formula,showEnv=FALSE)
		}
	if(!is.null(lags))
		lags.exp <- paste("l(",pollutant[[i]][1],", ",lags,")","+l(",pollutant[[i]][2],", ",lags,")",sep="")
	if(!is.null(mavs))
		mavs.exp <- paste("ma(",pollutant[[i]][1],", ",ma.base,", ",mavs,")","+ma(",pollutant[[i]][2],", ",ma.base,", ",mavs,")",sep="")
	exposure <- c(lags.exp,mavs.exp)
	estimate.temp1 <- NULL
	estimate.temp2 <- NULL
	for(j in 1:nrows)
		{	
		formula <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+",exposure[j]))
		if (inherits(model,"gam"))
			run.model <- gam(formula,family=model.family(overdispersion),data=model$data,control=model$control)
		else
			run.model <- glm(formula,family=model.family(overdispersion),data=model$data,control=model$control)
		index=unlist(strsplit(exposure[j],"+",fixed=TRUE)) # split dual exposure into vector
		# first pollutant
		beta <- coef(run.model)[index[1]]
		se <- sqrt(diag(vcov(run.model)))[index[1]]
		poll <- eval(parse(text=pollutant[[i]][1]))
		rr <- rr.eval(beta,se,unit[i],confidence.level)
		if(perc.rr)
			rr <- (rr-1)*100
		pvalue <- 2*(1-pnorm(abs(beta/se)))
		estimate.vector.p1 <- c(rr,pvalue)
		estimate.temp1 <- rbind(estimate.temp1,t(estimate.vector.p1))
		# second pollutant
		beta <- coef(run.model)[index[2]]
		se <- sqrt(diag(vcov(run.model)))[index[2]]
		poll <- eval(parse(text=pollutant[[i]][2]))
		rr <- rr.eval(beta,se,unit[i],confidence.level)
		if(perc.rr)
			rr <- (rr-1)*100
		pvalue <- 2*(1-pnorm(abs(beta/se)))
		estimate.vector.p2 <- c(rr,pvalue)
		estimate.temp2 <- rbind(estimate.temp2,t(estimate.vector.p2))
		}
	estimates[,,i] <- estimate.temp1
	estimates[,,i+npoll] <- estimate.temp2
	}
if (verbose)
	cat("\n")
attr(estimates,"unit") <- unit
attr(estimates,"labels") <- poll.labels
attr(estimates,"plot.labels") <- labels
attr(estimates,"lab.pairs") <- lab.pairs
attr(estimates,"core.model") <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3]))
attr(estimates,"perc.rr") <- perc.rr
class(estimates) <- c(class(estimates),"risk","dual.risk")

return(estimates)
}


#estimate.risks.pdlm.dual <- function(model,pollutant,pdlm.struc,labels,confidence.level,unit,perc.rr,interaction,overdispersion,verbose)
# estimate effects using the pdlm framework
#{
#	get.pairs <- function(x)
#	# generate labels given a list of pairs
#	{
#	l <- NULL
#	for (i in 1:length(x))
#		l <- c(l,paste(x[[i]][1],x[[i]][2],sep="+"))
#	return(l)
#	}
# if(exists("ares.active.dataset",envir=.aresEnv))
	# ares.active.dataset <- get("ares.active.dataset",envir=.aresEnv)
# else
	# stop("Setup() has not been used yet")
# for (p in pollutant)
	# assign(p,ares.active.dataset[p])
#estimates <- NULL
#formula.core <- as.character(model$formula)
#if(missing(pdlm.struc))
#	stop("Polynomial structure is missing")
#npoll <- length(pollutant)
#if(any(duplicated(pollutant)))
#	warning("It seems that at least one pollutant is duplicated")
#lags <- pdlm.struc$l
#if (length(pdlm.struc$d)==1)
#	degrees <- rep(pdlm.struc$d,npoll)
#else
#	degrees <- pdlm.struc$d
# need to redefine labels in order dimensions to agree in plot_risk()
#poll.labels <- labels
#if(!is.null(pdlm.struc$labels))
#	labels <- pdlm.struc$labels
#else
#	labels <- paste("Lag",seq(0,lags))
#overall <- ifelse(is.null(pdlm.struc$overall),TRUE,pdlm.struc$overall)
#if(overall)
#	{
# 	overall <- TRUE
# 	labels <- c(labels,"Overall")
# 	ovr <- 1
# 	}
# else
# 	{
# 	overall <- FALSE
# 	ovr <- 0
# 	}
# pollutant <- combn(pollutant,2,simplify=FALSE)
# #pollutant <- get.pairs(poll.pairs)
# lab.pairs <- combn(single.labels,2,simplify=FALSE)
# poll.labels <- get.pairs(lab.pairs)
# npoll <- length(pollutant)
# 
# estimates <- array(NA,dim=c(1+lags+ovr,4,2*npoll),dimnames=list(labels,c("RR","LBRR","UBRR","p.value"),poll.labels))
# beta.estimates <- array(NA,dim=c(1+lags+ovr,2,2*npoll),dimnames=list(labels,c("beta","se"),poll.labels))
# cat("Working...",npoll,"pollutants to estimate risks...\n")
# for (i in 1:npoll)
# 	{
# 	if (verbose)
# 		{
# 		cat(i,": ",sep="")
# 		print(as.formula(paste(formula.core[2],formula.core[1],formula.core[3],"+pdl(",pollutant[i],",",lags,",",degrees[[i]][1],")","+pdl(",pollutant[[i]][2],",",lags,",",degrees[i],")")),showEnv=FALSE)
# 		}
# 	run.model <- pdlm(model,pollutant[i],lags,degrees[i])
# 	if(overall)
# 		{
# 		beta <- c(run.model$beta$beta,run.model$beta$overall.beta)
# 		se <- c(run.model$beta$se,run.model$beta$overall.se)
# 		}
# 	else
# 		{
# 		beta <- run.model$beta$beta
# 		se <- run.model$beta$se
# 		}
# 	poll <- eval(parse(text=pollutant[i]))
# 	rr <- rr.eval(beta,se,unit[i],confidence.level)
# 	if(perc.rr)
# 		rr <- (rr-1)*100
# 	pvalue <- 2*(1-pnorm(abs(beta/se)))
# 	estimates[,,i] <- cbind(rr,pvalue)
# 	beta.estimates[,,i] <- cbind(beta,se)
# 	}
# attr(estimates,"unit") <- unit
# attr(estimates,"labels") <- poll.labels
# attr(estimates,"plot.labels") <- labels
# attr(estimates,"core.model") <- as.formula(paste(formula.core[2],formula.core[1],formula.core[3]))
# attr(estimates,"perc.rr") <- perc.rr
# attr(estimates,"beta") <- beta.estimates
# class(estimates) <- c(class(estimates),"risk","pdlm.dual.risk")
# 
# return(estimates)
# }
	
	
print_risk <- function(x,digits=getOption("digits"),...)
# print a risk object
{
if(missing(x))
	stop("Risk object is missing")
model.formula <- attr(x,"core.model")
poll.labels <- attr(x,"labels") 
lab.pairs <- attr(x,"lab.pairs")
npoll <- length(poll.labels)

cat("\nCore model formula:\n")
print(model.formula,showEnv=FALSE)
cat("\nEffect estimates:\n")
for (i in 1:npoll)
	{
	if(inherits(x,"dual.risk"))
		{
		cat("\nDual pollutant model:",poll.labels[i],"\n",sep=" ")
		cat("First pollutant effect:",lab.pairs[[i]][1],"\n",sep=" ")
		print(round(x[,,i],digits))
		cat("\nSecond pollutant effect:",lab.pairs[[i]][2],"\n",sep=" ")
		print(round(x[,,i+npoll],digits))
		}
	else 
		{
		cat("\nPollutant:" ,poll.labels[i],"\n",sep=" ")
		if(inherits(x,"lag.risk.interaction"))
			cat("Levels of the interaction term:",paste(attr(x,"interaction.levels"),collapse=" and "),"\n",sep=" ")
		print(round(x[,,i],digits))
		}
	}
}


plot_cook <- function(x,type="deviance",line=0.1,new=TRUE,...)
# plot Cook's distance
{
if (missing(x))
	stop("Fitted model is missing")
dat <- na.omit(cbind(res=resid(x),hat=hatvalues(x)))
r <- dat[,1]
h <- ifelse(dat[,2]<1,dat[,2],dat[,2]-0.001)
p <- x$rank+sum(x$nl.df)
D <- (r^2/p)*(h/(1-h)) # Cook and Weisberg (1982) p.117 # I think it can be improved, but I don't have time right now
n <- length(D)
order <- seq(1:n)
if (!is.null(line))
	max <- max(line,max(D))*1.1
else
	max <- max(D)*1.1

set_graph_window(new)
plot(order,D,type="p",ylim=c(0,max),xlab="Observation",ylab="Distance",bg="white")
if (!is.null(line))
	lines(order,rep(line,n),type="l",col="red",lwd=2)
title(main=paste("Cook\'s distance of observations of series",formula(x)[2]))
invisible(D)
}


fillup_hours <-
function(input,date.time,readings=24,cycle=seq(0,23,by=24/readings),start.date=NULL,end.date=NULL,date.format="%d/%m/%Y",hour.format="%H:%M:%S",offset=0,sep=" ")
# fill up missing hours on the environmental dataset
{
if(missing(input))
       stop("Missing input data")
if(missing(date.time))
       stop("Date/time column name is missing")
if(offset>=60)
 stop("Offset must be less than one hour")
# input
input <- as.data.frame(input)

dt.index <- match(table=colnames(input),x=date.time)
# colnames(input) <- c(date.time,colnames(input)[-dt.index])
input[dt.index] <-
as.character(strptime(input[,dt.index],format=paste(date.format,hour.format,sep=sep)))
# reference data frame
if (is.null(start.date))
       start.date <- as.Date(input[1,dt.index],format=date.format)
else
       start.date <- as.Date(start.date,date.format)
if (is.null(end.date))
       end.date <-
as.Date(input[length(input[,dt.index]),dt.index],format=date.format)
else
       end.date <- as.Date(end.date,date.format)
dates <- seq(start.date,end.date+1,by=as.difftime(24/readings,format="%
H",units="hours"))
nhours <- length(dates)
hours <- rep(cycle,nhours/readings)
hours.h <- paste(ifelse(hours<10,"0",""),hours,sep="")
hours.m <- paste(ifelse(offset<10,"0",""),offset,sep="")
hours <- format(paste(hours.h,":",hours.m,sep=""),format=hour.format)
date.time.ref <- as.character(strptime(paste(dates,hours,sep=sep),"%Y-%m-%d %H:%M"))
data.ref <- cbind.data.frame(date.time.ref)
colnames(data.ref) <- c(date.time)
# merging
merged <- merge(data.ref,input,by=date.time,all.x=TRUE)
return(merged)
}


daily_stats <- function(dataset,parameter,first.column=2,date=TRUE,samples=24,statistic="mean",daylight=c(6,19),date.format="%d/%m/%Y")
# compute daily statistics
{
if (missing(dataset))
	stop("Missing dataset")
if (missing(parameter))
	stop("Missing parameter")
dataset <- as.data.frame(dataset)
ndays <- trunc(dim(dataset)[1]/samples)*samples	# ensures that number of hourly obs is multiple of samples
if (substr(statistic,1,3) == "max")
	stat <- max
else if (substr(statistic,1,3) == "min")
	stat <- min
else stat <- mean

# TEMP or HUMID
# 24-hour mean
if (toupper(parameter) %in% c("TEMP","HUMID"))
	{
		dailymean <- function(parvec)
		{
		localmean <- NULL
		for (r in seq(1,ndays,by=samples))	# runs parameter vector submitted
			{
			day <- parvec[r:(r+(samples-1))]	# extracts one day
			localmean <- c(localmean,ifelse(sum(is.na(day))<samples,stat(day,na.rm=TRUE),NA))
			}
		return(localmean)
		}
	}
# PM10, SO2 and NO2
# 24-hour mean
else if (toupper(parameter) %in% c("PM10","SO2","NO2"))
	{
		dailymean <- function(parvec)
		{
		localmean <- NULL
		for (r in seq(1,ndays,by=samples))	# runs parameter vector submitted
			{
			day <- parvec[r:(r+(samples-1))]	# extracts one day
			valid <- sum(!(is.na(day)))	# counts non-missing hours
			localmean <- c(localmean,ifelse(sum(!is.na(day))>=(3/4)*samples,stat(day,na.rm=TRUE),NA))
			}
		return(localmean)
		}
	}
# O3
# daylight max
else if (toupper(parameter) == "O3MAX")
	{
		dailymean <- function(parvec)
		{
		localmean <- NULL
		
		for (r in seq(1,ndays,by=samples))	# runs parameter vector submitted
			{
			day <- parvec[r:(r+(samples-1))]	# extracts one day
			hours <- seq(1,samples,by=24/samples)		# gets hours vector
			#valid <- sum(!(is.na(day)))	# counts non-missing hours
			valid.allday <- !(is.na(day))	# counts non-missing hours
			valid.hours <- 0
			for (v in 1:length(hours))
				if ((valid.allday[v] == TRUE) && (hours[v] >= daylight[1]) && (hours[v] <= daylight[2]))
					valid.hours <- valid.hours + 1 # tests daylight hours
			localmean <- c(localmean,ifelse((valid.hours >= (3/4)*(daylight[2]-daylight[1])),stat(day,na.rm=TRUE),NA))
			}
		return(localmean)
		}
	}
# CO
# max 8-hour running mean 
else if (toupper(parameter) %in% c("CO","O3"))
	{
		if (samples != 24)
			stop("It is not possible to compute 8-hour running mean with less than 24 daily samples")
		dailymean <- function(parvec)
		{
		localmean <- NULL
		for (r in seq(1,ndays,by=24))	# runs parameter vector submitted
			{
			day <- parvec[r:(r+23)]	# extracts one day
			valid <- sum(!(is.na(day)))	# counts non-missing hours
			if (valid >= 18)
				{
				day <- rep(NA,18)	# creates a day vector with 18 positions
				for (m in 1:18)
					{
					local <- parvec[(r+m-1):(r+m-1+7)]	# jumps the hours within the day
					day[m] <- ifelse(sum(!is.na(local))>=6,mean(local,na.rm=TRUE),NA)	# averages and store in position m if it is ok
					}
				localmean <- c(localmean,ifelse(sum(is.na(day))<18,max(day,na.rm=TRUE),NA))
				}
			else
				localmean <- c(localmean,NA)
			}
		return(localmean)
		}
	}
else
	stop("Rule for this parameter is not implemented yet")

polldata <- dataset[,first.column:dim(dataset)[2]]
if (dim(as.data.frame(polldata))[2] == 1)
	meanmatrix <- dailymean(polldata)
else
	meanmatrix <- apply(polldata,2,dailymean)

if (date)
	{
	daysdate <- as.Date(as.character(dataset[seq(1,ndays,by=samples),1]),date.format)	# gets daily dates
	meanmatrix <- cbind.data.frame(daysdate,meanmatrix)
	colnames(meanmatrix) <- c("Date",colnames(polldata))
	}
else
	meanmatrix <- as.data.frame(meanmatrix)
return(meanmatrix)
}


setup <- function(dataset,date.var,weights=NULL,selection=NULL,date.format="%d/%m/%Y",weekday.ref="Sun",holidays=TRUE,...)
# initialize ares environment
{
if (missing(dataset))
	stop("Missing dataset")
if (missing(date.var)) 
	stop("Missing date variable")	

options(na.action="na.exclude",width=120)
d <- dim(dataset)

if (is.null(selection))
	ares.selection <- rep(1,d[1])
else
	{
	if(typeof(selection)=="character")
		ares.selection <- eval(parse(text=paste("dataset$",selection,sep="")))
	else
		ares.selection <- selection
	}
assign("ares.selection",ares.selection,envir=.aresEnv)

if (is.null(weights))
	{
	assign("ares.weights",rep(1,d[1])[ares.selection==1],envir=.aresEnv)
	}
else
	{
	assign("ares.weights",weights[ares.selection==1],envir=.aresEnv)
	}

if (!is.null(date.var))
	{
	lct <- Sys.getlocale("LC_TIME")
	Sys.setlocale("LC_TIME", "C")
	doe <- as.Date(eval(parse(text=paste("dataset$",date.var,sep=""))),format=date.format) 
	time <- seq(1,dim(dataset)[1])
	weekdays <- as.factor(weekdays(doe,abbreviate=TRUE))
	try(weekdays <- relevel(weekdays,weekday.ref))
	months <- as.factor(months(doe,abbreviate=TRUE))
	quarters <- as.factor(quarters(doe,abbreviate=TRUE))
	years <- as.factor(format(doe,format="%Y"))
	if(holidays)
		{
		holidays <- gen_holidays(doe,selection=FALSE,...)
		ares.active.dataset <- cbind.data.frame(dataset,doe,time,weekdays,months,quarters,years,holidays)
		}
	else
		{
		holidays <- NULL
		ares.active.dataset <- cbind.data.frame(dataset,doe,time,weekdays,months,quarters,years)
		}
	Sys.setlocale("LC_TIME",lct)
	}

attr(ares.active.dataset,"var.labels") <- toupper(names(ares.active.dataset))
assign("ares.active.dataset",ares.active.dataset,envir=.aresEnv)
# attach(.aresEnv$ares.active.dataset)
cat("Analysis environment initialized.\n")
cat("The original dataset has",d[1],"rows and",d[2],"columns \n",sep=" ")
# invisible(ares.active.dataset)
}


gen_holidays <- function(date,holidays=NULL,dates=NULL,selection=TRUE,country=NULL)
# generate holidays variables
{
if (missing(date))
	stop("Missing date variable")
if(selection)
	date <- date[.aresEnv$ares.selection==1]
# fixed international holidays
if(is.null(holidays) && is.null(dates))
	{
	intl.holidays <- c("newyear","christmas")
	intl.dates <- c("01/01","25/12")
		# country specific holidays
	if (!is.null(country) )
		{
		# try to find a holiday file in etc. File format is based upon International Standard ISO-3166-1993
		path <- file.path(path.package("ares2",quiet=TRUE)[1],"etc")
		file <- paste(path,"/",toupper(country),".hol",sep="")
		if(!file.exists(file))
			stop("It does not exist a file with holiday dates for this country")
		ctry <- try(read.table(file,sep=",",header=FALSE,colClasses="character",col.names=c("holidays","dates")),silent=TRUE)
			if(inherits(ctry,"try-error"))
				stop("There was a problem opening the holiday file")
		ctry.holidays <- ctry$holidays
		ctry.dates <- ctry$dates
		}
	else if(exists(".holidays",envir=.aresEnv) && exists(".dates",envir=.aresEnv))
		{
		# if file not found, try to find hidden vectors
		ctry.holidays <- get(".holidays",envir=.aresEnv)
		ctry.dates <- get(".dates",envir=.aresEnv)
		}
	else
		{
		# if still nothing, go on without it
		message("Country-specific holidays not found. Continuing without it.")
		ctry.holidays <- NULL
		ctry.dates <- NULL
		}
	holidays <- c(intl.holidays,ctry.holidays)
	dates <- c(intl.dates,ctry.dates)
	extra.holidays <- FALSE
	}
else
	{
	if (length(holidays)!=length(dates))
		stop("Lenghts of holidays and dates differ. It cannot continue.")
	extra.holidays <- TRUE
	}
holidays.matrix <- NULL
holidays.names <- NULL
nholidays <- length(holidays)
for (h in 1:nholidays)
	{
	if(nchar(dates[h])==5)
		hol.format <- "%d/%m"
	else if(nchar(dates[h])==10)
		hol.format <- "%d/%m/%Y"
	else
		stop("Something is wrong in the date format")
	hol.var <- ifelse(format(date,format=hol.format)==dates[h],1,0)
	if (sum(hol.var)>0)
		{
		holidays.matrix <- cbind(holidays.matrix,hol.var)
		holidays.names <- c(holidays.names,h)
		}
	}
if (is.null(holidays.matrix))
	cat("There are no holidays in the given period.\n")
else
	colnames(holidays.matrix) <- holidays[holidays.names]
if(extra.holidays)
	return(holidays.matrix)
# moving holidays
years <- as.integer(levels(as.factor(format(date,format="%Y"))))
nyears <- length(years)
mov.holidays <- moving.holidays(years)
n <- length(date)
mov.holidays.matrix <- NULL
# carnaval
carnaval <- integer(n)
mov.dates <- as.Date(mov.holidays$carnaval,format="%d/%m/%Y")
for (i in 1:nyears)
	carnaval <- carnaval+(date==mov.dates[i])
if(sum(carnaval)>0)
	mov.holidays.matrix <- cbind(mov.holidays.matrix,carnaval)
# thursdayst
thursdayst <- integer(n)
mov.dates <- as.Date(mov.holidays$thursdayst,format="%d/%m/%Y")
for (i in 1:nyears)
	thursdayst <- thursdayst+(date==mov.dates[i])
if(sum(thursdayst)>0)
	mov.holidays.matrix <- cbind(mov.holidays.matrix,thursdayst)
# passion
passion <- integer(n)
mov.dates <- as.Date(mov.holidays$passion,format="%d/%m/%Y")
for (i in 1:nyears)
	passion <- passion+(date==mov.dates[i])
if(sum(passion)>0)
	mov.holidays.matrix <- cbind(mov.holidays.matrix,passion)
# easter
easter <- integer(n)
mov.dates <- as.Date(mov.holidays$easter,format="%d/%m/%Y")
for (i in 1:nyears)
	easter <- easter+(date==mov.dates[i])
if(sum(easter)>0)
	mov.holidays.matrix <- cbind(mov.holidays.matrix,easter)
# corpus
corpus <- integer(n)
mov.dates <- as.Date(mov.holidays$corpus,format="%d/%m/%Y")
for (i in 1:nyears)
	corpus <- corpus+(date==mov.dates[i])
if(sum(corpus)>0)
	mov.holidays.matrix <- cbind(mov.holidays.matrix,corpus)

return(cbind(holidays.matrix,mov.holidays.matrix))
}


moving.holidays <- function(year)
# compute moving holidays (religious)
{
if(typeof(year)=="character")
	year <- as.integer(year)
# finding easter
n1 <- year%%19
n2 <- year%/%100
n3 <- year%%100
n4 <- n2%/%4
n5 <- n2%%4
n6 <- (n2+8)%/%25
n7 <- (n2-n6+1)%/%3
n8 <- (19*n1+n2-n4-n7+15)%%30
n9 <- n3%/%4
n10 <- n3%%4
n11 <- (32+2*n5+2*n9-n8-n10)%%7
n12 <- (n1+11*n8+22*n11)%/%451
month <- (n8+n11-7*n12+114)%/%31
day <- (n8+n11-7*n12+114)%%31+1

easter.day <- as.Date(ISOdate(year,month,day)) 
easter <- format(easter.day,format="%d/%m/%Y")
carnaval <- format(easter.day-47,format="%d/%m/%Y")
thursdayst <- format(easter.day-3,format="%d/%m/%Y")
passion <- format(easter.day-2,format="%d/%m/%Y")
corpus <- format(easter.day+60,format="%d/%m/%Y")

return(list(carnaval=carnaval,thursdayst=thursdayst,passion=passion,easter=easter,corpus=corpus))
}


import_data <- function(file,text.format="csv",...)
# import dataset
{
if(missing(file))
	stop("File name is missing")

flen <- nchar(file)
ext <- get.extension(file)
if(toupper(ext)=="DTA")
	try(dataobj <- read.dta(file))
else if(toupper(ext)=="DBF")
	try(dataobj <- read.dbf(file))
else if(toupper(ext)=="SAV")
	try(dataobj <- read.spss(file,to.data.frame=TRUE,reencode='latin1'))
else if(toupper(ext)=="RDA")
	dataobj <- try(get(load(file)))
else if(tolower(text.format)=="csv")
	dataobj <- read.csv(file,...)
else if (tolower(text.format)=="csv2")
	dataobj <- read.csv2(file,...)
else if (tolower(text.format)=="tab")
	dataobj <- read.delim(savefile,...)
else if (tolower(text.format)=="tab2")
	dataobj <- read.delim2(file,...)
else if (tolower(text.format)=="spc")
	dataobj <- read.table(file,header=TRUE,sep=" ",dec=".",...)
else if (tolower(text.format)=="spc2")
	dataobj <- read.table(file,header=TRUE,sep=" ",dec=",",...)
else
	stop("Sorry! File format not supported. You must import it manually.\nCheck library \'foreign\' for help.")
return(dataobj)
}


export_data <- function(data,file,text.format="csv",...)
# export some data file formats
{
if(missing(data))
	stop("Data object is missing")
if(missing(file))
	stop("File name is missing")

flen <- nchar(file)
ext <-get.extension(file)
if(toupper(ext)=="DTA")
	try(dataobj <- write.dta(data,file,version=7))
else if(toupper(ext)=="DBF")
	try(dataobj <- write.dbf(data,file))
else if(toupper(ext)=="RDA")
	try(save(data,file=file))
# text formats
else if (tolower(text.format)=="csv")
	write.table(data,file,row.names=TRUE,col.names=NA,na="",sep=",",dec=".",...)
else if (tolower(text.format)=="csv2")
	write.table(data,file,row.names=TRUE,col.names=NA,na="",sep=";",dec=",",...)
else if (tolower(text.format)=="tab")
	write.table(data,file,row.names=TRUE,col.names=NA,na="",sep="\t",dec=".",...)
else if (tolower(text.format)=="tab2")
	write.table(data,file,row.names=TRUE,col.names=NA,na="",sep="\t",dec=",",...)
else if (tolower(text.format)=="spc")
	write.table(data,file,row.names=TRUE,col.names=NA,na="",sep=" ",dec=".",...)
else if (tolower(text.format)=="spc2")
	write.table(data,file,row.names=TRUE,col.names=NA,na="",sep=" ",dec=",",...)
else
	stop("Sorry! File format not supported. You must export it manually.\nCheck library \'foreign\' for help.")
}


save_plot <- function(file,width=520,height=480)
# save active plot device
{
if(missing(file))
	stop("File name is missing")
cur <- dev.cur()
dev.set(dev.cur())
ext <- get.extension(file)
if(ext=="pdf")
	dev.copy(pdf,file=file,width=width/72,height=height/72)
else if(ext=="png")
	dev.copy(png,file=file,width=width,height=height,bg="white")
else if(ext=="jpg")
	dev.copy(jpeg,file=file,width=width,height=height)
else if(ext=="svg")
	dev.copy(devSVG,file=file,width=width/72,height=height/72)
else if(ext=="eps")
	{
	dev.copy2eps(file=file,width=width/72,height=height/72)
	return()
	}
else 
	stop("Graph format not supported")
dev.off()
}


get.extension <- function(path)
# get the extension part of a filename
{
parts <- try(strsplit(path, "\\.")[[1]],silent=TRUE)
if (length(parts) > 1)
		return(tolower(parts[length(parts)]))
else
		return("")
}


pgps <- function(model)
# estimate generalized Pearson statistics for Poisson models
{
gps <-sum(((model$y-model$fitted.values)^2)/model$fitted.values)
return(gps)
}


resdf <- function(model)
# extract residual degrees of freedom
{
#res.df <- model$df.null-sum(diag(model$edf))
res.df <- model$df.residual
return(res.df)
}


dispersion <- function(model)
# estimate consistent and asymptotic dispersion parameter based on generalized Pearson statistics for Poisson model
{
dispersion <- pgps(model)/resdf(model)
return(dispersion)
}


get_residuals <- function(model,type="adj_deviance",plot=FALSE,...)
# encapsulated residuals function
{
rd <- resid(model,type="deviance")

if (type == "adj_deviance")
	{
	resid <- rd+1/(6*sqrt(fitted(model)))
	attr(resid,"type") <- "Adjusted deviance residuals"
	}
else if (type == "std_deviance")
	{
	resid <- rd/sqrt(1-hatvalues(model))
	attr(resid,"type") <- "Standardized deviance residuals"
	}
else if (type == "std_scl_deviance")
	{
	resid <- rd/sqrt(dispersion(model)*(1-hatvalues(model)))
	attr(resid,"type") <- "Standardized scaled deviance residuals"
	}
else if (type == "deviance")
	{
	resid <- rd
	attr(resid,"type") <- "Deviance residuals"
	}
else
    stop("Residuals rules not implemented yet")
class(resid) <- c("residuals",class(resid))
attr(resid,"varname") <- formula(model)[2]
if (plot)	
	plot(resid)
return(resid)
}


diagnostics <- function(model,single.graph=TRUE)
# graph and diagnostics
{
if (single.graph)
	par(mfrow=c(2,3),cex.main=.80)
plot_fitted(model,new=!(single.graph))
plot_residuals(model,new=!(single.graph))
plot_cook(model,new=!(single.graph))
plot_pacf(model,lags=25,new=!(single.graph))
periodogram(model,rows=15,new=!(single.graph))
plot_qq(model,new=!(single.graph))
print_summary(model)
}


l <- function(var,k,selection=TRUE)
# lag variable
{
if (missing(k))
	stop("Missing lag index")
if (k!=0)
	{
	n <- length(var)
	na <- rep(NA,k)
	varl <- c(na,var)[1:n]
	}
else
	varl <- var
if(selection)
	return(varl[.aresEnv$ares.selection==1])
else
	return(varl)
}

ma <- function(var,begin,end,selection=TRUE)
# compute moving avarages
{
if(missing(begin)||missing(end))
	stop("Missing either start or end index")
data.matrix <- NULL
for (i in begin:end)
	data.matrix <- cbind(data.matrix,l(var,i,FALSE))
means <- apply(data.matrix,1,mean,na.rm=FALSE)
if (selection)
	return(means[.aresEnv$ares.selection==1])
else
	return(means)
}


rr.eval <- function(beta,se,unit=1,confidence.level=.95)
# compute rr and confidence.level interval
{
z <- abs(qnorm((1-confidence.level)/2))
rr <- exp(unit*beta)
lbrr <- exp(unit*beta - unit*z*se)
ubrr <- exp(unit*beta + unit*z*se)
return(cbind(rr,lbrr,ubrr))
}


plot_risk <- function(x,labels=rownames(x),new=TRUE,graph.scale=FALSE,...)
# plot relative risks in a fashion way
{
	# plotting RR
		doPlot <- function(rr,lbrr,ubrr,new,ref.line,labels,ticks,k,var.name,perc.rr,graph.scale)
		{
		nticks <- 20
		set_graph_window(new)
		stockplot(rr,lbrr,ubrr,ref.line=ref.line,xlabels=labels,ticks=nticks,new=new,main=paste("Relative risk for",k,ifelse(k>1,"units","unit"),"variation of the pollutant"),sub=ifelse(is.null(var.name),"",paste("Pollutant:",var.name)),xlab="Exposure",ylab=paste(ifelse(perc.rr,"%",""),"Relative risk"),graph.scale=graph.scale,...)
		}

if(missing(x))
	stop("Risks object is missing")
if(!inherits(x,"risk"))
	stop("Risks object is not of class risk")
unit <- attr(x,"unit")
#if (length(unique(unit))==1)
#	k <- unique(unit)
#else
#	k <- "k"
perc.rr <- attr(x,"perc.rr") 
ref.line <- ifelse(perc.rr,0,1)
# set ylim if required
# kept separated for the future
if(is.logical(graph.scale))
	{
	if(graph.scale)
		{
		ymin <- min(x[,2,])
		ymax <- max(x[,3,])
		graph.scale <- c(ymin,ymax)
		}
	else
		graph.scale <- NULL
	}
for(opt in 1:dim(x)[3])
	{
	k <- unit[opt]
	if(inherits(x,"pdlm.risk"))
		{
		if(opt==1)
			cat("\nThis is a polynomial distributed lag model.\n")
		var.name <- dimnames(x)[[3]][opt]
		xopt <- as.data.frame(x[,,opt])
		doPlot(xopt[,1],xopt[,2],xopt[,3],new,ref.line,labels,ticks,k,var.name,perc.rr,graph.scale)
		}
	else if (inherits(x,"lag.risk"))
		{
		if(opt==1)
			cat("\nThis is a single lag model.\n")
		var.name <- dimnames(x)[[3]][opt]
		xopt <- as.data.frame(x[,,opt])
		doPlot(xopt[,1],xopt[,2],xopt[,3],new,ref.line,labels,ticks,k,var.name,perc.rr,graph.scale)
		}
	else if (inherits(x,"lag.risk.interaction"))
		{
		if(opt==1)
			cat("\nThis is a single lag model with interaction.\n")
		var.name <- dimnames(x)[[3]][opt]
		xopt <- as.data.frame(x[,,opt])
		doPlot(xopt[,1],xopt[,2],xopt[,3],new,ref.line,labels,ticks,k,paste(var.name,"|",attr(x,"interaction"),"=",attr(x,"interaction.levels")[1],sep=" "),perc.rr,graph.scale)
		doPlot(xopt[,5],xopt[,6],xopt[,7],new,ref.line,labels,ticks,k,paste(var.name,"|",attr(x,"interaction"),"=",attr(x,"interaction.levels")[2],sep=" "),perc.rr,graph.scale)
		}
	else if(inherits(x,"dual.risk"))
		{
		cat("\nThis is a dual pollutant model.\n")
		lab.pairs <- attr(x,"lab.pairs")
		npoll <- length(lab.pairs)
		for (i in 1:npoll)
			cat(i,":",lab.pairs[[i]][1],"in",paste(lab.pairs[[i]][1],lab.pairs[[i]][2],sep="+"),"\n")
		for (i in 1:npoll)
			cat(i+npoll,":",lab.pairs[[i]][2],"in",paste(lab.pairs[[i]][1],lab.pairs[[i]][2],sep="+"),"\n")
		cat("0 :","Quit","\n")
		opt <- as.integer(readline("Select the pollutant to plot: "))
		if(opt==0)
			break
		else if (!(opt%in%seq(1,dim(x)[3])))
			return(cat("Option out of range.\n"))
		else
			{
			var.name <- dimnames(x)[[3]][opt]
			xopt <- as.data.frame(x[,,opt])
			doPlot(xopt[,1],xopt[,2],xopt[,3],new,ref.line,labels,ticks,k,var.name,perc.rr,graph.scale)
			}
		}
	else
		stop("Model not supported")	
	}
}


stockplot <- function(mid,low,high,ref.line=NULL,xlabels=seq(1:length(mid)),ticks=20,mid.pch=19,lim.pch=15,graph.scale=NULL,ylim=NULL,...)
# plot stock type lines
{
l <- length(mid)
if((length(low)!=l)|(length(high)!=l))
	stop("Mid points and end points must have the same size")
x <- seq(1:l)
y <- mid
ymin <- ifelse(is.null(graph.scale),min(low),graph.scale[1])
if (ymin > 1)
	ymin <- 1-(ymin-1)
ymax <- ifelse(is.null(graph.scale),max(high),graph.scale[2])
if (ymax < 1)
	ymax <- 1-(ymax-1)
yax <- as.character(round(seq(ymin,ymax,by=(ymax-ymin)/ticks),3))
if (is.null(ylim))
	ylim <- c(ymin,ymax)
plot(x,y,type="p",pch=mid.pch,xlim=c(1,l),ylim=ylim,col.axis="transparent",bg="white",xaxt="n",yaxt="n",...)
axis(1,at=x,labels=xlabels)
axis(2,at=yax,labels=yax)
segments(x,low,x,high)
points(x,high,type="p",pch=lim.pch)
points(x,low,type="p",pch=lim.pch)
if(!is.null(ref.line))
	{
	#axis(2,at=ref.line,labels=as.character(ref.line),col.axis="red",col="red")
	axis(4,at=ref.line,labels=as.character(ref.line),col.axis="red",col="red")
	abline(ref.line,0,col="red",lwd=2)
	}
}


bubbleplot <- function(x,y,z,bs=0.1,...)
# draw a buble plot
{
z <- z/max(z)
yrange <- max(y)-min(y)
xrange <- max(x)-min(x)
plot(x,y,type="n",...)

for (i in 1:length(x))
  {
  theta <- seq(0,2*pi,pi/200)
  yv <- z[i]*sin(theta)*bs*yrange
  xv <- z[i]*cos(theta)*bs*xrange
  lines(x[i]+xv,y[i]+yv,...)
  }
}


count_na <- function(var)
# output statiscs about missing data in var
{
arguments <- as.list(match.call())
var  <- eval(arguments$var, .aresEnv$ares.active.dataset)
n.total <- length(var)
n.missing <- sum(is.na(var))
n.valid <- n.total-n.missing
na.percent <- (n.missing/n.total)*100
retval <- list("n.total"=n.total,"na"=n.missing,"n.valid"=n.valid,"percent.na"=na.percent)
return(retval)
}


desc_data <- function(dataset=NULL)
# list variables in the ares active dataset
{
	lnames <- function(data)
		{
		d <- dim(data)
		obs <- d[1]
		index <-  seq(1,d[2])
		varnames <- colnames(data)
		if (!is.null(attr(data,"var.labels")))
			labels <- attr(data,"var.labels")
		else
			labels <- toupper(varnames)
		classes <- NULL
		nas <- NULL
		for(i in 1:d[2])
			{
			strvar <- paste("data$",varnames[i],sep="")
			myclass <- class(eval(parse(text=strvar)))
			classes <- c(classes,myclass)
			thisnas <- do.call("count_na",list(eval(parse(text=strvar))))
			nas <- c(nas,thisnas$na)
			}
		desc <- as.data.frame(list("Variable"=varnames,"Class"=classes,"Missing"=nas,"Label"=labels))
		cat("The dataset has",d[2],"variables and",d[1],"observations \n",sep=" ")
		print(desc)
		return(desc)
		}

if (is.null(dataset))
	{
	if (!exists("ares.active.dataset",envir=.aresEnv))
		stop("Neither initialized nor supplied dataset") 
	else
		dataset <- .aresEnv$ares.active.dataset
	}
else
	invisible(lnames(dataset))
}


desc_vars <- function(vars,by=NULL,stats=c("n","na","mean","sd","min","max","centiles"),probs=c(.25,.50,.75),labels=NULL,print=TRUE,digits=getOption("digits"),...)
# compute descriptives statistics of variables
{

if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
else
	dataset <- .aresEnv$ares.active.dataset

if (missing(vars))
	stop("Variables list is missing")
else
	if (class(vars)=="formula")
		vars <- all.vars(vars)
if (is.null(labels))
	labels <- vars
if (length(vars)!=length(labels))
	stop("Length of 'vars' and 'labels' differ")
if (any(!(stats %in% c("n","na","mean","sd","min","max","centiles"))))
	stop("There is at least one alien statistic in 'stats'")

	compute.stats <- function(data,vars,stats)
		{
		# initializing variables
		desc <- NULL
		n <- NULL
		na <- NULL
		mean <-NULL
		sd <- NULL
		min <- NULL
		max <- NULL
		centiles <- NULL
		
		for (i in 1:nvars)
			{
			if (any("n" %in% stats))
				n <- c(n,sum(!is.na(eval(parse(text=paste("data$",vars[i],sep=""))))))
			if (any("na" %in% stats))
				na <- c(na,sum(is.na(eval(parse(text=paste("data$",vars[i],sep=""))))))
			if (any("mean" %in% stats))
				mean <- c(mean,mean(eval(parse(text=paste("data$",vars[i],sep=""))),na.rm=TRUE),...)
			if (any("sd" %in% stats))
				sd <- c(sd,sd(eval(parse(text=paste("data$",vars[i],sep=""))),na.rm=TRUE))
			if (any("min" %in% stats))
				min <- c(min,min(eval(parse(text=paste("data$",vars[i],sep=""))),na.rm=TRUE))
			if (any("max" %in% stats))
				max <- c(max,max(eval(parse(text=paste("data$",vars[i],sep=""))),na.rm=TRUE))
			if (any("centiles" %in% stats))
				{
				centiles <- rbind(centiles,quantile(eval(parse(text=paste("data$",vars[i],sep=""))),probs=probs,na.rm=TRUE))
				colnames(centiles) <- paste("p",probs*100,sep="")
				}
			}
		if (any("n" %in% stats))
			desc <- cbind(desc,n)
		if (any("na" %in% stats))
			desc <- cbind(desc,na)
		if (any("mean" %in% stats))
			desc <- cbind(desc,mean)
		if (any("sd" %in% stats))
			desc <- cbind(desc,sd)
		if (any("min" %in% stats))
			desc <- cbind(desc,min)
		if (any("max" %in% stats))
			desc <- cbind(desc,max)
		if (any("centiles" %in% stats))
			desc <- cbind(desc,centiles)
		desc <- as.data.frame(desc)
		return(desc)
		}
stats <- tolower(stats)
nvars <- length(vars)

if(is.null(by))
	{
	desc <- compute.stats(dataset,vars,stats)
	# putting variable names
	rownames(desc) <- labels
	}
else
	{
	byfactor <- as.factor(dataset[,deparse(substitute(by))])
	bylevels <- levels(byfactor)
	desc <- NULL
	for (i in 1:nlevels(byfactor))
		{
		dataset.tmp <- subset(dataset,byfactor==bylevels[i]) 
		desc.tmp <- compute.stats(dataset.tmp,vars,stats)
		rownames(desc.tmp) <- paste(labels, " [",deparse(substitute(by))," == ",bylevels[i],"]",sep="")
		desc <- rbind.data.frame(desc,desc.tmp)
		}
	}

if(print)
	{
	header <- "Descriptive statistcs for variables:"
	cat(header,"\n")
	if(is.null(by))
		cat(vars,"\n",sep=" ")
	else
		cat(vars,"\nby",deparse(substitute(by)), ": levels = {",bylevels,"}\n",sep=" ")
	print(round(desc,digits))
	}
invisible(desc)
}


ljungbox_test <- function(x,k=25,...)
# perform Ljung and Box test on x
{
# under the null hypothesis:
# H0: rho(1)=rho(2)=...=rho(k)
# k<(n-1)
# Spanos, A. Probability Theory and Statistical Inference. p.749
x.acf <- acf(x,lag.max=k,type="correlation",plot=FALSE,na.action=na.pass,...)
rho <- x.acf$acf[1:k+1]
n <- x.acf$n.used
tau <- 1:k
data.name <- deparse(substitute(x))

LB <- sum(((n*(n+2))/(n-tau))*rho^2)
p.value <- 1-pchisq(LB,k)

# prints chi-square test with k df
cat("\n        Ljung-Box first order independence test\n\n")
cat("data: ",data.name,"\n",sep=" ")
cat("Ljung-Box statistic = ",LB,", df = ",k," p.value = ",format.pval(p.value),"\n",sep="")
cat("Null hypothesis: rho(1) = rho(2) = ... = rho(",k,") = 0\n\n",sep="")
retval <- list(statistic=LB,p.value=p.value,df=k,rho=rho,n.used=n,data.name=data.name)
class(retval) <- "htest"
invisible(retval)
}


whitenoise_test <- function(object,type="deviance",k=25,...)
# test if a vector is white noise
{
if(inherits(object,"lm"))
	resid <- na.exclude(get_residuals(object,type))
else
	resid <- object
n <- length(resid)
# Skewness e Kurtosis
m1 <- (1/n)*sum(resid)
m2 <- (1/n)*sum((resid-m1)^2)
m3 <- (1/n)*sum((resid-m1)^3)
m4 <- (1/n)*sum((resid-m1)^4)
S <- m3/sqrt(m2^3)
K <- m4/m2^2
	
cat("\nTest for normality of residuals\n")
cat("\n        Bowman and Shenton\n")
cat("Skewness = ",S,", N(0,",6/n,") distributed\n",sep="")
cat("Kurtosis = ",K,", N(3,",24/n,") distributed\n",sep="")

cat("\n        Jarque-Bera\n")
JB = n*((S^2/6)+((K-3)^2/24))
p.value = 1-pchisq(JB,2)
cat("\nJarque-Bera statistic = ",JB,", df = ",2,", p.value = ",p.value,"\n")

print(ks.test(resid,y="pnorm",alternative="two.sided"))
print(shapiro.test(resid))

cat("\nTest for independence of residuals\n")
ljungbox_test(resid,k=k,...)
}


pdf_report <- function(model,file,pollutants,method="pdlm",labels=toupper(pollutants),unit=10,outcome.label=NULL,city=NULL,df=0,...)
# print out a report given the core model
{
# assembling information
header <- paste("Ares Library Analysis Report - Version",packageDescription("ares2")["Version"])
doe <- .aresEnv$ares.active.data.frame$doe
period <- paste("From",min(doe),"to",max(doe))
formula <- as.character(model$formula)
outcome <- formula[2]
if (length(unit)==1)
	unit <- rep(unit,length(pollutants))
pdf(file,title=header,version="1.4",paper="a4",width=0,height=0)
textplot(paste(header,paste("Outcome:",outcome.label),paste("City:",city),paste("Period:",period),paste("Pollutants:",paste(labels,collapse=" ")),"Baseline model:",paste(formula[2],formula[1],formula[3],sep=" "),paste("On: ",date()),sep="\n"))
plot_event(outcome,df=df,new=FALSE)
textplot(capture.output(print_summary(model)))
par(mfrow=c(2,3),cex.main=.80)
plot_fitted(model,new=FALSE)
plot_residuals(model,new=FALSE)
plot_cook(model,new=FALSE)
plot_pacf(model,lags=25,new=FALSE)
frequencies <- c("Spectral frequencies ordered by intensity:",capture.output(periodogram(model,rows=15,new=FALSE)))
plot_qq(model,new=FALSE)
par(mfrow=c(1,1),cex.main=.80)
textplot(frequencies)
if(!(tolower(method) %in% c("pdlm","singlelag","both")))
	stop("Model for the exposure not implemeted yet")
# save and change warning options
ow <- options("warn")
options(warn=-1)
if(tolower(method)=="pdlm")
	{
	textplot("Effect estimates using distributed lag models")
	for (i in 1:length(pollutants))
		{
		risks <- capture.output(estimate_risks(model,pollutants[i],labels=labels[i],method="pdlm",unit=unit[i],digits=5,new=FALSE,...))
		textplot(risks)
		}
	}
else if(tolower(method)=="singlelag")
	{
	textplot("Effect estimates using single lag models")
	for (i in 1:length(pollutants))
		{
		risks <- capture.output(estimate_risks(model,pollutants[i],labels=labels[i],method="singlelag",unit=unit[i],digits=5,new=FALSE,...))
		textplot(risks)
		}
	}
else if(tolower(method)=="both")
	{
	textplot("Effect estimates using single lag models")
	for (i in 1:length(pollutants))
		{
		risks <- capture.output(estimate_risks(model,pollutants[i],labels=labels[i],method="singlelag",unit=unit[i],digits=5,new=FALSE,...))
		textplot(risks)
		}
	textplot("Effect estimates using distributed lag models")
	for (i in 1:length(pollutants))
		{
		risks <- capture.output(estimate_risks(model,pollutants[i],labels=labels[i],method="pdlm",unit=unit[i],digits=5,new=FALSE,...))
		textplot(risks)
		}
	}
options(ow) # reset warnings
while (!(is.null(dev.list())))
	dev.off()
}


unload <- function()
# unload ares and clear some things
{
if(exists(".aresEnv"))
	{
	rm(list=ls(all),envir=.aresEnv)
	gc(verbose=FALSE)  # memory clean up
	}
if (sum(search()=="package:ares2")>0)
	detach(package:ares2) 
else
	stop("The ares2 library is not attached")
}


# buildts <- function(date,x,unit="day")
# # build a times series aggregated by unit
# {
# names <- c(deparse(substitute(date)),deparse(substitute(x)))
# from <- min(date)
# to <- max(date)
# ref <- as.data.frame("date"=seq(from,to))
# data <- cbind.data.frame("date"=date,"x"=x)
# matched <- merge(data,ref,by="date",all=TRUE)
# out <- cbind.data.frame(matched$date,matched$x.x)
# colnames(out) <- names
# return(out)
# }
# 


model.family <- function(overdispersion)
# set the model family based on the option overdispersion
{
if(overdispersion)
	retval <- quasipoisson(link=log)
else
	retval <- poisson(link=log)
}


perc.range <- function(vars,min.perc=10,max.perc=90,digits=1)
# compute interpercentilic range given min and max
{
min <- sapply(vars,function(x)quantile(eval(parse(text=x)),prob=min.perc/100,na.rm=TRUE))
max <- sapply(vars,function(x)quantile(eval(parse(text=x)),prob=max.perc/100,na.rm=TRUE))
retval <- round(max-min,digits)
return(retval)
}


set_graph_window <- function(new_window)
# set graphic device
{
if (new_window & Sys.getenv("RSTUDIO")!="1") # workaround to fix buggy device in  RStudio
	getOption("device")()
}
