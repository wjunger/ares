# distributed lag models
lagvar <- function(var,k)
# generate lagged var. lag 0 is included
{
n <- length(var)
new.vars <- var   # matrix to hold lagged vars with zero lagged

for (i in 1:k)
    {
    lagged <- c(rep(NA,i),var[1:(n-i)])
    new.vars <- cbind(new.vars,lagged)
    }

# naming new variables
dimnames(new.vars)[[2]] <- paste(deparse(substitute(var)),".L",0:k,sep="")

return(new.vars)
}


pdl <- function(var,lags,degrees)
# generate a distributed lags basis for var
{
if (degrees < 1)
    stop("No polynomial structure specified.")

lagged <- lagvar(var,lags)  # put lagged variable into a matrix
# lagged <- na.omit(lagged)

new.vars <- apply(lagged,1,"sum")     # W_0=Z_0+Z_1+...+Z_q
for (i in 1:degrees)
    {
    var.temp <- 0
    for (j in 1:(lags))
        {
        var.temp <- var.temp + j^i*lagged[,j]
        }
    new.vars <- cbind(new.vars,var.temp)
    }

# naming new variables
dimnames(new.vars)[[2]] <- paste(deparse(substitute(var)),".P",0:degrees,sep="")

return(new.vars)
}

get.beta <- function(coeff,var.coeff,lags,degrees,prefix="beta")
# get beta and standard errors
{
# parameters
beta <- NULL    # room for beta
#npar <- length(diag(var.coeff))   # parameters in the model
#fetapar <- npar-degrees    # first eta parameter
eta <- coeff[substring(names(coeff),1,4)=="pdl("]  # get the coefficients of the pdl basis
#eta <- coeff[(length(coeff)-degrees):(length(coeff))]      # get eta parameters !!! must correct for missing in holidays!!!! some coefficients disappear in the covariance matrix 
names <- names(eta)

for (i in 0:(lags))
    {
    beta.temp <- 0
    for (j in 0:(degrees))
        {
        beta.temp <- beta.temp+eta[j+1]*i^j
        }
    beta <- c(beta,beta.temp)
    }
names(beta) <- paste(prefix,".L",0:lags,sep="")
# variances
var.beta <- NULL    # room for beta se
#var.eta <- diag(var.coeff[fetapar:npar,fetapar:npar])
#covar.eta <- var.coeff[fetapar:npar,fetapar:npar]
covar.eta <- var.coeff[names,names]
var.eta <- diag(covar.eta)
low.covar.eta <- lower.tri(covar.eta)*covar.eta  # lower.tri() returns TRUE/FALSE

for (i in 0:lags)
    {
    var.beta.temp <- 0
    var.beta.temp1 <- 0
    var.beta.temp2 <- 0
    for (j in 0:degrees)
        var.beta.temp1 <- var.beta.temp1+var.eta[j+1]*i^(2*j)
    for (j in 0:(degrees-1))
        for (s in 1:(j+1))
            var.beta.temp2 <- var.beta.temp2+low.covar.eta[j+2,s]*2*(i^j)*(i^s)
    var.beta.temp <- var.beta.temp1+var.beta.temp2
        
    var.beta <- c(var.beta,var.beta.temp)
    }
names(var.beta) <- paste(prefix,".L",0:lags,sep="")

# overall estimates
overall.beta <- sum(beta)
overall.var <- sum(covar.eta)

return(list(beta=beta,se=sqrt(var.beta),overall.beta=overall.beta,overall.se=sqrt(overall.var)))
}


pdlm <- function(model,var,lags=5,degrees=2,...)
# fit a polynomial distributed lags model 
{
called <- match.call()
# this is a bit messy -- check it later
if(!exists("ares.active.dataset",envir=.aresEnv))
	stop("Setup() has not been used yet")
if(exists("ares.selection",envir=.aresEnv))
	ares.selection <- get("ares.selection",envir=.aresEnv)
else
	stop("Setup() has not been used yet")
if(!is.character(var))
	stop("Argument x must be a literal")
#	var <- deparse(substitute(var))
#variate <- eval(parse(text=var),envir=list(model$data,.aresEnv,ares.active.dataset))
variate <- .aresEnv$ares.active.dataset[var]
# ------------------------
form <- as.formula(paste("~.+pdl(",var,", ",lags,", ",degrees,")",sep=""))
formula <- update(model$formula,form)
pdl.model <- update(model,formula,...)

coeff <- coef(pdl.model)
coeff.var <- vcov(pdl.model)
pars <- get.beta(coeff,coeff.var,lags,degrees)

retval <- list(cmodel=pdl.model,variate=variate,var.name=var,beta=pars,lags=lags,degrees=degrees,call=called)
class(retval) <- c("pdlm")
return(retval)
}


plot.pdlm <- function(x,unit=10,confidence.level=0.95,labels=NULL,new=TRUE,...)
# plot pdl coefficients or relative risks
{
if(missing(x))
	stop("Model object is missing")
if(!inherits(x,"pdlm"))
	stop("Model object is not of class pdlm")

beta <- x$beta$beta
se.beta <- x$beta$se
z <- abs(qnorm((1-confidence.level)/2))
# fixed unit
RRk <- exp(unit*beta)
LBRRk <- exp(unit*beta-unit*z*se.beta)
UBRRk <- exp(unit*beta+unit*z*se.beta)

if (is.null(labels))
	labels <- paste("Lag",0:(length(beta)-1),sep=" ")
else if (length(labels)==1)
	labels <- paste(labels,0:(length(beta)-1),sep=" ")

# plotting RR
nticks <- 20
if (new)
	getOption("device")()
stockplot(RRk,LBRRk,UBRRk,ref.line=1,xlabels=labels,ticks=nticks,new=new,main=paste("Relative risk for",round(unit,2),"variation of the pollutant"),sub=paste("Pollutant:",x$var.name),xlab="Exposure",ylab="Relative risk",...)
}


summary.pdlm <- function(object,...)
# get summary for pdl model
{
if(missing(object))
	stop("Model object is missing")
if(!inherits(object,"pdlm"))
	stop("Model object is not of class pdlm")

call <- object$call
beta <- object$beta$beta
se.beta <- object$beta$se
t <- beta/se.beta
pvalue <- 2*(1-pnorm(abs(t)))
coef.table <- cbind(beta,se.beta,t,pvalue)
dimnames(coef.table) <- list(names(beta),c("Estimate","Std. Error","t value","Pr(>|t|)"))
res <- resid(object$cmodel,type="deviance")
resid.stats <- quantile(res,probs=c(0.0,0.25,0.50,0.75,1.0),na.rm=TRUE)
names(resid.stats) <- c("Min","1Q","Median","3Q","Max")
deviance <- object$cmodel$deviance
null.deviance <- object$cmodel$null.deviance  
aic <- object$cmodel$aic  
iter <- object$cmodel$iter  
df.null <- object$cmodel$df.null
df.residual <- object$cmodel$df.residual  
na <- length(object$cmodel$na.action)

retval <- list(call=call,coef.table=coef.table,resid.stats=resid.stats,deviance=deviance,null.deviance=null.deviance,aic=aic,iter=iter,df.null=df.null,df.residual=df.residual,na=na)
class(retval) <- "summary.pdlm"
return(retval)
}


print.summary.pdlm <- function(x,digits=getOption("digits"),...)
# print method for summary
{
cat("\nCall:\n")
print(x$call)
cat("\nDeviance Residuals:\n")
print(round(x$resid.stat,digits))
cat("\nCoefficients:\n")
print(x$coef.table)
cat("\n    Null deviance: ",round(x$null.deviance,digits)," on ",round(x$df.null,digits)," degrees of freedom\n",sep="")
cat("\nResidual deviance: ",round(x$deviance,digits)," on ",round(x$df.residual,digits)," degrees of freedom\n",sep="")
cat("(",x$na," observations deleted due to missingness)\n",sep="")
cat("\nAIC: ",x$aic,"\n",sep="")
cat("\nNumber of Fisher Scoring iterations: ",x$iter,"\n",sep="")
}

print.pdlm <- function(x,digits=getOption("digits"),...)
# print method for pdlm
{
cat("\nCall:\n")
print(x$call)
cat("\n    Null deviance: ",round(x$cmodel$null.deviance,digits)," on ",round(x$cmodel$df.null,digits)," degrees of freedom",sep="")
cat("\nResidual deviance: ",round(x$cmodel$deviance,digits)," on ",round(x$cmodel$df.residual,digits)," degrees of freedom\n",sep="")
cat("(",length(x$cmodel$na.action)," observations deleted due to missingness)\n",sep="")
cat("\nAIC: ",x$cmodel$aic,"\n",sep="")
}
