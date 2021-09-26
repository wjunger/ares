# first and last

.aresEnv <- new.env(parent=emptyenv())
# workaround to attach thing
#.aresEnv$ares.active.dataset <- as.data.frame(NA)
#assign("ares.active.dataset",function(x)NA,envir=.aresEnv)
#assign("ares.selection",function(x)NA,envir=.aresEnv)
#assign("ares.weights",function(x)NA,envir=.aresEnv)

.onLoad <- function(libname,pkgname)
{
if(getRversion() >= "2.15.1") 
	utils::globalVariables(c(".aresEnv","ares.active.dataset","ares.selection","ares.weights","savefile","doe","ticks","devSVG"),package="ares")
#if(getRversion() >= "3.1.0") 
#	utils::suppressForeignCheck("localvariable")

}

.onAttach <- function(libname,pkgname)
# start up options
{
#ver <- utils::packageDescription(pkgname,fields="Version")
msg="Please, take a look at the help pages.\n"
packageStartupMessage(msg)
}

.onUnload <- function(libpath)
{
#library.dynam.unload("ares",libpath)
}


