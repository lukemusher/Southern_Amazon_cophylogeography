#  MOD.SUMMARY(model, white.adjust=c("none", "const", "HC3", "HC0", "HC1", "HC2", "HC4", "cluster"), cluster=NULL, std.coef=c("none", "beta", "semi", "robust"), bic=TRUE, odds=FALSE, CI=FALSE, CIlevel=0.95)

#note: originally based on a method written by John Fox
# white.adjust = specify the type of White's standard errors to use to adjust for clustering in the data (ex. for heteroscedasticity), the original is HC0, and the others are modifications 
	# specifying cluster gives cluster-robust standard errors in one-dimension - i.e. bases the correction on a variable that represents how the data might be clustered (ex. by nation, or by family, or repeated measurements on the same person)
# cluster = the variable that represents the groups that are creating the heteroscedasticity
# std.coef = the type of standardized regression coefficient to return - can be beta (fully standardized), semi-standardized, a mix of fully and semi-standardized where dichotomous predictors are semi-standardized, or robust (based on inter-quartile ranges, less affected by outliers)
# bic = logical; should the Bayesian Information Criterion be returned?
# odds = logical; should odds ratios be returned instead of unstandardized regression coefficients? - only works for logistic regressions run using glm()
# CI = logical; should confidence intervals be returned?  confidence intervals will display in the scale used for the coefficients (e.g. odds ratios if odds=TRUE) - confidence intervals will also adjust for heteroscedasticity - they do NOT calculate confidence intervals based on profiling (which is the way that confint.glm() does)
# CIlevel = if CI=TRUE, what confidence level should be used?  default is 0.95
#returnobj = logical, should the results be returned as a summary object?  if FALSE (the default) results are simply printed


#TO DO:  Fix confidence intervals - right now they don't do profiling for glms


mod.summary = function(model, white.adjust=c("none", "const", "HC3", "HC0", "HC1", "HC2", "HC4", "cluster"), cluster=NULL, std.coef=c("none", "beta", "semi", "mix", "robust"), bic=TRUE, odds=FALSE, CI=FALSE, CIlevel=0.95) {

#function to estimate the variance of the latent variable underlying a binary or ordinal logit outcome
estyvar=function(m) {
		#adapted from fitstat in Stata
		modtype=class(m)[1]
		if(modtype=="glm") modtype=family(m)$link
		modframe=apply(m$model,2, as.numeric)
		#are there any weights?
		if("(weights)" %in% colnames(modframe)){
			w = modframe[,"(weights)"]
			modframe=modframe[,-which(colnames(modframe)=="(weights)")]
		} else
		w = 1
		#this makes all data deviated from its mean
		modframe=apply(modframe,2,function(x) {x-mean(x)})
		new.n=dim(modframe)[1]-1
		#get just the X matrix, not the outcome y
		require(car)
		Vx = (wcrossprod(modframe, modframe, w)/new.n)[-1,]
		#adjust for single variable models where the matrix becomes a vector
		if(class(Vx)=="numeric") Vx = matrix(Vx, nrow=1)[,-1] else {
		Vx=Vx[,-1]
		}
		#explained variance, make sure there is no intercept
		ss.ex=switch(modtype, polr=coef(m) %*% Vx %*% coef(m), logit=coef(m)[-1] %*% Vx %*% coef(m)[-1])
		#error variance, assumed to be 3.29 in logit models
		sse = pi^2/3
		totalvar=ss.ex + sse
		return(totalvar)
}

#standardized, semi-standardized, and robust regression coefficients
coeffs=function(MOD, type)
	{
		#get the model matrix (this handles any categorical variable correctly)
    dat=as.data.frame(model.matrix(MOD$terms, data=MOD$model))
  
    b = switch(modtype, lm=summary(MOD)$coef[-1,1], logit=summary(MOD)$coef[-1,1], polr=coef(MOD))
		se = switch(modtype, lm=summary(MOD)$coef[-1,2], logit=summary(MOD)$coef[-1,2], polr=sqrt(diag(vcov(MOD)))[1:length(coef(MOD))])
		if(type!="semi")sx = sapply(dat[-1], sd)
		sy = switch(modtype, lm=sapply(MOD$model[1], sd), logit=sqrt(estyvar(MOD)), polr=sqrt(estyvar(MOD)))
		if(type=="beta") {beta = b * sx/sy; ses=se*(sx/sy)}
		if(type=="semi") {beta = b/sy; ses=se/sy}
    if(type=="mix") {
      get.uniques=sapply(dat[-1], unique, simplify=FALSE)
      lengths=sapply(get.uniques, length)
      max1=sapply(get.uniques, max)
      min0=sapply(get.uniques, min)
      is.binary=lengths==2 & max1==1 & min0==0
      beta = ifelse(is.binary, b/sy, b * sx/sy)
      ses = ifelse(is.binary, se/sy, se*(sx/sy))
    }
		if(type=="robust"){
			sx <- sapply(MOD$model[-1], IQR)
    		sy <- sapply(MOD$model[1], IQR)
    		beta <- b * sx/sy
    		ses=se*(sx/sy)
		}
		return(cbind(beta, ses))
	}

#Confidence intervals - I have to include this because the built in function confint() won't use corrected vcov matrices
confints=function(object, parm, level = 0.95, ...)
	{
    	cf <- object$coef[,1]
    	pnames <- names(cf)
    	if (missing(parm)) 
        	parm <- pnames
    	a <- (1 - level)/2
    	a <- c(a, 1 - a)
    	pct = paste(round(100 * a, 1), "%")
      	fac <- qt(a, object$df[2])
    	ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    	ses <- sqrt(diag(object$vcov))[parm]
    	ci[] <- cf[parm] + ses %o% fac
    	ci
	}

#Corrected vcov matrices based on clustering (rather than just assuming a generic approach like HC0 or HC3)
clx = function(fm, dfcw, cluster){
         # R-codes (www.r-project.org) for computing
         # clustered-standard errors. Mahmood Arai, Jan 26, 2008.
	 # The arguments of the function are:
         # fitted model, degree of freedom correction (not used in my implementation) and cluster) 
         # reweighting the var-cov matrix for the within model (i.e. the meat part of the sandwhich matrix)
         #first make sure you have the same observations in your cluster variable that were used in the fitted model
         if(length(cluster) != length(unlist(labels(fm$model)[1]))) cluster.mod = cluster[as.numeric(unlist(labels(fm$model)[1]))]
         else cluster.mod=cluster
         M <- length(unique(cluster.mod))   
         N <- length(cluster.mod)          
         K <- fm$rank                        
         dfc <- (M/(M-1))*((N-1)/(N-K))  #finite sample adjustment
         uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster.mod, sum));  #i.e. summing within each unique level of cluster.mod for each row in the matrix containing the empirical estimating functions of the fitted model
         #NOTE: the crossprod() is equivalent to sum(u*t(u)), i.e. the sum of the vector created by u*transpose(u)
         vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)*dfcw
	}

#Warnings and data checks	
	if(odds) 
		{
		if (class(model)[1] != "glm") stop("Fitted model must of class 'glm' in order to use odds ratios")
		if (model$family$link != "logit") stop("Cannot produce odds ratios, model is not a logistic regression")
		}
	if(!require(sandwich)) stop("Required sandwich package is missing.")
	if(!require(lmtest)) stop("Required lmtest package is missing.")
	if(match.arg(white.adjust)=="cluster") if(is.null(cluster)) stop("Must specify a variable representing hetero. clustering to use the 'cluster' option in white.adjust")
#drop weights from model matrix if there are any
if("(weights)" %in% names(model$model)) {
	model$model = model$model[-which(names(model$model)=="(weights)")]
}
#type of model
modtype=class(model)[1]
if(modtype=="glm") modtype=model$family$link

#Calculating vcov matrix for SEs
	adjust=match.arg(white.adjust)
	if(adjust=="none") V = vcov(model)
	if(adjust=="cluster") V = clx(model, 1, cluster)
	if((adjust != "none") & (adjust != "cluster")) V = vcovHC(model, type=adjust)
	
#Building the summary object
	sumry = summary(model)
	table = coeftest(model, vcov=V)
	sumry$vcov=V

#Update coefficients based on specified options	
	if(CI) {
#if you can just change the relevant features of a model object, you can use the built in functions confint() and not have to rewrite them - i.e. you need to change the vcov matrix that the confint() functions get
		table2 = confints(sumry, level=CIlevel)
		if(odds) table2=exp(table2)
		if (match.arg(std.coef)!="none")
		{
			cimod=model
			for(i in 1:2)
			{
				cimod$coefficients=table2[,i]
				table2[,i]=c(NA,coeffs(cimod, std.coef)[,1])
			}
		} 
		sumry$confint=table2
	} #end if statement for confident intervals

	if (match.arg(std.coef)!="none") {
		if(modtype %in% c("lm", "logit")) {
      table[,c(1,2)]= rbind(rep(NA, 2), coeffs(model, std.coef))
			table[1,] = rep(NA, 4)
			}
		if(modtype=="polr") table[,c(1,2)]= coeffs(model, std.coef)
		nms=colnames(table)
		nms[1]=switch(std.coef, beta="Std. coef", semi="Semi-std. coef", mix="Mixed semi/std. coef", robust="Robust coef")
		colnames(table)=nms
		}
	if(odds) {
		nms=colnames(table)
		nms=c(nms[1], "Odds ratio", nms[2:4])
		table=cbind(table[,1], exp(coef(model)), table[,2:4])		
		colnames(table)=nms
	}	
	sumry$coefficients = table
#additional information to attach to summary object
	sumry$bic = AIC(model, k=log(length(model$fitted.values)))
	#if(returnobj) return(sumry)
	#else print(sumry)

#information for printing
	cat("Information about displayed results", "\n")
	if(bic) {
		cat("BIC: ", sumry$bic, "\n")
	}
	if(odds) cat("Odds ratios used", "\n")
	if((adjust != "const") & (adjust != "none"))  cat("Note: Heteroscedasticity-consistent SEs using adjustment", adjust, "\n")
	if(match.arg(std.coef)=="beta") cat("Note: Standardized regression coefficients/SEs used", "\n")
	if(match.arg(std.coef)=="robust") cat("Note: Robust standardized regression coefficients/SEs used", "\n")
	if(match.arg(std.coef)=="semi") cat("Note: Semi-standardized regression coefficients/SEs used", "\n")
  if(match.arg(std.coef)=="mix") cat("Note: Standardized regression coefficients/SEs used for non-binary predictors", "\n", "Semi-standardized regression coefficients/SEs used for binary predictors", "\n")
	if(CI) {
		cat("\n")
		cat(CIlevel*100,"% confidence intervals", "\n")
		print(table2)
		}
			return(sumry)
	}
	