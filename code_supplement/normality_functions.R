#-------------------------------------------------------------------------------
# Common functions useful for any simulation study with HLM residuals.
# 
# March 2013
# Adam Loy
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# The "usual tests"
#-------------------------------------------------------------------------------

## testing normality for a list of residuals
normality.tests <- function(.list) {
	pvals <- lapply(.list, 
         function(r) {
			ad  <- try(ad.test(r))
			cvm <- try(cvm.test(r))
			ks  <- try(lillie.test(r))
           
			pv <- data.frame(AD = ifelse(class(ad) == "try-error", NA, ad$p.value), 
							 CVM = ifelse(class(cvm) == "try-error", NA, cvm$p.value), 
							 KS = ifelse(class(ks) == "try-error", NA, ks$p.value))
			return(pv)  
		})
	
	pvals <- do.call("rbind", pvals)	 
}

## testing normality for a list of residuals
normality.test.statistics <- function(.list, column = NULL) { 
	tstats <- ldply(.list, 
         function(r, column) {
         	if(is.null(column)) {
         		ad  <- try(ad.test.stat(r))
				cvm <- try(cvm.test.stat(r))
				ks  <- try(lillie.test.stat(r))
         	} else {
         		ad  <- try(ad.test.stat(r[,column]))
				cvm <- try(cvm.test.stat(r[,column]))
				ks  <- try(lillie.test.stat(r[,column]))

         	}
           
			ts <- data.frame(AD = ifelse(class(ad) == "try-error", NA, ad), 
							 CVM = ifelse(class(cvm) == "try-error", NA, cvm), 
							 KS = ifelse(class(ks) == "try-error", NA, ks))
			return(ts)  
		}, column = column)
	
	return(tstats)
}


## Summarizing tests individually
summarize.tests <- function(pvals, alpha) {
  apply(pvals, 2, 
        FUN = function(pvs) {sum(pvs <= alpha, na.rm = TRUE) / sum(!is.na(pvs))})
}

## Calcualte a Monte Carlo p-value
mc.pvalue <- function(stat, dsn) {
	dsn <- dsn[complete.cases(dsn)]
	pvalue <- sum(dsn > stat) / length(dsn)
	return(pvalue)
}

## Testing the normality via Monte Carlo
monte.carlo.results <- function(sims, settings, column = NULL, null.dsn) {
    test.stats <- lapply(sims, normality.test.statistics, column = column)

    p.values   <- lapply(test.stats, 
    	function(x, null.dsn){
    		ad   <- sapply(x[, "AD"], mc.pvalue, dsn = null.dsn[, "AD"])
    		cvm  <- sapply(x[, "CVM"], mc.pvalue, dsn = null.dsn[, "CVM"])
    		ks   <- sapply(x[, "KS"], mc.pvalue, dsn = null.dsn[, "KS"])
    		
    		return(data.frame(AD = ad, CVM = cvm, KS = ks))
    		
    	}, null.dsn = null.dsn)
    	
    a10 <- ldply(p.values, summarize.tests, alpha = .1)
    a10 <- data.frame(settings, alpha = rep(.10, nrow(a10)), a10)
    
    a05 <- ldply(p.values, summarize.tests, alpha = .05)
    a05 <- data.frame(settings, alpha = rep(.05, nrow(a05)), a05)
    
    RVAL <- rbind(a10, a05)
	RVAL <- arrange(RVAL, error, ranef, alpha)
	
	return(RVAL)
    
}

## Testing normality for the results of the simulations
test.simulation.results <- function(sims, settings, column = NULL) {

    if(!is.null(column)) {
    	sims <- lapply(sims, FUN = function(x) lapply(x, FUN = function(x) x[,column]))
    }
    
	pvals <- lapply(sims, normality.tests)
	
	a10 <- lapply(pvals, summarize.tests, alpha = .10)
	a10 <- round(do.call("rbind", a10), 3)
	a10 <- data.frame(settings, alpha = rep(.10, nrow(a10)), a10)
	
	a05 <- lapply(pvals, summarize.tests, alpha = .05)
	a05 <- round(do.call("rbind", a05), 3)
	a05 <- data.frame(settings, alpha = rep(.05, nrow(a05)), a05)
	
	RVAL <- rbind(a10, a05)
	RVAL <- arrange(RVAL, error, ranef, alpha)
	
	return(RVAL)
}

print.format <- function(x, ranef = FALSE) {
	if(ranef) x <- x[,c(2, 1, seq(ncol(x))[-c(1,2)])]
	x$error <- ordered(x$error, levels = levels(x$error)[c(2,3,1)])
	x$ranef <- ordered(x$ranef, levels = levels(x$ranef)[c(2,3,1)])

	if(ranef) RVAL <- arrange(x, ranef, error, alpha)
	else RVAL <- arrange(x, error, ranef, alpha)
	
	return(RVAL)
}


## Testing the normality of the Lange and Ryan residuals
test.statistic.langeryan <- function(.model, column) {
	
	r <- lev2.langeryan.resid(.model)[,column]
	w <- lev2.marginal.var(.model)[,column]
	
	stat <- try(lillie.test.stat(x = r, weights = w))
	stat <- ifelse(class(stat) == "try-error", NA, stat)
	
	return(stat)
}

## Carrying out the monte carlo test for the lange and ryan residuals
langeryan.test <- function(test.stats, settings, null.dsn) {
	p.values <- lapply(test.stats, FUN = function(x, null.dsn) {
		sapply(x, mc.pvalue, dsn = null.dsn)
	},  null.dsn = null.dsn)
	
	a10 <- sapply(p.values, FUN = function(pvs, alpha) sum(pvs <= alpha, na.rm = TRUE) / sum(!is.na(pvs)), alpha = .1)
    a10 <- data.frame(settings, alpha = rep(.10, length(a10)), KS = a10)
    
    a05 <- sapply(p.values, FUN = function(pvs, alpha) sum(pvs <= alpha, na.rm = TRUE) / sum(!is.na(pvs)), alpha = .05)
    a05 <- data.frame(settings, alpha = rep(.05, length(a05)), KS = a05)
    
    RVAL <- rbind(a10, a05)
	RVAL <- arrange(RVAL, error, ranef, alpha)
	
	return(RVAL)

}

#-------------------------------------------------------------------------------
# Using a Monte Carlo distribution for p-value
#-------------------------------------------------------------------------------

# weighted empirical CDF
wecdf <- function(x, weights) {
    stopifnot(length(x) == length(weights))
    sw <- sum(weights)
    if (length(x) < 1) 
        stop("'x' must have 1 or more non-missing values")
    stopifnot(all(weights >= 0))
    ox <- order(x)
    x  <- x[ox]
    w  <- weights[ox]
    vals <- sort(unique(x))
    xmatch <- factor(match(x, vals), levels = seq_along(vals))
    wmatch <- tapply(w, xmatch, sum)
    wmatch[is.na(wmatch)] <- 0
    rval <- approxfun(vals, cumsum(wmatch) / sw, method = "constant", 
        yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    return(rval)
}   

# Lilliefor's (K-S) test statistic
lillie.test.stat <- function(x, weights = NULL) {
	if(is.null(weights)) {
		x <- sort(x[complete.cases(x)])
    	n <- length(x)
    	if (n < 5) 
        	stop("sample size must be greater than 4")
    	p <- pnorm((x - mean(x))/sd(x))
    	Dplus <- max(seq(1:n)/n - p)
    	Dminus <- max(p - (seq(1:n) - 1)/n)
    	K <- max(Dplus, Dminus)
	} else{
		x <- x[complete.cases(x)]
    	w <- weights[complete.cases(x)]
    	ox <- order(x)
    	x <- x[ox]
    	w <- w[ox]
    	n <- length(x)
    	if (n < 5) 
        	stop("sample size must be greater than 4")
    	edf <- wecdf(x, w)
    	ep <- edf(x)
    	p <- pnorm((x - mean(x))/sd(x))  
    	Dplus <- max(ep - p)
    	Dminus <- max(p - ep)
    	K <- max(Dplus, Dminus)
	}
    
    return(K)
}

# Cramer-von Mises test statistic
cvm.test.stat <- function(x, weights = NULL) {
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if (n < 8) 
        stop("sample size must be greater than 7")
    p <- pnorm((x - mean(x))/sd(x))
    W <- (1/(12 * n) + sum((p - (2 * seq(1:n) - 1)/(2 * n))^2))
    
    return(W)
}

# Anderson-Darling test statistic
ad.test.stat <- function(x, weights = NULL) {
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if (n < 8) 
        stop("sample size must be greater than 7")
    p <- pnorm((x - mean(x))/sd(x))
    h <- (2 * seq(1:n) - 1) * (log(p) + log(1 - rev(p)))
    A <- -n - mean(h)
	
	return(A)
}
