library(tidyverse)
library(braidrm)

estimateAIC <- function(act,fit,npar) {
	df <- length(act)
	resid <- fit-act
	ssErr <- sum(resid^2)
	lLik <- -(df/2) * (log(2*pi)+log(ssErr)-log(df)+1)
	k <- npar
	aic <- 2*k - 2*lLik + 2*k*(k+1)/(df-k-1)
	return(aic)
}
estimateBIC <- function(act,fit,npar) {
	df <- length(act)
	resid <- fit-act
	ssErr <- sum(resid^2)
	lLik <- -(df/2) * (log(2*pi)+log(ssErr)-log(df)+1)
	k <- npar
	bic <- log(df)*k - 2*lLik 
	return(bic)
}

estimateParameters <- function(data,kprior) {
	result <- as.list(data %>% select(RowSid,RowName,ColSid,ColName) %>% slice(1))
	
	data$bliss <- blissDeviation(cbind(data$ColConc,data$RowConc),data$Value,range=c(100,0),clip="pre")
	posrel <- data$RowConc>0 & data$ColConc>0
	result$VBliss <- mean(data$bliss[posrel])
	result$MinBliss <- min(data$bliss)
	result$MaxBliss <- max(data$bliss)
	
	cdata <- data %>% select(RowConc,ColConc,Value) %>%
		mutate(Value=pmin(pmax(Value,0),100))
	start <- c(exp(mean(log(cdata$ColConc[cdata$ColConc>0]))),
				exp(mean(log(cdata$RowConc[cdata$RowConc>0]))),
				1,1,0,100,0,0,0)
	llim <- c(min(cdata$ColConc[cdata$ColConc>0])/10,
			   min(cdata$RowConc[cdata$RowConc>0])/10,
			   0.25,0.25,-1.96,0,0,0)
	ulim <- c(max(cdata$ColConc[cdata$ColConc>0])*10,
			   max(cdata$RowConc[cdata$RowConc>0])*10,
			   4,4,100,200,200,200)
	
	bfit <- braidrm(Value~ColConc+RowConc,cdata,model="kappa2",start=start,
					direction=-1,lower=llim,upper=ulim,prior=kprior,getCIs=FALSE)
	bpar <- coef(bfit)
	result <- c(result,as.list(bpar))
	
	cdata$fit <- fitted(bfit)
	if (var(cdata$fit)>0) { result$r2 <- cor(cdata$Value,cdata$fit)^2 }
	else { result$r2 <- 0 }
	if (bpar[[9]]>bpar[[6]]) {
		iaes <- as.list(c(1,1))
	} else {
		bpar[[9]] <- min(bpar[7:9])
		iaes <- as.list(estimateIAE(bpar,
									levels=c(50,10),
									limits=c(max(cdata$ColConc),max(cdata$RowConc))))
	}
	names(iaes) <- c("iae50","iae90")
	result <- c(result,iaes)
	
	result$AIC <- estimateAIC(cdata$Value,cdata$fit,8)
	result$BIC <- estimateBIC(cdata$Value,cdata$fit,8)
	
	return(result)
}

### ESTIMATE SURFACE PARAMETERS FRO FLUCONAZOLE SCREEN ###

tab <- read_csv("data/responses_flu.csv", show_col_types = FALSE)

measns <- tab %>% filter(RowConc==0 & ColConc==0) %>% pull(Value) %>% sd()
kprior <- kappaPrior(spread=measns)

blocks <- tab %>% 
	group_by(RowSid,RowName,ColSid,ColName) %>%
	group_split()

results <- map(blocks,~estimateParameters(.x,kprior),.progress=TRUE) %>% 
	bind_rows()

write_tsv(results,file="fits/surfaceParameters_flu.txt")

### ESTIMATE SURFACE PARAMETERS FRO CASPOFUNGIN SCREEN ###

tab <- read_csv("data/responses_cas.csv", show_col_types = FALSE)

measns <- tab %>% filter(RowConc==0 & ColConc==0) %>% pull(Value) %>% sd()
kprior <- kappaPrior(spread=measns)

blocks <- tab %>% 
	group_by(RowSid,RowName,ColSid,ColName) %>%
	group_split()

results <- map(blocks,~estimateParameters(.x,kprior),.progress=TRUE) %>% 
	bind_rows()

write_tsv(results,file="fits/surfaceParameters_cas.txt")