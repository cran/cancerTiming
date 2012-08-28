#for now, history is given as matrix 
#(1/ncopies, ..., 1); Must be length equal to the number of copies of chromosome for the history
#later will have function that takes encoding of history and creates this vector
#if exactAllele=TRUE, then assume that x/m is exactly the true allele frequency P -- not implemented yet

###Should code some EM steps to get good init point...
eventTiming<-function(x, m, history, totalCopy, method=c("fullMLE","partialMLE","Bayes"),type=c("gain","CNLOH"),
seqError=0, bootstrapCI=NULL, B=if(method=="Bayes") 10000 else 500,CILevel=0.95, alpha=1,tdf=4,normCont=0,verbose=TRUE,
returnAssignments=FALSE,coverageCutoff=1,minMutations=10,returnData=FALSE,init=NULL,maxiter=100, tol=0.0001){
	method<-match.arg(method)
	doCI<-!is.null(bootstrapCI)
	type<-match.arg(type)

	nMuts<-length(x)
	if(length(m)!=nMuts) stop("x and m must be of same length")

	A<-history #this will be a function later
	nCopies<-totalCopy
	if(is.null(dim(A))) stop("'history' should be a matrix of size (nEvents +1) x (nEvents +1)")
	if(ncol(A)!=nrow(A)) stop("'history' should be a square matrix (hint: do not include the allele frequency 1 except for CNLOH events)")
	K<-ncol(A)-1
	if(length(normCont)!=1 || normCont>1 || normCont<0) stop("normCont must be given and should be single value between 0 and 1.")

	#get final allele frequencies
	possAlleles<-allAF(nCopies,normCont=normCont,type="mutation")[[1]] 
	if(type=="gain") possAlleles<-head(possAlleles,-1) #can't have the allele frequency 1
	if(nrow(A)!=length(possAlleles)) stop(length(possAlleles),"possible alleles for this scenario, but history only has",nrow(A),"rows.")
	possAlleles<-errorAF(possAlleles,seqError=seqError)

	if(coverageCutoff<1) stop("coverageCutoff must be at least 1")
	nBelowCutoff<-length(which(m<coverageCutoff))
	nZero<-length(which(m==0))
	if(any(m<coverageCutoff)){		
		if(verbose) warning(nBelowCutoff, "mutations present below cutoff",coverageCutoff,", will be ignored")
		if(verbose & nZero>0) warning(nZero, "mutations have no reads at all")
		x<-x[which(m>=coverageCutoff)]
		m<-m[which(m>=coverageCutoff)]
	}
	summaryTable<-cbind(nMuts,length(m),nZero,nBelowCutoff)
	colnames(summaryTable)<-c("Total Mutations","Qualify for fitting","Zero Read Coverage","Coverage Below Cutoff")
	summaryTable<-t(summaryTable)
	colnames(summaryTable)<-"Number"
	call<-list(alleleSet=possAlleles,
		history=history, 
		totalCopy=totalCopy,
		type=type,
		method=method,
		seqError=seqError,
		normCont=normCont,
		coverageCutoff=coverageCutoff,
		minMutations=minMutations,
		B=B,
		init=init,
		maxiter=maxiter, 
		tol=tol)		
	if(returnData) call<-c(call,list(inputData=data.frame(x=x,m=m)))
	
	failReason<-NULL
	success<-TRUE
	if(length(m)<minMutations){
		failReason<-paste(failReason,"not enough mutations above cutoff that satisfy coverage criteria",sep=",")
		if(verbose) warning("less than ",minMutations," mutated locations meet coverage criteria; no estimate of pi will be calculated")
		success<-FALSE
	}
	# if(all(x==m)){
	# 	failReason<-paste(failReason,"all x equal to m",sep=",")
	# 	success<-FALSE
	# }
	dummyOutput<-list("pi"=rep(NA,length=ncol(A)),alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=NA,pathOfQ=NA),perLocationProb=NA,q=NA)
	names(dummyOutput$pi)<-paste("Stage",0:K,sep="")
	ci<-matrix(rep(NA,length=2*length(dummyOutput$pi)),ncol=2)
	row.names(ci)<-names(dummyOutput$pi)
	dummyOutput$piCI<-ci
	if(!is.null(failReason)){
		output<-dummyOutput
	}
	else{	
		rankA<-qr(A)$rank
		piId <- rankA==ncol(A)
		piEst<-rep(NA,ncol(A))
		if(piId){
			if(method%in%c("fullMLE","Bayes")){		
				### Estimate the q vector:
				output<-try(.estimateQ(x,m,possAlleles,history=history,init=init,maxiter=maxiter, tol=tol))
				if(output$optimDetails$nIter==maxiter){
					failReason<-paste(failReason,"hit maximum number of iterations",sep=",")
					success<-FALSE
				}
			}
			if(method=="partialMLE" ){
				mleAlleles<-mleAF(x,m,totalCopy=nCopies,seqError=seqError,normCont=normCont,maxCopy=if(type=="gain") totalCopy-1 else totalCopy)
				#q<-mleAlleles$alleleSet$Freq/sum(mleAlleles$alleleSet$Freq)
				output<-try(.estimateQ(x,m,possAlleles,alleleFreq=mleAlleles$alleleSet$Freq,history=history,init=init,maxiter=maxiter, tol=tol))				
#				output<-list(q=q,perLocationProb=mleAlleles$perLocationProb,alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=NA,pathOfQ=NA))
			}
			if(!inherits(output, "try-error")){	
				piEst<-solve(A)%*%output$q
				piEst<-as.vector(piEst)
				piEst<-piEst/sum(piEst)
				if(method=="Bayes"){
					initBayes<-piEst
					#initBayes<-seq(1,length(possAlleles))/length(possAlleles)
					bayesout<-.bayesEstimate(x,m,alleleSet=possAlleles,history=A,alpha=alpha,init=initBayes,tdf=tdf,nsim=B,CILevel=CILevel)
					piEst<-bayesout$pi
					q<-A%*%piEst
					q<-q/sum(q)
					output<-list(piCI=bayesout$piCI,q=q,perLocationProb=NA,alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=rbind(initBayes,bayesout$mode),pathOfQ=NA))					
				}
			}
			else{
				if(method=="Bayes"){
					initBayes<-seq(1,length(possAlleles))/length(possAlleles)
					bayesout<-.bayesEstimate(x,m,alleleSet=possAlleles,history=A,alpha=alpha,init=initBayes,tdf=tdf,nsim=B)
					piEst<-bayesout$pi
					q<-A%*%piEst
					q<-q/sum(q)
					output<-list(piCI=bayesout$piCI,q=q,perLocationProb=NA,alleleSet=possAlleles,optimDetails=list(nIter=NA,initial=rbind(initBayes,bayesout$mode),pathOfQ=NA))					
					
				}
				else{
					output<-dummyOutput
					failReason<-paste(failReason,"error returned in estimating q",sep=",")
					success<-FALSE	
				}
				
			}
			if(method=="Bayes") row.names(output$optimDetails$initial)<-c("mle","postMode")
		}
		else{
			output<-dummyOutput
			failReason<-paste(failReason,paste("history matrix not invertible, rank=",rankA,sep=""),sep=",")
			success<-FALSE
			if(verbose) warning(paste("history matrix not invertible, rank=",rankA,sep=""))
		}
		names(piEst)<-paste("Stage",0:K,sep="")
		output$pi<-piEst
		if(doCI & method!="Bayes"){
			if(!any(is.na(output$pi))){
				bootsample<-bootstrapEventTiming(B=B,x=x,m=m,call=call,type=bootstrapCI,pi=output$pi)
				ci<-t(apply(bootsample,2,quantile,c((1-CILevel)/2,1-(1-CILevel)/2)))
			}
			else {
				ci<-matrix(rep(NA,length=2*length(output$pi)),ncol=2)
				row.names(ci)<-names(output$pi)
			}
			#colnames(ci)<-c("lower","upper")
			output$piCI<-ci
			output<-output[c("pi","piCI","q","perLocationProb","optimDetails")]
		}
		else{
			if(method!="Bayes") output<-output[c("pi","q","perLocationProb","optimDetails")]
			else{
				if(!any(is.na(output$pi))){ row.names(output$piCI)<-names(output$pi)}
				output<-output[c("pi","piCI","q","perLocationProb","optimDetails")]
			}
		}
	}
	if(!returnAssignments) output<-output[-grep("perLocationProb",names(output))]
	return(c(output,list(summaryTable=summaryTable,success=success,failReason=if(!is.null(failReason)) gsub("^,","",failReason) else failReason,call=call)))
}

#simple EM algorithm to calculate \hat{q}, the proportions of the multinomial; can do this regardless of whether identifiable
#also gives probabilities of assignments to different alleles.
.estimateQ<-function(x,m,alleleSet,alleleFreq=NULL,history,init=NULL,maxiter=100, tol=0.0001){
	N<-length(x)
	possAlleles<-alleleSet
	lengthOfTheta<-length(possAlleles)-1
	if(nrow(history)!=length(possAlleles)) stop("length of alleleSet must equal the number of rows of history.")
	Ainv<-solve(history)
	if(is.null(init)){
		init<-history%*%rep(1/length(possAlleles),length=length(possAlleles))	
		init<-as.vector(init)/sum(init)
	} 
	##Break it into 2 because useful to have a function that assigns most likely allele to each
	EStep<-function(q){ ####Calculate P(Pi=a | Xi,q) each i and a
		PPerAllele<-mapply(possAlleles,q,FUN=function(a,qa){
			dbinom(x=x, size=m, prob=a)*qa 
		}) #returns a N x numb.Alleles Matrix; needs to be normalized
		PPerAllele<- prop.table(PPerAllele,1) #divide by the sum of each row
		return(PPerAllele)
	}
	
	###Constraints for the q
	#The feasible region is defined by ui %*% theta - ci >= 0. 
	uiL1<-diag(rep(-1,lengthOfTheta),nrow=lengthOfTheta,ncol=lengthOfTheta)
	ciL1<-rep(1,lengthOfTheta)
	uiGr0<-diag(rep(1,lengthOfTheta),nrow=lengthOfTheta,ncol=lengthOfTheta)
	ciGr0<-rep(0,lengthOfTheta)
	uiA<-Ainv%*%(rbind(diag(lengthOfTheta),rep(-1,lengthOfTheta)))
	ciA<- - Ainv%*%c(rep(0,lengthOfTheta),1)
	ui<-rbind(uiL1,uiGr0,uiA)
	ci<- c(-ciL1,-ciGr0,ciA)
	
	MStep<-function(PMat=NULL,qinit){ #estimate of q constrained so that solve(history)%*%q is all positive
		if(is.null(alleleFreq)) Na<-colSums(PMat) #from E-Step, the estimated values from the multinomial...basically just the maximum likelihood estimation if knew Pi exactly
		else Na<-alleleFreq
		f<-function(q){
			qS<-1-sum(q)
			#-dmultinom(Na, prob=c(q,qS), log = TRUE) #just to check; should be equivalent, but faster to not calculated normalizing constant
			- sum(Na*log(c(q,qS)))
		}

		gr<-function(q){
			qS<-1-sum(q)
			NS<-tail(Na,1)
			- (head(Na,-1)/q - NS/qS*q)
		}
		theta<-head(qinit,-1)
		if(lengthOfTheta==1){
			upper<-min((ci/ui)[ui<0])
			lower<-max((ci/ui)[ui>0])
			out<-optimize(f=f,interval=c(lower,upper))
			qout<-c(out$minimum,1-out$minimum)
		}
		else{
			out<-constrOptim(theta=theta,f = f, grad=NULL, ui=ui, ci=ci) #when use gr, get wrong answer -- returns me to init, doesn't iterate, etc. WHY????
			qout<-c(out$par,1-sum(out$par))
			}
		names(qout)<-names(possAlleles)
		return(qout)
	}
	TotalStep<-function(q){
	 	MStep(EStep(q),qinit=q)
	}
	if(is.null(alleleFreq)){
		qOld <- init;
		qNew<-TotalStep(qOld)
		allQ<-cbind(qOld,qNew)
		nIter<-1
		#convMessages<-c(firstOut$convergence)
		while(sum(abs(qNew - qOld)) >= tol & nIter <= maxiter){
			qOld<-qNew
			pMat<-EStep(qOld)
			qNew<-MStep(pMat,qinit=qOld)
			#qNew<-c(mstepOut$par,1-sum(mstepOut$par))
			#convMessages<-c(convMessages,mstepOut$convergence)
			nIter <- nIter + 1
			allQ<-cbind(allQ,qNew)
			#if(nIter %% 20 ==0) print(sprintf("finished iteration %s", nIter))
		}
		pMat<-EStep(qNew)
	}
	else{
		qNew<-MStep(qinit=init)
		pMat<-NA
		nIter<-NA
		allQ<-NA
	}
	names(qNew)<-names(possAlleles)
	return(list(q=qNew,perLocationProb=pMat,alleleSet=possAlleles,optimDetails=list(nIter=nIter,initial=init,pathOfQ=allQ)))
}
#run into problems when get close to the boundary -- conversion between theta and pi is tenuous
.bayesEstimate<-function(x,m,alleleSet,history,alpha=1,method=c("sir","imp","rej"),init,tdf=4,nsim=10000,CILevel=0.95){
	method<-match.arg(method)
	if(method!="sir") stop("Other methods than 'sir' are not yet operational")
	require(LearnBayes)
	pOfPi<-function(theta){
		(alpha-1)*(sum(theta-log(1+sum(exp(theta)))))
	}
	pOfX<-function(theta){
		thetaNum<-c(exp(theta),1)
		thetaDenom<-1+sum(exp(theta))
		Api<-history%*%thetaNum
		PPerX<-mapply(alleleSet,Api,FUN=function(a,qa){
			dbinom(x=x, size=m, prob=a)*qa 
		}) #returns a N x numb.Alleles Matrix
		PPerX<-rowSums(PPerX)
		cpi<-(1:length(alleleSet) +1) %*% thetaNum
		return(sum(log(PPerX))-length(x)*log(cpi))
	}
	
	
	#change of variables so that now theta[j]=log(pi[j-1]/pi[K])
	#assume alpha a single scalar
	chgOfVar<-function(theta){
		sum( theta+log(1+sum(exp(theta))-exp(theta))-2*log(1+sum(exp(theta))) )
	}
	logposterior<-function(theta,data){
		pOfPi(theta) + chgOfVar(theta) + pOfX(theta)
	}	
	pi2theta<-function(piV){
		head(log(piV/tail(piV,1)),-1)
	}
	theta2pi<-function(theta){
		piV<-exp(theta)/(1+sum(exp(theta)))
		return(c(piV,1-sum(piV)))
	}
	
	#use code in LearnBayes
	fit<-try(laplace(logposterior,pi2theta(init),x))
	if(inherits(fit, "try-error")){ #try with bad initial value that is in the interior
		fit<-try(laplace(logposterior,pi2theta(rep(1,length(init))/length(init)),x))		
	}
	tpar<-list(m=fit$mode,var=2*fit$var,df=tdf)
	if(method=="sir"){
		theta.s<-sir(logposterior,tpar,nsim,x)		
	}
	if(method=="imp"){ #not operational
		
	}
	if(method=="rej"){ #not operational
		Tdiff<-function(theta,datapar){
			data<-datapar$data
			fit<-datapar$fit
			d<-logposterior(theta,x)-dmt(theta,mean=fit$mode,S=2*fit$var,df=tdf,log=TRUE)
			return(d)
		}
		maxD<-laplace(Tdiff,fit$mode,list(data=x,fit=fit))
		d<-Tdiff(maxD$mode,list(data=x,fit=fit))
		theta.s<-rejectsampling(logposterior,tpar,d,nsim,x)		
	}
	if(is.null(dim(theta.s))){
		pis<-sapply(theta.s,theta2pi)
	}
	else{
		pis<-do.call(cbind,lapply(1:nrow(theta.s),function(i){theta2pi(theta.s[i,])}))
	} 
	out<-t(apply(pis,1,quantile,prob=c((1-CILevel)/2, 0.5,1-(1-CILevel)/2)))
	
	return(list(pi=out[,2],piCI=out[,c(1,3)], modeFit=theta2pi(fit$mode), nsim=nsim))
}

