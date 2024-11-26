
###install required R packages
library(nnet)
library(dfcrm)

bpocrm.imp<-function(p.skel,ttr,y,n,cs){
 
# if a single ordering is inputed as a vector, convert it to a matrix
	if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel));
	
	nord.tox = nrow(p.skel);
	mprior.tox = rep(1/nord.tox, nord.tox);  # prior for each toxicity ordering

asd<-1
bcrmh<-function(a,p,y,n){
	s2=asd
	lik=exp(-0.5*a*a/s2)
	for(j in 1:length(p)){
		pj=p[j]**exp(a)
		lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
		}
	return(lik);
    }

bcrmht<-function(a,p,y,n){
	s2=asd
	lik=a*exp(-0.5*a*a/s2)
	for(j in 1:length(p)){
		pj=p[j]**exp(a)
		lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
		}
	return(lik);
    }

bcrmht2<-function(a,p,y,n){
	s2=asd
	lik=a^2*exp(-0.5*a*a/s2)
	for(j in 1:length(p)){
		pj=p[j]**exp(a)
		lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
		}
	return(lik);
    }

### run a trial 	
    ncomb = ncol(p.skel);   #number of combos
    ptox.hat = numeric(ncomb); # estimate of toxicity prob
    stop=0; #indicate if trial stops early
    
		marginal.tox = est.tox=e2.tox=rep(0, nord.tox);
		for(k in 1:nord.tox)
		{
			marginal.tox[k] = integrate(bcrmh,lower=-Inf,upper=Inf, p=p.skel[k,], y=y,n=n,abs.tol = 0)$value;
			est.tox[k]=integrate(bcrmht,lower=-10,upper=10, p.skel[k,], y, n,abs.tol = 0)$value/marginal.tox[k]
			e2.tox[k]=integrate(bcrmht2,lower=-10,upper=10, p.skel[k,], y, n,abs.tol = 0)$value/marginal.tox[k]
		}		
		postprob.tox = (marginal.tox*mprior.tox)/sum(marginal.tox*mprior.tox);
		# toxicity model selection, identify the model with the highest posterior prob
		if(nord.tox>1){ 
			mtox.sel = which.is.max(postprob.tox); 
		} else{
			mtox.sel = 1;
		}
		ptox.hat=p.skel[mtox.sel,]**exp(est.tox[mtox.sel])
		post.var.tox=e2.tox[mtox.sel]-(est.tox[mtox.sel])^2
		crit.tox=qnorm(0.5+cs/2)
		lb.tox=est.tox[mtox.sel]-crit.tox*sqrt(post.var.tox)
		ub.tox=est.tox[mtox.sel]+crit.tox*sqrt(post.var.tox)
		ptoxL=p.skel[mtox.sel,]**exp(ub.tox)
		ptoxU=p.skel[mtox.sel,]**exp(lb.tox)
		
		if(ptoxL[1]>ttr){
			stop=1
			}	
	
		distance=abs(ptox.hat-ttr)
		comb.curr<-which.is.max(-distance)
		
		return(list(comb.select=comb.curr,DLTs.observed=y,pts.treated=n,ptox.est=round(ptox.hat,3),postmodp=round(postprob.tox,3),ptox1.lower.bound=round(ptoxL[1],3),stop=stop))
}
##########'bpocrm.imp' end here

#####Specify the total number of combinations
d<-4

#####Specify the number of possible toxicity orderings
s<-2

###Specifiy the possible toxicity orderings of the drug combinations
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,3,4)  
orders[2,]<-c(1,2,4,3)  

###Specify a set of toxicity skeleton values
skeleton<-c(0.122,0.204,0.30,0.402)

###The following will adjust the location of each skeleton value to correspond 
###to the 's' possible orderings specified above. 
###'p.skel' is the matrix of toxicity skeleton values
p.skel<-matrix(0,nrow=s,ncol=d)
for(j in 1:s){
	p.skel[j,]<-skeleton[order(orders[j,])]
}
#print(p.skel)

ttr=0.30  ##target toxicity rate 
cs=0.99	##confidence level for the confidence interval

y=c(0,1,0,0)   ##number of DLT's in each combo
n=c(3,3,0,0)   ##number of patients treated on each combo
bpocrm.imp(p.skel,ttr,y,n,cs)

