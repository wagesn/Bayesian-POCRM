
###install required R packages
library(nnet)
library(dfcrm)

###Load the function 'bpocrm' 
bpocrm<-function(p0,p.skel,ttr,cohortsize,ncohort,n.stop,start.comb,cs){
 
# if a single ordering is inputed as a vector, convert it to a matrix
	if(is.vector(p.skel)) p.skel=t(as.matrix(p.skel));
	
	nord.tox = nrow(p.skel);
	mprior.tox = rep(1/nord.tox, nord.tox);  # prior for each toxicity ordering

asd<-1.34
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
    y=rep(0,ncomb);  #number of toxicities at each dose level
    n=rep(0,ncomb);  #number of treated patients at each dose level
    comb.curr = start.comb;  #current dose level	 
    ptox.hat = numeric(ncomb); #estimate of toxicity prob
    comb.select=rep(0,ncomb); #a vector of indicators for dose selection
    stop=0; #indicate if trial stops early
    i=1	
while(i <= ncohort)
    {
	# generate data for a new cohort of patients
		y[comb.curr] = y[comb.curr] + rbinom(1,cohortsize,p0[comb.curr]);
		n[comb.curr] = n[comb.curr] + cohortsize;

		if(any(n>n.stop)){
			stop<-0
			break
		}

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
			break
			}	
	
		distance=abs(ptox.hat-ttr)
		comb.curr<-which.is.max(-distance)
		
		i=i+1
		}
	if(stop==0){
		comb.select[comb.curr]=comb.select[comb.curr]+1;
		}
	return(list(comb.select=comb.select,tox.data=y,pt.allocation=n,stop=stop))
}
##########'bpocrm' end here

###Load the function 'bpocrm.sim' 
bpocrm.sim<-function(p0,p.skel,ttr,cohortsize,ncohort,n.stop,ntrial,cs){
	ncomb=length(p0)
	
	comb.select<-y<-n<-matrix(nrow=ntrial,ncol=ncomb)
	nstop=0
	
	for(i in 1:ntrial){
		result<-bpocrm(p0,p.skel,ttr,cohortsize,ncohort,n.stop,start.comb,cs)
		comb.select[i,]=result$comb.select
		y[i,]=result$tox.data
		n[i,]=result$pt.allocation
		nstop=nstop+result$stop
	}
	cat("True tox probability: ", round(p0,3), sep="\t",  "\n");
	cat("selection percentage: ", formatC(colMeans(comb.select)*100, digits=1, format="f"), sep="\t",  "\n");
	cat("number of toxicities: ", formatC(colMeans(y), digits=1, format="f"), sep="\t",   "\n");
	cat("number of patients treated: ", formatC(colMeans(n), digits=1, format="f"), sep="\t",   "\n");
	cat("percentage of stop: ", nstop/ntrial*100, "\n");
}
##########'bpocrm.sim' end here


#####Specify the total number of combinations
d<-6

#####Specify the number of possible toxicity orderings
s<-6   

###Specifiy the possible toxicity orderings of the drug combinations
orders<-matrix(nrow=s,ncol=d)
orders[1,]<-c(1,2,3,4,5,6)  
orders[2,]<-c(1,2,3,5,4,6)  
orders[3,]<-c(1,2,3,5,6,4)  
orders[4,]<-c(1,2,5,3,4,6)  
orders[5,]<-c(1,2,5,3,6,4)  
orders[6,]<-c(1,2,5,6,3,4)

###Specify a set of toxicity skeleton values
skeleton<-c(0.07,0.18,0.33,0.49,0.63,0.75)

###The following will adjust the location of each skeleton value to correspond 
###to the 's' possible orderings specified above. 
###'p.skel' is the matrix of toxicity skeleton values
p.skel<-matrix(0,nrow=s,ncol=d)
for(j in 1:s){
	p.skel[j,]<-skeleton[order(orders[j,])]
}
#print(p.skel)

###True toxicity probability scenarios
p1<-c(0.26,0.33,0.51,0.62,0.78,0.86)
p2<-c(0.12,0.21,0.34,0.50,0.66,0.79)
p3<-c(0.04,0.07,0.20,0.33,0.55,0.70)
p4<-c(0.01,0.04,0.05,0.17,0.33,0.67)
p5<-c(0.01,0.02,0.05,0.15,0.20,0.33)
p6<-c(0.50,0.60,0.70,0.80,0.90,0.95)


####################################
#
#	Input
#
#	p0 = true toxicity probabilities 
#	p.skel = toxicity skeleton values
#	ttr = target toxicity rate
# cohortsize = size of each cohort inclusion
# nchort = number of cohorts in a sinlge trial
#	start.comb = starting combination
#	ntrial = number of simulated trials
#	cs = confidence level used in stopping rule
# 
####################################

ttr=0.33       ##target toxicity rate 
cohortsize=1   ##cohort size for each inclusion
ncohort=22     ##number of cohorts
start.comb=3   ##starting combination
n.stop=100     ##Number of patients needed on one combination to stop the trial
ntrial=100     ##number of simulated trials 
cs=0.90        ##confidence level for the confidence interval 

##simulate a single trial
#bpocrm(p0,p.skel,ttr,cohortsize,ncohort,n.stop,start.comb,cs)

p0<-p2 ##true toxicity probability scenario
set.seed(580)  ##random seed
bpocrm.sim(p0,p.skel,ttr,cohortsize,ncohort,n.stop,ntrial,cs)

