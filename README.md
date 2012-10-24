# Version 20120925
# Changelog:
# Test with JIT enabled
# Added ILich's initial arrangement
# Number of iteration before temperature reduced is changed to min(n.var*n.var + 1,rate)
# Max rate = n instead of 1000
# Initial t0 = f.old/n/n.var where f.old is the inital difference 
#Test for 50/n,
src <- '
arma::mat xy = Rcpp::as<arma::mat> (xys);
return wrap(cor(xy));
'
 fun3 <- cxxfunction(signature(xys="numeric"),body=src,plugin="RcppArmadillo")
 
 cor <- cmpfun(fun3)

sai = function(xy,cor.matrix,method="pearson",c.rate=0.85,maxitr=1e6,logs =0, plotxy=0,chart=0)
{
   require(psych)
	 round = 2
	 #Find the deviations of all correlation coefficents with weight <- 1 - Cost function
	 #(round(cor(update.xy(xy,cor.matrix)),2))
	 	update.xy <- function (xy,cor.matrix)
	{
		#set.seed(1)
		f.reg.mat = function(cor.mat,mean_std,nvar)
		{
			a = matrix(0,nrow=1,ncol=nvar)
			reg.matrix = matrix(0,nrow=(nvar-1),ncol=nvar)
			for (cnvar in 2:(nvar)) #start of calculating reg.matrix
			{
				b = solve(cor.mat[1:(cnvar-1),1:(cnvar-1)],cor.mat[1:(cnvar-1),(cnvar)])
				a0 = mean_std[1,cnvar]
				for (i in 1:(cnvar-1)) 
				{
					a[i] = mean_std[2,cnvar]/mean_std[2,i]*b[i]
					a0 = a0 - a[i]*mean_std[1,i]
				}
				reg.matrix[(cnvar-1),1]=a0
				reg.matrix[(cnvar-1),2:(i+1)]=a[1,1:(cnvar-1)]
			}#end of calculating reg.matrix
			return(reg.matrix)
		}
		#------------------------------- Calculating standard errors coefficients ----------------------------------
		f.e.mat = function(cor.mat,mean_std,nvar)
		{
			e = matrix(0,nrow=1,ncol=(nvar-1))
			for (i in 1:(nvar-1)) 
			{
				e[i]=mean_std[2,i+1]*sqrt(1-cor.mat[1,i+1]^2)
			}
			return(e)
		}
		n.var    = length(xy[1,])
		sample.size = length(xy[,1])
		sample   = matrix(nrow=sample.size,ncol=n.var)
		mean.std = matrix(nrow=2,ncol=n.var)
			{
			#----- Generate random number for the 1st distribution
			sample[,1] = xy[,1]
			for (i in 2:n.var) #J1
				{ 
						{	
									{   #J2
												sample[,i] = xy[,i]                           #Generate random for the ith variable
												for (j in 1:i)
													{
														mean.std[1,j] = mean(sample[,j])
														mean.std[2,j] = sd(sample[,j])	
													}
												reg.mat = matrix(0,nrow=(n.var-1),ncol=n.var)
												reg.mat = f.reg.mat(cor.matrix,mean.std,n.var)
												e.mat	= matrix(0,nrow=1,ncol=(n.var-1))
												e.mat	= f.e.mat(cor.matrix,mean.std,n.var)
												op.mat	= matrix(nrow=sample.size,ncol=i)
												t       = rnorm(sample.size,mean=0, sd = e.mat[1,(i-1)])        	
												for (j in 1:(i-1))
													{
														op.mat[,j] = sample[,j]
													}
												{
													y.d     = reg.mat[(i-1),1]								   											
													for (k in 1:(i-1))
															{
															  y.d = y.d + sample[,k]*reg.mat[(i-1),(k+1) ]	  
															}												                                        
												}
												op.mat[,i] = y.d + t*runif(1,0.75,1)
												op.mat	   = op.mat[order(op.mat[,(i)]),]
												op.mat[,i] = sort(sample[,i])													
												sample[,1:i] = op.mat[,1:i]												
									}
						}
				}
	return(sample)
			}
	} 
	
	 f.check <- function(n.var,matrix1,matrix2)		
	 {
		 return(sqrt(2*sqrt((sum((matrix1-matrix2)^2))/2)/n.var/(n.var-1)))
	 }
	#This function is used for plotting Temperature vs Error, Enabled if plotxy<-1
	plot.xy <- function(count,temperature,f.newr,f.oldr,f.bestr)
	{
		x  = 1:(count)
		y0 = temperature
		y1 = f.newr
		y2 = f.oldr
		y3 = f.bestr
		plot(x, y1,ylim=as.numeric(cbind(min(y1),max(y1))), axes=F, xlab="",ylab="Weighted Errors",col="dark grey", type="s",main="Errors vs Iterations")
		axis(2, pretty(as.numeric(cbind(min(y1) ,max(y1))),10))
		#par(new=T)
		axis(1,pretty(range(x),10))
		points(x, y2, type='l', col=" blue", xlab='x', ylab='y')
		points(x, y3, type='l', col="dark green", xlab='x', ylab='y')
		graphics::legend(max(x)*0.85,max(y1)*.95,c("New","Previous","Best"), lty=c(1),lwd=c(2.5),col=c("dark grey","blue","green"))
		#points(x, y0, type='l', col="red", xlab='x', ylab='y'
	}
	# initial parameters
	 n <- length(xy[,1])
	 n.var<- length(xy[1,])							#Count how many variables there are from the input
	 count 		 <- 0
	 xy.old		 <- update.xy(xy,cor.matrix)   #----UPDATED TEST
	 f.old		 <- f.check(n.var,cor(xy.old),cor.matrix) 
	 rate		 <- n
	 f.best 	 <- f.old
	 xy.best 	 <- xy.old
	 temperature <- NULL									#Used for storing the temperature of each iteration
	 f.bestr	 <- NULL									#Used for storing the sum of deviations of each iteration
	 f.newr		 <- NULL
	 f.oldr		 <- NULL
	 #loop 		 <- n*n.var	#n.var*n.var  #----UPDATED TEST
	 cols 	 <- (1:n.var)
	 ross	 <- (1:n)
	 t 		 <- 50/n  #init.cost			  #----UPDATED TEST							#Determine the inital temperature for SA. 
	 	 if(logs==1){
			 print(paste("Inital correlation"))
			 print(round(cor(xy.old),2))							
			 print(paste("Initial deviation",f.old))		
		}
	 repeat {
		s <- 0 
		#set.seed(0)
		k	<- 0
		loop	<- min(100 + n.var,rate)#loop,
			#Find the best solution in a neighbour region at t[j] temperature. Repeat Ntrials times before reducing the temperature 
			while ((k < loop) && (count < maxitr))
				{
					k 		<- k + 1
					count	<- count + 1
					column	<- sample(cols,1)								    #Randomly select n column
					rows	<- sample(ross,2)									#Randomly select 2 rows					
					swap1	<- xy.old[rows[1],column]
					swap2	<- xy.old[rows[2],column]
					xy.old[rows[1],column] <- swap2	#Swap the two indexes
					xy.old[rows[2],column] <- swap1
					cor.new	<- cor(xy.old)
					corl 			<- {cor.new - cor.matrix}
					f.new	<- sqrt(2*sqrt({sum(corl*corl)/2})/n.var/{n.var-1}) #Find the cost value
							if (plotxy == 1 ) {
								temperature	<-	cbind(temperature,t)	  	    #Record the current temperature
								f.newr		<-  cbind(f.newr,f.new)	   			#Record the current deviation		
								f.bestr		<-	cbind(f.bestr,f.best)
								f.oldr		<- 	cbind(f.oldr,f.old)
								#print(paste(count,"|fbest:",f.best,"| fold",f.old,"| fnew",f.new,"|t",t))
							}
					#Acceptance/Rejection Rule
					if (f.new<f.old) {						    #Is the new neighbour better than the older one? If YES than ACCEPT
							s		<- s + 1
							f.old   <- f.new
							if (f.new < f.best)				    #Is the new neighbour better than the current best solution?
								{
									k 			<- 1
									f.best		<- f.new
									xy.best		<- xy.old
									if (chart == 1) {chart.Correlation(xy.best)}
									if (all.equal(round(cor.new,round),round(cor.matrix,round),check.attributes=FALSE) == TRUE)
									{
										print(paste("Stop here"))
										s <- 0
										break
									}
								}
						} else
					if (exp(-(f.new-f.old)/t) >= runif(1))							#Else is the new neighbour be ACCEPTED with a Probability?
						{
							s		<- s + 1
							f.old   <- f.new
						} else
					{ 	xy.old[rows[1],column] <- swap1	#Return the two indexes
						xy.old[rows[2],column] <- swap2
					}
				}
		#Now reduce the temperature	
		t	<- t*c.rate    #rate of change ---UPDATED
		xy.old	<- xy.best
		f.old	<- f.best
		if (s==0) {break}
		if (count ==  maxitr) {break}
	}
	#close(pb)
	print(paste("Completed in",count,"iterations"))
	print(paste("The least deviation:",f.best))
	print(round(cor(xy.best),round))
	if (plotxy == 1) {plot.xy(count,temperature,f.newr,f.oldr,f.bestr)}
	return(list(cor(xy.best),count))
}

#### Test run
#### 2 variables
n = 10000
x = rnorm(n,0,1)
y = rpois(n,0.5)
xy = cbind(x,y)
cor.matrix<-diag(2)
cor.matrix[1,2] <- cor.matrix[2,1] <- 0.4
system.time(out<-sai(xy,cor.matrix))
sqrt(sum((cor(out)-cor.matrix)^2)/2)
system.time(out<-sai(xy,cor.matrix))

p = parallel(sai(xy,cor.matrix))
system.time(out<-collect(p))
#### Test run
#### 4 variables
 n = 1000
 x = rnorm (n,5,3)
 y = rpois (n,7)
 z = rexp(n,0.1)
 t = runif(n,0,100)
 xy = cbind(x,y,z,t)
		cor.matrix<-diag(4)
		cor.matrix[1,2] <- cor.matrix[2,1] <- 0.2
		cor.matrix[1,3] <- cor.matrix[3,1] <- 0.6
		cor.matrix[1,4] <- cor.matrix[4,1] <-  0
		cor.matrix[2,3] <- cor.matrix[3,2] <-  -0.4
		cor.matrix[2,4] <- cor.matrix[4,2] <-  0.5
		cor.matrix[3,4] <- cor.matrix[4,3] <- -0.7
system.time(out<-sai(xy,cor.matrix,maxitr=1e6,round=1))
sqrt(sum((cor(out)-cor.matrix)^2)/2)

#### Test run
#### 8 variables


library(evd)
library(PearsonDS)

n = 500
 var1 = rweibull(n,2.65,10.33)
 var2 = rextreme(n,qgamma,2.76,7.65)
 var3 = rlnorm(n,3.26,0.53)
 var4 = rbinom(n,19,0.46)
 var5 = rgamma(n,4.48,1.24)
 var6 = rpois(n,8.26)
 var7 = rpearsonV(n,shape=7.45,scale=60.15,location=1)
 var8 = rchisq(n,10)
 xy = cbind(var1,var2,var3,var4,var5,var6,var7,var8)
 cor.matrix<-diag(8)
		cor.matrix[1,2] <- cor.matrix[2,1] <- 0.901
		cor.matrix[1,3] <- cor.matrix[3,1] <- 0.684
		cor.matrix[1,4] <- cor.matrix[4,1] <-  0.567
		cor.matrix[1,5] <- cor.matrix[5,1] <-  -0.521
		cor.matrix[1,6] <- cor.matrix[6,1] <-  0.487
		cor.matrix[1,7] <- cor.matrix[7,1] <-  0.418
		cor.matrix[1,8] <- cor.matrix[8,1] <-  0.393
		cor.matrix[2,3] <- cor.matrix[3,2] <- 0.838
		cor.matrix[2,4] <- cor.matrix[4,2] <-  0.648
		cor.matrix[2,5] <- cor.matrix[5,2] <-  -0.570
		cor.matrix[2,6] <- cor.matrix[6,2] <-  0.577
		cor.matrix[2,7] <- cor.matrix[7,2] <-  0.519
		cor.matrix[2,8] <- cor.matrix[8,2] <-  0.483
		cor.matrix[3,4] <- cor.matrix[4,3] <-  0.866
		cor.matrix[3,5] <- cor.matrix[5,3] <-  -0.738
		cor.matrix[3,6] <- cor.matrix[6,3] <-  0.800
		cor.matrix[3,7] <- cor.matrix[7,3] <-  0.770
		cor.matrix[3,8] <- cor.matrix[8,3] <-  0.734
		cor.matrix[4,5] <- cor.matrix[5,4] <-  -0.910
		cor.matrix[4,6] <- cor.matrix[6,4] <-  0.938
		cor.matrix[4,7] <- cor.matrix[7,4] <-  0.857
		cor.matrix[4,8] <- cor.matrix[8,4] <-  0.877
		cor.matrix[5,6] <- cor.matrix[6,5] <-  -.919
		cor.matrix[5,7] <- cor.matrix[7,5] <-  -.788
		cor.matrix[5,8] <- cor.matrix[8,5] <-  -.822
		cor.matrix[6,7] <- cor.matrix[7,6] <-  0.926
		cor.matrix[6,8] <- cor.matrix[8,6] <-  0.940
		cor.matrix[7,8] <- cor.matrix[8,7] <-  0.942
system.time(out<-sai(xy,cor.matrix,maxitr=1e6))
sqrt(sum((cor(out)-cor.matrix)^2)/2)


round(cor(cornode(xy,target=cor.matrix),2))

system.time(for (i in 1:100000) cor(xy))
p = mcparallel(for ( i in 1:1000) cor(xy))
system.time(out<-collect(p))