library(statnet)
library(mgcv)
library(parallel)

num_iter <- 1000
time_steps <- 100

normalize_stochastic <- function(A) {
  vertices <- length(A[1,])
  for ( i in 1:vertices ) {		# stochasticize the matrix
    links <- sum(A[i,])
    if (links==0) {
      A[i,i]=1			        # 1s on diagonal if no links
    } else {
      A[i,] = A[i,] / links		# otherwise, normalize so sum is 1
    }
  }
  return(A)
}

reflected_appraisal <- function(A) {
  vertices<- length(A[1,])
  AT <- t(A)                  		# transpose for left eigenvector
  C <- A - diag(diag(A))		# remove the diagonal for later
  val <- Re(zapsmall(eigen(AT)$values)) # find all eigenvalues
  vec <- Re(zapsmall(eigen(AT)$vectors))# find all eigenvectors
  if (sum(vec < 0)>0) { vec <- -vec }
  newvec <- vector('double',vertices)
  newvec[] <- 0
  for (i in 1:vertices) {
    if (val[i]==1) {
      temp <- vec[,i]/sum(vec[,i])	# normalize to sum 1
      newvec <- newvec + temp		# add to social power vector
    }
  }
  W <- diag(newvec) + (diag(vertices) - diag(newvec))%*%C
  return(W)				# form W matrix
}

data(faux.mesa.high)			# load faux mesa high network and fit ERGM model

model <- ergm(faux.mesa.high ~ edges + nodematch('Sex',diff=TRUE) + nodematch('Race') + nodematch('Grade',diff=TRUE) + gwesp,control = control.ergm(MCMLE.maxit=100,parallel = detectCores()))

C <- as.matrix(faux.mesa.high)		# get adjacency matrix

vertices <- length(C[1,])

y0 <- runif(vertices,0,1) 		#rnorm(vertices,0,1) # initialize opinions at time 0

y <- matrix(nrow=vertices,ncol=time_steps)
y_hat <- matrix(0,nrow=vertices,ncol=time_steps)
y2 <- matrix(nrow=vertices,ncol=time_steps)

A <- normalize_stochastic(C)		# normalize adjacency to C row stochastic
W <- A
norm_diff <- 10
while (norm_diff > 0.00000001) {		# iterate reflected appraisal to convergence
  temp <- diag(W)
  W <- reflected_appraisal(W)
  norm_diff <- norm((diag(W)-temp),"2")
}

y[,1] <- W%*%y0
y2[,1] <- W%*%y0

for (t in 2:time_steps) {		# iterate opinion exchange over the real network. Also, do the new-network-each-time model
  y[,t] <- A%*%y[,t-1]
  y2[,t] <- normalize_stochastic(as.matrix(simulate(model)))%*%y2[,t-1]
}

for (i in 1:num_iter) {			# Markov Chain simulation of mean opinion
  B <- normalize_stochastic(as.matrix(simulate(model)))
  W <- B				# simulate new network and iterate 1. social power and 2. opinions
  norm_diff <- 10
  while (norm_diff < 0.0001) {
    temp <- diag(W)
    W <- reflected_appraisal(W)
    norm_diff <- norm((diag(W)-temp),"2")
  }
    y_temp <- W%*%y0
  y_hat[,1] <- y_hat[,1] + y_temp
  for (t in 2:time_steps) {
    y_temp <- W%*%y_temp
    y_hat[,t] <- y_hat[,t] + y_temp
  }
}

y_hat <- (1/num_iter)*y_hat		# average the resulting opinion vectors

resids <- y - y_hat
norm_resids <- vector(mode="double",length=10)

for (t in 1:time_steps) {
  norm_resids[t] <- norm(resids[,t],'2')
}

avg_op <- vector('double')
avg_op_hat <- vector('double')
avg_op_revised <- vector('double')
med_op <- vector('double')
med_op_hat <- vector('double')

for (t in 1:time_steps) {		# track and plot mean and median opinions over time
  avg_op[t] <- mean(y[,t])
  avg_op_hat[t] <- mean(y_hat[,t])
  avg_op_revised[t] <- mean(y2[,t])
  
  med_op[t] <- median(y[,t])
  med_op_hat[t] <- median(y_hat[,t])
}

m <- vector('double',length=time_steps)
m[] <- mean(y0)
act <- vector('double',length=time_steps)
act[] <- 0.5
plot(avg_op_revised,xlab="Time Step",ylab="Opinion",ylim = c(min(c(min(avg_op_hat),min(avg_op),min(avg_op_revised),min(med_op),min(med_op_hat))),max(c(max(avg_op_hat),max(avg_op),max(avg_op_revised),max(med_op),max(med_op_hat)))))
lines(m,col="red")
lines(m,col='yellow')
points(1:time_steps,avg_op_hat,col="blue",pch="+")
points(1:time_steps,avg_op,col="dark red",pch="*")
points(1:time_steps,med_op,col="purple")
points(1:time_steps,med_op_hat)
points(1:time_steps,avg_op_revised,col="green",pch="o")
legend("right",legend=c("Predicted Average","Actual Average","Predicted Median","Actual Median","DeGroot Reforming","Average of Initial States"),col=c("blue","red","black","purple","green","yellow"),pch=c("+","*","o","o","o","-"))


