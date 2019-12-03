library(statnet)
library(mgcv)
library(parallel)
library(ggplot2)

num_iter <- 1000
time_steps <- 100

normalize_stochastic <- function(A) {
  vertices <- length(A[1,])
  for ( i in 1:vertices ) {			#stochasticize the matrix
    links <- sum(A[i,])
    if (links==0) {
      A[i,i]=1			#1s on diagonal if no links
    } else {
      A[i,] = A[i,] / links		#otherwise, normalize so sum is 1
    }
  }
  return(A)
}

data(faux.mesa.high)

model <- ergm(faux.mesa.high ~ edges + nodematch('Sex',diff=TRUE) + nodematch('Race') + nodematch('Grade',diff=TRUE) + gwesp,control = control.ergm(MCMLE.maxit=100,parallel = detectCores()))

C <- as.matrix(faux.mesa.high)

vertices <- length(C[1,])

y0 <- rnorm(vertices,0,1)

y <- matrix(nrow=vertices,ncol=time_steps)
y_hat <- matrix(0,nrow=vertices,ncol=time_steps)
y2 <- matrix(nrow=vertices,ncol=time_steps)

A <- normalize_stochastic(C)

y[,1] <- A%*%y0
y2[,1] <- A%*%y0

for (t in 2:time_steps) {
  y[,t] <- A%*%y[,t-1]
  y2[,t] <- normalize_stochastic(as.matrix(simulate(model)))%*%y2[,t-1]
}

for (i in 1:num_iter) {
  B <- normalize_stochastic(as.matrix(simulate(model)))
  y_temp <- B%*%y0
  y_hat[,1] <- y_hat[,1] + y_temp
  for (t in 2:time_steps) {
    y_temp <- B%*%y_temp
    y_hat[,t] <- y_hat[,t] + y_temp
  }
}

y_hat <- (1/num_iter)*y_hat

resids <- y - y_hat
norm_resids <- vector(mode="double",length=10)

for (t in 1:time_steps) {
  norm_resids[t] <- norm(resids[,t],'2')
}

avg_op <- vector('double')
avg_op_hat <- vector('double')
avg_op_revised <- vector('double')

for (t in 1:time_steps) {
  avg_op[t] <- mean(y[,t])
  avg_op_hat[t] <- mean(y_hat[,t])
  avg_op_revised[t] <- mean(y2[,t])
}
m <- vector('double',length=time_steps)
m[] <- mean(y0)
plot(avg_op_revised,ylim = c(min(c(min(avg_op_hat),min(avg_op),min(avg_op_revised))),max(c(max(avg_op_hat),max(avg_op),max(avg_op_revised)))))
lines(m,col="red")
points(1:100,avg_op_hat,col="blue",pch="+")
points(1:100,avg_op,col="dark red",pch="*")
points(1:100,avg_op_revised,col="green",pch="o")