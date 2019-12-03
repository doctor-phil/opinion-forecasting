library(statnet)
library(mgcv)
library(igraph)

data(faux.mesa.high)

model <- ergm(faux.mesa.high ~ edges + nodematch('Sex',diff=TRUE) + nodematch('Race',diff=TRUE) + nodematch('Grade',diff=TRUE) + gwesp)

C <- as.matrix(faux.mesa.high)

G <- graph_from_adjacency_matrix(C)

num_components <- components(faux.mesa.high)
vertices <- length(C[1,])

subgraph

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
