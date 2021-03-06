\dontrun{

## Simulate Poisson cylinder system,
## generate a non-overlapping system by RSA,
## densify a cluster of cylinders

library(unfoldr)

## Unless MS-Windows platform	
# library(parallel)
# options(simLife.mc=2L)

lam <- 35
box <- list("xrange"=c(0,3),"yrange"=c(0,3),"zrange"=c(0,9))

## Spheroids of constant sizes
theta <- list("size"=list(.5),"shape"=list("radius"=0.1),
		"orientation"=list("kappa"=0.1))

S <- simPoissonSystem(theta,lam,size="const",type="cylinders",
		orientation="rbetaiso",box=box,pl=1,label="P")

## secondary phase: particles as spheres
F <- simPoissonSystem(list("size"=list(0.075)), lam=5, size="const",
		type="spheres",box=box, pl=1, label="F")

## apply RSA
S2 <- rsa(S,F,pl=101,verbose=TRUE)
####################################################################
## Optional: 3D visualization of cylinder projection areas
####################################################################

#require(rgl)
#id <- c(1,5,9,32,10)
#cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
#cylinders3d(S2[id], box, col=cols)	
#P <- getCylinderProjection(S2[id], B=c(0,1,0,1,1), draw=TRUE, conv = TRUE, np=20)
#P <- getCylinderProjection(S2[id], B=c(0,0,0,0,0), draw=TRUE, conv = TRUE, np=20)
#P <- getCylinderProjection(S2[id], B=c(1,1,1,1,1), draw=TRUE, conv = TRUE, np=20)

## construct clusters
CL <- simPoissonSystem(list("size"=list(0.35)), lam=0.1, size="const",
		type="spheres", box=box, pl=1, label="F")

CLUST <- simCluster(S2, CL, cond=list("eps"=1e-7,"minSize"=1L), verbose=TRUE, pl=1)
	
## densify
ctrl <- list(threshold.stop=0.01,max.call=5000,verbose=FALSE)
RET <- densifyCluster(S2, CLUST, ctrl, weight=10, cores = 1L, cl = NULL)	
G <- RET$cluster

####################################################################
## Optional: 3D visualization of densified sphere clusters
####################################################################

#open3d()
#lapply(CLUST,function(x) cylinders3d(S2[x$id],box=box,col=cols))	
#X <- do.call(rbind,lapply(CLUST, function(x) c(x$center,x$r)))
#invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))
#
### draw densified cluster
#open3d()
#invisible(lapply(G,function(x) { cylinders3d(x,box=box,col=cols) }))
#invisible(lapply(CLUST, function(x) rgl::spheres3d(X[,1:3],radius=X[,4],col="gray",alpha=0.2)))

}
