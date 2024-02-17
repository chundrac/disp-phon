library(phytools)
library(rstan)

n_cores = 16

options(mc.cores = parallel::detectCores())
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

rstan_options(threads_per_chain = as.integer(n_cores/4))

t = as.integer(commandArgs(trailingOnly = TRUE)[1])

set.seed(t)

data.df <- read.csv('character_data.tsv',sep='\t',row.names=1)

data.df <- data.df[,colSums(data.df)>1]

seg.feats <- read.csv('segmental_features.tsv',sep='\t',row.names=1)
J <- ncol(seg.feats)

trees <- read.nexus('proxy_matching/austronesian_proxy.nex')
i = sample(c(1:length(trees)),1)
tree <- trees[[i]]
tree <- keep.tip(tree,which(tree$tip.label %in% rownames(data.df)))
tree <- reorder.phylo(tree,'pruningwise')

data.df.i <- data.df[tree$tip.label,]

sim.data <- NULL

for (i in 1:10) {
  alpha <- rgamma(1,1,1)
  beta <- rgamma(1,1,1)
  simhist <- sim.history(tree,Q=rbind(c(-alpha,alpha),c(beta,-beta)))
  states.i <- simhist$states
  sim.data <- cbind(sim.data,as.numeric(states.i)-1)
}

rownames(sim.data) <- rownames(data.df.i)

print('done')

#simulate continuous coordinates

root.lon <- runif(1,0,300)
root.lat <- runif(1,-90,90)

sigma <- abs(rnorm(1,0,1))
phy.cov <- vcv.phylo(tree)*sigma

z.lon <- rnorm(nrow(data.df.i),0,1)
z.lat <- rnorm(nrow(data.df.i),0,1)

lon <- phy.cov %*% z.lon
lat <- phy.cov %*% z.lat

parent <- tree$edge[,1]
child <- tree$edge[,2]
b.lens <- tree$edge.length/1000
B <- length(b.lens)
N <- length(unique(c(parent,child)))
T <- length(child[which(!child %in% parent)])

orig.states <- data.df.i

states <- cbind(orig.states,1-orig.states)
node.states <- as.data.frame(matrix(1,nrow=tree$Nnode,ncol=ncol(states)))
colnames(node.states) <- colnames(states)
states <- rbind(states,node.states)
bin.states <- states[,c(rbind(c(1:(ncol(states)/2)),c(1:(ncol(states)/2))+(ncol(states)/2)))]
bin.states[is.na(bin.states)] <- 1
tip.lik <- bin.states

D <- ncol(tip.lik)/2

langs <- rownames(data.df.i)

lon <- lon[langs,1]
lat <- lat[langs,1]

mrca.matrix <- mrca(tree,full=T)

#get branch lengths leading to specific node
child.branch.lengths <- rep(0,nrow(mrca.matrix))
for (i in 1:length(child)) {
  child.branch.lengths[child[i]] <- b.lens[i]
}

#ancestral branches of each possible mrca; sum(child.branch.lengths*ancestors) gives the distance from root to some mrca node
ancestors <- matrix(0,nrow(mrca.matrix),ncol(mrca.matrix))
for (i in 1:nrow(mrca.matrix)) {
  for (j in unique(mrca.matrix[i,])) {
    ancestors[i,j] = 1
  }
}

lens.to.root <- ancestors*matrix(child.branch.lengths,nrow=nrow(ancestors),ncol=length(child.branch.lengths),byrow=T)
lens.to.root <- lens.to.root[child,child]

root.index <- T+1

missing.inds <- which(is.na(lon))
missing.inds <- missing.inds[which(missing.inds!=root.index)]
present.inds <- which(!is.na(lon))

X = length(missing.inds)
Y = length(present.inds)

Sigma22 = matrix(nrow=Y,ncol=Y)
for (x in 1:Y) {
  for (y in 1:Y) {
    Sigma22[x,y] <- mrca.matrix[present.inds[x],present.inds[y]]
  }
}

lon <- lon[!is.na(lon)]
lat <- lat[!is.na(lat)]

data.list <- list(
  N=N,
  M=X,
  P=Y,
  B=B,
  D=D,
  J=J,
  child=child,
  parent=parent,
  brlen=b.lens,
  tiplik=tip.lik,
  lon=lon,
  lat=lat,
  ancestor_lens=lens.to.root,
  mrca22=Sigma22
)

fit <- stan(file='model_joint.stan',data=data.list,cores=n_cores,control=list(adapt_delta=.99))

saveRDS(fit,file=paste('model_fits_sim/model_fit_sim_',t,'.RDS',sep=''))
