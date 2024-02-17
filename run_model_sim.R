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

model_code = "functions {
  vector gen_relaxed_dists(real[] rho, matrix ancestor_lens, int[] child, int B, int N) {
    vector[N] relaxed_dists = rep_vector(0,N);
    for (b in 1:B) {
      relaxed_dists[child[b]] = dot_product(ancestor_lens[b,],to_vector(rho));
    }
    return(relaxed_dists);
  }
  real geo_llik(int[,] mrca22, vector relaxed_dists, vector tip_vals, real root_val, int T, real sigma) {
    matrix[T,T] cov22;
    vector[T] mu2;for (i in 1:(T-1)) {
      cov22[i,i] = relaxed_dists[mrca22[i,i]];
      for (j in (i+1):T) {
        cov22[i,j] = relaxed_dists[mrca22[i,j]];
        cov22[j,i] = cov22[i,j];
      }
    }
    cov22[T,T] = relaxed_dists[mrca22[T,T]];
    mu2 = rep_vector(root_val,T);
    return(multi_normal_cholesky_lpdf(tip_vals|mu2,cholesky_decompose(sigma*cov22)));
  }
  matrix evprob(real z, real alpha, real beta) {
    matrix[2,2] P;
    P[1,1] = (beta/(alpha+beta)) + (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[1,2] = (alpha/(alpha+beta)) - (alpha/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,1] = (beta/(alpha+beta)) - (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    P[2,2] = (alpha/(alpha+beta)) + (beta/(alpha+beta)*exp(-(alpha+beta)*z));
    return P;
  }
  //compute likelihood via Felsenstein's Pruning Algorithm
  real pruning_vec(int N, int B, int[] child, int[] parent, real[] brlen, matrix tiplik, vector alpha, vector beta) {
    vector[2] pi;                         //stationary probability
    matrix[N,2] lambda = log(tiplik);     //likelihoods at tips+nodes
    real llik;                       //log likelihoods for each family for feature d
    for (b in 1:B) {
      matrix[2,2] P = evprob(brlen[b], alpha[b], beta[b]); //via matrix exponentiation
      for (d in 1:2) {
        lambda[parent[b],d] += log(dot_product(P[d],exp(lambda[child[b]])));
      }
    }
    pi[1] = -log(2) + lambda[parent[B],1];
    pi[2] = -log(2) + lambda[parent[B],2];
    llik = log_sum_exp(pi);
    return(llik);
  }
}
data {
  int<lower=1> N;                           //number of nodes
  int<lower=1> P;                           //number of present geo values
  int<lower=1> B;                           //number of branches
  int<lower=1> D;                           //number of segments
  int<lower=1> J;                           //number of features
  int<lower=1> child[B];                    //child of each branch
  int<lower=1> parent[B];                   //parent of each branch
  real<lower=0> brlen[B];                   //length of each branch
  matrix[N,D*2] tiplik;                     //likelihoods for data at tips+internal nodes in of each tree
  vector[P] lon;
  vector[P] lat;
  matrix[B,B] ancestor_lens;
  int mrca22[P,P];
}
parameters {
  real<lower=0> sigma;
  real<lower=0> rho[B];
  real lon_root;
  real lat_root;
  real<lower=0> tau;
  real alpha_s;
  real alpha_p;
  vector[D] alpha_seg_s;
  vector[D] alpha_seg_p;
  real beta_s;
  real beta_p;
  vector[D] beta_seg_s;
  vector[D] beta_seg_p;
  real delta_s[B];
  real delta_p[B];
  vector[B] delta_seg_s[D];
  vector[B] delta_seg_p[D];
  real<lower=0> std_delta_s;
  real<lower=0> std_delta_p;
  real<lower=0> std_delta_seg_s[D];
  real<lower=0> std_delta_seg_p[D];
}
transformed parameters {
  vector[D] log_lik;
  vector[N] relaxed_dists = gen_relaxed_dists(rho,ancestor_lens,child,B,N);
  real tip_lon_loglik = geo_llik(mrca22,relaxed_dists,lon,lon_root,P,sigma);
  real tip_lat_loglik = geo_llik(mrca22,relaxed_dists,lat,lat_root,P,sigma);
  {vector[B] s[D];
   vector[B] p[D];
  for (d in 1:D) {
    for (b in 1:B) {
      s[d,b] = inv_logit(alpha_s + alpha_seg_s[d] + (beta_s + beta_seg_s[d])*rho[b] + delta_s[b]*std_delta_s + delta_seg_s[d,b]*std_delta_seg_s[d])*tau;
      p[d,b] = inv_logit(alpha_p + alpha_seg_p[d] + (beta_p + beta_seg_p[d])*rho[b] + delta_p[b]*std_delta_p + delta_seg_p[d,b]*std_delta_seg_p[d]);
    }
    log_lik[d] = pruning_vec(N,B,child,parent,brlen,tiplik[,((2*d)-1):(2*d)],p[d,].*s[d,],(1-p[d,]).*s[d,]);
  }
  }
}
model {
  rho ~ normal(0,1);
  sigma ~ normal(0,1);
  lon_root ~ uniform(0,360);
  lat_root ~ uniform(-90,90);
  tau ~ uniform(0,10);
  alpha_s ~ normal(0,1);
  alpha_p ~ normal(0,1);
  alpha_seg_s ~ normal(0,1);
  alpha_seg_p ~ normal(0,1);
  beta_s ~ normal(0,1);
  beta_p ~ normal(0,1);
  beta_seg_s ~ normal(0,1);
  beta_seg_p ~ normal(0,1);
  delta_s ~ normal(0,1);
  delta_p ~ normal(0,1);
  for (d in 1:D) {
    delta_seg_s[d] ~ normal(0,1);
    delta_seg_p[d] ~ normal(0,1);
  }
  std_delta_s ~ normal(0,1);
  std_delta_p ~ normal(0,1);
  std_delta_seg_s ~ normal(0,1);
  std_delta_seg_p ~ normal(0,1);
  target += tip_lon_loglik;
  target += tip_lat_loglik;
  target += sum(log_lik);
}"

data.list <- list(
  N=N,
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

fit <- stan(model_code=model_code,data=data.list,cores=n_cores,control=list(adapt_delta=.99))

saveRDS(fit,file=paste('fitted_models_sim/model_fit_sim_',t,'.RDS',sep=''))
