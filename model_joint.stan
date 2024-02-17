functions {
  real gc_dist(real lon1, real lat1, real lon2, real lat2) {
    real lon1r = (lon1/180)*pi();
    real lat1r = (lat1/180)*pi();
    real lon2r = (lon2/180)*pi();
    real lat2r = (lat2/180)*pi();
    real dlon = lon2r - lon1r;
    real dlat = lat2r - lat1r;
    real a = square(sin(dlat/2)) + cos(lat1r) * cos(lat2r) * square(sin(dlon/2));
    real c = 2 * asin(sqrt(a));
    real r = 6371;
    return(c * r);
  }
  vector gen_relaxed_dists(real[] rho, matrix ancestor_lens, int[] child, int B, int N) {
    vector[N] relaxed_dists = rep_vector(0,N);
    for (b in 1:B) {
      relaxed_dists[child[b]] = dot_product(ancestor_lens[b,],to_vector(rho));
    }
    return(relaxed_dists);
  }
  matrix make_sigma_square(int[,] mrca, vector relaxed_dists, real sigma, int T) {
    matrix[T,T] cov;
    for (i in 1:(T-1)) {
      cov[i,i] = relaxed_dists[mrca[i,i]];
      for (j in (i+1):T) {
        cov[i,j] = relaxed_dists[mrca[i,j]];
        cov[j,i] = cov[i,j];
      }
    }
    cov[T,T] = relaxed_dists[mrca[T,T]];
    return(cov*sigma);
  }
  matrix make_sigma_nonsquare(int[,] mrca, vector relaxed_dists, real sigma, int N, int T) {
    matrix[N,T] cov;
    for (i in 1:N) {
      for (j in 1:T) {
        cov[i,j] = relaxed_dists[mrca[i,j]];
      }
    }
    return(cov*sigma);
  }
  vector latent_geo(matrix cov11, matrix cov22, matrix cov12, vector tip_vals, real root_val, int T, int N, vector eta) {
    matrix[T,T] invcov22;
    vector[N] f;
    vector[N] mu1;
    vector[T] mu2;
    vector[N] mu;
    matrix[N,N] Sigma;
    mu1 = rep_vector(root_val,N);
    mu2 = rep_vector(root_val,T);
    invcov22 = inverse(cov22);
    mu = mu1 + cov12*invcov22*(tip_vals-mu2);
    Sigma = cov11 - cov12*invcov22*cov12';
    f = mu + cholesky_decompose(Sigma + diag_matrix(rep_vector(1e-10,N)))*eta;
    return(f);
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
  int<lower=0> M;                           //number of missing geo values
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
  int mrca11[M,M];
  int mrca12[M,P];
}
parameters {
  real<lower=0> sigma_lon;
  real<lower=0> sigma_lat;
  real<lower=0> rho_lon[B];
  real<lower=0> rho_lat[B];
  real lon_root;
  real lat_root;
  vector[M] eps_lon;
  vector[M] eps_lat;
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
  vector[B] eps;
  real<lower=0> std_delta_s;
  real<lower=0> std_delta_p;
  real<lower=0> std_delta_seg_s[D];
  real<lower=0> std_delta_seg_p[D];
}
transformed parameters {
  vector[D] log_lik;
  vector[B] dist;
  vector[N] all_lon;
  vector[N] all_lat;
  real tip_lon_loglik;
  real tip_lat_loglik;
  {
    vector[N] relaxed_dists_lon = gen_relaxed_dists(rho_lon,ancestor_lens,child,B,N);
    vector[N] relaxed_dists_lat = gen_relaxed_dists(rho_lat,ancestor_lens,child,B,N);
    matrix[P,P] cov22_lon = make_sigma_square(mrca22, relaxed_dists_lon, sigma_lon, P);
    matrix[P,P] cov22_lat = make_sigma_square(mrca22, relaxed_dists_lat, sigma_lat, P);
    matrix[M,M] cov11_lon = make_sigma_square(mrca11, relaxed_dists_lon, sigma_lon, M);
    matrix[M,M] cov11_lat = make_sigma_square(mrca11, relaxed_dists_lat, sigma_lat, M);
    matrix[M,P] cov12_lon = make_sigma_nonsquare(mrca12, relaxed_dists_lon, sigma_lon, M, P);
    matrix[M,P] cov12_lat = make_sigma_nonsquare(mrca12, relaxed_dists_lat, sigma_lat, M, P);
    tip_lon_loglik = multi_normal_cholesky_lpdf(lon|rep_vector(lon_root,P),cholesky_decompose(cov22_lon));
    tip_lat_loglik = multi_normal_cholesky_lpdf(lat|rep_vector(lat_root,P),cholesky_decompose(cov22_lat));
    all_lon[1:P] = lon;
    all_lon[P+1] = lon_root;
    all_lon[(P+2):N] = latent_geo(cov11_lon, cov22_lon, cov12_lon, lon, lon_root, P, M, eps_lon);
    all_lat[1:P] = lat;
    all_lat[P+1] = lat_root;
    all_lat[(P+2):N] = latent_geo(cov11_lat, cov22_lat, cov12_lat, lat, lat_root, P, M, eps_lat);
  }
  for (b in 1:B) {
    dist[b] = log(gc_dist(lon_root, lat_root, all_lon[child[b]], all_lat[child[b]]));
  }
  {vector[B] s[D];
   vector[B] p[D];
  for (d in 1:D) {
    for (b in 1:B) {
      s[d,b] = inv_logit(alpha_s + alpha_seg_s[d] + (beta_s + beta_seg_s[d])*dist[b] + delta_s[b]*std_delta_s + delta_seg_s[d,b]*std_delta_seg_s[d])*tau;
      p[d,b] = inv_logit(alpha_p + alpha_seg_p[d] + (beta_p + beta_seg_p[d])*dist[b] + delta_p[b]*std_delta_p + delta_seg_p[d,b]*std_delta_seg_p[d]);
    }
    log_lik[d] = pruning_vec(N,B,child,parent,brlen,tiplik[,((2*d)-1):(2*d)],p[d,].*s[d,],(1-p[d,]).*s[d,]);
  }
  }
}
model {
  rho_lon ~ normal(0,1);
  rho_lat ~ normal(0,1);
  sigma_lon ~ normal(0,1);
  sigma_lat ~ normal(0,1);
  lon_root ~ uniform(0,360);
  lat_root ~ uniform(-90,90);
  eps_lon ~ normal(0,1);
  eps_lat ~ normal(0,1);
  tau ~ uniform(0,20);
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
}
