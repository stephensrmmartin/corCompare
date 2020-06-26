
functions {
#include fns/tbeta.stan

  matrix compare_all(vector par) {
    int P = rows(par);
    matrix[P, P] out;
    for(i in 1:P) {
      for(j in 1:i) {
	out[i, j] = par[i] - par[j];
	out[j, i] = -out[i, j];
      }
    }
    return(out);
  }

  matrix to_cov(vector sigma, real rho) {
    matrix[2, 2] cor = diag_matrix(rep_vector(1.0, 2));
    matrix[2, 2] out;
    cor[1,2] = rho;
    cor[2,1] = rho;
    out = quad_form_diag(cor, sigma);
    return(out);
  }

}

data {
  // Data required
  int N; // Total number of rows (long format)
  int group[N]; // "Grouping" variable (i.e., a data-subset of interest).
  vector[N] y; // Outcome variable with which all X's will be correlated.
  vector[N] x; // Variables to correlate with y

  // Settings
  int prior_only; // Whether to sample from prior or not.
  vector[max(group)] tbeta_mu; // Transformed beta-mean(s) for each group.
  vector[max(group)] tbeta_psi; // Transformed beta-psi(s) for each group.
}

transformed data {
  int G = max(group);
  matrix[N, 2] yx = append_col(y, x);
  row_vector[2] yx_arr[N];
  for(n in 1:N) {
    yx_arr[n] = yx[n];
  }
}

parameters {
  vector[2] mu[G];
  vector<lower = 0>[2] sigma[G];
  vector<lower = -1, upper = 1>[G] rho;
}

transformed parameters {
  matrix[2, 2] cov_g[G];
  for(g in 1:G) {
    cov_g[g] = to_cov(sigma[g], rho[g]);
  }
  
}

model {
  for(g in 1:G) {
    rho[g] ~ tbeta(tbeta_mu[g], tbeta_psi[g]);
    // Uniform mean; uniform SD; can check this later
  }

  if(!prior_only) {
    for(n in 1:N) {
      yx_arr[n] ~ multi_normal(mu[group[n]], cov_g[group[n]]);
    }
  }
  
}

generated quantities {
  /* Disabling for now; needs to unpack the array-of-vec to vec */
  // matrix[G, G] mu_diff = compare_all(mu);
  // matrix[G, G] sigma_diff = compare_all(sigma);
  matrix[G, G] rho_diff = compare_all(rho);
}
