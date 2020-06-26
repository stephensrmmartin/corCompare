

  /*
    Beta density (mu/psi parameterization)
    psi = (a + b)
    mu = a / psi
    a = mu(psi)
    b = psi - mu(psi)
      = psi(1 - mu)
      = psi - a
   
    E(x) = mu
    V(x) = mu(1 - mu) / (psi + 1)
    psi= [mu(1 - mu) / V(x)] - 1

   */ 
  real mupsibeta_lpdf(real x, real mu, real psi) {
    real a = mu * psi;
    real b = psi - a;
    real out = beta_lpdf(x | a, b);
    return(out);
  }

  /*
    Beta density (mu/variance parameterization)
    E(x) = mu
    V(x) = v
   */
  real muvarbeta_lpdf(real x, real mu, real v) {
    real psi = mu*(1 - mu) / (v) - 1;
    return(mupsibeta_lpdf(x | mu, psi));
  }

  /*
    Transformed-beta density (mu/psi parameterization)
    tbeta(x, mu, psi) = beta(x / 2 + .5, mu / 2 + .5, psi)
    E(x) = mu
    V(x) = (1 + mu)(1 - mu) / (psi + 1)
   */
  real tbeta_lpdf(real x, real mu, real psi) {
    real mu_beta = (mu + 1.0)*.5;
    real x_beta = (x + 1) * .5;
    real out = mupsibeta_lpdf(x_beta | mu_beta, psi) + log(.5);

    return(out);
  }
