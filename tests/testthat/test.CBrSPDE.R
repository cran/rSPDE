context("CBrSPDE.operators")

test_that("Checking covariances of CBrSPDE",{
  set.seed(123)
  nobs <- 100
  s <- seq(from = 0, to = 1, length.out = nobs)
  fem <- rSPDE.fem1d(s)
  kappa <- 40
  sigma <- 1
  d <- 1
  for (nu in c(0.8,1.7,2.6)){
  op2 <- matern.operators(C=fem$C, G=fem$G,nu = nu,kappa = kappa,sigma = sigma,
                                  d=1,m = 2)
  v <- t(rSPDE.A1d(s,0.5))
  c.true <- matern.covariance(abs(s - 0.5), kappa, nu, sigma)
  Q <- rspde.matern.precision(kappa=kappa,nu=nu,sigma=sigma,rspde_order=2,dim=1,
                              fem_mesh_matrices = op2$fem_mesh_matrices)
  A <- Diagonal(nobs)
  Abar <- cbind(A,A,A)
  w <- rbind(v,v,v)
  c.approx2 <- (Abar)%*%solve(Q,w)
  res <- sum((c.approx2-c.true)^2)
  expect_equal(res, 0, tolerance = 1e-4)
}
})
  
test_that("Checking loglike of CBrSPDE", {
  set.seed(123)
  nobs <- 100
  s <- seq(from = 0, to = 1, length.out = nobs)
  fem <- rSPDE.fem1d(s)
  kappa <- 40
  sigma <- 1
  d <- 1
  for (nu in c(0.8,1.7,2.6)){
    op2 <- matern.operators(C=fem$C, G=fem$G,nu = nu,kappa = kappa,sigma = sigma,
                                    d=1,m = 2)
    A <- Diagonal(nobs)
    sim_data = A %*% simulate(op2) + rnorm(dim(A)[1], sd=0.1)
    loglike2 = rSPDE.matern.loglike(object=op2, Y = sim_data, A=A, 
                                      sigma.e=0.1,user_nu = nu, 
                                      user_kappa = kappa, 
                                      user_sigma=sigma)
    
    op1 <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                            G = fem$G, C = fem$C, d = 1, type="operator")
    loglike1 = rSPDE.matern.loglike(op1, sim_data, A, sigma.e=0.1)
    
    expect_equal(loglike1,loglike2, tol=0.3)
  }
})

test_that("Checking Predict of CBrSPDE", {
  set.seed(123)
  nobs <- 100
  s <- seq(from = 0, to = 1, length.out = nobs)
  fem <- rSPDE.fem1d(s)
  kappa <- 40
  sigma <- 1
  d <- 1
  for (nu in c(0.8,1.7,2.6)){
    Aprd <- rSPDE.A1d(s,0.5)
    op2 <- matern.operators(C=fem$C, G=fem$G,nu = nu,kappa = kappa,sigma = sigma,
                                    d=1,m = 2)
    A <- Diagonal(nobs)
    sim_data = A %*% simulate(op2) + rnorm(dim(A)[1], sd=0.1)
    predict2 = predict(object=op2, Y = sim_data, A=A, sigma.e=0.1,
                       Aprd = Aprd,
                       compute.variances=TRUE)
    
    op1 <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                            G = fem$G, C = fem$C, d = 1,
                            type = "operator")
    predict1 = predict(object=op1, Y=sim_data, A=A, Aprd=Aprd,  sigma.e=0.1,
                       compute.variances = TRUE)
    
    expect_equal(as.double(predict1$mean),as.double(predict2$mean), tol=0.02)
    expect_equal(as.double(predict1$variance),as.double(predict2$variance), tol=0.002)
  }
})

test_that("Checking loglike of CBrSPDE with replicates", {
  set.seed(123)
  nobs <- 100
  s <- seq(from = 0, to = 1, length.out = nobs)
  fem <- rSPDE.fem1d(s)
  kappa <- 40
  sigma <- 1
  d <- 1
  for (nu in c(0.8,1.7,2.6)){
    op2 <- matern.operators(C=fem$C, G=fem$G,nu = nu,kappa = kappa,sigma = sigma,
                                    d=1,m = 2)
    A <- Diagonal(nobs)
    sim_data1 = A %*% simulate(op2) + rnorm(dim(A)[1], sd=0.1)
    sim_data2 = A%*%simulate(op2) + rnorm(dim(A)[1], sd=0.1)
    sim_data = cbind(sim_data1,sim_data2)
    loglike2 = rSPDE.matern.loglike(object=op2, Y = sim_data, A=A, sigma.e=0.1)
    
    op1 <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                            G = fem$G, C = fem$C, d = 1,
                            type="operator")
    loglike1 = rSPDE.matern.loglike(op1, sim_data, A, sigma.e=0.1)
    
    expect_equal(loglike1,loglike2, tol=0.3)
  }
})