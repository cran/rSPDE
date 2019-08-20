## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(rSPDE)
library(excursions)
require.nowarnings(INLA)
set.seed(1)

## ------------------------------------------------------------------------
s <- seq(from = 0, to = 1, length.out = 101)

## ------------------------------------------------------------------------
fem <- rSPDE.fem1d(s)

## ------------------------------------------------------------------------
kappa <- 20
sigma <- 2
nu <- 0.8
op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                       G = fem$G, C = fem$C, d = 1, m = 1)

## ------------------------------------------------------------------------
v <- t(rSPDE.A1d(s,0.5))
c.approx <- op$Pr %*% solve(op$Q, op$Pr %*% v)
c.true <- matern.covariance(abs(s - 0.5), kappa, nu, sigma) 

## ---- fig.show='hold',fig.height = 2.5, fig.width = 7, fig.align = "center",echo=FALSE----
opar <- par(mfrow = c(1,2),mgp = c(1.3, 0.5, 0),
            mar = c(2,2,0.5,0.5) + 0.1)
plot(s, c.true, type = "l", ylab = "C(|s-0.5|)", xlab = "s",ylim=c(0,5),
     cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8)
lines(s, c.approx, col = 2)
legend("topright", bty = "n",
       legend = c("Matérn", "m=1 rSPDE"),
       col = c("black", "red"),
       lty = rep(1,2), ncol = 1,
       cex = 0.8)

plot(s, c.true - c.approx, type = "l", ylab = "Error", xlab = "s",
     cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8)
par(opar)

## ------------------------------------------------------------------------
op2 <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                       G = fem$G, C = fem$C, d = 1, m = 2)
c.approx2 <- op2$Pr %*% solve(op2$Q, op2$Pr %*% v)

s2 <- seq(from = 0, to = 1, length.out = 501)
fem2 <- rSPDE.fem1d(s2)
op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                       G = fem2$G, C = fem2$C, d = 1, m=1)
A <- rSPDE.A1d(s2,s)
v <- t(rSPDE.A1d(s2,0.5))
c.approx3 <- A%*%op$Pr %*% solve(op$Q, op$Pr %*% v)
op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                       G = fem2$G, C = fem2$C, d = 1, m=2)
c.approx4 <- A%*%op$Pr %*% solve(op$Q, op$Pr %*% v)

## ---- fig.show='hold',fig.height = 3, fig.width = 7, fig.align = "center",echo=FALSE----
opar <- par(mgp = c(1.3, 0.5, 0), mar = c(2,2,0.5,0.5) + 0.1)
plot(s, c.true - c.approx, type = "l", ylab = "Error", xlab = "s", col = 1,
     cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8)
lines(s, c.true - c.approx2, col = 2)
lines(s, c.true - c.approx3, col = 3)
lines(s, c.true - c.approx4, col = 4)
legend("bottomright", bty = "n",
       legend = c("m=1 coarse mesh", "m=2 coarse mesh", "m=1 fine mesh", "m=2 fine mesh"),
       col = c(1,2,3,4),
       lty = rep(1,2), ncol = 1,
       cex = 0.8)
par(opar)

## ------------------------------------------------------------------------
errors <- rep(0,4)
for(i in 1:4){
  op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                       G = fem2$G, C = fem2$C, d = 1, m = i)
  c.app <- A%*%op$Pr %*% solve(op$Q, op$Pr %*% v)
  errors[i] <- norm(c.true-c.app)
}
print(errors)

## ------------------------------------------------------------------------
errors2 <- rep(0,4)
for(i in 1:4){
  op <- matern.operators(kappa = kappa, sigma = sigma, nu = nu,
                       G = fem2$G, C = fem2$C, d = 1, m = i)
  c.app <- A%*%Sigma.mult(op, v)
  errors2[i] <- norm(c.true-c.app)
}
print(errors2)

## ------------------------------------------------------------------------
s <- seq(from = 0, to = 1, length.out = 501)
fem <- rSPDE.fem1d(s)
kappa <-  10*(1+2*s^2)
tau <-  0.1*(1 - 0.7*s^2)
op <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu, 
                            G = fem$G, C = fem$C, d = 1, m=1)

## ------------------------------------------------------------------------
v <- t(rSPDE.A1d(s, c(0.1,0.5,0.9)))
covs <- Sigma.mult(op, v)

## ---- fig.show='hold',fig.height = 3, fig.width = 7, fig.align = "center",echo=FALSE----
opar <- par(mgp = c(1.3, 0.5, 0), mar = c(2,2,0.5,0.5) + 0.1)
plot(s, covs[,1], type = "l", ylab = "C(s,s_i)", xlab = "s",
     cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8)
lines(s, covs[,2], col = 2)
lines(s, covs[,3], col = 3)
par(opar)

## ------------------------------------------------------------------------
c = min(kappa)^2
L = fem$G + fem$C %*% Diagonal(501, kappa^2)

## ------------------------------------------------------------------------
op <- fractional.operators(L = L, beta = (nu + 1/2)/2, C = fem$C, 
                           scale.factor = c, tau = tau, m = 1)

## ------------------------------------------------------------------------
covs2 <- Sigma.mult(op,v)
norm(covs-covs2)

## ------------------------------------------------------------------------
u <- simulate(op)

## ------------------------------------------------------------------------
n.obs <- 20
obs.loc <- runif(n = n.obs, min = 0, max = 1)
A <- rSPDE.A1d(s, obs.loc)

## ------------------------------------------------------------------------
sigma.e <- 0.3
Y <- as.vector(A %*% u + sigma.e * rnorm(n.obs))

## ---- fig.show='hold',fig.height = 4, fig.width = 5, fig.align = "center"----
A.krig <- rSPDE.A1d(s, s)
u.krig <- predict(op, A = A, Aprd = A.krig, Y = Y, sigma.e = sigma.e)

## ---- fig.show='hold',fig.height = 2.5, fig.width = 7, fig.align = "center",echo=FALSE----
opar <- par(mgp = c(1.3, 0.5, 0), mar = c(2,2,0.5,0.5) + 0.1)
plot(obs.loc, Y, ylab = "u(s)", xlab = "s", 
     ylim = c(min(c(min(u), min(Y))), max(c(max(u), max(Y)))),
          cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8)
lines(s, u)
lines(s, u.krig$mean, col = 2)
par(opar)

## ---- eval = require.nowarnings(INLA)------------------------------------
x <- seq(from = 0, to = 10, length.out = 70)
mesh <- inla.mesh.create(lattice = inla.mesh.lattice(x = x, y = x), 
                         extend = FALSE, refine = FALSE)
fem <- inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2,
                         output = list("c0", "c1", "g1"))
C <- fem$c0
G <- fem$g1

## ---- eval = require.nowarnings(INLA)------------------------------------
kappa <- 0.5
tau <- 1
nu <- 0.5
op <- spde.matern.operators(kappa = kappa, tau = tau, nu = nu, G = G, C = C, d = 2, m = 1)

## ---- eval = require.nowarnings(INLA)------------------------------------
u <- simulate(op)
n.obs <- 4000
obs.loc <- cbind(runif(n = n.obs, min = 0, max = 1), runif(n = n.obs, min = 0, max = 1))
A <- inla.spde.make.A(mesh, loc = obs.loc)
sigma.<- 0.1
Y = as.vector(A %*% u + sigma.e*rnorm(n.obs))

## ---- fig.show='hold',fig.height = 3.5, fig.width = 7, fig.align = "center",echo=FALSE, eval = require.nowarnings(INLA)----
opar <- par(mfrow = c(1,2),mgp = c(1.2, 0.5, 0), mar = c(2,2,0.5,0.5) + 0.1)
proj <- inla.mesh.projector(mesh, dims = c(70, 70))
image(inla.mesh.project(proj, field = as.vector(u)),xlab="",ylab="",
      cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8)
plot(obs.loc[,1],obs.loc[,2],cex=0.2,pch=16,xlab="",ylab="",
     cex.main = 0.8, cex.axis = 0.8, cex.lab = 0.8)
par(opar)

## ------------------------------------------------------------------------
mlik <- function(theta, Y, G, C, A) {
  return(-spde.matern.loglike(exp(theta[1]), exp(theta[2]), exp(theta[3]), exp(theta[4]),
                              Y = Y, G = G, C = C, A = A, d = 2, m=1))
}

## ---- eval = require.nowarnings(INLA)------------------------------------
theta0 = log(c(2, sqrt(var(Y)), 1,0.1*sqrt(var(Y))))
pars <- optim(theta0, mlik, Y = Y, G = G, C = C, A = A, method = "L-BFGS-B")
results <- data.frame(kappa = c(kappa, exp(pars$par[1])), 
                      tau = c(tau, exp(pars$par[2])),
                      nu = c(nu, exp(pars$par[3])),
                      sigma.e = c(sigma.e, exp(pars$par[4])),
                      row.names = c("True", "Estimate"))
print(results)

