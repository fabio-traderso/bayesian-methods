rbiTrNormGibbs <- function (initx = 2, inity = -2, rho, burnin = 100, R = 10000) 
{
  kernel = function(x, mu, rooti) {
    z = as.vector(t(rooti) %*% (x - mu))
    (z %*% z)
  }
  if (missing(rho)) {
    pandterm("Requires rho argument ")
  }
  cat("Bivariate Truncated Normal Gibbs Sampler", fill = TRUE)
  cat("rho= ", rho, fill = TRUE)
  cat("initial x,y coordinates= (", initx, ",", inity, ")", 
      fill = TRUE)
  cat("burn-in= ", burnin, " R= ", R, fill = TRUE)
  cat(" ", fill = TRUE)
  cat(" ", fill = TRUE)
  sd = (1 - rho^2)^(0.5)
  sigma = matrix(c(1, rho, rho, 1), ncol = 2)
  rooti = backsolve(chol(sigma), diag(2))
  mu = c(0, 0)
  x = seq(-3.5, 3.5, length = 100)
  y = x
  z = matrix(double(100 * 100), ncol = 100)
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      z[i, j] = kernel(c(x[i], y[j]), mu, rooti)
    }
  }
  prob = c(0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  lev = qchisq(prob, 2)
  par(mfrow = c(1, 1))
  contour(x, y, z, levels = lev, labels = prob, xlab = "theta1", 
          ylab = "theta2", drawlabels = TRUE, col = "green", labcex = 1.3, 
          lwd = 2)
  title(paste("Gibbs Truncated Sampler with Intermediate Moves: Rho =", 
              rho))
  points(initx, inity, pch = "B", cex = 1.5)
  oldx = initx
  oldy = inity
  continue = "y"
  r = 0
  draws = matrix(double(R * 2), ncol = 2)
  draws[1, ] = c(initx, inity)
  cat(" ")
  cat("Starting Gibbs Sampler ....", fill = TRUE)
  cat("(hit enter or y to display moves one-at-a-time)", fill = TRUE)
  cat("('go' to paint all moves without stopping to prompt)", 
      fill = TRUE)
  cat(" ", fill = TRUE)
  while (continue != "n" && r < R) {
    if (continue != "go") 
      continue = readline("cont?")
    newy = sd * rtrun(0,1,-Inf,0) + rho * oldx
    lines(c(oldx, oldx), c(oldy, newy), col = "magenta", 
          lwd = 1.5)
    newx = sd * rtrun(0,1,0,Inf) + rho * newy
    lines(c(oldx, newx), c(newy, newy), col = "magenta", 
          lwd = 1.5)
    oldy = newy
    oldx = newx
    r = r + 1
    draws[r, ] = c(newx, newy)
  }
  continue = readline("Show Comparison to iid Sampler?")
  if (continue != "n" & continue != "No" & continue != "no") {
    par(mfrow = c(1, 2))
    contour(x, y, z, levels = lev, xlab = "theta1", ylab = "theta2", 
            drawlabels = TRUE, labels = prob, labcex = 1.1, col = "green", 
            lwd = 2)
    title(paste("Gibbs Draws: Rho =", rho))
    points(draws[(burnin + 1):R, ], pch = 20, col = "magenta", 
           cex = 0.7)
    idraws = t(chol(sigma)) %*% matrix(rtrun(2 * (R - burnin)), 
                                       nrow = 2)
    idraws = t(idraws)
    contour(x, y, z, levels = lev, xlab = "theta1", ylab = "theta2", 
            drawlabels = TRUE, labels = prob, labcex = 1.1, col = "green", 
            lwd = 2)
    title(paste("IID draws: Rho =", rho))
    points(idraws, pch = 20, col = "magenta", cex = 0.7)
  }
  attributes(draws)$class = c("bayesm.mat", "mcmc")
  attributes(draws)$mcpar = c(1, R, 1)
  return(draws)
}
