## Functions to simulate Matérn GRF (based on Krainski et al. (2018))

# Matern correlation
cMatern <- function(h, nu, kappa) {
  ifelse(h > 0, besselK(h * kappa, nu) * (h * kappa)^nu / 
           (gamma(nu) * 2^(nu - 1)), 1)
}

# Function to sample from zero mean multivariate normal
rmvnorm0 <- function(n, cov, R = NULL){ 
  if (is.null(R))
    R <- chol(cov)
  
  return(crossprod(R, matrix(rnorm(n * ncol(R)), ncol(R))))
}

# function to simulate Matérn GRF in a grid (x,y) with parameter nu, range, sigma2u
sim_GF_Matern <- function(loc_x, nu, range, sigma2u){
  # Define locations and distance matrix
  mdist <- as.matrix(dist(loc_x[,c('x','y')]))
  
  
  # Covariance parameter scenarios
  params <- cbind(nu = rep(nu, length(range)), 
                  range = rep(range, each = length(nu)))
  # Sample error
  z <- matrix(rnorm(nrow(mdist)))
  # Compute the correlated samples
  yy <- lapply(1:nrow(params), function(j) { 
    param <- c(params[j, 1],kappa = sqrt(8 * params[j, 1]) / params[j, 2],
               params[j, 2])
    names(param)[2] <- "kappa"
    v <- cMatern(mdist, param[1], param[2])
    v <- sigma2u * v 
    # Parameter scenario and computed sample
    return(list(params = param, y = crossprod(chol(v), z)))
  })
  return(yy)
}

# ## Example
# # Define parameters
# nu <- 1
# range <- sqrt(prod(grid_dim))/(5)*4  # Range ~ 2*Scale
# sigma2u <- 1
# 
# # simulate data
# loc_x$rn <- sim_GF_Matern(loc_x, nu, range, sigma2u)[[1]]$y
# 
# ggplot(loc_x)+
#   geom_point(aes(x=x,y=y,col=rn))+
#   theme_bw()


# RFfit(RMgauss(var=NA, scale=NA),data=yy[[1]]$y,x = loc_x$x, y = loc_x$y)