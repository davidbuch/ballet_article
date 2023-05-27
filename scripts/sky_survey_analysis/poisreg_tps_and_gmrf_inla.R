library(rstanarm)
data <- read.csv("data/intermediate_data/edsgc.csv")

res <- 2^4
breaks_x <- with(data, seq(min(x), max(x), length.out = res + 1))
breaks_y <- with(data, seq(min(y), max(y), length.out = res + 1))


midpoints_x <- (breaks_x[1:res] + breaks_x[2:(res + 1)]) / 2
midpoints_y <- (breaks_y[1:res] + breaks_y[2:(res + 1)]) / 2


binned_x <- midpoints_x[with(data, cut(x, breaks = breaks_x, labels = FALSE))]
binned_y <- midpoints_y[with(data, cut(y, breaks = breaks_y, labels = FALSE))]

y <- as.integer(table(binned_x, binned_y))
X <- expand.grid(midpoints_x, midpoints_y)

# Create Regression Dataframe
df <- 5
bs <- splines::bs(1:res, df)
Z <- apply(expand.grid(1:df,1:df), 1, \(ids) c(outer(bs[,ids[1]], bs[,ids[2]])))
rdata <- data.frame(Z)
rdata$y <- y

stan_poisreg <- stan_glm(y ~ . - 1, 
                         family = poisson, 
                         data = rdata, chains = 2, iter = 1000)

stan_samps <- rstan::extract(stan_poisreg$stanfit)
pred_samps <-  stan_samps$beta %*% t(Z)
f_pe <- colMeans(pred_samps)

bench::bench_time(poisreg <- glm(y ~ Z, family = poisson))
f_pe <- poisreg$fitted.values

N <- nrow(X)
Mesh <- inla.mesh.2d(as.matrix(X), max.edge = c(20, 40))
HostsA <- inla.spde.make.A(Mesh, loc = as.matrix(X)) # Making A matrix
Hosts.spde = inla.spde2.pcmatern(mesh = Mesh, prior.range = c(10, 0.5), prior.sigma = c(.5, .5)) # Making SPDE
w.Host <- inla.spde.make.index('w', n.spde = Hosts.spde$n.spde) # making the w
StackHost <- inla.stack(
  data = list(y = y), # specify the response variable
  A = list(1, HostsA), # Vector of Multiplication factors for random and fixed effects 
  effects = list(
    Intercept = rep(1, N),
    w = w.Host)) 


poisr_inla <- inla(y ~ -1 + Intercept +
                     f(w, model = Hosts.spde), 
                   family = "poisson", 
                   data = inla.stack.data(StackHost), 
                   control.compute = list(dic = TRUE), 
                   control.predictor = list(A = inla.stack.A(StackHost)))

ggplot() + 
  geom_contour_filled(aes(x = Var1, y = Var2, z = poisr_inla$summary.fitted.values[1:4096,1]), data = X) + 
  geom_point(aes(x = x, y = y), alpha = 0.2, size = 0.1, data = data) + 
  guides(fill = 'none')
