library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

T <- 143
N <- 150
T0 <- 7
D <- 100

infprofile <- read.csv("serial_interval.csv")$fit

uk_cases <- read.csv("uk_cases.csv")

C <- uk_cases[1:150,3:145]

cori_dat <- list(N = N, T = T, T0 = T0, D = D, C = C, infprofile = infprofile)

fit <- stan(file = 'cori-simple.stan', data = cori_dat)
# fit <- stan(file = 'cori-gp.stan', data = cori_dat)
print(fit)

s <- summary(fit, pars="Rt", probs=c(0.025, .5, .975))$summary
Rt <- s[,"50%"]
Rt <- t(t(Rt))
df <- data.frame(area = uk_cases[1:150,2], Rt = Rt)

write.csv(df,"Rt.csv",row.names=FALSE)

