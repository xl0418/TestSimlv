library(ggplot2)
library(reshape2)

simrc <-
  function(n0,
           r0,
           mu,
           k,
           v,
           m,
           rin,
           tmax,
           mean_d_n = 0,
           sd_d_n = 0.1,
           ht = 1,
           sto_mode = "multiplicative",
           growth_mode = "linear") {
    no_species <- length(n0)
    no_resources <- length(r0)
    nt <- data.frame(matrix(NA, nrow = tmax, ncol = no_species))
    rt <- data.frame(matrix(NA, nrow = tmax, ncol = no_resources))
    nt[1,] <- n0
    rt[1,] <- r0
    # simulate the lotka-volterra equations
    for (tt in 1:(tmax - 1)) {
      temp_nt <- nt[tt,]
      temp_rt <- rt[tt,]
      # if(tt == 100){
      #   browser()
      # }
      if (sto_mode == "additive") {
        for (tti in 1:ht) {
          delta_nt <- rep(0, no_species)
          for (ii in 1:no_species) {
            if(growth_mode == "linear"){
              # additive error
              delta_nt[ii] <- (mu[ii] * temp_rt - m) / ht
            }else{
              # additive error
              delta_nt[ii] <- (mu[ii] * temp_rt / (k[ii] + temp_rt) - m) / ht
            }

            temp_nt[ii] <- max(temp_nt[ii] + rnorm(1, mean_d_n, sd_d_n) / ht, 0 )

          }
          temp_rt <-
            temp_rt + ((rin - temp_rt) - sum(v * temp_nt * temp_rt)) / ht
        }
      } else{
        for (tti in 1:ht) {
          delta_nt <- rep(0, no_species)
          for (ii in 1:no_species) {
            if(growth_mode == "lienar"){
              # multiplicative error
              delta_nt[ii] <- (mu[ii] * temp_rt - m)/ ht
            }else{
              # multiplicative error
              delta_nt[ii] <-(mu[ii] * temp_rt / (k[ii] + temp_rt) - m) / ht
            }
            temp_nt[ii] <-
              temp_nt[ii] * exp(rnorm(1, delta_nt[ii], sd_d_n))

          }
          temp_rt <-
            temp_rt + ((rin - temp_rt) - sum(v * temp_nt * temp_rt)) / ht
        }
      }

      # max of temp_nt and 0
      temp_nt[temp_nt < 0] <- 0
      temp_rt[temp_rt < 0] <- 0

      nt[tt + 1,] <- temp_nt
      rt[tt + 1,] <- temp_rt
    }
    df <- data.frame(time = 1:tmax, nt, rt)
    colnames(df) <-
      c("time",
        paste0("n", 1:no_species),
        paste0("r", 1:no_resources))
    return(df)

  }

# Simulate data
# test
n0 <- c(0.1, 0.1)
r <- c(1)
mu <- c(1, 0.5)
k <- c(1, 1)
v <- 1.3 * mu
m <- 0.1
rin <- 1
tmax <- 200
sd_d_n <- 0.1

growth_mode <- "mm"

simresult_multiplicative <-
  simrc(
    n0 = n0,
    r0 = r,
    mu = mu,
    k = k,
    v = v,
    m = m,
    rin = rin,
    tmax = tmax,
    sto_mode = "multiplicative",
    growth_mode = growth_mode,
    sd_d_n = sd_d_n,
    ht = 10
  )
simresult_multiplicative$model <- "multiplicative"

cov_n1r <-
  cov(simresult_multiplicative$n1, simresult_multiplicative$r1)

simresult_additive <-
  simrc(
    n0 = n0,
    r0 = r,
    mu = mu,
    k = k,
    v = v,
    m = m,
    rin = rin,
    tmax = tmax,
    sto_mode = "additive",
    growth_mode = growth_mode,
    sd_d_n = sd_d_n,
    ht = 10
  )
simresult_additive$model <- "additive"

simresult <- rbind(simresult_multiplicative, simresult_additive)
# long data frame
df <- melt(simresult, id.vars = c("time", "model"))

r1_star <- m / mu[1]
r2_star <- m / mu[2]

# plot the result
p <-
  ggplot(df, aes(
    x = time,
    y = value,
    group = variable,
    color = factor(variable)
  )) + geom_line(linewidth = 2) +
  theme_minimal() + xlab("Time") + ylab("Population size") + ggtitle("Simulated population dynamics") + scale_color_discrete(name = "Species") +
  facet_wrap( ~ model) + geom_hline(yintercept = r1_star,
                                    linetype = "dashed",
                                    color = "red") + geom_hline(yintercept = r2_star,
                                                                linetype = "dashed",
                                                                color = "green") + theme(legend.position = "top")

print(p)



#####################  periodic growth rates ############################


growth_rate <- function(mu, temp, temp_opt, sigma){
  return(mu * exp( -(temp_opt - temp)^2 / sigma^2))
}


Sim_rc <- function(tt, mu1, mu2,
                   v1, v2,
                   r1in,
                   N10, N20, R10,
                   m,
                   temp, temp_opt1, temp_opt2, sigma, tau
) {
  N1 <- c(N10)
  N2 <- c(N20)
  R1 <- c(R10)
  mu1_vec <- c(mu1)
  mu2_vec <- c(mu2)
  for(ti in 1:tt){
    temp_N1 <- N1[ti]
    temp_N2 <- N2[ti]
    temp_R1 <- R1[ti]
    for(tii in 1:100){
      temp <- sin(2 * pi * (ti + tii / 100) / tau)
      mu1_g <- growth_rate(mu1, temp, temp_opt1, sigma)
      mu2_g <- growth_rate(mu2, temp, temp_opt2, sigma)
      dN1 <- (mu1_g * temp_R1 - m) * temp_N1 / 100
      dN2 <- (mu2_g * temp_R1 - m) * temp_N2 / 100
      dR1 <- ((r1in - temp_R1) - v1 * mu1_g * temp_N1* temp_R1  - v2 * mu2_g * temp_N2* temp_R1 ) / 100

      temp_N1 <- max(temp_N1 + dN1, 0)
      temp_N2 <- max(temp_N2 + dN2, 0)
      temp_R1 <- max(temp_R1 + dR1, 0)

    }
    mu1_vec <- c(mu1_vec, mu1_g)
    mu2_vec <- c(mu2_vec, mu2_g)

    N1 <- c(N1, temp_N1)
    N2 <- c(N2, temp_N2)

    R1 <- c(R1, temp_R1)
  }
  return(data.frame(N1 = N1, N2 = N2, R1 = R1, mu1 = mu1_vec, mu2 = mu2_vec, cc = 1:(tt+1)))

}

### sim ###

mu1 <- 1
mu2 <- 1
v1 <- 1.3
v2 <- 1.3
r1in <- 10
N10 <- 0.1
N20 <- 0.1
R10 <- 1
m <- 0.1
temp <- 0
temp_opt1 <- -0.25
temp_opt2 <- 0.3
sigma <- 1
tau <- 20

simresult <-
  Sim_rc(2000, mu1, mu2, v1, v2, r1in, N10, N20, R10, m, temp, temp_opt1, temp_opt2, sigma, tau)

plot_df <- simresult[, c("cc", "N1", "N2", "R1")]

### plot ###
p <-
  ggplot(melt(plot_df, id.vars = "cc"), aes(x = cc, y = value, group = variable, color = variable)) +
  geom_line(linewidth = 2) +
  theme_minimal() +
  xlab("Time") +
  ylab("Population size") +
  ggtitle("Simulated population dynamics") +
  scale_color_discrete(name = "Species") +
  theme(legend.position = "top")

p


#####

storage1 <- cov(simresult$mu1, simresult$R1)
storage2 <- cov(simresult$mu2, simresult$R1)

R1_star <- (m - storage1) / mean(simresult$mu1)
R2_star <- (m - storage2) / mean(simresult$mu2)
R1_mean <- mean(simresult$R1)


