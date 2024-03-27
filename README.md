# How to install
Install the package from my [Github repo](https://github.com/xl0418/SimLV) by running this command in R console:

```
library(devtools)
install_github("xl0418/SimLV")
```
# How to use

Try this:

```r
library(SimLV)
library(ggplot2)
library(reshape2)

# Simulate data
# test
n0 <- c(1, 1, 1)
r <- c(1, 1, 1)
alpha <- matrix(c(0.2, 0.1, 0.3, 0.1, 0.5, 0.3, 0.1, 0.12, 0.3), nrow = 3, byrow = TRUE)
tmax <- 100
simresult <- simlv(n0 = n0, r = r, alphaij = alpha, tmax = tmax)

# long data frame
df <- melt(simresult, id.vars=c("time"))

# plot the result
ggplot(df, aes(x = time, y = value, group = variable, color = factor(variable))) + geom_line(linewidth = 2) + theme_minimal() + xlab("Time") + ylab("Population size") + ggtitle("Simulated population dynamics") + scale_color_discrete(name = "Species")

```
