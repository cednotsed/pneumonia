require(MKpower)
require(tidyverse)

# Calculate power
rx <- function(n) rnorm(n, mean = 2.31 - 2, sd = 1.51)
ry <- function(n) rnorm(n, mean = 2.31, sd = 1.51) 

## two-sample
set.seed(66)
pwr <- sim.ssize.wilcox.test(rx = rx, ry = ry, 
                      n.min = 1, 
                      n.max = 60, 
                      step.size = 1,
                      iter = 5000,
                      alternative = "less",
                      BREAK = F)

tibble(n = pwr$n, power = pwr$emp.power) %>%
  ggplot(aes(x = n, y = power)) +
  geom_point() +
  geom_line() +
  geom_text(aes(label = power),
            vjust = 2) +
  labs(x = "Per-group sample size", y = "Empirical power", title = "Monte carlo simulation of Mann-Whitney U test")


pwr.t.test(d = 2,
           n = 5,
           sig.level = 0.05,
           type = "paired",
           alternative = "two.sided")
