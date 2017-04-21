library(ggplot2)
setwd("K:/FORSK-Projekt/Projekter/Scientific Projects/141_Jorne_phd/Scripts/RobustStatistics")

set.seed(516516)

x <- seq(0, 1, 0.1)
y <- rnorm(length(x), mean = x, sd = 0.1)
outlier <- c(0, 1, rep(0, length(x) - 2))

lm1 <- lm(y ~ x)

plot1 <- ggplot(data.frame(Covariate = x, Outcome = y, outlier), aes(Covariate, Outcome, col = outlier)) + 
  geom_point(size = 3) +  
  stat_smooth(method = "lm", se = F) +
  ggtitle("Original data") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave("./Figures/linearNoOutlier.pdf", plot1, width =  5, height = 5)


xOut <- c(0, 1, seq(0.2, 1, 0.1))  
plot2 <- ggplot(data.frame(Covariate = xOut, Outcome = y, outlier), aes(Covariate, Outcome, col = outlier)) + 
  geom_point(size = 3) + 
  stat_smooth(method = "lm", se = F) + 
  ggtitle("Contaminated data") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave("./Figures/linearOutlier.pdf", plot2, width =  5, height = 5)