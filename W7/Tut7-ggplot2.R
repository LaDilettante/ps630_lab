### Tutorial 7: ggplot2
gdp_vals=seq(1000,10000,by=1000)

nd = data.frame(l1gdp_pc=rep(gdp_vals, 3),l1fdi=c(rep(0.0269332,10),rep(0.6382053,10),rep(1.9903931,10)))
summary(intmodel)

nd_predict <- predict(intmodel, type = "response", newdata = nd, se.fit = T)

nd_predict=cbind(nd, nd_predict)

# Confidence Intervals
nd_predict$low  <- nd_predict$fit - 2* nd_predict$se.fit 
nd_predict$high  <- nd_predict$fit + 2* nd_predict$se.fit 

# Plot this in ggplot2
library(ggplot2)
ggplot(data = nd_predict, aes ( x = l1gdp_pc, y = fit, color = factor(l1fdi))) +
  scale_colour_discrete( name = "FDI", 
                         labels = c("1st Quartile", "2nd Quartile", "3rd Quartile")) +
  geom_line() +
  labs( x = "GDP per capita",
        y = "Marginal effect") + 
  geom_ribbon(data = nd_predict , aes(ymin = low , ymax = high), alpha = 0.2)