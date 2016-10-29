x <- rnorm(500)
u <- rnorm(500, sd = sqrt(abs(x)))
y <- 0.01 * x + u

m_white <- lm(y ~ x)
squared_uhat <- resid(m_white)**2 ; yhat <- predict(m_white)
m_white_stage2 <- lm(squared_uhat ~ yhat + I(yhat^2))
R_squared <- summary(m_white_stage2)$r.squared

n <- length(y) ; k <- 2 # k = 1 because we have 2 regressors, yhat and yhat squared
(F_statistic <- (R_squared / k) / ((1 - R_squared) / (n - k - 1)))

1- pf(F_statistic, df1 = 2, df2 = n - 3)

summary(m_white_stage2)$fstatistic["value"]
bptest(m_white, ~ yhat + I(yhat^2))$statistic
