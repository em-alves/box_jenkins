# Econ 112 Macroeconomic Data Analysis

# Writing Assignment 1

# Box-Jenkins Method - Time Series Model of UK Quarterly Unemployment Rate

# load data into workspace
data.unem = read.table("C:/Users/mlylv/OneDrive/Desktop/R/UK_Unem_Rate1.csv", header = T, sep = ",")

# turn data into time series
unemr = ts(data.unem[,2], start = 1960, frequency = 4)
plot(unemr, main = "UK Unemployment Rate", ylab = "Unemployment Rate", type = 'l')

# compute annualized quarterly unemployment rate
dlogunemr = 100 * diff(log(unemr))

# rename variable
u = dlogunemr
u = ts(dlogunemr, start = 1960, frequency = 4)

# Augmented Dickey Fuller test - Tests the null hypothesis that a unit root is present
# in a time series sample. Rejecting the null means series is stationary.

install.packages('tseries')
library(tseries)
adf.test(u, alternative = c("stationary", "explosive"), k = 12)


# IDENTIFICATION STAGE
# plot data, ACF, PACF

plot(u, main = "UK Quarterly Unemployment Rate", xlab = "Time", ylab = "Unemployment Rate", type = "l")

(Acf(u, main = "Autocorrelation Unemployment Rate, UK"))

(Pacf(u, main = "Partial Autocorrelation Unemployment Rate, UK"))

# ESTIMATION STAGE
library(dyn)
library(forecast)
library(sandwich)
library(lmtest)

# possible candidates models: AR(1), ARMA(1,1)

# estimate models
AR.1.MA.0.u <- arima(u, order = c(1,0,0))
coeftest(AR.1.MA.0.u, vcv = NeweyWest)
summary(AR.1.MA.0.u)
(fit <- arima(u, order = c(1,0,0)))
autoplot(fit)

AR.2.MA.0.u <- arima(u, order = c(2,0,0))
coeftest(AR.2.MA.0.u, vcv = NeweyWest)

AR.1.MA.1.u <- arima(u, order = c(1,0,1))
coeftest(AR.1.MA.1.u, vcv = NeweyWest)

AR.2.MA.1.u <- arima(u, order = c(2,0,1))
coeftest(AR.2.MA.1.u, vcv = NeweyWest)

AR.1.MA.2.u <- arima(u, order = c(1,0,2))
coeftest(AR.1.MA.2.u, vcv = NeweyWest)

plot(arima.errors(AR.1.MA.0.u))

AIC(AR.1.MA.0.u,AR.2.MA.0.u,AR.1.MA.1.u,AR.2.MA.1.u,AR.1.MA.2.u)
BIC(AR.1.MA.0.u,AR.2.MA.0.u,AR.1.MA.1.u,AR.2.MA.1.u,AR.1.MA.2.u)

# DIAGNOSTIC CHECK

# chosen model: AR(1)

# re-estimate using OLS

AR.1.MA.0.u <- dyn$lm(u ~ lag(u,-1))
coeftest(AR.1.u, vcv = NeweyWest)

AR.1.u.fitted <- fitted(AR.1.MA.0.u)
plot(AR.1.u.fitted)

plot(u, main = "Fitted AR(1) Estimates - UK UnemRate", xlab = "Time", type = "l", col = "blue")
lines(AR.1.u.fitted, col = "red")

# compute residuals

AR.1.u.residuals <- residuals(AR.1.MA.0.u)

plot(AR.1.u.residuals, main = "Residuals of AR(1) Estimates", xlab = "Time", type = "l", col = "blue")

acf(AR.1.u.residuals, lag.max = 100, main = "Autocorrelation of Residuals for AR(1) of UK Unemployment Rate")

# plot histogram of residuals and contrast with Gaussian distribution
# to check if residuals are normally distributed

x <- AR.1.u.residuals
h <- hist(x, breaks = 40, col = "red", xlab = "UK Unemployment Rate", 
        main = "Histogram with Normal Curve") 
xfit <- seq(min(x), max(x), length = 40)
yfit <- dnorm(xfit, mean = mean(x), sd = sd(x)) 
yfit <- yfit * diff(h$mids[1:2]) * length(x) 
lines(xfit, yfit, col = "blue", lwd = 2)

checkresiduals(AR.1.u.residuals)

summary(AR.1.u)
summary(AR.1.u.residuals)
summary(u)

adf.test(AR.1.u.fitted, alternative = c("stationary", "explosive"), k = 12)


qqnorm(AR.1.u.residuals); qqline(AR.1.u.residuals)


autoplot(AR.1.u.fitted)

Box.test(AR.1.u.residuals, lag=24, fitdf=4, type="Ljung")

Box.Ljung.Test(AR.1.u.residuals, lag = NULL, main = NULL)


var(AR.1.u.residuals)
