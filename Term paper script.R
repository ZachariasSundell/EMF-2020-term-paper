# Download packages:
library("quantmod")
library('psych')
library('TSA')
library('tsoutliers')
library('urca')
library(aTSA)
library(forecast)
library('FinTS')
library('rugarch')
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(extrafont)
library(zoo)
library('rmgarch')
library('WeightedPortTest')
library(writexl)
library('tidyquant')


#library('plyr')
#library('matrixcalc')
#library(data.table)
#library(plyr)
#library(tsbox)
#library(TSstudio)
#library(dygraphs)
#library(plotly)
#library('MTS')
#library('mgarchBEKK')
#library(vars)

### Data manipulation ####
startDate <- as.Date("2016-01-01")
endDate <- as.Date("2020-12-11")
getSymbols("^GSPC", from = startDate, to=endDate) #S&P 500
getSymbols("^N100", from = startDate, to=endDate) #Euronext 100
getSymbols("^N225", from = startDate, to=endDate) # Nikkei 225
getSymbols("000001.SS", from = startDate, to=endDate) # SSEC
SSEC <- `000001.SS`

par(mfrow = c(2,2))
plot.xts(GSPC$'GSPC.Adjusted')
plot.xts(N100$'N100.Adjusted')
plot.xts(N225$'N225.Adjusted')
plot.xts(SSEC$'000001.SS.Adjusted')
#SSEC <- xts(SSE[,-1], order.by = SSE$Date)
# Log returns
logGSPC <- diff(log(GSPC$GSPC.Adjusted)*100)
logGSPC
logN100 <- diff(log(N100$N100.Adjusted)*100)
logN100
logN225 <- diff(log(N225$N225.Adjusted)*100)
logN225
logSSEC <- diff(log(SSEC$`000001.SS.Adjusted`)*100)
logSSEC
#Merge
logTS <- merge(logGSPC,logN100,logN225,logSSEC)
names(logTS) <- c("GSPC","N100","N225","SSEC")
logTS <- logTS[date(logTS) >= "2016-03-01" & date(logTS) <= "2020-12-10",]
# Missing variable transformation
#logTS <- na.omit(logTS)
logTS[is.na(logTS)] <- 0
length(logTS$GSPC)
#### Visualization of volatility ####
# Visualization
par(mfrow= c(2,2))
plot.xts(logTS$GSPC, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "GSPC")
plot.xts(logTS$N100, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "N100")
plot.xts(logTS$N225, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "N225")
plot.xts(logTS$SSEC, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "SSEC")
par(mfrow= c(1,1))

### Dividing data into 2 periods ####
period1 <- logTS[date(logTS) >= "2016-03-01" & date(logTS) < "2020-01-01",]
length(period1$GSPC)
period2 <- logTS[date(logTS) >= "2020-01-01",]
length(period2$GSPC)
### Descriptive statistics: Period 1 ####
# Descriptive Statistics:
nrow(period1)
desc_stat1 <- describe(period1, quant = c(0.25,0.75))
round(desc_stat1$mean, 5)
round(desc_stat1$median, 5)
round(desc_stat1$sd, 5)
round(desc_stat1$skew, 4)
round(desc_stat1$kurtosis, 4)
# Visualization
par(mfrow= c(2,2))
plot.xts(period1$GSPC, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "GSPC")
plot.xts(period1$N100, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "N100")
plot.xts(period1$N225, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "N225")
plot.xts(period1$SSEC, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "SSEC")
par(mfrow= c(1,1))

### MA,AR, ARMA or ARIMA model estimation Period 1 ####
# Visualization
# ACF
par(mfrow=c(2,2))
acf(period1$GSPC, lag.max = 30, main = "GSPC")
acf(period1$N100, lag.max = 30, main = "N100")
acf(period1$N225, lag.max = 30, main = "N225")
acf(period1$SSEC, lag.max = 30, main = "SSEC")
par(mfrow=c(1,1))
# PACF
par(mfrow=c(2,2))
pacf(period1$GSPC, lag.max = 30, main = "GSPC")
pacf(period1$N100, lag.max = 30, main = "N100")
pacf(period1$N225, lag.max = 30, main = "N225")
pacf(period1$SSEC, lag.max = 30, main = "SSEC")
par(mfrow=c(1,1))
# Jarque Bera tests for normality 
JB1_GSPC <- JarqueBera.test(period1$GSPC)
JB1_GSPC
round(JB1_GSPC[[1]][["statistic"]][["X-squared"]], 4)
JB1_N100 <- JarqueBera.test(period1$N100)
JB1_N100
round(JB1_N100[[1]][["statistic"]][["X-squared"]], 4)
JB1_N225 <- JarqueBera.test(period1$N225)
JB1_N225
round(JB1_N225[[1]][["statistic"]][["X-squared"]], 4)
JB1_SSEC <- JarqueBera.test(period1$SSEC)
JB1_SSEC
round(JB1_SSEC[[1]][["statistic"]][["X-squared"]], 4)
# None of the series are normally distributed
# STATIONARITY: 
# ADF unit root test for stationarity of time-series 
adf1GSPC <- ur.df(period1$GSPC, type = "none", lags = trunc((length(period1$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf1GSPC)
adf1N100 <- ur.df(period1$N100, type = "none", lags = trunc((length(period1$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf1N100)
adf1N225 <- ur.df(period1$N225, type = "none", lags = trunc((length(period1$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf1N225)
adf1SSEC <- ur.df(period1$SSEC, type = "none", lags = trunc((length(period1$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf1SSEC)
# KPSS test as complement to ADF test
kpss1GSPC <- kpss.test(period1$GSPC, lag.short = TRUE, output = TRUE)
kpss1GSPC
kpss1N100 <- kpss.test(period1$N100, lag.short = TRUE, output = TRUE)
kpss1N100
kpss1N225 <- kpss.test(period1$N225, lag.short = TRUE, output = TRUE)
kpss1N225
kpss1SSEC <- kpss.test(period1$SSEC, lag.short = TRUE, output = TRUE)
kpss1SSEC
# Ljung-Box test for autocorrelation. 
LB1_GSPC <- Box.test(period1$GSPC, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
LB1_GSPC
round(LB1_GSPC[["statistic"]][["X-squared"]], 4)
LB1_N100 <- Box.test(period1$N100, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
LB1_N100
round(LB1_N100[["statistic"]][["X-squared"]], 4)
LB1_N225 <- Box.test(period1$N225, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
LB1_N225
round(LB1_N225[["statistic"]][["X-squared"]], 4)
LB1_SSEC <- Box.test(period1$SSEC, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
LB1_SSEC
round(LB1_SSEC[["statistic"]][["X-squared"]], 4)
# Find best fitting ARMA model and fitting it with arima()
GSPC1_aic <- auto.arima(period1$GSPC,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
GSPC1_aic
arma_GSPC1 <- arima(period1$GSPC, order = c(1,0,0))
arma_GSPC1
N1001_aic <- auto.arima(period1$N100,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
N1001_aic
arma_N1001 <- arima(period1$N100, order = c(1,0,0))
arma_N1001
N2251_aic <- auto.arima(period1$N225,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
N2251_aic
arma_N2251 <- arima(period1$N225, order = c(1,0,0))
arma_N2251
SSEC1_aic <- auto.arima(period1$SSEC,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
SSEC1_aic
arma_SSEC1 <- arima(period1$SSEC, order = c(1,0,0))
arma_SSEC1
# LM-ARCH test for conditional heteroskedacticity 
ARCH1_GSPC <- ArchTest(period1$GSPC,lags = 5)
ARCH1_GSPC
ARCH1_N100 <- ArchTest(period1$N100,lags = 5)
ARCH1_N100
ARCH1_N225 <- ArchTest(period1$N225,lags = 5)
ARCH1_N225
ARCH1_SSEC <- ArchTest(period1$SSEC,lags = 5)
ARCH1_SSEC
#arch.test(arma_GSPC1, output = T)
#arch.test(arma_N1001, output = T)
#arch.test(arma_N2251, output = T)
#arch.test(arma_SSEC1, output = T)

# ARCH-GARCH estimations: Period 1 
# GSPC ####
ar1_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),variance.model = list(garchOrder = c(1,1), model = "eGARCH"), distribution.model = "sstd")
argarch1_GSPC <- ugarchfit(spec = ar1_garch11_spec, data = period1$GSPC)
argarch1_GSPC
# Specs #
df <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df %>% add_row(
  type_garch = argarch1_GSPC@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_GSPC@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_GSPC)[1],
  BIC = infocriteria(argarch1_GSPC)[2],
  LogLikelihood = argarch1_GSPC@fit[["LLH"]]
) -> df
df %>%
  mutate(model_name = paste0("AR(1,0)-",type_garch,"(1,1)")) -> mydf
df_long <- tidyr::pivot_longer(data = mydf %>%
                                 select(model_name,
                                        type_garch,
                                        type_dist,
                                        AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long$model_name)
# Visualization
p1 <- ggplot(df_long %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (GSPC)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p1
# N100 ####
ar4_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),variance.model = list(garchOrder = c(1,1), model = "eGARCH"), distribution.model = "std")
argarch1_N100 <- ugarchfit(spec = ar4_garch11_spec, data = period1$N100)
argarch1_N100
# Specs #
df1 <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df1 %>% add_row(
  type_garch = argarch1_N100@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_N100@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_N100)[1],
  BIC = infocriteria(argarch1_N100)[2],
  LogLikelihood = argarch1_N100@fit[["LLH"]]
) -> df1
df1 %>%
  mutate(model_name = paste0("AR(4,0)-",type_garch,"(1,1)")) -> mydf1
df_long1 <- tidyr::pivot_longer(data = mydf1 %>%
                                 select(model_name,
                                        type_garch,
                                        type_dist,
                                        AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long1$model_name)
# Visualization
p2 <- ggplot(df_long1 %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (N100)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p2
# N225 ####
ar0_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),variance.model = list(garchOrder = c(1,1), model = "eGARCH"), distribution.model = "sstd")
argarch1_N225 <- ugarchfit(spec = ar0_garch11_spec, data = period1$N225)
argarch1_N225
# Specs #
df2 <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df2 %>% add_row(
  type_garch = argarch1_N225@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_N225@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_N225)[1],
  BIC = infocriteria(argarch1_N225)[2],
  LogLikelihood = argarch1_N225@fit[["LLH"]]
) -> df2
df2 %>%
  mutate(model_name = paste0("AR(0,0)-",type_garch,"(1,1)")) -> mydf2
df_long2 <- tidyr::pivot_longer(data = mydf2 %>%
                                  select(model_name,
                                         type_garch,
                                         type_dist,
                                         AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long2$model_name)
# Visualization
p3 <- ggplot(df_long2 %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (N225)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p3
# SSEC ####
ar5_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),variance.model = list(garchOrder = c(1,1), model = "sGARCH"), distribution.model = "std")
argarch1_SSEC <- ugarchfit(spec = ar5_garch11_spec, data = period1$SSEC)
argarch1_SSEC
# Specs #
df3_5 <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df3_5 %>% add_row(
  type_garch = argarch1_SSEC@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_SSEC@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_SSEC)[1],
  BIC = infocriteria(argarch1_SSEC)[2],
  LogLikelihood = argarch1_SSEC@fit[["LLH"]]
) -> df3_5
df3_1 %>%
  mutate(model_name = paste0("AR(0,0)-",type_garch,"(1,1)")) -> mydf3_1
df_long3_1 <- tidyr::pivot_longer(data = mydf3_1 %>%
                                  select(model_name,
                                         type_garch,
                                         type_dist,
                                         AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long3_5$model_name)
merged <- do.call("rbind", list(df_long3,df_long3_1,df_long3_2, df_long3_3, df_long3_4, df_long3_5))
merged %>%
  arrange(model_name) -> merged
# Visualization
p4 <- ggplot(merged %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (SSEC)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p4

### DCC-GARCH estimation: Period 1 ####
# Combining univariate GARCH spec with multispec
mspec <- multispec(c(ar1_garch11_spec,ar4_garch11_spec,ar0_garch11_spec,ar5_garch11_spec))
# DCC specification - GARCH(1,1) for conditional correlations
dccgarch11_spec1 <- dccspec(uspec = mspec,lag =1, lag.criterion = "AIC", lag.max= 20, dccOrder = c(1,1), distribution = "mvt")
dccgarch11_spec1
# Fitting DCC model
dccfit1 <- dccfit(dccgarch11_spec1, data = period1, fit.control = list(eval.se=TRUE))
dccfit1
# Extracting time varying variance-covariance and correlation matrices 
cov1.garch11 <- rcov(dccfit1)
cov1.garch11
cor1.garch11 <- rcor(dccfit1)
cor1.garch11
# Visualization
# Time-varying correlations
par(mfrow=c(3,2))
plot(as.xts(cor1.garch11[1,2,]),main="GSPC and N100", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor1.garch11[1,3,]),main="GSPC and N225", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor1.garch11[1,4,]),main="GSPC and SSEC", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor1.garch11[2,3,]),main="N100 and N225", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor1.garch11[2,4,]),main="N100 and SSEC", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor1.garch11[3,4,]),main="N225 and SSEC", yaxis.right = FALSE, ylim = c(-0.1,1))
par(mfrow=c(1,1))
# Time-varying covariances
par(mfrow=c(3,2))
plot(as.xts(cov1.garch11[1,2,]),main="GSPC and N100", yaxis.right = FALSE, ylim = c(0,7))
plot(as.xts(cov1.garch11[1,3,]),main="GSPC and N225", yaxis.right = FALSE, ylim = c(0,7))
plot(as.xts(cov1.garch11[1,4,]),main="GSPC and SSEC", yaxis.right = FALSE, ylim = c(0,7))
plot(as.xts(cov1.garch11[2,3,]),main="N100 and N225", yaxis.right = FALSE, ylim = c(0,7))
plot(as.xts(cov1.garch11[2,4,]),main="N100 and SSEC", yaxis.right = FALSE, ylim = c(0,7))
plot(as.xts(cov1.garch11[3,4,]),main="N225 and SSEC", yaxis.right = FALSE, ylim = c(0,7))
par(mfrow=c(1,1))
# DCC model diagnostics ####
# Li-Mak test for presence of autocorrelation in the residuals 
dccLMtest_GSPC <- Weighted.LM.test(dccfit1@model$residuals[,1], dccfit1@model$sigma[,1]^2, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
dccLMtest_GSPC
dccLMtest_N100 <- Weighted.LM.test(dccfit1@model$residuals[,2], dccfit1@model$sigma[,2]^2, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
dccLMtest_N100
dccLMtest_N225 <- Weighted.LM.test(dccfit1@model$residuals[,3], dccfit1@model$sigma[,3]^2, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
dccLMtest_N225
dccLMtest_SSEC <- Weighted.LM.test(dccfit1@model$residuals[,4], dccfit1@model$sigma[,4]^2, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
dccLMtest_SSEC
# LM-ARCH test on standardized residuals (=all good)
DCCARCH1_GSPC <- ArchTest(dccfit1@mfit$stdresid[,1],lags = 5)
DCCARCH1_GSPC
DCCARCH1_N100 <- ArchTest(dccfit1@mfit$stdresid[,2],lags = 5)
DCCARCH1_N100
DCCARCH1_N225 <- ArchTest(dccfit1@mfit$stdresid[,3],lags = 5)
DCCARCH1_N225
DCCARCH1_SSEC <- ArchTest(dccfit1@mfit$stdresid[,4],lags = 5)
DCCARCH1_SSEC
# Ljung-Box test for autocorrelation in squared standardized residuals. 
dccLB1_GSPC <- Box.test(dccfit1@mfit$stdresid[,1]^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
dccLB1_GSPC
round(dccLB1_GSPC[["statistic"]][["X-squared"]], 4)
dccLB1_N100 <- Box.test(dccfit1@mfit$stdresid[,2]^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
dccLB1_N100
round(dccLB1_N100[["statistic"]][["X-squared"]], 4)
dccLB1_N225 <- Box.test(dccfit1@mfit$stdresid[,3]^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
dccLB1_N225
round(dccLB1_N225[["statistic"]][["X-squared"]], 4)
dccLB1_SSEC <- Box.test(dccfit1@mfit$stdresid[,4]^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
dccLB1_SSEC
round(dccLB1_SSEC[["statistic"]][["X-squared"]], 4)
#Residuals visualized
par(mfrow= c(2,2))
plot(dccfit1@mfit$stdresid[,1],type= "l", main = "GSPC", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit1@mfit$stdresid[,1]), b=0, col = "blue")
plot(dccfit1@mfit$stdresid[,2],type= "l", main = "N100", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit1@mfit$stdresid[,2]), b=0, col = "blue")
plot(dccfit1@mfit$stdresid[,3],type= "l", main = "N225", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit1@mfit$stdresid[,3]), b=0, col = "blue")
plot(dccfit1@mfit$stdresid[,4],type= "l", main = "SSEC", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit1@mfit$stdresid[,2]), b=0, col = "blue")
par(mfrow= c(1,1))
### BEKK-GARCH estimation: Period 1 ####
# Convert period 1 to data frame and export to excel for use in Oxmetrics
#tmp <- "C:\\Users\\Sundell\\Documents\\EMF 2020\\Term Paper\\dataframe3"
#write.zoo(period1,sep=",",file=tmp)
period1_df <- as.data.frame(period1)
write_xlsx(period1_df,"C:\\Users\\Sundell\\Documents\\EMF 2020\\dataframe2.xlsx")
#test2 <- period1_df[,c(1,3)]
#test3 <- period1_df[,c(1,4)]
#test4 <- period1_df[,c(2,3)]
#test5 <- period1_df[,c(2,4)]
#test6 <- period1_df[,c(3,4)]
#Bekk1GSPC_N100 <- BEKK11(test1, include.mean = T, cond.dist = "normal", ini.estimates = NULL)
#Bekk1GSPC_N225 <- BEKK11(test2, include.mean = T, cond.dist = "normal", ini.estimates = NULL)
#Bekk1GSPC_SSEC <- BEKK11(test3, include.mean = T, cond.dist = "normal", ini.estimates = NULL)
#Bekk1N100_N225 <- BEKK11(test4, include.mean = T, cond.dist = "normal", ini.estimates = NULL)
#Bekk1N100_SSEC <- BEKK11(test5, include.mean = T, cond.dist = "normal", ini.estimates = NULL)
#Bekk1N225_SSEC <- BEKK11(test6, include.mean = T, cond.dist = "normal", ini.estimates = NULL)

# BEKK GARCH (1,1) Model Diagnostics: Period 1 (dataframe1 extracted from Oxmetrics)
# Li-Mak test for presence of autocorrelation in the residuals 
BEKKLMtest_GSPC <- Weighted.LM.test(dataframe1$Res_GSPC, dataframe1$CondVar_GSPC, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKKLMtest_GSPC
BEKKLMtest_N100 <- Weighted.LM.test(dataframe1$Res_N100, dataframe1$CondVar_N100, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKKLMtest_N100
BEKKLMtest_N225 <- Weighted.LM.test(dataframe1$Res_N225, dataframe1$CondVar_N225, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKKLMtest_N225
BEKKLMtest_SSEC <- Weighted.LM.test(dataframe1$Res_SSEC, dataframe1$CondVar_SSEC, lag = trunc((length(period1)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKKLMtest_SSEC
# LM-ARCH test on standardized residuals 
BEKKARCH1_GSPC <- ArchTest(dataframe1$StdRes_GSPC,lags = 5)
BEKKARCH1_GSPC
BEKKARCH1_N100 <- ArchTest(dataframe1$StdRes_N100,lags = 5)
BEKKARCH1_N100
BEKKARCH1_N225 <- ArchTest(dataframe1$StdRes_N225,lags = 5)
BEKKARCH1_N225
BEKKARCH1_SSEC <- ArchTest(dataframe1$StdRes_SSEC,lags = 5)
BEKKARCH1_SSEC
# Ljung-Box test for autocorrelation in squared standardized residuals. 
BEKKLB1_GSPC <- Box.test(dataframe1$StdRes_GSPC^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
BEKKLB1_GSPC
round(BEKKLB1_GSPC[["statistic"]][["X-squared"]], 4)
BEKKLB1_N100 <- Box.test(dataframe1$StdRes_N100^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
BEKKLB1_N100
round(BEKKLB1_N100[["statistic"]][["X-squared"]], 4)
BEKKLB1_N225 <- Box.test(dataframe1$StdRes_N225^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
BEKKLB1_N225
round(BEKKLB1_N225[["statistic"]][["X-squared"]], 4)
BEKKLB1_SSEC <- Box.test(dataframe1$StdRes_SSEC^2, lag=sqrt(length(period1$GSPC)), type = "Ljung-Box")
BEKKLB1_SSEC
round(BEKKLB1_SSEC[["statistic"]][["X-squared"]], 4)

#Residuals visualized
par(mfrow= c(2,2))
plot(dataframe1$StdRes_GSPC,type= "l", main = "GSPC", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe1$StdRes_GSPC), b=0, col = "blue")
plot(dataframe1$StdRes_N100,type= "l", main = "N100", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe1$StdRes_N100), b=0, col = "blue")
plot(dataframe1$StdRes_N225,type= "l", main = "N225", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe1$StdRes_N225), b=0, col = "blue")
plot(dataframe1$StdRes_SSEC,type= "l", main = "SSEC", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe1$StdRes_SSEC), b=0, col = "blue")
par(mfrow= c(1,1))

# Period 2:
### Descriptive statistics: Period 2 ####
# Descriptive Statistics:
nrow(period2)
desc_stat2 <- describe(period2, quant = c(0.25,0.75))
round(desc_stat2$mean, 5)
round(desc_stat2$median, 5)
round(desc_stat2$sd, 5)
round(desc_stat2$skew, 4)
round(desc_stat2$kurtosis, 4)
# Visualization
par(mfrow= c(2,2))
plot.xts(period2$GSPC, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "GSPC")
plot.xts(period2$N100, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "N100")
plot.xts(period2$N225, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "N225")
plot.xts(period2$SSEC, col = "blue", ylim = c(-15,15), yaxis.right = FALSE, main = "SSEC")
par(mfrow= c(1,1))

### MA,AR, ARMA or ARIMA model estimation Period 2 ####
# Visualization
# ACF
par(mfrow=c(2,2))
acf(period2$GSPC, lag.max = 30, main = "GSPC")
acf(period2$N100, lag.max = 30, main = "N100")
acf(period2$N225, lag.max = 30, main = "M225")
acf(period2$SSEC, lag.max = 30, main = "SSEC")
par(mfrow=c(1,1))
# PACF
par(mfrow=c(2,2))
pacf(period2$GSPC, lag.max = 30, main = "GSPC")
pacf(period2$N100, lag.max = 30, main = "N100")
pacf(period2$N225, lag.max = 30, main = "N225")
pacf(period2$SSEC, lag.max = 30, main = "SSEC")
par(mfrow=c(1,1))
# Jarque Bera tests for normality 
JB2_GSPC <- JarqueBera.test(period2$GSPC)
JB2_GSPC
round(JB2_GSPC[[1]][["statistic"]][["X-squared"]], 4)
JB2_N100 <- JarqueBera.test(period2$N100)
JB2_N100
round(JB2_N100[[1]][["statistic"]][["X-squared"]], 4)
JB2_N225 <- JarqueBera.test(period2$N225)
JB2_N225
round(JB2_N225[[1]][["statistic"]][["X-squared"]], 4)
JB2_SSEC <- JarqueBera.test(period2$SSEC)
JB2_SSEC
round(JB2_SSEC[[1]][["statistic"]][["X-squared"]], 4)
# None of the series are normally distributed
# STATIONARITY: 
# ADF unit root test for stationarity of time-series 
adf2GSPC <- ur.df(period2$GSPC, type = "none", lags = trunc((length(period2$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf2GSPC)
adf2N100 <- ur.df(period2$N100, type = "none", lags = trunc((length(period2$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf2N100)
adf2N225 <- ur.df(period2$N225, type = "none", lags = trunc((length(period2$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf2N225)
adf2SSEC <- ur.df(period2$SSEC, type = "none", lags = trunc((length(period2$GSPC)-1)^(1/3)), selectlags = "AIC")
summary(adf2SSEC)
# KPSS test as complement to ADF test
kpss2GSPC <- kpss.test(period2$GSPC, lag.short = TRUE, output = TRUE)
kpss2GSPC
kpss2N100 <- kpss.test(period2$N100, lag.short = TRUE, output = TRUE)
kpss2N100
kpss2N225 <- kpss.test(period2$N225, lag.short = TRUE, output = TRUE)
kpss2N225
kpss2SSEC <- kpss.test(period2$SSEC, lag.short = TRUE, output = TRUE)
kpss2SSEC
# Ljung-Box test for autocorrelation.
LB2_GSPC <- Box.test(period2$GSPC, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
LB2_GSPC
round(LB2_GSPC[["statistic"]][["X-squared"]], 4)
LB2_N100 <- Box.test(period2$N100, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
LB2_N100
round(LB2_N100[["statistic"]][["X-squared"]], 4)
LB2_N225 <- Box.test(period2$N225, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
LB2_N225
round(LB2_N225[["statistic"]][["X-squared"]], 4)
LB2_SSEC <- Box.test(period2$SSEC, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
LB2_SSEC
round(LB2_SSEC[["statistic"]][["X-squared"]], 4)
# Find best fitting ARMA model and fitting it with arima()
GSPC2_aic <- auto.arima(period2$GSPC,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
GSPC2_aic
arma_GSPC2 <- arima(period2$GSPC, order = c(1,0,0))
arma_GSPC2
N1002_aic <- auto.arima(period2$N100,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
N1002_aic
arma_N1002 <- arima(period2$N100, order = c(1,0,0))
arma_N1002
N2252_aic <- auto.arima(period2$N225,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
N2252_aic
arma_N2252 <- arima(period2$N225, order = c(1,0,0))
arma_N2252
SSEC2_aic <- auto.arima(period2$SSEC,max.p = 7, max.q = 7, start.p = 0, start.q = 0, stationary = F, seasonal = FALSE, ic = "aic", trace = TRUE)
SSEC2_aic
arma_SSEC2 <- arima(period2$SSEC, order = c(1,0,0))
arma_SSEC2
# LM-ARCH test for conditional heteroskedacticity 
ARCH2_GSPC <- ArchTest(period2$GSPC,lags = 5)
ARCH2_GSPC
ARCH2_N100 <- ArchTest(period2$N100,lags = 5)
ARCH2_N100
ARCH2_N225 <- ArchTest(period2$N225,lags = 5)
ARCH2_N225
ARCH2_SSEC <- ArchTest(period2$SSEC,lags = 5)
ARCH2_SSEC
#arch.test(arma_GSPC2, output = T)
#arch.test(arma_N1002, output = T)
#arch.test(arma_N2252, output = T)
#arch.test(arma_SSEC2, output = T)

# ARCH-GARCH estimations: Period 2
# GSPC ####
ar21_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(3,1)),variance.model = list(garchOrder = c(1,1), model = "eGARCH"), distribution.model = "sstd")
argarch2_GSPC <- ugarchfit(spec = ar21_garch11_spec, data = period2$GSPC)
argarch2_GSPC
# Specs #
df <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df %>% add_row(
  type_garch = argarch1_GSPC@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_GSPC@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_GSPC)[1],
  BIC = infocriteria(argarch1_GSPC)[2],
  LogLikelihood = argarch1_GSPC@fit[["LLH"]]
) -> df
df %>%
  mutate(model_name = paste0("AR(1,0)-",type_garch,"(1,1)")) -> mydf
df_long <- tidyr::pivot_longer(data = mydf %>%
                                 select(model_name,
                                        type_garch,
                                        type_dist,
                                        AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long$model_name)
# Visualization
p1 <- ggplot(df_long %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (GSPC)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p1
# N100 ####
ar22_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),variance.model = list(garchOrder = c(1,1), model = "eGARCH"), distribution.model = "sstd")
argarch2_N100 <- ugarchfit(spec = ar22_garch11_spec, data = period2$N100)
argarch2_N100
# Specs #
df1 <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df1 %>% add_row(
  type_garch = argarch1_N100@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_N100@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_N100)[1],
  BIC = infocriteria(argarch1_N100)[2],
  LogLikelihood = argarch1_N100@fit[["LLH"]]
) -> df1
df1 %>%
  mutate(model_name = paste0("AR(4,0)-",type_garch,"(1,1)")) -> mydf1
df_long1 <- tidyr::pivot_longer(data = mydf1 %>%
                                  select(model_name,
                                         type_garch,
                                         type_dist,
                                         AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long1$model_name)
# Visualization
p2 <- ggplot(df_long1 %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (N100)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p2
# N225 ####
ar23_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),variance.model = list(garchOrder = c(1,1), model = "sGARCH"), distribution.model = "std")
argarch2_N225 <- ugarchfit(spec = ar23_garch11_spec, data = period2$N225)
argarch2_N225
# Specs #
df2 <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df2 %>% add_row(
  type_garch = argarch1_N225@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_N225@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_N225)[1],
  BIC = infocriteria(argarch1_N225)[2],
  LogLikelihood = argarch1_N225@fit[["LLH"]]
) -> df2
df2 %>%
  mutate(model_name = paste0("AR(0,0)-",type_garch,"(1,1)")) -> mydf2
df_long2 <- tidyr::pivot_longer(data = mydf2 %>%
                                  select(model_name,
                                         type_garch,
                                         type_dist,
                                         AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long2$model_name)
# Visualization
p3 <- ggplot(df_long2 %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (N225)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p3
# SSEC ####
ar24_garch11_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),variance.model = list(garchOrder = c(1,1), model = "eGARCH"), distribution.model = "sstd")
argarch2_SSEC <- ugarchfit(spec = ar24_garch11_spec, data = period2$SSEC)
argarch2_SSEC
# Specs #
df3_5 <- data.frame(
  type_garch=character(), 
  type_dist=character(), 
  stringsAsFactors=TRUE,
  AIC = numeric(),
  BIC = numeric(),
  LogLikelihood = numeric()) 
# Data gathering, manipluation & visualization 
df3_5 %>% add_row(
  type_garch = argarch1_SSEC@model[["modeldesc"]][["vmodel"]],
  type_dist = argarch1_SSEC@model[["modeldesc"]][["distribution"]],
  AIC = infocriteria(argarch1_SSEC)[1],
  BIC = infocriteria(argarch1_SSEC)[2],
  LogLikelihood = argarch1_SSEC@fit[["LLH"]]
) -> df3_5
df3_1 %>%
  mutate(model_name = paste0("AR(0,0)-",type_garch,"(1,1)")) -> mydf3_1
df_long3_1 <- tidyr::pivot_longer(data = mydf3_1 %>%
                                    select(model_name,
                                           type_garch,
                                           type_dist,
                                           AIC, BIC, LogLikelihood),  cols = c('AIC', 'BIC','LogLikelihood'))
model_name <- unique(df_long3_5$model_name)
merged <- do.call("rbind", list(df_long3,df_long3_1,df_long3_2, df_long3_3, df_long3_4, df_long3_5))
merged %>%
  arrange(model_name) -> merged
# Visualization
p4 <- ggplot(merged %>%
               arrange(type_garch), 
             aes(x = reorder(model_name, 
                             order(type_garch)),
                 y = value, 
                 shape = type_dist,
                 color = type_garch)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  scale_shape_manual(values = c(15,16,17,18))+
  theme_bw(base_family = "Bahnschrift") + 
  facet_wrap(~name, scales = 'free_x') +
  labs(title = 'GARCH model selection (SSEC)', 
       subtitle = 'The best model is the one with lowest AIC, BIC and lowest LogLikelihood in absolute terms',
       x = '',
       y = 'Value of Criteria',
       shape = 'Type of Distribution',
       color = 'Type of GARCH model') + 
  theme(legend.position = "right")
p4

### DCC-GARCH estimation: Period 1 ####
# Combining univariate GARCH spec with multispec
mspec2 <- multispec(c(ar21_garch11_spec,ar22_garch11_spec,ar23_garch11_spec,ar24_garch11_spec))
# DCC specification - GARCH(1,1) for conditional correlations
dccgarch11_spec2 <- dccspec(uspec = mspec2,lag =1, lag.criterion = "AIC", lag.max= 20, dccOrder = c(1,1), distribution = "mvnorm")
dccgarch11_spec2
# Fitting DCC model
dccfit2 <- dccfit(dccgarch11_spec2, data = period2, fit.control = list(eval.se=TRUE))
dccfit2
# Extracting time varying variance-covariance and correlation matrices 
cov2.garch11 <- rcov(dccfit2)
cov2.garch11
cor2.garch11 <- rcor(dccfit2)
cor2.garch11
# Visualization
# Time-varying correlations
par(mfrow=c(3,2))
plot(as.xts(cor2.garch11[1,2,]),main="GSPC and N100", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor2.garch11[1,3,]),main="GSPC and N225", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor2.garch11[1,4,]),main="GSPC and SSEC", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor2.garch11[2,3,]),main="N100 and N225", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor2.garch11[2,4,]),main="N100 and SSEC", yaxis.right = FALSE, ylim = c(-0.1,1))
plot(as.xts(cor2.garch11[3,4,]),main="N225 and SSEC", yaxis.right = FALSE, ylim = c(-0.1,1))
par(mfrow=c(1,1))
# Time-varying covariances
par(mfrow=c(3,2))
plot(as.xts(cov2.garch11[1,2,]),main="GSPC and N100", yaxis.right = FALSE, ylim = c(0,50))
plot(as.xts(cov2.garch11[1,3,]),main="GSPC and N225", yaxis.right = FALSE, ylim = c(0,50))
plot(as.xts(cov2.garch11[1,4,]),main="GSPC and SSEC", yaxis.right = FALSE, ylim = c(0,50))
plot(as.xts(cov2.garch11[2,3,]),main="N100 and N225", yaxis.right = FALSE, ylim = c(0,50))
plot(as.xts(cov2.garch11[2,4,]),main="N100 and SSEC", yaxis.right = FALSE, ylim = c(0,50))
plot(as.xts(cov2.garch11[3,4,]),main="N225 and SSEC", yaxis.right = FALSE, ylim = c(0,50))
par(mfrow=c(1,1))
# DCC model diagnostics ####
# Li-Mak test for presence of autocorrelation in the residuals 
dccLMtest_GSPC2 <- Weighted.LM.test(dccfit2@model$residuals[,1], dccfit2@model$sigma[,1]^2, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 6, weighted = FALSE)
dccLMtest_GSPC2
dccLMtest_N1002 <- Weighted.LM.test(dccfit2@model$residuals[,2], dccfit2@model$sigma[,2]^2, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
dccLMtest_N1002
dccLMtest_N2252 <- Weighted.LM.test(dccfit2@model$residuals[,3], dccfit2@model$sigma[,3]^2, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
dccLMtest_N2252
dccLMtest_SSEC2 <- Weighted.LM.test(dccfit2@model$residuals[,4], dccfit2@model$sigma[,4]^2, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
dccLMtest_SSEC2
# LM-ARCH test on standardized residuals 
DCCARCH2_GSPC <- ArchTest(dccfit2@mfit$stdresid[,1],lags = 5)
DCCARCH2_GSPC
DCCARCH2_N100 <- ArchTest(dccfit2@mfit$stdresid[,2],lags = 5)
DCCARCH2_N100
DCCARCH2_N225 <- ArchTest(dccfit2@mfit$stdresid[,3],lags = 5)
DCCARCH2_N225
DCCARCH2_SSEC <- ArchTest(dccfit2@mfit$stdresid[,4],lags = 5)
DCCARCH2_SSEC
# Ljung-Box test for autocorrelation in squared standardized residuals. 
dccLB2_GSPC <- Box.test(dccfit2@mfit$stdresid[,1]^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
dccLB2_GSPC
round(dccLB2_GSPC[["statistic"]][["X-squared"]], 4)
dccLB2_N100 <- Box.test(dccfit2@mfit$stdresid[,2]^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
dccLB2_N100
round(dccLB2_N100[["statistic"]][["X-squared"]], 4)
dccLB2_N225 <- Box.test(dccfit2@mfit$stdresid[,3]^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
dccLB2_N225
round(dccLB2_N225[["statistic"]][["X-squared"]], 4)
dccLB2_SSEC <- Box.test(dccfit2@mfit$stdresid[,4]^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
dccLB2_SSEC
round(dccLB2_SSEC[["statistic"]][["X-squared"]], 4)
# Visualization
par(mfrow= c(2,2))
plot(dccfit2@mfit$stdresid[,1],type= "l", main = "GSPC", ylab = "Standardised residuals", ylim =c(-4,4))
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit2@mfit$stdresid[,1]), b=0, col = "blue")
plot(dccfit2@mfit$stdresid[,2],type= "l", main = "N100", ylab = "Standardised residuals",ylim =c(-4,4))
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit2@mfit$stdresid[,2]), b=0, col = "blue")
plot(dccfit2@mfit$stdresid[,3],type= "l", main = "N225", ylab = "Standardised residuals",ylim =c(-4,4))
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit2@mfit$stdresid[,3]), b=0, col = "blue")
plot(dccfit2@mfit$stdresid[,4],type= "l", main = "SSEC", ylab = "Standardised residuals",ylim =c(-8,4))
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dccfit2@mfit$stdresid[,2]), b=0, col = "blue")
par(mfrow= c(1,1))

### BEKK-GARCH estimation: Period 1 ####
# Convert period 1 to data frame and export to excel for use in Oxmetrics
#tmp <- "C:\\Users\\Sundell\\Documents\\EMF 2020\\Term Paper\\dataframe3"
#write.zoo(period1,sep=",",file=tmp)
period2_df <- as.data.frame(period2)
write_xlsx(period2_df,"C:\\Users\\Sundell\\Documents\\EMF 2020\\dataframe2.xlsx")
# Model Diagnostics
BEKK2LMtest_GSPC <- Weighted.LM.test(dataframe2$Res_GSPC, dataframe2$CondVar_GSPC, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKK2LMtest_GSPC
round(BEKK2LMtest_GSPC[["statistic"]][["X-squared on Squared Residuals for fitted ARCH process"]], 4)
BEKK2LMtest_N100 <- Weighted.LM.test(dataframe2$Res_N100, dataframe2$CondVar_N100, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKK2LMtest_N100
round(BEKK2LMtest_N100[["statistic"]][["X-squared on Squared Residuals for fitted ARCH process"]], 4)
BEKK2LMtest_N225 <- Weighted.LM.test(dataframe2$Res_N225, dataframe2$CondVar_N225, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKK2LMtest_N225
round(BEKK2LMtest_N225[["statistic"]][["X-squared on Squared Residuals for fitted ARCH process"]], 4)
BEKK2LMtest_SSEC <- Weighted.LM.test(dataframe2$Res_SSEC, dataframe2$CondVar_SSEC, lag = trunc((length(period2)-1)^(1/3)), type = c("partial"), fitdf = 2, weighted = FALSE)
BEKK2LMtest_SSEC
round(BEKK2LMtest_SSEC[["statistic"]][["X-squared on Squared Residuals for fitted ARCH process"]], 4)
# LM-ARCH test on standardized residuals 
BEKKARCH2_GSPC <- ArchTest(dataframe2$StdRes_GSPC,lags = 5)
BEKKARCH2_GSPC
round(BEKKARCH2_GSPC[["statistic"]][["Chi-squared"]], 4)
BEKKARCH2_N100 <- ArchTest(dataframe2$StdRes_N100,lags = 5)
BEKKARCH2_N100
round(BEKKARCH2_N100[["statistic"]][["Chi-squared"]], 4)
BEKKARCH2_N225 <- ArchTest(dataframe2$StdRes_N225,lags = 5)
BEKKARCH2_N225
round(BEKKARCH2_N225[["statistic"]][["Chi-squared"]], 4)
BEKKARCH2_SSEC <- ArchTest(dataframe2$StdRes_SSEC,lags = 5)
BEKKARCH2_SSEC
round(BEKKARCH2_SSEC[["statistic"]][["Chi-squared"]], 4)
# Ljung-Box test for autocorrelation in squared standardized residuals. 
BEKKLB2_GSPC <- Box.test(dataframe2$StdRes_GSPC^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
BEKKLB2_GSPC
round(BEKKLB2_GSPC[["statistic"]][["X-squared"]], 4)
JB1_GSPC <- JarqueBera.test(period1$GSPC)
JB1_GSPC
round(JB1_GSPC[[1]][["statistic"]][["X-squared"]], 4)
BEKKLB2_N100 <- Box.test(dataframe2$StdRes_N100^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
BEKKLB2_N100
round(BEKKLB2_N100[["statistic"]][["X-squared"]], 4)
BEKKLB2_N225 <- Box.test(dataframe2$StdRes_N225^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
BEKKLB2_N225
round(BEKKLB2_N225[["statistic"]][["X-squared"]], 4)
BEKKLB2_SSEC <- Box.test(dataframe2$StdRes_SSEC^2, lag=sqrt(length(period2$GSPC)), type = "Ljung-Box")
BEKKLB2_SSEC
round(BEKKLB2_SSEC[["statistic"]][["X-squared"]], 4)
#Residuals visualized
par(mfrow= c(2,2))
plot(dataframe2$StdRes_GSPC,type= "l", main = "GSPC", ylab = "Standardised residuals", ylim = c(-4,4))
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe2$StdRes_GSPC), b=0, col = "blue")
plot(dataframe2$StdRes_N100,type= "l", main = "N100", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe2$StdRes_N100), b=0, col = "blue")
plot(dataframe2$StdRes_N225,type= "l", main = "N225", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe2$StdRes_N225), b=0, col = "blue")
plot(dataframe2$StdRes_SSEC,type= "l", main = "SSEC", ylab = "Standardised residuals")
abline(a=-3, b=0, col = "red")
abline(a=3, b=0, col = "red")
abline(a = mean(dataframe2$StdRes_SSEC), b=0, col = "blue")
par(mfrow= c(1,1))

# Mean return shocks period 1 & 2 ####
GSPC_meanshock1 <- mean(dataframe1$Res_GSPC)
N100_meanshock1 <- mean(dataframe1$Res_N100)
N225_meanshock1 <- mean(dataframe1$Res_N225)
SSEC_meanshock1 <- mean(dataframe1$Res_SSEC)
GSPC_meanshock2 <- mean(dataframe2$Res_GSPC)
N100_meanshock2 <- mean(dataframe2$Res_N100)
N225_meanshock2 <- mean(dataframe2$Res_N225)
SSEC_meanshock2 <- mean(dataframe2$Res_SSEC)
 
GSPC_meanshock1 
N100_meanshock1
N225_meanshock1 
SSEC_meanshock1 
GSPC_meanshock2 
N100_meanshock2 
N225_meanshock2 
SSEC_meanshock2 