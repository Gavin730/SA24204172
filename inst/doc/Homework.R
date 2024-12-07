## -----------------------------------------------------------------------------
n <- 10000
u <- runif(n)#generate random number
Rayleigh <- function(sigma){
  x <- sigma * sqrt(-2 * log(1-u)) #Inverse transform algorithm
  hist(x,breaks=100,prob = TRUE)#plot density 
  
  #title
  usr <- par("usr")
  text(x = usr[2] - 0.3 * (usr[2] - usr[1]), 
     y = usr[4] - 0.3 * (usr[4] - usr[3]), 
     expression(f(x) == frac(x, sigma^2) * e^{-frac(x^2, 2 * sigma^2)} * "," ~ x >= 0 ~ ","))
  text(x = usr[2] - 0.13 * (usr[2] - usr[1]), 
     y = usr[4] - 0.3 * (usr[4] - usr[3]),  substitute(sigma == sigma_val, list(sigma_val = sigma)))
  
  y <- seq(0, 100, .01) 
  lines(y, (y / sigma^2) * exp(-y^2 / (2 * sigma^2))) # plot oracle density function
}
Rayleigh(1)
Rayleigh(2)
Rayleigh(3)

## -----------------------------------------------------------------------------
p1 <- 0.25
n <- 100000

mixgauss <-function(p,n){
  x <- rep(0, n)
  for(i in 1:n){
    a <- runif(1) #generate random number
    if (a < p){  #check if a < p1
      x[i] <- rnorm(1)
    }
    else{
      x[i] <- 3+ rnorm(1)
    }
  }
  return(x)
}
#par(mfrow = c(2, 2))
hist(mixgauss(0.65,n),breaks=150,prob = TRUE,main = expression("Mixed Gaussian Distribution, " ~ p[1] == 0.65))
hist(mixgauss(0.7,n),breaks=150,prob = TRUE,main = expression("Mixed Gaussian Distribution, " ~ p[1] == 0.7))
hist(mixgauss(0.75,n),breaks=150,prob = TRUE,main = expression("Mixed Gaussian Distribution, " ~ p[1] == 0.75))
hist(mixgauss(0.8,n),breaks=150,prob = TRUE,main = expression("Mixed Gaussian Distribution, " ~ p[1] == 0.8))

## -----------------------------------------------------------------------------
lambda = 0.5
alpha <-2
beta <- 2
X<- function(t){#generate compound Poisson process
  n = rpois(1,lambda*t)
  result = 0
  for (i in 1:n){
    result = result + rgamma(1,shape = alpha,rate = beta)
  }
  return(result)
}

n <- 10000
data <- rep(0, n)
for (i in 1:n){
  data[i] = X(10)
}
{
cat("When lambda = ",lambda,"alpha = ",alpha,"beta = ",beta,"\n")
cat("The mean of X(10) is", mean(data),"\n")
cat("The variance of X(10) is", var(data),"\n")
cat("λtE[Y_1] = ", lambda*10*alpha/beta,"\n")
cat("λtE[Y_1^2] = ", lambda*10*alpha*(beta+1)/beta^2,"\n\n")}

#The second set of parameters 
lambda = 1
alpha <-1
beta <- 2
data <- rep(0, n)
for (i in 1:n){
  data[i] = X(10)
}
{
cat("When lambda = ",lambda,"alpha = ",alpha,"beta = ",beta,"\n")
cat("The mean of X(10) is", mean(data),"\n")
cat("The variance of X(10) is", var(data),"\n")
cat("λtE[Y_1] = ", lambda*10*alpha/beta,"\n")
cat("λtE[Y_1^2] = ", lambda*10*alpha*(alpha+1)/beta^2,"\n\n")}

#The third set of parameters
lambda = 1
alpha <-2
beta <- 0.5
data <- rep(0, n)
for (i in 1:n){
  data[i] = X(10)
}
{
cat("When lambda = ",lambda,"alpha = ",alpha,"beta = ",beta,"\n")
cat("The mean of X(10) is", mean(data),"\n")
cat("The variance of X(10) is", var(data),"\n")
cat("λtE[Y_1] = ", lambda*10*alpha/beta,"\n")
cat("λtE[Y_1^2] = ", lambda*10*alpha*(alpha+1)/beta^2,"\n\n")}

## -----------------------------------------------------------------------------
beta_cdf <- function(x,a,b){#模特卡罗模拟beta分布cdf
  m <- 1e4; u <- runif(m, min=0, max=x)
  estimation <- mean(x*u^2*(1-u)^2/beta(a,b)) 
  return (estimation)
}
for (x in seq(0.1, 0.9, 0.1)){
  cat("x=",x,"\n")
  cat("蒙特卡洛模拟累计分布函数值：",beta_cdf(x,3,3),"\n")
  cat("真实累计分布函数值：",pbeta(x,3,3),"\n")
}

## -----------------------------------------------------------------------------
n <- 10000
u <- runif(n)#generate random number

Rayleigh_1 <- function(sigma){
  x <- sigma * sqrt(-2 * log(1-u)) #Inverse transform algorithm
}
Rayleigh_antithetic <- function(sigma){ #generated by antithetic variables
  x <- sigma * sqrt(-2 * log(u)) 
}
Rayleigh_2 <- function(sigma){
  u <- runif(n)
  x <- sigma * sqrt(-2 * log(1-u))
}

data <- matrix(1:9, nrow = 3, ncol = 3)
rownames(data) <- c('var((X + X\')/ 2)' , 'var((X_1 + X_2)/ 2)', 'reduction')
colnames(data) <- c("σ=1", "σ=2", "σ=3")
for (i in 1:nrow(data)){
  for (j in 1:ncol(data)){
    if (i == 1){
    data[i,j] <- var(Rayleigh_1(j)+Rayleigh_antithetic(j))/4}
    if (i == 2){
    data[i,j] <- var(Rayleigh_1(j)+Rayleigh_2(j))/4}
    if (i == 3){
    data[i,j] <- 1-data[1,j]/data[2,j]}
  }
}
data

## -----------------------------------------------------------------------------
f_1 <- function(x){#指数分布
  return (exp(-(x-1)))
}
f_2 <- function (x){#半正态分布
  return(2/sqrt(2*pi)*exp(-(x-1)^2/2))
}
g <- function (x){
  return (x^2/sqrt(2*pi)*exp(-x^2/2))
}


n <- 1000000
x_1 <- rexp(n, 1) +1
x_2 <- abs(rnorm(n,0)) +1
{
  cat("通过f_1的估计",mean(g(x_1)/f_1(x_1)),'\n')
  cat("通过f_2的估计",mean(g(x_2)/f_2(x_2)),'\n')
  cat('f_1/‖f_1‖-g/‖g‖的L1范数',integrate(function(x) abs(f_1(x)-g(x)/integrate(function(x) g(x), lower = 1, upper = Inf)$value), lower = 1, upper = Inf)$value,'\n')
  cat('f_2/‖f_2‖-g/‖g‖的L1范数',integrate(function(x) abs(f_2(x)-g(x)/integrate(function(x) g(x), lower = 1, upper = Inf)$value), lower = 1, upper = Inf)$value,'\n')
  cat("通过f_1的方差",var (g(x_1)/f_1(x_1)),'\n')
  cat("通过f_2的方差",var (g(x_2)/f_1(x_2)),'\n')
}

## -----------------------------------------------------------------------------

fast_sorting_algorithm <- function(x){#快排算法
  if (length(x) <= 1){
    return (x)}
  a=sample(1:length(x),1)
  pivot <- x [a]
  left <- fast_sorting_algorithm(x[x<pivot])
  middle <- x[x==pivot]
  right <- fast_sorting_algorithm(x[x>pivot])
  return (c(left,middle,right))
}
library(ggplot2)
n_values <- c(1e4, 2 * 1e4, 4 * 1e4, 6 * 1e4, 8 * 1e4)
a_n <- rep(0,5)
count <- 1

for(n in n_values) {
  random_nums <- sample(1:n, n)
  
  cat("n =", n, "\n")
  temp = 0
  for (i in 1:10){
    start_time <- Sys.time()  # 记录开始时间
    sorted_nums <- fast_sorting_algorithm(random_nums)  # 快速排序
    end_time <- Sys.time()  # 记录结束时间
    temp <- temp + end_time - start_time
  }
  a_n[count] = temp / 10
  cat("排序用时:",a_n[count] , "秒\n")
  cat('n*logn=',n*log(n),"\n")
  cat("-----------------------\n")
  count = count+1
}
t_n <- n_values * log(n_values)
model <- lm(a_n ~ t_n)
summary(model)

ggplot(data = data.frame(t_n, a_n), aes(x = t_n, y = a_n)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "n log(n)",
       y = "a_n")

## -----------------------------------------------------------------------------
n <- 1000  # 样本大小
t <- 10000  # 模拟次数

#蒙特卡洛模拟
skewness <- function(n,t) {
  skewness<- numeric(t)
  for (i in 1:t){
    x <- rnorm(n)
    mean_x <- mean(x)
    s <- sd(x)
    skewness[i] <- sum((x - mean_x)^3) / (n * s^3)
  }
  return(skewness)
}

b1_hat <- skewness(n,t)

# 计算偏度分位数
quant <- c(0.025, 0.05, 0.95, 0.975)
b1_hat_quant <- quantile(b1_hat, probs = quant)
#计算标准差
b1_density <- dnorm(b1_hat_quant, mean = 0, sd = sqrt(6/n))
b1_se <- sqrt(quant * (1 - quant) / (n * b1_density^2))
#大样本近似分位数
b1_quant <- qnorm(quant, mean = 0, sd = sqrt(6/n))

result <- matrix(1:12, nrow = 3, ncol = 4)
rownames(result) <- c('蒙特卡洛偏度分位数估计' , '标准误差', '大样本近似估计')
colnames(result) <- c("2.5%", "5%", "95%","97.5%")
for (i in 1:nrow(result)){
    if (i == 1){
    result[i,] <- b1_hat_quant}
    if (i == 2){
    result[i,] <- b1_se}
    if (i == 3){
    result[i,] <- b1_quant}
}
print(result)

## -----------------------------------------------------------------------------
library(MASS) 
n <- 100
t <- 1000
alpha <- 0.05
Sigma <- matrix(c(1, 0.3, 0.3, 1), ncol=2)#协方差阵
p_values <- matrix(NA, ncol=3, nrow=t)
#对多元正态分布
for (i in 1:t) {
    data <- mvrnorm(n=n, mu=c(0,0), Sigma=Sigma)
    # Pearson相关性检验
    p_values[i, 1] <- cor.test(data[,1], data[,2], method="pearson")$p.value
    # Spearman相关性检验
    p_values[i, 2] <- cor.test(data[,1], data[,2], method="spearman")$p.value
    # Kendall相关性检验
    p_values[i, 3] <- cor.test(data[,1], data[,2], method="kendall")$p.value
}
#计算检验功效
power <- colMeans(p_values < alpha)
names(power) <- c("Pearson", "Spearman", "Kendall")
print(power)

## -----------------------------------------------------------------------------
#对非正态分布
#x为指数分布，y为x/3加上指数分布随机扰动
for (i in 1:t) {
  x <- rexp(n,1)
  y <- x/3 + rexp(n,1)
  data <- cbind(x,y)
  # Pearson相关性检验
  p_values[i, 1] <- cor.test(data[,1], data[,2], method="pearson")$p.value
  # Spearman相关性检验
  p_values[i, 2] <- cor.test(data[,1], data[,2], method="spearman")$p.value
  # Kendall相关性检验
  p_values[i, 3] <- cor.test(data[,1], data[,2], method="kendall")$p.value
}
power <- colMeans(p_values < alpha)
names(power) <- c("Pearson", "Spearman", "Kendall")
print(power)

## -----------------------------------------------------------------------------
N <- 1000  
n_null <- 950  
n_alt <- 50  
alpha <- 0.1  
m <- 10000  
set.seed(123)
results <- matrix(0, nrow = 3, ncol = 2)
colnames(results) <- c("Bonferroni correction", "B-H correction")
rownames(results) <- c("FWER", "FDR", "TPR")

#校正函数，返回FWER，FDR,TPR
correction <- function(p,method,alpha,n_null,n_alt){
  N=n_null+n_alt
  p_correction <- p.adjust(p,method=method)
  reject <- p_correction<alpha
  false_pos <- sum(reject[1:n_null])
  true_pos <- sum(reject[(n_null + 1):N])
  if (false_pos>0){FWER=1}
  else FWER=0
  FDR <- false_pos/sum(reject)
  TPR <- true_pos/n_alt
  return (c(FWER,FDR,TPR))
}

for (i in 1:m) {
  p_null <- runif(n_null)  
  p_alt <- rbeta(n_alt, 0.1, 1)  
  p_values <- c(p_null, p_alt)
  #bonferroni校正
  p_bf <- correction(p_values,"bonferroni",alpha,n_null,n_alt)
  results["FWER", "Bonferroni correction"] <- results["FWER", "Bonferroni correction"]+p_bf[1]
  results["FDR", "Bonferroni correction"] <- results["FDR", "Bonferroni correction"]+p_bf[2]
  results["TPR", "Bonferroni correction"] <- results["TPR", "Bonferroni correction"]+p_bf[3]
  
  # B-H校正
  p_bh <- correction(p_values,"BH",alpha,n_null,n_alt)
  results["FWER", "B-H correction"] <- results["FWER", "B-H correction"]+p_bh[1]
  results["FDR", "B-H correction"] <- results["FDR", "B-H correction"]+p_bh[2]
  results["TPR", "B-H correction"] <- results["TPR", "B-H correction"]+p_bh[3]
}

# 计算平均值
results["FWER", ] <- results["FWER", ] / m
results["FDR", ] <- results["FDR", ] / m
results["TPR", ] <- results["TPR", ] / m
print(results)

## -----------------------------------------------------------------------------

library(boot)
data <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
#mle估计
mle_lambda <- 1 / mean(data)

lambda <- function(data, i) {
  return(1 / mean(data[i]))  
}
# bootstrap估计
boot_results <- boot(data, lambda, R = 1000)
boot_bias <- mean(boot_results$t) - mle_lambda
boot_se <- sd(boot_results$t)  
{
cat("最大似然估计 : ", mle_lambda, "\n")
cat("估计偏差: ", boot_bias, "\n")
cat("估计标准误差: ", boot_se, "\n")}



## -----------------------------------------------------------------------------
datamean <- function(x, i) mean(x[i])

de <- boot(data=data, statistic=datamean, R=1000)

ci <- boot.ci(de, type=c("norm", "basic", "perc", "bca"))
{
cat("standard normal CI: ", ci$norm[2:3], "\n")
cat("basic CI: ", ci$basic[4:5], "\n")
cat("percentile CI: ", ci$percent[4:5], "\n")
cat("BCa CI: ", ci$bca[4:5], "\n")}

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())

library(boot)
library(bootstrap)
data(scor)

cov_matrix <- cov(scor)
eigenvalues <- eigen(cov_matrix)$values
theta_hat <- eigenvalues[1] / sum(eigenvalues)
#jacknife估计
jackknife <- function(data) {
  n <- nrow(data)
  theta_jack <- numeric(n)
  for (i in 1:n) {
    eigen_jack <- eigen(cov(data[-i, ]))$values
    theta_jack[i] <- eigen_jack[1] / sum(eigen_jack)
  }
  
  bias_jack <- (n - 1) * (mean(theta_jack) - theta_hat)
  se_jack <- sqrt((n - 1) * mean((theta_jack - theta_hat)^2))
  
   round(c(original=theta_hat,bias_jack=bias_jack,
 se_jack=se_jack),3)
}

jackknife_results <- jackknife(scor)
jackknife_results

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
library(DAAG)
data(ironslag)


magnetic <- ironslag$magnetic
chemical <- ironslag$chemical
n <- length(magnetic)

e1 <- e2 <- e3 <- e4 <- numeric(n)
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  # 线性模型
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  # 二次多项式模型
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  # 对数-线性模型
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  # 三次多项式模型
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + 
           J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}

# 预测误差
pred_error <- array(c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)),dim = c(1, 4))
colnames(pred_error) <- c("Linear", "Quadratic", "Log-Linear", "Cubic")
rownames(pred_error)<-c('prediction error')
pred_error
# 调整后的R^2
r2 <- array(c(summary(lm(magnetic ~ chemical))$adj.r.squared,
            summary(lm(magnetic ~ chemical + I(chemical^2)))$adj.r.squared,
            summary(lm(log(magnetic) ~ chemical))$adj.r.squared,
            summary(lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3)))$adj.r.squared),dim=c(1,4))
colnames(r2) <- c("Linear", "Quadratic", "Log-Linear", "Cubic")
rownames(r2) <- c("adjusted R2")
r2


## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
set.seed(123)
library(DAAG)
data(chickwts)
set.seed(123)
x <- sort(as.vector(chickwts$weight[chickwts$feed == "soybean"]))
y <- sort(as.vector(chickwts$weight[chickwts$feed == "linseed"]))

data <- c(x, y)
n_x <- length(x)
n_y <- length(y)
n <- n_x + n_y
R <- 1000  
cv_stat <- numeric(R)  

#Cramer-von Mises统计量
cvm <- function(data, indices) {
  x <- data[indices[1:n_x]]
  y <- data[indices[(n_x + 1):(n_x + n_y)]]
  ecdf_x <- ecdf(x)
  ecdf_y <- ecdf(y)
  W <- (1 / (n_x * n_y)) * sum((ecdf_x(data) - ecdf_y(data))^2)
  return(W)
}

boot_result <- boot(data,cvm, R)
W <- boot_result$t0
p_value <- mean(c(boot_result$t0,boot_result$t) >= W)


{cat("Cramer-von Mises统计量 W=", W, "\n")
cat("p值:", p_value, "\n")}
hist(boot_result$t, main = "Cramer-von Mises统计量的分布", 
     xlab = "W", freq = FALSE, breaks = "scott")
points(W, 0, cex = 1, pch = 16, col = "red")

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
library(boot)
# 生成数据
set.seed(123)
n <- 30
x <- rnorm(n)
y <- 0.1*x + rnorm(n)

# Spearman秩相关系数
spearman <- function(data, indices) {
  d <- data[indices, ]
  return(cor(d[, 1], d[, 2], method = "spearman"))
}

data <- data.frame(x, y)

R<- 10000
boot_results <- boot(data, spearman, R)
spearman_stat <- boot_results$t0
p_perm <- mean(abs(c(boot_results$t0,boot_results$t)) >= abs(boot_results$t0))
cor_test <- cor.test(x, y, method = "spearman")
{cat("置换检验的p值:", p_perm, "\n")
cat("cor.test 报告的p值:", cor_test$p.value, "\n")}

# 绘制分布直方图
hist(boot_results$t, main = "Spearman统计量的置换分布", xlab = "Spearman 统计量", breaks = "scott")
points(spearman_stat, 0, cex = 1.5, pch = 16, col = "red")

## -----------------------------------------------------------------------------
rm(list = ls())

target_pdf <- function(x) (1 / (pi * (1 + x^2)))#目标分布
n <- 1000000
set.seed(123)
# Metropolis-Hastings算法，采用正态作为提议分布
m_h <- function( n) {
  samples <- numeric(n+100000)
  samples[1] <- 0
  
  for (i in 2:n) {
    proposal <- rnorm(1, mean = samples[i - 1], sd = 1)
    accept_ratio <- target_pdf(proposal) / target_pdf(samples[i - 1])
    
    if (runif(1) < accept_ratio)  samples[i] <- proposal
    else  samples[i] <- samples[i - 1]
  }
  return(samples[-(1:100000)])#舍弃前1000个样本
}

samples <- m_h(n)


# 比较十分位数
{cat("生成样本的十分位数:\n")
print(quantile(samples,seq(0, 1, 0.1)))
cat("\n标准柯西分布的十分位数:\n")
print(qcauchy(seq(0, 1, 0.1)))}

#Gelman-Rubin
chains <- list()
k=3
for (i in 1:k) {
  chains[[i]] <- m_h(n)
}

# 计算样本均值和方差
chain_means <- sapply(chains, mean)
chain_vars <- sapply(chains, var)
overall_mean <- mean(chain_means)

# 链间方差
B <- (n / (k - 1)) * sum((chain_means - overall_mean)^2)
# 链内方差
W <- mean(chain_vars)
# 总体方差
V <- ((n - 1) / n) * W + (1 / n) * B
R_hat <- sqrt(V / W)

if (R_hat < 1.2) {
  cat("链已收敛，R_hat=",R_hat,"\n")
} else {
  cat("链未收敛，R_hat=",R_hat,"\n")
}

## -----------------------------------------------------------------------------
rm(list = ls())
n <- 100  
a <- 2    
b <- 3    
iter<- 5000  
#初始化
x <- 1  
y <- runif(1)     

gibbs <- function(iter){
  results <- matrix(0, nrow = iter, ncol = 2)
  colnames(results) <- c("x", "y")
  for (i in 1:iter) {
    # 生成 x
    x <- rbinom(1, n, y)
    #生成 y
    results[i, ] <- c(x, y)
    y <- rbeta(1, x + a, n - x + b)
  }
  return(results[-(1:1000),])
}

#Gelman-Rubin
R_hat<-function(j){
  chains <- list()
  k=3
  for (i in 1:k) {
    chains[[i]] <- gibbs(iter)[ ,j]
  }
  chain_means <- sapply(chains, mean)
  chain_vars <- sapply(chains, var)
  overall_mean <- mean(chain_means)
  B <- (n / (k - 1)) * sum((chain_means - overall_mean)^2)
  W <- mean(chain_vars)
  V <- ((n - 1) / n) * W + (1 / n) * B
  R_hat <- sqrt(V / W)
  return(R_hat)
}


if (R_hat(1) < 1.2) {
  cat("链已收敛，R_hat(x),R_hat(y)=",c(R_hat(1),R_hat(2)),"\n")
} else {
  cat("链未收敛，R_hat(x),R_hat(y)=",c(R_hat(1),R_hat(2)),"\n")
}

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
k_term <- function(k, a, d) {
  #取对数
  term_log <-  - (lgamma(k + 1) + k * log(2)) + 
              (2 * k + 2) * log(sqrt(sum(a^2))) - log(2 * k + 1) - log(2 * k + 2) +
              lgamma((d + 1) / 2) + lgamma(k + 1.5) - lgamma(k + d / 2 + 1)
  term <- (-1)^k*exp(term_log)
  return(term)
}

## -----------------------------------------------------------------------------
function_sum <- function(a, d) {
  sum_result <- 0
  k<- 0:100
  term_log <-  - (lgamma(k + 1) + k * log(2)) + 
              (2 * k + 2) * log(sqrt(sum(a^2))) - log(2 * k + 1) - log(2 * k + 2) +
              lgamma((d + 1) / 2) + lgamma(k + 1.5) - lgamma(k + d / 2 + 1)
  term <- (-1)^k*exp(term_log)
  return(sum(term))
}

## -----------------------------------------------------------------------------

a <- c(1, 2)
d <- length(a)
result <- function_sum(a, d)
print(result)

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
c_k <- function(a,k) {
  sqrt(a^2 * k / (k + 1 - a^2))
}

lhs <- function(a,k){
   2 * exp(lgamma(k / 2)-lgamma((k - 1) / 2)) / sqrt(pi * (k - 1))  * integrate(function(u) (1 + u^2 / (k - 1))^(-k / 2), 0, c_k(a,k-1))$value
}

rhs <- function(a,k){
  2 * exp(lgamma((k+1) / 2)-lgamma(k/ 2)) / sqrt(pi * k) * integrate(function(u) (1 + u^2 / k)^(-(k + 1) / 2), 0, c_k(a,k))$value
}



solve_a <- function(k) {
  # 方程的解
  f <- function(a) {
    lhs(a,k)-rhs(a,k)
  }
  
  # 使用 uniroot 求解
  solution_a <- uniroot(f, lower = 0.1, upper = min(sqrt(k),2))$root
  return(solution_a)
}

k_values <- c(4:25, 100, 500, 1000)
solution <- sapply(k_values, solve_a)
print(rbind(k_values,solution))

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
S_k1 <- function(a, k) {
  1 - pt(sqrt(a^2 * (k - 1) / (k - a^2)), df = k - 1)
}
S_k <- function(a, k) {
  1 - pt(sqrt(a^2 * k / (k + 1 - a^2)), df = k)
}

solve_a <- function(k) {
  # 方程的解
  f <- function(a) {
    S_k1(a,k)-S_k(a,k)
  }
  
  # 使用 uniroot 求解
  solution_a <- uniroot(f, lower = 0.1, upper = min(sqrt(k),2))$root
  return(solution_a)
}

k_values <- c(4:25, 100, 500, 1000)
solution <- sapply(k_values, solve_a)
print(rbind(k_values,solution))

## -----------------------------------------------------------------------------
observed_data <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1  
n <- length(observed_data)

# E-M 算法
em <- function(data, tau, tol = 1e-6, max_iter = 1000) {
  # 初始化 lambda
  lambda <- 1 / mean(data)  # 初始值
  
  for (i in 1:max_iter) {
    # E步
    truncated_data <- data[data == tau]  
    untruncated_data <- data[data < tau] 
    expected_truncated <- sum(truncated_data + 1 / lambda)
    
    # M步
    lambda_new <- n / (sum(untruncated_data) + expected_truncated)
    
    if (abs(lambda - lambda_new) < tol) break
    lambda <- lambda_new
  }
  return(lambda)
}

lambda_em <- em(observed_data, tau)

# 使用观测数据直接计算 MLE  
log_likelihood <- function(lambda) {
  -sum(ifelse(observed_data < tau, log(lambda) - lambda * observed_data, -lambda * tau))
}
mle_result <- optim(par = 1, log_likelihood, method = "Brent", lower = 1e-5, upper = 10)
lambda_mle<- mle_result$par


# 输出结果
cat("EM 算法估计的 λ:", lambda_em, "\n")
cat("MLE 估计的 λ:", lambda_mle, "\n")

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
library(boot) 
# 定义约束
A1 <- rbind(c(2, 1, 1),  # 2x + y + z <= 2
            c(1, -1, 3)) # x - y + 3z <= 3
b1 <- c(2, 3)
a <- c(4, 2, 9)

# simplex
result <- simplex(a = a, A1 = A1, b1 = b1, maxi = FALSE)

cat("最优解的值为:\n")
print(result$soln)
cat("目标函数的最小值为:\n")
print(result$value)

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())

formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
model_1 <- lapply(formulas, function(f) lm(f, data = mtcars))
print(model_1)


## -----------------------------------------------------------------------------

bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
fit_model <- function(data) {
  lm(mpg ~ disp, data = data)
}

model_2 <- lapply(bootstraps, fit_model)
print(model_2)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
r2_lapply <- sapply(c(model_1,model_2), rsq)

# 输出 R^2 值
print(r2_lapply)

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

p_values <- sapply(trials, function(x) x$p.value)
print(p_values)

## -----------------------------------------------------------------------------
p_values <- sapply(trials, `[[`, "p.value")
print(p_values)

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
parallel_lapply <- function(FUN, x,y, output_type) {
  results <- Map(FUN, x,y)
  vapply(results, identity, output_type)
}

a <- c(1, 2, 3)
b <- c(4, 5, 6)
output <- parallel_lapply(function(x, y) x + y, a, b, numeric(1))
print(output)

## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
fast_chisq_test <- function(x, y) {

  #列联表
  contingency_table <- table(x, y)
  row_sum <- rowSums(contingency_table)
  col_sum <- colSums(contingency_table)
  total_sample_size <- sum(contingency_table)
  
  #计算
  theoretical_frequency <- outer(row_sum, col_sum) / total_sample_size
  chisq_statistic <- sum((contingency_table - theoretical_frequency)^2 / theoretical_frequency)
  df <- (nrow(contingency_table) - 1) * (ncol(contingency_table) - 1)
  p_value <- pchisq(chisq_statistic, df, lower.tail = FALSE)

  return(list(p_value = p_value))
}

x <- c(1,2,3,1,2,3)
y <- c(3,2,1,3,2,1)
result <- fast_chisq_test(x, y)
print(result)


## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
fast_table <- function(x, y) {
  
  #获得观测值范围
  x_vals <- sort(unique(x))
  y_vals <- sort(unique(y))
  contingency_matrix <- matrix(0, nrow = length(x_vals), ncol = length(y_vals),
                               dimnames = list(x_vals, y_vals))
  
  x_indices <- match(x, x_vals)
  y_indices <- match(y, y_vals)
  for (i in seq_along(x)) {
    contingency_matrix[x_indices[i], y_indices[i]] <- contingency_matrix[x_indices[i], y_indices[i]] + 1
  }
  
  return(contingency_matrix)
}

fast_chisq_test <- function(x, y) {

  #列联表
  contingency_table <- fast_table(x, y)
  row_sum <- rowSums(contingency_table)
  col_sum <- colSums(contingency_table)
  total_sample_size <- sum(contingency_table)
  
  #计算
  theoretical_frequency <- outer(row_sum, col_sum) / total_sample_size
  chisq_statistic <- sum((contingency_table - theoretical_frequency)^2 / theoretical_frequency)
  df <- (nrow(contingency_table) - 1) * (ncol(contingency_table) - 1)
  p_value <- pchisq(chisq_statistic, df, lower.tail = FALSE)

  return(list(p_value = p_value))
}

x <- c(1,2,3,1,2,3)
y <- c(3,2,1,3,2,1)
result <- fast_chisq_test(x, y)
print(result)


## -----------------------------------------------------------------------------
#清空内存
rm(list = ls())
library(Rcpp)
library(ggplot2)

#Rcpp
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix gibbs_sampler(int n, int a, int b,  int N) {
    NumericMatrix samples(N, 2); // 存储结果
    double x = 0; // 初始化 x
    double y = 1; // 初始化 y

    for (int i = 0; i < N; ++i) {// Gibbs采样
        x = R::rbinom(n, y);
        y = R::rbeta(x + a, n - x + b);
        samples(i, 0) = x;
        samples(i, 1) = y;
    }
    return samples;
}
')
n <- 10       
a <- 2        
b <- 3        
N <- 1000     
#生成数据与绘图
set.seed(123)
samples <- gibbs_sampler(n, a, b, N)
plot(samples,xlab="x", ylab="y")

## -----------------------------------------------------------------------------

gibbs_sampler_r <- function(n, a, b, num_steps) {
  x <- 0
  y <- 1
  samples <- matrix(0, nrow = num_steps, ncol = 2)
  # Gibbs采样
  for (i in 1:num_steps) {
    x <- rbinom(1, size = n, prob = y)
    y <- rbeta(1, shape1 = x + a, shape2 = n - x + b)
    samples[i, ] <- c(x, y)
  }
  samples
}
set.seed(123)
samples_rcpp <- gibbs_sampler(n, a, b, N)
set.seed(123)
samples_r <- gibbs_sampler_r(n, a, b, N)

## -----------------------------------------------------------------------------
# 比较QQ图
qqplot(
  samples_r[,1], samples_rcpp[,1],
  main = "x QQ Plot",
  xlab = "OracleR_x",
  ylab = "Rcpp_x"
)
abline(0, 1, col = "blue")

# QQ plot comparison for y
qqplot(
  samples_r[,2], samples_rcpp[,2],
  main = "y QQ Plot",
  xlab = "OracleR_y",
  ylab = "Rcpp_y"
)
abline(0, 1, col = "blue")


## -----------------------------------------------------------------------------
library(microbenchmark)
#比较运行时间
benchmark_results <- microbenchmark(
  samples_rcpp <- gibbs_sampler(n, a, b, N),
  samples_r <- gibbs_sampler_r(n, a, b, N),
  times = 100 
)
print(benchmark_results)

