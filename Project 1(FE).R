# 0. 环境设置与数据加载
##########################################
rm(list=ls())
require(forecast)
require(urca)
require(zoo)
require(tseries)
require(knitr)
require(ggplot2)

load("Block1Data.RData")

# 选择变量和国家
vnm <- "W"   # 名义工资 (Nominal Wages)
cnt <- "CZ"  # 捷克 (Czechia)
data <- dt[[vnm]]
Y <- data[,cnt]

# 转换: 4 -- q/q% (季度环比百分比变化率)
# 这是我们用来建模的平稳序列 x
x <- 100 * diff(log(Y), 1) 
vNam <- paste0(vnm, " in ", cnt, " (QoQ %)")

# 数据对齐与清理
start.date <- as.yearqtr("1980 Q1") 
x <- na.omit(x)
x <- na.omit(window(x, start=start.date)) 
Y_orig <- na.omit(window(Y, start=start.date)) # 用于原始序列的平稳性检验

cat("数据准备完成。建模序列名称:", vNam, "\n")


#a
# 1. 确保必要的包已加载
require(zoo)

# 变量 Y_orig 应该已经在数据准备部分被定义为：
# Y_orig <- na.omit(window(Y, start=as.yearqtr("1980 Q1")))

# 2. 绘制原始序列图
# 设置绘图参数，以确保图表清晰专业
par(mfrow=c(1,1), cex = 0.8, bty="l") 

# 绘制时间序列图
plot(Y_orig, 
     main="Original Series: Nominal Wages in CZ (Level)", # 主标题
     sub="", 
     xlab="Year", 
     ylab="Nominal Wages Level",
     col="blue", # 设置线条颜色
     lwd=2)      # 设置线条粗细

# 可选：添加一条趋势线来强调非平稳性
# T <- 1:length(Y_orig)
# fit <- lm(Y_orig ~ T)
# abline(fit, col="red", lty=2) # 绘制线性趋势线

# 恢复默认绘图布局（如果需要继续绘制多图）
# par(mfrow=c(2,2)) 

# 打印信息
cat("原始序列 Y_orig 的时间序列图已生成。\n")


print(start(x))
print(length(x))


# 2.1. 图形化分析
par(mfrow=c(2,2), cex = 0.8, bty="l")
# 时间序列图：观察均值和波动
plot(x, main=vNam, sub="", xlab="", ylab=""); abline(h=mean(x),col="red")
# 直方图：观察分布形状
hist(x, main="Histogram", ylab="", sub="", xlab="", probability = TRUE)
lines(density(x), col = "blue", lwd = 2) 
# ACF和PACF：判断自相关结构和潜在的AR/MA阶数
Acf(x, 16, main="ACF (QoQ %)", ylab="", sub="", xlab="")
Pacf(x, 16, main="PACF (QoQ %)", ylab ="", sub="", xlab="")

# 关闭图形设备以准备下一个 ggplot 输出
dev.off()

cat("图形分析完成。请查看绘制的四张图。\n")

# 2.2. 单位根检验

#stationary check（需要输出在excel）
library(urca)

# ADF（Augmented Dickey-Fuller）单位根检验
summary(ur.df(x, type = "drift", lags = 4, selectlags = "BIC"))  #diff之后的数据
summary(ur.df(na.omit(Y), type = "drift", lags = 4, selectlags = "BIC"))   #原数据

# PP（Phillips-Perron）检验
summary(ur.pp(x, type = "Z-tau", model = "constant"))
summary(ur.pp(na.omit(Y), type = "Z-tau", model = "constant"))

# KPSS（Kwiatkowski–Phillips–Schmidt–Shin）检验
summary(ur.kpss(x, type = "mu", lags = "short"))
summary(ur.kpss(na.omit(Y), type = "mu", lags = "short"))


#b
# 1. 模型选择 (Model Selection)
cat("\n--- B.1 模型选择: BIC & AIC 准则 ---\n")

# LagSel 函数用于计算不同 p, q 组合下的信息准则
LagSel <- function(x, Pmax=4, Qmax=2, d=0, crit="BIC"){
  IC <- matrix(NA, Pmax+1, Qmax+1)
  for(p in 0:Pmax){
    for(q in 0:Qmax){
      tryCatch({
        # 使用 ML (Maximum Likelihood) 方法估计
        model <- Arima(x, order=c(p,d,q), method="ML")
        n <- length(x)
        if(crit == "AIC"){ IC[p+1,q+1] <- AIC(model, k=2) }
        if(crit == "BIC"){ IC[p+1,q+1] <- AIC(model, k=log(n)) }
      }, error = function(e){ IC[p+1,q+1] <- NA })
    }
  }
  rownames(IC) <- paste('ar', 0:Pmax, sep=""); colnames(IC) <- paste('ma', 0:Qmax, sep="")
  return(IC)
}

# 输出 BIC 表格 (用于演示文稿截图)
BIC_table <- LagSel(x, Pmax=4, Qmax=2, crit="BIC")
cat("\n[B 部分输出 1: BIC 准则表格]\n")
print(kable(BIC_table, format = "pipe", digits = 2))

# 输出 AIC 表格 (用于参考)
AIC_table <- LagSel(x, Pmax=4, Qmax=2, crit="AIC")
cat("\n[B 部分输出 2: AIC 准则表格]\n")
print(kable(AIC_table, format = "pipe", digits = 2))

# 确定最优模型（请根据您的实际结果替换 c(p, 0, q)）
# 假设 BIC 选择了 ARMA(1, 0, 1) 为最优模型（示例）
pdq_best <- c(1, 0, 1) 
arma_best <- Arima(x, order = pdq_best, method = "ML")
cat(paste("\n**已选定最优模型为 ARMA(", pdq_best[1], ", 0, ", pdq_best[3], ")**\n"))

# 打印模型系数和统计量 (用于参考和诊断)
summary(arma_best)


# 2. 模型诊断: Ljung-Box Test (残差白噪声检验)
cat("\n--- B.2 模型诊断: Ljung-Box Test ---\n")

# Ljung-Box Test for Autocorrelation (White Noise Check)
k <- length(arma_best[["coef"]]) - 1 # 自由度 = 模型参数个数 - 1
Jmax <- 12 # 检验到 12 阶滞后 (3年)
LB <- matrix(NA, Jmax-k, 2)

for(h in (k+1):Jmax){
  # H0: 残差是白噪声 (White Noise)
  z <- Box.test(residuals(arma_best), lag = h, type = "Ljung-Box", fitdf = k)
  LB[h-k,] <- c(z[["statistic"]], z[["p.value"]])
}
rownames(LB) <- paste('h=', (k+1):Jmax, sep=""); colnames(LB) <- c("LB Stat.", "P-value")

# 打印 Ljung-Box 表格 (用于演示文稿截图)
cat("\n[B 部分输出 3: Ljung-Box 检验表格 (残差白噪声)]\n")
print(kable(LB, format = "pipe", digits = 3)) 
# 决策：如果所有 P-value > 0.05，则不拒绝 H0，残差为白噪声。

# 确保 ggplot2 包已加载
require(ggplot2)

# --- 1. 数据准备 ---
# 计算标准化残差 (Standardized Residuals)
# e: (残差 - 残差均值) / 残差标准差
# 在 R 中，残差均值通常接近于零，所以我们使用 residuals(arma_best) / sd(residuals(arma_best))
e <- residuals(arma_best) / sd(residuals(arma_best))

# 定义直方图的 bin 宽度 (可根据您的数据调整，0.2 是一个好的起点)
bwdth <- 0.2

# --- 2. 绘制图形 ---
# 注意：stat_function 中的 dnorm(x) 默认是标准正态分布 (均值0, 方差1)
# 乘以 length(e) * bwdth 是为了将理论密度曲线的高度缩放到与直方图的计数相匹配
Normality_Check_Plot <- ggplot(data.frame(e), aes(x = e)) +
  theme_light() +
  
  # 绘制直方图
  geom_histogram(binwidth = bwdth, colour = "white", fill = "orange", alpha = 0.8) +
  
  # 叠加正态分布曲线（理论频率）
  stat_function(fun = function(x) dnorm(x) * length(e) * bwdth, 
                color = "red", 
                linewidth = 1.2) +
  
  # 标签和标题
  labs(title="Model Residuals vs. Normal Distribution (Normality Check)", 
       y="Frequency", 
       x="Standardized Residuals")

# 打印图表（用于 PPT 演示文稿）
print(Normality_Check_Plot)
cat("残差直方图已生成，用于视觉检查正态性。\n")

# 确保 tseries 包已加载
require(tseries)

cat("\n--- Jarque-Bera Test for Residual Normality ---\n")

# 提取 ARMA 模型 (arma_best) 的残差
residuals_arma <- residuals(arma_best)

# 执行 Jarque-Bera 检验
jb_test_result <- jarque.bera.test(residuals_arma)

# 打印关键结果
print(jb_test_result)

# 提取关键值
JB_stat <- round(jb_test_result$statistic, 3)
JB_p_value <- round(jb_test_result$p.value, 3)

cat(paste("\nJB Statistic:", JB_stat, "\n"))
cat(paste("P-value:", JB_p_value, "\n"))

# 3. 脉冲响应函数 (IRF)
ar_coef <- arma_best[["coef"]][grep("^ar", names(arma_best[["coef"]]))]
ma_coef <- arma_best[["coef"]][grep("^ma", names(arma_best[["coef"]]))]
Kmax = 16 # 16 个季度 (4 年)

# ARMAtoMA 计算 MA 系数 (即 IRF)
IRF_coeffs <- ARMAtoMA(ar = ar_coef, ma = ma_coef, Kmax)
IRF_data <- data.frame(
  H = 0:Kmax,
  IRF = c(1, IRF_coeffs) # IRF at H=0 is always 1
)

# 绘制 IRF 图 (使用 ggplot)
IRF_plot <- ggplot(IRF_data, aes(x = H, y = IRF)) +
  geom_line(color = "orange", linewidth = 1) +
  geom_point(color = "orange") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title=paste0("IRF from ARMA(", pdq_best[1], ", 0, ", pdq_best[3], ") for ", cnt), 
       x="Quarters after Shock (H)", y="Response of QoQ % Wage Growth")

# 打印 IRF 图 (用于演示文稿)
print(IRF_plot)
cat("\n[B 部分输出 4: ARMA IRF 图已生成]\n")

#c

# 确保已安装并加载 vars 包
require(vars) 
require(zoo) 
require(forecast) 
require(ggplot2)

# --- 1. 数据准备 (Data Preparation for Euro Area) ---
# 变量: 欧元区名义工资 (假设列名为 EA20)
vnm_euro <- "W"
Y_euro <- dt[[vnm_euro]][,"EA20"] # 假设 EA20 是欧元区

# 应用相同的 QoQ % 转换到欧元区数据
x_euro <- 100 * diff(log(Y_euro), 1)

# 数据清理，确保从同一日期开始
start.date <- as.yearqtr("1980 Q1")
#start.date <- start(x) # 使用 CZ 序列的起始日期对齐
x_euro <- na.omit(window(x_euro, start=start.date, end=end(x))) # 对齐时间范围

# --- 2. 数据对齐与向量创建 (Alignment and Vector Creation) ---

# 合并两个平稳序列，确保它们的时间索引完全一致
Z <- merge(x_euro, x) 
Z <- na.omit(Z) # 清理任何合并后出现的 NA

# 重新命名列，以在输出中清晰识别 (y*_t, y_t)
colnames(Z) <- c("W_EA", "W_CZ") 

cat("VAR 建模数据已准备。样本量:", nrow(Z), "\n")

# --- 3. 确定最优滞后阶数 (Optimal Lag Selection) ---

# 使用 VARselect 确定最优滞后阶数 (p)
# Pmax=8 通常适用于季度数据 (2年)
var_select <- VARselect(Z, lag.max = 8, type = "const") 

# 打印信息准则，用于确定最佳 P 阶数
cat("\n--- C.1 滞后阶数选择结果 (VARselect) ---\n")
print(var_select$criteria)

# 假设 BIC 选择了最优滞后阶数 p_opt (请根据您的实际结果替换)
p_opt <- var_select$selection["BIC"] 
cat(paste("\nBIC 选定的最优滞后阶数 p_opt =", p_opt, "\n"))

# 1. 设定滞后阶数 (已经完成)
p_opt <- 3

# 2. 估计 VAR 模型 (缺失的步骤，必须运行)
# 假设 Z 变量（合并后的数据）存在
var_model <- VAR(Z, p = p_opt, type = "const") 

# 3. 脉冲响应函数 (IRF)
var_irf <- irf(var_model, impulse = c("W_EA", "W_CZ"), response = c("W_EA", "W_CZ"), 
               n.ahead = 16, ortho = TRUE, cumulative = FALSE)

# 绘图（检查是否成功）
plot(var_irf)

# --- 1. 估计 VAR 模型 (VAR Estimation) ---

# 手动设定最优滞后阶数 p_opt 为 3
p_opt <- 3 

# 假设 Z 变量（包含 W_EA 和 W_CZ 的合并数据）已经定义
# var_model <- VAR(Z, p = p_opt, type = "const") 
# 如果模型已运行成功，跳过此行，使用已有的 var_model 对象

# --- 2. 结构分析 (Structural Analysis) ---

# 2.1. 脉冲响应函数 (IRF)
# ortho=TRUE 使用 Cholesky 分解
var_irf <- irf(var_model, impulse = c("W_EA", "W_CZ"), response = c("W_EA", "W_CZ"), 
               n.ahead = 16, ortho = TRUE, cumulative = FALSE)
cat("\n--- C.3 VAR 脉冲响应函数 (IRF) --- (图表输出)\n")
plot(var_irf) 


# 2.2. 预测误差方差分解 (FEVD)
var_fevd <- fevd(var_model, n.ahead = 16)
cat("\n--- C.4 FEVD 输出 (用于 PPT 表格) ---\n")
print(var_fevd) 

# -----------------------------------------------------
# C. VAR 模型估计与结构分析 (无冗余精简版)
# -----------------------------------------------------

# 1. 设定最优滞后阶数 p = 3
p_opt <- 3 

# 2. 估计 VAR 模型 (只需运行一次)
cat("\n--- C.2 估计 VAR(3) 模型 ---\n")
var_model <- VAR(Z, p = p_opt, type = "const") 
# 确保模型对象 var_model 被创建并存储

# 3. 脉冲响应函数 (IRF) - Cholesky 分解
cat("\n--- C.3 VAR 脉冲响应函数 (IRF) --- (图表输出)\n")
var_irf <- irf(var_model, 
               impulse = c("W_EA", "W_CZ"), 
               response = c("W_EA", "W_CZ"), 
               n.ahead = 16, 
               ortho = TRUE, 
               cumulative = FALSE)
plot(var_irf) # 输出四张 IRF 图


# 4. 预测误差方差分解 (FEVD)
cat("\n--- C.4 FEVD 输出 (用于 PPT 表格) ---\n")
var_fevd <- fevd(var_model, n.ahead = 16)
print(var_fevd) 

# --- (接续 HD 和 Part D) ---

# 假设 var_model 对象已成功创建
# -----------------------------------------------------

# 1. 初始化变量和参数
i <- 2                             # 研究变量 W_CZ 在 Z 中的列索引 (第二列)
VARE <- var_model                  # 统一模型名称
e <- residuals(VARE)               # 提取 VAR 残差 e (reduced-form residuals)

# 2. 获取结构性冲击 (Structural Shocks, u)
# 使用 Cholesky 分解的 B 矩阵。
B <- t(chol(cov(e)))             # B 是下三角矩阵，从残差协方差矩阵 Sigma 得到
u <- t(solve(B) %*% t(e))        # 计算结构性冲击 u = B^{-1} * e
T_obs <- dim(u)[1]               # 观测数量
colnames(u) <- c("u_WEA", "u_WCZ") # 命名冲击：u_WEA (Euro Area Shock), u_WCZ (CZ Shock)

# 3. 计算 IRFs (结构性 SVMA 矩阵)
# 使用 vars 包的 Phi 函数来计算 IRF (即 SVMA 矩阵)
SVMA <- Phi(VARE, nstep=T_obs)
SIRF <- t(SVMA[i,,])             # SIRF 针对变量 i (W_CZ)

# 4. 历史分解 (HistDec) 核心循环
HistDec <- matrix(0, T_obs, 2)
colnames(HistDec) <- colnames(u)

for(t in 1:T_obs){
  # tmp1: 过去的冲击 (从 t 到 1)
  tmp1 <- as.matrix(u[t:1,]) 
  # tmp2: 对应的 IRF 响应 (从 t 到 1)
  tmp2 <- as.matrix(SIRF[1:t,])
  # 卷积求和 (Cumulative effect of all past shocks)
  HistDec[t,] <- colSums(tmp1 * tmp2)
}

# 5. 初始条件和均值项 (InitCond)
# InitCond 包含了初始条件和确定性项（截距）的影响
z_i <- tail(Z[,i], T_obs) 
InitCond <- z_i - rowSums(HistDec)
# 分离均值 (mu)
mu <- mean(z_i - rowSums(HistDec))
InitCond <- InitCond - mu

# 6. 最终 HD 数据整合
HistDec <- zoo(HistDec, order.by = index(InitCond))
HistDec <- merge(HistDec, InitCond = InitCond, Mean = mu) # 加上初始条件和均值项

# 7. 绘图 (Output)
par(mfrow=c(1,1), cex = 0.7, bty="l")
# 仅绘制冲击项 (u_WEA, u_WCZ, InitCond)
tmp <- barplot(t(HistDec[,1:3]), # 不包含 Mean
               main="Historical Decomposition of W_CZ (QoQ %)", 
               legend.text=c("Euro Area Shock", "CZ Shock", "Initial Conditions"), 
               xlab="", 
               col = c("blue","red", "grey"), 
               border=NA, # 移除柱子边框
               space = 0, # 移除柱子间隙
               las=2)

# 叠加实际值 (Actual Value) 的曲线
lines(x=tmp, y=z_i, lwd=3, col="black")

cat("\n✅ 历史分解图已生成。请检查图表以完成 Part C 分析。\n")

#d
# ==============================================================================
# 0. 数据加载与预处理 (Data Loading & Setup)
# ==============================================================================

# --- 配置区域 ---
vnm <- "W"        # 变量: 名义工资 (Nominal Wages)
cnt_target <- "CZ"   # 目标国家: 捷克
cnt_comp   <- "EA20" # 对比区域: 欧元区 (作为 VAR 的第二个变量)

# --- 数据提取 ---
# 从 W 数据集中分别提取 CZ 和 EA20
# 注意: dt[[vnm]] 是一个 mts (多变量时间序列) 对象
raw_target <- dt[[vnm]][, cnt_target]
raw_comp   <- dt[[vnm]][, cnt_comp]

# --- 数据转换 (QoQ Growth Rate) ---
# 转换为对数差分 (近似百分比增长率)
# 公式: 100 * (log(Y_t) - log(Y_{t-1}))
growth_target <- 100 * diff(log(raw_target), 1)
growth_comp   <- 100 * diff(log(raw_comp), 1)

# 合并为 VAR 所需的矩阵 (去除 NA 以对齐)
data_matrix <- cbind(growth_target, growth_comp)
colnames(data_matrix) <- c(cnt_target, cnt_comp)
data_matrix <- na.omit(data_matrix)

# 定义我们的主要分析对象 x (捷克的增长率)
x <- data_matrix[, cnt_target]

# 获取原始水平值的最后部分，用于后续还原
last_obs_level <- tail(raw_target, 1)
last_date_raw  <- time(raw_target)[length(raw_target)]

# ==============================================================================
# Task d: 模型精度比较 (RW, ARMA, VAR)
# Compare accuracy: MFE, RMSFE, DM test, sequential forecasts graph
# ==============================================================================

h_test <- 8 # 测试集长度 (2年)
n_total <- nrow(data_matrix)
n_train <- n_total - h_test

# 存储结果
actuals <- x[(n_train + 1):n_total]
fc_rw <- numeric(h_test)
fc_arma <- numeric(h_test)
fc_var <- numeric(h_test)

# 滚动预测循环 (Rolling Forecast)
for(i in 1:h_test) {
  curr_end <- n_train + i - 1
  
  # 当前训练窗口
  train_mat <- window(data_matrix, end = time(data_matrix)[curr_end])
  train_x   <- train_mat[, cnt_target]
  
  # 1. RW (Random Walk with Drift)
  # 因为名义工资通常有正增长趋势，带漂移项的 RW 更合理
  fit_rw <- rwf(train_x, h=1, drift=TRUE) 
  fc_rw[i] <- fit_rw$mean[1]
  
  # 2. ARMA
  # 自动选择最佳 ARIMA(p,0,q)
  fit_arma <- auto.arima(train_x, d=0, stationary=TRUE, seasonal=FALSE)
  fc_arma[i] <- forecast(fit_arma, h=1)$mean[1]
  
  # 3. VAR (CZ 和 EA20 的双变量模型)
  # 选择滞后阶数 (AIC)
  lag_sel <- VARselect(train_mat, lag.max=4, type="const")
  k_opt <- lag_sel$selection["AIC(n)"]
  if(is.na(k_opt)) k_opt <- 1
  
  fit_var <- VAR(train_mat, p=k_opt, type="const")
  pred_var <- predict(fit_var, n.ahead=1)
  fc_var[i] <- pred_var$fcst[[cnt_target]][1, "fcst"]
}

# --- 4.1 计算指标 (MFE, RMSFE) ---
calc_metrics <- function(act, pred) {
  err <- act - pred
  c(MFE = mean(err), RMSFE = sqrt(mean(err^2)))
}

# 1. 计算指标
res_table <- rbind(
  RW   = calc_metrics(actuals, fc_rw),
  ARMA = calc_metrics(actuals, fc_arma),
  VAR  = calc_metrics(actuals, fc_var)
)

# 2. 确保是 data.frame
res_table <- as.data.frame(res_table)

# 3. 计算相对 RMSFE
res_table$Rel_RMSFE <- round(res_table$RMSFE / res_table["RW", "RMSFE"], 4)

# 4. 打印最终结果
print("=== Task d: Forecast Accuracy Metrics (with Relative RMSFE) ===")
print(res_table)

# --- 4.2 DM Test (Diebold-Mariano) ---
e_rw   <- actuals - fc_rw
e_arma <- actuals - fc_arma
e_var  <- actuals - fc_var

cat("\n=== Task d: DM Test P-values ===\n")
cat("RW vs ARMA: ", dm.test(e_rw, e_arma, h=1, power=2)$p.value, "\n")
cat("RW vs VAR:  ", dm.test(e_rw, e_var, h=1, power=2)$p.value, "\n")
cat("ARMA vs VAR:", dm.test(e_arma, e_var, h=1, power=2)$p.value, "\n")

# ==============================================================================
# 修正后的 Task d 绘图代码
# ==============================================================================

# 1. 提取纯数值，避免索引冲突
y_actual <- as.numeric(actuals)
y_rw     <- as.numeric(fc_rw)
y_arma   <- as.numeric(fc_arma)
y_var    <- as.numeric(fc_var)

# 2. 获取时间标签 (用于 x 轴)

date_labels <- time(actuals) 

# 3. 绘图
# ylim 使用纯数值计算范围，不仅不会报错，还更准确
plot(date_labels, y_actual, type="o", lwd=2, pch=16, 
     ylim=range(c(y_actual, y_rw, y_arma, y_var)),
     main=paste("Task d: Sequential Forecasts (", cnt_target, " W Growth %)", sep=""),
     ylab="QoQ Growth %", xlab="Time")

# 4. 添加预测线
lines(date_labels, y_rw,   col="green", lty=2, lwd=2, type="o", pch=1)
lines(date_labels, y_arma, col="blue",  lty=2, lwd=2, type="o", pch=2)
lines(date_labels, y_var,  col="red",   lty=2, lwd=2, type="o", pch=3)

# 5. 添加图例
legend("top", legend=c("Actual", "RW", "ARMA", "VAR"),
       col=c("black", "green", "blue", "red"), 
       lty=c(1, 2, 2, 2), pch=c(16, 1, 2, 3), lwd=2)

grid() 
# ==============================================================================
# Task e: 未来两年预测 (修正版: 精确匹配到 2027 Q2)
# ==============================================================================

h_future <- 8 # 预测步长: 8个季度 (2025 Q3 - 2027 Q2)

# --- 1. 模型拟合与预测 (增长率) ---
# ARMA (全样本)
fit_arma_full <- auto.arima(x, d=0, stationary=TRUE, seasonal=FALSE)
fc_arma_growth <- forecast(fit_arma_full, h=h_future)$mean

# VAR (全样本)
lag_sel_full <- VARselect(data_matrix, lag.max=4, type="const")
k_opt_full <- lag_sel_full$selection["AIC(n)"]
if(is.na(k_opt_full)) k_opt_full <- 1 # 兜底防止NA
fit_var_full <- VAR(data_matrix, p=k_opt_full, type="const")
pred_var_full <- predict(fit_var_full, n.ahead=h_future)
fc_var_growth <- pred_var_full$fcst[[cnt_target]][, "fcst"]

# --- 2. 还原为水平值 (Levels) ---
# 辅助函数: 从最后一个观测值开始，累乘增长率
reconstruct_levels <- function(last_val, rates) {
  lvls <- numeric(length(rates))
  curr <- last_val
  for(k in 1:length(rates)) {
    # rate 是 100 * log_diff -> exp(rate/100) 是乘数
    curr <- curr * exp(rates[k]/100)
    lvls[k] <- curr
  }
  # 创建时间序列对象 (自动推算开始时间为历史数据结束后的下一个季度)
  start_time <- last_date_raw + 0.25 
  ts(lvls, start=start_time, frequency=4)
}

lvl_arma <- reconstruct_levels(last_obs_level, fc_arma_growth)
lvl_var  <- reconstruct_levels(last_obs_level, fc_var_growth)

# --- 3. European Commission (EC) Forecast (精确年份匹配) ---
# 数据来源: European Commission (Nominal wages growth)
# 2025: 5.9%, 2026: 5.4%, 2027: 4.8%
# 转换为季度复合增长率: (1 + r)^(1/4)

r_2025 <- 0.059
r_2026 <- 0.054
r_2027 <- 0.048

g_q_2025 <- (1 + r_2025)^(1/4)
g_q_2026 <- (1 + r_2026)^(1/4)
g_q_2027 <- (1 + r_2027)^(1/4)

# *** 关键修正 ***
# 历史数据截止于 2025 Q2。
# 预测期覆盖: 
# 2025: Q3, Q4 (2个季度)
# 2026: Q1, Q2, Q3, Q4 (4个季度)
# 2027: Q1，Q2 (2个季度)
factors <- c(rep(g_q_2025, 2), rep(g_q_2026, 4), rep(g_q_2027, 2))

lvl_ec <- numeric(h_future)
curr_ec <- last_obs_level
for(k in 1:h_future) {
  curr_ec <- curr_ec * factors[k]
  lvl_ec[k] <- curr_ec
}
lvl_ec <- ts(lvl_ec, start=start(lvl_arma), frequency=4)

# --- 4. 结果汇总表 ---
# 生成对应的时间标签 (例如 "2025 Q2" 等)
future_dates <- as.yearqtr(time(lvl_arma))

df_future <- data.frame(
  Date = future_dates,
  ARMA = as.numeric(lvl_arma),
  VAR  = as.numeric(lvl_var),
  EC   = as.numeric(lvl_ec)
)

print("=== Task e: Future Forecasts Table (Levels) ===")
print(df_future)

# --- 5. 最终绘图 ---

# 1. 截取最近 3 年的历史数据
hist_plot <- window(raw_target, start=time(raw_target)[length(raw_target)-11])

# 2. 计算范围 (全部转为纯数字，防止报错)
y_vals_all <- c(as.numeric(hist_plot), as.numeric(lvl_arma), as.numeric(lvl_var), as.numeric(lvl_ec))
y_lims <- range(y_vals_all)

# 3. 确定 X 轴范围 (关键修改：加 as.numeric)
# 将时间对象 (yearqtr) 强制转换为纯小数年份 
x_start <- as.numeric(time(hist_plot)[1])
x_end   <- as.numeric(time(lvl_arma)[h_future])

# 右侧留白 0.5 以便放下 "2027"
x_limits <- c(x_start, x_end + 0.5) 

# 4. 绘图
plot(hist_plot, type="o", pch=20, col="black", lwd=2,
     xlim=x_limits, ylim=y_lims, xaxt="n", # 禁止自动 X 轴
     main=paste("Task e: Forecast of Nominal Wages (", cnt_target, ") to 2027 Q2", sep=""),
     ylab="Index Level", xlab="Year")

# 5. 手动添加年份标签
years_to_label <- seq(floor(x_start), ceiling(x_end))

# 画主刻度 (年份)
axis(1, at=years_to_label, labels=years_to_label, las=1)
# 画次刻度 (季度)
axis(1, at=seq(floor(x_start), ceiling(x_end), by=0.25), labels=FALSE, tcl=-0.2)

# 6. 添加预测线
lines(lvl_arma, col="blue", lwd=2, type="o", pch=1)
lines(lvl_var,  col="red",  lwd=2, type="o", pch=2)
lines(lvl_ec,   col="darkgreen", lwd=2, type="o", pch=3)

# 7. 图例
legend("topleft", 
       legend=c("History", "ARMA", "VAR (w/ EA20)", "EC Forecast"),
       col=c("black", "blue", "red", "darkgreen"),
       lty=1, pch=c(20, 1, 2, 3), lwd=2)

grid()

