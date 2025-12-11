# Time Series Analysis: Nominal Wages in the Czech Republic 🇨🇿

**Authors:** Mengzhen Shang, Liyuan Cao  
**Project:** Time Series Analysis (Project 1)

## 📌 Project Overview
This project performs a comprehensive time series analysis of Nominal Wages (W) in the Czech Republic (CZ). The analysis aims to model wage dynamics, identify structural dependencies on the Euro Area (EA20), and forecast future wage growth up to 2027.

The study compares univariate models (Random Walk, ARMA) against multivariate models (VAR) to determine the impact of external economic shocks from the Euro Area on the Czech labor market.

## 📂 Repository Structure
* **`Project 1(FE).R`**: The main R script containing all data preprocessing, statistical testing, modeling (ARMA/VAR), and visualization code.
* **`TIME SERIES ANALYSIS.pdf`**: The final presentation report summarizing findings, plots, and economic interpretations.

## 🛠 Methodology

### 1. Data Preprocessing & Stationarity
* **Transformation:** Log-difference transformation was applied to convert nominal wage levels into Quarter-over-Quarter (QoQ) percentage change.
* **Stationarity Checks:** Verified using ADF (Augmented Dickey-Fuller), PP (Phillips-Perron), and KPSS tests. The transformed series was confirmed as $I(0)$.

### 2. Univariate Modeling (ARMA)
* **Model Selection:** Based on AIC and BIC criteria.
* **Optimal Model:** ARMA(1,1) was selected.
* **Diagnostics:** Ljung-Box test confirmed residuals are White Noise; Jarque-Bera test checked for normality.

### 3. Multivariate Modeling (VAR)
* **External Variable:** Euro Area Nominal Wage Growth (EA20).
* **Lag Selection:** $p=3$ selected based on BIC/AIC/SC.
* **Structural Analysis:**
    * **Impulse Response Functions (IRF):** Analyzed the delayed impact of Euro Area shocks on CZ wages.
    * **FEVD:** Found that in the long run, ~25% of Czech wage volatility is driven by Euro Area shocks.
    * **Historical Decomposition:** Quantified the contribution of domestic vs. external shocks during the 2020-2025 period.

### 4. Forecasting
* **Evaluation:** Rolling window approach (last 8 quarters).
* **Metrics:** MFE, RMSFE, and Diebold-Mariano (DM) test.
* **Future Outlook (2025-2027):** Comparison of ARMA and VAR forecasts against European Commission (EC) projections.

## 📊 Key Results
1.  **Stationarity:** The original wage series is non-stationary; the QoQ% change is stationary.
2.  **Euro Area Impact:** There is a significant "catch-up" effect. Shocks from the Euro Area positively impact Czech wages with a lag of approx. 3 quarters.
3.  **Forecast Accuracy:** The Random Walk (RW) model performed best for short-term accuracy (lowest RMSFE), suggesting strong trend dominance in the recent volatile period.
4.  **Future Scenario:** The VAR model predicts higher wage growth ("optimistic scenario") compared to the more conservative ARMA and EC forecasts, driven by projected growth in the Euro Area.

## 💻 How to Run
1.  Ensure you have the required R packages installed:
    ```r
    install.packages(c("forecast", "urca", "zoo", "tseries", "knitr", "ggplot2", "vars"))
    ```
2.  **Prerequisite:** The script requires the data file `Block1Data.RData` to be in the working directory.
3.  Run `Project 1(FE).R` to generate all statistics and plots used in the report.

## 📜 License
This project is for educational and analytical purposes.
