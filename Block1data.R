# Creating our database
require(eurostat)
require(zoo)
require(lubridate)
library(openxlsx)
require(seasonal)
require(seasonalview)

# House price index
# browseURL("https://ec.europa.eu/eurostat/databrowser/view/prc_hpi_q__custom_14945099/default/table")
data <- get_eurostat("prc_hpi_q",
                     filters = list(purchase = "TOTAL",
                                    unit     = "I15_Q"
                     ))

HPI <- do.call(merge,lapply(unique(data$geo),function(cnt){
  ind  <- data$geo ==cnt
  zoo(data$values[ind],as.yearqtr(data$time[ind]))
}))
names(HPI)<-unique(data$geo) 

# expludes countries with more than 20 missing observations
colSums(is.na(HPI))<=20
HPI <- HPI[,colSums(is.na(HPI))<=20]

# these data require seasonal adjustment
plot(diff(log(HPI)))

HPISA <- do.call(merge,lapply(1:dim(HPI)[2],function(i){
  Y0    <- HPI[,i]
  m     <- seas(transform.function = "log", na.omit(ts(Y0,start=c(2005,1),freq=4)))
  Y0[!is.na(Y0)]   <- as.numeric(predict(m))
  Y0
}))
names(HPISA) <- names(HPI)
plot(diff(log(HPISA)))
HPI <- HPISA

# Consumption deflator and real consumption
# browseURL("https://ec.europa.eu/eurostat/databrowser/view/namq_10_gdp__custom_14945120/default/table")
data <- get_eurostat("namq_10_gdp",
                     filters = list(na_item = "P31_S14_S15",
                                    unit   = c("CLV20_MNAC","PD20_NAC"),
                                    s_adj  = "SCA"
                     ))


C <- do.call(merge,lapply(unique(data$geo),function(cnt){
  ind  <- data$geo ==cnt & data$unit == "CLV20_MNAC"
  zoo(data$values[ind],as.yearqtr(data$time[ind]))
}))
names(C)<-unique(data$geo) 
C <- window(C, start=as.yearqtr("1995 Q1"))

P <- do.call(merge,lapply(unique(data$geo),function(cnt){
  ind  <- data$geo ==cnt & data$unit == "PD20_NAC"
  zoo(data$values[ind],as.yearqtr(data$time[ind]))
}))
names(P)<-unique(data$geo) 
P <- window(P, start=as.yearqtr("1995 Q1"))

# Disposable income
# browseURL("https://ec.europa.eu/eurostat/databrowser/view/nasq_10_nf_tr__custom_14945133/default/table")
data <- get_eurostat("nasq_10_nf_tr",
                     filters = list(na_item = "B6G",
                                    unit   = c("CP_MNAC"),
                                    s_adj  = "SCA",
                                    sector = "S14_S15",
                                    direct = "RECV"
                     ))

DI <- do.call(merge,lapply(unique(data$geo),function(cnt){
  ind  <- data$geo ==cnt 
  zoo(data$values[ind],as.yearqtr(data$time[ind]))
}))
names(DI)<-unique(data$geo)
DI <- window(DI, start=as.yearqtr("1999 Q1"))
DI <- DI[,colSums(is.na(DI))<=40]


# LABOUR COSTS -- WAGES
# browseURL("https://ec.europa.eu/eurostat/databrowser/view/lc_lci_r2_q__custom_14945168/default/table")
data <- get_eurostat("lc_lci_r2_q",
                     filters = list(lcstruct = "D11",
                                    unit   = c("I20"),
                                    s_adj  = "SCA",
                                    nace_r2 = "B-S"
                     ))

W <- do.call(merge,lapply(unique(data$geo),function(cnt){
  ind  <- data$geo ==cnt 
  zoo(data$values[ind],as.yearqtr(data$time[ind]))
}))
names(W)<-unique(data$geo)
W <- window(W, start=as.yearqtr("1999 Q1"))
W <- W[,colSums(is.na(W))<=40]

# should be seasonally adjusted, but...
plot(diff(log(W)))


# Saving all data in R
# countries with all data
cntList <- Reduce(intersect, list(names(HPI),names(C),names(P),names(DI),names(W)))
dt <-list(HPI  = HPI[,cntList],
          RHPI = HPI[,cntList]/P[,cntList],
          C    = C[,cntList],
          P    = P[,cntList],
          DI   = DI[,cntList],
          RDI  = DI[,cntList]/P[,cntList],
          W    = W[,cntList],
          RW   = W[,cntList]/P[,cntList])

save(dt, file = "Block1Data.RData")

# Saving to excel data -- the same time sample

idx <- Reduce(intersect, list(index(HPI),index(C),index(P),index(DI),index(W)))
idx <- as.yearqtr(idx)

wb <- createWorkbook("Block1")
addWorksheet(wb, "C");   writeData(wb, sheet = "C", data.frame(dates=as.character(idx),C[idx,cntList]))
addWorksheet(wb, "HPI"); writeData(wb, sheet = "HPI", data.frame(dates=as.character(idx),HPI[idx,cntList]))
addWorksheet(wb, "DI");  writeData(wb, sheet = "DI", data.frame(dates=as.character(idx),DI[idx,cntList]))
addWorksheet(wb, "W");   writeData(wb, sheet = "W", data.frame(dates=as.character(idx),W[idx,cntList]))
addWorksheet(wb, "P");   writeData(wb, sheet = "P", data.frame(dates=as.character(idx),P[idx,cntList]))
saveWorkbook(wb, "Block1data.xlsx", overwrite = TRUE)





