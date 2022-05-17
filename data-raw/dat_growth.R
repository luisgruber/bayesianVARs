## code to prepare `dat_growth` dataset goes here
library(tidyverse)
library(xts)
library(lubridate)
dat <- readr::read_csv("https://files.stlouisfed.org/files/htdocs/fred-md/quarterly/2021-07.csv")

type <- "growth" # "stationary"

# remove rows that only contain NAs
dat <- dat %>% filter(rowSums(is.na(dat)) != ncol(dat))
#dat <- dat[rowSums(is.na(dat)) != ncol(dat),] the same in base R
dates <- mdy(as.matrix(dat[-(1:2),1]))
data <- as.matrix(dat[-1, c("GDPC1","PCECC96","GPDIC1", "PRFIx", "INDPRO",
                            "CUMFNS","SRVPRD", "CE16OV", "AWHMAN", "PCECTPI",
                            "GDPCTPI", "GPDICTPI", "CPIAUCSL", "CES2000000008x",
                            "FEDFUNDS", "GS1", "GS10", "M2REAL", "EXUSUKx",
                            "UMCSENTx", "S&P 500")])

data_raw <- xts(data[-1,], dates)

# (approx) stationary data ---------------------------------------------------------
if(type=="stationary"){
# Transform variables according to McCracken & NG (2020): FRED-QD: A Quarterly Database for Macroeconomic Research (https://doi.org/10.20955/wp.2020.005)
data_trans <- matrix(NA, nrow(data)-1,ncol(data))
colnames(data_trans) <- colnames(data)
unique(data[1,]) # show transformations
#1: raw; 2: diff; 5: diff(log); 6: diff^2(log)
for (i in 1:ncol(data)) {
  if(data[1,i]==1) data_trans[,i] <- data[-1,i]
  if(data[1,i]==2) data_trans[2:nrow(data_trans),i] <- diff(data[-1,i])
  if(data[1,i]==5) data_trans[2:nrow(data_trans),i] <- diff(log(data[-1,i]))
  if(data[1,i]==6) data_trans[3:nrow(data_trans),i] <- diff(log(data[-1,i]), differences = 2)
}
# create xts object
data_trans <- xts(data_trans[-(1:3),], order.by = dates[-(1:3)])
small <- c("GDPC1", "CPIAUCSL", "FEDFUNDS")
medium <- c("GDPC1","PCECC96", "GPDIC1", "AWHMAN",  "CPIAUCSL", "CES2000000008x", "FEDFUNDS")
dat <- list(small=data_trans[,small], medium=data_trans[,medium], large=data_trans)
}

# Growth data: diff(log) --------------------------------------------------
if(type=="growth"){
data_growth <- matrix(NA, nrow(data)-1,ncol(data))
colnames(data_growth) <- colnames(data)
data[4,"UMCSENTx"] <- data[3,"UMCSENTx"]
for (i in colnames(data)) {
  if(i %in% c("FEDFUNDS","GS1","GS10")) data_growth[,i] <- data[-1,i]/100
  else data_growth[2:nrow(data_growth),i] <- diff(log(data[-1,i]))
}
# create xts object
data_growth <- xts(data_growth[-(1:3),], order.by = dates[-(1:3)])
small <- c("GDPC1", "CPIAUCSL", "FEDFUNDS")
medium <- c("GDPC1","PCECC96", "GPDIC1", "AWHMAN",  "CPIAUCSL", "CES2000000008x", "FEDFUNDS")
dat_growth <- list(small=data_growth[,small], medium=data_growth[,medium], large=data_growth)
}

usethis::use_data(dat_growth, overwrite = TRUE)
