## code to prepare `dat_growth` dataset goes here
library(tidyverse)
library(xts)
library(lubridate)
dir <- paste0("")
dat <- readr::read_csv(paste0(dir,"2021-07.csv"))

# remove rows that only contain NAs
dat <- dat %>% filter(rowSums(is.na(dat)) != ncol(dat))
#dat <- dat[rowSums(is.na(dat)) != ncol(dat),] the same in base R
dates <- mdy(as.matrix(dat[-(1:2),1]))
data <- as.matrix(dat[-1, c("GDPC1","PCECC96","GPDIC1", "PRFIx", "INDPRO",
                            "CUMFNS","SRVPRD", "CE16OV", "AWHMAN", "PCECTPI",
                            "GDPCTPI", "GPDICTPI", "CPIAUCSL", "CES2000000008x",
                            "FEDFUNDS", "GS1", "GS10", "M2REAL", "EXUSUKx",
                            "UMCSENTx", "S&P 500")])


# (approx) stationary data ---------------------------------------------------------

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
data_final <- list(small=data_trans[,small], medium=data_trans[,medium], large=data_trans)
saveRDS(data_final, file = paste0(dir,"data_transformed_new.RData"))

# Growth data: diff(log) --------------------------------------------------

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
data_final <- list(small=data_growth[,small], medium=data_growth[,medium], large=data_growth)
saveRDS(data_final, file = paste0(dir,"data_growth.RData"))


# Data in log-levels ------------------------------------------------------

# (except interest rates)
data_llevel <- matrix(NA, nrow(data)-1,ncol(data))
colnames(data_llevel) <- colnames(data)
for (i in colnames(data)) {
  if(i == "FEDFUNDS") data_llevel[,i] <- data[-1,i]
  else if(i == "GS1") data_llevel[,i] <- data[-1,i]
  else if(i == "GS10") data_llevel[,i] <- data[-1,i]
  else data_llevel[,i] <- log(data[-1,i])
}
data_llevel <- xts(data_llevel[-(1:3),], order.by = dates[-(1:3)])
small <- c("GDPC1", "CPIAUCSL", "FEDFUNDS")
medium <- c("GDPC1","PCECC96", "GPDIC1", "AWHMAN",  "CPIAUCSL", "CES2000000008x", "FEDFUNDS")
datall_final <- list(small=data_llevel[,small], medium=data_llevel[,medium], large=data_llevel)
saveRDS(datall_final, file = paste0(dir,"data_loglevel.RData"))

usethis::use_data(dat_growth, overwrite = TRUE)
