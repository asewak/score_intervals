require(data.table)
data("GBSG2", package = "TH.data")

data <- as.data.table(GBSG2)
stratum <- NULL
t <- "time"
event <- "cens"
# Must only have two levels
treatment <- "horTh"
treat_lvl <- "yes"
b <- -0.3

if(is.null(stratum)){
  data[,strata:=1]
  stratum <- "strata"
}

# Set column names to avoid using get
setnames(data, c(t, treatment, event), c("time","treatment","event"))

# Convert treatment variable to indicator
data[,X:=fifelse(treatment==treat_lvl,1,0)]

# Add exp vars
data[,`:=`(E_exp=exp(b*X), E_exp_x=exp(b*X)*X)]

# Remove ties
keys <- c("strata","time")
setorderv(data, keys, c(1,-1))
data_cum <- data[,.(.N, D=sum(event), E_exp_x=sum(E_exp_x), E_exp=sum(E_exp)), by = keys][
  ,.(strata, time, D, R=cumsum(N), E=cumsum(E_exp_x)/cumsum(E_exp))]

# Scoring calculation
data_sc <- data[data_cum, on = keys, .(time, U=event*(X-E), I=event*(1-i.E)*i.E*(i.R-i.D)/(i.R-1))]
S_b <- data_sc[!is.na(I)][,.(U=sum(U), I=sum(I), sc=sum(U)^2/sum(I))]




