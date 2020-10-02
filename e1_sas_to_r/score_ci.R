require(data.table)
data("GBSG2", package = "TH.data")

data <- as.data.table(GBSG2)
stratum <- NULL
time <- "time"
event <- "cens"
# Must only have two levels
treatment <- "horTh"
treat_lvl <- "yes"
b <- -0.3

# Convert treatment variable to indicator
data[,X:=fifelse(get(treatment)==treat_lvl,1,0)]

# Add exp vars
data[,`:=`(E_exp=exp(b*X), E_exp_x=exp(b*X)*X)]

# Remove ties
setorder(data, -time)
data_cum <- data[,.(.N, D=sum(get(event)), E_exp_x=sum(E_exp_x), E_exp=sum(E_exp)), by = time][
  ,.(time, D, R=cumsum(N), E=cumsum(E_exp_x)/cumsum(E_exp))]

# Scoring calculation
data_sc <- data[data_cum, on = "time", .(time, U=get(event)*(X-E), I=get(event)*(1-i.E)*i.E*(i.R-i.D)/(i.R-1))]
S_b <- data_sc[,.(U=sum(U), I=sum(I, na.rm = T), sc=sum(U)^2/sum(I, na.rm = T))]




