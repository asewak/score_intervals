require(data.table)
require(survival)
data("GBSG2", package = "TH.data")

# data <- as.data.table(GBSG2)
# stratum <- NULL
# t <- "time"
# event <- "cens"
# treatment <- "horTh" # Must only have two levels
# treat_lvl <- "yes"
# b <- 0

fn_score <- function(dt, stratum=NULL, t, event, treatment, treat_lvl, b){
  
  # Convert to a data.table and create strata column
  data <- as.data.table(dt)
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
  return(S_b[,sc])
}


fn_cox <- function(dt, t, event, treatment){
  stime <- with(dt, Surv(time, cens))
  form <- as.formula(paste0("stime~",treatment))
  fit <- coxph(form, data = dt)
  sum_fit <- summary(fit)
  return(sum_fit$coefficients[,c("coef", "se(coef)")])
}

# fn_score(dt=GBSG2, stratum=NULL, t="time", event="cens", treatment="horTh", treat_lvl="yes", b=0)

est_cox <- fn_cox(dt=GBSG2, t="time", event="cens", treatment="horTh")

grid_b_low <- seq(est_cox["coef"] - 2.96*est_cox["se(coef)"], est_cox["coef"], length.out = 5)
grid_b_up <- seq(est_cox["coef"], est_cox["coef"] + 2.96*est_cox["se(coef)"], length.out = 5)


# Interpolate to find the intersection with the relevant
# Note: on the sqrt scale for speed
fn_interpolate <- function(g, alpha, dt, stratum=NULL, t, event, treatment, treat_lvl){
  S_grid <- sapply(g, function(b) fn_score(dt=dt, stratum=stratum, t=t, event=event, treatment=treatment, treat_lvl=treat_lvl, b=b))
  s <- spline(g, sqrt(S_grid))
  int <- approx(x = s$y, y = s$x, xout = qnorm(1-alpha/2))
  return(int)
}

fn_interpolate(g=grid_b_low, alpha=0.05, dt=GBSG2, stratum=NULL, t="time", event="cens", treatment="horTh", treat_lvl="yes")
fn_interpolate(g=grid_b_high, alpha=0.05, dt=GBSG2, stratum=NULL, t="time", event="cens", treatment="horTh", treat_lvl="yes")

