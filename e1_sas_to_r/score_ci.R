require(data.table)
require(survival)
data("GBSG2", package = "TH.data")

# Score as a function of beta - S(b)
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

# Cox PH estimate
fn_cox <- function(dt, t, event, treatment){
  stime <- with(dt, Surv(time, cens))
  form <- as.formula(paste0("stime~",treatment))
  fit <- coxph(form, data = dt)
  sum_fit <- summary(fit)
  return(sum_fit$coefficients[,c("coef", "se(coef)")])
}

# Interpolate to find the intersection with the relevant quantile
# Note: on the sqrt scale for speed
fn_interpolate <- function(g, alpha, dt, stratum=NULL, t, event, treatment, treat_lvl){
  S_grid <- sapply(g, function(b) fn_score(dt=dt, stratum=stratum, t=t, event=event, treatment=treatment, treat_lvl=treat_lvl, b=b))
  s <- spline(g, sqrt(S_grid))
  int <- approx(x = s$y, y = s$x, xout = qnorm(1-alpha/2))
  return(int$y)
}

# Combines the score function and interpolation to output a score CI
# Note: Assumes a 'wide' 99.7% Wald interval to determine range to search for the bounds
fn_score_ci <- function(alpha=0.05, dt, stratum=NULL, t, event, treatment, treat_lvl){
  
  # Estimate the cox coefficient
  est_cox <- fn_cox(dt=dt, t=t, event=event, treatment=treatment)
  
  # Create a 99% Wilcoxon CI for the range to search
  grid_b_low <- seq(est_cox["coef"] - 2.96*est_cox["se(coef)"], est_cox["coef"], length.out = 5)
  grid_b_up <- seq(est_cox["coef"], est_cox["coef"] + 2.96*est_cox["se(coef)"], length.out = 5)
  
  # Interpolate along the grid to find the intersection with relevant quantile
  ci_low <- fn_interpolate(g=grid_b_low, alpha=alpha, dt=dt, stratum=stratum, t=t, event=event, treatment=treatment, treat_lvl=treat_lvl)
  ci_up <- fn_interpolate(g=grid_b_up, alpha=alpha, dt=dt, stratum=stratum, t=t, event=event, treatment=treatment, treat_lvl=treat_lvl)
  ci <- cbind(ci_low, ci_up)
  colnames(ci) <- paste0(round(c(alpha/2,1-alpha/2) * 100,1)," %")
  return(ci)
}

# Example usage
fn_score_ci(alpha=0.05, dt=GBSG2, stratum=NULL, t="time", event="cens", treatment="horTh", treat_lvl="yes")
