rm(list = ls())
# for lalonde data
library(LalRUtils)
# libraries for fast computation
libreq(Matrix, matrixStats, data.table, fixest, CVXR)
source("R/ebweights.R")

# %%
############################################################
# ATT
############################################################

data(lalonde.psid); setDT(lalonde.psid); df = lalonde.psid
yn = 're78'; wn = 'treat'; Xn = setdiff(colnames(df), c(yn, wn))

# rest of the code is generic
X = df[get(wn) == 0, ..Xn] |> as.matrix()
# treatment covariate means (for ATT)
target = colMeans2(as.matrix(df[get(wn) == 1, ..Xn]))

# solve for balancing weights
ω = entrBal(X, target)

# ATT estimate
regdf = rbind(
  # treated obs : uniform weights
  data.table(y = df[get(wn) == 1, get(yn)], w = 1, wt = 1 / sum(df[[wn]])),
  # untreated obs : solved weights
  data.table(y = df[get(wn) == 0, get(yn)], w = 0, wt = ω)
)

feols(y ~ w, weights = ~wt, data = regdf, vcov = "HC1")

# OLS estimation, Dep. Var.: y
# Observations: 2,675
# Standard-errors: Heteroskedasticity-robust
#             Estimate Std. Error t value        Pr(>|t|)
# (Intercept)     3924      659.8   5.948 0.0000000030644 ***
# w               2425      876.5   2.766 0.0057107947524 **
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# RMSE: 205.2   Adj. R2: 0.025064

# %% check that primal and dual solutions are identical
ω_primal = eb_solve_primal(c(1, target), cbind(1, X))
ωs = cbind(ω, ω_primal)
# should lie along 45° line
plot(ωs[, 1], ωs[,2]); abline(0, 1)

# %%
library(microbenchmark)
microbenchmark(
  primal = {
    eb_solve_primal(c(1, target), cbind(1, X))
  },
  dual = {
    eb_solve_dual(  c(1, target), cbind(1, X))
  },
  times = 100L
)

# Unit: milliseconds
#    expr    min    lq   mean median     uq    max neval cld
#  primal 570.58 634.1 740.75 731.61 774.22 1826.5   100   b
#    dual  15.33  33.6  60.46  53.03  73.82  322.3   100  a

# %%
############################################################
# ATE
############################################################

data(lalonde.exp); setDT(lalonde.exp); df = lalonde.exp
yn = 're78'; wn = 'treat'; Xn = setdiff(colnames(df), c(yn, wn))

X0 = df[get(wn) == 0, ..Xn] |> as.matrix()
X1 = df[get(wn) == 1, ..Xn] |> as.matrix()

# overall covariate means (for ATE)
target2 = colMeans2(as.matrix(df[, ..Xn]))

# solve for balancing weights - treatment obs
ω1 = entrBal(X1, target2)
ω0 = entrBal(X0, target2)

regdf = rbind(
	# treated obs : uniform weights of 1
	data.table(y = df[get(wn)== 1, get(yn)], w = 1, wt = ω1),
	data.table(y = df[get(wn)== 0, get(yn)], w = 0, wt = ω0)
)

feols(y ~ w, weights = ~wt, data = regdf, vcov = "HC1")

# OLS estimation, Dep. Var.: y
# Observations: 445
# Standard-errors: Heteroskedasticity-robust
#             Estimate Std. Error t value  Pr(>|t|)
# (Intercept)     4591      353.3   12.99 < 2.2e-16 ***
# w               1572      689.3    2.28  0.023057 *
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# RMSE: 453.8   Adj. R2: 0.011074


