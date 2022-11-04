rm(list = ls())
# for lalonde data
library(LalRUtils)
# libraries for fast computation
libreq(Matrix, matrixStats, data.table, fixest, CVXR)

source("R/ebweights.R")
source("R/qrweights.R")
# %% ############################################################
data(lalonde.psid); setDT(lalonde.psid); df = lalonde.psid
yn = 're78'; wn = 'treat'; Xn = setdiff(colnames(df), c(yn, wn))
# rest of the code is generic
X = df[get(wn) == 0, ..Xn] |> as.matrix()
# treatment covariate means (for ATT)
target = colMeans2(as.matrix(df[get(wn) == 1, ..Xn]))
# %%
wts_eb = eb_solve_primal(c(1, target), cbind(1, X))
wts_qr = qr_solve_primal(c(1, target), cbind(1, X))

# %%
mean(df[get(wn)==1, get(yn)]) - weighted.mean(df[get(wn) == 0, get(yn)], wts_eb)
mean(df[get(wn)==1, get(yn)]) - weighted.mean(df[get(wn) == 0, get(yn)], wts_qr)

# %%
