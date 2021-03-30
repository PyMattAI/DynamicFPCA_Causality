# DynamicFPCA_Causality
# Function DFPCA.test(X, Y, num, ndx, ndy, L, q)
#
# Non-causality test in the sense of Granger between functional observations X and Y if X cause Y.
#
# Parameters:
# num: number of observations by curves 
# X: data matrix num*n with n sample size
# Y: data matrix num*n with n sample size
# ndx: number of Dynamic functional principal components used for X
# ndy: number of Dynamic functional principal components used for Y
# L: non negative integer, DFPCA filters coefficients at lags L. By default L=30.
# q: positive integer, window size for the kernel estimator of the spectral density estimator.By default, q=max(1,floor(n^0.33)). 
