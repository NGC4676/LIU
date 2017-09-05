#==============================================================================
#  A simple fit with PyMC can be accomplished as follows. Here we will fit 
# the mean of a distribution—perhaps an overly simplistic example for MCMC,
# but useful as an introductory example:
#==============================================================================
import numpy as np
import pymc
from matplotlib import pyplot as plt
# generate random Gaussian data with mu=0, sigma=1
N = 10000
x = np.random.normal(size=N)

# define the MCMC model: uniform prior on mu , fixed (known) sigma
mu = pymc.Uniform('mu', -5, 5)               #输入先验几率
sigma = pymc.Uniform('sigma', -5, 5)               #输入先验几率
M = pymc.Normal('M', mu , sigma ,
					  observed=True , value=x)     #value处输入数据
model = dict(M=M, mu=mu, sigma=sigma)                #这一步把mu和M两个模型映射到字符串？

# run the model , and get the trace of mu
S = pymc.MCMC(model)
S.sample(50000, burn=10000)
sigma_sample = S.trace('sigma')[:]
mu_sample = S.trace('mu')[:]

ax1=plt.subplot(231)
ax2=plt.subplot(232)
ax1.hist(x,bins=100, histtype='stepfilled')
ax2.hist(sigma_sample,bins=100, histtype='stepfilled')

# print the MCMC estimate
print ''
print("Bayesian (MCMC): %.3f +/- %.3f" 
% (np.mean(sigma_sample), np.std(sigma_sample)))

# compare to the frequentist estimate
print("Frequentist: %.3f +/- %.3f"
 % (np.mean(x), np.std(x, ddof=1) / np.sqrt(N)))
