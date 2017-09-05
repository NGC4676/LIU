#==============================================================================
# https://pymc-devs.github.io/pymc/tutorial.html
#==============================================================================
#Pymc的几个重要属性——parents，children，value， logp

#所谓的probability model就是指linked collection of variables
import matplotlib.pyplot as plt
from pymc import DiscreteUniform, Exponential, deterministic, Poisson
import numpy as np
import pymc
#==============================================================================
# Building models
#==============================================================================
#input data
disasters_array =  np.array(
						[ 4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6,
                   3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
                   2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0,
                   1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
                   0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2,
                   3, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
                   0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1])

#离散均匀分布先验密度
switchpoint = DiscreteUniform('switchpoint', lower=0, upper=110, doc='Switchpoint[year]')

#指数分布先验密度
early_mean = Exponential('early_mean', beta=1.)
late_mean = Exponential('late_mean', beta=1.)

#定义随机变量rate，在转折点改变
@deterministic(plot=False)  #这一步把rate转化成deterministic类
def rate(s=switchpoint, e=early_mean, l=late_mean):
    ''' Concatenate Poisson means '''
    out = np.empty(len(disasters_array))
    out[:s] = e
    out[s:] = l
    return out

#定义随机变量disaster表示事故发生次数，observed=True表示值不变（默认False）
disasters = Poisson('disasters', mu=rate, value=disasters_array, observed=True)

#==============================================================================
# Fit with MCMC
#==============================================================================
model=dict(disasters=disasters,switchpoint=switchpoint,
											early_mean=early_mean,
											late_mean=late_mean)
#M will expose variables switchpoint, early_mean, late_mean and disasters as attributes
M=pymc.MCMC(model)  

#characterizing its posterior distribution (Metropolis-Hastings algorithm)
M.sample(iter=11000, burn=1000, thin=10)  #thin——retaining every k th sample

#Histogram of the marginal posterior probability
plt.hist(M.trace('late_mean')[:])

#For each variable in the model, plot generates a composite figure
pymc.Matplot.plot(M)