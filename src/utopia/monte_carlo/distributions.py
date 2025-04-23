from scipy.stats import uniform, lognorm, randint, rv_discrete
import numpy as np


# Define your custom distributions (SciPy version)
def MP_Uniform(low, high):
    """Uniform distribution for microplastic params (SciPy version)"""
    return uniform(loc=low, scale=high-low)  # SciPy uniform(loc=min, width)

def MP_Lognormal(mean, std):
    """Lognormal distribution (SciPy parameterization)"""
    sigma = np.sqrt(np.log(1 + (std**2)/(mean**2)))
    mu = np.log(mean) - 0.5 * sigma**2
    return lognorm(s=sigma, scale=np.exp(mu))

def MP_Discrete(values, weights=None):
    """Custom discrete distribution (e.g., for scenario selection)"""
    if weights is None:
        weights = [1/len(values)] * len(values)  # Uniform weights if unspecified
    return rv_discrete(values=(range(len(values)), weights))