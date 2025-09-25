import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform, triang, lognorm


def plot_distribution_from_dict(name, param_distributions, n_samples=10000, bins=50):
    """
    Plot histogram and PDF from a param_distributions dictionary.

    Parameters
    ----------
    name : str
        Name of the parameter (key in param_distributions)
    param_distributions : dict
        Dictionary like { "param": ("dist_name", {"param": value, ...}), ... }
    n_samples : int
        Number of random samples
    bins : int
        Number of bins for histogram
    """

    # --- Look up distribution ---
    if name not in param_distributions:
        raise ValueError(f"{name} not found in dictionary")

    dist_name, params = param_distributions[name]

    # --- Build distribution object ---
    if dist_name == "uniform":
        dist = uniform(**params)
    elif dist_name == "triangular":  # in case you add triangular later
        dist = triang(**params)
    elif dist_name in ["lognorm", "lognormal"]:
        dist = lognorm(**params)
    else:
        raise ValueError(f"Unsupported distribution: {dist_name}")

    # --- Generate samples ---
    samples = dist.rvs(size=n_samples)

    # --- Plot ---
    x = np.linspace(min(samples), max(samples), 500)
    pdf = dist.pdf(x)

    plt.figure(figsize=(8, 5))
    plt.hist(
        samples, bins=bins, density=True, alpha=0.6, color="skyblue", label="Samples"
    )
    plt.plot(x, pdf, "r-", lw=2, label=f"{dist_name} PDF")
    plt.title(f"{name} ~ {dist_name}")
    plt.xlabel(name)
    plt.ylabel("Density")
    plt.legend()
    plt.show()
