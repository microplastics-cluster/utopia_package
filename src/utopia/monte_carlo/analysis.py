import monaco as mc
from typing import Dict, Any
import numpy as np

def run_monte_carlo(model_func, param_distributions, n_simulations=1000, **kwargs):
    """
    Core Monte Carlo runner compatible with utopiaModel outputs.
    """
    sim = mc.Simulation(
        name='UTOPIA_Sensitivity',
        model=model_func,
        first_parameter=param_distributions,
        n_samples=n_simulations,
        **kwargs
    )
    sim.run()
    sim.compute_sensitivities()
    return sim

def summarize_sensitivity(sim, output_key='outputs'):
    """
    Extract key sensitivity metrics in a utopia-friendly format.
    """
    return {
        'first_order': sim.sensitivities.first_order[output_key],
        'total_order': sim.sensitivities.total_order[output_key],
        'output_distributions': {
            k: {'mean': np.mean(v), 'std': np.std(v)}
            for k, v in sim.outvals[output_key].items()
        }
    }