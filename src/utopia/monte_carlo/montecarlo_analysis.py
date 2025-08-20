import monaco as mc
import scipy.stats as st
from scipy.stats import randint, rv_discrete, lognorm, uniform
import pandas as pd
from utopia.utopia import utopiaModel
from utopia.results_processing.process_results import *


def run_mc_analysis(base_config, base_data, n_cases=10, param_distributions=None):
    """
    Run Monte Carlo uncertainty & sensitivity analysis on the UTOPIA model.

    Parameters
    ----------
    base_config : dict of config parameters.
    base_data : dict of input data.

    n_cases : int
        Number of Monte Carlo cases.
    param_distributions : dict
        Dictionary: param_name -> (scipy_dist_name, kwargs dict)

    Returns (not implemented yet)
    -------
    results_df : pd.DataFrame
        DataFrame with all inputs and outputs for each case.
    sa_df : pd.DataFrame
        DataFrame with sensitivity analysis metrics from Monaco.
    """

    # Define preprocess for Monaco
    def preprocess(case):
        # Variables that will be passed to the run function and will be sampled
        inputs = copy.deepcopy(base_data)
        # Replace only the sampled ones with numbers
        for param in param_distributions.keys():
            inputs[param] = float(case.invals[param].val)

        return (inputs,)

    # Define run function
    def run_fn(inputs):

        model = utopiaModel(config=base_config, data=inputs)
        model.run()
        processor = ResultsProcessor(model)  # Custom processor to handle results
        # Process results to obtain other outputs such as overall residence time and persistence

        processor.estimate_flows()
        processor.generate_flows_dict()
        processor.process_results()

        # Extract results by compartment and include the concnetration in mass and number for each compartment in the result dictionary (only used one compartment here as example, to be considered if relevant for all or fix a different mechanism to store these results)

        processor.extract_results_by_compartment()
        df = processor.results_by_comp

        # result={"mass_g":processor.Results["mass_g"]}

        processor.estimate_exposure_indicators()
        result = {
            "residence_time_mass": processor.processed_results[
                "Overall_exposure_indicators"
            ]["Overall residence time (years)"][0],
            "residence_time_number": processor.processed_results[
                "Overall_exposure_indicators"
            ]["Overall residence time (years)"][1],
            "persistence_mass": processor.processed_results[
                "Overall_exposure_indicators"
            ]["Overall persistence (years)"][0],
            "persistence_number": processor.processed_results[
                "Overall_exposure_indicators"
            ]["Overall persistence (years)"][1],
            "C_g_m3_Ocean_Surface_Water": df.loc[
                df["Compartments"] == "Ocean_Surface_Water", "Concentration_g_m3"
            ][0],
        }
        return (result,)

    # Define postprocess for Monaco
    def postprocess(case, result):
        for k, v in result.items():
            case.addOutVal(k, v)

    # Create simulation
    sim = mc.Sim(
        name="UTOPIA_MC_simulation",
        ndraws=n_cases,
        fcns={"preprocess": preprocess, "run": run_fn, "postprocess": postprocess},
        firstcaseismedian=False,
        seed=12362398,
        singlethreaded=True,
        savecasedata=False,
        savesimdata=False,
        verbose=True,
        debug=True,
    )

    # Add input distributions
    for param, (dist_name, dist_kwargs) in param_distributions.items():
        sim.addInVar(param, getattr(st, dist_name), dist_kwargs)

    # Run simulation
    sim.runSim()

    # Process results (to be added)

    return sim
