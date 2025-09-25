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
        df2 = processor.Results_extended

        # Process model outputs into log-relative abundances per size fraction.
        comp_relative_abundance_number = {}
        for comp in df["Compartments"].unique():
            comp_relative_abundance_number[comp] = extract_log_rel_abundance(
                processor.Results_extended,
                compartment=comp,
                value_col="number_fraction",
            )

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
            "log_relative_abundance_num_5um_Ocean_Surface_Water": comp_relative_abundance_number[
                "Ocean_Surface_Water"
            ]
            .loc[
                comp_relative_abundance_number["Ocean_Surface_Water"][
                    "Size_Fraction_um"
                ]
                == 5,
                "log_rel_abundance",
            ]
            .values[0],
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


def extract_log_rel_abundance(model_df, compartment, value_col="fraction"):
    """
    Process model outputs into log-relative abundances per size fraction.

    Parameters
    ----------
    model_df : pd.DataFrame
        Must contain columns: 'compartment', 'aggregation_state', 'size_fraction', and a value column (e.g. 'fraction').
    compartment : str
        The compartment to filter on (e.g., 'water').
    value_col : str
        The column with the fraction values (e.g., 'fraction_mass' or 'fraction_number').

    Returns
    -------
    pd.DataFrame
        Columns: ['size_fraction', 'relative_abundance', 'log_size', 'log_rel_abundance']
    """

    # 1. Filter for compartment of interest
    df = model_df[model_df["Compartment"] == compartment].copy()

    # 2. Aggregate over aggregation states (sum fractions within size fractions)
    df_grouped = df.groupby("Size_Fraction_um")[value_col].sum().reset_index()

    # 3. Normalize to relative abundance across all size fractions
    df_grouped["relative_abundance"] = (
        df_grouped[value_col] / df_grouped[value_col].sum()
    )

    # 4. Convert to log10 scale
    df_grouped["log_size"] = np.log10(df_grouped["Size_Fraction_um"])
    df_grouped["log_rel_abundance"] = np.log10(df_grouped["relative_abundance"])

    return df_grouped[
        ["Size_Fraction_um", "relative_abundance", "log_size", "log_rel_abundance"]
    ]


def build_result_dict(
    processor, df, comp_relative_abundance_number, compartments, sizes
):
    result = {
        "residence_time_mass": processor.processed_results[
            "Overall_exposure_indicators"
        ]["Overall residence time (years)"][0],
        "residence_time_number": processor.processed_results[
            "Overall_exposure_indicators"
        ]["Overall residence time (years)"][1],
        "persistence_mass": processor.processed_results["Overall_exposure_indicators"][
            "Overall persistence (years)"
        ][0],
        "persistence_number": processor.processed_results[
            "Overall_exposure_indicators"
        ]["Overall persistence (years)"][1],
    }

    # Add concentrations per compartment
    for comp in compartments:
        key = f"C_g_m3_{comp}"
        value = df.loc[df["Compartments"] == comp, "Concentration_g_m3"].values[0]
        result[key] = value

    # Add relative abundances per (size, compartment)
    for comp in compartments:
        for size in sizes:
            key = f"log_relative_abundance_num_{size}um_{comp}"
            value = (
                comp_relative_abundance_number[comp]
                .loc[
                    comp_relative_abundance_number[comp]["Size_Fraction_um"] == size,
                    "log_rel_abundance",
                ]
                .values[0]
            )
            result[key] = value

    return result


def set_emission(emissions_dict, compartment, fraction, value):
    """
    Set emission value in one compartment/fraction, zero elsewhere.

    Parameters
    ----------
    emissions_dict : dict
        Nested dictionary {compartment: {fraction: value}}
    compartment : str
        The compartment to receive the emission (key of outer dict)
    fraction : str
        The fraction (key of inner dict) to receive the emission
    value : float or int
        The emission value to set

    Returns
    -------
    dict
        New dictionary with emission localized
    """
    # Make a deep copy so original dict isn’t modified
    from copy import deepcopy

    new_dict = deepcopy(emissions_dict)

    # Zero everything
    for comp in new_dict:
        for frac in new_dict[comp]:
            new_dict[comp][frac] = 0

    # Assign chosen emission
    if compartment not in new_dict:
        raise ValueError(f"Compartment '{compartment}' not found")
    if fraction not in new_dict[compartment]:
        raise ValueError(f"Fraction '{fraction}' not found in {compartment}")

    new_dict[compartment][fraction] = value
    return new_dict
