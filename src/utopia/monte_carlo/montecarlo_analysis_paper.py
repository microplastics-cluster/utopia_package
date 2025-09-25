import monaco as mc
import scipy.stats as st
from scipy.stats import randint, rv_discrete, lognorm, uniform
import pandas as pd
from utopia.utopia import utopiaModel
from utopia.results_processing.process_results import *
from scipy.stats import spearmanr


def run_mc_analysis_paper(
    base_config,
    base_data,
    art,
    n_cases,
    emission_comp,
    param_distributions=None,
    target_comp=None,
):
    """
    Run Monte Carlo uncertainty & sensitivity analysis on the UTOPIA model and perform comparison with observed data.

    Parameters
    ----------
    base_config : dict of config parameters.
    base_data : dict of input data.
    overlay_compartments_mapping : dict
        Mapping of observed datasets to model compartments.
    observed_df : pd.DataFrame
        DataFrame with observed data, must contain columns: 'Article', 'log_Size', 'log_Abundance'.
    n_cases : int
        Number of Monte Carlo cases.
    param_distributions : dict
        Dictionary: param_name -> (scipy_dist_name, kwargs dict)

    Returns (not implemented yet)
    -------
    results_df : pd.DataFrame
        DataFrame with all inputs and outputs for each case.
    sa_df : pd.DataFrame
        DataFrame with sensitivity analysis metrics from Monaco.?
    """

    # LOAD OBSERVED DATA from Kooi and Koelmans 2019
    # (DOI: https://doi.org/10.1021/acs.estlett.9b00379)

    observed_file = "../src/utopia/data/observed_data_long.xlsx"

    observed_df = pd.read_excel(observed_file)

    observed_df = observed_df.rename(
        columns={
            "dataset identifier": "Dataset",
            "Article name": "Article",
            "log‑transformed size": "log_Size",
            "log‑transformed abundance": "log_Abundance",
        }
    )

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
        processor.estimate_flows()
        processor.generate_flows_dict()
        processor.process_results()

        # Extract results by compartment and include the concentration in mass and number for each compartment in the result dictionary (only used one compartment here as example, to be considered if relevant for all or fix a different mechanism to store these results)

        processor.extract_results_by_compartment()
        df = processor.results_by_comp
        df2 = processor.Results_extended
        # processor.estimate_exposure_indicators()

        # Process model outputs into log-relative abundances per size fraction and compare with observed data

        obs_subset = observed_df[observed_df["Article"] == art]
        x_obs_all = obs_subset["log_Size"].astype(float).values
        y_obs_all = obs_subset["log_Abundance"].astype(float).values

        if len(x_obs_all) > 2:
            slope_obs, intercept_obs = np.polyfit(x_obs_all, y_obs_all, 1)
        else:
            slope_obs, intercept_obs = np.nan, np.nan

        result = {}

        # comp_relative_abundance_resultsaccording to target compartment set in the input parameters

        # Exclude size bins that are not withing the observed data range (0.5 and 5 um)
        excluded_sizes = [0.5, 5]
        filtered_Results = df2[~df2["Size_Fraction_um"].isin(excluded_sizes)]
        rel_abun = extract_log_rel_abundance(
            filtered_Results,
            compartment=target_comp,
            value_col="number_of_particles",
        )

        if rel_abun.empty:
            # fill with NaNs
            result.update(
                {
                    "Spearman_r": np.nan,
                    "Spearman_p": np.nan,
                    "Pass_Spearman": False,
                    "RMSE": np.nan,
                    "R_squared": np.nan,
                    "Model_slope": np.nan,
                    "Model_intercept": np.nan,
                    "Target_compartment": target_comp,
                    "Emission_compartment": emission_comp,
                }
            )

        else:

            x_mod_all = rel_abun["log_size"].values
            y_mod_all = rel_abun["log_rel_abundance"].values

            # Spearman correlation
            r_value, p_value = spearmanr(x_mod_all, y_mod_all)

            spearman_threshold = -1

            # Initialize metrics
            rmse, r2, slope_mod, intercept_mod = np.nan, np.nan, np.nan, np.nan
            if np.isnan(r_value):
                pass_spearman = False
            elif (r_value) == spearman_threshold:
                pass_spearman = True
            else:
                pass_spearman = False

            # Compute RMSE, R², slope, intercept only if threshold passed

            if pass_spearman == True and len(x_mod_all) > 2 and len(x_obs_all) > 2:
                slope_mod, intercept_mod = np.polyfit(x_mod_all, y_mod_all, 1)

                # Align with observed range
                x_min = max(x_obs_all.min(), x_mod_all.min())
                x_max = min(x_obs_all.max(), x_mod_all.max())
                mask_obs = (x_obs_all >= x_min) & (x_obs_all <= x_max)
                x_obs_filt = x_obs_all[mask_obs]
                y_obs_filt = y_obs_all[mask_obs]

                if len(x_obs_filt) > 2:
                    pred_y = slope_mod * x_obs_filt + intercept_mod
                    ss_res = np.sum((y_obs_filt - pred_y) ** 2)
                    ss_tot = np.sum((y_obs_filt - np.mean(y_obs_filt)) ** 2)
                    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
                    rmse = np.sqrt(ss_res / len(x_obs_filt))

                    # Plot observed vs model with metrics only if spearman test is true
                    # plot_obs_vs_model_with_metrics_True(
                    #     x_obs_all,
                    #     y_obs_all,
                    #     x_mod_all,
                    #     y_mod_all,
                    #     comp,
                    #     art,
                    #     rmse,
                    #     r2,
                    #     MP_density=model.MPdensity_kg_m3,
                    #     FI=model.FI,
                    #     t_half_deg=model.t_half_deg_free,
                    #     t_frag=model.t_frag_gen_FreeSurfaceWater,
                    # )

                    # Store all metrics in result dict
                    result.update(
                        {
                            "Spearman_r": r_value,
                            "Spearman_p": p_value,
                            "Pass_Spearman": pass_spearman,
                            "RMSE": rmse,
                            "R_squared": r2,
                            "Model_slope": slope_mod,
                            "Model_intercept": intercept_mod,
                            "Target_compartment": target_comp,
                            "Emission_compartment": emission_comp,
                        }
                    )
                else:
                    # If not enough points to compute metrics, store NaNs
                    result.update(
                        {
                            "Spearman_r": np.nan,
                            "Spearman_p": np.nan,
                            "Pass_Spearman": False,
                            "RMSE": np.nan,
                            "R_squared": np.nan,
                            "Model_slope": np.nan,
                            "Model_intercept": np.nan,
                            "Target_compartment": target_comp,
                            "Emission_compartment": emission_comp,
                        }
                    )
            else:
                # Store all metrics in result dict
                result.update(
                    {
                        "Spearman_r": np.nan,
                        "Spearman_p": np.nan,
                        "Pass_Spearman": False,
                        "RMSE": np.nan,
                        "R_squared": np.nan,
                        "Model_slope": np.nan,
                        "Model_intercept": np.nan,
                        "Target_compartment": target_comp,
                        "Emission_compartment": emission_comp,
                    }
                )

            # Add global exposure indicators
            # result.update(
            #     {
            #         "residence_time_mass": processor.processed_results[
            #             "Overall_exposure_indicators"
            #         ]["Overall residence time (years)"][0],
            #         "residence_time_number": processor.processed_results[
            #             "Overall_exposure_indicators"
            #         ]["Overall residence time (years)"][1],
            #         "persistence_mass": processor.processed_results[
            #             "Overall_exposure_indicators"
            #         ]["Overall persistence (years)"][0],
            #         "persistence_number": processor.processed_results[
            #             "Overall_exposure_indicators"
            #         ]["Overall persistence (years)"][1],
            #     }
            # )

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
        samplemethod="sobol_random",
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


def extract_log_rel_abundance(model_df, compartment, value_col):
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
    ) * 100  # in percentage

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


def sim_to_dataframe(sim):
    """Convert Monaco simulation results to a pandas DataFrame."""
    records = []
    for case in sim.cases:
        row = {}
        # Inputs
        for name, val in case.invals.items():
            row[name] = val.val
        # Outputs
        for name, val in case.outvals.items():
            row[name] = val.val
        records.append(row)
    return pd.DataFrame(records)


def reshape_mc_results(df, input_params, dataset_name):
    """
    Reshape MC results dataframe into long format for selected compartments.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame from sim_to_dataframe(sim).
    compartments : list of str
        List of target compartments, e.g. ["Sediment_Coast", "Sediment_Ocean"].
    input_params : list of str
        Names of input parameters to keep (from Monte Carlo inputs).
    dataset_name : str
        Observed dataset name (article).
    emission_comp : str
        Emission compartment name (scenario).

    Returns
    -------
    pd.DataFrame
        Long-format dataframe with metrics per target compartment.
    """
    rows = []
    for idx, row in df.iterrows():

        rows.append(
            {
                "Case": idx,
                "Observed_dataset": dataset_name,
                "Emission_Compartment": row.get("Emission_compartment"),
                "Target_Compartment": row.get("Target_compartment"),
                "RMSE": row.get("RMSE"),
                "R2": row.get("R_squared"),
                "Slope": row.get("Model_slope"),
                "Spearman_r": row.get("Spearman_r"),
                "Pass_Spearman": row.get("Pass_Spearman"),
                **{param: row.get(param) for param in input_params},
            }
        )
    return pd.DataFrame(rows)


def plot_obs_vs_model_with_metrics(obs_df, model_df, compartment, art, input_dict):
    """Plot observed vs modeled log-abundance by size with regression and metrics.

    Parameters
    ----------
    obs_df : pd.DataFrame
        Observed dataset with columns x_col and y_obs_col
    model_df : pd.DataFrame
        Modeled dataset with columns x_col and y_mod_col
    compartment : str
        Name of the compartment for title
    art : str
        Name of the article for title"""
    # Observed dataset
    x_obs_all = obs_df["log_Size"].astype(float).values
    y_obs_all = obs_df["log_Abundance"].astype(float).values

    if len(x_obs_all) > 2:
        slope_obs, intercept_obs = np.polyfit(x_obs_all, y_obs_all, 1)
    else:
        slope_obs, intercept_obs = np.nan, np.nan

    # Plot observed  dataset
    x_fit_obs = np.linspace(x_obs_all.min(), x_obs_all.max(), 100)
    y_fit_obs = slope_obs * x_fit_obs + intercept_obs
    plt.figure(figsize=(8, 6))
    plt.scatter(
        x_obs_all, y_obs_all, color="blue", marker="o", label=f"Observed: {art}"
    )

    if x_fit_obs is not None:
        plt.plot(x_fit_obs, y_fit_obs, color="blue", linestyle="--", label="_nolegend_")

    # Plot modelled dataset
    x_mod_all = model_df["log_size"].values
    y_mod_all = model_df["log_rel_abundance"].values

    slope_mod, intercept_mod = np.polyfit(x_mod_all, y_mod_all, 1)
    x_fit_mod = np.linspace(x_mod_all.min(), x_mod_all.max(), 100)
    y_fit_mod = slope_mod * x_fit_mod + intercept_mod
    plt.scatter(x_mod_all, y_mod_all, marker="x", label=f"Model data: {compartment}")
    plt.plot(x_fit_mod, y_fit_mod, label="_nolegend_")

    # Align modelled with observed range to copute metrics
    x_min = max(x_obs_all.min(), x_mod_all.min())
    x_max = min(x_obs_all.max(), x_mod_all.max())
    mask_obs = (x_obs_all >= x_min) & (x_obs_all <= x_max)
    x_obs_filt = x_obs_all[mask_obs]
    y_obs_filt = y_obs_all[mask_obs]

    # Compute metrics
    pred_y = slope_mod * x_obs_filt + intercept_mod
    ss_res = np.sum((y_obs_filt - pred_y) ** 2)
    ss_tot = np.sum((y_obs_filt - np.mean(y_obs_filt)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = np.sqrt(ss_res / len(x_obs_filt))

    # Plotting with metrics

    plt.xlabel("log(Size [µm])")
    plt.ylabel("log(Relative Abundance [%])")
    plt.title(
        f"Overlay: {art} - {compartment}-{input_dict['MPdensity_kg_m3']}-{input_dict['FI']}-{input_dict['t_half_deg_free']}-{input_dict['t_frag_gen_FreeSurfaceWater']}"
    )
    # Annotate metrics
    plt.text(
        0.05,
        0.95,
        f"RMSE = {rmse:.2f}\nR² = {r2:.2f}",
        transform=plt.gca().transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
    )
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_obs_vs_model_with_metrics_True(
    x_obs_all,
    y_obs_all,
    x_mod_all,
    y_mod_all,
    comp,
    art,
    rmse,
    r2,
    MP_density,
    FI,
    t_half_deg,
    t_frag,
):
    """Plot observed vs modeled log-abundance by size with regression and metrics when spearman test is true.

    Parameters
    ----------
    obs_df : pd.DataFrame
        Observed dataset with columns x_col and y_obs_col
    model_df : pd.DataFrame
        Modeled dataset with columns x_col and y_mod_col
    compartment : str
        Name of the compartment for title
    art : str
        Name of the article for title"""

    slope_obs, intercept_obs = np.polyfit(x_obs_all, y_obs_all, 1)
    x_fit_obs = np.linspace(x_obs_all.min(), x_obs_all.max(), 100)
    y_fit_obs = slope_obs * x_fit_obs + intercept_obs
    plt.figure(figsize=(8, 6))
    plt.scatter(
        x_obs_all, y_obs_all, color="blue", marker="o", label=f"Observed: {art}"
    )

    if x_fit_obs is not None:
        plt.plot(x_fit_obs, y_fit_obs, color="blue", linestyle="--", label="_nolegend_")

    # Plot modelled dataset (filtered for excluded sizes)

    slope_mod, intercept_mod = np.polyfit(x_mod_all, y_mod_all, 1)
    x_fit_mod = np.linspace(x_mod_all.min(), x_mod_all.max(), 100)
    y_fit_mod = slope_mod * x_fit_mod + intercept_mod
    plt.scatter(x_mod_all, y_mod_all, marker="x", label=f"Model data: {comp}")
    plt.plot(x_fit_mod, y_fit_mod, label="_nolegend_")

    # Compute metrics

    plt.xlabel("log(Size [µm])")
    plt.ylabel("log(Relative Abundance [%])")
    plt.title(f"Overlay: {art} - {comp}-{MP_density}-{FI}-{t_half_deg}-{t_frag}")
    # Annotate metrics
    plt.text(
        0.05,
        0.95,
        f"RMSE = {rmse:.2f}\nR² = {r2:.2f}",
        transform=plt.gca().transAxes,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
    )
    plt.legend()
    plt.grid(True)
    plt.show()


from utopia.utopia import utopiaModel
from utopia.results_processing.process_results import *
from utopia.monte_carlo.montecarlo_analysis_paper import extract_log_rel_abundance


def run_and_plot_top_results(top_results, data_data, config_data):
    observed_file = "../src/utopia/data/observed_data_long.xlsx"

    observed_df = pd.read_excel(observed_file)

    observed_df = observed_df.rename(
        columns={
            "dataset identifier": "Dataset",
            "Article name": "Article",
            "log‑transformed size": "log_Size",
            "log‑transformed abundance": "log_Abundance",
        }
    )
    for i in range(len(top_results)):
        art = top_results.iloc[i]["Observed_dataset"]
        comp_E = top_results.iloc[i]["Emission_Compartment"]
        comp_T = top_results.iloc[i]["Target_Compartment"]
        new_data = data_data.copy()
        new_data["emiss_dict_g_s"] = set_emission(
            new_data["emiss_dict_g_s"].copy(), comp_E, "e", 100
        )
        new_data["MPdensity_kg_m3"] = top_results.iloc[i]["MPdensity_kg_m3"]
        new_data["FI"] = top_results.iloc[i]["FI"]
        new_data["t_half_deg_free"] = top_results.iloc[i]["t_half_deg_free"]
        new_data["t_frag_gen_FreeSurfaceWater"] = top_results.iloc[i][
            "t_frag_gen_FreeSurfaceWater"
        ]
        model_top = utopiaModel(config=config_data, data=new_data)
        model_top.run()
        processor_top = ResultsProcessor(model_top)
        processor_top.estimate_flows()
        processor_top.generate_flows_dict()
        processor_top.process_results()
        processor_top.extract_results_by_compartment()
        df = processor_top.results_by_comp
        df2 = processor_top.Results_extended
        processor_top.estimate_exposure_indicators()
        model_df = df2.copy()
        excluded_sizes = [0.5, 5]
        filtered_Results = model_df[~model_df["Size_Fraction_um"].isin(excluded_sizes)]
        rel_abun = extract_log_rel_abundance(
            filtered_Results,
            compartment=comp_T,
            value_col="number_of_particles",
        )

        input_dict = {
            "MPdensity_kg_m3": new_data["MPdensity_kg_m3"],
            "FI": new_data["FI"],
            "t_half_deg_free": new_data["t_half_deg_free"],
            "t_frag_gen_FreeSurfaceWater": new_data["t_frag_gen_FreeSurfaceWater"],
        }

        obs_subset = observed_df[observed_df["Article"] == art]

        plot_obs_vs_model_with_metrics(
            obs_df=obs_subset,
            model_df=rel_abun,
            compartment=comp_T,
            art=art,
            input_dict=input_dict,
        )
