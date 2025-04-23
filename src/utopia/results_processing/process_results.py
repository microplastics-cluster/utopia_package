import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LogNorm
import pandas as pd
from utopia.helpers import *
from utopia.preprocessing.fill_interactions_dictionaries import *
from utopia.results_processing.exposure_indicators_calculation import *
from utopia.solver_steady_state import *
from utopia.results_processing.emission_fractions_calculation import *
#from utopia.results_processing.pdf_reporting import *


class ResultsProcessor:
    """Provides functionalities for restructuring, analysing and plotting the UTOPIA model results."""

    def __init__(self, model):
        self.processed_results = {}  # empty dictionary to store results
        self.model = model
        self.R = model.R
        self.Results_extended = None

    def process_results(self):
        """Reformat results dataframe for easier analysis by specifying size fractions, MP forms and compartments and deriving mass and number fractions, input and outup flows."""
        # Reformat results (R) dataframe
        median_size_of_bin=[self.model.size_dict[x[0]] for x in self.R.index]
        self.R["Size_Fraction_um"] = [self.model.size_distribution_dict[str(x)] for x in median_size_of_bin]
        self.R["MP_Form"] = [
            self.model.MP_form_dict_reverse[x[1]] for x in self.R.index
        ]
        self.R["Compartment"] = [
            self.model.comp_dict_inverse[float(x[2:-7])] for x in self.R.index
        ]

        Results = self.R[
            [
                "Compartment",
                "MP_Form",
                "Size_Fraction_um",
                "mass_g",
                "number_of_particles",
                "concentration_g_m3",
                "concentration_num_m3",
            ]
        ]
        # Calculate mass and number fractions relative to the total mass and number of particles and store in new dataframe "Results_extended" (needed for plotting results heatmaps)
        total_mass = sum(Results["mass_g"])
        total_number = sum(Results["number_of_particles"])
        Results_extended = Results.copy()
        Results_extended.loc[:, "mass_fraction"] = [
            x / total_mass for x in Results["mass_g"]
        ]
        Results_extended.loc[:, "number_fraction"] = [
            x / total_number for x in Results["number_of_particles"]
        ]

        mass_fraction_df = Results_extended.loc[
            :, ["Compartment", "MP_Form", "Size_Fraction_um", "mass_fraction"]
        ]

        number_fraction_df = Results_extended.loc[
            :, ["Compartment", "MP_Form", "Size_Fraction_um", "number_fraction"]
        ]

        """ Add input and output flows dict to results extended dataframe"""

        Results_extended2 = self.addFlows_to_results_df(Results_extended)

        """ Fix input flows dict to results extended dataframe"""

        interactions_pp_df = fillInteractions_fun_OOP_dict(
            self.model.system_particle_object_list,
            self.model.SpeciesList,
            self.surfComp_list,
        )
        # Estimate Pnum_SS (particle number at steady state) for each particle object in the system
        for p in self.model.system_particle_object_list:
            p.Pnum_SS = mass_to_num(p.Pmass_g_SS, p.Pvolume_m3, p.Pdensity_kg_m3)
        # Create a dictionary of recieving inflows per particle taking the values from the interactions matrix
        particle_inflows_dict_mass = {}
        particle_inflows_dict_number = {}
        for p in self.model.system_particle_object_list:
            inflows_p_mass = []
            inflows_p_num = []
            emission_rate_g_s = self.model.emiss_dict_g_s[p.Pcompartment.Cname][
                p.Pcode[0]
            ]
            emission_rate_num_s = mass_to_num(
                emission_rate_g_s, p.Pvolume_m3, p.Pdensity_kg_m3
            )
            for p2 in self.model.system_particle_object_list:
                interaction_rate = interactions_pp_df[p2.Pcode][p.Pcode]
                if type(interaction_rate) == dict:
                    inflow = {k: v * p2.Pmass_g_SS for k, v in interaction_rate.items()}
                    inflows_p_mass.append(inflow)
                    inflows_p_num.append(
                        {k: v * p2.Pnum_SS for k, v in interaction_rate.items()}
                    )
                else:
                    inflows_p_mass.append(interaction_rate * p2.Pmass_g_SS)
                    inflows_p_num.append(interaction_rate * p2.Pnum_SS)
            dict_list = [item for item in inflows_p_mass if isinstance(item, dict)]
            dict_list_num = [item for item in inflows_p_num if isinstance(item, dict)]
            merged_dict = {}
            merged_dict_num = {}
            for d in dict_list:
                for k, v in d.items():
                    if k in merged_dict:
                        merged_dict[k] += v
                    else:
                        merged_dict[k] = v
            for d in dict_list_num:
                for k, v in d.items():
                    if k in merged_dict_num:
                        merged_dict_num[k] += v
                    else:
                        merged_dict_num[k] = v

            particle_inflows_dict_mass[p.Pcode] = merged_dict
            particle_inflows_dict_number[p.Pcode] = merged_dict_num
            # Add the emission rate to the inflow dictionary
            merged_dict["Emission_flow"] = emission_rate_g_s
            merged_dict_num["Emission_flow"] = emission_rate_num_s

        # Substitute the inputflow values in the results_extended dataframe:
        for ele in particle_inflows_dict_mass:
            Results_extended2.at[ele, "inflows_g_s"] = particle_inflows_dict_mass[ele]
        for ele in particle_inflows_dict_number:
            Results_extended.at[ele, "inflows_num_s"] = particle_inflows_dict_number[
                ele
            ]
        # Add total input and putput flows to Results extended dataframe
        Results_extended2["Total_inflows_g_s"] = [
            sum(Results_extended2.iloc[i].inflows_g_s.values())
            for i in range(len(Results_extended2))
        ]

        Results_extended2["Total_outflows_g_s"] = [
            sum(Results_extended2.iloc[i].outflows_g_s.values())
            for i in range(len(Results_extended2))
        ]

        Results_extended2["Total_inflows_num_s"] = [
            sum(Results_extended2.iloc[i].inflows_num_s.values())
            for i in range(len(Results_extended2))
        ]

        Results_extended2["Total_outflows_num_s"] = [
            sum(Results_extended2.iloc[i].outflows_num_s.values())
            for i in range(len(Results_extended2))
        ]
        self.Results_extended = Results_extended2
        self.processed_results["Results_extended"] = Results_extended2

    def addFlows_to_results_df(self, Results_extended):
        """Calculate inflows and outflows (mass and number) and update Results_extended."""
        inflows_mass_list = []
        inflows_num_list = []
        outflows_mass_list = []
        outflows_num_list = []

        for n in range(len(Results_extended)):
            compartment = Results_extended.iloc[n]["Compartment"]
            size_fraction = Results_extended.iloc[n]["Size_Fraction_um"]
            mp_form = Results_extended.iloc[n]["MP_Form"]

            # Calculate inflows and outflows for mass
            inflows_mass = process_flows(
                compartment, size_fraction, mp_form, "input_flows", self.flows_dict_mass
            )
            outflows_mass = process_flows(
                compartment,
                size_fraction,
                mp_form,
                "output_flows",
                self.flows_dict_mass,
            )
            inflows_mass_list.append(inflows_mass)
            outflows_mass_list.append(outflows_mass)

            # Calculate inflows and outflows for number
            inflows_num = process_flows(
                compartment,
                size_fraction,
                mp_form,
                "input_flows",
                self.flows_dict_number,
            )
            outflows_num = process_flows(
                compartment,
                size_fraction,
                mp_form,
                "output_flows",
                self.flows_dict_number,
            )
            inflows_num_list.append(inflows_num)
            outflows_num_list.append(outflows_num)

        # Update the Results_extended DataFrame with the calculated flows
        Results_extended["inflows_g_s"] = inflows_mass_list
        Results_extended["inflows_num_s"] = inflows_num_list
        Results_extended["outflows_g_s"] = outflows_mass_list
        Results_extended["outflows_num_s"] = outflows_num_list

        return Results_extended

    def plot_fractionDistribution_heatmaps(self, fraction):
        """Plots the mass and number fractions after they have been extracted to the Results_extended df."""
        if self.Results_extended is None:
            raise ValueError(
                "Mass and particle number fractions not extracted. Call process_results() first."
            )

        # Define the order for the MP_Form labels
        mp_form_order = [
            "freeMP",
            "heterMP",
            "biofMP",
            "heterBiofMP",
        ]  # Replace with your desired order

        # Define the order for the Compartment labels
        compartment_order = [
            "Ocean_Surface_Water",
            "Ocean_Mixed_Water",
            "Ocean_Column_Water",
            "Coast_Surface_Water",
            "Coast_Column_Water",
            "Surface_Freshwater",
            "Bulk_Freshwater",
            "Sediment_Freshwater",
            "Sediment_Ocean",
            "Sediment_Coast",
            "Beaches_Soil_Surface",
            "Beaches_Deep_Soil",
            "Background_Soil_Surface",
            "Background_Soil",
            "Impacted_Soil_Surface",
            "Impacted_Soil",
            "Air",
        ]  # Replace with your desired order

        # Pivot the DataFrame to have one row per combination of MP_Form, Compartment, and Size_Fraction_um
        pivot_table = self.Results_extended.pivot_table(
            index=["MP_Form", "Size_Fraction_um"],
            columns="Compartment",
            values=fraction,
            aggfunc="mean",
        )

        # Reorder the rows based on mp_form_order and columns based on compartment_order
        pivot_table = pivot_table.loc[mp_form_order, compartment_order]

        # Apply log scale to the pivot table
        pivot_table_log = np.log10(pivot_table)

        # Replace -inf values with NaN
        pivot_table_log.replace(-np.inf, np.nan, inplace=True)

        # Stablish a lower limit
        # Set the lower limit for the values
        lower_limit = -14
        upper_limit = np.nanmax(pivot_table_log)

        # Replace values below the lower limit with NaN
        pivot_table_log = pivot_table_log.applymap(
            lambda x: np.nan if x < lower_limit else x
        )

        # Define a custom colormap with grey color for NaN values
        cmap = sns.color_palette("viridis", as_cmap=True)
        cmap.set_bad("white")

        # Plot the heatmap with logarithmic scale and custom colormap
        plt.figure(figsize=(12, 8))
        sns.heatmap(
            pivot_table_log,
            cmap=cmap,
            cbar=True,
            cbar_kws={"label": "log10 (" + fraction + ") "},
            annot=False,
            linewidths=0.5,
            linecolor="grey",
            vmin=lower_limit,
            vmax=upper_limit,
        )

        # Set compartment labels to cover all size fractions underneath
        compartment_labels = pivot_table.columns
        compartment_label_positions = np.arange(len(compartment_labels)) + 0.5
        plt.xticks(
            ticks=compartment_label_positions, labels=compartment_labels, rotation=90
        )

        # Set MP_Form and Size_Fraction_um labels
        row_labels = [
            f"{mp_form} - {size_frac_um}" for mp_form, size_frac_um in pivot_table.index
        ]
        row_label_positions = np.arange(len(pivot_table.index)) + 0.5
        plt.yticks(ticks=row_label_positions, labels=row_labels, rotation=0)
        titlename = (
            "Heatmap of log10 ("
            + fraction
            + " by MP_Form, Compartment, and Size_Fraction_um"
        )
        plt.title(titlename)
        plt.xlabel("Compartment", fontsize=14)
        plt.ylabel("MP_Form - Size_Fraction_um", fontsize=14)
        plt.tight_layout()

        fig = plt.gcf()
        plt.show()

        return fig  # , titlename

    def generate_flows_dict(self):
        for unit in ["mass", "number"]:
            if unit == "mass":
                tables_inputFlows = self.tables_inputFlows_mass
                tables_outputFlows = self.tables_outputFlows_mass
            elif unit == "number":
                tables_inputFlows = self.tables_inputFlows_number
                tables_outputFlows = self.tables_outputFlows_number
            else:
                raise ValueError("Unit must be 'mass' or 'number'.")
            flows_dict = dict()
            flows_dict["input_flows"] = {}
            flows_dict["output_flows"] = {}

            # Decode index in input and output flow tables
            for comp in tables_outputFlows.keys():
                df1 = tables_outputFlows[comp].copy()
                MP_size_df1 = []
                MP_form_df1 = []
                for x in df1.index:
                    MP_size_df1.append(self.model.size_dict[x[0]])
                    MP_form_df1.append(self.model.MP_form_dict_reverse[x[1:2]])

                df1.insert(0, "MP_size", MP_size_df1)
                df1.insert(1, "MP_form", MP_form_df1)
                flows_dict["output_flows"][comp] = df1

            for comp in tables_inputFlows:
                df2 = tables_inputFlows[comp].copy()
                MP_size_df2 = []
                MP_form_df2 = []
                for y in df2.index:
                    MP_size_df2.append(self.model.size_dict[y[0]])
                    MP_form_df2.append(self.model.MP_form_dict_reverse[y[1:2]])
                df2.insert(0, "MP_size", MP_size_df2)
                df2.insert(1, "MP_form", MP_form_df2)
                flows_dict["input_flows"][comp] = df2
            if unit == "mass":
                self.flows_dict_mass = flows_dict
            else:
                self.flows_dict_number = flows_dict

    def estimate_flows(self):

        self.surfComp_list = [c for c in self.model.dict_comp if "Surface" in c]
        """Estimate flows corresponding to each mode process based on the model results."""
        # Outflows ( in mass and particle number)
        for p in self.model.system_particle_object_list:
            p.outFlow_mass_g_s = {}
            p.outFlow_number_g_s = {}
            for c in p.RateConstants:
                if type(p.RateConstants[c]) == list:
                    p.outFlow_mass_g_s[c] = [
                        R * p.Pmass_g_SS for R in p.RateConstants[c]
                    ]
                    p.outFlow_number_g_s[c] = [
                        R * p.Pnum_SS for R in p.RateConstants[c]
                    ]
                else:
                    p.outFlow_mass_g_s[c] = p.RateConstants[c] * p.Pmass_g_SS
                    p.outFlow_number_g_s[c] = p.RateConstants[c] * p.Pnum_SS

        # Tables of output flows per compartmet
        tables_outputFlows_mass = {}
        tables_outputFlows_number = {}
        for c in list(self.model.dict_comp.keys()):
            part_dic_mass = {}
            part_dic_number = {}
            for p in self.model.system_particle_object_list:
                if p.Pcompartment.Cname == c:
                    part_dic_mass[p.Pcode] = pd.DataFrame.from_dict(
                        p.outFlow_mass_g_s, orient="index"
                    )
                    part_dic_number[p.Pcode] = pd.DataFrame.from_dict(
                        p.outFlow_number_g_s, orient="index"
                    )
            tables_outputFlows_mass[c] = pd.concat(part_dic_mass, axis=1).transpose()
            tables_outputFlows_number[c] = pd.concat(
                part_dic_number, axis=1
            ).transpose()

        for k in tables_outputFlows_mass:
            tables_outputFlows_mass[k] = (
                tables_outputFlows_mass[k].reset_index(level=1).drop("level_1", axis=1)
            )
            tables_outputFlows_number[k] = (
                tables_outputFlows_number[k]
                .reset_index(level=1)
                .drop("level_1", axis=1)
            )
        self.tables_outputFlows_mass = tables_outputFlows_mass
        self.tables_outputFlows_number = tables_outputFlows_number

        # Inflows: Tables of recieving flows through transport from other compartments
        tables_inputFlows_mass = {}
        tables_inputFlows_number = {}
        for comp in list(self.model.dict_comp.keys()):
            comp_input_flows_mass = []
            comp_input_flows_num = []
            for e_comp in self.model.dict_comp:
                if comp in self.model.dict_comp[e_comp].connexions:
                    inpProc = self.model.dict_comp[e_comp].connexions[comp]
                    if (
                        type(inpProc) == list
                    ):  # When there is more than one process of inflow into the compartment
                        df_inflows = tables_outputFlows_mass[e_comp].loc[
                            :, ["k_" + ele for ele in inpProc]
                        ]
                        df_inflows_num = self.tables_outputFlows_number[e_comp].loc[
                            :, ["k_" + ele for ele in inpProc]
                        ]

                        for proc in inpProc:
                            if proc == "dry_deposition":
                                position = self.surfComp_list.index(comp)
                                df_inflows["k_" + proc] = df_inflows["k_" + proc].apply(
                                    lambda x: x[position] if isinstance(x, list) else x
                                )
                                df_inflows_num["k_" + proc] = df_inflows_num[
                                    "k_" + proc
                                ].apply(
                                    lambda x: x[position] if isinstance(x, list) else x
                                )

                            elif proc == "mixing":

                                if (
                                    e_comp == "Ocean_Mixed_Water"
                                    and comp == "Ocean_Surface_Water"
                                ):
                                    df_inflows["k_" + proc] = df_inflows[
                                        "k_" + proc
                                    ].apply(
                                        lambda x: x[0] if isinstance(x, list) else x
                                    )
                                    df_inflows_num["k_" + proc] = df_inflows_num[
                                        "k_" + proc
                                    ].apply(
                                        lambda x: x[0] if isinstance(x, list) else x
                                    )

                                elif (
                                    e_comp == "Ocean_Mixed_Water"
                                    and comp == "Ocean_Column_Water"
                                ):
                                    df_inflows["k_" + proc] = df_inflows[
                                        "k_" + proc
                                    ].apply(
                                        lambda x: x[1] if isinstance(x, list) else x
                                    )
                                    df_inflows_num["k_" + proc] = df_inflows_num[
                                        "k_" + proc
                                    ].apply(
                                        lambda x: x[1] if isinstance(x, list) else x
                                    )
                                else:
                                    pass
                                # Revisit for percollation and tillage
                            else:
                                pass
                        comp_input_flows_mass.append(df_inflows)
                        comp_input_flows_num.append(df_inflows_num)

                    else:
                        df_inflows = (
                            tables_outputFlows_mass[e_comp]
                            .loc[:, "k_" + inpProc]
                            .to_frame()
                        )
                        df_inflows_num = (
                            tables_outputFlows_number[e_comp]
                            .loc[:, "k_" + inpProc]
                            .to_frame()
                        )
                        for ele in df_inflows["k_" + inpProc]:
                            if type(ele) == list:
                                connecting_comp = {
                                    key: value
                                    for key, value in self.model.dict_comp[
                                        e_comp
                                    ].connexions.items()
                                    if value == inpProc
                                }
                                poss_dict = {
                                    key: index
                                    for index, key in enumerate(connecting_comp.keys())
                                }
                                possition = poss_dict[comp]
                                df_inflows["k_" + inpProc] = df_inflows[
                                    "k_" + inpProc
                                ].apply(
                                    lambda x: x[possition] if isinstance(x, list) else x
                                )
                                df_inflows_num["k_" + inpProc] = df_inflows_num[
                                    "k_" + inpProc
                                ].apply(
                                    lambda x: x[possition] if isinstance(x, list) else x
                                )

                            else:
                                pass
                        comp_input_flows_mass.append(df_inflows)
                        comp_input_flows_num.append(df_inflows_num)
                else:
                    pass

            tables_inputFlows_mass[comp] = pd.concat(comp_input_flows_mass).fillna(0)
            tables_inputFlows_number[comp] = pd.concat(comp_input_flows_num).fillna(0)

        self.tables_inputFlows_mass = tables_inputFlows_mass
        self.tables_inputFlows_number = tables_inputFlows_number

    def extract_results_by_compartment(self):
        if self.Results_extended is None:
            raise ValueError(
                "Mass and particle number fractions not extracted. Call process_results() first."
            )
        mass_g = []
        particle_number = []
        mass_frac_100 = []
        num_frac_100 = []
        mass_conc_g_m3 = []
        num_conc = []
        for comp in list(self.model.dict_comp.keys()):
            mass_g.append(
                sum(
                    self.Results_extended[self.Results_extended["Compartment"] == comp][
                        "mass_g"
                    ]
                )
            )
            particle_number.append(
                sum(
                    self.Results_extended[self.Results_extended["Compartment"] == comp][
                        "number_of_particles"
                    ]
                )
            )
            mass_frac_100.append(
                sum(
                    self.Results_extended[self.Results_extended["Compartment"] == comp][
                        "mass_fraction"
                    ]
                )
                * 100
            )
            num_frac_100.append(
                sum(
                    self.Results_extended[self.Results_extended["Compartment"] == comp][
                        "number_fraction"
                    ]
                )
                * 100
            )
            mass_conc_g_m3.append(
                sum(
                    self.Results_extended[self.Results_extended["Compartment"] == comp][
                        "concentration_g_m3"
                    ]
                )
            )
            num_conc.append(
                sum(
                    self.Results_extended[self.Results_extended["Compartment"] == comp][
                        "concentration_num_m3"
                    ]
                )
            )

        results_by_comp = pd.DataFrame(columns=["Compartments"])
        results_by_comp["Compartments"] = list(self.model.dict_comp.keys())
        results_by_comp["mass_g"] = mass_g
        results_by_comp["number_of_particles"] = particle_number
        results_by_comp["%_mass"] = mass_frac_100
        results_by_comp["%_number"] = num_frac_100
        results_by_comp["Concentration_g_m3"] = mass_conc_g_m3
        results_by_comp["Concentration_num_m3"] = num_conc

        # self.results_by_comp = results_by_comp
        # return results_by_comp
        """Calculate inflows and outflows (mass and number) by compartment and update results_by_comp."""
        inflows_mass_list = []
        inflows_num_list = []
        outflows_mass_list = []
        outflows_num_list = []

        for n in range(len(results_by_comp)):
            compartment = results_by_comp.iloc[n]["Compartments"]

            # Calculate inflows and outflows for mass
            inflows_mass = process_flows_comp(
                compartment, "input_flows", self.flows_dict_mass
            )
            outflows_mass = process_flows_comp(
                compartment, "output_flows", self.flows_dict_mass
            )
            inflows_mass_list.append(inflows_mass)
            outflows_mass_list.append(outflows_mass)

            # Calculate inflows and outflows for number
            inflows_num = process_flows_comp(
                compartment, "input_flows", self.flows_dict_number
            )
            outflows_num = process_flows_comp(
                compartment, "output_flows", self.flows_dict_number
            )
            inflows_num_list.append(inflows_num)
            outflows_num_list.append(outflows_num)

        # Update the Results_extended DataFrame with the calculated flows
        results_by_comp["inflows_g_s"] = inflows_mass_list
        results_by_comp["inflows_num_s"] = inflows_num_list
        results_by_comp["outflows_g_s"] = outflows_mass_list
        results_by_comp["outflows_num_s"] = outflows_num_list
        results_by_comp["Total_inflows_g_s"] = [
            sum(results_by_comp.iloc[i].inflows_g_s.values())
            for i in range(len(results_by_comp))
        ]
        results_by_comp["Total_inflows_num_s"] = [
            sum(results_by_comp.iloc[i].inflows_num_s.values())
            for i in range(len(results_by_comp))
        ]
        results_by_comp["Total_outflows_g_s"] = [
            sum(results_by_comp.iloc[i].outflows_g_s.values())
            for i in range(len(results_by_comp))
        ]
        results_by_comp["Total_outflows_num_s"] = [
            sum(results_by_comp.iloc[i].outflows_num_s.values())
            for i in range(len(results_by_comp))
        ]

        self.results_by_comp = results_by_comp
        self.processed_results["results_by_comp"] = results_by_comp

    def create_rateConstants_table(self):
        df_dict = {
            "Compartment": [],
            "MP_form": [],
            "Size_Bin": [],
            "Rate_Constants": [],
        }

        for p in self.model.system_particle_object_list:
            df_dict["Compartment"].append(p.Pcompartment.Cname)
            df_dict["MP_form"].append(p.Pform)
            df_dict["Size_Bin"].append(p.Pname[:3])
            df_dict["Rate_Constants"].append(p.RateConstants)

        df = pd.DataFrame(df_dict)
        df2 = df["Rate_Constants"].apply(pd.Series)
        df = df.drop(columns="Rate_Constants")
        df3 = pd.concat([df, df2], axis=1)
        self.RC_df = df3
        self.processed_results["RateConstants_df"] = df3

    def plot_rateConstants(self):
        def sum_if_list(value):
            """Returns the sum of a list if the input is a list, otherwise returns the value itself."""
            return sum(value) if isinstance(value, list) else value

        """(FIX RC for wet deposition, now its given as a list of rate constants per surface compartment only for dry deposition and wet depossition is turned off)This needs to be fixed also for the matrix of interactions and estimation of flows"""
        rateConstants_df = self.RC_df.fillna(0)
        selected_columns = rateConstants_df.columns[3:]
        data_raw = self.RC_df[selected_columns]
        selected_data = data_raw.applymap(sum_if_list)
        log_data = selected_data.applymap(lambda x: np.log10(x) if x > 0 else np.nan)

        # Violin Plot
        plt.figure(figsize=(10, 6))
        sns.violinplot(data=log_data)
        # plt.yscale('log')
        plt.xticks(rotation=90)
        plt.title("Distribution of rate constants as log(k_s-1)")
        plt.show()
        fig = plt.gcf()
        self.processed_results["RC_violin_plot"] = fig

    def plot_compartment_distribution(
        self, mass_or_number
    ):  # mass_or_number: "%_mass" or ""%_number""
        """Bar chart plot of the mass or particle number distribution of particles by compartment."""
        compartment_colors = {
            "Ocean_Surface_Water": "#756bb1",
            "Ocean_Mixed_Water": "#756bb1",
            "Ocean_Column_Water": "#756bb1",
            "Coast_Surface_Water": "#2c7fb8",
            "Coast_Column_Water": "#2c7fb8",
            "Surface_Freshwater": "#9ebcda",
            "Bulk_Freshwater": "#9ebcda",
            "Sediment_Freshwater": "#fdae6b",
            "Sediment_Ocean": "#fdae6b",
            "Sediment_Coast": "#fdae6b",
            "Beaches_Soil_Surface": "#ffeda0",
            "Beaches_Deep_Soil": "#ffeda0",
            "Background_Soil_Surface": "#e5f5e0",
            "Background_Soil": "#e5f5e0",
            "Impacted_Soil_Surface": "#d95f0e",
            "Impacted_Soil": "#d95f0e",
            "Air": "#deebf7",
        }

        # Add your dataframe sorting and rounding logic
        df = (
            self.results_by_comp[["Compartments", mass_or_number]]
            .round(2)
            .sort_values(by=mass_or_number, ascending=False)
        )

        # Get the list of colors based on the Compartments in the df
        bar_colors = df["Compartments"].map(compartment_colors)

        # Plot the horizontal bar chart
        plt.barh(
            df["Compartments"], df[mass_or_number], color=bar_colors
        )  # Apply specific colors to each bar
        plt.xlabel(mass_or_number)
        plt.ylabel("Compartments")
        plt.title(f"{mass_or_number} Distribution by Compartment")
        plt.gca().invert_yaxis()  # Invert y-axis to match sorting order
        plt.show()
        fig = plt.gcf()
        return fig

    def process_all(self):
        """Runs all processing steps in order automatically and stores the results in the processed_results dictionary of the the class."""
        self.create_rateConstants_table()
        self.plot_rateConstants()
        self.estimate_flows()
        self.generate_flows_dict()
        self.process_results()
        fraction_heatmaps = dict()
        for fraction in ["mass_fraction", "number_fraction"]:
            fraction_heatmaps[fraction] = self.plot_fractionDistribution_heatmaps(
                fraction
            )
        self.processed_results["fraction_heatmaps"] = fraction_heatmaps
        self.extract_results_by_compartment()
        comp_distribution_barcharts = dict()
        for fraction in ["%_mass", "%_number"]:
            comp_distribution_barcharts[fraction] = self.plot_compartment_distribution(
                fraction
            )
        self.processed_results["comp_distribution_barcharts"] = (
            comp_distribution_barcharts
        )

        # Calculate exposure indicators
        self.estimate_exposure_indicators()
        self.estimate_emission_fractions()

    def estimate_exposure_indicators(self):
        """Estimate overall size dependent exposure indicators"""
        (
            self.processed_results["Overall_exposure_indicators"],
            self.processed_results["size_fraction_indicators"],
        ) = Exposure_indicators_calculation(self)

    def estimate_emission_fractions(self):
        """Estimate mass emission fractions:
        - Environmentally Dispersed Fraction (ϕ1): quantifies the relative extent to which the pollutants (MPs) can reach remote regions.
        - Remotely transferred fraction of mass (ϕ2) expresses the relative extent to which the MPs are (net) transferred to the target remote compartment following environmental dispersion to the remote region.
        """
        (
            self.processed_results["emission_fractions_mass_data"],
            self.processed_results["emission_fractions_mass_figure"],
        ) = estimate_emission_fractions(self)

    # def generate_pdf_report(self):
    #     # Create and populate the PDF
    #     pdf = PDFReport()
    #     pdf.add_page()

    #     for name, result in self.processed_results.items():
    #         pdf.set_font("Arial", "B", 12)
    #         pdf.cell(0, 10, name, ln=True)
    #         if isinstance(result, dict):
    #             pdf.add_dict(result)
    #         elif isinstance(result, pd.DataFrame):
    #             pdf.add_dataframe(result)
    #         elif isinstance(result, plt.Figure):
    #             pdf.add_figure(result)
    #         else:
    #             pdf.add_text(str(result))

    #     pdf.output("results_report.pdf")
