import copy
import string
import json
import pandas as pd
import os
from datetime import datetime
import math
import numpy as np
from pathlib import Path
from utopia.preprocessing.objects_generation import *
from utopia.preprocessing.generate_rate_constants import *
from utopia.preprocessing.fill_interactions_df import *
from utopia.solver_steady_state import *

import json


class utopiaModel:
    """The class that controls usage of the UTOPIA model
    Parameters
    ----------
    config : dict
        Configuration dictionary with model configuration options regarding to parameterization of the unit world compartments
    data : dict
        Data dictionary with input data required to run the model. This contains the selected options from the imput menues of the user interface including microplastics properties and emission scenario selected (in the future to also include compartment properties or  the selected compartment parameterization scenario).
    validate : bool. deafult=True
        Validate input data and configuration options
    """

    def __init__(self, config: None, data: None, validate=True):
        """Initialize the UTOPIA model"""
        self.base_path = Path(__file__).resolve().parent / "data"

        # Load default config and data if not provided
        if config is None:
            config = self.load_json_file(self.base_path / "default_config.json")
        if data is None:
            data = self.load_json_file(self.base_path / "default_data.json")

        # Validate the config and data inputs
        if validate:
            self.validate_inputs(config, data)

        # Assign attributes
        self.config = config
        self.data = data

        # Load parameters
        self.load_parameters()
        self.particles_df = self.generate_particles_dataframe()
        self.generate_coding_dictionaries()

    @staticmethod
    def load_json_file(filepath):
        base_path = os.path.dirname(__file__)
        full_path = os.path.join(base_path, filepath)

        with open(full_path, "r", encoding="utf-8") as f:
            return json.load(f)

    @staticmethod
    def load_csv_column(filename, column_name):
        """Load a column from input CSV file: Reads a single column from a CSV file and returns it as a list"""
        base_path = Path(__file__).resolve().parent / "data"
        file_path = base_path / filename
        df = pd.read_csv(file_path, usecols=[column_name])
        return df[column_name].tolist()

    @staticmethod
    def check_required_keys(dictionary, required_keys, dict_name):
        """Checks if all required keys are present in a dictionary."""
        missing_keys = [key for key in required_keys if key not in dictionary]
        if missing_keys:
            raise KeyError(f"Missing keys in {dict_name}: {', '.join(missing_keys)}")

    def validate_inputs(self, config, data):
        """Validates input dictionaries to ensure all required keys exist."""
        required_data_keys = [
            "MPdensity_kg_m3",
            "MP_composition",
            "shape",
            "FI",
            "MP_form",
            "t_half_deg_free",
            "heter_deg_factor",
            "biof_deg_factor",
            "factor_deepWater_soilSurface",
            "factor_sediment",
            "t_frag_gen_FreeSurfaceWater",
            "biof_frag_factor",
            "heter_frag_factor",
            "emiss_dict_g_s",
        ]

        required_config_keys = [
            
            
            "vol_algal_cell_m3",
            "spm_density_kg_m3",
            "comp_input_file_name",
            "comp_interactFile_name",
            "boxName",
            "MPforms_list",
            "solver",
            "compartment_types",
        ]

        self.check_required_keys(data, required_data_keys, "data")
        self.check_required_keys(config, required_config_keys, "config")
        
        # Type and value checks
        if not isinstance(data["MPdensity_kg_m3"], (int, float)):
            raise TypeError("MPdensity_kg_m3 must be a number.")

        if data["MPdensity_kg_m3"] <= 0:
            raise ValueError("MPdensity_kg_m3 must be positive.")

    # Add more checks as needed (TO BE ADDED!)

    
    def modify_and_save_data(self, data, modifications, filename):
        """
        Modify the provided data dictionary with the given modifications and save to a JSON file.

        Parameters:
        - data: Original data dictionary.
        - modifications: Dictionary containing keys and new values to update in the data.
        - filename: Name of the JSON file to save the modified data.
        """
        # Apply modifications
        for key, value in modifications.items():
            if key in data:
                data[key] = value
            else:
                raise KeyError(f"Invalid key in modifications: {key}")

        # Save the modified data
        self.save_json_file(data, filename)
    
    def save_json_file(self, data, filename):
        output_path = self.base_path / filename
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)
        print(f"Modified data saved to {output_path}")
    
    def load_parameters(self):
        """Loads required parameters from config and data dictionaries."""

        # Microplastics physical properties
        self.MPdensity_kg_m3 = self.data["MPdensity_kg_m3"]
        self.MP_composition = self.data["MP_composition"]
        self.shape = self.data["shape"]
        self.MP_form = self.data["MP_form"]
        #self.big_bin_diameter_um = self.config["big_bin_diameter_um"]
        #self.N_sizeBins = self.config["N_sizeBins"]
        self.FI = self.data["FI"]
        self.t_half_deg_free = self.data["t_half_deg_free"]
        self.t_frag_gen_FreeSurfaceWater = self.data["t_frag_gen_FreeSurfaceWater"]

        self.heter_deg_factor = self.data["heter_deg_factor"]
        self.biof_deg_factor = self.data["biof_deg_factor"]
        self.factor_deepWater_soilSurface = self.data["factor_deepWater_soilSurface"]
        self.factor_sediment = self.data["factor_sediment"]

        self.biof_frag_factor = self.data["biof_frag_factor"]
        self.heter_frag_factor = self.data["heter_frag_factor"]

        # Environmental characteristics

        self.vol_algal_cell_m3 = self.config["vol_algal_cell_m3"]
        self.spm_density_kg_m3 = self.config["spm_density_kg_m3"]
        self.comp_input_file_name = self.config["comp_input_file_name"]
        self.comp_interactFile_name = self.config["comp_interactFile_name"]
        self.boxName = self.config["boxName"]

        self.MPforms_list = self.config["MPforms_list"]

        # Load parameters from config and data dictionaries
        self.MPforms_list = self.config["MPforms_list"]
        self.compartments_list = self.load_csv_column(
            self.comp_input_file_name, "Cname"
        )
        self.solver = self.config["solver"]
        self.compartment_types = self.config["compartment_types"]

        # Derived environmental parameters
        self.radius_algae_m = ((3.0 / 4.0) * (self.vol_algal_cell_m3 / math.pi)) ** (
            1.0 / 3.0
        )
        self.spm_radius_um = self.radius_algae_m * 1e6

        # Emission scenario
        self.emiss_dict_g_s = self.data["emiss_dict_g_s"]

    def generate_particles_dataframe(self):
        """Generates the microplastics input DataFrame from Utopia model attributes."""
        MPdensity_kg_m3 = self.MPdensity_kg_m3
        shape = self.shape
        # Number of size bins and sizes are now fixed ando should not be part of the input data files
        self.N_sizeBins = 5
        N_sizeBins = self.N_sizeBins
        self.size_distribution=[0.3, 3,30,750,15000] # median values of the particle diameter of the size bin in nanometers
        self.size_distribution_dict={
            "0.3": "Nanoplastic",
            "3.0": "< 10 µm",
            "30.0": "< 100 µm",
            "750.0": "< 5000 µm",
            "15000.0": "Mesoplastic",
        }
        size_distribution = self.size_distribution
        #big_bin_diameter_um = self.big_bin_diameter_um

        # Generate size distribution
        # size_distribution = [big_bin_diameter_um]
        # for _ in range(N_sizeBins - 1):
        #     size_distribution.append(size_distribution[-1] / 10)
        # size_distribution.reverse()

        # Only supports spherical particles for now
        if shape == "sphere":
            data = {
                "Name": [f"mp{i+1}" for i in range(N_sizeBins)],
                "form": ["freeMP"] * N_sizeBins,
                "shape": [shape] * N_sizeBins,
                "composition": [self.MP_composition] * N_sizeBins,
                "density_kg_m3": [MPdensity_kg_m3] * N_sizeBins,
                "dimensionX_um": [d / 2 for d in size_distribution],
                "dimensionY_um": [d / 2 for d in size_distribution],
                "dimensionZ_um": [d / 2 for d in size_distribution],
            }
            return pd.DataFrame(data)
        else:
            raise ValueError("Shape not supported yet")

    def generate_coding_dictionaries(self):
        """Generates Mp form, size and compartment coding dictionaries as attributes."""
        if self.particles_df is None:
            raise ValueError("Particles DataFrame has not been generated.")

        # Dictionary mapping particle names to sizes
        self.dict_size_coding = dict(
            zip(self.particles_df["Name"], self.particles_df["dimensionX_um"] * 2)
        )

        # Generate size codes (a-z based on number of bins)
        self.size_codes = list(string.ascii_lowercase[: self.N_sizeBins])

        # Dictionary mapping size codes to sizes
        self.size_dict = dict(zip(self.size_codes, self.dict_size_coding.values()))

        # Dictionary mapping MP form codes to MP forms
        self.particle_forms_coding = dict(zip(self.MPforms_list, ["A", "B", "C", "D"]))
        self.MP_form_dict_reverse = {
            v: k for k, v in self.particle_forms_coding.items()
        }

        # Dictionary mapping compartment names to compartment codes
        self.particle_compartmentCoding = dict(
            zip(
                self.compartments_list,
                list(range(len(self.compartments_list))),
            )
        )
        self.comp_dict_inverse = {
            v: k for k, v in self.particle_compartmentCoding.items()
        }

    def run(self):
        """Runs the UTOPIA model with the configured parameters."""
        # Generate model objects based on model configuration and input data
        print("Running UTOPIA model with configured parameters...")
        (
            self.system_particle_object_list,
            self.SpeciesList,
            self.spm,
            self.dict_comp,
            self.particles_properties_df,
        ) = generate_objects(self)
        print("Generated model objects.")

        # Estimate rate contants for all processess for each particle in the system
        generate_rate_constants(self)
        print("Generated rate constants for model particles.")

        # Build matrix of interactions
        self.interactions_df = fillInteractions_fun_OOP(
            system_particle_object_list=self.system_particle_object_list,
            SpeciesList=self.SpeciesList,
            dict_comp=self.dict_comp,
        )
        print("Built matrix of interactions.")
        # Solve system of ODEs
        if self.solver == "SteadyState":

            (self.R, self.PartMass_t0, self.input_flows_g_s, self.input_flows_num_s) = (
                solver_SS(self)
            )
            print("Solved system of ODEs for steady state.")
        else:
            raise ValueError("Solver not implemented yet")

        # Test that there are no negative results
        for i, idx in zip(self.R["mass_g"], self.R.index):
            if i < 0:
                print("negative values in the solution for " + idx)
            else:
                pass

    def summarize(self):
        """Prints a summary of the model's key parameters."""
        print(f"Model: UTOPIA")
        # print(f"Box Name: {self.boxName}")
        print(f"Microplastic Density (kg/m3): {self.MPdensity_kg_m3}")
        print("MP shape: ", self.shape)
        print("Emissions made to MP form: ", self.MP_form)

        def print_fragmentation_style(F):
            if F == 0:
                print("Fragmentation style: Erosive")
            elif F == 1:
                print("Fragmentation style: Sequential")
            else:
                print(f"Fragmentation style: Mixed (F = {F})")

        print_fragmentation_style(self.FI)
        print("Fragmetation timescale (days): ", self.t_frag_gen_FreeSurfaceWater)
        print("Discorporation timescale (days): ", self.t_half_deg_free)
        for compartment, size_fractions in self.emiss_dict_g_s.items():
            for fraction, value in size_fractions.items():
                if value > 0:
                    print(
                        f"Emissions to {compartment} for size fraction {self.size_dict[fraction]} {chr(181)}m: {value} g/s"
                    )
    def run_sensitivity_analysis(self, param_distributions, n_simulations=1000, **monaco_kwargs):
        """
        Run Monte Carlo sensitivity analysis using monaco.
        
        Args:
            param_distributions (dict): {param_name: distribution_object}
            n_simulations (int): Number of Monte Carlo runs
            **monaco_kwargs: Passed to monaco.Simulation
            
        Returns:
            monaco.Simulation: Full results object
        """
        from utopia.monte_carlo.analysis import run_monte_carlo
        from utopia.monte_carlo.distributions import MP_Gaussian, Emission_Uniform  # Your custom distributions
        
        # Create a wrapper that preserves model state
        def model_wrapper(**params):
            # Create a copy of the model with new params
            modified_model = self.copy()
            modified_model.data.update(params)  # Update input data
            modified_model.run()
            
            return {
                'outputs': modified_model.R,
                'indicators': modified_model.processed_results
            }
        
        
        
        return run_monte_carlo(
            model_func=model_wrapper,
            param_distributions=param_distributions,
            n_simulations=n_simulations,
            **monaco_kwargs
        )