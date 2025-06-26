    """
    This module contains the input and configuration data examples for the Utopia package.
    """
    import numpy as np
    
    #Example model config with all available variables
    full_config={
        "big_bin_diameter_um": 5000,
        "N_sizeBins": 5,
        "vol_algal_cell_m3": 2.0e-16,# REF: Kooi et al. (2017)
        "spm_density_kg_m3": 1388, # REF: Kooi et al. (2017),
        "boxName": "Utopia",
        "compartment_names": [],
        "comp_input_file_name": "inputs_compartments.csv",
        "comp_interactions_file_name": "compartment_interactions.csv",
        "alpha_heter_filename": "alpha_heter.csv",
        "compartment_types":{ "UTOPIA_surfaceSea_water_compartments":["Ocean_Surface_Water", "Coast_Surface_Water"],"UTOPIA_water_compartments":[
    "Ocean_Mixed_Water",
    "Ocean_Column_Water",
    "Coast_Column_Water",
    "Surface_Freshwater",
    "Bulk_Freshwater",
],"UTOPIA_deep_soil_compartments":[
    "Beaches_Deep_Soil",
    "Background_Soil",
    "Impacted_Soil",
],"UTOPIA_soil_surface_compartments":[
    "Beaches_Soil_Surface",
    "Background_Soil_Surface",
    "Impacted_Soil_Surface",
],"UTOPIA_sediment_compartment":[
    "Sediment_Freshwater",
    "Sediment_Ocean",
    "Sediment_Coast",
],"UTOPIA_air_compartments":["Air"]}
        }
    
    #: Example model data
    full_data={
        "MPdensity_kg_m3": 980,
        "MP_composition": "PE",
        "shape": "sphere",
        "FI": 0.5,
        "t_half_deg_free":66000,  # in days (10 times slower than the rate of degradation (to form dissolved organics) shown in Pfohl et al. 2023 for TPU-arom),
        "heter_deg_factor": 10, #model assumption
        "biof_deg_factor":1/2, #model assumption
        "factor_deepWater_soilSurface": 10, #model assumption
        "factor_sediment": 100, #model assumption
        "t_frag_gen_FreeSurfaceWater": 36.5, #in days
        "biof_frag_factor":2, #model assumption
        "heter_frag_factor":100, #model assumption,
        "emiss_dict_g_s": {'Ocean_Surface_Water': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 250000}, 'Ocean_Mixed_Water': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Ocean_Column_Water': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Coast_Surface_Water': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Coast_Column_Water': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Surface_Freshwater': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Bulk_Freshwater': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Sediment_Freshwater': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Sediment_Ocean': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Sediment_Coast': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Beaches_Soil_Surface': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Beaches_Deep_Soil': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Background_Soil_Surface': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Background_Soil': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Impacted_Soil_Surface': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Impacted_Soil': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}, 'Air': {'a': 0, 'b': 0, 'c': 0, 'd': 0, 'e': 0}},}
    
    
    