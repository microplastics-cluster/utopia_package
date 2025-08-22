# Standardized procedure for integration of sub-models into UTOPIA 

Utopia currently integrates 21 sub-process models described in src\utopia\preprocessing\RC_generator.py  

For each of these processes UTOPIA derives sub-proces reate constants for each model compartment, size fraction and aggregation state combination (17x5x4=340).

Currently UTOPIA considers the following input parameters that determine the calculation of the process sub-model rate contats claculation:

- Size
- Aggregation state
- Shape
- Density
- Fragmentation profile (given by FI)
- Fragmentation/discorporation timescale

New process sub-models to be incorporated in UTOPIA should be able to integrate these properties or be compatible with these inputs.

In order to incorporate a new process sub-model go trougth the following steps:

# STEP 1: Process sub-models procedural check-list

Go through the procedural checklist to ensure you have all information needed:

 - From the model inputs parameters (MP properties) what is the coverage of the process sub-model? Define limits of applicability of model
  
 - Provide dictionary of input parameters per compartment and MP form and size fraction when relevant
  
 - Indicate assumptions followed and provide references for set values and or process formulation
  
 - Provide uncertainty ranges for process-sub model parameters (type of distribution and range)

# STEP 2: Generate new entry for the model sub-process in the src\utopia\preprocessing\RC_generator.py file following the template code:

def new_process_name(particle):
    # describe the process providing references and indicating the assumption followed

    # The process code will calculate a particle specific rate constant for the process (size and aggregation state). Which will be also specific to the compartment where the particle is in

    # write process code to return rate conatnt value in s-1
    
    return k_new_process_name

# STEP 3: Assign process to the compartment where it takes place

Any new process defined in the model has to be assigned to the compartment or compartments where it will take place. This is done by adding the process name to the list of processes of the sub-compartment class in the compartment_clases.py function under the object folder of the model. 

As example in the function compartment_clases.py we can add further processes to the different compartment classes by including the process name into the attribute processes list:

class compartment_deep_soil(Compartment):
    def __init__(
        self,
        Cname,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.processess = [
            "discorporation",
            "fragmentation",
            "sequestration_deep_soils",
            "soil_convection", “new_process_name”
        ]


# STEP 4: Include the process in the compartments interactions matrix (compartment_interactions.csv)

This file specifies the interactions between compartments and is limited to transport processes. If the new process submodel is not a transport process this step can be skipped.

# STEP 5: In case a new input parameter is needed for the new process-sub model, this parameters shoild be added in the configuration and data input files needed to run the model and the parameter would have to be added as attibute to be asigned to the particle object as well as in the UTOPIAModel class validation functions.

--> Add input to config_data.json file as new entry
--> Assing new attribute in particle object class
--> Add new entry to the class validation function

