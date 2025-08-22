# Standardized Procedure for Integration of Sub-Models into UTOPIA  

UTOPIA currently integrates **21 sub-process models** described in  
`src/utopia/preprocessing/RC_generator.py`.  

For each of these processes, UTOPIA derives sub-process rate constants for each model compartment, size fraction, and aggregation state combination (**17 × 5 × 4 = 340**).  

Currently, UTOPIA considers the following input parameters that determine the calculation of the process sub-model rate constants:

- Size
- Aggregation state
- Shape
- Density
- Fragmentation profile (given by **FI**)
- Fragmentation/discorporation timescale

New process sub-models to be incorporated into UTOPIA should be able to integrate these properties or be compatible with these inputs.  

---

## STEP 1: Process Sub-Models Procedural Checklist  

Go through the checklist to ensure you have all information needed:

- From the model input parameters (MP properties), what is the coverage of the process sub-model? Define limits of applicability of the model.  
- Provide a dictionary of input parameters per compartment and MP form and size fraction (when relevant).  
- Indicate assumptions followed and provide references for set values and/or process formulation.  
- Provide uncertainty ranges for process sub-model parameters (type of distribution and range).  

---

## STEP 2: Generate a New Entry for the Model Sub-Process  

Create a new entry in  
`src/utopia/preprocessing/RC_generator.py`  
following this template:  

```python
def new_process_name(particle):
    """
    Describe the process, provide references, and indicate the assumptions followed.

    The process code should calculate a particle-specific rate constant
    (based on size and aggregation state). It will also be specific to the
    compartment where the particle is located.
    """

    # Example calculation of rate constant
    k_new_process_name = ...  

    return k_new_process_name

```
## STEP 3: Assign the Process to the Relevant Compartment

Any new process defined in the model has to be assigned to the compartment(s) where it takes place.
This is done by adding the process name to the list of processes of the sub-compartment class in
compartment_clases.py (under the object folder).

Example:

```python
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
        self.processes = [
            "discorporation",
            "fragmentation",
            "sequestration_deep_soils",
            "soil_convection",
            "new_process_name"
        ]
```
## STEP 4: Update the Compartment Interactions Matrix

Include the process in the `compartment_interactions.csv` file.
This file specifies the interactions between compartments and is limited to transport processes.
If the new process sub-model is not a transport process, this step can be skipped.

## STEP 5: Add Any New Input Parameters

If a new input parameter is required, it must be added in the configuration and data input files needed to run the model. The parameter should also be added as an attribute to the particle object and included in the UTOPIAModel class validation functions.

- Add input to config_data.json as a new entry.

- Assign the new attribute in the particle object class.

- Add the new entry to the class validation function.