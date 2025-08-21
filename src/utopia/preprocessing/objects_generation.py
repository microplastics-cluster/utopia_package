"""Contains all the functions to generate the model objects: model box, model compartments and the model particles (MP ans SPM)"""

import os
import math
import copy
from pathlib import Path
from utopia.globalConstants import *
from utopia.objects.compartment_classes import *
from utopia.objects.particulate_classes import *
from utopia.objects.box_class import *

from utopia.preprocessing.readinputs_from_csv import *


def generate_objects(model):
    """Function for generating the UTOPIA model objects: model box, model compartments and the model particles"""
    # Boxes
    UTOPIA = Box(model.boxName)
    # print(f"The model box {boxName} has been created")

    modelBoxes = [UTOPIA]
    # modelBoxes=instantiateBoxes_from_csv(boxFile)
    boxNames_list = [b.Bname for b in modelBoxes]

    # Compartmets
    """Call read imput file function for compartments"""
    # base_path = Path(__file__).resolve().parent / "data"
    inputs_path_file = model.base_path / model.comp_input_file_name

    # inputs_path_file = f"data/{model.comp_input_file_name}"
    compartments = instantiate_compartments(
        inputs_path_file=inputs_path_file, compartment_types=model.compartment_types
    )

    # Establish connexions between compartments defining their interaction mechanism: only listed those compartments wich will recieve particles from the define compartment. i.e. the ocean surface water compartment transports particles to the ocean mix layer through settling and to air through sea spray resuspension
    connexions_path_file = model.base_path / model.comp_interactFile_name
    set_interactions(compartments, connexions_path_file)

    # Assign modelling code to compartments
    for c in range(len(compartments)):
        compartments[c].Ccode = c + 1

    ##Calculate compartments volume
    for c in compartments:
        c.calc_volume()

    ## Dictionary of compartments
    dict_comp = {
        item.Cname: item for item in compartments
    }  # Before the compartments association to RS...migth need to come after to also reflect the CBox connexion

    compartmentNames_list = [item.Cname for item in compartments]

    # PARTICLES

    ##Free microplastics (freeMP)

    # MP_freeParticles = instantiateParticles_from_csv(inputs_path + mp_imputFile_name)

    MP_freeParticles = generate_particles_from_df(model.particles_df)

    dict_size_coding = dict(
        zip(
            [p.Pname for p in MP_freeParticles],
            [p.diameter_um for p in MP_freeParticles],
        )
    )

    ###Calculate freeMP volume
    for i in MP_freeParticles:
        i.calc_volume()
        # print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")

    ##Biofouled microplastics (biofMP)
    spm = Particulates(
        Pname="spm1",
        Pform="suspendedParticulates",
        Pcomposition="Mixed",
        Pdensity_kg_m3=model.spm_density_kg_m3,
        Pshape="sphere",
        PdimensionX_um=model.spm_radius_um*2,
        PdimensionY_um=0,
        PdimensionZ_um=0,
    )
    spm.calc_volume()
    # print(f"spm Volume: {spm.Pvolume_m3} m3")
    # print(f"Density of spm: {spm.Pdensity_kg_m3} kg_m3")

    MP_biofouledParticles = []
    for i in MP_freeParticles:
        MP_biofouledParticles.append(ParticulatesBF(parentMP=i, spm=spm))
    # print(
    #     f"The biofouled MP particles {[p.Pname for p in MP_biofouledParticles]} have been generated"
    # )

    ###Calculate biofMP volume
    for i in MP_biofouledParticles:
        i.calc_volume()
        # print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")

    ##Heteroaggregated microplastics (heterMP)

    MP_heteroaggregatedParticles = []
    for i in MP_freeParticles:
        MP_heteroaggregatedParticles.append(ParticulatesSPM(parentMP=i, parentSPM=spm))
    # print(
    #     f"The heteroaggregated MP particles {[p.Pname for p in MP_heteroaggregatedParticles]} have been generated"
    # )

    ###Calculate heterMP volume
    for i in MP_heteroaggregatedParticles:
        i.calc_volume_heter(i.parentMP, spm)
        # print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")

    ##Biofouled and Heteroaggregated microplastics (biofHeterMP)
    MP_biofHeter = []
    for i in MP_biofouledParticles:
        MP_biofHeter.append(ParticulatesSPM(parentMP=i, parentSPM=spm))
    # for i in MP_biofHeter:
    #     print(f"Density of {i.Pname}: {i.Pdensity_kg_m3} kg_m3")
    # print(
    #     f"The biofouled and heteroaggregated MP particles {[p.Pname for p in MP_biofHeter]} have been generated"
    # )

    ###Calculate biofHeterMP volume
    for i in MP_biofHeter:
        i.calc_volume_heter(i.parentMP, spm)

    particles = (
        MP_freeParticles
        + MP_biofouledParticles
        + MP_heteroaggregatedParticles
        + MP_biofHeter
    )

    particles_properties = {
        "Particle": ([p.Pname for p in particles]),
        "Radius_m": ([p.radius_m for p in particles]),
        "Volume_m3": ([p.Pvolume_m3 for p in particles]),
        "Density_kg_m3": ([p.Pdensity_kg_m3 for p in particles]),
        "Corey Shape Factor": ([p.CSF for p in particles]),
    }

    particles_properties_df = pd.DataFrame(data=particles_properties)

    # Assign compartmets to UTOPIA

    for comp in compartments:
        UTOPIA.add_compartment(
            copy.deepcopy(comp)
        )  # Check if the use of copy is correct!!

    # print(
    #     f"The compartments {[comp.Cname for comp in UTOPIA.compartments]} have been assigned to {UTOPIA.Bname } model box"
    # )

    # Estimate volume of UTOPIA box by adding volumes of the compartments addedd
    # UTOPIA.calc_Bvolume_m3() #currently volume of soil and air boxess are missing, to be added to csv file

    # Add particles to compartments
    for b in modelBoxes:
        for c in b.compartments:
            for p in particles:
                c.add_particles(copy.deepcopy(p))
        # print(f"The particles have been added to the compartments of {b.Bname}")

    # List of particle objects in the system:
    system_particle_object_list = []

    for b in modelBoxes:
        for c in b.compartments:
            for freeMP in c.particles["freeMP"]:
                system_particle_object_list.append(freeMP)
            for heterMP in c.particles["heterMP"]:
                system_particle_object_list.append(heterMP)
            for biofMP in c.particles["biofMP"]:
                system_particle_object_list.append(biofMP)
            for heterBiofMP in c.particles["heterBiofMP"]:
                system_particle_object_list.append(heterBiofMP)

    # Generate list of species names and add code name to object
    SpeciesList = generate_system_species_list(
        system_particle_object_list,
        model.MPforms_list,
        compartmentNames_list,
        boxNames_list,
    )

    return (
        system_particle_object_list,
        SpeciesList,
        spm,
        dict_comp,
        particles_properties_df,
    )


###CONTINUE HERE###
