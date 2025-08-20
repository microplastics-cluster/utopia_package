import math
import numpy as np


class Particulates:
    """Class Particulates generates particulate objects, especifically microplastic particle objects. The class defines a particle object by its composition, shape and dimensions"""

    # constructor
    def __init__(
        self,
        Pname,
        Pform,
        Pcomposition,
        Pdensity_kg_m3,
        Pshape,
        PdimensionX_um,
        PdimensionY_um,
        PdimensionZ_um,
        t_half_d=5000,
        Pnumber_t0=None,
    ):
        self.Pname = Pname
        self.Pform = Pform  # Pform has to be in the particles type list: ["freeMP",""heterMP","biofMP","heterBiofMP"]
        self.Pcomposition = Pcomposition
        self.Pdensity_kg_m3 = Pdensity_kg_m3
        self.Pshape = Pshape
        self.PdimensionX_um = PdimensionX_um  # shortest size
        self.PdimensionY_um = PdimensionY_um  # longest size
        self.PdimensionZ_um = PdimensionZ_um  # intermediate size
        self.PdimensionX_m = PdimensionX_um / 1000000  # shortest size
        self.PdimensionY_m = PdimensionY_um / 1000000  # longest size
        self.PdimensionZ_m = PdimensionZ_um / 1000000  # intermediate size
        self.Pnumber_t0 = Pnumber_t0  # number of particles at time 0. to be objetained from emissions and background concentration of the compartment
        
        # extract radius (for spherical particles) or calculate equivalent radius (for fibers)
        if self.Pshape == "sphere":
            self.radius_m = (
                self.PdimensionX_m / 2
            )  # In spherical particles from MP radius (x dimension)
        elif self.Pshape in {"fiber", "fibre", "cylinder"}:
            self.radius_m = (
                (3 / 2) * self.PdimensionX_m * self.PdimensionY_m * self.PdimensionZ_m
            ) ** (1 / 3) / 2
        else:
            print("Error: shape not supported yet")
            # print error message for shapes other than spheres
            
        self.diameter_m = self.radius_m * 2
        self.diameter_um = self.diameter_m * 1e6
        self.Pemiss_t_y = 0  # set as 0
        self.t_half_d = t_half_d

    def __repr__(self):
        return (
            "{"
            + self.Pname
            + ", "
            + self.Pform
            + ", "
            + self.Pcomposition
            + ", "
            + self.Pshape
            + ", "
            + str(self.Pdensity_kg_m3)
            + ", "
            + str(self.radius_m)
            + "}"
        )

    # methods

    def calc_volume(self):
        """Particle volume calculation. Different formulas for different particle shapes, currently defined for spheres, fibres, cylinders, pellets and irregular fragments"""

        if self.Pshape == "sphere":
            self.Pvolume_m3 = 4 / 3 * math.pi * (self.radius_m) ** 3
            # calculates volume (in m3) of spherical particles from MP radius (x dimension)
            self.CSF = 1
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)
            # print(
            #     "Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3"
            # )
            # print("Calculated Corey Shape Factor: " + str(self.CSF))

        elif (
            self.Pshape == "fibre"
            or self.Pshape == "fiber"
            or self.Pshape == "cylinder"
        ):
            self.Pvolume_m3 = math.pi * (self.PdimensionX_m) ** 2 * (self.PdimensionY_m)
            # calculates volume (in m3) of fibres or cylinders from diameter and
            # length assuming cylindrical shape where X is the shorterst size (radius) ans Y the longest (heigth)
            self.CSF = (self.PdimensionX_m) / math.sqrt(self.PdimensionY_m * self.PdimensionZ_m)
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)
            # print(
            #     "Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3"
            # )
            # print("Calculated Corey Shape Factor: " + str(self.CSF))

        elif self.Pshape == "pellet" or self.Pshape == "fragment":
            self.Pvolume_m3 = (
                self.PdimensionX_m * self.PdimensionY_m * self.PdimensionZ_m
            )
            # approximate volume calculation for irregular fragments
            # approximated as a cuboid using longest, intermediate and shortest length
            #!! Note: not sure if pellets fits best here or rather as sphere/cylinder
            # might adjust later!!
            self.CSF = self.PdimensionX_m / math.sqrt(
                self.PdimensionY_m * self.PdimensionZ_m
            )
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)
            # print(
            #     "Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3"
            # )
            # print("Calculated Corey Shape Factor: " + str(self.CSF))

        else:
            print("Error: unknown shape")
            # print error message for shapes other than spheres
            # (to be removed when other volume calculations are implemented)

    def calc_numConc(self, concMass_mg_L, concNum_part_L):

        if concNum_part_L == 0:
            self.concNum_part_m3 = (
                concMass_mg_L / 1000 / self.Pdensity_kg_m3 / self.Pvolume_m3
            )
            # if mass concentration is given, it is converted to number concentration
        else:
            self.concNum_part_m3 = concNum_part_L * 1000
            # if number concentration is given, it is converted from part/L to part/m3

    def assign_compartment(self, comp):
        self.Pcompartment = comp


class ParticulatesBF(Particulates):
    "This is a class to create ParticulatesBIOFILM objects"

    # class attribute
    species = "particulate"

    # constructor
    def __init__(self, parentMP, spm):

        self.Pname = parentMP.Pname + "_BF"
        self.Pcomposition = parentMP.Pcomposition
        self.Pform = "biofMP"
        self.parentMP = parentMP
        self.BF_density_kg_m3 = spm.Pdensity_kg_m3
        self.BF_thickness_um = spm.PdimensionX_um
        self.radius_m = parentMP.radius_m + (
            self.BF_thickness_um / 1e6
        )  # In spherical particles from MP radius (x dimension)
        self.diameter_m = self.radius_m * 2
        self.diameter_um = self.diameter_m * 1e6
        self.t_half_d = 25000  # As per The Full Multi parameterization
        if parentMP.PdimensionY_um == 0:
            self.PdimensionY_um = 0
        else:
            self.PdimensionY_um = parentMP.PdimensionY_um + self.BF_thickness_um * 2

        if parentMP.PdimensionZ_um == 0:
            self.PdimensionZ_um = 0
        else:
            self.PdimensionZ_um = parentMP.PdimensionZ_um + self.BF_thickness_um * 2

        if parentMP.PdimensionX_um == 0:
            self.PdimensionX_um = 0
        else:
            self.PdimensionX_um = parentMP.PdimensionX_um + self.BF_thickness_um * 2

        self.Pshape = (
            parentMP.Pshape
        )  # to be updated for biofilm, could argue that shape is retained (unlike for SPM-bound)
        self.Pdensity_kg_m3 = (
            self.parentMP.radius_m**3 * self.parentMP.Pdensity_kg_m3
            + (
                (self.parentMP.radius_m + (self.BF_thickness_um / 1e6)) ** 3
                - self.parentMP.radius_m**3
            )
            * self.BF_density_kg_m3
        ) / ((self.parentMP.radius_m + (self.BF_thickness_um / 1e6)) ** 3)
        # equation from Kooi et al for density

        self.PdimensionX_m = self.PdimensionX_um / 1000000  # shortest size
        self.PdimensionY_m = self.PdimensionY_um / 1000000  # longest size
        self.PdimensionZ_m = self.PdimensionZ_um / 1000000  # intermediate size


class ParticulatesSPM(Particulates):
    "This is a class to create ParticulatesSPM objects"

    # class attribute
    species = "particulate"

    # constructor
    def __init__(self, parentSPM, parentMP):

        self.Pname = parentMP.Pname + "_SPM"
        self.Pcomposition = parentMP.Pcomposition
        if parentMP.Pform == "biofMP":
            self.Pform = "heterBiofMP"
            self.t_half_d = 50000  # As per The Full multi parameterization
        else:
            self.Pform = "heterMP"
            self.t_half_d = 100000  # As per The Full multi parameterizatio
        self.parentMP = parentMP
        self.parentSPM = parentSPM
        self.Pdensity_kg_m3 = parentMP.Pdensity_kg_m3 * (
            parentMP.Pvolume_m3 / (parentMP.Pvolume_m3 + parentSPM.Pvolume_m3)
        ) + parentSPM.Pdensity_kg_m3 * (
            parentSPM.Pvolume_m3 / (parentMP.Pvolume_m3 + parentSPM.Pvolume_m3)
        )
        self.radius_m = (
            3 * (parentMP.Pvolume_m3 + parentSPM.Pvolume_m3) / (4 * math.pi)
        ) ** (
            1 / 3
        )  # Note: this is an equivalent radius. MP-SPM most likely not truly spherical
        self.diameter_m = self.radius_m * 2
        self.diameter_um = self.diameter_m * 1e6
        self.Pshape = (
            parentMP.Pshape
        )  # to be updated for biofilm, could argue that shape is retained (unlike for SPM-bound)
        
        # add dimensions
        if parentMP.PdimensionY_um == 0:
            self.PdimensionY_um = 0
        else:
            self.PdimensionY_um = parentMP.PdimensionY_um + parentSPM.diameter_um

        if parentMP.PdimensionZ_um == 0:
            self.PdimensionZ_um = 0
        else:
            self.PdimensionZ_um = parentMP.PdimensionZ_um + parentSPM.diameter_um

        if parentMP.PdimensionX_um == 0:
            self.PdimensionX_um = 0
        else:
            self.PdimensionX_um = parentMP.PdimensionX_um + parentSPM.diameter_um
            
        self.PdimensionX_m = self.PdimensionX_um / 1000000  # shortest size
        self.PdimensionY_m = self.PdimensionY_um / 1000000  # longest size
        self.PdimensionZ_m = self.PdimensionZ_um / 1000000  # intermediate size

    # methods

    # volume calculation - currently simple version.
    # more complexity to be added later:
    # different formulas for different particle shapes.
    # currently defined for spheres, fibres, cylinders, pellets and irregular fragments
    def calc_volume_heter(self, parentMP, parentSPM):
        if self.Pshape == "sphere":
            self.Pvolume_m3 = parentMP.Pvolume_m3 + parentSPM.Pvolume_m3
            # calculates volume (in m3) of spherical particles from MP radius (x dimension)
            self.CSF = 1
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)

        elif (
            self.Pshape == "fibre"
            or self.Pshape == "fiber"
            or self.Pshape == "cylinder"
        ):
            self.Pvolume_m3 = parentMP.Pvolume_m3 + parentSPM.Pvolume_m3
            # calculates volume (in m3) of fibres or cylinders from diameter and
            # length assuming cylindrical shape where X is the shorterst size (radius) ans Y the longest (heigth)
            self.CSF = (self.radius_m) / math.sqrt(self.PdimensionY_m * self.radius_m)
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)

        elif self.Pshape == "pellet" or self.Pshape == "fragment":
            self.Pvolume_m3 = parentMP.Pvolume_m3 + parentSPM.Pvolume_m3
            # approximate volume calculation for irregular fragments
            # approximated as a cuboid using longest, intermediate and shortest length
            #!! Note: not sure if pellets fits best here or rather as sphere/cylinder
            # might adjust later!!
            self.CSF = self.PdimensionX_m / math.sqrt(
                self.PdimensionY_m * self.PdimensionZ_m
            )
            # calculate corey shape factor (CSF)
            # (Waldschlaeger 2019, doi:10.1021/acs.est.8b06794)

        else:
            print("Error: unknown shape")

        # print("Calculated " + self.Pname + " volume: " + str(self.Pvolume_m3) + " m3")
