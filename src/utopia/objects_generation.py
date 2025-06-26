""" Contains all the functions to generate the model objects: model box, model compartments and the model particles (MP ans SPM)"""

import math

from helpers.globalConstants import *


class Box:
    """Class box generates one object box representing the unit world in the case of the UTOPIA parameterization that can contain 17 compartments and that can have conexions to other boxes (for example if conecting several UTOPIA boxes to give spatial resolution)"""

    # class attributes
    description = "Generic Box class"

    def __init__(
        self,
        Bname,
        Bdepth_m=None,
        Blength_m=None,
        Bwidth_m=None,
        Bvolume_m3=None,
        Bconexions=None,
    ):
        # Assign attributes to self (instance attributes). Those set up as None are optional attributes
        self.Bname = Bname  # name of the box
        self.Bdepth_m = Bdepth_m  # depth of the box
        self.Blength_m = Blength_m  # length of the box
        self.Bwidth_m = Bwidth_m  # width of the box
        self.Bvolume_m3 = Bvolume_m3  # volume of the box
        self.compartments = []  # list of compartments in the box
        self.Bconexions = Bconexions  # conexions to other model boxes

    def __repr__(self):
        return (
            "{"
            + self.Bname
            + ", "
            + str(self.Bdepth_m)
            + ", "
            + str(self.Blength_m)
            + ", "
            + str(self.Bwidth_m)
            + "}"
        )

    def add_compartment(self, comp):
        self.compartments.append(comp)
        comp.assign_box(self)

    def calc_Bvolume_m3(self):
        if self.Bvolume_m3 is None:
            if any(
                attr is None for attr in [self.Bdepth_m, self.Blength_m, self.Bwidth_m]
            ):
                print(
                    "Missing parameters needded to calculate Box volume --> calculating based on compartments volume"
                )
                if len(self.compartments) == 0:
                    print(
                        "No compartments assigned to this model box --> use add_compartment(comp)"
                    )
                else:
                    vol = []
                    for c in range(len(self.compartments)):
                        if self.compartments[c].Cvolume_m3 is None:
                            print(
                                "Volume of compartment "
                                + self.compartments[c].Cname
                                + " is missing"
                            )
                            continue
                        else:
                            vol.append(self.compartments[c].Cvolume_m3)
                    self.Bvolume_m3 = sum(vol)
            else:
                self.Bvolume_m3 = self.Bdepth_m * self.Blength_m * self.Bwidth_m
                # print("Box volume: " + str(self.Bvolume_m3)+" m3")
        else:
            print("Box volume already assigned: " + str(self.Bvolume_m3) + " m3")


class Compartment:
    """Class Compartment (parent class) generates compartment objects that belong by default to an assigned model box (Cbox). Each compartment contains four different particle objects corresponding to the 4 described aggregation states of UTOPIA (freeMP, heterMP, biofMP, heterBiofMP) and the processes that can occur in the compartment are listed under the processess attribute. Each compartment has a set of connexions withing the UTOPIA box listed in the conexions attribute wich will be asigned by reading on the conexions input file of the model."""

    def __init__(
        self,
        Cname,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
    ):
        self.Cname = Cname
        self.Cdepth_m = Cdepth_m
        self.Clength_m = Clength_m
        self.Cwidth_m = Cwidth_m
        self.Cvolume_m3 = Cvolume_m3
        self.CsurfaceArea_m2 = CsurfaceArea_m2
        self.particles = {
            "freeMP": [],
            "heterMP": [],
            "biofMP": [],
            "heterBiofMP": [],
        }  # dictionary of particles in the compartment
        self.processess = [
            "degradation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
        ]
        self.connexions = []

    def assign_box(self, Box):
        self.CBox = Box

    def add_particles(self, particle):
        self.particles[particle.Pform].append(particle)
        particle.assign_compartment(self)

    def calc_volume(self):
        if self.Cvolume_m3 is None:
            if any(
                attr is None for attr in [self.Cdepth_m, self.Clength_m, self.Cwidth_m]
            ):
                print(
                    "Missing parameters needded to calculate compartment volume --> Try calc_vol_fromBox or add missing values to compartment dimensions"
                )

            else:
                self.Cvolume_m3 = self.Cdepth_m * self.Clength_m * self.Cwidth_m
                # print(
                #     "Calculated "
                #     + self.Cname
                #     + " volume: "
                #     + str(self.Cvolume_m3)
                #     + " m3"
                # )
        else:
            pass
            # print("Assigned " + self.Cname + " volume: " + str(self.Cvolume_m3) + " m3")

    def calc_vol_fromBox(self):
        self.Cvolume_m3 = (
            self.CBox.Bvolume_m3 * self.CBox.CvolFractionBox[self.Cname.lower()]
        )

    def calc_particleConcentration_Nm3_initial(self):
        for p in self.particles:
            for s in self.particles[p]:
                self.particles[p][s].initial_conc_Nm3 = (
                    self.particles[p][s].Pnumber / self.Cvolume_m3
                )


"""Compartment Subclasses (inheritances) of the class compartment add extra attributes to the compatment that define the type of compartment """


class compartment_water(Compartment):
    def __init__(
        self,
        Cname,
        SPM_mgL,
        waterFlow_m3_s,
        T_K,
        G,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
        flowVelocity_m_s=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.SPM_mgL = SPM_mgL
        self.flowVelocity_m_s = flowVelocity_m_s
        self.waterFlow_m3_s = waterFlow_m3_s
        self.T_K = T_K
        self.G = G  # Shear rate (G, in s−1)
        self.processess = [
            "discorporation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
            "mixing",
        ]


class compartment_surfaceSea_water(Compartment):
    # added new processess to the list of processess. new attributes that migth be needed to this processess should be added here
    def __init__(
        self,
        Cname,
        SPM_mgL,
        waterFlow_m3_s,
        T_K,
        G,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
        flowVelocity_m_s=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.SPM_mgL = SPM_mgL
        self.flowVelocity_m_s = flowVelocity_m_s
        self.waterFlow_m3_s = waterFlow_m3_s
        self.T_K = T_K
        self.G = G  # Shear rate (G, in s−1)
        self.processess = [
            "discorporation",
            "fragmentation",
            "heteroaggregation",
            "heteroaggregate_breackup",
            "biofouling",
            "defouling",
            "advective_transport",
            "settling",
            "rising",
            "mixing",
            "sea_spray_aerosol",
            "beaching",
        ]


class compartment_sediment(Compartment):
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
            "sediment_resuspension",
            "burial",
        ]


class compartment_soil_surface(Compartment):
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
            "runoff_transport",
            "percolation",
            "soil_air_resuspension",
            "soil_convection",
        ]
        # self.earthworm_density_in_m3 = earthworm_density_in_m3
        # self.Qrunoff_m3 = Qrunoff_m3


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
            "soil_convection",
        ]

    # retention_in_soil (straining?) of the particles in soil following heteroaggregation with geocolloids?
    # shall we also include heteroaggregation/heteroaggegrate break-up processess in the soil compartment? In SimpleBox for Nano they do account for aggregation and attachment

    # Difference between retention in soil and sequestration deep soil: sequestrations deep soil is like burial in deep sediments (elemination process-->out of the system)


class compartment_air(Compartment):
    def __init__(
        self,
        Cname,
        T_K=None,
        wind_speed_m_s=None,
        I_rainfall_mm=None,
        Cdepth_m=None,
        Clength_m=None,
        Cwidth_m=None,
        Cvolume_m3=None,
        CsurfaceArea_m2=None,
        flowVelocity_m_s=None,
    ):
        super().__init__(
            Cname, Cdepth_m, Clength_m, Cwidth_m, Cvolume_m3, CsurfaceArea_m2
        )
        self.T_K = T_K
        self.wind_speed_m_s = wind_speed_m_s
        self.I_rainfall_mm = I_rainfall_mm
        self.flowVelocity_m_s = flowVelocity_m_s
        self.processess = [
            "discorporation",
            "fragmentation",
            "wind_trasport",
            "dry_deposition",
            "wet_deposition",
        ]
        # shall we also include heteroaggregation/heteroaggegrate break-up processess in the air compartment?


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
        self.radius_m = (
            self.PdimensionX_um / 1e6
        )  # In spherical particles from MP radius (x dimension)
        self.diameter_m = self.radius_m * 2
        self.diameter_um = self.diameter_m * 1e6
        self.Pemiss_t_y = 0  # set as 0 until the inputs_emissions.csv file is used
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
            self.Pvolume_m3 = math.pi * (self.radius_m) ** 2 * (self.PdimensionY_m)
            # calculates volume (in m3) of fibres or cylinders from diameter and
            # length assuming cylindrical shape where X is the shorterst size (radius) ans Y the longest (heigth)
            self.CSF = (self.radius_m) / math.sqrt(self.PdimensionY_m * self.radius_m)
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


def generate_objects():
    """Function for generating the UTOPIA model objects: model box, model compartments and the model particles"""


###CONTINUE HERE###
