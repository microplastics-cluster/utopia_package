import numpy as np

# ---- Constants dependants of Temperature, Pressure, altitud. Have to be adapted to specific regions ------------
muAir = 1.789e-5  # [kg−1.s−1], dynamic viscosity of the air (Heigth, Temp dependent)
rhoAir = 1.1438  # [kg.m-3], density of air (Heigth, Temp dependent)
g0 = 9.81  # [m.s-2] gravitational acceleration on earth (Heigth dependent)


# Constant for Cunningham correction, here from Jennings, S, 1988
A1 = 1.257
A2 = 0.400
A3 = 1.100
mfpAir = 6.635e-8  # [m] Mean free path on dry air, from Jennings, S, 1988


# Calculate the Reynolds Number (based on settling velocity of Stokes,
# may be optimised with self consistancy)
def ReynoldsNumberFromStokes(d, rho):
    vg = np.multiply(d, d)
    vg = vg * g0 / (18.0 * muAir)
    vg = np.multiply(vg, (rho - rhoAir))

    Rep = np.multiply(vg, rhoAir)
    Rep = np.multiply(Rep, d)
    Rep = np.multiply(Rep, (1.0 / muAir))

    # print("XRe= ", Rep, Rep / 3.0)

    # Set Rep to a specific value for test (e.g., 13.4) !!!!!!!!!!! to remoove after test
    # Rep_value = 13.35
    # Rep = np.full_like(d, Rep_value)

    return Rep


# Placeholder function to calculate Reynolds number
def ReynoldsNumberFromVg(d, rho, vg):
    # Placeholder implementation
    # This function should calculate the Reynolds number based on the settling velocity
    Rep = np.multiply(vg, rhoAir)
    Rep = np.multiply(Rep, d)
    Rep = np.multiply(Rep, (1.0 / muAir))

    return Rep


# Calculate the Drag coefficient number
def dragCoefficient(d, rho, Rep):
    # Rep = ReynoldsNumberFromStokes (d, rho)
    # Set Rep to a specific value for test (e.g., 13.4) !!!!!!!!!!! to remoove after test

    # Rep = 13.35
    Cd = np.zeros_like(Rep)

    # Cd[Rep <= 1] = 24.0 / Rep[Rep <= 1]
    # Cd[(1 < Rep) & (Rep <= 1000.0)] = (
    #     24.0 / Rep[(1 < Rep) & (Rep <= 750.0)]
    #     + 5.0 / np.power(Rep[(1 < Rep) & (Rep <= 750.0)], 0.6)
    #     + 0.44
    # )
    # Cd[(1000.0 < Rep) & (Rep <= 2.0e5)] = 0.44
    # Cd[Rep >= 2.0e5] = 0.10

    Rep_conditions = [
        Rep <= 1,
        (1 < Rep) & (Rep <= 1000.0),
        (1000.0 < Rep) & (Rep <= 2.0e5),
        Rep > 2.0e5,
    ]

    Cd_functions = [
        lambda x: 24.0 / x,
        lambda x: 24.0 / x + 5.0 / np.power(x, 0.6) + 0.44,
        0.44,
        0.10,
    ]

    Cd = np.piecewise(Rep, Rep_conditions, Cd_functions)
    
    # Cd = 24.0 * (1.0 + (0.15 * Rep**0.687)) / Rep
    # Cd = Cd + 0.42 / (1.0 + (42500.0/ (Rep**1.16)))
    return Cd


# Calculate rate constants of dry settling based on Newton regime
# (for big particles, generating a turbulent flow)
def kineticCstdrySettlingNewtonSphere(d, rho, Rep):
    d = np.array(d)
    rho = np.array(rho)

    # Calculate the drag coefficient Cd
    Cd = dragCoefficient(d, rho, Rep)
    # print("Cd sphere=", Cd)
    # Calculate the settling velocity
    v = np.sqrt(4.0 * d * g0 * (rho - rhoAir) / (3.0 * Cd * rhoAir))
    # v = np.sqrt(4.0 * d * g0 * (rho) / (3.0 * Cd * rhoAir))

    # elt = v / shape  # [m.s-1] / [m] The Boundary layer heigth
    elt = v
    return elt


# Placeholder function to calculate settling velocity
def calculate_settling(reynolds_number):
    # Placeholder implementation
    # This function should calculate a new estimate for the settling velocity based on the Reynolds number
    return new_settling_velocity


# Subroutine to calculate settling velocity iteratively
def get_settling(initial_Settling, d, rho, initial_Rep):

    settling_old = initial_Settling

    # Set convergence threshold
    tolerance = 0.001

    # Maximum number of iterations
    max_iterations = 20

    # Initialize iteration counter
    iteration = 0

    # Iterate until convergence or maximum iterations reached
    while iteration < max_iterations:
        # Calculate Reynolds number using current settling velocity
        reynolds = ReynoldsNumberFromVg(d, rho, settling_old)

        # Use Reynolds number to calculate new estimate for settling velocity
        settling_new = kineticCstdrySettlingNewtonSphere(d, rho, reynolds)

        # Check for convergence
        cvg = abs((settling_new - settling_old) / settling_new)
        # print("convergence <", tolerance, "?=", cvg)
        if np.all(cvg < tolerance):
            break

        # Update settling velocity for next iteration
        settling_old = settling_new

        # Increment iteration counter
        iteration += 1

    # Use final settling velocity for further calculations
    final_settling_velocity = settling_new
    reynolds = ReynoldsNumberFromVg(d, rho, final_settling_velocity)

    return final_settling_velocity
