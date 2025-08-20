import math


def calculate_settling_velocity(particle, d_p, rho_p, rho_f, mu, g=9.81):
    """
    Estimates the settling velocity of a particle in water based on its size.
    Automatically selects the correct equation depending on Reynolds number.

    Parameters:
    d_p   : Particle equivalent diameter (m)
    rho_p : Particle density (kg/m³)
    rho_f : Fluid density (kg/m³) (typically ~1000 kg/m³ for water)
    mu    : Dynamic viscosity of water (Pa·s or kg/(m·s)) (e.g., ~0.001 for water at 20°C)
    g     : Gravitational acceleration (m/s²) (default: 9.81 m/s²)

    Returns:
    v_s   : Settling velocity (m/s)
    regime: Flow regime used (Stokes, Intermediate, or Newton)
    """

    if particle.Pshape == "sphere":

        # Stokes' Law (Re < 0.1, viscous-dominated, d_p < ~100 µm)
        v_s_stokes = g * (1 / 18) * ((rho_p - rho_f) / mu) * (d_p**2)
        Re_stokes = (rho_f * v_s_stokes * d_p) / mu  # Calculate Reynolds number

        if Re_stokes < 0.1:
            return v_s_stokes  # , "Stokes' Law (laminar flow)"

    # Intermediate regime (0.1 < Re < 1000)
    # Iterative approach to solve for velocity since Cd depends on Re
    v_s = v_s_stokes  # Initial guess
    for _ in range(10):  # Iterate for better accuracy
        Re = (rho_f * v_s * d_p) / mu
        Cd = (24 / Re) * (
            1 + 0.15 * Re**0.687
        )  # Empirical drag coefficient (empirical Schiller–Naumann correlation)
        v_s = math.sqrt((4 * g * d_p * (rho_p - rho_f)) / (3 * Cd * rho_f))
        # Intermediate regime (0.1 < Re < 1000)
        # Iterative approach to solve for velocity since Cd depends on Re
        v_s = v_s_stokes  # Initial guess
        for _ in range(10):  # Iterate for better accuracy
            Re = (rho_f * v_s * d_p) / mu
            Cd = max(
                (24 / Re) * (1 + 0.15 * Re**0.687), 0.44
            )  # Empirical drag coefficient
            v_s = math.sqrt((4 * g * d_p * (rho_p - rho_f)) / (3 * Cd * rho_f))

        if 0.1 <= Re < 1000:
            return v_s  # , "Intermediate regime (empirical drag correction)"

        # Newton's Law (Re > 1000, inertial-dominated)
        Cd_newton = 0.44  # Approximate constant drag coefficient
        v_s_newton = math.sqrt(
            (4 * g * d_p * (rho_p - rho_f)) / (3 * Cd_newton * rho_f)
        )

        return v_s_newton  # , "Newton's Law (turbulent flow)"

    if particle.Pshape == "fiber" or "fibre" or "cylinder":

        # Stokes' Law (Re < 0.1, viscous-dominated, d_p < ~100 µm)
        v_s_stokes = g * (1 / 18) * ((rho_p - rho_f) / mu) * (d_p**2)
        Re_stokes = (rho_f * v_s_stokes * d_p) / mu  # Calculate Reynolds number

        if Re_stokes < 1:
            return v_s_stokes  # , "Stokes' Law (laminar flow)"

        # Intermediate regime (1 < Re < 1000)
        # Iterative approach to solve for velocity since Cd depends on Re
        v_s = v_s_stokes  # Initial guess
        for _ in range(10):  # Iterate for better accuracy
            Re = (rho_f * v_s * d_p) / mu
            Cd = max(
                19 * (Re ** (-0.6)), 0.86
            )  # Empirical drag coefficient for fibers (ref. DOI: 10.1016/j.envres.2023.115783)
            v_s = math.sqrt((4 * g * d_p * (rho_p - rho_f)) / (3 * Cd * rho_f))

        if 1 <= Re < 1000:
            return v_s  # , "Intermediate regime (empirical drag correction)"

        # Newton's Law (Re > 1000, inertial-dominated)
        Cd_newton = 0.86  # Approximate constant drag coefficient for fibers (ref. DOI: 10.1016/j.envres.2023.115783)
        v_s_newton = math.sqrt(
            (4 * g * d_p * (rho_p - rho_f)) / (3 * Cd_newton * rho_f)
        )

        return v_s_newton  # , "Newton's Law (turbulent flow)"

    else:
        raise ValueError("Unsupported particle shape: {}".format(particle.Pshape))


def calculate_rising_velocity(particle, d_p, rho_p, rho_f, mu, g=9.81):
    """
    Estimates the rising velocity of a particle in water based on its size.
    Automatically selects the correct equation depending on Reynolds number.

    Parameters:
    d_p   : Particle equivalent diameter (m)
    rho_p : Particle density (kg/m³)
    rho_f : Fluid density (kg/m³) (typically ~1000 kg/m³ for water)
    mu    : Dynamic viscosity of water (Pa·s or kg/(m·s)) (e.g., ~0.001 for water at 20°C)
    g     : Gravitational acceleration (m/s²) (default: 9.81 m/s²)

    Returns:
    v_r   : Rising velocity (m/s)
    regime: Flow regime used (Stokes, Intermediate, or Newton)
    """
    if particle.Pshape == "sphere":

        # Stokes' Law (Re < 0.1, viscous-dominated, d_p < ~100 µm)
        v_s_stokes = g * (1 / 18) * ((rho_p - rho_f) / mu) * (d_p**2)
        Re_stokes = (rho_f * v_s_stokes * d_p) / mu  # Calculate Reynolds number

        if Re_stokes < 0.1:
            return -v_s_stokes
            # "Stokes' Law (laminar flow)" # Negative for rising particles

        # Intermediate regime (0.1 < Re < 1000)
        # Iterative approach to solve for velocity since Cd depends on Re
        v_s = v_s_stokes  # Initial guess
        for _ in range(10):  # Iterate for better accuracy
            Re = (rho_f * v_s * d_p) / mu
            Cd = max(
                (24 / Re) * (1 + 0.15 * Re**0.687), 0.44
            )  # Empirical drag coefficient
            v_s = math.sqrt((4 * g * d_p * (rho_f - rho_p)) / (3 * Cd * rho_f))

        if 0.1 <= Re < 1000:
            return -v_s
        # ,"Intermediate regime (empirical drag correction)",
        # Negative for rising particles

        # Newton's Law (Re > 1000, inertial-dominated)
        Cd_newton = 0.44  # Approximate constant drag coefficient
        v_s_newton = math.sqrt(
            (4 * g * d_p * (rho_f - rho_p)) / (3 * Cd_newton * rho_f)
        )

        return -v_s_newton

    elif particle.Pshape == "fiber" or "fibre" or "cylinder":

        # Stokes' Law (Re < 0.1, viscous-dominated, d_p < ~100 µm)
        v_s_stokes = g * (1 / 18) * ((rho_p - rho_f) / mu) * (d_p**2)
        Re_stokes = (rho_f * v_s_stokes * d_p) / mu  # Calculate Reynolds number

        if Re_stokes < 1:
            return -v_s_stokes  # , "Stokes' Law (laminar flow)"

        # Intermediate regime (1 < Re < 1000)
        # Iterative approach to solve for velocity since Cd depends on Re
        v_s = v_s_stokes  # Initial guess
        for _ in range(10):  # Iterate for better accuracy
            Re = (rho_f * v_s * d_p) / mu
            Cd = max(
                19 * (Re ** (-0.6)), 0.86
            )  # Empirical drag coefficient for fibers (ref. DOI: 10.1016/j.envres.2023.115783)
            v_s = math.sqrt((4 * g * d_p * (rho_f - rho_p)) / (3 * Cd * rho_f))

        if 1 <= Re < 1000:
            return -v_s  # , "Intermediate regime (empirical drag correction)"

        # Newton's Law (Re > 1000, inertial-dominated)
        Cd_newton = 0.86  # Approximate constant drag coefficient for fibers (ref. DOI: 10.1016/j.envres.2023.115783)
        v_s_newton = math.sqrt(
            (4 * g * d_p * (rho_f - rho_p)) / (3 * Cd_newton * rho_f)
        )

        return (
            -v_s_newton
        )  # , "Newton's Law (turbulent flow)" # Negative for rising particles

    else:
        raise ValueError("Unsupported particle shape: {}".format(particle.Pshape))
