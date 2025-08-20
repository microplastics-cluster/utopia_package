'''
Author: Xiaoyu Zhang (xiaoyu.zhang@aces.su.se)
Date: 2025-07-23 13:42:13
LastEditTime: 2025-07-23 16:13:28
Description: this script is for testing the settling velocity function
             and comparing the settling velocity, Reynolds number and drag coefficient with measured data.
             (test data extracted from the paper: https://doi.org/10.1016/j.envres.2023.115783)
'''
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np

def calculate_settling_velocity(shape, d_p, rho_p, rho_f, mu, g=9.81):
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
    Re    : Reynolds number
    Cd    : Drag coefficient
    """
    # Shape normalization
    shape = shape.lower()
    
    if 'sphere' in shape:
        # Stokes' Law (Re < 0.1, viscous-dominated, d_p < ~100 µm)
        v_s_stokes = g * (1 / 18) * ((rho_p - rho_f) / mu) * (d_p**2)
        Re_stokes = (rho_f * v_s_stokes * d_p) / mu  # Calculate Reynolds number
        Cd_stokes = 24 / Re_stokes  # Drag coefficient for Stokes regime

        if Re_stokes < 0.1:
            return v_s_stokes, Re_stokes, Cd_stokes

        # Intermediate regime (0.1 < Re < 1000)
        # Iterative approach to solve for velocity since Cd depends on Re
        v_s = v_s_stokes  # Initial guess
        Re = Re_stokes
        for _ in range(10):  # Iterate for better accuracy
            Re = (rho_f * v_s * d_p) / mu
            Cd = max((24 / Re) * (1 + 0.15 * Re**0.687), 0.44)  # Empirical drag coefficient
            v_s = math.sqrt((4 * g * d_p * (rho_p - rho_f)) / (3 * Cd * rho_f))

        if 0.1 <= Re < 1000:
            return v_s, Re, Cd

        # Newton's Law (Re > 1000, inertial-dominated)
        Cd_newton = 0.44  # Approximate constant drag coefficient
        v_s_newton = math.sqrt((4 * g * d_p * (rho_p - rho_f)) / (3 * Cd_newton * rho_f))
        Re_newton = (rho_f * v_s_newton * d_p) / mu

        return v_s_newton, Re_newton, Cd_newton

    elif "fiber" in shape or "fibre" in shape or "cylinder" in shape:

        # Stokes' Law (Re < 0.1, viscous-dominated, d_p < ~100 µm)
        v_s_stokes = g * (1 / 18) * ((rho_p - rho_f) / mu) * (d_p**2)
        Re_stokes = (rho_f * v_s_stokes * d_p) / mu  # Calculate Reynolds number
        Cd_stokes = 24 / Re_stokes  # Drag coefficient for Stokes regime

        if Re_stokes < 1:
            return v_s_stokes, Re_stokes, Cd_stokes

        # Intermediate regime (1 < Re < 1000)
        # Iterative approach to solve for velocity since Cd depends on Re
        v_s = v_s_stokes  # Initial guess
        Re = Re_stokes
        for _ in range(10):  # Iterate for better accuracy
            Re = (rho_f * v_s * d_p) / mu
            Cd = max(19 * (Re ** (-0.6)), 0.86)  # Empirical drag coefficient for fibers
            v_s = math.sqrt((4 * g * d_p * (rho_p - rho_f)) / (3 * Cd * rho_f))

        if 1 <= Re < 1000:
            return v_s, Re, Cd

        # Newton's Law (Re > 1000, inertial-dominated)
        Cd_newton = 0.86  # Approximate constant drag coefficient for fibers
        v_s_newton = math.sqrt((4 * g * d_p * (rho_p - rho_f)) / (3 * Cd_newton * rho_f))
        Re_newton = (rho_f * v_s_newton * d_p) / mu

        return v_s_newton, Re_newton, Cd_newton


# Read CSV data
df = pd.read_csv('settlingvdata.csv')

# Convert units and simulations
results = []
for i, row in df.iterrows():
    # Convert units (mm to m)
    d_p = row['dp_mm'] * 0.001  # diameter in meters
    
    # parameter extraction
    shape = row['Shape']
    rho_p = row['rho_p_kgm3']
    rho_f = row['rho_w_kgm3']
    mu_dynamic = 9.820e-4  # kg/(m·s)

    # calculate settling velocity, Reynolds number and drag coefficient
    try:
        ws_calc, Re_calc, Cd_calc = calculate_settling_velocity(shape, d_p, rho_p, rho_f, mu_dynamic)
        ws_meas = row['ws_m_s']
        Rep_meas = row['Rep']
        CD_meas = row['CD']
        
        # Calculate errors
        ws_error = abs(ws_calc - ws_meas) / ws_meas * 100
        Re_error = abs(Re_calc - Rep_meas) / Rep_meas * 100
        Cd_error = abs(Cd_calc - CD_meas) / CD_meas * 100
        
        results.append({
            'No': row['No'],
            'Shape': shape,
            'd_p (mm)': row['dp_mm'],
            
            # Settling velocity
            'Measured ws (m/s)': ws_meas,
            'Calculated ws (m/s)': ws_calc,
            'ws Error (%)': ws_error,
            
            # Reynolds number
            'Measured Rep': Rep_meas,
            'Calculated Re': Re_calc,
            'Re Error (%)': Re_error,
            
            # Drag coefficient
            'Measured CD': CD_meas,
            'Calculated Cd': Cd_calc,
            'Cd Error (%)': Cd_error
        })
    except Exception as e:
        print(f"Error calculating row {i}: {e}")
        results.append({
            'No': row['No'],
            'Shape': shape,
            'd_p (mm)': row['dp_mm'],
            'Measured ws (m/s)': row['ws_m_s'],
            'Calculated ws (m/s)': None,
            'ws Error (%)': None,
            'Measured Rep': row['Rep'],
            'Calculated Re': None,
            'Re Error (%)': None,
            'Measured CD': row['CD'],
            'Calculated Cd': None,
            'Cd Error (%)': None
        })

# Create results dataframe
results_df = pd.DataFrame(results)

# Print results with all parameters
print(results_df[['No', 'Shape', 'd_p (mm)', 
                 'Measured ws (m/s)', 'Calculated ws (m/s)', 'ws Error (%)',
                 'Measured Rep', 'Calculated Re', 'Re Error (%)',
                 'Measured CD', 'Calculated Cd', 'Cd Error (%)']])

# Visualization - 3x1 subplots for ws, Re and Cd
plt.figure(figsize=(18, 6))

# Define colors for different shapes
shape_colors = {
    'Sphere': '#1f78b4',
    'Fiber': "#5a9755",
    'Circular_cylinder': "#fdb454"
}

# Settling velocity comparison
plt.subplot(1, 3, 1)
for shape in results_df['Shape'].unique():
    df_shape = results_df[results_df['Shape'] == shape]
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Measured ws (m/s)'],
        edgecolors=shape_colors.get(shape), 
        facecolors='none',  
        label=f'{shape} (Measured)',
        marker='o',
    )
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Calculated ws (m/s)'],
        color=shape_colors.get(shape),
        label=f'{shape} (Calculated)',
        marker='x'
    )

plt.xlabel('Equivalent Diameter (mm)')
plt.ylabel('Settling Velocity (m/s)')
plt.title('Measured vs Calculated Settling Velocity')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.yscale('log')
plt.xscale('log')

# Reynolds number comparison
plt.subplot(1, 3, 2)
for shape in results_df['Shape'].unique():
    df_shape = results_df[results_df['Shape'] == shape]
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Measured Rep'],
        edgecolors=shape_colors.get(shape), 
        facecolors='none',  
        label=f'{shape} (Measured)',
        marker='o',
    )
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Calculated Re'],
        color=shape_colors.get(shape),
        label=f'{shape} (Calculated)',
        marker='x'
    )

plt.xlabel('Equivalent Diameter (mm)')
plt.ylabel('Reynolds Number')
plt.title('Measured vs Calculated Reynolds Number')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.yscale('log')
plt.xscale('log')

# Drag coefficient comparison
plt.subplot(1, 3, 3)
for shape in results_df['Shape'].unique():
    df_shape = results_df[results_df['Shape'] == shape]
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Measured CD'],
        edgecolors=shape_colors.get(shape), 
        facecolors='none',  
        label=f'{shape} (Measured)',
        marker='o',
    )
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Calculated Cd'],
        color=shape_colors.get(shape),
        label=f'{shape} (Calculated)',
        marker='x'
    )

plt.xlabel('Equivalent Diameter (mm)')
plt.ylabel('Drag Coefficient (Cd)')
plt.title('Measured vs Calculated Drag Coefficient')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.yscale('log')
plt.xscale('log')

plt.tight_layout()
plt.savefig('settling_velocity_re_cd_comparison.png', dpi=300)
plt.show()


# Error analysis - 3x1 subplots for errors
plt.figure(figsize=(18, 6))

# Settling velocity error
plt.subplot(1, 3, 1)
for shape in results_df['Shape'].unique():
    df_shape = results_df[results_df['Shape'] == shape].dropna()
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['ws Error (%)'],
        color=shape_colors.get(shape),
        label=shape,
        s=80
    )

plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('Equivalent Diameter (mm)')
plt.ylabel('Settling Velocity Error (%)')
plt.title('Settling Velocity Error by Particle Shape')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.ylim(-50, 100)

# Reynolds number error
plt.subplot(1, 3, 2)
for shape in results_df['Shape'].unique():
    df_shape = results_df[results_df['Shape'] == shape].dropna()
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Re Error (%)'],
        color=shape_colors.get(shape),
        label=shape,
        s=80
    )

plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('Equivalent Diameter (mm)')
plt.ylabel('Reynolds Number Error (%)')
plt.title('Reynolds Number Error by Particle Shape')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.ylim(-50, 100)

# Drag coefficient error
plt.subplot(1, 3, 3)
for shape in results_df['Shape'].unique():
    df_shape = results_df[results_df['Shape'] == shape].dropna()
    plt.scatter(
        df_shape['d_p (mm)'], 
        df_shape['Cd Error (%)'],
        color=shape_colors.get(shape),
        label=shape,
        s=80
    )

plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.xlabel('Equivalent Diameter (mm)')
plt.ylabel('Drag Coefficient Error (%)')
plt.title('Drag Coefficient Error by Particle Shape')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.ylim(-50, 100)

plt.tight_layout()
plt.savefig('settling_velocity_re_cd_errors.png', dpi=300)
plt.show()

# Save detailed results to CSV file
results_df.to_csv('settling_velocity_results_comparison.csv', index=False)
print("Detailed results have been saved to csv file")