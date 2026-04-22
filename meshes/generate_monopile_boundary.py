import math


def generate_dealii_config(V_avg_ref, H_s, output_filename="bc_config.txt"):
    """
    Generates a boundary condition configuration file for deal.II.

    Parameters:
    -----------
    V_avg_ref : float
        Reference mean wind speed at 10m height (m/s).
    H_s : float
        Significant wave height (m).
    output_filename : str
        Name of the output text file.
    """

    # ==========================================
    # 1. PHYSICAL CONSTANTS & GEOMETRY
    # ==========================================
    # Geometry
    D_pile = 5.5
    t_pile = 0.060
    z_hub = 87.0
    D_top = 3.9
    D_base = 5.5
    t_tower = 0.060

    # Environment & Material
    rho_air = 1.225
    rho_water = 1025.0
    g = 9.81

    # Coefficients
    C_D_tower = 0.7
    C_D_pile = 0.7
    C_M_pile = 2.0
    alpha = 0.12
    z_ref = 10.0

    # RNA Properties
    M_RNA = 350000.0  # kg
    D_rotor = 126.0  # m
    R_rotor = D_rotor / 2.0
    A_rotor = math.pi * R_rotor**2
    V_rated = 11.4  # m/s

    # Loading Parameters
    freq = 1.0  # Hz
    T_wave = 10.0  # s (Fixed per previous instructions)

    # ==========================================
    # 2. HELPER CALCULATION FUNCTIONS
    # ==========================================
    def get_wind_speed(z):
        if z <= 0:
            return 0.0
        return V_avg_ref * (z / z_ref) ** alpha

    def get_thrust_coeff(V_hub):
        if V_hub <= V_rated:
            return min(7.0 / V_rated, 1.0)
        else:
            val = (7.0 * V_rated**2) / (V_hub**3)
            return max(val, 0.005)

    def calc_tower_traction_max(z_top):
        """Returns PEAK traction (Pa) for windward tower segment."""
        V_z = get_wind_speed(z_top)
        # Normalized by semi-circumference (pi*D/2)
        traction = (rho_air * C_D_tower / math.pi) * (V_z**2)
        return traction

    def calc_wave_traction_max():
        """Returns PEAK traction (Pa) for windward pile."""
        u_surf = (math.pi * H_s) / T_wave
        u_dot_surf = (2 * math.pi**2 * H_s) / (T_wave**2)

        drag_term = 0.5 * rho_water * C_D_pile * D_pile * (u_surf**2)
        inertia_term = rho_water * C_M_pile * (math.pi * D_pile**2 / 4.0) * u_dot_surf
        f_surf = drag_term + inertia_term

        # Normalized by semi-circumference
        t_wave = (2.0 * f_surf) / (math.pi * D_pile)
        return t_wave

    def calc_rna_max():
        """Returns PEAK Thrust (Pa) and CONSTANT Weight (Pa)."""
        V_hub = get_wind_speed(z_hub)
        C_T = get_thrust_coeff(V_hub)

        F_thrust = 0.5 * rho_air * A_rotor * C_T * (V_hub**2)

        r_ext = D_top / 2.0
        r_int = r_ext - t_tower
        A_top = math.pi * (r_ext**2 - r_int**2)

        tr_thrust = F_thrust / A_top
        tr_weight = (M_RNA * g) / A_top
        return tr_thrust, tr_weight

    # ==========================================
    # 3. GENERATE CONFIG ENTRIES
    # ==========================================
    lines = []

    # --- Surf-1: Bottom Pinned (Dirichlet) ---
    # Fix u_x (0), u_y (1), u_z (2) to 0.0
    lines.append("1 dirichlet 0 0.0")
    lines.append("1 dirichlet 1 0.0")
    lines.append("1 dirichlet 2 0.0")

    # --- Surf-2: RNA Loads (Top) ---
    thrust_max, weight_const = calc_rna_max()

    # 1. Constant Weight (Neumann)
    # Direction: (0, 0, -1). Magnitude: weight_const.
    # Vector applied: (0, 0, -weight_const)
    lines.append(f"2 neumann 0.0 0.0 {-weight_const:.6e}")

    # 2. Cyclic Thrust (Triangular Neumann)
    # R=0 means Mean = Max/2, Amp = Max/2
    # Direction: (1, 0, 0)
    thrust_mean = thrust_max / 2.0
    thrust_amp = thrust_max / 2.0
    lines.append(f"2 triangularneumann 1 0 0 {freq} {thrust_mean:.6e} {thrust_amp:.6e}")

    # --- Surf-3: Wave Load (Pile) ---
    wave_max = calc_wave_traction_max()
    wave_mean = wave_max / 2.0
    wave_amp = wave_max / 2.0
    # Direction: (1, 0, 0)
    lines.append(f"3 triangularneumann 1 0 0 {freq} {wave_mean:.6e} {wave_amp:.6e}")

    # --- Surf-4 to Surf-12: Wind Load Segments ---
    # Tower levels: 0, 10, 20, ..., 80, 87
    z_levels = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 87.0]

    # Surfaces 4 through 12 correspond to segments 0-1, 1-2, etc.
    for i in range(len(z_levels) - 1):
        surf_id = 4 + i
        z_top = z_levels[i + 1]  # Use top z for conservative calc

        wind_max = calc_tower_traction_max(z_top)
        wind_mean = wind_max / 2.0
        wind_amp = wind_max / 2.0

        lines.append(
            f"{surf_id} triangularneumann 1 0 0 {freq} {wind_mean:.6e} {wind_amp:.6e}"
        )

    # ==========================================
    # 4. WRITE TO FILE
    # ==========================================
    with open(output_filename, "w") as f:
        f.write("\n".join(lines))

    print(f"Generated '{output_filename}' for V_avg={V_avg_ref} m/s, Hs={H_s} m.")


# ==========================================
# MAIN EXECUTION
# ==========================================
if __name__ == "__main__":
    # Example: Rated Load Case
    # V_avg = 8.8 m/s, Hs = 5.0 m
    generate_dealii_config(
        V_avg_ref=8.8, H_s=5.0, output_filename="monopile_mid_boundary.txt"
    )
    generate_dealii_config(
        V_avg_ref=4.4, H_s=5.0, output_filename="monopile_low_boundary.txt"
    )
    generate_dealii_config(
        V_avg_ref=17.6, H_s=8.0, output_filename="monopile_high_boundary.txt"
    )

    # Uncomment to generate others:
    # generate_dealii_config(V_avg_ref=4.4, H_s=5.0, output_filename="config_low.txt")
    # generate_dealii_config(V_avg_ref=17.6, H_s=8.0, output_filename="config_overload.txt")
