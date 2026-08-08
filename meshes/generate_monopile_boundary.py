import math


def generate_dealii_config(V_avg_ref, H_s, output_filename="bc_config.txt",
                           amp_ratio=0.2, half_model=True,
                           half_symmetry_ids=(13, 14)):
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
    amp_ratio : float
        Cyclic amplitude as a fraction of the mean load. The loads computed
        from V_avg_ref and H_s are the MEAN of each cycle; the triangular wave
        then spans mean*(1 -/+ amp_ratio). amp_ratio=0.2 gives a peak of
        1.2*mean, a trough of 0.8*mean, and hence a load ratio
        R = min/max = (1-amp_ratio)/(1+amp_ratio) = 0.667.
    half_model : bool
        True when the mesh is the y>0 half produced by revolving the profile
        through Pi (the current monopile_gmsh.geo). Adds a u_y=0 roller on the
        exposed y=0 faces. Set False only for a full 2*Pi mesh, where those
        faces do not exist and the extra ids would silently match nothing.
    half_symmetry_ids : tuple of int
        Boundary ids of the two exposed y=0 faces. These come from the Abaqus
        round-trip (*Surface, name=Surf-N -> id N), NOT from the .geo, which
        defines no Physical entities. Update these to whatever the export
        actually assigns -- an id that does not exist in the mesh applies no
        constraint and the model will drift out of plane.
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

    def calc_tower_traction_mean(z_top):
        """Returns PEAK traction (Pa) for windward tower segment."""
        V_z = get_wind_speed(z_top)
        # Normalized by semi-circumference (pi*D/2)
        traction = (rho_air * C_D_tower / math.pi) * (V_z**2)
        return traction

    def calc_wave_traction_mean():
        """Returns PEAK traction (Pa) for windward pile."""
        u_surf = (math.pi * H_s) / T_wave
        u_dot_surf = (2 * math.pi**2 * H_s) / (T_wave**2)

        drag_term = 0.5 * rho_water * C_D_pile * D_pile * (u_surf**2)
        inertia_term = rho_water * C_M_pile * (math.pi * D_pile**2 / 4.0) * u_dot_surf
        f_surf = drag_term + inertia_term

        # Normalized by semi-circumference
        t_wave = (2.0 * f_surf) / (math.pi * D_pile)
        return t_wave

    def calc_rna_mean():
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

    # --- Surf-13/14: y=0 symmetry plane of the HALF model (roller) ---
    # monopile_gmsh.geo revolves the profile by Pi, not 2*Pi, so only the
    # y>0 half is meshed and the y=0 plane is left as two flat exposed faces
    # (sweep start on the +X side, sweep end on the -X side).
    #
    # These get a SYMMETRY condition: u_y = 0, u_x and u_z FREE. One line per
    # face, component 1 only. AbstractField::setup_dirichlet_boundary_condition
    # passes fields.component_masks[..._1] to interpolate_boundary_values, so
    # exactly one component is constrained.
    #
    # Emitting all three components here would be an ENCASTRE, not a roller:
    # it would clamp the crack faces shut and suppress the very opening the
    # phase field is meant to resolve.
    #
    # No load rescaling is needed anywhere else in this file. The tractions
    # below are per-unit-AREA, and the half model's surfaces have half the
    # area, so each resultant force halves automatically -- which is correct
    # for a structure carrying half the domain.
    if half_model:
        for surf_id in half_symmetry_ids:
            lines.append(f"{surf_id} dirichlet 1 0.0")

    # --- Surf-2: RNA Loads (Top) ---
    thrust_mean_val, weight_const = calc_rna_mean()

    # 1. Constant Weight (Neumann)
    # Direction: (0, 0, -1). Magnitude: weight_const.
    # Vector applied: (0, 0, -weight_const)
    lines.append(f"2 neumann 0.0 0.0 {-weight_const:.6e}")

    # 2. Cyclic Thrust (Triangular Neumann)
    # The value from V_avg_ref/H_s is the MEAN of the cycle; the amplitude is
    # a prescribed fraction of it, so the wave spans mean*(1 -/+ amp_ratio).
    # Direction: (1, 0, 0)
    thrust_mean = thrust_mean_val
    thrust_amp = thrust_mean_val * amp_ratio
    lines.append(f"2 triangularneumann 1 0 0 {freq} {thrust_mean:.6e} {thrust_amp:.6e}")

    # --- Surf-3: Wave Load (Pile) ---
    wave_mean = calc_wave_traction_mean()
    wave_amp = wave_mean * amp_ratio
    # Direction: (1, 0, 0)
    lines.append(f"3 triangularneumann 1 0 0 {freq} {wave_mean:.6e} {wave_amp:.6e}")

    # --- Surf-4 to Surf-12: Wind Load Segments ---
    # Tower levels: 0, 10, 20, ..., 80, 87
    z_levels = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 87.0]

    # Surfaces 4 through 12 correspond to segments 0-1, 1-2, etc.
    for i in range(len(z_levels) - 1):
        surf_id = 4 + i
        z_top = z_levels[i + 1]  # Use top z for conservative calc

        wind_mean = calc_tower_traction_mean(z_top)
        wind_amp = wind_mean * amp_ratio

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
    # V_avg = 8.8 m/s, Hs = 2.0 m
    generate_dealii_config(
        V_avg_ref=8.8, H_s=2.0, output_filename="monopile_mid_boundary.txt"
    )
    # Uncomment to generate others:
    # generate_dealii_config(V_avg_ref=4.4, H_s=5.0, output_filename="config_low.txt")
    # generate_dealii_config(V_avg_ref=17.6, H_s=8.0, output_filename="config_overload.txt")
