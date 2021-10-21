import numpy as np

namelist_items = {
    "gridlist": [
        ("geometry", str),
        ("x_start", (int, np.integer, float)),
        ("x_end", (int, np.integer, float)),
        ("gridpoints", (int, np.integer)),
        ("mesh_accumulation", bool),
        ("ev_1", (int, np.integer, float)),
        ("ev_2", (int, np.integer, float)),
        ("sigma_1", (int, np.integer, float)),
        ("sigma_2", (int, np.integer, float)),
        ("force_r0", bool),
        ("coaxial", bool),
    ],
    "equilibriumlist": [
        ("equilibrium_type", str),
        ("boundary_type", str),
        ("use_defaults", bool),
        ("remove_spurious_eigenvalues", bool),
        ("nb_spurious_eigenvalues", (int, np.integer)),
    ],
    "savelist": [
        ("write_matrices", bool),
        ("write_eigenfunctions", bool),
        ("write_derived_eigenfunctions", bool),
        ("show_results", bool),
        ("basename_datfile", str),
        ("basename_logfile", str),
        ("output_folder", str),
        ("logging_level", (int, np.integer)),
        ("dry_run", bool),
        ("write_eigenfunction_subset", bool),
        ("eigenfunction_subset_center", complex),
        ("eigenfunction_subset_radius", (int, np.integer, float)),
    ],
    "physicslist": [
        ("mhd_gamma", float),
        ("flow", bool),
        ("radiative_cooling", bool),
        ("ncool", (int, np.integer)),
        ("cooling_curve", str),
        ("external_gravity", bool),
        ("thermal_conduction", bool),
        ("use_fixed_tc_para", bool),
        ("fixed_tc_para_value", (int, np.integer, float)),
        ("use_fixed_tc_perp", bool),
        ("fixed_tc_perp_value", (int, np.integer, float)),
        ("resistivity", bool),
        ("use_fixed_resistivity", bool),
        ("fixed_eta_value", (int, np.integer, float)),
        ("use_eta_dropoff", bool),
        ("dropoff_edge_dist", (int, np.integer, float)),
        ("dropoff_width", (int, np.integer, float)),
        ("viscosity", bool),
        ("viscosity_value", (int, np.integer, float)),
        ("hall_mhd", bool),
        ("hall_dropoff", bool),
        ("elec_pressure", bool),
        ("elec_inertia", bool),
        ("inertia_dropoff", bool),
        ("incompressible", bool),
    ],
    "unitslist": [
        ("cgs_units", bool),
        ("unit_density", (int, np.integer, float)),
        ("unit_temperature", (int, np.integer, float)),
        ("unit_magneticfield", (int, np.integer, float)),
        ("unit_length", (int, np.integer, float)),
        ("mean_molecular_weight", (int, np.integer, float)),
    ],
    "solvelist": [
        ("solver", str),
        ("arpack_mode", str),
        ("number_of_eigenvalues", (int, np.integer)),
        ("which_eigenvalues", str),
        ("maxiter", (int, np.integer)),
        ("sigma", (int, np.integer, float, complex)),
    ],
}
