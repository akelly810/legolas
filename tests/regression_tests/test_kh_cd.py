kh_cd_setup = {
    "name": "kh_cd",
    "config": {
        # geometry is hard-coded for this equilibrium
        "gridpoints": 51,
        "parameters": {
            "k2": -1.0,
            "V": 1.63,
            "cte_rho0": 1.0,
            "cte_p0": 1.0,
            "Bz0": 0.25,
            "rc": 0.5,
            "rj": 1.0,
        },
        "flow": True,
        "equilibrium_type": "kelvin_helmholtz_cd",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-200, 200), "ylims": (-1, 1), "RMS_TOLERANCE": 2.1},
        {"xlims": (-20, 20), "ylims": (-1, 1)},
        {"xlims": (-4, 4), "ylims": (-1, 1)},
    ],
}
