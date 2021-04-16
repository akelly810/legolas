khi_setup = {
    "name": "kelvin_helmholtz",
    "config": {
        "gridpoints": 51,
        "parameters": {
            "k2": 0.0,
            "k3": 1.0,
            "cte_rho0": 1.0,
            "cte_p0": 10.0,
            "delta": 0.0,
            "g": 0.0,
            "alpha": 0.0,
            "theta": 0.0,
            "p1": 0.0,
            "p2": 0.0,
            "p3": 1.0,
            "p4": 0.0,
            "tau": 11.0,
        },
        "flow": True,
        "external_gravity": False,
        "equilibrium_type": "kelvin_helmholtz",
        "logging_level": 0,
        "show_results": False,
        "write_eigenfunctions": False,
        "write_matrices": False,
    },
    "image_limits": [
        {"xlims": (-300, 300), "ylims": (-0.2, 0.2)},
        {"xlims": (-30, 30), "ylims": (-0.2, 0.2)},
        {"xlims": (-1.2, 1.2), "ylims": (-0.2, 0.2)}
    ]
}
