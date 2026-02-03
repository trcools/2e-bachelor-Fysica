"""
`parameters.py` is used to implement al parameters for the Rössler system.
"""
# ------ General parameters -----
def get_parameters():
    
    # Rössler system parameters
    a = 0.2
    b = 0.2
    # c will be varied as the bifurcation parameter (similar to rho in Lorenz)
        
    # --- Time settings ---
    t_min, t_max = 0.0, 200.0  # Longer time for Rössler to show attractor
    dt = 0.01
        
    # --- Initial conditions ---
    x0, y0, z0 = 1, 1, 1
    initial_condition = (x0, y0, z0)
    
    initial_conditions = [(1, 1, 1),
                          (2, 1, 0),
                          (-1, 2, 1),]
    return a, b, t_min, t_max, dt, initial_condition, initial_conditions
