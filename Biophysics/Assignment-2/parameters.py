"""
`parameters.py` is used to implement al parameters.
"""
# ------ General parameters -----
def get_parameters():
    
    sigma = 10.0
    beta = 8/3
    #rho will be varied.
        
    # --- Time settingsv---
    t_min, t_max = 0.0, 50.0
    dt = 0.01
        
    # --- Initial conditions ---
    x0, y0, z0 = 2, 1, 1
    initial_condition = (x0, y0, z0)
    
    initial_conditions = [(3, 1, 0),
                          (-2, 2, 1),
                          (1, 1, 5),]
    return sigma, beta, t_min, t_max, dt, initial_condition, initial_conditions
