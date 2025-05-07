# initialize_inputs.py

import numpy as np
import pandas as pd

def initialize_inputs(file_path, s_fc, soil_depth):
    precip = pd.read_csv(file_path, index_col=0)
    n_steps = len(precip)

    # Pre-allocate arrays
    interception = np.zeros(n_steps)
    infil = np.zeros(n_steps)
    perco = np.zeros(n_steps)
    sm = np.zeros(n_steps)
    E_stress_tc = np.zeros(n_steps)
    E_stress_sc = np.zeros(n_steps)
    E_stress_bs = np.zeros(n_steps)
    evap_actual_tc = np.zeros(n_steps)
    evap_actual_sc = np.zeros(n_steps)
    evap_actual_bs = np.zeros(n_steps)
    total_evap = np.zeros(n_steps)
    AWC = np.zeros(n_steps)
    run_off = np.zeros(n_steps)
    y = np.zeros(n_steps)

    # Initial conditions
    sm[0] = s_fc
    E_stress_tc[0] = 1
    E_stress_sc[0] = 1
    E_stress_bs[0] = 1
    infil[0] = 0.5 * s_fc
    perco[0] = 0
    total_evap[0] = 0
    evap_actual_tc[0] = 0.5
    evap_actual_sc[0] = 0.5
    evap_actual_bs[0] = 0.5
    interception[0] = 0.1
    AWC[0] = 0.5

    return {
        "precip": precip,
        "n_steps": n_steps,
        "interception": interception,
        "infil": infil,
        "perco": perco,
        "sm": sm,
        "E_stress_tc": E_stress_tc,
        "E_stress_sc": E_stress_sc,
        "E_stress_bs": E_stress_bs,
        "evap_actual_tc": evap_actual_tc,
        "evap_actual_sc": evap_actual_sc,
        "evap_actual_bs": evap_actual_bs,
        "total_evap": total_evap,
        "AWC": AWC,
        "run_off": run_off,
        "y": y
    }
