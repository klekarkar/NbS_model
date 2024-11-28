import numpy as np
from compute_interception import daily_interception_threshold
from compute_ET import *
import parameters
from parameters import *
import initialize_parameter_arrays
from data_processing import *
import importlib
import os
from infiltration_models import *
#==============================================================================
"""Soil Moisture Balance Model

Initial arrays of soil moisture balance variables
"""
#==============================================================================
# Number of time steps
src=r"W:\VUB\_main_research\scripts\wetland_model\conceptual_model\data"
precip = precip=pd.read_csv(os.path.join(src,"Boechout_precip_ETo.csv"), index_col=0)
n_steps = len(precip)

# Pre-allocate arrays with zeros
interception = np.zeros(n_steps)  # Interception array
infil = np.zeros(n_steps)         # Infiltration array
perco = np.zeros(n_steps)         # Percolation array
sm = np.zeros(n_steps) 
E_stress_tc = np.zeros(n_steps) # Evaporation stress tall canopy
E_stress_sc = np.zeros(n_steps) # Evaporation stress short canopy
E_stress_bs = np.zeros(n_steps) # Evaporation stress
evap_actual_tc = np.zeros(n_steps)
evap_actual_sc = np.zeros(n_steps)
evap_actual_bs = np.zeros(n_steps)
total_evap = np.zeros(n_steps)
AWC = np.zeros(n_steps)
run_off = np.zeros(n_steps)
y = np.zeros(n_steps)

# Initial conditions
sm[0] = s_fc  # Initial soil moisture is at field capacity
E_stress_tc[0] = 1  # Initial evaporation stress is 1
E_stress_sc[0] = 1  
E_stress_bs[0] = 1  
infil[0] = 0.5 * s_fc  # max infiltration at the start
perco[0] = 0  # No percolation at the start
total_evap[0] = 0
evap_actual_tc[0] =  0.5 #initial evaporation
evap_actual_sc[0] =  0.5
evap_actual_bs[0] =  0.5
interception[0] = 0.1  #
AWC[0] = 0.5
soil_depth = 300  # mm

def simulate_soil_water_balance(
    time, sm, pET_k, precip, LAI, s_fc, s_wp, tau, soil_depth, Ks, S, dt, 
    frac_tall_canopy, frac_short_canopy, frac_bare_soil, alpha_tall_canopy, 
    alpha_short_canopy, interception, total_evap, evap_actual_tc, evap_actual_sc, 
    evap_actual_bs, AWC, infil, perco
):
    """
    Simulate the water balance for a single time step.
    """
    # Step 1: Calculate interception
    interception[time] = min(precip['precipitation'].iloc[time], 0.05 * LAI)

    # Step 2: Evaporate interception first
    if interception[time] >= pET_k.iloc[time]:
        total_evap[time] = pET_k.iloc[time] # If interception alone satisfies PET, total ET is limited to PET
        evap_actual_tc[time] = 0  # No evaporation from tall canopy
        evap_actual_sc[time] = 0  # No evaporation from short canopy
        evap_actual_bs[time] = 0  # No evaporation from bare soil
    else:
        remaining_PET = pET_k.iloc[time] - interception[time]

        # Step 3: Calculate evaporation stress for different components
        E_stress_tc = max(0, 1 - ((s_fc - sm[time - 1]) / (s_fc - s_wp))**2)
        E_stress_sc = max(0, 0.5 * (1 - np.sqrt((s_fc - sm[time - 1]) / (s_fc - s_wp)) + tau / 0.8))
        E_stress_bs = max(0, 1 - np.sqrt((s_fc - sm[time - 1]) / (s_fc - s_wp)))

        # Calculate potential evaporation for each component
        pot_evap_tc = remaining_PET * frac_tall_canopy * alpha_tall_canopy
        pot_evap_sc = remaining_PET * frac_short_canopy * alpha_short_canopy
        pot_evap_bs = remaining_PET * frac_bare_soil * alpha_short_canopy

        # Calculate actual evaporation for each component
        evap_actual_tc[time] = pot_evap_tc * E_stress_tc
        evap_actual_sc[time] = pot_evap_sc * E_stress_sc
        evap_actual_bs[time] = pot_evap_bs * E_stress_bs

        # Step 4: Calculate total ET
        total_evap[time] = (
            interception[time] + evap_actual_tc[time] + evap_actual_sc[time] + evap_actual_bs[time]
        )

    # Step 5: Calculate infiltration
    # Calculate available water capacity
    AWC[time] = max(0, (s_fc - sm[time - 1]) * soil_depth)

    # Create an instance of the infiltration model
    infiltration_model = InfiltrationModel(Ks, S)

    # Calculate cumulative daily infiltration in mm
    daily_infiltration = infiltration_model.valiantzas_model(dt) * 24  # Convert hourly to daily
    infil[time] = min(daily_infiltration, precip['precipitation'].iloc[time] - interception[time])

    # Step 6: Update soil moisture (before percolation)
    sm[time] = sm[time - 1] + (infil[time] - total_evap[time]) / soil_depth

    # Step 7: Calculate percolation if soil moisture exceeds field capacity
    if sm[time] >= s_fc:
        perco[time] = (sm[time] - s_fc) * soil_depth  # Excess water becomes percolation
        sm[time] = s_fc  # Soil moisture is capped at field capacity
    else:
        perco[time] = 0  # No percolation if below field capacity

    # Step 8: Enforce lower boundary for soil moisture
    sm[time] = max(s_wp, sm[time])  # Prevent negative soil moisture

    return interception, total_evap, evap_actual_tc, evap_actual_sc, evap_actual_bs, AWC, infil, perco, sm
