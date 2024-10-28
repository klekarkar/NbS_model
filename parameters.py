"""Model parameters and initial conditions for the water balance model"""
import numpy as np

# Model parameters
A1 = 38755
Kh_GWlocal = 200
Kh_GWout = 1800
Kh_GWreg = 1
L_GWlocal = 3300
L_GWout = 3000
L_GWreg = 3e5
M = 0.9
MeanRootDepth = 2000
kevap_BS = 0.8
kv_WL = 1.5
kveg = 1.1
kvegsat = 0.9
LAI = 4.5
n = 0.20
phi_GW = n
phi_GWlocal = n
phi_GWout = n
ref_elev = 13.42 * 1000 #elevation of the catchment in mm
ref_elev_WL = 10.9 * 1000 #elevation of the wetland in mm
r_P = 0.2
sfc = 0.50
s_initial = sfc * 0.5
s_molecular_suction = 0.10
s_sce = 0.10
s_scs = 0.50

# Interception parameters
a = 0.50
k_extinct = 0.463

# Initialize water balance variables and variable catchment properties empty arrays.
arrays = {
    "A_ratio": np.array([ ]),
    "A2": np.array([ ]),
    "CALCdeltaH_GW": np.array([ ]),
    "CALCdeltaH_WL": np.array([ ]),
    "day": np.array([ ]),
    "deltaH_groundwater": np.array([ ]),
    "deltaH_wetland": np.array([ ]),
    "evap_BS": np.array([ ]),
    "evap_baresoil": np.array([ ]),
    "evap_baresoil_cont_to_s": np.array([ ]),
    "evap_BS_mm": np.array([ ]),
    "ep_WL": np.array([ ]),
    "eveg_sat": np.array([ ]),
    "eveg_sat_cont_to_y": np.array([ ]),
    "eveg_us": np.array([ ]),
    "eveg_us_cont_to_s": np.array([ ]),
    "gradient_local": np.array([ ]),
    "gradient_outflow": np.array([ ]),
    "Inf_cont_to_s_cm": np.array([ ]),
    "Inf_cont_to_s_moist": np.array([ ]),
    "interception_threshold": np.array([ ]),
    "net_precipitation": np.array([ ]),
    "moisture_dep_factor": np.array([ ]),
    "plant_stress_factor": np.array([ ]),
    "Q_GWreg": np.array([ ]),
    "Qss_GW": np.array([ ]),
    "qGW_Local": np.array([ ]),
    "qGW_out": np.array([ ]),
    "qLOCAL_GW": np.array([ ]),
    "qLOCAL_WL": np.array([ ]),
    "Recharge": np.array([ ]),
    "recharge_cont_to_y": np.array([ ]),
    "Rlinear": np.array([ ]),
    "R_us": np.array([ ]),
    "s": np.array([ ]),
    "water_table_elevation": np.array([ ]),
    "wetland_wl_elev": np.array([ ]),
    "tR_local": np.array([ ]),
    "tR_out": np.array([ ]),
    "water_level_wetland": np.array([ ]),
    "water_table_depth": np.array([ ]),
    "WSC": np.array([ ]),
    "y": np.array([ ]),
    "y_rech": np.array([ ]),
    "y_WL": np.array([ ]),
    "net_p_input": np.array([ ]),
    "run_off": np.array([ ]),
    "interception": np.array([ ]),
    "s_max": np.array([ ]),
    "delta_s": np.array([ ]),
    "cumulative_infiltration": np.array([ ])}

# model/parameters.py

# Initialize water balance variables and variable catchment properties as np.array([ ])
A_ratio = np.array([])
A2 = np.array([ ])
CALCdeltaH_GW = np.array([ ])
CALCdeltaH_WL = np.array([ ])
day = np.array([ ])
deltaH_groundwater = np.array([ ])
deltaH_wetland = np.array([ ])
evap_BS = np.array([ ])
eb_US = np.array([ ])
eb_US_cont_to_s = np.array([ ])
evap_BS_mm = np.array([ ])
ep_WL = np.array([ ])
eveg_sat = np.array([ ])
eveg_sat_cont_to_y = np.array([ ])
eveg_us = np.array([ ])
eveg_us_cont_to_s = np.array([ ])
gradient_local = np.array([ ])
gradient_outflow = np.array([ ])
Inf_cont_to_s_cm = np.array([ ])
Inf_cont_to_s_moist = np.array([ ])
net_precipitation = np.array([ ])
moisture_dep_factor = np.array([ ])
plant_stress_factor = np.array([ ])
Q_GWreg = np.array([ ])
Qss_GW = np.array([ ])
qGW_Local = np.array([ ])
qGW_out = np.array([ ])
qLOCAL_GW = np.array([ ])
qLOCAL_WL = np.array([ ])
Recharge = np.array([ ])
recharge_cont_to_y = np.array([ ])
Rlinear = np.array([ ])
R_us = np.array([ ])
s = np.array([ ])
water_table_elevation = np.array([ ])
wetland_wl_elev = np.array([ ])
tR_local = np.array([ ])
tR_out = np.array([ ])
water_level_wetland = np.array([ ])
water_table_depth = np.array([ ])
WSC = np.array([ ])
y = np.array([ ])
y_rech = np.array([ ])
y_WL = np.array([ ])
net_p_input = np.array([ ])
run_off = np.array([ ])
interception = np.array([ ])
s_max = np.array([ ])
delta_s = np.array([ ])
interception_threshold = np.array([ ])
evap_baresoil = np.array([ ])
evap_baresoil_cont_to_s = np.array([ ])
cumulative_infiltration = np.array([ ])
