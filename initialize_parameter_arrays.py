import parameters
import numpy as np


#calculate lenght of simulation
def calculate_simulation_length(input_data):
    """Returns length of simulation based on length of input data
    Input: input dataframe of daily weather data
    Output: length of simulation
    """
    return len(input_data)

def initialized_arrays(input_data):
    # Calculate length of simulation
    len_sim = calculate_simulation_length(input_data)
    # initialize arrays as zeros arrays based on the length of input data
    for key in parameters.arrays.keys():
        parameters.arrays[key] = np.zeros(len_sim).reshape(len_sim, 1)

    return parameters.arrays
