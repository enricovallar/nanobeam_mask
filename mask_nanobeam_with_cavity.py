#________________________________________________________________
#       MASK FOR NANOBEMS WITH CAVITY
#________________________________________________________________
# Let's import basic stuff
import samplemaker.layout as smlay # used for layout 
import samplemaker.makers as sm # used for drawing
import samplemaker.devices as smdev # used for device function
# Let's use numpy arrays
import numpy as np
import mydevices
from mydevices import *


themask= smlay.Mask("Nanobeams_with_cavity")
geom = sm.GeomGroup()
nanobeam = smdev.Device.build_registered("NANOBEAM")
nanobeam.set_relevant_params(["taper_length", 
                              "top_length", 
                              "top_width", 
                              "top_tip_width", 
                              "theter_length", 
                              "theter_width", 
                              "theters_distance",  
                              "shrink_um", 
                              "N_atoms_left",
                              "N_atoms_right",
                              "atom_radius",
                              "gap",
                              "lattice_constant",])

# Set the parameters for the nanobeam
nanobeam.set_param("theter_length", 1)
nanobeam.set_param("top_length", 10)
nanobeam.set_param("taper_length", 8)
nanobeam.set_param("top_width", 1)
nanobeam.set_param("m", 0.85)
nanobeam.set_param("offset_from_tip", 1)
nanobeam.set_param("shrink_um", 0.012)

nanobeam.set_param("atom_radius", 0.080)
nanobeam.set_param("one_side_theters", 1)
nanobeam.set_param("theters_distance", 8)
nanobeam.set_param("N_atoms_left", 14)
nanobeam.set_param("theter_width", 0.2)




from itertools import product
group_index = 0
element_row_index = 0
def generate_param_dict(**kwargs):
    # Initialize the dictionary that will store lists for each parameter
    param_dict = {key: [] for key in kwargs}
    
    # Get the cartesian product of all parameter values
    param_combinations = product(*kwargs.values())
    
    # Populate the param_dict with the combinations
    for combination in param_combinations:
        for i, key in enumerate(kwargs):
            param_dict[key].append(combination[i])
    
    return param_dict, len(param_dict[key])


# Define the parameters for the device table
theters_distance = [5,10]
lattice_constant = [-0.040, -0.030, -0.020,-0.010, 0, 0.010, 0.020, 0.03, 0.040]
lattice_constant = [0.350 + x for x in lattice_constant]
print(lattice_constant)
gap = [-0.040, -0.030, -0.020,-0.010, 0, 0.010, 0.020, 0.03, 0.040]
gap = [0.880 + x for x in gap]
print(gap)
N_atoms_right = [14, 4]

# Generate param dictionary with arbitrary parameters
param_cols, N_cols= generate_param_dict(
    lattice_constant = lattice_constant,
    N_atoms_right = N_atoms_right,
)

param_rows, N_rows = generate_param_dict(
    gap = gap
)

row_id = [i for i in range(N_rows)]
column_id = [i for i in range(N_cols)]

param_rows["row_id"] = row_id
param_cols["col_id"] = column_id


# Device Table with waveguide with 30 um length
nanobeam.set_param("top_length", 30)
nanobeam.set_param("group_id", 0)   
print("start writing first table")
tab30um = smlay.DeviceTable(nanobeam,
                        N_rows, N_cols,
                        param_rows,
                        param_cols
                        )


# Device Table with devices with 20 um length
print("start writing second table")
nanobeam.set_param("top_length", 20)
nanobeam.set_param("group_id", 1)   
tab20um = smlay.DeviceTable(nanobeam,
                        N_rows, N_cols,
                        param_rows,
                        param_cols
                        )


tab30um.auto_align(20, 20, numkey=5)
tab20um.auto_align(20, 20, numkey=5)


themask.addDeviceTable(tab30um, 0, 0)
themask.addDeviceTable(tab20um, 0, 600)
themask.exportGDS()
