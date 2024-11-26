##______________________________________________________________
#       MASK FOR ONLY NANOBEAMS
#_______________________________________________________________

# Let's import basic stuff
import samplemaker.layout as smlay # used for layout 
import samplemaker.makers as sm # used for drawing
import samplemaker.devices as smdev # used for device function
# Let's use numpy arrays
import numpy as np
import mydevices
from mydevices import *

# Create a simple mask layout
themask_only_nanobeams= smlay.Mask("NanoBeams_no_cavity")

geom = sm.GeomGroup()
nanobeam = smdev.Device.build_registered("NANOBEAM")
nanobeam.set_relevant_params(["taper_length", 
                              "top_length", 
                              "top_width", 
                              "top_tip_width", 
                              "theter_length", 
                              "theter_width", 
                              "theters_distance",  
                              "shrink_um"])


nanobeam.set_param("theter_length", 1)
nanobeam.set_param("top_length", 10)
nanobeam.set_param("taper_length", 8)
nanobeam.set_param("m", 0.85)
nanobeam.set_param("offset_from_tip", 1)
nanobeam.set_param("shrink_um", 0.012)
nanobeam.set_param("N_atoms_left", 0)
nanobeam.set_param("N_atoms_right", 0)


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
theter_width = np.arange(0.150, 0.35, 0.05).tolist() 
theters_distance = [4, 5, 8, 10]
one_side_theters = [0, 0, 1, 1]


# Generate param dictionary with arbitrary parameters
param_cols, N_cols= generate_param_dict(
    theters_distance = theters_distance,
    one_side_theters = one_side_theters,
)

param_rows, N_rows = generate_param_dict(
    theter_width = theter_width,

)

row_id = [i for i in range(N_rows)]
column_id = [i for i in range(N_cols)]

param_rows["row_id"] = row_id
param_cols["col_id"] = column_id


# Device Table with waveguide with 10 um length
nanobeam.set_param("top_length", 10)
nanobeam.set_param("group_id", 0)  
tab10um = smlay.DeviceTable(nanobeam,
                        N_rows, N_cols,
                        param_rows,
                        param_cols
                        )

# Device Table with devices with 20 um length
nanobeam.set_param("top_length", 20)
nanobeam.set_param("group_id", 1)   
tab20um = smlay.DeviceTable(nanobeam,
                        N_rows, N_cols,
                        param_rows,
                        param_cols
                        )


tab10um.auto_align(20, 20, numkey=5)
tab20um.auto_align(20, 20, numkey=5)


themask_only_nanobeams.addDeviceTable(tab10um, 0, 0)
themask_only_nanobeams.addDeviceTable(tab20um, 0, 400)



themask_only_nanobeams.exportGDS()

