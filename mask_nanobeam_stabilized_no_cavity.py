#________________________________________________________________
#       MASK FOR STABILIZED NANOBEMS 
#________________________________________________________________
# Let's import basic stuff
import samplemaker.layout as smlay # used for layout 
import samplemaker.makers as sm # used for drawing
import samplemaker.devices as smdev # used for device function
# Let's use numpy arrays
import numpy as np
import mydevices
from mydevices import *

themask_stabilized_nanobeams= smlay.Mask("Nanobeams_stabilized_no_cavity")
StabilizedNanobeam = smdev.Device.build_registered("FULLDEVICE")
StabilizedNanobeam.get_params()


from itertools import product
geom = sm.GeomGroup()
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

# Example usage:
theter_width_pad = np.arange(0.150, 0.35, 0.05).tolist()
theters_distance_pad = [4, 5, 8, 10, 20]
top_length = [10, 10, 20, 20]

StabilizedNanobeam.set_relevant_params(["taper_length", 
                              "top_length", 
                              "top_width", 
                              "top_tip_width", 
                              "theter_length", 
                              "theter_width", 
                              "theters_distance", 
                              "theter_width_pad",
                              "theter_length_pad",
                              "theters_distance_pad",
                              "aperture_size",
                              "holes_distance",
                              "pad_edge",
                              "shrink_um"])


StabilizedNanobeam.set_param("theter_width", 0.3)
StabilizedNanobeam.set_param("theter_length", 0.5)
StabilizedNanobeam.set_param("theter_length_pad", 1)
StabilizedNanobeam.set_param("theter_width_pad", 0.2)
StabilizedNanobeam.set_param("aperture_size", 0.5)
StabilizedNanobeam.set_param("taper_length", 8)
StabilizedNanobeam.set_param("m", 0.85)
StabilizedNanobeam.set_param("offset_from_tip", 1)
StabilizedNanobeam.set_param("shrink_um", 0.012)
StabilizedNanobeam.set_param("N_atoms_left", 0)
StabilizedNanobeam.set_param("N_atoms_right", 0)





# Generate param dictionary with arbitrary parameters
param_cols, N_cols= generate_param_dict(
    theter_width_pad = theter_width_pad,
    top_length = top_length
)

param_rows, N_rows = generate_param_dict(
    theters_distance_pad = theters_distance_pad
)

row_id = [i for i in range(N_rows)]
column_id = [i for i in range(N_cols)]


param_rows["row_id"] = row_id
param_cols["col_id"] = column_id

print("start writing")
StabilizedNanobeam.set_param("group_id", 2)
tabPad = smlay.DeviceTable(StabilizedNanobeam,
                        N_rows, N_cols,
                        param_rows,
                        param_cols
                        )


tabPad.auto_align(50, 50,  numkey=5)


#GeomView(tab.get_geometries().flatten(), fix_aspect_ratio=False)




themask_stabilized_nanobeams.addDeviceTable(tabPad, 3000, 0)
themask_stabilized_nanobeams.exportGDS()