# Basic classes

import samplemaker.layout as smlay
import samplemaker.makers as sm
from samplemaker.devices import Device, registerDevicesInModule
from samplemaker.plotly_viewers import DeviceInspect, GeomView
from samplemaker.baselib.waveguides import BaseWaveguideSequencer, BaseWaveguidePort
from samplemaker.routers import WaveguideConnect
from samplemaker.shapes import GeomGroup
# Create a simple mask layout
import numpy as np
import math
themask = smlay.Mask("Adiabatic Vertical Taper Mask")
import samplemaker.devices


#_____________________________________________________________________________________________
#
#  Nanobeam
#_____________________________________________________________________________________________

class Nanobeam(Device):

    def __init__(self):
        super().__init__()
        self.relevant_params = []

    def initialize(self):
        self.set_name("NANOBEAM")
        self.set_description("Nanobeam with theters")
    
    def write_notes(self): 
        p = self.get_params()
        t = ""
        for key in self.relevant_params:
            t += f"{key}: {p[key]}\n"
        text = sm.make_text(0, -p["theter_length"]-p["top_width"], t, 1, 0.2, to_poly=False, numkey=8, layer=10)
        return text
    
    def set_relevant_params(self, keys):
        p = self.get_params()
        for key in keys:
            if key in p.keys():
                self.relevant_params.append(key)

            

    
    
    def parameters(self):
        self.addparameter("taper_length", 8, "Length of each taper")
        self.addparameter("top_length", 20, "Length of top waveguide")
        self.addparameter("top_width", 0.5, "Width of top waveguide")
        self.addparameter("top_tip_width", 0.05, "Width of top waveguide tip")
        self.addparameter("n_points", 1000, "Number of points used in the taper", param_type=int)
        self.addparameter("m", 0.85, "Top taper shape factor")
        self.addparameter("theters_distance", 5, "Distance between theters")
        self.addparameter("theter_length", 1, "Length of theters")   
        self.addparameter("theter_width", 0.2, "Width of theters")
        self.addparameter("offset_from_tip", 1, "Offset from the tip of the top waveguide as a fraction of the taper length", param_range=[0,1.5])
        self.addparameter("group_id", 0, "Group ID for the device", param_type=int,  param_range=(0, 200))
        self.addparameter("row_id", 0, "Row ID for the device", param_type=int, param_range=(0, 200))
        self.addparameter("col_id", 0, "Column ID for the device", param_type=int, param_range=(0, 200))
        self.addparameter("shrink_um", 0.015, "Shrinkage in microns", param_range=(0.0, 0.5))
        self.addparameter("one_side_theters", 0, "If 0, the theters are only on one side of the top waveguide, if 1, the theters are on both sides", param_type=int, param_range=(0, 1))
        self.addparameter("atom_radius", 0.080, "Radius of the atom in um", param_range=(0.0, 1.0))
        self.addparameter("lattice_constant", 0.350, "Lattice constant in um", param_range=(0.0, 1.0))
        self.addparameter("N_atoms_left", 10, "Number of atoms on the left side of the defect", param_type=int, param_range=(0, 100))
        self.addparameter("N_atoms_right", 10, "Number of atoms on the right side of the defect", param_type=int, param_range=(0, 100))
        self.addparameter("gap", 0.880, "Gap between the centers of the two atoms that form the cavity", param_range=(0.0, 2.0))

    def geom(self):
        p = self.get_params()
        length_total = p["top_length"] + 2*p["taper_length"] 
        
        #top waveguide
        seq = [
            ["S", p["top_length"]/2],
        ]

        sequencer = BaseWaveguideSequencer(seq)
        sequencer.reset()
        sequencer.options["defaultWidth"] = p["top_width"]
        sequencer.options["wgLayer"] = 1

        top = sm.GeomGroup()
        top += sequencer.run()
        
       
        taper_length = p["taper_length"]
        n_points = p["n_points"]
        x = np.linspace(0,p["taper_length"], p["n_points"])
        m = p["m"]
        w1 = p["top_width"]
        w2 = p["top_tip_width"]
        a = (w1-w2)/(taper_length**m)
        w = a*(taper_length-x)**m + w2

        top += sm.make_tapered_path(x, np.zeros(n_points), w,1).translate(p["top_length"]/2, 0)        
        top += top.copy().mirrorX(0)

        #bottom waveguide
        seq = [
            ["S", p["top_length"]/2],
        ]       
       
        XP  = length_total /2
        YP  = 0
        bending_radius = 10


        def WaveguideConnector(port1, port2, width):
            res = WaveguideConnect(port1, port2, bending_radius)
            if res[0] == True: 
                so = BaseWaveguideSequencer(res[1])
                so.options = sequencer.options
                so.options["defaultWidth"] = width
                g = so.run()
                g.rotate_translate(port1.x0, port1.y0, math.degrees(port1.angle()))
                return g
            else: 
                return GeomGroup()

        p1 = BaseWaveguidePort(-XP, YP, "west",p["top_tip_width"], "p1")
        p2 = BaseWaveguidePort(XP, YP, "east", p["top_tip_width"], "p2")

        p1.connector_function = WaveguideConnector
        p2.connector_function = WaveguideConnector
        self.addlocalport(p1)
        self.addlocalport(p2)


        x_centers_left = []
        x = 0
        for i in range(p["N_atoms_left"]):
            x_centers_left.append(x)
            x += p["lattice_constant"]
        x_centers_left = [x - x_centers_left[-1] -p["gap"]/2 for x in x_centers_left]
        

        x_centers_right = []
        x = 0
        for i in range(p["N_atoms_right"]):
            x_centers_right.append(x)
            x += p["lattice_constant"]
        x_centers_right = [x + p["gap"]/2 for x in x_centers_right]
        
        
        x_atoms_min = x_centers_left[0] - p["atom_radius"] if p["N_atoms_left"] > 0 else 0
        x_atoms_max = x_centers_right[-1] + p["atom_radius"] if p["N_atoms_right"] > 0 else 0

        atom = sm.make_circle(0, 0, p["atom_radius"], layer=2, to_poly=True, vertices=64)
        atoms_left = GeomGroup()
        for i,x in enumerate(x_centers_left):
            atoms_left += atom.copy().translate(x, 0)
        
        atoms_right = GeomGroup()
        for i,x in enumerate(x_centers_right):
            atoms_right += atom.copy().translate(x, 0)
        
        atoms = atoms_left + atoms_right

        

        
        # We now focus on the theters
        theter = sm.make_rect(0, 0, p["theter_width"], p["top_width"]  + 2*p["theter_length"]  , layer=1, numkey=5)
        theters = GeomGroup()

        x_valid_max = p["top_length"]/2  + p["taper_length"] - p["offset_from_tip"]*taper_length
        x = -x_valid_max
        x_centers_tapers = []
        while x <= x_valid_max:
            x_centers_tapers.append(x)
            x += p["theters_distance"]
        x_last = x_centers_tapers[-1]
        x_rest =  x_valid_max - x_last
        if x_rest > 0:
            x_centers_tapers = [x + x_rest/2 for x in x_centers_tapers]
        for x in x_centers_tapers:
            if p["N_atoms_left"] == 0 and p["N_atoms_right"] == 0:
                theters += theter.copy().translate(x, 0)
            elif x + p["theter_width"] < x_atoms_min or x - p["theter_width"] > x_atoms_max:
                theters += theter.copy().translate(x, 0)
        if p["one_side_theters"]:
            rect = sm.make_rect(0, 0, p["top_length"]+2*p["taper_length"], p["top_width"]/2 + p["theter_length"], layer=2, numkey=8)
            theters = theters.boolean_difference(rect, 1, 2)
        


        # Device shape before shrinking
        device_shape = top + theters
        device_shape.boolean_union(1)
        
         
            
        

        atoms_no_exp = atoms.copy().poly_resize(p["shrink_um"], layer = 1).set_layer(2)
        top_no_exp = top.copy().poly_resize(p["shrink_um"], layer = 1)
        theters_no_exp = theters.copy().poly_resize(p["shrink_um"], layer = 1)
        no_exp = top_no_exp + theters_no_exp
        no_exp.boolean_union(1)
        no_exp = no_exp.boolean_difference(atoms_no_exp, 1, 2)
       
              
        
        no_exp_binding_box = no_exp.bounding_box().toRect().set_layer(7)
        if p["one_side_theters"]:
            no_exp_binding_box += sm.make_rect(0, 0, p["top_length"]+2*p["taper_length"] + p["shrink_um"]*2, p["top_width"]/2 + p["theter_length"]*0.75, layer=7, numkey=8)
        
         


        low_current = top_no_exp.copy().poly_resize(0.3, layer = 1).set_layer(2) 
        low_current += theters_no_exp.copy().poly_resize( 0.3, layer = 1).set_layer(2)
        low_current = low_current.boolean_difference(no_exp, 2, 1)
        low_current.set_layer(1)
        low_current = low_current.boolean_difference(no_exp_binding_box.copy().poly_outlining(10, layer = 7), 1, 7)

        

        
        
        high_current = top_no_exp.copy().poly_resize(0.150, layer = 1).set_layer(2)
        high_current += theters_no_exp.copy().poly_resize(0.150, layer = 1).set_layer(2)
        high_current = no_exp_binding_box.copy().boolean_difference(high_current, 7, 2).set_layer(2)
        
        
        str = f"{p['group_id']}_{p['row_id']}_{p['col_id']}"
        text = sm.make_text(0, 4*(p["theter_length"] + p["top_width"]), str, 2, 0.2, to_poly=True, numkey=2, layer=4)
        
        
        length_total = no_exp_binding_box.bounding_box().width 
        width_total = no_exp_binding_box.bounding_box().height    

        rect_left_low  = sm.make_rect(0, 0, 0.3, p["top_tip_width"]+p["shrink_um"]*2+0.300*2, layer=1, numkey=6)  
        rect_right_low = sm.make_rect(0, 0, 0.3, p["top_tip_width"]+p["shrink_um"]*2+0.300*2, layer=1, numkey=4)
    
        
        rect_left_high  = sm.make_rect(0, p["top_width"]/2 + p["theter_length"]+p["shrink_um"], 1, width_total, layer=2, numkey=9)
        rect_right_high = sm.make_rect(0, p["top_width"]/2 + p["theter_length"]+p["shrink_um"], 1, width_total, layer=2, numkey=7)

        rect_left_sub = sm.make_rect(0, 0, 0.15, p["top_tip_width"]+p["shrink_um"]*2+0.15*2, layer=2, numkey=6)  
        rect_right_sub= sm.make_rect(0, 0, 0.15, p["top_tip_width"]+p["shrink_um"]*2+0.15*2, layer=2, numkey=4)
        
        rect_left_high.boolean_difference(rect_left_sub, 2, 2)
        rect_right_high.boolean_difference(rect_right_sub, 2, 2)


         
        low_current += rect_left_low.translate(-0.5*length_total, 0)
        low_current += rect_right_low.translate(0.5*length_total, 0)
        low_current.boolean_union(1)

        high_current += rect_right_high.translate(0.5*length_total-0.01, 0)
        high_current += rect_left_high.translate(-0.5*length_total+0.01, 0)
        high_current.boolean_union(2)

         

        
        nanobeam = low_current + high_current + text  
        nanobeam += self.write_notes()
        
        #avc += bot  
        return nanobeam
    








coupon = Nanobeam.build()
print(coupon.get_params())
print(coupon.relevant_params)
coupon.set_relevant_params(["taper_length", "top_length", "top_width", "top_tip_width"])
print(coupon.relevant_params)

#DeviceInspect(coupon, fix_aspect_ratio=True)
registerDevicesInModule(__name__)
geom = coupon.geom()
#GeomView(geom, fix_aspect_ratio=True, plot_height = 600)






#_____________________________________________________________________________________________
#
#  Nanobeam with Pads
#_____________________________________________________________________________________________
class FullDevice(Nanobeam):
    def initialize(self):
        self.set_name("FULLDEVICE")
        self.set_description("full device")
        
    def parameters(self):
        super().parameters()
        self.addparameter("pad_edge", 100, "Edge of the pad")
        self.addparameter("offset_from_corner_um", 5, "Offset from the corner of the pad")
        self.addparameter("theters_distance_pad", 5, "Distance between theters in the pad")
        self.addparameter("theter_length_pad", 1, "Length of theters in the pad")
        self.addparameter("theter_width_pad", 1, "Width of theters in the pad")
        self.addparameter("aperture_size", 2, "Size of the aperture")
        self.addparameter("holes_distance_from_edge", 5, "Distance of the holes from the edge of the pad")
        self.addparameter("holes_distance", 5, "Distance between holes") 


        

    def geom(self):
        p = self.get_params()
        
        
        
        nanobeam_dev = Nanobeam.build()
        p2 = nanobeam_dev.get_params()
        for key in p2.keys():
            nanobeam_dev.set_param(key, p[key])
            
        
        nanobeam = nanobeam_dev.geom()
        
        nanobeam_width_original = 2*p["theter_length"]+p["top_width"]
        pads = sm.make_rect(0, 0, p["pad_edge"], 2*p["pad_edge"]+ nanobeam_width_original , layer = 3, numkey=5)
        theters = GeomGroup()

        
        theter = sm.make_rect(0, p["pad_edge"]+nanobeam_width_original/2, p["theter_width_pad"], p["theter_length_pad"]  , layer=3, numkey=2)
        theters_x = GeomGroup()
        
        x_valid_max = p["pad_edge"]/2  - p["offset_from_corner_um"]
        x = -x_valid_max 
        while x <= x_valid_max:
            theters_x += theter.copy().translate(x, 0)
            x_last = x
            x += p["theters_distance"]
        x_rest =  x_valid_max - x_last
        if x_rest > 0:
            theters_x.translate(x_rest/2, 0)
        theters_x += theters_x.copy().rotate(0, 0, 180)

    

        theter = sm.make_rect(-p["pad_edge"]/2, 0, p["theter_length_pad"], p["theter_width_pad"], layer=3, numkey=6)
        y_valid_max = p["pad_edge"]+ nanobeam_width_original/2 - p["offset_from_corner_um"]
        y = -y_valid_max
        theters_y = GeomGroup()
        while y <= y_valid_max:
            theters_y += theter.copy().translate(0, y)
            y_last = y
            y += p["theters_distance"]
        y_rest =  y_valid_max - y_last
        if y_rest > 0:
            theters_y.translate(0, y_rest/2)
        theters_y += theters_y.copy().rotate(0, 0, 180)

        pads+= theters_x + theters_y
        pads.boolean_union(3)


        hole = sm.make_rect(0, 0, p["aperture_size"], p["aperture_size"], layer=2, numkey=5)
        global themask
        themask.addCell("HOLE", hole)
        holes = GeomGroup()

        distance = p["holes_distance"]+p["aperture_size"]
        x_valid_max = p["pad_edge"] / 2 - p["holes_distance_from_edge"] 
        y_valid_max = p["pad_edge"] + nanobeam_width_original/2 - p["holes_distance_from_edge"] 
        x = -x_valid_max
        x_centers = []
        while x <= x_valid_max:
            x_centers.append(x)
            x_last = x
            x += distance
        x_rest = x_valid_max - x_last
        x_centers = [xc + x_rest/2 for xc in x_centers]

        y = - y_valid_max
        y_centers = []
        while y <= y_valid_max:
            y_centers.append(y)
            y_last = y
            y += distance
        y_rest = y_valid_max-y_last
        y_centers = [yc + y_rest/2 for yc in y_centers]
       

        # Exclude points from the centers grid
        nanobeam_length_no_shrink = nanobeam.select_layers([1,2]).bounding_box().width
        nanobeam_width_no_shrink = nanobeam.select_layers([1,2]).bounding_box().height
        holes = GeomGroup()
        for xc in x_centers:
            for yc in y_centers:
                if abs(xc) <= nanobeam_length_no_shrink / 2 and abs(yc) <= nanobeam_width_no_shrink/2:
                    pass
                else:
                    holes+=hole.copy().translate(xc, yc)
        pads = pads.boolean_difference(holes, 3, 2)
        pads.boolean_union(3)




        bb = pads.bounding_box()
        rbox = bb.toRect()
        rbox.set_layer(4)
        pads = rbox.boolean_difference(pads, 4, 3)
        pads.boolean_union(4)
        pads.poly_resize(-p["shrink_um"], layer=4)
        
        nanobeam_bb = nanobeam.select_layers([1,2]).bounding_box()
        rbox2 = nanobeam_bb.toRect()
        rbox2.set_layer(5)
        pads = pads.boolean_difference(rbox2, 4, 5)
        pads = pads.set_layer(3)
        
        nanobeam = nanobeam.select_layers([1,2,4])
        FullDevice =  nanobeam + pads 
        text = FullDevice.select_layer(4)

        text.translate(0, p["pad_edge"]*1.1)
        FullDevice+= self.write_notes()
        
        
        return FullDevice

FullDevice.build()
#DeviceInspect(FullDevice, fix_aspect_ratio=True)
registerDevicesInModule(__name__)
geom = FullDevice.build().geom()
#GeomView(geom, fix_aspect_ratio=True, plot_height = 600)   

# Basic classes

import samplemaker.layout as smlay
import samplemaker.makers as sm
from samplemaker.devices import Device, registerDevicesInModule
from samplemaker.plotly_viewers import DeviceInspect, GeomView
from samplemaker.baselib.waveguides import BaseWaveguideSequencer, BaseWaveguidePort
from samplemaker.routers import WaveguideConnect
from samplemaker.shapes import GeomGroup
# Create a simple mask layout
import numpy as np
import math
themask = smlay.Mask("Adiabatic Vertical Taper Mask")
import samplemaker.devices

from mydevices import *

#__________________________________________________
#  PHOTONIC CRYSTAL
#__________________________________________________
import samplemaker.phc as sphc

class PhCNanobeam(Nanobeam):
    def initialize(self):
        self.set_name("PHC_NANOBEAM")
        self.set_description("full device")
        
    def parameters(self):
        super().parameters()
        self.addparameter("pad_edge", 100, "Edge of the pad")
        self.addparameter("offset_from_corner_um", 5, "Offset from the corner of the pad")
        self.addparameter("theters_distance_pad", 5, "Distance between theters in the pad")
        self.addparameter("theter_length_pad", 1, "Length of theters in the pad")
        self.addparameter("theter_width_pad", 1, "Width of theters in the pad")
        self.addparameter("aperture_size", 2, "Size of the aperture")
        self.addparameter("holes_distance_from_edge", 5, "Distance of the holes from the edge of the pad")
        self.addparameter("holes_distance", 5, "Distance between holes") 
        self.addparameter("lattice_mirror", 0.45, "Lattice constant of the mirror")
        self.addparameter("radius_mirror", 0.1,  "Mirror radius")


        

    def geom(self):
        p = self.get_params()
        
        
        
        nanobeam_dev = Nanobeam.build()
        p2 = nanobeam_dev.get_params()
        for key in p2.keys():
            nanobeam_dev.set_param(key, p[key])

        
        nanobeam_dev.set_param("theter_length", 0)
        p["theter_length"] = 0
        
        nanobeam = nanobeam_dev.geom()



        
        nanobeam_width_original = 2*p["theter_length"]+p["top_width"]
        pads = sm.make_rect(0, 0, p["pad_edge"], 2*p["pad_edge"]+ nanobeam_width_original , layer = 3, numkey=5)
        theters = GeomGroup()

        
        theter = sm.make_rect(0, p["pad_edge"]+nanobeam_width_original/2, p["theter_width_pad"], p["theter_length_pad"]  , layer=3, numkey=2)
        theters_x = GeomGroup()
        
        x_valid_max = p["pad_edge"]/2  - p["offset_from_corner_um"]
        x = -x_valid_max 
        while x <= x_valid_max:
            theters_x += theter.copy().translate(x, 0)
            x_last = x
            x += p["theters_distance"]
        x_rest =  x_valid_max - x_last
        if x_rest > 0:
            theters_x.translate(x_rest/2, 0)
        theters_x += theters_x.copy().rotate(0, 0, 180)

    

        theter = sm.make_rect(-p["pad_edge"]/2, 0, p["theter_length_pad"], p["theter_width_pad"], layer=3, numkey=6)
        y_valid_max = p["pad_edge"]+ nanobeam_width_original/2 - p["offset_from_corner_um"]
        y = -y_valid_max
        theters_y = GeomGroup()
        while y <= y_valid_max:
            theters_y += theter.copy().translate(0, y)
            y_last = y
            y += p["theters_distance"]
        y_rest =  y_valid_max - y_last
        if y_rest > 0:
            theters_y.translate(0, y_rest/2)
        theters_to_remove =theter.copy()
        theters_y = theters_y.boolean_difference(theters_to_remove, 3, 3)
        theters_y += theters_y.copy().rotate(0, 0, 180)

        pads+= theters_x + theters_y
        pads.boolean_union(3)
        
        def circular_atom_cell(x,y,params):
            return sm.make_circle(x,y, params[0], to_poly=True, vertices=16).set_layer(1)
        
        crystal = sphc.Crystal()
        Nx = int(p["pad_edge"]//p["lattice_mirror"])
        Ny = 9
        r = p["radius_mirror"]

        crystal =  crystal.triangular_box(Nx=int((Nx-1)//2), Ny=int((Ny-1)//2), Nparams=1)
        print(crystal.xpts)
        print(crystal.ypts)
    

        
        phc = sphc.make_phc(crystal=crystal, scaling = p["lattice_mirror"],  cellparams = [r], x0=0, y0=0)
        phc_bb = phc.bounding_box()
        phc = phc.copy().translate(0, phc_bb.height/2 + p["top_width"]/2 + p["radius_mirror"]) + phc.copy().translate(0, - phc_bb.height/2 - p["top_width"]/2 - p["radius_mirror"])
        phc.set_layer(1)
        phc.all_to_poly(32)
        phc_bb2 = phc.bounding_box()
        pads = pads.boolean_difference(phc, 3, 1)
        


         



        hole = sm.make_rect(0, 0, p["aperture_size"], p["aperture_size"], layer=2, numkey=5)
        global themask
        themask.addCell("HOLE", hole)
        holes = GeomGroup()

        distance = p["holes_distance"]+p["aperture_size"]
        x_valid_max = p["pad_edge"] / 2 - p["holes_distance_from_edge"] 
        y_valid_max = p["pad_edge"] + nanobeam_width_original/2 - p["holes_distance_from_edge"] 
        x = -x_valid_max
        x_centers = []
        while x <= x_valid_max:
            x_centers.append(x)
            x_last = x
            x += distance
        x_rest = x_valid_max - x_last
        x_centers = [xc + x_rest/2 for xc in x_centers]

        y = - y_valid_max
        y_centers = []
        while y <= y_valid_max:
            y_centers.append(y)
            y_last = y
            y += distance
        y_rest = y_valid_max-y_last
        y_centers = [yc + y_rest/2 for yc in y_centers]
       


        
        # Exclude points from the centers grid
        nanobeam_length_no_shrink = nanobeam.select_layers([1,2]).bounding_box().width
        nanobeam_width_no_shrink = nanobeam.select_layers([1,2]).bounding_box().height
        holes = GeomGroup()
        for xc in x_centers:
            for yc in y_centers:
                if abs(xc) <= phc_bb2.width / 2 and abs(yc) <= phc_bb2.height/2:
                    pass
                else:
                    holes+=hole.copy().translate(xc, yc)
        pads = pads.boolean_difference(holes, 3, 2)
        pads.boolean_union(3)


        bb = pads.bounding_box()
        rbox = bb.toRect()
        rbox.set_layer(4)
        pads = rbox.boolean_difference(pads, 4, 3)
        pads.boolean_union(4)
        pads.poly_resize(-p["shrink_um"], layer=4)
        
        nanobeam_bb = nanobeam.select_layers([1,2]).bounding_box()
        rbox2 = nanobeam_bb.toRect()
        rbox2.set_layer(5)
        pads = pads.boolean_difference(rbox2, 4, 5)
        pads = pads.set_layer(3)

        separator = sm.make_rect(0, 0, p["pad_edge"], p["top_width"], layer=4, numkey=5)
        separator.poly_resize(-p["shrink_um"], layer=4) 
        separator = separator.boolean_difference(rbox2, 4, 5)
        pads+=separator.set_layer(3)

        
        nanobeam = nanobeam.select_layers([1,2,4])
        FullDevice =  nanobeam + pads 
        text = FullDevice.select_layer(4)

        text.translate(0, p["pad_edge"]*1.1)
        FullDevice+= self.write_notes()
        
        
        return FullDevice

PhCNanobeam.build()
#DeviceInspect(PhCNanobeam, fix_aspect_ratio=True)

registerDevicesInModule(__name__)
geom = PhCNanobeam.build().geom()
#GeomView(geom, fix_aspect_ratio=True, plot_height = 600)   






