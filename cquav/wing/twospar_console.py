from math import sin, cos, tan, sqrt, pi, exp, acos, asin, atan
from copy import deepcopy
from time import time

import numpy as np
from scipy.optimize import minimize_scalar, shgo, brentq, direct
import cadquery as cq
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
import splinecloud_scipy as scloud
import matplotlib.pyplot as plt

from .profiles.airfoil_section import AirfoilSection
from .meta.geom import Wing
from ..constants import *


class TwoSparWingConsole(Wing):
    """ 
    Wing console of the tapered swept wing, made of foam, and reinforced with two simple thin-walled spars. 
    """
    
    def __init__(self, ccls_surface, span_efficiency, root_airfoil, 
                 tip_airfoil=None, span=2000, root_twist_angle=0, tip_twist_angle=0, dihedral=0,
                 sweep_angle=None, chord_align_factor=0.0, aspect_ratio=9, taper_ratio=0.7, 
                 front_spar_rel_position=0.25, front_spar_thickness=None, 
                 rear_spar_rel_position=0.7, rear_spar_thickness=None,
                 shell_thickness=None, materials=None):

        """
            Initialization parametrs:
            ========================
            @root_airfoil - Airfoil object for the root chord
            @tip_airfoil - Airfoil object for the tip chord
            @root_twist_angle - root chord twist angle (along the vector normal to the chord plane)
            @tip_twist_angle - tip chord twist angle (along the vector normal to the chord plane)
            @span - full win span(two consoles), [mm]
            @sweep_angle - front edge sweep angle, [degrees]
            @aspect_ratio - full wing aspect ratio = (root chord)**2/(wing area)
            @taper_ratio - wing taper ratio = (tip chord)/(root chord)
            @dihedral - dihedral angle, [degrees]
            @chord_align_factor - parameter defining the position of the tip chord relative to the root chod 
                                  (used todetermine sweep angle if sweep_angle is None)
            @front_spar_rel_position - relative (to the chord length) position of the font spar
            @front_spar_thickness - front spar thickness
            @rear_spar_rel_position - relative (to the chord length) position of the rear spar
            @rear_spar_thickness - rear spar thickness
            @shell_thickness - thikness of the wing console shell, [mm]. 
                               If shell_thickness == None - no shell will be constructed
            @materials - dictionary with materials for each part of the wing structure in the next form:
                         {
                            "shell": <shell material instance>,
                            "foam": <foam material instance>,
                            "spar": <spar material instance"
                         }
        """

        super().__init__(ccls_surface, span_efficiency, sweep_angle, aspect_ratio, span, root_airfoil, 
                         taper_ratio=taper_ratio, dihedral=dihedral)

        tip_airfoil = tip_airfoil or root_airfoil
        self.length = 0.5 * span
        self.root_twist_angle = root_twist_angle
        self.tip_twist_angle = tip_twist_angle
        self.tip_vertical_offset = self.length*sin(dihedral*pi/180)

        if sweep_angle is None:           
            self.tip_chord_setback = (self.root_chord - self.tip_chord) * chord_align_factor
            self.sweep_angle = atan(2*self.tip_chord_setback/self.span) * 180 / pi
        else:
            self.tip_chord_setback = tan(self.sweep_angle * pi / 180)*self.length

        self.root_section = AirfoilSection(
            root_airfoil, chord=self.root_chord, twist_angle=self.root_twist_angle,
        )
        self.tip_section = AirfoilSection(
            tip_airfoil, chord=self.tip_chord, twist_angle=self.tip_twist_angle,
            x_offset=self.tip_chord_setback, y_offset=self.tip_vertical_offset, z_offset=self.length,
        )

        self.front_spar_rel_position = front_spar_rel_position
        self.rear_spar_rel_position = rear_spar_rel_position
        self.front_spar_thickness = front_spar_thickness or 0.01*self.root_chord
        self.rear_spar_thickness = rear_spar_thickness or 0.01*self.root_chord
        self.shell_thickness = shell_thickness

        self.materials = materials or {}
        self.spar_material = None
        self.shell_material = None
        self.foam_material = None

        if materials:
            self.assign_materials(materials)

        self.build_geometry()
   
    def assign_materials(self, materials: dict):
        self.materials = materials
        if self.materials:
            if "spar" in materials.keys():
                self.spar_material = materials["spar"]
            
            if "shell" in materials.keys():
                self.shell_material = materials["shell"]
                
            if "foam" in materials.keys():
                self.foam_material = materials["foam"]

    def update_geom_params(self):
        self.update_dimensions()
        self.update_aero_params()

        self.length = 0.5 * self.span
        self.tip_chord_setback = tan(self.sweep_angle * pi / 180)*self.length

        root_airfoil = self.root_section.airfoil
        tip_airfoil = self.tip_section.airfoil

        self.root_section = AirfoilSection(
            root_airfoil, chord=self.root_chord, twist_angle=self.root_twist_angle,
        )
        self.tip_section = AirfoilSection(
            tip_airfoil, chord=self.tip_chord, twist_angle=self.tip_twist_angle,
            x_offset=self.tip_chord_setback, y_offset=self.tip_vertical_offset, z_offset=self.length,
        )

    def set_new_aspect_ratio(self, ar):
        print(f"Setting new aspect ratio. AR = {ar}")
        self.ar = ar
        self.update_geom_params()

    def set_new_span(self, span):
        print(f"Setting new span. Span = {span}")
        self.span = span
        self.update_geom_params()

    def build_geometry(self, length=None):         
        if length:
            self.length = length
        
        self._airfoil_body = self.build_airfoil_body()
        self.front_spar = self.build_spar_body(self.front_spar_rel_position, self.front_spar_thickness)
        self.rear_spar = self.build_spar_body(self.rear_spar_rel_position, self.rear_spar_thickness)
        self.foam = self.build_foam_body()
        
        x_ac, z_ac = self.console_lift_force_position()
        self.aero_center = (
            cq.Workplane()
            .workplane(offset=z_ac)
            .move(x_ac, 0)
            .sphere(5)
        )

        self.original_parts = {
            "_airfoil_body": self._airfoil_body,
            "front_spar": self.front_spar,
            "rear_spar": self.rear_spar,
            "foam": self.foam,
        }

        self.transformed_parts = [
            "front_spar",
            "rear_spar",
            "foam",
            "_airfoil_body",
            "aero_center"
        ]

        if self.shell_thickness:
            self.shell = self.build_outset_shell()
            self.transformed_parts.append("shell")
            self.original_parts['shell'] = self.shell
  
    def build_airfoil_body(self):
        root_airfoil_sketch = self.root_section.sketch.copy()
        tip_airfoil_sketch = self.tip_section.sketch.copy()

        airfoil_body = (
            cq.Workplane()
            .placeSketch(
                root_airfoil_sketch.moved(self.root_section.location), 
                tip_airfoil_sketch.moved(self.tip_section.location))
            .loft()
        )

        return airfoil_body

    def build_spar_body(self, spar_rel_position, spar_thickness):
        spar_root_x_coord = spar_rel_position * self.root_chord
        spar_tip_x_coord = spar_rel_position * self.tip_chord
        tip_chord_shift = cq.Vector(self.tip_chord_setback, 0, self.length)

        spar_root_section = (
            cq.Sketch()
            .push([cq.Vector(spar_root_x_coord, 0)])
            .rect(spar_thickness, self.root_chord*2)
        )

        spar_tip_section = (
            cq.Sketch()
            .push([cq.Vector(spar_tip_x_coord, 0)])
            .rect(spar_thickness, self.root_chord*2)
        )

        spar_box = (
            cq.Workplane()
            .placeSketch(
                spar_root_section, 
                spar_tip_section.moved(cq.Location(tip_chord_shift)))
            .loft()
        )

        spar_body = spar_box.intersect(self._airfoil_body)

        return spar_body

    def build_foam_body(self):
        foam_body = (
            cq.Workplane()
            .add(self._airfoil_body)
            .cut(self.front_spar)
            .cut(self.rear_spar)
        )

        return foam_body

    def build_outset_shell(self):
        shell = (
            cq.Workplane()
            .add(self._airfoil_body)
            .faces("%PLANE")
            .shell(self.shell_thickness)
        )

        return shell

    def translate(self, x_offset, y_offset, z_offset):
        for part_name in self.transformed_parts:
            part = getattr(self, part_name)
            translated = part.translate((x_offset, y_offset, z_offset))
            setattr(self, part_name, translated)

        ## outset shell is error prone and may not transform correctly
        if self.shell_thickness:
            self.shell = self.build_outset_shell()

    def rotate_xyz(self, axis, angle):
        p1 = (0, 0, 0)
        if axis in ("X", "x"):
            p2 = (1, 0, 0)
        elif axis in ("Y", "y"):
            p2 = (0, 1, 0)
        elif axis in ("Z", "z"):
            p2 = (0, 0, 1)

        for part_name in self.transformed_parts:
            part = getattr(self, part_name)
            translated = part.rotate(p1, p2, angle)
            setattr(self, part_name, translated)

        ## outset shell is error prone and may not transform correctly
        if self.shell_thickness:
            self.shell = self.build_outset_shell()

    def mirror(self, plane):
        for part_name in self.transformed_parts:
            part = getattr(self, part_name)
            translated = part.mirror(plane)
            setattr(self, part_name, translated)

        ## outset shell is error prone and may not transform correctly
        if self.shell_thickness:
            self.shell = self.build_outset_shell()

    @property
    def moment_Yx(self):
        ## returns moment about absolute CS center in vertical plane
        front_spar_mass = self.get_components_mass([self.front_spar], self.spar_material)
        rear_spar_mass = self.get_components_mass([self.rear_spar], self.spar_material)
        foam_mass = self.get_foam_mass()
        shell_mass = self.get_shell_mass()

        front_spar_moment = self.front_spar.vals().Center().x * 1e-3 * front_spar_mass
        rear_spar_moment = self.front_spar.vals().Center().x * 1e-3 * rear_spar_mass
        foam_moment = self.foam.vals()[0].Center().x * 1e-3 * foam_mass
        shell_moment = self._airfoil_body.vals()[0].Center().x * 1e-3 * shell_mass

        return front_spar_moment + rear_spar_moment + shell_moment

    @property
    def moment_Yz(self):
        ## returns moment about absolute CS center in vertical plane
        front_spar_mass = self.get_components_mass([self.front_spar], self.spar_material)
        rear_spar_mass = self.get_components_mass([self.rear_spar], self.spar_material)
        foam_mass = self.get_foam_mass()
        shell_mass = self.get_shell_mass()

        front_spar_moment = self.front_spar.vals().Center().z * 1e-3 * front_spar_mass
        rear_spar_moment = self.front_spar.vals().Center().z * 1e-3 * rear_spar_mass
        foam_moment = self.foam.vals()[0].Center().z * 1e-3 * foam_mass
        shell_moment = self._airfoil_body.vals()[0].Center().z * 1e-3 * shell_mass

        return front_spar_moment + rear_spar_moment + shell_moment

    def get_aero_center(self):
        ac = self.aero_center.vals()[0].Center()
        return ac.x, ac.y, ac.z

    def xa(self):
        return self.get_aero_center()[0]

    ## Mass properties
    def get_components_mass(self, components, material):
        volume = 0
        
        for c in components:
            volume += c.solids().vals()[0].Volume()
            
        return volume * (1e-3)**3 * material.density
    
    def get_spars_mass(self):
        return self.get_components_mass([self.front_spar, self.rear_spar], self.spar_material)
        
    def get_foam_mass(self):
        return self.get_components_mass([self.foam], self.foam_material)
                 
    def get_shell_mass(self, airfoil_body=None):
        ## computation of thin BREP volume is error-prone - in some cases produces wrong or negative results
        ## here a thin shell mass can be computed using outer area of the airoil body

        if not self.shell_thickness:
            return 0

        airfoil_body = airfoil_body or self._airfoil_body
        wetted_area = 0
        for f in airfoil_body.faces().vals():
            wetted_area += f.Area()

        for f in airfoil_body.faces("<Z or >Z").vals():
            wetted_area -= f.Area()

        shell_volume = wetted_area * self.shell_thickness
        shell_mass = 1e-9 * shell_volume * self.shell_material.density  
                  
        return shell_mass

    def get_console_mass(self):
        spars_mass = self.get_spars_mass()
        foam_mass = self.get_foam_mass()
        shell_mass = self.get_shell_mass()
        
        return spars_mass + foam_mass + shell_mass   

    @property
    def mass(self):
        return self.get_console_mass()

    def length_cut(self, offset, components="all", merge_result=False, group_on_workplane=True):
        if offset <= 0:
            raise ValueError("Offset should be positive")
            
        if components == "all":
            components = list(self.original_parts.values())
        elif type(components) == str:
            components = [self.original_parts[components]]

        max_chord = max(self.root_section.chord, self.tip_section.chord)

        max_x = max([shape.vals()[0].BoundingBox().xlen for shape in components])
        max_y = max([shape.vals()[0].BoundingBox().ylen for shape in components])

        cutbox = (
            cq.Workplane()
            .workplane(offset=offset)
            .rect(max_x*2.5, max_y*2.5) ## slighly increase required size of the cut box
            .extrude(self.length)
        )
        
        wp = cq.Workplane()
        sliced_components = []
        for component in components:
            component_cut = component.cut(cutbox)
            
            if merge_result:
                wp = wp.union(component_cut)
            else:
                wp = wp.add(component_cut.solids())
            
        if group_on_workplane or merge_result:
            return wp
            
        else:
            return sliced_components
        
    @classmethod
    def eval_section_inertia_moments(self, body, selector="start"):
        """
        Evaluates moments of inartia of the component section
        
        @body: cq.Workplane containing solid
        """
        
        ## convert verbose selector to proper CadQuery selector
        if selector == "start":
            selector = "<Z"
        elif selector == "end":
            selector = ">Z"
            
        section_occ = body.faces(selector).toOCC()
        
        properties = GProp_GProps()
        BRepGProp.SurfaceProperties_s(section_occ, properties)
        matrix_of_inertia = properties.MatrixOfInertia()
        
        Ixx, Iyy, Izz = matrix_of_inertia.Value(1,1), matrix_of_inertia.Value(2,2), matrix_of_inertia.Value(3,3)
        
        return Ixx, Iyy, Izz

    @classmethod
    def eval_body_slice_inertia_moments(self, body, height):
        """
        Evaluates moments of inartia of the component section
        
        @body: cq.Workplane containing solid
        """
                   
        section_occ = body.section(height).toOCC()
        properties = GProp_GProps()
        BRepGProp.SurfaceProperties_s(section_occ, properties)
        matrix_of_inertia = properties.MatrixOfInertia()
        
        Ixx, Iyy, Izz = matrix_of_inertia.Value(1,1), matrix_of_inertia.Value(2,2), matrix_of_inertia.Value(3,3)
        
        return Ixx, Iyy, Izz
    
    def eval_aero_midspan(self, alpha, velocity, delta=0, n=100):
        """Calculates span offset where the total console lift force is applied"""
        
        root_airfoil = self.root_section.airfoil
        tip_airfoil = self.tip_section.airfoil
        
        alpha_root = self._alpha_actual(alpha, delta+self.root_twist_angle)
        alpha_tip = self._alpha_actual(alpha, delta+self.tip_twist_angle)

        reynolds_root = self.eval_reynolds(velocity, self.root_chord)
        reynolds_tip = self.eval_reynolds(velocity, self.tip_chord)

        cl_root = root_airfoil.eval_cl(alpha_root, reynolds_root)
        cl_tip = tip_airfoil.eval_cl(alpha_tip, reynolds_tip)
        cl_prev = cl_root

        dl = self.length / n
        Cl_area = 0
        Cl_mom = 0
        for i in range(1, n+1):
            ns = (dl*i - dl*0.5) / self.length
            ccls_i = self.ccls_surface(self.ar, ns)

            cl_i = (cl_root - cl_tip) * (self.length - dl*i) / self.length + cl_tip
            cl_i_avg = (cl_i + cl_prev) / 2

            Cl_area += cl_i_avg * ccls_i * dl
            Cl_mom +=  cl_i_avg * ccls_i * dl * (dl*i - dl*0.5)

        midspan = Cl_mom / Cl_area
        
        return midspan
        
    def eval_lift_force(self, alpha, velocity, delta=0, n=20, load_container=None,
            compute_weight_load=False, load_factor=1, g=STANDARD_GRAVITY):

        length = self.length * 1e-3
        dl = length / n
        root_airfoil = self.root_section.airfoil
        tip_airfoil = self.tip_section.airfoil

        alpha_root = self._alpha_actual(alpha, delta+self.root_twist_angle)
        alpha_tip = self._alpha_actual(alpha, delta+self.tip_twist_angle)

        reynolds_root = self.eval_reynolds(velocity, self.root_chord)
        reynolds_tip = self.eval_reynolds(velocity, self.tip_chord)

        c_dihedral = cos(self.dihedral * pi / 180)
        dyn_airpressure = 0.5 * AIR_DENSITY * velocity**2

        cl_root = root_airfoil.eval_cl(alpha_root, reynolds_root)
        cl_tip = tip_airfoil.eval_cl(alpha_tip, reynolds_tip)

        cl_prev = cl_root
        root_chord = self.root_chord*1e-3

        ## compute aerodynamic load acting on each console slice along its length 
        Lift_force = 0  
        prev_cut_mass = 0

        if load_container is None:
            load_container = {}
        
        for i in range(1, n+1):
            ns = (dl*i - 0.5*dl) / length
            ccls_i = self.ccls_surface(self.ar, ns)

            cl_i = (cl_root - cl_tip) * (length - dl*i) / length + cl_tip
            cl_i_avg = (cl_i + cl_prev) / 2
            cl_prev = cl_i

            w_i = dyn_airpressure * cl_i_avg * ccls_i * c_dihedral * root_chord
            
            if compute_weight_load:
                current_cut_mass = self._eval_components_cut_masses(dl*i*1e3)
                console_section_specific_weight = (current_cut_mass - prev_cut_mass) * g * load_factor / dl
                w_i -= console_section_specific_weight
                prev_cut_mass = current_cut_mass
                
            load_container[i] = w_i
            Lift_force += w_i * dl
        
        return Lift_force
 
    def cl_slope(self, alpha1=0, alpha2=5, alpha_unit="degree"):
        ## accounts to possile variable airfoil

        velocity = 35 ## does not influense result
        dyn_airpressure = 0.5 * AIR_DENSITY * velocity**2
        k_dyn = dyn_airpressure * self.root_chord * self.length * 1e-6

        lift1 = self.eval_lift_force(alpha1, velocity, n=100)
        lift2 = self.eval_lift_force(alpha2, velocity, n=100)

        scale = 1 if alpha_unit == "degree" else 180/pi
        
        return scale * (lift2-lift1) / (k_dyn * (alpha2-alpha1))

    def _eval_components_cut_masses(self, root_section_offset, front_spar_cut=None, rear_spar_cut=None, 
                                    foam_cut=None, airfoil_body_cut=None):       
        
        front_spar_cut = front_spar_cut or self.length_cut(root_section_offset, components='front_spar')
        rear_spar_cut = rear_spar_cut or self.length_cut(root_section_offset, components='rear_spar')
        spars_cut_mass = self.get_components_mass([front_spar_cut, rear_spar_cut], self.spar_material)

        foam_cut = foam_cut or self.length_cut(root_section_offset, components='foam')
        foam_cut_mass = self.get_components_mass([foam_cut], self.foam_material)

        if self.shell_thickness:
            airfoil_body_cut = airfoil_body_cut or self.length_cut(root_section_offset, components='_airfoil_body')
            shell_cut_mass = self.get_shell_mass(airfoil_body_cut)
        else:
            shell_cut_mass = 0

        return spars_cut_mass + foam_cut_mass + shell_cut_mass

    def eval_drag_force(self, alpha, velocity, delta=0, n=20):
        ## calculate net drag force

        root_airfoil = self.root_section.airfoil
        tip_airfoil = self.tip_section.airfoil

        alpha_root = self._alpha_actual(alpha, delta+self.root_twist_angle)
        alpha_tip = self._alpha_actual(alpha, delta+self.tip_twist_angle)

        reynolds_root = self.eval_reynolds(velocity, self.root_chord)
        reynolds_tip = self.eval_reynolds(velocity, self.tip_chord)

        c_dihedral = cos(self.dihedral * pi / 180)
        dyn_airpressure = 0.5 * AIR_DENSITY * velocity**2

        cd0_root = root_airfoil.eval_cd(alpha_root, reynolds_root)
        cd0_tip = tip_airfoil.eval_cd(alpha_tip, reynolds_tip)
        cd0_prev = cd0_root

        root_chord = self.root_chord * 1e-3
        tip_chord = self.tip_chord * 1e-3
        chord_prev = root_chord

        length = self.length * 1e-3
        dl = length / n

        zero_drag_force = 0

        for i in range(1, n+1):
            cd0_i = (cd0_root - cd0_tip) * (length - dl*i) / length + cd0_tip
            cd0_i_avg = (cd0_i + cd0_prev) / 2
            cd0_prev = cd0_i

            chord_i = (root_chord - tip_chord) * (length - dl*i) / length + tip_chord
            chord_i_avg = (chord_i + chord_prev) / 2
            chord_prev = chord_i

            drag0_i = dyn_airpressure * cd0_i_avg * chord_i_avg * dl
            zero_drag_force += drag0_i

        cdi = self.induced_drag_coefficient(alpha_root)
        induced_drag_force = dyn_airpressure * cdi * self.area*0.5
        
        return (zero_drag_force + induced_drag_force) * INCREASED_DRAG_FACTOR

    def eval_bend_deflection(self, alpha, velocity, relative=False,
            delta=0, n=20, load_factor=1, compute_weight_load=True, g=STANDARD_GRAVITY, verbose=False,
            extra_tip_lift_force=0, extra_tip_lift_force_arm=0, extra_tip_g_force=0, extra_tip_g_force_arm=0):

        length = self.length * 1e-3
        extra_tip_lift_force_arm = extra_tip_lift_force_arm*1e-3
        extra_tip_g_force_arm = extra_tip_g_force_arm*1e-3

        ## calculate aerodynamic loads on the wingspan sections
        lift_force_eval_time = time()
        W = {}
        self.eval_lift_force(alpha, velocity, delta=delta, n=n, load_container=W)
        W[n+1] = 0
        if verbose:
            print("lift_force_eval_time:", time() - lift_force_eval_time)

        inertia_moments_eval_time = time()
        ## compute console longeron section moments of ineria and console specific weights
        I = {}
        E = self.spar_material.tensile_modulus
        RMC = {}
        Ymax = {}
        prev_cut_mass = 0

        I[1] = (self.eval_section_inertia_moments(self.original_parts['front_spar'], selector="start") + 
                self.eval_section_inertia_moments(self.original_parts['rear_spar'], selector="start"))
        
        Ymax[1] = max(self._get_section_max_offsets(self.front_spar, selector="start")[1], 
                      self._get_section_max_offsets(self.rear_spar, selector="start")[1])

        I_foam = self.eval_section_inertia_moments(self.original_parts['foam'], selector="start")
        EI_foam = self.foam_material.tensile_modulus * I_foam[0]

        if self.shell_thickness:
            I_shell = self.eval_section_inertia_moments(self.original_parts['shell'], selector="start")
            EI_shell = self.shell_material.tensile_modulus * I_shell[0]
        else:
            I_shell = 0
            EI_shell = 0
        
        RMC[1] = E*I[1][0]/(EI_foam + EI_shell + E*I[1][0])

        dl = length / n
        for i in range(1, n+1):
            front_spar_cut = self.length_cut(dl*i*1e3, components='front_spar')
            rear_spar_cut = self.length_cut(dl*i*1e3, components='rear_spar')

            I[i+1] = (self.eval_section_inertia_moments(front_spar_cut, selector="end") + 
                      self.eval_section_inertia_moments(rear_spar_cut, selector="end"))

            # evaluate reduced moment coefficient
            foam_cut = self.length_cut(dl*i*1e3, components='foam')
            I_foam = self.eval_section_inertia_moments(foam_cut, selector="end")
            EI_foam = self.foam_material.tensile_modulus * I_foam[0]

            if self.shell_thickness:
                shell_cut = self.length_cut(dl*i*1e3, components='shell')
                I_shell = self.eval_section_inertia_moments(shell_cut, selector="end")
                EI_shell = self.shell_material.tensile_modulus * I_shell[0]
                
            else:
                I_shell = 0
                EI_shell = 0

            RMC[i+1] = E*I[i+1][0]/(EI_foam + EI_shell + E*I[i+1][0])
            Ymax[i+1] = max(self._get_section_max_offsets(self.front_spar, selector="end")[1], 
                            self._get_section_max_offsets(self.rear_spar, selector="end")[1])
            
            if compute_weight_load:
                current_cut_mass = self._eval_components_cut_masses(dl*i*1e3, 
                    front_spar_cut=front_spar_cut, rear_spar_cut=rear_spar_cut, foam_cut=foam_cut
                )

                console_section_specific_weight = (current_cut_mass - prev_cut_mass) * g * load_factor / dl
                W[i] -= console_section_specific_weight
                prev_cut_mass = current_cut_mass
                
        if verbose:
            print("inertia_moments_eval_time:", time() - inertia_moments_eval_time)

        delta_eval_time = time()
        ## integrate bend deflection caused by aerodynamic and weight loads, acting in opposite direction
        delta_sum = 0
        theta_prev = 0
        delta_prev = 0
        max_bend_stress = 0
        for i in range(1, n+1):
            Mi_extra_lift_start = (length - (i-1)*dl + extra_tip_lift_force_arm) * extra_tip_lift_force
            Mi_extra_lift_end = (length - (i)*dl + extra_tip_lift_force_arm) * extra_tip_lift_force
            Mi_extra_g_start = (length - (i-1)*dl + extra_tip_g_force_arm) * extra_tip_g_force
            Mi_extra_g_end = (length - (i)*dl + extra_tip_g_force_arm) * extra_tip_g_force

            Mi_start = Mi_extra_lift_start - Mi_extra_g_start
            Mi_end = Mi_extra_lift_end - Mi_extra_g_end

            for j in range(i, n+1):
                ## sum all momments starting from the begining of the current section to the tip
                Mi_start += (j-i) * W[j] * dl**2 + W[j] * dl**2 / 2
                ## sum all momments starting from the end of the current section to the tip
                Mi_end += (j-i) * W[j+1] * dl**2 + W[j+1] * dl**2 / 2

            Ixx_start = I[i][0] * (1e-3)**4
            Ixx_end = I[i+1][0] * (1e-3)**4
            RMC_avg = (RMC[i] + RMC[i+1])/2

            theta_i = dl * RMC_avg * (Mi_start/(E * Ixx_start) + Mi_end/(E * Ixx_end)) / 2
            delta_i = delta_prev + dl * (theta_prev + theta_i) / 2

            bend_stress = Mi_start * RMC[i] * Ymax[i]*(1e-3) / Ixx_start ## Pa
            max_bend_stress = max(max_bend_stress, bend_stress)

            theta_prev = theta_i
            delta_prev = delta_i

            delta_sum += delta_i

        if verbose:  
            print("delta_eval_time:", time() - delta_eval_time)
            print()

        if relative:
            return delta_sum / length, max_bend_stress*1e-6

        return delta_sum * 1e3, max_bend_stress*1e-6 ## [mm], [MPa]
    
    def _get_section_max_offsets(self, body, selector="end"):
        """
        Evaluates body section max offsets from the center point
        
        @body: cq.Workplane containing solid
        """

        ## convert verbose selector to proper CadQuery selector
        if selector == "start":
            selector = "<Z"
        elif selector == "end":
            selector = ">Z"
            
        section = body.faces(selector)
        center = section.val().Center()
        section_bb = section.val().BoundingBox()

        max_x_offset = max(
            abs(section_bb.xmin - center.x),
            abs(section_bb.xmax - center.x),
        )

        max_y_offset = max(
            abs(section_bb.ymin - center.y),
            abs(section_bb.ymax - center.y),
        )

        return max_x_offset, max_y_offset
   
    def fit_span_to_required_lift(self, alpha, velocity, required_lift_force, 
            g=STANDARD_GRAVITY, load_factor=1, n=20, delta=0, adopt_result=True, compute_weight_load=True,
            span_min=200, span_max=8000, method="minimize", verbose=True
        ):
 
        def excess_lift_force(span):
            span = float(span)
            if span <= span_min or span >= span_max:
                return 1e6
            
            self.set_new_span(span)
            self.build_geometry()
            
            lift_force = self.eval_lift_force(alpha, velocity, 
                load_factor=load_factor, g=g, delta=delta, compute_weight_load=compute_weight_load,
            )

            error = lift_force*2 - required_lift_force
            if verbose:
                print(f"Current span: {span}")
                print(f"Mass: {self.get_console_mass()}, [kg]")
                print(f"Excess console lift force: {lift_force}, [N]")
                print(f"Lift force error: {error}, [N]")
                print()

            return abs(error)

        print(f"Finding span for the required lift force = {required_lift_force}")
        if method == "minimize":
            result = minimize_scalar(
                excess_lift_force, 
                bounds=[span_min, span_max], 
                method="bounded",
                options={"xatol": 1}
            )

            solved_span = result.x
        
        elif method == "shgo":
            result = shgo(
                excess_lift_force, 
                bounds=[(span_min, span_max)], 
                # options={"xatol": 0.001}
            )

            solved_span = result.x[0]

        if adopt_result:
            self.set_new_span(solved_span)
            self.build_geometry()
 
        print("==============================================")
        print("Wing console span selection result:\n", solved_span)
        
        return solved_span

    def print_summary(self, summary_dict):
        for section, params in summary_dict.items():
            print("--------------")
            print(f"{section}")
            print("--------------")
            for label, (value, unit) in params.items():
                print(f"{label}: {value} {unit}")
        print()

    def flatten_summary_dict_wide(self, summary_dict):
        flat_dict = {}
        for section, params in summary_dict.items():
            for label, (value, unit) in params.items():
                key = f"{label}{', ' + unit if unit else ''}"
                flat_dict[key] = value
        
        return flat_dict

    def stats(self, alpha, velocity, load_factor=1, g=STANDARD_GRAVITY, verbose=True):
        lift_force = self.eval_lift_force(alpha, velocity, load_factor=load_factor, g=g)
        drag_force = self.eval_drag_force(alpha, velocity)

        bend_deflection_rel, bend_stress = self.eval_bend_deflection(alpha, velocity, load_factor=load_factor, g=g, relative=True)
        center_of_pressure = self.console_lift_force_position()
        
        spars_mass = self.get_spars_mass()
        foam_mass = self.get_foam_mass()
        shell_mass = self.get_shell_mass()

        mass = spars_mass + foam_mass + shell_mass
        weight = mass * g * load_factor

        stats_dict = {
            "Flight conditions": {
                "AoA": (alpha, "[degrees]"),
                "Velocity": (velocity, "[m/s]"),
                "Load factor": (load_factor, "")
            },
            "Wing geometry": {
                "Root airfoil": (self.root_section.airfoil.name, ""),
                "Root chord": (self.root_section.chord, "[mm]"),
                "Root section twist angle": (self.root_section.twist_angle, "[degrees]"),

                "Tip airfoil": (self.tip_section.airfoil.name, ""),
                "Tip chord": (self.tip_section.chord, "[mm]"),
                "Tip section twist angle": (self.tip_section.twist_angle, "[degrees]"),

                "Wing sweep angle": (self.sweep_angle, "[degrees]"),
                "Wing span": (self.span, "[mm]"),
                "Wing aspect ratio": (self.ar, ""),
                "Wing area": (self.area, "[m^2]")
            },
            "Mass breakdown": {
                "Console mass": (mass, "[kg]"),
                "Spars mass": (spars_mass, "[kg]"),
                "Relative spars mass": (f"{100*spars_mass/mass:.2f}", "[%]"),
                "Foam mass": (foam_mass, "[kg]"),
                "Relative foam mass": (f"{100*foam_mass/mass:.2f}", "[%]"),
                "Shell mass": (shell_mass, "[kg]"),
                "Relative shell mass": (f"{100*shell_mass/mass:.2f}", "[%]"),

                "Front spar thickness": (self.front_spar_thickness, "[mm]"),
                "Rear spar thickness": (self.rear_spar_thickness, "[mm]"),
                "Shell thickness": (self.shell_thickness, "[mm]"),
            },
            "Strength:": {
                "Relative bend deflection": (bend_deflection_rel * 100, "[%]"),
                "Max bend stress": (bend_stress, "[MPa]")
            },
            "Aerodynamics": {
                "Excess console lift force": (lift_force - weight, "[N]"),
                "Total console lift force": (lift_force, "[N]"),
                "Total console drag force": (drag_force, "[N]"),
                "Aerodynamic quality": (lift_force / drag_force, ""),
                "Lift to weight ratio": (lift_force / weight, ""),
                "Center of aerodynamic pressure (chord, span)": (center_of_pressure, "[mm, mm]")
            }
        }

        if verbose:
            self.print_summary(stats_dict)

        return stats_dict

    def plot_aero_forces_by_alpha(self, velocity, n=100):
        alpha_vals = np.linspace(0, self.alpha_max_lift()-self.root_twist_angle, n)
        LF_vals = [self.eval_lift_force(aoa, velocity)*2 for aoa in alpha_vals]
        DF_vals = [self.eval_drag_force(aoa, velocity)*2 for aoa in alpha_vals]
        LF_to_DF_vals = [LF_vals[i]/DF_vals[i] for i in range(len(alpha_vals))]

        fig, axs = plt.subplots(1, 3)
#         fig.tight_layout()

        axs[0].plot(alpha_vals, LF_vals)
        axs[0].set_title('Lift force')
        axs[0].grid()
        axs[1].plot(alpha_vals, DF_vals, 'tab:orange')
        axs[1].set_title('Drag force')
        axs[1].grid()
        axs[2].plot(alpha_vals, LF_to_DF_vals, 'tab:green')
        axs[2].set_title('Aerodynamic quality')
        axs[2].grid()
