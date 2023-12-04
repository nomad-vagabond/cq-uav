import math

import numpy as np
from scipy.optimize import minimize_scalar, shgo, brentq, direct
import cadquery as cq

from cquav.constants import STANDARD_GRAVITY
from .profile import AirfoilSection, ThreeChamberBoxedWingSection


class RectangularWingConsole:
    """ 
    Aerodynamic wing console reinforced with tree-pockets thin-walled longeron box and constant section along its length.
    Units:
        Length: mm
        Mass: kg
        Volume: mm^3
    """
    
    def __init__(self, root_section, length=1000, box_thickness=None, materials=None, shell_thickness=1, 
                 min_length=50, max_length=5000, min_chord=50, max_chord=600):
        
        if isinstance(root_section, ThreeChamberBoxedWingSection):
            self.box_section = root_section
            self.airfoil_section = root_section.airfoil_section  
        elif isinstance(root_section, AirfoilSection):
            self.airfoil_section = root_section
            self.box_section = ThreeChamberBoxedWingSection(root_section, box_thickness=box_thickness, shell_thickness=shell_thickness)
        else:
            raise AttributeError("'root_section' should be an istance of ThreeChamberBoxedWingSection or AirfoilSection classes")
    
        self.chord = self.box_section.chord
        self.length = length
        self.shell_thickness = shell_thickness
        self.box_thickness = self.box_section.box_thickness
        self.box_material = None
        self.shell_material = None
        self.foam_material = None
        
        self.min_length = min_length
        self.max_length = max_length
        self.min_chord = min_chord
        self.max_chord = max_chord
        
        if materials:
            self.assign_materials(materials)

        self.build_geometry()
    
    def build_geometry(self, length=None):         
        if length:
            self.length = length
        
        self.inner_body = self.build_airfoil_body()    
        self.foam = self.build_foam_body()
        self.front_box = self.build_box_compartment(self.box_section.front_section)
        self.central_box = self.build_box_compartment(self.box_section.central_section)
        self.rear_box = self.build_box_compartment(self.box_section.rear_section)
        self.shell = self.build_external_shell()
        self.linearized_body = self.build_linearized_body()

    def build_airfoil_body(self): ## causes irregularities in external shape
        body = (
            cq.Workplane()
            .placeSketch(self.airfoil_section.sketch.copy())
            .extrude(self.length)
        )

        return body

    def build_external_shell(self): ## causes irregularities in external shape
        shell = (
            cq.Workplane()
            .placeSketch(self.airfoil_section.sketch.copy())
            .extrude(self.length)
            .faces(">Z or <Z")
            .shell(self.shell_thickness)
        )

        return shell
       
    def build_linearized_body(self):
        linearized_body = (
            cq.Workplane()
            .sketch()
            .polygon(self.airfoil_section.profile_points)
            .finalize()
            .extrude(self.length)
        )
        
        return linearized_body    
          
    def build_inner_body(self):
        profile_points = self.airfoil_section.profile_points
        if not self.airfoil_section.sharp_tail:
            profile_points = profile_points + [profile_points[0]]
            
        inner_section = (
            cq.Sketch()
            .polygon(profile_points)
            .wires()
            .offset(-self.shell_thickness, "i")
        )
        
        inner_body = (
            cq.Workplane()
            .placeSketch(inner_section)
            .extrude(self.length)
        )
        
        return inner_body
    
    def build_foam_body(self): ## internal foam body - medium between shall and box
        foam = (
            cq.Workplane()
            .add(self.inner_body)
            .center(0,0)
            .placeSketch(
                self.box_section.front_section,
                self.box_section.central_section,
                self.box_section.rear_section
            )
            .cutThruAll()
        )
        
        return foam
      
    def build_box_compartment(self, box_section):
        box = (
            cq.Workplane()
            .placeSketch(box_section)
            .extrude(self.length)
            .faces(">Z")
            .wires()
            .toPending()
            .offset2D(-self.box_thickness)
            .cutThruAll()
            .clean()
        )
        
        return box
        
    def assign_materials(self, materials: dict):
        if "box" in materials.keys():
            self.box_material = materials["box"]
        
        if "shell" in materials.keys():
            self.shell_material = materials["shell"]
            
        if "foam" in materials.keys():
            self.foam_material = materials["foam"]
            
    def get_max_bend_displacement(self, load):
        Ixx, Iyy, Ixy = self.box_section.inertia_moments
        nu_max = (load * (self.length*1e-3)**4)/(8 * self.box_material.tensile_modulus * Ixx * (1e-3)**4)
        
        return nu_max * 1e3 ## [mm]

    def get_bend_displacement(self, load, dist):
        Ixx, Iyy, Ixy = self.box_section.inertia_moments
        l = self.length*1e-3
        d = dist*1e-3
        E = self.box_material.tensile_modulus
        nu = load * (2*l*d**3 - 3*(l**2)*(d**2) - (d**4)/2) / (12 * E * Ixx * (1e-3)**4)
        
        return -nu * 1e3 ## [mm]

    
    ## Mass properties
    def get_components_mass(self, components, material):
        volume = 0
        
        for c in components:
            volume += c.solids().vals()[0].Volume()
            
        return volume * (1e-3)**3 * material.density
    
    def get_box_mass(self):
        return self.get_components_mass([self.front_box, self.central_box, self.rear_box], self.box_material)
        
    def get_foam_mass(self):
        return self.get_components_mass([self.foam], self.foam_material)
              
    def get_shell_mass(self):
        ## computation of BREP volume is error-prone - in some cases produces wrong or negative results
        ## thus we will use linearized shell body to estimate its mass
        
        shell_mass = self.get_components_mass([self.shell], self.shell_material)
        if shell_mass <= 0:
            print("BREP shell mass is negative. Realculating approximte value usng linearized bodies")
            body_mass = self.get_components_mass([self.linearized_body], self.shell_material)
            inner_body_mass = self.get_components_mass([self.inner_body], self.shell_material)
                  
            shell_mass = body_mass - inner_body_mass
                  
            if shell_mass <= 0:
                raise ValueError("Can not calculate shell mass correcly")     
                  
        return shell_mass
    
    def get_console_mass(self):
        box_mass = self.get_box_mass()
        foam_mass = self.get_foam_mass()
        shell_mass = self.get_shell_mass()
                  
        if shell_mass > box_mass:
            print(f"Warning: poor design. Shell mass is higher than box mass ({shell_mass} > {box_mass})")
        
        return box_mass + foam_mass + shell_mass    
      
    def compute_lift_force(self, alpha, fluid_props, compute_weight_load=True, 
                           load_factor=1, g=STANDARD_GRAVITY):
        
        console_weight = self.get_console_mass() * g * load_factor
        reynolds = self.airfoil_section.eval_reynolds(fluid_props)
        cl = self.airfoil_section.airfoil.eval_cl(alpha, reynolds)
        lift_force = fluid_props.dynamic_pressure * cl * (self.chord * self.length * 1e-6)
        
        if compute_weight_load:
            lift_force = lift_force - console_weight
            
        return lift_force, self.length / 2

    def compute_drag_force(self, alpha, fluid_props):
        reynolds = self.airfoil_section.eval_reynolds(fluid_props)
        cd = self.airfoil_section.airfoil.eval_cd(alpha, reynolds)
        drag_force = fluid_props.dynamic_pressure * cd * (self.chord * self.length * 1e-6)
        
        return drag_force

    def compute_bend_force(self, alpha, lift_force, drag_force, mass, g=STANDARD_GRAVITY):
        alpha_rad = alpha * math.pi / 180
        bend_force_up = lift_force * math.cos(alpha_rad) + drag_force * math.sin(alpha_rad)
        bend_force_down = mass * g * math.cos(alpha_rad)

        return bend_force_up - bend_force_down

    def compute_bend_coefficient(self, alpha, cl, cd):
        alpha_rad = alpha * math.pi / 180
        cb = cl * math.cos(alpha_rad) + cd * math.sin(alpha_rad)

        return cb

    def fit_length_to_required_tip_deflection(self, alpha, fluid_props, tip_delta_max=0.1,
                                              g=STANDARD_GRAVITY, load_factor=1, verbose=True, 
                                              adopt_result=True):
        """
        Find console length that ensures tip deflection equal to tip_delta_max
        """
        self.tip_delta_max = tip_delta_max

        reynolds = self.airfoil_section.eval_reynolds(fluid_props)
        alpha_rad = alpha * math.pi / 180

        cl = self.airfoil_section.airfoil.eval_cl(alpha, reynolds)
        cd = self.airfoil_section.airfoil.eval_cd(alpha, reynolds)
        cb = self.compute_bend_coefficient(alpha, cl, cd)

        specific_aerodynamic_load = fluid_props.dynamic_pressure * cb * self.chord * 1e-3

        def eval_displacement(length):
            length = float(length)
            if length < self.min_length or length > self.max_length:
                return 1e6

            self.build_geometry(length)
            mass = self.get_console_mass()
            specific_weight = mass * g * math.cos(alpha_rad) * load_factor / (length * 1e-3)
            total_specific_load = specific_aerodynamic_load - specific_weight

            tip_deflection = self.get_max_bend_displacement(total_specific_load)
            relative_deflection = tip_deflection / length
            deflection_error = abs(relative_deflection) - self.tip_delta_max
            
            if verbose:
                print(f"Console length: {length} [mm]")
                print(f"Сonsole mass: {mass} [kg]")
                print(f"Console spect ratio: {length / self.chord}")
                print(f"Lift force: {specific_aerodynamic_load * length * 1e-3} [N]")
                print(f"Absolute tip deflection: {tip_deflection} [mm]")
                print(f"Relative tip deflection: {100*relative_deflection} %")
                print(f"Box thickness: {self.box_thickness}, [mm]")
                print(f"Tip deflection error: {deflection_error}")
                print()

            return abs(deflection_error)**2

        print(f"Finding wing console length for the chord = {self.chord}")
        result = minimize_scalar(
            eval_displacement, 
            bounds=[self.min_length, self.max_length], 
            method="bounded",
            options={"xatol": 10}
        )
        if verbose:
            print("=============================================")
            print("Wing console length selection result:\n", result)
 
        solved_length = result.x
        
        if adopt_result:
            self.build_geometry(solved_length)
        
        return solved_length
    
    def set_new_chord(self, chord):
        airfoil = self.airfoil_section.airfoil
        self.airfoil_section = AirfoilSection(airfoil, chord=chord)
        self.box_section = ThreeChamberBoxedWingSection(self.airfoil_section)
        self.chord = self.box_section.chord
        self.box_thickness = self.box_section.box_thickness
        
    def fit_chord_to_required_lift_force(self, alpha, fluid_props, required_lift_force, 
            load_factor=1, g=STANDARD_GRAVITY, method="local", adopt_result=True, verbose=True):
        """
        Find proper chord length to fit required lift force 
        and ensure tip deflection is not greater than self.tip_delta_max
        """
        
        def excess_lift_force(chord):
            chord = float(chord)
            if chord < self.min_chord or chord > self.max_chord:
                return 1e6

            self.set_new_chord(chord)
            solved_length = self.fit_length_to_required_tip_deflection(
                alpha, fluid_props, load_factor=load_factor, verbose=False
            )
            
            lift_force, lift_force_arm = self.compute_lift_force(
                alpha, fluid_props, load_factor=load_factor, compute_weight_load=True
            )

            error = lift_force - required_lift_force

            if verbose:
                print(f"Current chord: {chord}, [mm]")
                print(f"Console length: {solved_length}, [mm]")
                print(f"Console excess lift force: {lift_force}, [N]")
                print(f"Required lift force: {required_lift_force}, [N]")
                print(f"Console mass: {self.get_console_mass()}, [kg]")
                print(f"Box thickness: {self.box_thickness}, [mm]")
                print(f"Absolute error: {error}, [N]")
                print()

            if method == "brentq":
                return error
            else:
                return abs(error)**2
                
        if method == "local":
            result = minimize_scalar(
                excess_lift_force, 
                bounds=[self.min_chord, self.max_chord], 
                method="bounded"
            )

            solved_chord = result.x
        
        elif method == "shgo":
            result = shgo(
                excess_lift_force, 
                bounds=[(self.min_chord, self.max_chord)], 
            )

            solved_chord = result.x[0]

        elif method == "direct":
            result = direct(
                excess_lift_force, 
                bounds=[(self.min_chord, self.max_chord)], 
            )

            solved_chord = result.x[0]

        elif method == "brentq":
            result = brentq(excess_lift_force, self.min_chord, self.max_chord)
            solved_chord = result

        if adopt_result:
            self.set_new_chord(solved_chord)
            self.fit_length_to_required_tip_deflection(
                alpha, fluid_props, load_factor=load_factor, verbose=False
            )
            
        print("==============================================")
        print("Wing console chord selection result:\n", solved_chord)
        
        return solved_chord

    def stats(self, alpha, fluid_props, load_factor=1, g=STANDARD_GRAVITY):
        lift_force, lift_force_arm = self.compute_lift_force(
            alpha, fluid_props, load_factor=load_factor, compute_weight_load=True
        )
        center_of_pressure = (0.25*self.chord, lift_force_arm)
        
        drag_force = self.compute_drag_force(alpha, fluid_props)
        
        box_mass = self.get_box_mass()
        foam_mass = self.get_foam_mass()
        shell_mass = self.get_shell_mass()
        console_mass = box_mass + foam_mass + shell_mass
        weight = console_mass * g * load_factor
        
        print(f"Length: {self.length}, [mm]")
        print(f"Chord: {self.chord}, [mm]")
        print(f"Mass: {console_mass}, [kg] (box: {box_mass}, foam: {foam_mass}, shell: {shell_mass})")
        print(f"Excess lift force: {lift_force}, [N]")
        print(f"Drag force: {drag_force}, [N]")
        print(f"Aspect ratio: {self.length / self.chord}")
        print(f"Area: {self.length*self.chord*1e-6}, [m^2]")
        print(f"Angle of attack: {alpha}, [degrees]")
        print(f"Reinforcement box thickness: {self.box_thickness}, [mm]")
        print(f"Shell thickness: {self.shell_thickness}, [mm]")
        print(f"Lift to weight ratio: {(lift_force + weight) / weight}")
        print(f"Center of aerodynamic pressure (lift): {center_of_pressure}, [mm, mm]")
        print()
