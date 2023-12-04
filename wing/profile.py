from functools import cached_property

import numpy as np
from scipy.interpolate import splprep, UnivariateSpline
import cadquery as cq
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCP.gp import gp_Pnt
from OCP.TColgp import TColgp_Array1OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.Geom import Geom_BSplineCurve


def build_occ_spline(knot_vector, control_points, degree):
    ## assume there are no repeating knots inside 
    knot_mults = list(map(lambda k: degree+1 if k in (0,1) else 1, knot_vector))

    poles = TColgp_Array1OfPnt(1, len(control_points))
    for i, cp in enumerate(control_points):
        poles.SetValue(i+1, gp_Pnt(cp[0], cp[1], 0))
        
    knots = TColStd_Array1OfReal(1, len(knot_vector))
    for i, k in enumerate(knot_vector):
        knots.SetValue(i+1, k)
        
    mults = TColStd_Array1OfInteger(1, len(knot_mults))
    for i, m in enumerate(knot_mults):
        mults.SetValue(i+1, m)
        
    bspl = Geom_BSplineCurve(poles, knots, mults, degree)
    
    return bspl     


class AirfoilSection:
    def __init__(self, airfoil, chord=1000, thicken_sharp_tail=True):
    
        self.airfoil = airfoil
        self.chord = chord
        self._setup_profile(thicken_sharp_tail)
        self.sketch = self.build_sketch()
        self.cm = self.get_profile_center_of_mass()
    
    def _setup_profile(self, thicken_sharp_tail):
        profile_points = self.airfoil.profile.to_numpy()
        self._patch_tail(profile_points, thicken_sharp_tail)

        t, c, k = self._bspl_profile_approx(profile_points)
        self.profile_tck = t, c*self.chord, k

        self.profile_points = [tuple(p) for p in profile_points*self.chord]
        self._setup_profile_curves()

    def _patch_tail(self, profile_points, thicken_sharp_tail):
        """
        Patch tail points to avoid malformed geometry construction
        and enable shell creation in CadQuery
        """
        thicken_ratio = 0.3 # share of neighboring point offsets

        self.sharp_tail = (
            profile_points[0][0] == profile_points[-1][0]
            and 
            profile_points[0][1] == profile_points[-1][1]
        )
        
        ## height at the first point pair after the tail point
        h1 = profile_points[1][1] - profile_points[-2][1] 
        if self.sharp_tail and thicken_sharp_tail:
            h0 = min(h1 * thicken_ratio, 0.001)

            profile_points[0][1] = profile_points[0][1] + h0/2
            profile_points[-1][1] = profile_points[-1][1] - h0/2

            self.sharp_tail = False

        # ## calculate x-offset of the first point pair after the tail point
        # v1 = np.array(profile_points[1]) - np.array(profile_points[0])
        # v2 = np.array(profile_points[-2]) - np.array(profile_points[-1])
        # tail_length = min(np.linalg.norm(v1), np.linalg.norm(v2))  
        # short_tail = tail_length < 0.01
        # thin_tail = h1 < 0.05

        # if short_tail or thin_tail: 
        #     ## reduce point density at the tail region
        #     profile_points.pop(1)
        #     profile_points.pop(-2)

    def _setup_profile_curves(self):
        """
        Inerpolate upper and butom profile countours with spline functions
        to enable geometry analysis (these are not used in geometry consruction).
        """
        profile_points = np.array(self.profile_points)

        nose_pnt_ind = None
        min_x_val = np.min(profile_points[:, 0])
        for i, p in enumerate(profile_points):
            if p[0] == min_x_val:
                nose_pnt_ind = i
                break

        upper_points = np.flip(profile_points[:nose_pnt_ind+1], axis=0)
        x_top, y_top = upper_points[:, 0], upper_points[:, 1]
        self.spl_top = UnivariateSpline(x_top, y_top, s=0.0)

        bottom_points = np.array(profile_points[nose_pnt_ind:])
        x_bot, y_bot = bottom_points[:, 0], bottom_points[:, 1]
        self.spl_bottom = UnivariateSpline(x_bot, y_bot, s=0.0)

    def _bspl_profile_approx(self, profile_points, s=2e-5, k=3):
        """
        B-Spline approximation of the discrete profile data, represented in Selig format
        
        Paameters:
        'pp' - array of [x,y] coords
        's' - smoothing factor
        'k' - spline degree
        
        Returns:
        't', 'cx',  - B-Spline represntation, instances of the SplineCloud ParametricUnivariateSpline
        """
        xvals = [p[0] for p in profile_points]
        yvals = [p[1] for p in profile_points]
        
        ## set weights to critical profile points
        weight_mod = lambda x: 20 if x==0.0 else (5 if x==1.0 else 1)
        weights = list(map(weight_mod, xvals))
        
        tck, u = splprep([xvals, yvals], s=s, k=k, w=weights)
        t, c, k = tck
        cx, cy = c
        
        ## adjust spline start and end points to match corrsponding profile points
        cx[0], cy[0] = profile_points[0]
        cx[-1], cy[-1] = profile_points[-1]

        print("number of profile data points:", len(profile_points))
        print("number of brep control points:", len(cx))
        
        return t[k:-k], np.array(list(zip(cx, cy))), k

    def build_sketch(self, smooth_profile=True):
        if smooth_profile:
            bspl_occ = build_occ_spline(*self.profile_tck)
            profile_edge = cq.Edge(BRepBuilderAPI_MakeEdge(bspl_occ).Edge())
            sketch = cq.Sketch().edge(profile_edge)
        else:
            sketch = cq.Sketch().spline(self.profile_points)

        if self.sharp_tail:
            return sketch.assemble()
        else:
            return sketch.close().assemble()
        
    def get_profile_center_of_mass(self):
        profile_face = cq.Workplane().placeSketch(self.sketch.copy()).toOCC()
        return profile_face.Center()
    
    @cached_property
    def profile_max_height(self):
        x_vals, y_vals = zip(*self.profile_points)
        y_max = max(y_vals)
        y_min = min(y_vals)

        return y_max - y_min

    @cached_property
    def profile_height(self):
        """
        Profile height at 25% of chord length
        """
        y_top = float(self.spl_top(0.25 * self.chord))
        y_bot = float(self.spl_bottom(0.25 * self.chord))

        return y_top - y_bot
    
    def eval_reynolds(self, fluid_props):
        reynolds = fluid_props.velocity * (self.chord*1e-3) / fluid_props.kinematic_viscosity

        return reynolds


class ThreeChamberBoxedWingSection:
    def __init__(self, airfoil_section, wall_positions=None, box_thickness=None, shell_thickness=1, thickness_tol=3):
        self.airfoil_section = airfoil_section
        self.chord = airfoil_section.chord
        
        if self.airfoil_section.sharp_tail:
            self.profile_points = list(self.airfoil_section.profile_points)
        else:
            self.profile_points = [*self.airfoil_section.profile_points, self.airfoil_section.profile_points[0]]
        
        if box_thickness:
            self.box_thickness = box_thickness
        else:
            delta_by_height = max(self.airfoil_section.profile_height * 0.08, 1)
            delta_by_chord = max(self.chord * 0.01, 1)
            self.box_thickness = round(min(delta_by_height, delta_by_chord), thickness_tol)
        
        self.min_foam_shell_thickness = min(self.airfoil_section.profile_height * 0.12, self.box_thickness * 2)
        self.min_foam_shell_thickness = max(self.min_foam_shell_thickness, 2 * shell_thickness)
        self.min_wall_height = self.box_thickness * 3
        self.wall_positions = wall_positions or [0.03, 0.25, 0.5, 0.75]
        self.adjust_wall_positions()
        
        self.front_section = self.build_box_section(self.wall_positions[:2])
        self.central_section = self.build_box_section(self.wall_positions[1:3])
        self.rear_section = self.build_box_section(self.wall_positions[2:])
        
        self.inertia_moments = self.eval_inertia_moments()
        
    def build_box_section(self, walls):
        profile_cut = self._profile_cut(self.chord * walls[0], self.chord * walls[1])
        
        return self._box_compartment_hull(profile_cut)
    
    def _profile_cut(self, x1, x2, cut_corners=True):
        xc = (x1 + x2) / 2
        section_cut = (
            cq.Sketch()
            .polygon(self.profile_points)
            .push([cq.Vector(xc, 0)])
            .rect(x2-x1, self.chord, mode="i")
        )
        
        if cut_corners:
            section_cornercut = (
                section_cut
                .edges().vertices(f">({xc}, 0, 0)").circle(self.min_foam_shell_thickness, mode="s")
                .edges().vertices(f"<({xc}, 0, 0)").circle(self.min_foam_shell_thickness, mode="s")
                .edges()
            )
            return section_cornercut
        
        else:
            return section_cut
            
    def _box_compartment_hull(self, section):
        verts = section.copy().edges('|Y').vertices().vals()
        points = [(v.Center().x, v.Center().y) for v in verts]
        
        hull = (
            cq.Sketch()
            .segment(points[0],points[1])
            .segment(points[2],points[3])
            .hull()
        )
        
        return hull
    
    def adjust_wall_positions(self):
        self._adjust_wall_position(0, increment=0.02, max_val=self.wall_positions[1]*0.9)
        self._adjust_wall_position(3, increment=-0.05, min_val=self.wall_positions[2]*1.1)
        
    def _adjust_wall_position(self, wall_ind, increment=0.05, max_val=1, min_val=0):
        correct_offset = False
        while not correct_offset:
            if self.wall_positions[wall_ind] >= max_val or self.wall_positions[wall_ind] <= min_val:
                raise ValueError("Critical wall size have reached before finding solution")
            
            wall_abs_position = self.chord * self.wall_positions[wall_ind]
            
            wall_vertices = (
                cq.Sketch()
                .polygon(self.profile_points)
                .rect(wall_abs_position*2, self.chord, mode="i")
                .edges("|Y")
                .vertices()
                .vals()
            )
            
            correct_offset = self._validate_section_wall_height(wall_vertices)
            
            if not correct_offset:
                self.wall_positions[wall_ind] += increment
                
                wall_name = {0: 'front', 1: 'second', 2: 'third', 3: 'rear'}
                print(f"new {wall_name[wall_ind]} wall position:", self.wall_positions[wall_ind])
       
    def _validate_section_wall_height(self, wall_vertices):
        wall_height = abs(wall_vertices[0].Center().y - wall_vertices[1].Center().y)
        return wall_height-self.min_foam_shell_thickness*2 >= self.min_wall_height
    
    def eval_inertia_moments(self):
        front_box = (
            cq.Workplane()
            .placeSketch(self.front_section.copy())
            .extrude(1)
            .faces(">Z or <Z")
            .shell(-self.box_thickness, kind='intersection')
        )
        
        central_box = (
            cq.Workplane()
            .placeSketch(self.central_section.copy())
            .extrude(1)
            .faces(">Z or <Z")
            .shell(-self.box_thickness, kind='intersection')
        )
        
        rear_box = (
            cq.Workplane()
            .placeSketch(self.rear_section.copy())
            .extrude(1)
            .faces(">Z or <Z")
            .shell(-self.box_thickness, kind='intersection')
        )
        
        box_section_occ = (
            cq.Workplane()
            .union(front_box)
            .union(central_box)
            .union(rear_box)
            .faces("<Z")
            .toOCC()
        )
        
        properties = GProp_GProps()
        BRepGProp.SurfaceProperties_s(box_section_occ, properties)
        matrix_of_inertia = properties.MatrixOfInertia()
        
        Ixx, Iyy, Izz = matrix_of_inertia.Value(1,1), matrix_of_inertia.Value(2,2), matrix_of_inertia.Value(3,3)
        
        return Ixx, Iyy, Izz
