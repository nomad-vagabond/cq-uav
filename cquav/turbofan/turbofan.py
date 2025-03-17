import json
from functools import cache, cached_property
from pathlib import Path

from OCP.gp import gp_Pnt
from OCP.TColgp import TColgp_Array1OfPnt
from OCP.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCP.Geom import Geom_BSplineCurve
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge

import cadquery as cq
import numpy as np
import splinecloud_scipy as sc


AIRFOILS_COLLECTION_JSON_FILE = Path(__file__).parent.parent/"wing"/"airfoil"/"airfoils_collection.json"

@cache
def get_refined_airfoils_collection(data = None, data_file = AIRFOILS_COLLECTION_JSON_FILE):
    if not data:
        with open(data_file) as file:
            data = json.load(file)

    curves_dict = {}
    for category, airfoils in data.items():
        for airfoil_name, airfoil_data in airfoils.items():
            if "profile_curve" in airfoil_data:
                if profile_curve := airfoil_data["profile_curve"]:
                    curves_dict[airfoil_name] = profile_curve[0]

    return curves_dict


class AirfoilSection:
    def __init__(self, profile_curve, chord, offset, twist_angle):
        self.profile_curve = profile_curve
        self.chord = chord
        self.offset = cq.Vector(0, 0, offset)
        self.twist_axis = cq.Vector(0, 0, 1)  # twist around Z axis
        self.twist_angle = twist_angle

    @cached_property
    def spline(self):
        if isinstance(self.profile_curve, sc.ParametricUnivariateSpline):
            return self.profile_curve
        return sc.load_spline(self.profile_curve)

    @staticmethod
    def _build_spline(knot_vector, control_points, degree):
        knot_mults = list(map(lambda k: degree + 1 if k in (0, 1) else 1, knot_vector))

        poles = TColgp_Array1OfPnt(1, len(control_points))
        for i, cp in enumerate(control_points):
            poles.SetValue(i + 1, gp_Pnt(cp[0], cp[1], 0))

        knots = TColStd_Array1OfReal(1, len(knot_vector))
        for i, k in enumerate(knot_vector):
            knots.SetValue(i + 1, k)

        mults = TColStd_Array1OfInteger(1, len(knot_mults))
        for i, m in enumerate(knot_mults):
            mults.SetValue(i + 1, m)

        return Geom_BSplineCurve(poles, knots, mults, degree)

    @cached_property
    def bspl(self):
        k = self.spline.k
        t = self.spline.knots[k:-k]
        cp = np.array([self.spline.coeffs_x, self.spline.coeffs_y]).T * self.chord
        return self._build_spline(t, cp, k)

    def build_sketch(self):
        edge = cq.Edge(BRepBuilderAPI_MakeEdge(self.bspl).Edge())
        return cq.Sketch().edge(edge).close().assemble()

    @cached_property
    def location(self):
        return cq.Location(self.offset, self.twist_axis, self.twist_angle)


class Turbofan:
    def __init__(self, sections: list, vanes_count: int, center_hole_diameter: int, hub_diameter: int):
        self.sections = sections
        self.vanes_count = vanes_count
        self.center_hole_diameter = center_hole_diameter
        self.hub_diameter = hub_diameter

    @cached_property
    def _vane_trailing_edge_max_coordinate(self):
        edges = self.sections[0].build_sketch().moved(self.sections[0].location).edges()
        vertices = [v for e in edges for v in e.Vertices()]
        vertex_coords = [v.toTuple() for v in vertices]
        trailing_edge_tip_coords = max(vertex_coords, key=lambda v: v[0])
        return max(trailing_edge_tip_coords)

    @cache
    def build_vane(self):
        if not self.sections:
            raise ValueError(f"The number of sections has to be > `0`, but you provided `{len(self.sections)}`")
        return cq.Workplane().placeSketch(
            *(section.build_sketch().moved(section.location) for section in self.sections)
        ).loft()

    @cached_property
    def hub_height(self):
        bbox = self.build_vane().val().BoundingBox()
        return bbox.ymax - bbox.ymin

    @cache
    def build_turbofan(self):
        vane_offset = ((self.hub_diameter / 2) ** 2 - self._vane_trailing_edge_max_coordinate ** 2) ** 0.5
        vanes = [
            self.build_vane().translate((0, 0, vane_offset)).rotate((0, 0, 0), (0, 1, 0), i * 360 / self.vanes_count)
            for i in range(self.vanes_count)]

        hub = cq.Workplane("XZ").circle(self.hub_diameter / 2).extrude(self.hub_height)
        hub = hub.faces("XZ").workplane().circle(self.hub_diameter / 2).extrude(-self.hub_height / 2)
        hub = hub.faces("+Y").edges().fillet(self.hub_height * 0.25)

        for vane in vanes:
            hub = hub.union(vane)

        return hub.faces("XZ").circle(self.center_hole_diameter / 2).cutThruAll()
