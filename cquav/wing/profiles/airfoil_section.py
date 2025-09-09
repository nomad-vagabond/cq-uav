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


class AirfoilSection:
    def __init__(self, airfoil, chord=1000, x_offset=0, y_offset=0, z_offset=0, twist_angle=0):
        self.airfoil = airfoil
        self.profile_curve = airfoil.profile_curve # spline curve
        self.chord = chord
        self.offset_vector = cq.Vector(x_offset, y_offset, z_offset)
        self.twist_axis = cq.Vector(0, 0, -1)  # twist around Z axis
        self.twist_angle = twist_angle
        self.sketch = self.build_sketch()

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
        k = self.profile_curve.k
        t = self.profile_curve.knots[k:-k]
        cp = np.array([self.profile_curve.coeffs_x, self.profile_curve.coeffs_y]).T * self.chord
        
        return self._build_spline(t, cp, k)

    def build_sketch(self):
        edge = cq.Edge(BRepBuilderAPI_MakeEdge(self.bspl).Edge())
        return cq.Sketch().edge(edge).close().assemble()

    @cached_property
    def location(self):
        return cq.Location(self.offset_vector, self.twist_axis, self.twist_angle)

    def get_profile_center_of_mass(self):
        profile_face = cq.Workplane().placeSketch(self.sketch.copy()).toOCC()
        return profile_face.Center()
