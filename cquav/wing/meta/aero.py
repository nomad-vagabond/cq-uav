from math import sin, cos, tan, sqrt, pi, exp, acos
from copy import deepcopy
from time import time
import re, os, csv, json
from collections import namedtuple

import pandas as pd
import numpy as np
from scipy.optimize import direct, minimize_scalar
import cadquery as cq
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
import splinecloud_scipy as scloud
import matplotlib.pyplot as plt


STANDARD_GRAVITY = 9.80665
## Physical constants

AIR_DENSITY = 1.225 #[kg/m^3]
STANDARD_GRAVITY = 9.80665 ## [m/s]
AIR_KINEMATIC_VISCOSITY_SEA_LEVEL = 1.460*1e-5 ## [m^2/s]

DELTA_MAX = 0.1 ## relative wing tip displacement
G = 9.81 ## [m/s**2]
LOAD_FACTOR = 1.0

## Wing Geometry

CHORD_MIN = 50
CHORD_MAX = 1000
CHORD_DEFAULT = 260

SPAN_MIN = 50
SPAN_MAX = 5000
SPAN_DEFAULT = 900

VELOCITY = 42 #[m/s]
SPEED_OF_SOUND = 340 #[m/s]
MACH = VELOCITY/SPEED_OF_SOUND

DYN_AIR_PRESSURE = 0.5 * AIR_DENSITY * VELOCITY**2 ## dynamic pressure

is_int = lambda s: s.isdigit() or (s[0] == '-' and s[1:].isdigit())


XFOIL_LIFT_OVERESTIMATE_ERROR = 0.13
XFOIL_DRAG_UNDERESTIMATE_ERROR = 0.14


def inverse(v, scloud_func):
    res = scloud_func.eval(v)
    if res is np.nan:
        return 1e6
    
    return -res


def neighbor_ar(ar_range, ar):
    ar_range_ = deepcopy(ar_range)
    ar_range_.append(ar)
    ar_range_ = sorted(ar_range_)
    ind = ar_range_.index(ar)
    if ind == 0:
        return ar_range[0], ar_range[1]
    elif ind >= len(ar_range_)-1:
        return ar_range[-1], ar_range[-1]

    return ar_range_[ind-1], ar_range_[ind+1]


class NormalizedCoeffDistribution:
    """
    Chord scaled normalized aerodynamic coefficient distribution over the normalized semi-span
    
    Parameters
    ccl_collection - nested dict with curve ids for the coefficient distributions
    
    """
    
    def __init__(self, ccl_collection):
        self.curves = self.load_curves(ccl_collection)
        self.ar_range = sorted(ccl_collection.keys())
        self.ccls_by_ar_ns = self.get_linear_sweep_by_ar()
    
    def load_curves(self, collection):
        return {ar: scloud.load_spline(curve_id) for ar, curve_id in collection.items()}
    
    def reduce_ar_range(self):
        self.ar_range = [self.ar_range[0], self.ar_range[-1]]
    
    def get_linear_sweep_by_ar(self):
        def interp_ar(ar, ns):
            ar_left, ar_right = neighbor_ar(self.ar_range, ar)
            
            ccls_left = self.curves[ar_left].eval(ns, extrapolate=True)
            ccls_right = self.curves[ar_right].eval(ns, extrapolate=True)
            ccls = (ar-ar_left)*(ccls_right-ccls_left)/(ar_right-ar_left) + ccls_left

            return ccls
        
        return interp_ar
    
    def surfplot(self, ar_num=40, ns_num=40, ar_min=None, ar_max=None):
        if ar_min is None:
            ar_min = min(self.ar_range)

        if ar_max is None:
            ar_max = max(self.ar_range)

        ar_range = np.linspace(ar_min, ar_max, ar_num)
        ns_range = np.linspace(0, 1, ns_num)
        ar_vals, ns_vals = np.meshgrid(ar_range, ns_range)
        ccls_vals = np.array([[self.ccls_by_ar_ns(ar_, ns_) for ar_ in ar_range] for ns_ in ns_range])
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(ar_vals, ns_vals, ccls_vals, cmap='viridis')
        ax.set_title('CCLS(AR, NS)')
        ax.set_xlabel('Aspect Ratio')
        ax.set_ylabel('Normalized Span')
        ax.set_zlabel('CCLS')
        ax.view_init(elev=15., azim=45.0)
        plt.show() 


class WingSpanEfficiency:
    """
    Span Efficiency Factor
    
    Parameters
    se_collection - dict with curve ids for span efficiency factor functions
    
    TODO: write tests
    """
    
    def __init__(self, se_collection):
        self.curves = self.load_curves(se_collection)
        self.ar_range = sorted(self.curves.keys())
        self.eval = self.get_linear_sweep_by_ar()
    
    def load_curves(self, collection):
        curves = {}
        for ar, curve_id in collection.items():
            curves[ar] = scloud.load_spline(curve_id)
        
        return curves
    
    def reduce_ar_range(self):
        self.ar_range = [self.ar_range[0], self.ar_range[-1]]
    
    def get_linear_sweep_by_ar(self):
        def interp_ar(ar, tr):
            ar_left, ar_right = neighbor_ar(self.ar_range, ar)
            
            delta_left = self.curves[ar_left].eval(tr, extrapolate=True)
            delta_right = self.curves[ar_right].eval(tr, extrapolate=True)
            delta = (ar-ar_left)*(delta_right-delta_left)/(ar_right-ar_left) + delta_left

            return delta
        
        return interp_ar

    def surfplot(self, ar_num=40, tr_num=40):
        ar_range = np.linspace(4, 10, ar_num)
        tr_range = np.linspace(0, 1, tr_num)
        ar_vals, tr_vals = np.meshgrid(ar_range, tr_range)
        delta_vals = np.array([[self.eval(ar_, tr_) for ar_ in ar_range] for tr_ in tr_range])
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(ar_vals, tr_vals, delta_vals, cmap='viridis')
        ax.set_title('Span Eficiency(AR, TR)')
        ax.set_xlabel('Aspect Ratio')
        ax.set_ylabel('Taper Ratio')
        ax.set_zlabel('delta')
        ax.view_init(elev=15., azim=45.0)
        plt.show()  
