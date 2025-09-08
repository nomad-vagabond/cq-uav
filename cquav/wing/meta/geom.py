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

from ...constants import *


class Wing:
    
    k_alpha_curve = "spl_yEjpFluF2H1i"
    dk_alpha_curve = "spl_EY51y783x087"
    
    def __init__(self, ccls_surface, span_efficiency, sweep_angle, aspect_ratio, span, airfoil, 
                 taper_ratio=0.75, fuselage_diam=0, fuselage_nose_offset=0, dihedral=0, velocity=VELOCITY):
        self.ccls_surface = ccls_surface
        self.se = span_efficiency
        self._load_curves()
        
        self.airfoil = airfoil
        self.ar = aspect_ratio
        self.tr = taper_ratio
        self.sweep_angle = sweep_angle # [degrees]
        self.span = span # [mm]
        self.dihedral = dihedral
        
        self.fuselage_diam = fuselage_diam
        self.fuselage_nose_offset = fuselage_nose_offset
        
        self.update_dimensions()
        self.update_aero_params()
        self.reynolds = self.eval_reynolds(velocity)
        
    def update_dimensions(self):
        self.area = 1e-6 * self.span**2 / self.ar # [m^2]
        self.root_chord = 2 * self.area*1e6 / ((1+self.tr) * self.span) # root chord
        self.tip_chord = self.root_chord * self.tr

    def update_aero_params(self):
        self.CCLS_avg = self.average_scaled_lift_coefficient()
        self.CCLS_avg_console = self.average_scaled_lift_coefficient(offset=self.fuselage_diam/self.span)
        self.rel_aero_midspan = self.eval_relative_aero_midspan()

    def _load_curves(self):
        self.k_alpha = scloud.load_spline(self.k_alpha_curve)
        self.dk_alpha = scloud.load_spline(self.dk_alpha_curve)
    
    def _alpha_actual(self, alpha, delta=0):
        if self.dihedral == 0:
            return alpha + delta
        
        vc = self.actual_velocity_coeff(alpha)
        alpha_ = acos(cos(alpha*pi/180)/vc)*180/pi
        
        return alpha_ + delta
    
    def actual_velocity_coeff(self, alpha):
        return sqrt(1 - (sin(alpha*pi/180)*sin(self.dihedral*pi/180))**2)
        
    def eval_reynolds(self, velocity):
        reynolds = velocity * self.root_chord*1e-3 / AIR_KINEMATIC_VISCOSITY_SEA_LEVEL
        
        return reynolds
    
    def lift_force(self, alpha, delta=0, velocity=VELOCITY):
        ## calculate net lift force
        alpha_ = self._alpha_actual(alpha, delta)
        cl = self.airfoil.eval_cl(alpha_, self.reynolds)
        c_dihedral = cos(self.dihedral * pi/180)

        dyn_airpressure = 0.5 * AIR_DENSITY * velocity**2 ## dynamic pressure
        rect_wing_area = self.root_chord * self.span * 1e-6

        return dyn_airpressure * rect_wing_area * cl * self.CCLS_avg * c_dihedral
    
    def mean_lift_coefficient(self, alpha):
        cl = self.airfoil.eval_cl(alpha, self.reynolds)
        c_dihedral = cos(self.dihedral * pi/180)
        
        return cl * self.CCLS_avg * c_dihedral
    
    def cl_slope(self, alpha1=0, alpha2=5, alpha_unit="degree"):
        cl1 = self.mean_lift_coefficient(alpha1)
        cl2 = self.mean_lift_coefficient(alpha2)

        scale = 1 if alpha_unit == "degree" else 180/pi
        
        return scale * (cl2-cl1) / (alpha2-alpha1)
        
    def drag_force(self, alpha, delta=0, velocity=VELOCITY):
        ## calculate net grag force
        alpha_ = self._alpha_actual(alpha, delta)
        cd0 = self.airfoil.eval_cd(alpha_, self.reynolds)
        cdi = self.induced_drag_coefficient(alpha_)

        dyn_airpressure = 0.5 * AIR_DENSITY * velocity**2 ## dynamic pressure
        
        return dyn_airpressure * self.area * (cd0+cdi)
      
    def induced_drag_coefficient(self, alpha):
        delta = self.se.eval(self.ar, self.tr)
        cl = self.CCLS_avg*self.airfoil.eval_cl(alpha, self.reynolds)
        
        return (cl**2)*(1+delta)/(pi*self.ar)
        
    def eval_relative_aero_midspan(self, n=100):
        """Calculates relative span coordinate where the total console lift force is applied"""
        ns_vals = np.linspace(0, 1, n) ## normalized half span coords
        dy = 1/len(ns_vals)

        ccls_area = 0
        ccls_mom = 0
        for ns in ns_vals:
            ccls_i = self.ccls_surface(self.ar, ns+0.5*dy)
            ccls_area += ccls_i * dy
            ccls_mom += ccls_i * dy * (ns + 0.5*dy)

        rel_midspan = ccls_mom / ccls_area
        
        return rel_midspan
    
    def average_scaled_lift_coefficient(self, offset=0, n=100):
        half_area_norm = (1+self.tr)/2
        ns_vals = np.linspace(offset, 1, n) ## normalized half span coords
        dy = 1/len(ns_vals)
        CCLS_area = sum([self.ccls_surface(self.ar, ns+0.5*dy)*dy for ns in ns_vals])
        
        return CCLS_area/1.0 # average CCLS
        
    def console_lift_force_position(self):
        midspan = 0.5 * self.span * self.rel_aero_midspan
        
        return tan(self.sweep_angle * pi/180)*midspan + 0.25*self.aero_mid_chord, midspan
    
    @property
    def aero_mid_chord(self):
        ## mid aero chord computed based on the lift distribution
        tr_midspan = 1 - self.rel_aero_midspan*(1-self.tr)
        
        return self.root_chord * tr_midspan
    
    @property
    def mean_aero_chord(self):
        mac = self.root_chord*(2/3)*(1+self.tr+self.tr**2)/(1+self.tr)
        
        return mac    
        
    @property
    def fuselage_chord(self):
        rel_zerospan = self.fuselage_diam/self.span
        tr_zerospan = 1 - rel_zerospan*(1-self.tr)
        
        return self.root_chord * tr_zerospan
      
    @property    
    def console_area(self):
        fuselage_chord = self.fuselage_chord
        
        return 0.5*(self.span-self.fuselage_diam)*0.5*(fuselage_chord + self.tip_chord)

    def alpha_max_lift(self):
        def eval_cl_inverse(alpha):
            cl = self.mean_lift_coefficient(alpha)
            if np.isnan(cl):
                return 2
            
            return -float(cl)

        res = direct(eval_cl_inverse, bounds=[(0.0, 30.0)])
        return float(res.x[0])

    def fuselage_interference_coeff(self, verbose=False):
        if self.fuselage_diam == 0:
            return 1
        
        df = self.fuselage_diam/self.span
        k_alpha = self.k_alpha.eval(df)
        dk_alpha = self.dk_alpha.eval(df)

        taper_coeff = (1+3*df - df*(1-df)*self.tip_chord/self.fuselage_chord)/(1+df)**2
        
        offset_ratio = self.fuselage_nose_offset/self.fuselage_diam
        reynolds = VELOCITY * self.fuselage_nose_offset / AIR_KINEMATIC_VISCOSITY_SEA_LEVEL
        displacement_thickness =(0.093/(reynolds)**(1/5))*offset_ratio*(1+0.4*MACH+0.147*MACH**2-0.006*MACH**3)
        eta = 1/self.tr
        boundary_layer_coeff = 1-df*(1+df**2)*(df+eta*(1+eta)-1)*(displacement_thickness/(1+eta**2))
        
        offset_coeff = 0.6 + 0.4*(1-exp(-0.5*offset_ratio))
        
        if verbose:
            print("k_alpha:", k_alpha)
            print("dk_alpha:", dk_alpha)
            print("taper_coeff:", taper_coeff)
            print("boundary_layer_coeff:", boundary_layer_coeff)
            print("offset_coeff:", offset_coeff)
        
        kaa = (k_alpha+dk_alpha)*taper_coeff*boundary_layer_coeff*offset_coeff
        
        return kaa
        
    def plot_wing_geom(self):
        tangam = tan(self.sweep_angle * pi/180)
        semispan = self.span/2
        aero_midspan = 0.5 * self.span * self.rel_aero_midspan
               
        console_poly = [
            (0, 0),
            (tangam*semispan, semispan),
            (tangam*semispan + self.tip_chord, semispan),
            (self.root_chord, 0),
            (0, 0)
        ]
        
        midaxis = [
            (self.root_chord/2, 0), 
            (tangam*semispan + self.tip_chord/2, semispan)
        ]
        
        aero_mid_chord_line = [
            (tangam*aero_midspan, aero_midspan),
            (tangam*aero_midspan + self.aero_mid_chord, aero_midspan)
        ]
        
        x_poly, y_poly = zip(*console_poly)
        x_midax, y_midax = zip(*midaxis)
        x_midchord, y_midchord = zip(*aero_mid_chord_line)
        xac, yac = self.console_lift_force_position()
        
        fig = plt.figure()
        plt.plot(x_poly, y_poly)
        plt.plot(x_midax, y_midax, ls='--')
        plt.plot(x_midchord, y_midchord, ls='--')
        plt.plot(xac, yac, 'go')
        plt.axis('equal')
        plt.show()

