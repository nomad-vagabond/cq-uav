# -*- coding: utf-8 -*-
import math
import requests, json
import numpy as np
import pandas as pd
import scipy.interpolate as si
from scipy.optimize import brentq
from functools import partial


_error_msg = {

    1: 'Parameter t must be a list or an array that represents knot vector.',

    2: 'Method parameter must be one of: "interp", "smooth", "lsq".'

}


def LoadSubset(subset_id_or_url):
    url_split = subset_id_or_url.split("/")
    if len(url_split) > 1:
        url = subset_id_or_url
    else:
        subset_id = url_split[-1]
        if "sbt_" not in subset_id or len(subset_id) != 16:
            raise ValueError("Wrong subset id was specified")
        url = "https://splinecloud.com/api/subsets/{}".format(subset_id)
    
    response = requests.get(url)
    subset = json.loads(response.content)['table']
    
    return pd.DataFrame.from_dict(subset)
        
        
        

def LoadSpline(curve_id_or_url):
    url_split = curve_id_or_url.split("/")
    if len(url_split) > 1:
        url = curve_id_or_url
    else:
        curve_id = url_split[-1]
        if "spl_" not in curve_id or len(curve_id) != 16:
            raise ValueError("Wrong curve id was specified")
        url = "https://splinecloud.com/api/curves/{}".format(curve_id)

    response = requests.get(url)
    curve = json.loads(response.content)

    curve_params = curve['spline']
    t = np.array(curve_params['t'])
    c = np.array(curve_params['c'])
    w = curve_params['w']
    tcck = t, c[:, 0], c[:, 1], curve_params['k']

    return ParametricUnivariateSpline.from_tcck(tcck)
        

class PPolyInvertible(si.PPoly):
    """Piecewise polynomial with ability to evaluate inverse dependency x(y)"""

    @classmethod
    def construct_fast(cls, c, x, extrapolate=None, axis=0):
        # self = super(PPolyInvertible, cls).construct_fast(c, x, extrapolate=extrapolate, axis=axis)
        self = super(PPolyInvertible, cls).construct_fast(
            c, x, extrapolate=extrapolate)
        self.k = len(self.c) - 1
        self.powers = np.arange(self.k, -1, -1)
        self.intervals = self._form_intervals(self.x)
        # self.powers = np.arange(self.k, 0, -1)
        return self

    @classmethod
    def from_splinefunc(cls, spline):
        self = cls.from_spline(spline.tck)
        self.project_intervals(spline)
        
        return self

    def eval_oninterval(self, n, numpoints=50):
        coeffs = self.c.T[n + self.k]
        tbreak = self.x[n + self.k]

        a = self.intervals[n][0]
        b = self.intervals[n][1]

        tpoints = np.linspace(a, b, numpoints)
        ppoints = np.zeros(len(tpoints))
        i = 0
        for t in tpoints:
            ppoints[i] = self.eval_poly(t, coeffs, tbreak)
            i += 1
        return ppoints

    def _get_poly(self, t, n, xvalue=0):
        coeffs = self.c.T[n + self.k]
        tbreak = self.x[n + self.k]
        poly = self.eval_poly(t, coeffs, tbreak)
        error = poly - xvalue

        if abs(error) < 1e-12:
            return 0.0

        return error

    def eval_poly(self, t, coeffs, tbreak):
        # poly = coeffs[0]*(t - tbreak)**3  + coeffs[1]*(t - tbreak)**2 + coeffs[2]*(t - tbreak) + coeffs[3]
        poly = 0
        for c, p in zip(coeffs, self.powers):
            poly += c*(t - tbreak)**p
        
        return poly

    def _get_interval(self, coord, intervals):
        i = 0
        for interval in intervals:
            if coord >= interval[0] and coord <= interval[1]:
                return i
            else:
                i += 1
        
        return None

    def _form_intervals(self, breaks):
        # n = len(breaks) - 1
        n = len(breaks) - 2*self.k - 1
        intervals = np.zeros((n, 2))
        i = self.k
        for interval in intervals:
            interval[0], interval[1] = breaks[i], breaks[i + 1]
            i += 1
        
        return intervals

    def project_intervals(self, sf):
        breaks = sf(self.x)
        self.pintervals = self._form_intervals(breaks)

    def _check_monotonous(self, intervals):
        check = True
        for interval in intervals:
            if interval[1] < interval[0]:
                check = False
                break
        
        return check

    def evalinv(self, x):
        pinterval = self._get_interval(x, self.pintervals)
        if pinterval is not None:
            interval = self.intervals[pinterval]
            t = brentq(partial(self._get_poly, n=pinterval, xvalue=x), 
                       interval[0], interval[1])
            return t
        else:
            return None


class ParametricUnivariateSpline(object):
    """
    One-dimensional parametric spline fit to a given set of data points.

    Fits a spline x, y = spl(t) of degree `k` to the provided `x`, `y` data.

    If fitting method is set to "interp" spline will interpolate 
    through all data points.

    If fitting method is set to "smooth" then normalized smoothing 
    factor sn will be used to choose the number of knots.
    Regular smoothing factor s used by underlying spline functions is evaluated as:
    s = sn*sum((y_data[i])**2)

    If fitting method is set to "lsq" and internal knot vector t is not specified 
    then uniform knot vector of length nk will be used to create least squares 
    spline approximation.

    """
    def __init__(self, x_data, y_data, t=None, method="interp", 
                 sn=None, k=3, w=None, nk=3, bbox=[None]*2):
        
        self.x_data = x_data
        self.y_data = y_data
        self.k = k
        self.data_len = len(self.x_data)
        self.xmax, self.xmin = max(x_data), min(x_data)
        self.nk = nk

        if w is None:
            w_ = np.ones(self.data_len)
        else:
            w_ = np.array(w)
        ## sum((w[i] * (y[i]-spl(x[i])))**2, axis=0)
        # sscale = sum( (y_data[i])**2 for i in range(self.data_len))
        if sn is not None:
            spl_smax = si.LSQUnivariateSpline(
                x_data, y_data, [], k=k, w=w, bbox=bbox)
            s_data = [spl_smax(d) for d in x_data]
            smax = sum((w_[i]*(y_data[i] - s_data[i]))**2 
                for i in range(self.data_len))
            s = sn*smax
        else:
            s = None

        if method == "interp":
            # self.splinefunc = si.InterpolatedUnivariateSpline(x_data, y_data)
            self.splinefunc = si.UnivariateSpline(
                x_data, y_data, k=k , s=0.0, w=w, bbox=bbox)
        elif method == "smooth":
            self.splinefunc = si.UnivariateSpline(
                x_data, y_data, k=k , s=s, w=w, bbox=bbox)
        elif method == "lsq":
            if t is None:
                knots = self._uniform_knotvector(self.nk)
            elif len(t) > 0:
                knots = t
            else:
                raise ValueError(_error_msg[0])
            self.splinefunc = si.LSQUnivariateSpline(
                x_data, y_data, knots, k=k, w=w, bbox=bbox)
        else:
            raise ValueError(_error_msg[1])

        knots = self.splinefunc.get_knots()
        self.knots = self._form_knotvector(knots)
        self.knots_norm = self._normalize_knotvector( d=self.data_len)
        # self.knots_norm = self._normalize_knotvector(self.knots) #newfix

        self.coeffs_x = self._get_controlpoints(self.knots)
        self.coeffs_y = self.splinefunc.get_coeffs()
        self.coeffs_t = self._get_controlpoints(self.knots_norm)

        self._build_splines(self.coeffs_x, self.coeffs_y)
        self._get_ppolyrep()


    @classmethod
    def from_tcck(cls, tcck):
        """Construct a parametric spline object from given tcck"""
        self = cls.__new__(cls)
        t, cx, cy, k = tcck
        self.k = k
        self.knots = t
        self.knots_norm = self._normalize_knotvector()
        self.coeffs_x = cx
        self.coeffs_y = cy
        self._build_splines(self.coeffs_x, self.coeffs_y)
        self._get_ppolyrep()
        
        return self

    def __call__(self, tpoints):
        x_points = self.spline_x(tpoints)
        y_points = self.spline_y(tpoints)
        
        return x_points, y_points

    def eval(self, x):
        if hasattr(x, '__iter__'):
            t = np.array([self.spline_x.ppoly.evalinv(xi) for xi in x])
            return self.spline_y.ppoly(t)
        
        else:
            t = self.spline_x.ppoly.evalinv(x)
            return self.spline_y.ppoly(t)

    def get_polypoints(self, n):
        xpoints = self.spline_x.ppoly.eval_oninterval(n)
        ypoints = self.spline_y.ppoly.eval_oninterval(n)
        
        return xpoints, ypoints

    def _get_ppolyrep(self):
        self.spline_x.ppoly = PPolyInvertible.from_splinefunc(self.spline_x)
        self.spline_y.ppoly = PPolyInvertible.from_splinefunc(self.spline_y)

    def polyrep(self, tpoints):
        return self.spline_x.ppoly(tpoints), self.spline_y.ppoly(tpoints)

    def _build_splines(self, coeffs_x, coeffs_y):
        tck_x = self.knots_norm, coeffs_x, self.k
        tck_y = self.knots_norm, coeffs_y, self.k
        self.spline_x = si.UnivariateSpline._from_tck(tck_x)
        self.spline_y = si.UnivariateSpline._from_tck(tck_y)
        self.spline_x.tck = tck_x
        self.spline_y.tck = tck_y

    def _form_knotvector(self, knots):
        knots_full = np.concatenate(
            ([knots[0]]*self.k, knots, [knots[-1]]*self.k ))
        return knots_full

    def _normalize_knotvector(self, knots=None, d=1.0):
        if knots is None: knots = self.knots
        num_knots = len(knots)
        ka = (knots[-1] - knots[0]) / d
        knots_norm = np.empty(num_knots)
        for i in range(num_knots):
            knots_norm[i] = d - ((knots[-1] - knots[i])) / ka
        
        return knots_norm

    def _get_controlpoints(self, knots):
        n = len(knots) - 1 - self.k
        cpoints = np.empty(n)
        for i in range(n):
            tsum = 0
            for j in range(1, self.k + 1):
                tsum += knots[i + j]
            cpoints[i] = tsum/float(self.k)
        
        return cpoints

    def _uniform_knotvector(self, nk):
        if nk == 0:
            return []
        elif nk == 1:
            return [(self.xmax - self.xmin) / 2.0 + self.xmin]
        else:
            knot_offset = float(self.xmax - self.xmin) / nk
            # ks = self.xmin + knotdist
            # ke = self.xmax - knotdist
            knots = np.linspace(self.xmin, self.xmax, nk+2)
            knots = knots[1:-1]
            # knots = np.linspace(knot_offset, self.xmax-knot_offset, nk-2)
            return knots
