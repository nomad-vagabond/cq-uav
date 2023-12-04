import numpy as np
from scipy.optimize import direct
import cadquery as cq

import cquav.splinecloud_scipy as scsp


class Airfoil:

    def __init__(self, profile_data):
        self.name = profile_data["name"]
        self.profile = scsp.LoadSubset(profile_data["profile"])
        self.profile_data = profile_data
        self.load_curves()
   
    def __repr__(self):
        return f"Airfoil-{self.name}"
    
    def get_reynolds_interval(self, reynolds_round):
        reynolds_range = set([int(key[3:]) for key in self.profile_data["curves"]])
        reynolds_values = sorted(reynolds_range)

        if reynolds_round in reynolds_values:
            return reynolds_round, reynolds_round 

        reynolds_range.add(reynolds_round)
        re_left = sorted([r for r in reynolds_range if r < reynolds_round])
        re_right = sorted([r for r in reynolds_range if r > reynolds_round])

        re_left = reynolds_values[0] if not len(re_left) else re_left[-1]
        re_right = reynolds_values[-1] if not len(re_right) else re_right[0]

        return re_left, re_right

    def interpolated_coefficient(self, alpha, reynolds, coeff):
        reynolds_round = int(reynolds*1e-3)*1000
        re0, re1 = self.get_reynolds_interval(reynolds_round)

        c0 = float(self.curves[f"{coeff}_{re0}"].eval(alpha))
        if re0 == re1:
            return c0

        c1 = float(self.curves[f"{coeff}_{re1}"].eval(alpha))
        c = c0 + (reynolds-re0)*(c1-c0)/(re1-re0)

        return c

    def eval_cl(self, alpha, reynolds):
        return self.interpolated_coefficient(alpha, reynolds, "CL")

    def eval_cd(self, alpha, reynolds):
        return self.interpolated_coefficient(alpha, reynolds, "CD")

    def eval_cm(self, alpha, reynolds):
        return self.interpolated_coefficient(alpha, reynolds, "CM")

    def load_curves(self):
        self.curves = {}
        for curve_name, curve_uid in self.profile_data["curves"].items():
            self.curves[curve_name] = scsp.LoadSpline(curve_uid)

    def cl_to_cd(self, alpha, reynolds):
        return float(self.eval_cl(alpha, reynolds)) / float(self.eval_cd(alpha, reynolds))
    
    def cd_to_cl(self, alpha, reynolds):
        return 1 / self.cl_to_cd(alpha, reynolds)

    def alpha_optimal(self, reynolds):
        def cd_to_cl(alpha):
            q = self.cd_to_cl(alpha, reynolds)
            if np.isnan(q):
                return 1e3

            return q
        
        res = direct(cd_to_cl, bounds=[(0.0, 30.0)])
        return float(res.x[0])
        
    def alpha_max_lift(self, reynolds):
        def eval_cl_inverse(alpha):
            cl = self.eval_cl(alpha, reynolds)
            if np.isnan(cl):
                return 2
            
            return -float(cl)

        res = direct(eval_cl_inverse, bounds=[(0.0, 30.0)])
        return float(res.x[0])

    def alpha_min_drag(self, reynolds):
        def eval_cd(alpha):
            cd = self.eval_cd(alpha, reynolds)
            if np.isnan(cd):
                return 1e3

            return cd

        res = direct(eval_cd, bounds=[(0.0, 30.0)])
        return float(res.x[0])
