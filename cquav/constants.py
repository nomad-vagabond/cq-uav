import os
import math

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

