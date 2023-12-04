from dataclasses import dataclass
from functools import cached_property


@dataclass
class IsotropicMaterial:
    density: float = 0.0 ## kg/m**3
    tensile_strength: float = 0.0 ## Pa
    tensile_modulus: float = 0.0 ## Pa
    shear_modulus: float = 0.0 ## Pa
    

@dataclass
class FluidProperties:
    density: float = 0.0 ## kg/m**3
    velocity: float = 0.0 ## m/s
    kinematic_viscosity: float = 0.0 ## [m^2/s]

    @cached_property
    def dynamic_pressure(self):
        return 0.5 * self.density * self.velocity**2