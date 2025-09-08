import yaml
from typing import Any, Dict, Optional, Tuple, List


class Experiment:
    def __init__(self, metadata: Dict, data: Dict):
        self.metadata = metadata  # e.g., ROOT_AIRFOIL, TIP_AIRFOIL, etc.
        self.data = data          # nested dictionary TR -> AR -> curve -> AoA

    def _retrieve(self, dictionary: Dict, key: Any, description: str) -> Any:
        """Helper to safely retrieve a value from a dictionary or raise KeyError with a message."""
        try:
            return dictionary[key]
        except KeyError:
            raise KeyError(f"{description} '{key}' not found.")

    def get_data(self, tr: float, ar: int, aoa: float, curve: str = "ccl*c/cref"):
        tr_dict = self._retrieve(self.data, f"TR={tr}", "TR")
        ar_dict = self._retrieve(tr_dict, f"AR={ar}", "AR")
        curve_dict = self._retrieve(ar_dict, curve, "Curve")
        
        return self._retrieve(curve_dict, aoa, "AoA")


    def list_tr_values(self) -> List[float]:
        return [float(k.split("=")[1]) for k in self.data.keys()]

    def list_ar_values(self, tr: float) -> List[int]:
        tr_dict = self._retrieve(self.data, f"TR={tr}", "TR")

        return [int(k.split("=")[1]) for k in tr_dict.keys()]

    def list_curves(self, tr: float, ar: int) -> List[str]:
        tr_dict = self._retrieve(self.data, f"TR={tr}", "TR")
        ar_dict = self._retrieve(tr_dict, f"AR={ar}", "AR")

        return list(ar_dict.keys())

    def list_aoa_values(self, tr: float, ar: int, curve: str = "ccl*c/cref") -> List[float]:
        tr_dict = self._retrieve(self.data, f"TR={tr}", "TR")
        ar_dict = self._retrieve(tr_dict, f"AR={ar}", "AR")
        curve_dict = self._retrieve(ar_dict, curve, "Curve")

        return sorted([float(k) for k in curve_dict.keys()])

    def get_curves_for_tr_aoa(
        self, tr: float, aoa: float, curve: str = "ccl*c/cref"
    ) -> Dict[int, Any]:
        """
        Return a dictionary {AR: curve_id} for given TR, AoA, and curve type.
        Example:
            {6: 'spl_xxx', 7: 'spl_yyy', ...}
        """
        tr_dict = self._retrieve(self.data, f"TR={tr}", "TR")
        result = {}
        for ar_key, ar_dict in tr_dict.items():
            ar_val = int(ar_key.split("=")[1])
            curve_dict = ar_dict.get(curve, {})
            if aoa in curve_dict:
                result[ar_val] = curve_dict[aoa]
        
        return result


class AeroLoadCollection:
    def __init__(self, yaml_file: str):
        with open(yaml_file, "r") as f:
            yaml_data = yaml.safe_load(f)

        self.experiments_list = yaml_data.get("experiments", [])

        # Build a dictionary keyed by metadata tuple for fast lookup
        self._lookup: Dict[Tuple, Experiment] = {}
        for exp in self.experiments_list:
            key = (
                exp.get("ROOT_AIRFOIL"),
                exp.get("TIP_AIRFOIL"),
                exp.get("TIP_TWIST"),
                exp.get("SA_LE"),
                exp.get("SA_RE"),
            )
            self._lookup[key] = Experiment(metadata=exp, data=exp.get("data", {}))

    def get_experiment(
        self,
        root_airfoil: str,
        tip_airfoil: str,
        tip_twist: Optional[float] = 0,
        sa_le: Optional[Any] = None,
        sa_re: Optional[Any] = 0,
    ) -> Experiment:
        """Return the Experiment object for given metadata."""
        key = (root_airfoil, tip_airfoil, tip_twist, sa_le, sa_re)
        if key not in self._lookup:
            raise KeyError(f"Experiment with metadata {key} not found.")
        return self._lookup[key]


# ------------------------------
# Example usage
# ------------------------------
if __name__ == "__main__":
    collection = AeroLoadCollection("airfoil_experiments.yaml")

    # Retrieve a specific experiment
    experiment = collection.get_experiment(
        root_airfoil="NACA4412",
        tip_airfoil="NACA4412",
        tip_twist=0,
        sa_le=None,
        sa_re=0,
    )

    # Retrieve curve value
    curve = experiment.get_data(tr=0.5, ar=6, aoa=1.0)
    print("Curve ID:", curve)
