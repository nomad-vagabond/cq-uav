import os, json
from .airfoil import Airfoil

current_dir = os.path.dirname(os.path.abspath(__file__))

def load_airfoils_collection():
    collecion_path = os.path.join(current_dir, "airfoils_collection.json")
    with open(collecion_path) as ac:
        airfoils_data = json.loads(ac.read())
    
    return airfoils_data
