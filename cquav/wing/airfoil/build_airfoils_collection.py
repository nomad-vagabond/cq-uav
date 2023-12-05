import os
import re
import json
import requests

import cquav.splinecloud_scipy as scsp


# API_BASE_URL = "http://127.0.0.1:8000/api/repositories/"
API_BASE_URL = "https://splinecloud.com/api/repositories/"

AIRFOIL_REPOS = [
    "rpo_Qo8SUrqwUoYM", # NACA airfoil database
    "rpo_M0syMhGxDgHA", # Symmetrical airfoils database
    "rpo_a9pWAT9ndbAa", # Airfoil database A-list
    "rpo_4CEeTxhHRPCX", # Airfoil database B-list
    "rpo_Q505mu93qDyv", # Airfoil database C-list
    "rpo_38wOPSs9ChyP", # Airfoil database I-list
    "rpo_Ax35m27Elm2L", # Airfoil database L-list
    "rpo_JB6jMwMSHQN4", # Airfoil database N-list
]


def extract_profile_name(datafile_path):
    pathsplit =  datafile_path.split('/')
    if len(pathsplit) < 2:
        return
    
    return pathsplit[-2]


def try_append_profile(profile_name, airfoils_collection):
    airfoil_data = airfoils_collection.get(profile_name)
    if not airfoil_data:
        airfoils_collection[profile_name] = {
            "name": profile_name,
            "curves": {}
        }


def process_relation(relation_data, airfoils_collection):
    name = relation_data["name"]

    reynolds_pattern = re.compile(r'Re(\d+)') # matches Reynolds number in this pattern: # Cl v Alpha nacam12-il Re200000
    profile_pattern = re.compile(r'Alpha\s(.*?)\sRe') # matches profile name in this pattern: # Cl v Alpha nacam12-il Re200000

    reynolds = reynolds_pattern.findall(name)
    profile_name = profile_pattern.findall(name)
    
    if not (reynolds and profile_name):
        return

    profile_name = extract_profile_name(relation_data['subset']['dataset']['datafile']['path'])
    if not profile_name:
        return

    reynolds = int(reynolds[0])

    if "Cl" in name:
        curve_type = f"CL_{reynolds}"
    elif "Cd" in name:
        curve_type = f"CD_{reynolds}"
    elif "Cm" in name:
        curve_type = f"CM_{reynolds}"
    else:
        return

    try_append_profile(profile_name, airfoils_collection)
    curve_uid = relation_data["curves"][0]["uid"]
    airfoils_collection[profile_name]["curves"][curve_type] = curve_uid


def process_dataset(dataset_data, airfoils_collection):
    datafile_name = dataset_data['datafile']['name']
    name, ext = os.path.splitext(datafile_name)
    if not ext == ".dat":
        return

    profile_name = extract_profile_name(dataset_data['datafile']['path'])
    if not profile_name:
        return

    try_append_profile(profile_name, airfoils_collection)
    subset_uid = dataset_data['subsets'][0]['uid']
    airfoils_collection[profile_name]["profile"] = subset_uid

    pathsplit = dataset_data["datafile"]["path"].split("/")
    if len(pathsplit) > 2:
        airfoils_collection[profile_name]["group"] = pathsplit[0]
    else:
        airfoils_collection[profile_name]["group"] = dataset_data["datafile"]["repo"]["name"]


def retreive_preformance_curves(repo_id, coefficient, airfoils_collection):
    url = f"{API_BASE_URL}{repo_id}/relations/?keyword={coefficient}"

    go_next = True
    while go_next:
        response = requests.get(url)
        data = response.json()
        for relation_data in data["results"]:
            process_relation(relation_data, airfoils_collection)
        
        go_next = data['next'] != 'None'
        if go_next:
            print("Next page found, downloading", data['next'])
            url = data['next']


def retreive_profile_points(repo_id, airfoils_collection):
    url = f"{API_BASE_URL}{repo_id}/datasets/"

    go_next = True
    while go_next:
        response = requests.get(url)
        data = response.json()
        for dataset_data in data["results"]:
            process_dataset(dataset_data, airfoils_collection)
        
        go_next = data['next'] != 'None'
        if go_next:
            print("Next page found, downloading", data['next'])
            url = data['next']


def missing_curve_warning(airfoil_name, coefficient, reynolds):
    print(f"Skipping airfoil {airfoil_name}. Missing {coefficient} curve for Re{reynolds}. \n")


def check_reynolds(airfoil_data, reynolds):
    has_cl = airfoil_data['curves'].get(f'CL_{reynolds}')
    if not has_cl:
        missing_curve_warning(airfoil_data['name'], "CL", reynolds)

    has_cd = airfoil_data['curves'].get(f'CD_{reynolds}')
    if not has_cd:
        missing_curve_warning(airfoil_data['name'], "CD", reynolds)
    
    has_cm = airfoil_data['curves'].get(f'CM_{reynolds}')
    if not has_cm:
        missing_curve_warning(airfoil_data['name'], "CM", reynolds)

    return has_cl and has_cd and has_cm


def process_repo(repo_uid, airfoils_collection):
    retreive_preformance_curves(repo_uid, "Cl", airfoils_collection)
    retreive_preformance_curves(repo_uid, "Cd", airfoils_collection)
    retreive_preformance_curves(repo_uid, "Cm", airfoils_collection)
    retreive_profile_points(repo_uid, airfoils_collection)


if __name__ == "__main__":
    airfoils_collection = {}
    for repo_uid in AIRFOIL_REPOS:
        process_repo(repo_uid, airfoils_collection)

    airfoils_group_collection = {}
    n_airfoils = 0
    for airfoil_name, airfoil_data in airfoils_collection.items():
        profile_subset = airfoil_data.get('profile')
        has_curves = airfoil_data.get('curves')
        if not profile_subset:
            print(f"Skipping airfoil {airfoil_name}. Mssing profile data.")
        if not has_curves:
            print(f"Skipping airfoil {airfoil_name}. Mssing or poor curves data.")

        if not (profile_subset and has_curves):
            continue

        reynolds_check = [
            check_reynolds(airfoil_data, 50000),
            check_reynolds(airfoil_data, 100000),
            check_reynolds(airfoil_data, 200000),
            check_reynolds(airfoil_data, 500000),
            check_reynolds(airfoil_data, 1000000)
        ]

        if not all(reynolds_check):
            continue

        group = airfoil_data["group"]
        if not airfoils_group_collection.get(group):
            airfoils_group_collection[group] = {}

        collection_group = airfoils_group_collection[group]
        collection_group[airfoil_name] = airfoil_data

        n_airfoils += 1

    print(f"Collected data on {n_airfoils} airfoils")

    with open('airfoils_collection.json', 'w') as ac:
        json.dump(airfoils_group_collection, ac, indent=2)
