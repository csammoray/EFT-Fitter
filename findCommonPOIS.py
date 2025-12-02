import sys
import argparse

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('--params', dest='params', nargs='+', required=True, 
                        help="List of params files to find common POIs from")
    return parser.parse_args()

def load_pois_from_module(module_name):
    sys.path.append("./params")
    pois_module = __import__(module_name, globals(), locals(), ["pois"], 0)
    if hasattr(pois_module, 'pois'):
        return pois_module.pois
    return {}

def find_common_pois(param_modules):
    all_pois_sets = []
    
    for module_name in param_modules:
        pois = load_pois_from_module(module_name)
        all_pois_sets.append(set(pois.keys()))
    
    if all_pois_sets:
        common_poi_names = set.intersection(*all_pois_sets)
        return sorted(common_poi_names)
    else:
        return []

def main():
    opt = get_options()
    
    common_pois = find_common_pois(opt.params)
    
    if common_pois:
        print("Common POIs:")
        for poi in common_pois:
            print(poi)
    else:
        print("No common POIs found.")

if __name__ == "__main__":
    main()
