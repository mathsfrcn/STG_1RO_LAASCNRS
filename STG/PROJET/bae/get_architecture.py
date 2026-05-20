import os

base_path = os.path.dirname(__file__)

def print_directory_structure(base_path, indent=0):
    for entry in os.listdir(base_path):
        path = os.path.join(base_path, entry)
        print("  " * indent + entry)
        if os.path.isdir(path):
            print_directory_structure(path, indent + 1)

print_directory_structure(base_path)