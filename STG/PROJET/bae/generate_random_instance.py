import random
import os
from os import listdir
from os.path import isfile, join

def generate_random_instance(filename, nb_periods, nb_items, demand_prob=0.3):    
    with open(filename, 'w') as f:

        f.write(f"{nb_periods}\n")
        f.write(f"{nb_items}\n")
        
        for _ in range(nb_items):
            demand_line = ["1" if random.random() <= demand_prob else "0" for _ in range(nb_periods)]
            f.write(" ".join(demand_line) + "\n")
            
        # Stock cost
        f.write("10\n")
        
        # Setup costs
        for i in range(nb_items):
            setup_line = []
            for j in range(nb_items):
                if i == j:
                    setup_line.append("0")
                else:
                    setup_line.append(str(random.randint(100, 200)))
            f.write(" ".join(setup_line) + "\n")
            
        # Optimal bound
        f.write(f"{random.randint(1000, 5000)}\n")

# ==========================================
# GENERATION
# ==========================================

outdir = "toy_instances"
nb_instances = 10
nb_periods = 20
nb_items = 10
demand_prob = 0.3

file_list = [f for f in listdir(outdir) if isfile(join(outdir, f))]

for file in file_list:
    os.remove(f"toy_instances/{file}")

if not os.path.exists(outdir):
    os.makedirs(outdir)

for i in range(1, nb_instances+1):
    file_path = f"{outdir}/toy_instance_{i}.txt"
    generate_random_instance(file_path, nb_periods, nb_items, demand_prob)

print(f"========== Generation complete ==========")