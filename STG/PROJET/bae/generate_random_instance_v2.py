import random
import os
from os import listdir
from os.path import isfile, join

def generate_random_instance_dense(filename, nb_periods, nb_items):    
    with open(filename, 'w') as f:
        f.write(f"{nb_periods}\n")
        f.write(f"{nb_items}\n")
        
        for _ in range(nb_items):
            # On remplace les probabilités binaires par des entiers strictement positifs.
            # Entre 5 et 20 garantit que la demande marginale d_t est toujours >= 5
            # ce qui laisse largement la place pour des Delta_t de 1 ou 2.
            demand_line = [str(random.randint(5, 20)) for _ in range(nb_periods)]
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
nb_periods = 10
nb_items = 10

file_list = [f for f in listdir(outdir) if isfile(join(outdir, f))]

for file in file_list:
    os.remove(f"toy_instances/{file}")

if not os.path.exists(outdir):
    os.makedirs(outdir)

for i in range(1, nb_instances+1):
    file_path = f"{outdir}/toy_instance_{i}.txt"
    generate_random_instance_dense(file_path, nb_periods, nb_items)

print(f"========== Generation complete ==========")