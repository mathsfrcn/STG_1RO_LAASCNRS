import random
import os
from os import listdir
from os.path import isfile, join
import numpy as np

def generate_random_instance(file_path: str, nb_periods: int, nb_items: int, margin_settings: int, demand_prob: float, use_periodicity: bool, start_in_period: bool, timespan_period: int, read_instance_rd_lb: float, read_instance_rd_ub: float):    
    with open(file_path, 'w') as f:
        # Parameters
        f.write(f"{nb_periods}\n")
        f.write(f"{nb_items}\n")

        dt_matrix = np.zeros((nb_items, nb_periods), dtype=int)
        
        if use_periodicity:
            if start_in_period:                             # If you want to start during a period of activity
                demand_periodicity: float = demand_prob
            else:
                demand_periodicity: float = 1-demand_prob

            for i in range(nb_items):
                for p in range(nb_periods):
                    if p%timespan_period == 0 and p != 0:   # We flip at the beginning of each period
                        demand_periodicity: float = 1-demand_periodicity
                    
                    base_demand = 1 if random.random() <= demand_periodicity else 0
                    noise: int = random.randint(0, 1)
                    dt_matrix[i, p] = base_demand + noise
                    
                f.write(" ".join(map(str, dt_matrix[i])) + "\n")
        else:
            for i in range(nb_items):
                for p in range(nb_periods):
                    base_demand: Literal[0, 1] = 1 if random.random() <= demand_prob else 0
                    noise: int = random.randint(0, 1)
                    dt_matrix[i, p] = base_demand + noise
                
                f.write(" ".join(map(str, dt_matrix[i])) + "\n")

        # dt, Dt
        dt = dt_matrix.sum(axis=0)
        Dt = np.cumsum(dt)

        # Costs
        cI: int = int(random.uniform(1, 10) + margin_settings)
        cB: int = int(random.uniform(cI/2, (3*cI)/2))
        cP: int = int(random.uniform(cI/2, (3*cI)/2))
        bP: int = int(random.uniform(cI/2, (3*cI)/2))
        f.write(f"{cI} {cB} {cP} {bP}\n")

        X = np.zeros(nb_periods, dtype=int)
        for t in range(nb_periods):
            if Dt[t] == 0:
                X[t] = 0
            else:
                A = int(read_instance_rd_lb * Dt[t])
                B = read_instance_rd_ub * Dt[t]
                X[t] = random.randint(int(B), int(B) + A)   # Generate int between B and B + A
        
        f.write(" ".join(map(str, X)) + "\n")

def main(outdir: str, nb_instances: int, nb_periods: int, nb_items: int, margin_settings: int, demand_prob: float, use_periodicity: bool, start_in_period: bool, timespan_period: int, read_instance_rd_lb: float, read_instance_rd_ub: float):
    if (not outdir
            or (type(use_periodicity) is not bool)
            or (type(start_in_period)  is not bool)
            or (nb_instances <= 0) 
            or (nb_periods <= 0) 
            or (nb_items <= 0)
            or (margin_settings <= 0)
            or not (0 < read_instance_rd_lb < read_instance_rd_ub < 1)
            or not (0 < demand_prob <= 1)
            or not (0 < timespan_period < nb_periods)):
        print("Error: Invalid configuration")
    else:
        file_list: list[str] = [f for f in listdir(outdir) if isfile(join(outdir, f))]

        for file in file_list:
            os.remove(f"toy_instances/{file}")

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        for i in range(1, nb_instances+1):
            file_path: str = f"{outdir}/toy_instance_{i}.txt"
            generate_random_instance(file_path, nb_periods, nb_items, margin_settings, demand_prob, use_periodicity, start_in_period, timespan_period, read_instance_rd_lb, read_instance_rd_ub)

        print("Succes: Generation complete")

##############################
# Generation
##############################

outdir = "toy_instances"
use_periodicity: bool      = True
start_in_period: bool      = True
timespan_period: int       = 13     # ]0, nb_periods[
nb_instances: int          = 50
nb_periods: int            = 52
nb_items: int              = 20
margin_settings: int       = 10    # Costs will be in [x, x+10[
demand_prob: float         = 0.7
read_instance_rd_lb: float = 0.4   # The production plan will be between lb% and ub% of the cumulative demand
read_instance_rd_ub: float = 0.8

main(outdir, nb_instances, nb_periods, nb_items, margin_settings, demand_prob, use_periodicity, start_in_period, timespan_period, read_instance_rd_lb, read_instance_rd_ub)