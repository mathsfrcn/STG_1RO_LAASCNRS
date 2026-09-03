#include <cstddef>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <string>
#include <bits/stdc++.h> 
#include <ctime>
#include <algorithm>
#include <stdlib.h>
#include <cstdlib>
#include <iomanip>
#include <dirent.h>
#include <cmath>
#include <sys/types.h>
#include <thread>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <queue>

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

//=========================================== structures

struct Instance{
	int T, cI, cB, bP;
	float Gamma;
	vector<float> deltat;
	vector<float> Dt;
	vector<float> dt;
	vector<float> X;
};

struct Solution{
	Instance inst;
	vector<float> Xt;
	vector<float> xt;
	float obj_val;
};

struct Instance_ADV{
	int T, cI, cB, bP;
	float Gamma;
	vector<float> deltat;
	vector<float> Xt;
	vector<float> xt;
};

struct Solution_ADV{
	Instance_ADV inst;
	vector<float> Dt;	// Cumulative demand
	vector<float> dt;	// Marginal demand

};

struct Benders_Result{
    int iter;
    float time;
	float time_master;
	float time_subproblem;
    float obj_value;
	Solution final_solution;
};

struct Arc_Decision{
    int t;  	// Time
    int i;  	// Budget at the start
    int j;  	// Budget at the end
    int type;   // 0 or 1 (overstock / stockout)

    // Definition of two identical arcs
    bool operator == (const Arc_Decision& other) const{
        return (t == other.t && i == other.i && j == other.j && type == other.type);
    }
};

typedef vector<Arc_Decision> Path;

enum class KC_Method{
	KC,
	RDK,
	RDKL,
	Unique,
	UniqueDual,
	HOG,
	HOGL,
}; 

struct MonteCarlo_Result{
	float mean_cost;
	float std_dev;
	float ci_lower;	// Lower bound for the CI
	float ci_upper;
	float worst_case_simulated;
};


// Structure for A* exploration in the budget graph
struct BFSNode {
    float potential_cost;
    float actual_cost;    // Cumulative cost
    int t;
    int j;
    vector<float> scenario;

    bool operator<(const BFSNode& other) const {
        return potential_cost < other.potential_cost; 
    }
};

// =========================================== Calculate the Manhattan distance for the orthogonality heuristic

int calculate_L1_distance(const Path& pathA, const Path& pathB){
    if(pathA.size() != pathB.size()){
        return 0.0f;
    }

    int sum_diff = 0;

    for(size_t k = 0; k < pathA.size(); k++){
        int delta_A = pathA[k].j - pathA[k].i;
        int delta_B = pathB[k].j - pathB[k].i;

        sum_diff += abs(delta_A - delta_B); 	// Manhattan local distance
    }

    return sum_diff;
}

// =========================================== Recursive extraction of worst-case scenarios following a Depth-First Search

void extract_paths_dfs(
            int t,                                                  
            int j,                                                  // Current node in the backtrack
            const vector<vector<float> >& pi_value,                 // Dynamic programming matrix
            const vector<vector<vector<vector<float> > > >& costs,  // Original costs
            const Solution& sol,
            Path& current_path,
            vector<Path>& all_paths,
            int limit_number_paths,
            float eps){
    
    if(all_paths.size() >= limit_number_paths){
        return;
    }

    // If we went back to the beginning we stop
    if(t == 0){
        Path reversed_path = current_path;  	// Because we start at the end, we have to reverse the path of this branch
        reverse(reversed_path.begin(), reversed_path.end());
        all_paths.push_back(reversed_path);
        return;
    }

    // We are trying to figure out where we came from, its like: what was the budget i? to arrive at j at step t
    for(int i = 0; i <= j; i++){
        if(j <= i + sol.inst.deltat[t-1] && (t != 1 || i == 0)){
            if((j-i) <= sol.inst.dt[t-1]){
                // Check for arc of type 0
                if(abs(pi_value[t][j] - (pi_value[t-1][i] + costs[t][i][j][0])) < eps){
                    Arc_Decision arc = {t, i, j, 0};
                    current_path.push_back(arc);
                    extract_paths_dfs(t-1, i, pi_value, costs, sol, current_path, all_paths, limit_number_paths, eps);
                    current_path.pop_back();
                }
            }
            
            // Check for arc of type 1
            if(abs(pi_value[t][j] - (pi_value[t-1][i] + costs[t][i][j][1])) < eps){
                Arc_Decision arc = {t, i, j, 1};
                current_path.push_back(arc);
                extract_paths_dfs(t-1, i, pi_value, costs, sol, current_path, all_paths, limit_number_paths, eps);
                current_path.pop_back();
            }
        }
    }
}

// =========================================== Function to export the complete graph for visualization (NEEDS CORRECTED)

void export_budget_graph_json(
        const Solution& sol, 
        const vector<vector<float>>& pi_value, 
        const vector<vector<vector<vector<float>>>>& costs, 
        const vector<vector<vector<vector<int>>>>& arcbool){

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    ostringstream oss;
    oss << "graph_visualization/benders_graph_output_" << put_time(&tm, "%Y-%m-%d_%H%M") << ".json";

    ofstream output(oss.str());

    if(!output.is_open()){
        cerr << "Error: The file 'graph_visualization' doesn't exist" << endl;
        return;
    }

    output << "{\n";
    output << "  \"T\": " << sol.inst.T << ",\n";
    output << "  \"Gamma\": " << sol.inst.Gamma << ",\n";
    
    // Nodes export
    output << "  \"nodes\": [\n";
    bool first_node = true;
    for(int t = 0; t < sol.inst.T+2; t++) {
        for(int b = 0; b < sol.inst.Gamma+1; b++) {
            if(!first_node) output << ",\n";
            output << "    {\"id\": \"" << t << "_" << b << "\", \"t\": " << t << ", \"b\": " << b << ", \"pi\": " << pi_value[t][b] << "}";
            first_node = false;
        }
    }

    output << "\n  ],\n";
    
    // Arcs export
    output << "  \"links\": [\n";
    bool first_arc = true;
    for(int t = 1; t < sol.inst.T+2; t++) {
        for(int i = 0; i < sol.inst.Gamma+1; i++) {
            for(int j = i; j < sol.inst.Gamma+1; j++){
                if(t <= sol.inst.T && j > i+sol.inst.deltat[t-1]) continue;
                
                for(int type = 0; type < 2; type++){
                    if(t == sol.inst.T + 1 && j != 0) continue;             // Do not export the arcs of the last layer if it is not the well
                    
                    if(!first_arc) output << ",\n";

                    float cost = (t == sol.inst.T+1) ? 0: costs[t][i][j][type];
                    int is_worst = arcbool[t][i][j][type];
                    
                    output << "    {\"source\": \"" << t-1 << "_" << i << "\", \"target\": \"" << t << "_" << j 
                        << "\", \"type\": " << type << ", \"cost\": " << cost << ", \"is_worst\": " << is_worst << "}";
                    first_arc = false;
                }
            }
        }
    }

    output << "\n  ]\n}\n";
    output.close();
}

//=========================================== Misc.

// Display an float vector, for debugging
void display_vector_float(vector<float> v){
	for (int i = 0; i < v.size(); ++i){
		cout<<v[i]<<" ";
	}

	cout<<endl;
}

// Display an int vector, for debugging
void display_vector_int(vector<int> v){
	for (int i = 0; i < v.size(); ++i){
		cout<<v[i]<<" ";
	}

	cout<<endl;
}

string BoolToString(bool b){
  return b ? "1": "0";
}

// =========================================== Instance related code

vector<float> standardToCumul(vector<float> data){
	vector<float> cumul;
	float tmp = 0;
	for(int i = 0; i < data.size(); i++){
		tmp += data[i];
		cumul.push_back(tmp);
	}

	return cumul;
}

vector<float> cumulToStandard(vector<float> data){
	vector<float> stand;
	stand.push_back(data[0]);
	for(int i = 1; i < data.size(); i++){
		stand.push_back(data[i] - data[i-1]);
	}

	return stand;
}

// Read an instance from psplib problem 58 pspInstance
Instance read_instance_psplib(string filename, int budget){
	Instance inst;
	int nbProd;
	int tmp;

	ifstream file(filename.c_str());
		if(!file){
			cout << "problem with file" << endl;
			exit(-1);
		} 

	file >> inst.T;
	file >> nbProd;
	inst.dt.resize(inst.T);
	for(int i = 0; i < nbProd; i++){
		for(int j = 0; j<inst.T; j++){
			file >> tmp;
			inst.dt[j] += tmp;
		}
	}

	inst.Dt = standardToCumul(inst.dt);

	inst.cI = 3; 	//stock cost
	inst.cB = 6;	//backorder cost
	inst.bP = 10; 	//selling price
	//inst.Gamma = int(inst.Dt[inst.Dt.size()-1]);
	inst.Gamma = budget;
	inst.deltat.resize(inst.T);
	for(int t = 0; t < inst.T;t++){
		inst.deltat[t] = int(inst.Dt[t] / float(2));
		// cout<<inst.deltat[t]<<endl;
	}

	file.close();

	return inst;
}

// Use this if you have to fix parameters
Instance read_hand_instance(string filename, int budget){
	Instance inst;
	int nbProd;
	int tmp;

	ifstream file(filename.c_str());

	if(!file){
		cout << "Error: Problem with file" << endl;
		exit(-1);
	} 

	file >> inst.T;
	file >> nbProd;
	inst.dt.resize(inst.T);

	for(int i = 0; i < nbProd; i++){
		for(int j = 0; j < inst.T; j++){
			file >> tmp;
			inst.dt[j] += tmp;
		}
	}
	
	inst.Dt = standardToCumul(inst.dt);
	inst.cI = 1; 	// Stock cost
	inst.cB = 2; 	// Backorder cost
	inst.bP = 10; 	// Selling price
	inst.Gamma = budget;
	inst.deltat.resize(inst.T);
	inst.X.resize(inst.T);

	for(int t = 0; t < inst.T; t++){
		if(t == 0){
			inst.deltat[t] = 1;
			inst.X[t] = inst.Dt[t];
		} else{
			inst.deltat[t] = 2;
			inst.X[t] = inst.Dt[t];
		}
	}

	file.close();

	return inst;
}

// Use this to read an instance generated by Python.
Instance read_instance_py(string filename, int budget, float adv_margin){
	Instance inst;
    int nbProd;
    int tmp;

    ifstream file(filename.c_str());
    if(!file){
        cout << "Error: Problem with file" << endl;
        exit(-1);
    } 

    file >> inst.T;
    file >> nbProd;
    
    inst.dt.resize(inst.T, 0);
    for(int i = 0; i < nbProd; i++){
        for(int j = 0; j < inst.T; j++){
            file >> tmp;
            inst.dt[j] += tmp; 
        }
    }
    
    // Deterministic calculation of Dt
    inst.Dt = standardToCumul(inst.dt);
    
    // Read parameters
    file >> inst.cB >> inst.cI >> inst.bP;
    inst.Gamma = budget;
    inst.deltat.resize(inst.T);

    // Deterministic calculation of deltat
    for(int t = 0; t<inst.T;t++){
        if(inst.Dt[t] == 0){
            inst.deltat[t] = 0;
        } else {
            float proportion = inst.dt[t] / inst.Dt[inst.T-1];
            inst.deltat[t] = ceil((inst.Gamma * adv_margin) * proportion);

            if(inst.deltat[t] > inst.Dt[t]){
                inst.deltat[t] = inst.Dt[t];
            }
        }
    }

    // Read X
    inst.X.resize(inst.T);
    for(int t = 0; t < inst.T; t++){
        file >> inst.X[t];
    }

    file.close();
    
    return inst;
}

Instance read_instance_randomized(string filename, int budget, float read_instance_rd_lb, float read_instance_rd_ub, float adv_margin){
	Instance inst;
	int nbProd;
	int tmp;

	ifstream file(filename.c_str());
	if(!file){
		cout << "Error: Problem with file" << endl;
		exit(-1);
	} 

	file >> inst.T;
	file >> nbProd;
	inst.dt.resize(inst.T);
	for(int i = 0; i < nbProd; i++){				// For each references
		for(int j = 0; j < inst.T; j++){			// During the period
			file >> tmp;
			inst.dt[j] += tmp + int(rand() % 2);	// According to the production plan, demand is added up to time t
		}
	}
	
	inst.Dt = standardToCumul(inst.dt);
	inst.cB = rand() % 10 + 10; // Backorder cost
	inst.cI = rand() % 10 + 10; // Stock cost
	inst.bP = rand() % 10 + 10; // Selling price
	//inst.Gamma = int(inst.Dt[inst.Dt.size()-1]);
	inst.Gamma = budget;
	inst.Dt = standardToCumul(inst.dt);
	inst.deltat.resize(inst.T);

	//cout << "Display Dt:\n";
	//display_vector_float(inst.Dt);
	for(int t = 0; t<inst.T;t++){
		if(inst.Dt[t] == 0){												// Ofc, if there is no demand we can't set up uncertainty
			inst.deltat[t] = 0;
		} else{
			float proportion = inst.dt[t] / inst.Dt[inst.T-1];				// We calculate the share of demande in périod t relative to total demand
			inst.deltat[t] = ceil((inst.Gamma * adv_margin) * proportion);	// This portion of the total budget is associated with period t

			if(inst.deltat[t] > inst.Dt[t]){								// Security: we cannot cancel more requests than there are
				inst.deltat[t] = inst.Dt[t];
			}
		}
	}

	inst.X.resize(inst.T);

	// Needs corrected
	for(int t = 0; t < inst.T; t++){
		if(inst.Dt[t] == 0){
			inst.X[t] = 0;
		} else{		// Cumulative production may vary between 80% and 120% of cumulative demand.
			inst.X[t] = int(rand() % (int(read_instance_rd_lb*inst.Dt[t])+1) + read_instance_rd_ub*inst.Dt[t]);
		}
	}

	file.close();
	
	return inst;
}

Instance duplicate_instance(Instance inst){
	return inst;
}

// Returns the graph with all the costs for KC_subproblem
vector<vector<vector<vector<float> > > > budget_graph_cost(Solution sol){
	vector<vector<vector<vector<float> > > > costs;

	costs.resize(sol.inst.T+2);

	for(int t = 1; t < sol.inst.T+2; t++){
		costs[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			costs[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				costs[t][i][j].resize(2);
				if(j <= i+sol.inst.deltat[t-1] and j >= i){
					if(t < sol.inst.T){
						costs[t][i][j][0] = sol.inst.cI*(sol.Xt[t-1]- (sol.inst.Dt[t-1] - (j-i)));	// Cost of Inventory
						
						costs[t][i][j][1] = sol.inst.cB*(sol.inst.Dt[t-1] + (j-i) - sol.Xt[t-1]);	// Cost of Backorders
					} else if(t == sol.inst.T){						
						costs[t][i][j][0] = sol.inst.cI*(sol.Xt[t-1]- (sol.inst.Dt[t-1] - (j-i))) - sol.inst.bP*(sol.inst.Dt[t-1] - (j-i));
						
						costs[t][i][j][1] = sol.inst.cB*(sol.inst.Dt[t-1] + (j-i) - sol.Xt[t-1]) - sol.inst.bP*sol.Xt[t-1];
						}
				}
			}

			if(t == sol.inst.T+1){
				costs[t][i][0][0] = 0;
				costs[t][i][0][1] = 0;
			}
		}
	}

	// cout<<"nominal scenario cost:"<<endl;
	// for(int t = 1; t<sol.inst.T+2;t++){
	// 	cout<<costs[t][0][0][0]<<" ";
	// }
	// cout<<endl;
	// for(int t = 1; t<sol.inst.T+2;t++){
	// 	cout<<costs[t][0][0][1]<<" ";
	// }
	// cout<<endl;

	return costs;
}

// Initialize budget graph  with nominal scenario
vector<vector<vector<vector<int> > > > init_graph(Instance inst){
	vector<vector<vector<vector<int> > > > arcbool;
	arcbool.resize(inst.T+2);

	// for t = 0 and t = T+1 only one variable in the array
	for(int t = 0; t < inst.T+2; t++){
		arcbool[t].resize(inst.Gamma+1);
		for(int i = 0; i < inst.Gamma+1; i++){
			arcbool[t][i].resize(inst.Gamma+1);
			for(int j = 0; j < inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}

			if(t >= 1) arcbool[t][0][0][0] = 1;
			if(t >= 1) arcbool[t][0][0][1] = 1;
		}
	}

	return arcbool;
}

// Initialize budget graph with ALL scenarios
vector<vector<vector<vector<int> > > > init_graph_full(Instance inst){

	vector<vector<vector<vector<int> > > > arcbool;
	arcbool.resize(inst.T+2);

	// for t = 0 and t = T+1 only one variable in the array
	for(int t = 0; t<inst.T+2;t++){
		arcbool[t].resize(inst.Gamma+1);
		for(int i = 0; i<inst.Gamma+1; i++){
			arcbool[t][i].resize(inst.Gamma+1);
			for(int j = 0; j<inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
				arcbool[t][i][j][0] = 1;
				arcbool[t][i][j][1] = 1;
			}
			
		}
	}

	return arcbool;
}

// Return the union of budget graphs
vector<vector<vector<vector<int> > > > merge_budget_graph(vector<vector<vector<vector<int> > > > arcbool1, vector<vector<vector<vector<int> > > > arcbool2){
	vector<vector<vector<vector<int> > > > merge_arcbool;

	merge_arcbool.resize(arcbool1.size());

	for(int t = 0; t < arcbool1.size(); t++){
		merge_arcbool[t].resize(arcbool1[t].size());
		for(int i = 0; i < arcbool1[t].size(); i++){
			merge_arcbool[t][i].resize(arcbool1[t][i].size());
			for(int j = 0; j < arcbool1[t][i].size(); j++){
				merge_arcbool[t][i][j].resize(2);
				
				if(arcbool1[t][i][j][0] == 1 or arcbool2[t][i][j][0] == 1){
					merge_arcbool[t][i][j][0] = 1;
				}

				if(arcbool1[t][i][j][1] == 1 or arcbool2[t][i][j][1] == 1){
					merge_arcbool[t][i][j][1] = 1;
				}
			}
		}
	}

	return merge_arcbool;
}

//=========================================== Solution related code

float objective_value(Solution sol, vector<float> Dt){
	//cout<<"============= obj_value: "<<endl;
	float obj = 0;
	for(int i = 0; i < sol.inst.T; i++){
		obj += max(sol.inst.cI*(sol.Xt[i] - Dt[i]), sol.inst.cB*(Dt[i] - sol.Xt[i]));
	}

	obj -=  sol.inst.bP*min(Dt[sol.inst.T-1], sol.Xt[sol.inst.T-1]); 

	return obj;
}

MonteCarlo_Result run_monte_carlo(const Solution& sol, int num_scenarios) {
    vector<float> simulated_costs(num_scenarios);
    float sum_costs = 0.0;
    float worst_cost = -1.0;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);

    for(int k = 0; k < num_scenarios; k++) {
        vector<float> Dt_random(sol.inst.T, 0.0);
        float budget_consomme = 0.0;
        
        // The deviation is generated directly based on the cumulative demand.
        for(int t = 0; t < sol.inst.T; t++) {
            float variation = 0.0;
            
            if(sol.inst.deltat[t] > 0) {
                float max_var = sol.inst.deltat[t];
                
                // Gamma safety
                if(budget_consomme + max_var > sol.inst.Gamma) {
                    max_var = sol.inst.Gamma - budget_consomme;
                    if(max_var < 0) max_var = 0;
                }
                
                std::uniform_real_distribution<float> dist(-max_var, max_var);
                variation = dist(rng);
                
                budget_consomme += abs(variation);
            }
            
            // Application of the variation based on the CUMULATIVE nominal base
            Dt_random[t] = sol.inst.Dt[t] + variation;
            
            // Physique safety
            if(t > 0 && Dt_random[t] < Dt_random[t-1]) {
                Dt_random[t] = Dt_random[t-1];
            }
        }
        
        // Evaluation of the production plan
        float cost = objective_value(sol, Dt_random);
        simulated_costs[k] = cost;
        sum_costs += cost;
        if(cost > worst_cost) worst_cost = cost;
    }
    
    // Statistics (LGN/TCL)
    float mean = sum_costs / num_scenarios;
    float variance = 0.0;
    for(int k = 0; k < num_scenarios; k++) {
        variance += pow(simulated_costs[k] - mean, 2);
    }
    variance /= num_scenarios;
    float std_dev = sqrt(variance);
    
    float margin_of_error = 1.96f * (std_dev / sqrt(num_scenarios));
    
    return {mean, std_dev, mean - margin_of_error, mean + margin_of_error, worst_cost};
}

// ========================================================================================================================================================================================================
// ========================================================================================================================================================================================================
// ================================================================================= Knowledge Compilation Benders Decomposition code =====================================================================
// ========================================================================================================================================================================================================
// ========================================================================================================================================================================================================

Solution KC_benders_Master(Instance inst, vector<vector<vector<vector<int> > > > arcsol){
	Solution sol;
	IloEnv env;
	IloModel model(env);

	// Used to retrieve relevant pi after the solve
	vector<vector<int> > pibool;
	pibool.resize(inst.T+2);

	// Vars
	IloArray<IloNumVarArray> pi(env,inst.T+2);
	// for t = 0 and t = T+1 only one variable in the array
	for(int t = 0; t < inst.T+2; t++){
		pi[t] = IloNumVarArray(env, inst.Gamma+1);
		pibool[t].resize(inst.Gamma+1);
		for(int i = 0; i < inst.Gamma+1; i++){
			char name[80];
			pi[t][i] = IloNumVar(env, -IloInfinity, IloInfinity);
			// cout<<pi[t][i].getLB()<<endl;
			sprintf(name,"pi_%d_%d",t,i);
			pi[t][i].setName(name);
		}
	}

	IloNumVarArray X(env, inst.T);
	for(int t = 0; t < inst.T; t++){
		X[t] = IloNumVar(env);
		char name[80];
		sprintf(name,"X_%d",t);
		X[t].setName(name);
	}

	// Upper bound of cumulative production
	for(int t = 1; t < inst.T+1; t++){
		model.add(X[t-1] <= inst.X[t-1]);
	}

	// Lower bound of cumulative production
	for(int t = 1; t < inst.T; t++){
    	model.add(X[t] >= X[t-1]);
	}

	// for t in 1...T-1 
	for(int t = 1; t < inst.T; t++){
		for(int i = 0; i < inst.Gamma+1; i++){
			for(int j = i; j < inst.Gamma+1; j++){
				if(j <= i+inst.deltat[t-1] and (arcsol[t][i][j][0] or arcsol[t][i][j][1])){
					IloExpr expr(env);
					if(t == 1) {
						expr = pi[0][0];
						pibool[0][0] = 1;
					} else{ expr = pi[t-1][i];}
					model.add(pi[t][j] - expr >= inst.cI*(X[t-1] - (inst.Dt[t-1] - (j-i))));
					model.add(pi[t][j] - expr >= inst.cB*(inst.Dt[t-1] + (j-i) - X[t-1]));
					pibool[t][j] = 1;
				}
			}
		}
	}

	// for t = T (T-1 -> T)
	int t = inst.T;
	for(int i = 0; i < inst.Gamma+1; i++){
		for(int j = i; j < inst.Gamma+1; j++){
			if(j <= i+inst.deltat[t-1] and (arcsol[t][i][j][0] or arcsol[t][i][j][1])){
				model.add(pi[t][j] - pi[t-1][i] >= inst.cI*(X[t-1] - (inst.Dt[t-1] - (j-i))) - inst.bP*(inst.Dt[t-1] - (j-i)));
				model.add(pi[t][j] - pi[t-1][i] >= inst.cB*(inst.Dt[t-1] + (j-i) - X[t-1]) - inst.bP*X[t-1]);
				pibool[t][j] = 1;
			}
		}
	}

	// Modeling the arcs of the last (T -> T+1) layer
	for(int i = 0; i < inst.Gamma+1; i++){
		if(pibool[t][i] == 1){
			model.add(pi[t+1][0] - pi[t][i] >= 0);
		}
	}
	
	model.add(pi[0][0]==0);

	pibool[0][0] = 1;
	pibool[inst.T+1][0] = 1;

	// Objective Value

	model.add(IloMinimize(env, pi[inst.T+1][0]));

	IloCplex cplex(model);
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
	cplex.setOut(env.getNullStream());
    if(!cplex.solve()){
    	env.error() << "Failed to optimize LP." << endl;
    	throw(-1);
	}

	// cout<<"budget graph:"<<endl;
	// for(int i = inst.Gamma; i>=0; i--){
	// 	stringstream buff;
	// 	for(int t = 0; t<inst.T+2; t++){
	// 		buff<< " ";
	// 		buff<<pibool[t][i];
	// 	}
	// 	cout<<buff.str()<<endl;
	// }

	sol.obj_val = cplex.getValue(pi[inst.T+1][0]);
	sol.Xt.resize(inst.T);
	for(int t = 0; t < inst.T; t++){
		sol.Xt[t] = cplex.getValue(X[t]);
	}

	sol.inst = inst;

	env.end();

	return sol;
}

vector<vector<vector<vector<int> > > > KC_benders_Subproblem(Solution sol, float approx_coeff, bool use_graph_export, float& ub_cost, float eps, float p_few){
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t < sol.inst.T+2; t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;	// Start at period 0 cost 0
	for(int t = 1; t < sol.inst.T+1; t++){
		for(int j = 0; j < sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j] + costs[t][j][j][0];
			for(int i = 0; i <= j; i++){
				if(j <= i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i] + costs[t][i][j][0] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][0];
					}

					if(pi_value[t-1][i] + costs[t][i][j][1] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					}
				}
			}

			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] > tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
	}

	pi_value[sol.inst.T+1][0] = tmp;
	
	// Its the value of the longest path (i.e. the obj_value for the Adv)
	ub_cost = pi_value[sol.inst.T+1][0];

	// ========================== Now the backtrack
	
	float sub_OPT;

	if(ub_cost>=0){
		sub_OPT = approx_coeff*ub_cost;
	} else{
		sub_OPT = (1-approx_coeff)*ub_cost+ub_cost;
	}
	
	// t = T+1
	pi_subopt_bool[sol.inst.T+1][0] = true;
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] >= sub_OPT){
			arcbool[sol.inst.T+1][i][0][0] = 1;
			arcbool[sol.inst.T+1][i][0][1] = 1;
			pi_subopt_bool[sol.inst.T][i] = true;
		}
	}

	// Matrix for storing worst-cost paths
    vector<vector<bool>> is_elite_node(sol.inst.T + 2, vector<bool>(sol.inst.Gamma + 1, false));
    
    // Initialization of elite nodes at time T
    for(int i = 0; i < sol.inst.Gamma+1; i++){
        if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
            is_elite_node[sol.inst.T][i] = true;
        }
    }

	for(int t = sol.inst.T; t > 0; t--){
        for(int j = 0; j < sol.inst.Gamma+1; j++){
            if(pi_subopt_bool[t][j]){
                bool has_incoming_arc = false;
                int last_valid_i = -1;
                int last_valid_type = -1;

                for(int i = 0; i <= j; i++){
                    if(j <= i+sol.inst.deltat[t-1] and (t != 1 or i == 0)){
                        if(abs(pi_value[t][j] - (pi_value[t-1][i]+costs[t][i][j][0])) < eps){	// Arc 0
                            last_valid_i = i; last_valid_type = 0;
                            
                            if(is_elite_node[t][j]){											// If we are on the path to the worst-case scenario, we approve it automatically
                                arcbool[t][i][j][0] = 1;
                                pi_subopt_bool[t-1][i] = true;
                                is_elite_node[t-1][i] = true;
                                has_incoming_arc = true;
                            }
                        }

                        if(abs(pi_value[t][j] - (pi_value[t-1][i]+costs[t][i][j][1])) < eps){	// Arc 1
                            last_valid_i = i; last_valid_type = 1;
                            
                            if(is_elite_node[t][j]){
                                arcbool[t][i][j][1] = 1;
                                pi_subopt_bool[t-1][i] = true;
                                is_elite_node[t-1][i] = true;
                                has_incoming_arc = true;
                            } else if((float)rand() / RAND_MAX < p_few){
                                arcbool[t][i][j][1] = 1;
                                pi_subopt_bool[t-1][i] = true;
                                has_incoming_arc = true;
                            }
                        }
                    }
                }

				// We take the last valid arc if no arc has been selected for this node, to ensure connectivity in the subgraph
                if(!has_incoming_arc && last_valid_i != -1){
                    arcbool[t][last_valid_i][j][last_valid_type] = 1;
                    pi_subopt_bool[t-1][last_valid_i] = true;
                    if(is_elite_node[t][j]) is_elite_node[t-1][last_valid_i] = true;
                }
            }
        }
    }

	if(use_graph_export){
        export_budget_graph_json(sol, pi_value, costs, arcbool);
    }

	// Display the subgraph
	// cout<<"subgraph:"<<endl;
	// stringstream bufft;
	// for(int t = 0; t<sol.inst.T+2; t++){
	// 	bufft<<t;
	// 	bufft<<" ";
	// }
	// cout<<bufft.str()<<endl;
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	string buff = "";
	// 	for(int t = 0; t<sol.inst.T+2; t++){
	// 		buff += BoolToString(pi_subopt_bool[t][i])+" " ;
	// 	}
	// 	cout<<buff<<endl;
	// }

	// cout<<endl;
	// //display the subgraph
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	stringstream buff;
	// 	for(int t = 0; t<sol.inst.T+2; t++){
	// 		buff<< " ";
	// 		buff<<pi_value[t][i];
	// 	}
	// 	cout<<buff.str()<<endl;
	// }

	return arcbool;
}

// KC with random K
vector<vector<vector<vector<int> > > > KC_benders_Subproblem_RDK(Solution sol, float approx_coeff, int nb_path_to_select, bool use_graph_export, int limit_number_paths, float& ub_cost, float eps){
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t < sol.inst.T+2; t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;									// Start at period 0 cost 0
	for(int t = 1; t < sol.inst.T+1; t++){
		for(int j = 0; j < sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j] + costs[t][j][j][0];	// It's the value of the dual problem that will store the value of the longest path from the start to t, having consumed j units of budget
			for(int i = 0; i <= j; i++){
				if(j <= i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i] + costs[t][i][j][0] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][0];
					}

					if(pi_value[t-1][i] + costs[t][i][j][1] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					} 
				}
			}

			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] > tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
	}

	pi_value[sol.inst.T+1][0] = tmp;
	ub_cost = pi_value[sol.inst.T+1][0];

	// ========================== Now the backtrack

	float sub_OPT;

	if(ub_cost>=0){
		sub_OPT = approx_coeff*ub_cost;
	} else{
		sub_OPT = (1-approx_coeff)*ub_cost+ub_cost;
	}

	vector<Path> candidates_path;
	Path current_path_buffer;
	vector<Path> optimal_paths;
	vector<Path> selected_paths;
	auto rng = std::default_random_engine(std::random_device{}());

	for(int i = 0; i < sol.inst.Gamma+1; i++){	// We extract optimal paths only
		if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, optimal_paths, limit_number_paths, eps);
		}
	}

	if(!optimal_paths.empty()){					// Optimal convergence
		selected_paths.push_back(optimal_paths[0]);
	}

	for(int i = 0; i < sol.inst.Gamma+1; i++){	// We extract sub-optimal paths
		if(pi_value[sol.inst.T][i] >= sub_OPT){
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, candidates_path, limit_number_paths, eps);
		}
	}
	
	if(!candidates_path.empty()){	// Candidate paths are mixed in order to randomly select some from the relaxed set.
		if(!optimal_paths.empty()){	// The selected path is removed to avoid duplicates.
			candidates_path.erase(
				std::remove(candidates_path.begin(), candidates_path.end(), optimal_paths[0]), candidates_path.end()
			);
			
			std::shuffle(std::begin(candidates_path), std::end(candidates_path), rng);

			int paths_needed = nb_path_to_select - selected_paths.size();
			for(int i = 0; i < paths_needed && i < candidates_path.size(); i++){
				selected_paths.push_back(candidates_path[i]);
			}
		}
	}

	// Connexion
	for(int t = 0; t < sol.inst.T+2; t++){
		for(int i = 0 ; i < sol.inst.Gamma+1; i++){
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j][0] = 0;
				arcbool[t][i][j][1] = 0;
			}
		}
	}

	for(size_t p = 0; p < selected_paths.size(); p++){
		for(size_t a = 0; a < selected_paths[p].size(); a++){
			Arc_Decision arc = selected_paths[p][a];
			arcbool[arc.t][arc.i][arc.j][arc.type] = 1;
		}

		int final_budget = selected_paths[p].back().j;
		arcbool[sol.inst.T+1][final_budget][0][0] = 1;
		arcbool[sol.inst.T+1][final_budget][0][1] = 1;
	}

	if(use_graph_export){
        export_budget_graph_json(sol, pi_value, costs, arcbool);
    }

	return arcbool;
}

// KC with random K
vector<vector<vector<vector<int> > > > KC_benders_Subproblem_RDKL(Solution sol, float approx_coeff, int nb_path_to_select, bool use_graph_export, int limit_number_paths, float& ub_cost, float eps){
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t < sol.inst.T+2; t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;	// Start at period 0 cost 0
	for(int t = 1; t < sol.inst.T+1; t++){
		for(int j = 0; j < sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j] + costs[t][j][j][0];
			for(int i = 0; i <= j; i++){
				if(j <= i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i] + costs[t][i][j][0] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][0];
					}

					if(pi_value[t-1][i] + costs[t][i][j][1] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					} 
				}
			}

			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] > tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
	}

	pi_value[sol.inst.T+1][0] = tmp;
	ub_cost = pi_value[sol.inst.T+1][0];
	
	// ========================== Now the backtrack

	float sub_OPT;

	if(ub_cost >= 0){
		sub_OPT = approx_coeff * ub_cost;
	} else{
		sub_OPT = (1-approx_coeff) * ub_cost + ub_cost;
	}

	vector<Path> candidates_path;
	Path current_path_buffer;
	vector<Path> optimal_paths;
	vector<Path> selected_paths;
	auto rng = std::default_random_engine(std::random_device{}());

	for(int i = 0; i < sol.inst.Gamma+1; i++){									// We extract optimal paths
		if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, optimal_paths, limit_number_paths, eps);
		}
	}

	if(!optimal_paths.empty()){
		std::shuffle(std::begin(optimal_paths), std::end(optimal_paths), rng); 	// We choose N random optimal paths

		int paths_to_take = std::min((int)optimal_paths.size(), nb_path_to_select);
		for(int i = 0; i < paths_to_take; i++){
			selected_paths.push_back(optimal_paths[i]);
		}
	}

	if(selected_paths.size() < nb_path_to_select){								// If we dont have enough paths, we complete by random sub-optimal paths
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			if(pi_value[sol.inst.T][i] >= sub_OPT && abs(pi_value[sol.inst.T][i] - ub_cost) >= eps){
				extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, candidates_path, limit_number_paths, eps);
			}
		}

		if(!candidates_path.empty()){
			std::shuffle(std::begin(candidates_path), std::end(candidates_path), rng);

			int paths_needed = nb_path_to_select - selected_paths.size();
			for(int i = 0; i < paths_needed && i < candidates_path.size(); i++){
				selected_paths.push_back(candidates_path[i]);
			}
		}
	}

	// Connexion
	for(int t = 0; t < sol.inst.T+2; t++){
		for(int i = 0 ; i < sol.inst.Gamma+1; i++){
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j][0] = 0;
				arcbool[t][i][j][1] = 0;
			}
		}
	}

	for(size_t p = 0; p < selected_paths.size(); p++){
		for(size_t a = 0; a < selected_paths[p].size(); a++){
			Arc_Decision arc = selected_paths[p][a];
			arcbool[arc.t][arc.i][arc.j][arc.type] = 1;
		}

		int final_budget = selected_paths[p].back().j;
		arcbool[sol.inst.T+1][final_budget][0][0] = 1;
		arcbool[sol.inst.T+1][final_budget][0][1] = 1;
	}

	if(use_graph_export){
        export_budget_graph_json(sol, pi_value, costs, arcbool);
    }

	return arcbool;
}

// KC with only one optimal path pushed back to te master at each iteration
vector<vector<vector<vector<int> > > > KC_benders_Subproblem_Unique(Solution sol, float eps, float& ub_cost){
    vector<vector<vector<vector<int> > > > arcbool; 
    vector<vector<float> > pi_value; 				
    vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); 

    pi_value.resize(sol.inst.T+2);
    arcbool.resize(sol.inst.T+2);
    for(int t = 0; t < sol.inst.T+2; t++){
        pi_value[t].resize(sol.inst.Gamma+1);
        arcbool[t].resize(sol.inst.Gamma+1);
        for(int i = 0; i < sol.inst.Gamma+1; i++){
            arcbool[t][i].resize(sol.inst.Gamma+1);
            for(int j = 0; j < sol.inst.Gamma+1; j++){
                arcbool[t][i][j].resize(2, 0);
            }
        }
    }

    // Dynamic prog. for longest path

    float tmp;
    pi_value[0][0] = 0;
    for(int t = 1; t < sol.inst.T+1; t++){
        for(int j = 0; j < sol.inst.Gamma+1; j++){
            tmp = pi_value[t-1][j] + costs[t][j][j][0];
            for(int i = 0; i <= j; i++){
                if(j <= i+sol.inst.deltat[t-1]){
                    if(pi_value[t-1][i] + costs[t][i][j][0] > tmp){
                        tmp = pi_value[t-1][i]+costs[t][i][j][0];
                    }

                    if(pi_value[t-1][i] + costs[t][i][j][1] > tmp){
                        tmp = pi_value[t-1][i]+costs[t][i][j][1];
                    } 
                }
            }

            pi_value[t][j] = tmp;
        }
    }

    // We calculate the worst cost
    tmp = pi_value[sol.inst.T][0];
    for(int i = 0; i < sol.inst.Gamma+1; i++){
        if(pi_value[sol.inst.T][i] > tmp){
            tmp = pi_value[sol.inst.T][i];
        }
    }

    pi_value[sol.inst.T+1][0] = tmp;
    ub_cost = tmp;

    // ========================== Now the backtrack

    int current_j = -1;
    
    for(int i = 0; i < sol.inst.Gamma+1; i++){
        if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
            current_j = i;
            break; 									// We stop when we have the first scenario
        }
    }

    if(current_j != -1){
        arcbool[sol.inst.T+1][current_j][0][0] = 1;	// Connexion T -> T+1
        arcbool[sol.inst.T+1][current_j][0][1] = 1;

        for(int t = sol.inst.T; t > 0; t--){
            bool found_arc = false;
            for(int i = 0; i <= current_j; i++){
                if(current_j <= i+sol.inst.deltat[t-1] && (t != 1 || i == 0)){
                    // Arc 0
                    if(abs(pi_value[t][current_j] - (pi_value[t-1][i]+costs[t][i][current_j][0])) < eps){
                        arcbool[t][i][current_j][0] = 1;
                        current_j = i;
                        found_arc = true;
                        break; 
                    }

                    // Arc 1
                    if(abs(pi_value[t][current_j] - (pi_value[t-1][i]+costs[t][i][current_j][1])) < eps){
                        arcbool[t][i][current_j][1] = 1;
                        current_j = i;
                        found_arc = true;
                        break;
                    }
                }
            }

            if(!found_arc) break;
        }
    }

    return arcbool;
}

// KCU with two optimal paths pushed back to te master at each iteration
vector<vector<vector<vector<int> > > > KC_benders_Subproblem_Unique_Dual(Solution sol, float eps, float& ub_cost){
    vector<vector<vector<vector<int> > > > arcbool; 
    vector<vector<float> > pi_value; 				
    vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); 

    pi_value.resize(sol.inst.T+2);
    arcbool.resize(sol.inst.T+2);
    for(int t = 0; t < sol.inst.T+2; t++){
        pi_value[t].resize(sol.inst.Gamma+1);
        arcbool[t].resize(sol.inst.Gamma+1);
        for(int i = 0; i < sol.inst.Gamma+1; i++){
            arcbool[t][i].resize(sol.inst.Gamma+1);
            for(int j = 0; j < sol.inst.Gamma+1; j++){
                arcbool[t][i][j].resize(2, 0); // Init to 0
            }
        }
    }

    // Dynamic prog. for longest path

    float tmp;
    pi_value[0][0] = 0;
    for(int t = 1; t < sol.inst.T+1; t++){
        for(int j = 0; j < sol.inst.Gamma+1; j++){
            tmp = pi_value[t-1][j] + costs[t][j][j][0];
            for(int i = 0; i <= j; i++){
                if(j <= i+sol.inst.deltat[t-1]){
                    if(pi_value[t-1][i] + costs[t][i][j][0] > tmp){
                        tmp = pi_value[t-1][i]+costs[t][i][j][0];
                    }

                    if(pi_value[t-1][i] + costs[t][i][j][1] > tmp){
                        tmp = pi_value[t-1][i]+costs[t][i][j][1];
                    } 
                }
            }

            pi_value[t][j] = tmp;
        }
    }

    // Same as before
    tmp = pi_value[sol.inst.T][0];
    for(int i = 0; i < sol.inst.Gamma+1; i++){
        if(pi_value[sol.inst.T][i] > tmp){
            tmp = pi_value[sol.inst.T][i];
        }
    }

    pi_value[sol.inst.T+1][0] = tmp;
    ub_cost = tmp;

    // ========================== Now the backtrack
    
	int start_j_late = -1;
	int start_j_early = -1;

	// We are looking for the optimal nodes
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
			if(start_j_late == -1) start_j_late = i;
			start_j_early = i;
		}
	}
	
	if(start_j_late != -1){
		// Connexion to the well
		arcbool[sol.inst.T+1][start_j_late][0][0] = 1;
		arcbool[sol.inst.T+1][start_j_late][0][1] = 1;
		arcbool[sol.inst.T+1][start_j_early][0][0] = 1;
		arcbool[sol.inst.T+1][start_j_early][0][1] = 1;

		// Late Consumption
		int current_j = start_j_late;
		for(int t = sol.inst.T; t > 0; t--){
			bool found = false;
			for(int i = 0; i <= current_j; i++){
				if(current_j <= i+sol.inst.deltat[t-1] && (t != 1 || i == 0)){
					if(abs(pi_value[t][current_j] - (pi_value[t-1][i]+costs[t][i][current_j][0])) < eps){
						arcbool[t][i][current_j][0] = 1; current_j = i; found = true; break; 
					}
					if(abs(pi_value[t][current_j] - (pi_value[t-1][i]+costs[t][i][current_j][1])) < eps){
						arcbool[t][i][current_j][1] = 1; current_j = i; found = true; break;
					}
				}
			}
			if(!found) break;
		}

		// Early Consumption
		current_j = start_j_early;
		for(int t = sol.inst.T; t > 0; t--){
			bool found = false;
			for(int i = current_j; i >= 0; i--){
				if(current_j <= i+sol.inst.deltat[t-1] && (t != 1 || i == 0)){
					if(abs(pi_value[t][current_j] - (pi_value[t-1][i]+costs[t][i][current_j][0])) < eps){
						arcbool[t][i][current_j][0] = 1; 
						current_j = i; 
						found = true; 
						break; 
					}

					if(abs(pi_value[t][current_j] - (pi_value[t-1][i]+costs[t][i][current_j][1])) < eps){
						arcbool[t][i][current_j][1] = 1; 
						current_j = i; 
						found = true; 
						break;
					}
				}
			}

			if(!found) break;
		}
	}

    return arcbool;
}

// KCOG
vector<vector<vector<vector<int> > > > KC_benders_Subproblem_OG(Solution sol, float approx_coeff, int nb_path_to_select, bool use_graph_export, int limit_number_paths, float& ub_cost, float eps){
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t < sol.inst.T+2; t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;	// Start at period 0 cost 0
	for(int t = 1; t < sol.inst.T+1; t++){
		for(int j = 0; j < sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j] + costs[t][j][j][0];	// It's the value of the dual problem that will store the value of the longest path from the start to t, having consumed j units of budget
			for(int i = 0; i <= j; i++){
				if(j <= i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i] + costs[t][i][j][0] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][0];
					}

					if(pi_value[t-1][i] + costs[t][i][j][1] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					} 
				}
			}

			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] > tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
	}

	pi_value[sol.inst.T+1][0] = tmp;
	ub_cost = pi_value[sol.inst.T+1][0];

	// ========================== Now the backtrack

	float sub_OPT;

	if(ub_cost>=0){
		sub_OPT = approx_coeff*ub_cost;
	} else{
		sub_OPT = (1-approx_coeff)*ub_cost+ub_cost;
	}
	
	vector<Path> candidates_path;				// Store the worst-case scenarios
	Path current_path_buffer; 					// Use for recursion
	vector<Path> optimal_paths;
	vector<Path> selected_paths;				// Scenarios we push up to the master

	for(int i = 0; i < sol.inst.Gamma+1; i++){	// Optimal scenarios
		if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, optimal_paths, limit_number_paths, eps);
		}
	}

	if(!optimal_paths.empty()){					// We take the first optimal path we check
		selected_paths.push_back(optimal_paths[0]);
	}

	for(int i = 0; i < sol.inst.Gamma+1; i++){	// Sub-optimal paths
		if(pi_value[sol.inst.T][i] >= sub_OPT){
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, candidates_path, limit_number_paths, eps);
		}
	}

	if(!optimal_paths.empty()){ 				// We delete the first we took
		candidates_path.erase(
			std::remove(candidates_path.begin(), candidates_path.end(), optimal_paths[0]), candidates_path.end()
		);
	}

	// MaxMin L1 distance
	if(!candidates_path.empty()){
		while(selected_paths.size() < nb_path_to_select && !candidates_path.empty()){
			float best_max_min_distance = -1.0;
			int best_candidate_index = -1;
            
			for(size_t c = 0; c < candidates_path.size(); c++){
				float min_distance_selected = 1e9;

				for(size_t s = 0; s < selected_paths.size(); s++){
					float dist = calculate_L1_distance(candidates_path[c], selected_paths[s]);
					if(dist < min_distance_selected){
						min_distance_selected = dist;
					}
				}

				if(min_distance_selected > best_max_min_distance){
					best_max_min_distance = min_distance_selected;
					best_candidate_index = c;
				}
			}

			selected_paths.push_back(candidates_path[best_candidate_index]);
			candidates_path.erase(candidates_path.begin() + best_candidate_index);
		}
	}

	// Reset arcbool
	for(int t = 0; t < sol.inst.T+2; t++){
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j][0] = 0;
				arcbool[t][i][j][1] = 0;
			}
		}
	}

	// Flip only the arcs that are on the paths
	for(size_t p = 0; p < selected_paths.size(); p++){
		for(size_t a = 0; a < selected_paths[p].size(); a++){
			Arc_Decision arc = selected_paths[p][a];
			arcbool[arc.t][arc.i][arc.j][arc.type] = 1;
		}

		// Reconnected the end of the path to node T+1 like in the previous code
		int final_budget = selected_paths[p].back().j;
		arcbool[sol.inst.T+1][final_budget][0][0] = 1;
		arcbool[sol.inst.T+1][final_budget][0][1] = 1;
	}

	if(use_graph_export){
        export_budget_graph_json(sol, pi_value, costs, arcbool);
    }

	// Display the subgraph
	// cout<<"subgraph:"<<endl;
	// stringstream bufft;
	// for(int t = 0; t<sol.inst.T+2; t++){
	// 	bufft<<t;
	// 	bufft<<" ";
	// }
	// cout<<bufft.str()<<endl;
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	string buff = "";
	// 	for(int t = 0; t<sol.inst.T+2; t++){
	// 		buff += BoolToString(pi_subopt_bool[t][i])+" " ;
	// 	}
	// 	cout<<buff<<endl;
	// }

	// cout<<endl;
	// //display the subgraph
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	stringstream buff;
	// 	for(int t = 0; t<sol.inst.T+2; t++){
	// 		buff<< " ";
	// 		buff<<pi_value[t][i];
	// 	}
	// 	cout<<buff.str()<<endl;
	// }

	return arcbool;
}

// KCOGL (KCOG Lexicographical)
vector<vector<vector<vector<int> > > > KC_benders_Subproblem_OGL(Solution sol, float approx_coeff, int nb_path_to_select, bool use_graph_export, int limit_number_paths, float& ub_cost, float eps){
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t < sol.inst.T+2; t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;	// Start at period 0 cost 0
	for(int t = 1; t < sol.inst.T+1; t++){
		for(int j = 0; j < sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j] + costs[t][j][j][0];
			for(int i = 0; i <= j; i++){
				if(j <= i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i] + costs[t][i][j][0] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][0];
					}

					if(pi_value[t-1][i] + costs[t][i][j][1] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					} 
				}
			}

			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] > tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
	}

	pi_value[sol.inst.T+1][0] = tmp;
	ub_cost = pi_value[sol.inst.T+1][0];

	// ========================== Now the backtrack

	float sub_OPT;

	if(ub_cost>=0){
		sub_OPT = approx_coeff*ub_cost;
	} else{
		sub_OPT = (1-approx_coeff)*ub_cost+ub_cost;
	}
	
	vector<Path> candidates_path;				// Store the worst-case scenarios
	Path current_path_buffer; 					// Use for recursion
	vector<Path> optimal_paths;
	vector<Path> selected_paths;				// Scenarios we push up to the master

	for(int i = 0; i < sol.inst.Gamma+1; i++){	// Optimal scenarios
		if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, optimal_paths, limit_number_paths, eps);
		}
	}

	// MaMin on the optimal set
	// maybe we can optimize this section by "memoïsation" 
	if(!optimal_paths.empty()){
		selected_paths.push_back(optimal_paths[0]);
        optimal_paths.erase(optimal_paths.begin());
		
		while(selected_paths.size() < nb_path_to_select && !optimal_paths.empty()){
			float best_max_min_distance = -1.0;
			int best_candidate_index = -1;

			for(int o = 0; o < optimal_paths.size(); o++){
				float min_distance_selected = 1e9;
				
				for(size_t s = 0;  s < selected_paths.size(); s++){
					float dist = calculate_L1_distance(optimal_paths[o], selected_paths[s]);
					if(dist < min_distance_selected){
						min_distance_selected = dist;
					}
				}

				if(min_distance_selected > best_max_min_distance){
					best_max_min_distance = min_distance_selected;
					best_candidate_index = o;
				}
			}

			selected_paths.push_back(optimal_paths[best_candidate_index]);
			optimal_paths.erase(optimal_paths.begin() + best_candidate_index);
		}
	}

	// If we dont have enough N
	if(selected_paths.size() < nb_path_to_select){
		for(int i = 0; i < sol.inst.Gamma+1; i++){	// Sub-optimal paths
			if(pi_value[sol.inst.T][i] >= sub_OPT && abs(pi_value[sol.inst.T][i] - ub_cost) >= eps){
				extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, candidates_path, limit_number_paths, eps);
			}
		}


		// MaxMin L1 distance on the sub-optimal set
		if(!candidates_path.empty()){
			while(selected_paths.size() < nb_path_to_select && !candidates_path.empty()){
				float best_max_min_distance = -1.0;
				int best_candidate_index = -1;
				
				for(size_t c = 0; c < candidates_path.size(); c++){
					float min_distance_selected = 1e9;

					for(size_t s = 0; s < selected_paths.size(); s++){
						float dist = calculate_L1_distance(candidates_path[c], selected_paths[s]);
						if(dist < min_distance_selected){
							min_distance_selected = dist;
						}
					}

					if(min_distance_selected > best_max_min_distance){
						best_max_min_distance = min_distance_selected;
						best_candidate_index = c;
					}
				}

				selected_paths.push_back(candidates_path[best_candidate_index]);
				candidates_path.erase(candidates_path.begin() + best_candidate_index);
			}
		}
	}

	// Reset arcbool
	for(int t = 0; t < sol.inst.T+2; t++){
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j][0] = 0;
				arcbool[t][i][j][1] = 0;
			}
		}
	}

	// Flip only the arcs that are on the paths
	for(size_t p = 0; p < selected_paths.size(); p++){
		for(size_t a = 0; a < selected_paths[p].size(); a++){
			Arc_Decision arc = selected_paths[p][a];
			arcbool[arc.t][arc.i][arc.j][arc.type] = 1;
		}

		// Reconnected the end of the path to node T+1 like in the previous code
		int final_budget = selected_paths[p].back().j;
		arcbool[sol.inst.T+1][final_budget][0][0] = 1;
		arcbool[sol.inst.T+1][final_budget][0][1] = 1;
	}

	if(use_graph_export){
        export_budget_graph_json(sol, pi_value, costs, arcbool);
    }

	// Display the subgraph
	// cout<<"subgraph:"<<endl;
	// stringstream bufft;
	// for(int t = 0; t<sol.inst.T+2; t++){
	// 	bufft<<t;
	// 	bufft<<" ";
	// }
	// cout<<bufft.str()<<endl;
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	string buff = "";
	// 	for(int t = 0; t<sol.inst.T+2; t++){
	// 		buff += BoolToString(pi_subopt_bool[t][i])+" " ;
	// 	}
	// 	cout<<buff<<endl;
	// }

	// cout<<endl;
	// //display the subgraph
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	stringstream buff;
	// 	for(int t = 0; t<sol.inst.T+2; t++){
	// 		buff<< " ";
	// 		buff<<pi_value[t][i];
	// 	}
	// 	cout<<buff.str()<<endl;
	// }

	return arcbool;
}

Benders_Result KC_benders_Main(Instance inst, float approx_coeff, KC_Method method, bool use_graph_export, int limit_number_paths, int nb_path_to_select, float eps, int max_iter, int max_time_s, float p_few){
	auto start = high_resolution_clock::now();
	Solution sol;
	Solution new_sol;
	Solution_ADV sol_adv;
	float total_time_master = 0.0;
	float total_time_subproblem = 0.0;
	float proc_time;
	int i = 0;
	bool stopCriterion = false;
	float ub_cost;
	vector<vector<vector<vector<int> > > > arcsol = init_graph(inst);
	vector<vector<vector<vector<int> > > > arcsol_new;

	auto start_m = high_resolution_clock::now();
	sol = KC_benders_Master(inst, arcsol);
	auto stop_m = high_resolution_clock::now();
	total_time_master += duration_cast<microseconds>(stop_m - start_m).count() * 1e-6;	

	while(!stopCriterion){
		auto start_s = high_resolution_clock::now();

		switch(method){
			case KC_Method::KC:
				arcsol_new = KC_benders_Subproblem(sol, approx_coeff, use_graph_export, ub_cost, eps, p_few);
				break;

			case KC_Method::RDK:
				arcsol_new = KC_benders_Subproblem_RDK(sol, approx_coeff, nb_path_to_select, use_graph_export, limit_number_paths, ub_cost, eps);
				break;

			case KC_Method::RDKL:
				arcsol_new = KC_benders_Subproblem_RDKL(sol, approx_coeff, nb_path_to_select, use_graph_export, limit_number_paths, ub_cost, eps);
				break;

			case KC_Method::Unique:
				arcsol_new = KC_benders_Subproblem_Unique(sol, eps, ub_cost);
				break;
			
			case KC_Method::UniqueDual:
				arcsol_new = KC_benders_Subproblem_Unique_Dual(sol, eps, ub_cost);
				break;

			case KC_Method::HOG:
				arcsol_new = KC_benders_Subproblem_OG(sol, approx_coeff, nb_path_to_select, use_graph_export, limit_number_paths, ub_cost, eps);
				break;

			case KC_Method::HOGL:
				arcsol_new = KC_benders_Subproblem_OGL(sol, approx_coeff, nb_path_to_select, use_graph_export, limit_number_paths, ub_cost, eps);
				break;
		}

		auto stop_s = high_resolution_clock::now();
		total_time_subproblem += duration_cast<microseconds>(stop_s - start_s).count() * 1e-6;

		// If the cost increased, the opponent has found a computer breach and we continue in the loop
		if(ub_cost <= sol.obj_val + eps){
			stopCriterion = true;
			break;
		}

		arcsol = merge_budget_graph(arcsol, arcsol_new);	// Merge the worst solution with the current solution
		auto start_m2 = high_resolution_clock::now();
		new_sol = KC_benders_Master(inst, arcsol);			// Proposes a new solution according to the merge
		auto stop_m2 = high_resolution_clock::now();
		total_time_master += duration_cast<microseconds>(stop_m2 - start_m2).count() * 1e-6;

		auto current_time = high_resolution_clock::now();
		long long elapsed_s = duration_cast<seconds>(current_time - start).count();

		if(elapsed_s >= max_time_s){
			cout << "TIMEOUT: benders stop." << endl;
			break;
		}

		if(i >= max_iter){
			cout << "MAXITER: benders stop." << endl;
			break;
		}

		i++;
		
		// display_vector_float(new_sol.Xt);

		//cout<<"new sol value KC_benders_main: " << new_sol.obj_val << endl;

		//cout<<"============ "<< new_sol.obj_val << " " << sol.obj_val<<endl;
		
		//cout << sol.obj_val << endl;
		
		// cout<<"worst case: ";
		// display_vector_float(sol_adv.Dt);
		
		cout << "Iteration " << i << " - Master Obj (LB): " << sol.obj_val << " | Subproblem Cost (UB): " << ub_cost << endl;

		sol = new_sol;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	proc_time = duration.count() * 1e-6;

	return {i, proc_time,total_time_master, total_time_subproblem, sol.obj_val, sol};
}

// ========================================================================================================================================================================================================
// ========================================================================================================================================================================================================
// ================================================================================= Standard Benders Decomposition code ==================================================================================
// ========================================================================================================================================================================================================
// ========================================================================================================================================================================================================

Solution BA_benders_Master(Instance inst, vector<vector<float> > scenarios){
	Solution sol;
	sol.inst = inst;
	IloEnv env;
	IloModel model(env);

	// Vars
	IloNumVar z(env, -IloInfinity, IloInfinity);
	//IloNumVar z(env, -1000, 70);
	z.setName("z");
	IloNumVarArray  X(env, inst.T);
	IloArray<IloNumVarArray> s(env,scenarios.size());
	IloArray<IloNumVarArray> B(env,scenarios.size());
	IloArray<IloNumVarArray> I(env,scenarios.size());
	for(int o = 0; o < scenarios.size(); o++){
		s[o] = IloNumVarArray(env, inst.T);
		B[o] = IloNumVarArray(env, inst.T);
		I[o] = IloNumVarArray(env, inst.T);
		for(int t = 0; t < inst.T; t++){
			char name[80];
			s[o][t] = IloNumVar(env);
			sprintf(name,"s_%d_%d",o,t);
			s[o][t].setName(name);

			B[o][t] = IloNumVar(env);
			sprintf(name,"B_%d_%d",o,t);
			B[o][t].setName(name);
	
			I[o][t] = IloNumVar(env);
			sprintf(name,"I_%d_%d",o,t);
			I[o][t].setName(name);

			X[t] = IloNumVar(env, 0, IloInfinity,  IloNumVar::Float);
			sprintf(name,"X_%d",t);
			X[t].setName(name);
		}
	}

	// Consts
	for(int t = 1; t < inst.T+1; t++){
		model.add(X[t-1] <= inst.X[t-1]);	// inst.X[t] is the maximum possible production of the factory at time t
	}

	for(int t = 1; t < inst.T; t++){	
    	model.add(X[t] >= X[t-1]);
	}

	for(int t = 0; t < inst.T; t++){
		for(int o = 0; o < scenarios.size(); o++){
			model.add(B[o][t] - I[o][t] == scenarios[o][t] - X[t]); //(2)
			IloExpr expr(env);
			for(int i = 0; i <= t; i++){
				expr += s[o][i];
			}

			model.add(expr == scenarios[o][t]-B[o][t]);				//(3)
		}
	}

	for(int o = 0; o < scenarios.size(); o++){
		IloExpr expr(env);
		for(int t = 0; t < inst.T; t++){
			expr += (inst.cI*I[o][t] + inst.cB*B[o][t] - inst.bP*s[o][t]);
		}

		model.add(z >= expr);
	}

	// Objective value
	model.add(IloMinimize(env, z));

	// Solve
	IloCplex cplex(model);
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
	cplex.setOut(env.getNullStream());
    if(!cplex.solve()){
    	env.error() << "Failed to optimize LP." << endl;
    	throw(-1);
	}

	vector<float> Xt;
	Xt.resize(inst.T);
	for(int t = 0; t < inst.T; t++){
		Xt[t] = cplex.getValue(X[t]);
	}

	sol.Xt = Xt;
	sol.xt = cumulToStandard(sol.Xt);
	sol.obj_val = cplex.getObjValue();

	env.end();

	return sol;
}

// BAO
vector<Solution_ADV> BAO_benders_Subproblem_DFS(Solution sol, float approx_coeff, int nb_path_to_select, float eps){
	vector<vector<float> > pi_value;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol);

	pi_value.resize(sol.inst.T+2);
	for(int t = 0; t < sol.inst.T+2; t++){
		pi_value[t].resize(sol.inst.Gamma+1, -1e9);
	}
	
	// Dynamic prog. for longest path

	float tmp;
    pi_value[0][0] = 0;
    for(int t = 1; t < sol.inst.T+1; t++){
        for(int j = 0; j < sol.inst.Gamma+1; j++){
            tmp = pi_value[t-1][j] + costs[t][j][j][0]; 
            for(int i = 0; i <= j; i++){
                if(j <= i+sol.inst.deltat[t-1]){
                    if(pi_value[t-1][i] + costs[t][i][j][0] > tmp) tmp = pi_value[t-1][i]+costs[t][i][j][0];
                    if(pi_value[t-1][i] + costs[t][i][j][1] > tmp) tmp = pi_value[t-1][i]+costs[t][i][j][1];
                }
            }
            pi_value[t][j] = tmp;
        }
    }

    float ub_cost = pi_value[sol.inst.T][0];
    for(int i = 0; i < sol.inst.Gamma+1; i++){
        if(pi_value[sol.inst.T][i] > ub_cost){
            ub_cost = pi_value[sol.inst.T][i];
        }
    }

    pi_value[sol.inst.T+1][0] = ub_cost;
	
	// ========================== Now the backtrack
	
	vector<Path> collected_paths;
	Path current_path_buffer;

	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(abs(pi_value[sol.inst.T][i] - ub_cost) < eps){
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, collected_paths, nb_path_to_select, eps);
			if(collected_paths.size() >= nb_path_to_select) break;
		}
	}
	
	if(collected_paths.size() < nb_path_to_select){
		float sub_OPT;
		if(ub_cost >= 0){
			sub_OPT = approx_coeff * ub_cost;
		} else{
			sub_OPT = (1-approx_coeff) * ub_cost + ub_cost;
		}

		for(int i = 0; i < sol.inst.Gamma+1; i++){
			if(pi_value[sol.inst.T][i] >= sub_OPT && abs(pi_value[sol.inst.T][i] - ub_cost) >= eps){
				extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, collected_paths, nb_path_to_select, eps);
				if(collected_paths.size() >= nb_path_to_select) break;
			}
		}
	}
	
	

	vector<Solution_ADV> adv_scenarios;

	for(size_t p = 0; p < collected_paths.size(); p++){
		Solution_ADV adv;
		adv.Dt.resize(sol.inst.T, 0.0);

		for(size_t a = 0; a < collected_paths[p].size(); a++){
			Arc_Decision arc = collected_paths[p][a];
			int t = arc.t;
			int i = arc.i;
			int j = arc.j;

			if(arc.type == 0){
				adv.Dt[t-1] = sol.inst.Dt[t-1] - (j-i);
			} else {
				adv.Dt[t-1] = sol.inst.Dt[t-1] + (j-i);
			}
		}
		
		adv_scenarios.push_back(adv);
	}

	
	return adv_scenarios;
}

Benders_Result BAO_benders_Main(Instance inst, float approx_coeff, int nb_path_to_select, float eps, int max_iter, int max_time_s){
    auto start = chrono::high_resolution_clock::now();
    Solution sol;
    float total_time_master = 0.0;
    float total_time_subproblem = 0.0;
    bool stopCriterion = false;
    float proc_time;
    int i = 1;
    vector<vector<float> > scenarios;
    
    scenarios.push_back(inst.Dt);

    auto start_m = chrono::high_resolution_clock::now();
    sol = BA_benders_Master(inst, scenarios);
    auto stop_m = chrono::high_resolution_clock::now();
    total_time_master += chrono::duration_cast<chrono::microseconds>(stop_m - start_m).count() * 1e-6;

    while(!stopCriterion){
        auto start_s = chrono::high_resolution_clock::now();
        vector<Solution_ADV> adv_scenarios = BAO_benders_Subproblem_DFS(sol, eps, nb_path_to_select, approx_coeff);
        auto stop_s = chrono::high_resolution_clock::now();
        total_time_subproblem += chrono::duration_cast<chrono::microseconds>(stop_s - start_s).count() * 1e-6;

        float upper_bound_cost = objective_value(sol, adv_scenarios[0].Dt);

        if(upper_bound_cost <= sol.obj_val + eps){
            stopCriterion = true;
            break;
        }

        for(const auto& adv : adv_scenarios){
            scenarios.push_back(adv.Dt);
        }

        auto start_m2 = chrono::high_resolution_clock::now();
        sol = BA_benders_Master(inst, scenarios); //[cite: 2]
        auto stop_m2 = chrono::high_resolution_clock::now();
        total_time_master += chrono::duration_cast<chrono::microseconds>(stop_m2 - start_m2).count() * 1e-6;

        auto current_time = chrono::high_resolution_clock::now();
        long long elapsed_s = chrono::duration_cast<chrono::seconds>(current_time - start).count();

        if(elapsed_s >= max_time_s){
            cout << "TIMEOUT: benders stop." << endl;
            break;
        }

        if(i >= max_iter){
            cout << "MAXITER: benders stop." << endl;
            break;
        }

        i++;
        cout << "Iteration " << i << " - Master Obj (LB): " << sol.obj_val << " | Subproblem Cost (UB): " << upper_bound_cost << endl;
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    proc_time = duration.count() * 1e-6;
    
    return {i, proc_time, total_time_master, total_time_subproblem, sol.obj_val, sol};
}

// BA
Solution_ADV BA_benders_Subproblem_DP(Solution sol, float eps){
	Solution_ADV sol_adv;
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs of the budget graph

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t < sol.inst.T+2; t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i < sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j < sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;
	for(int t=1; t < sol.inst.T+1; t++){
		for(int j = 0; j < sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j]+costs[t][j][j][0];	// Init of pi_value[t][j]
			for(int i = 0; i <= j; i++){
				if(j <= i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i]+costs[t][i][j][0] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][0];
					}

					if(pi_value[t-1][i]+costs[t][i][j][1] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					} 
				}
			}

			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] > tmp){
			tmp = pi_value[sol.inst.T][i];
		}
	}

	pi_value[sol.inst.T+1][0] = tmp;
	
	// ========================== Now the backtrack

	vector<float> scenario;
	scenario.resize(sol.inst.T);

	// t = T+1
	pi_subopt_bool[sol.inst.T+1][0] = true;
	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(abs(pi_value[sol.inst.T][i] - pi_value[sol.inst.T+1][0]) < eps){
			arcbool[sol.inst.T+1][i][0][0] = 1;
			arcbool[sol.inst.T+1][i][0][1] = 1;
			pi_subopt_bool[sol.inst.T][i] = true;
		}
	}

	for(int t=sol.inst.T; t > 0; t--){
		for(int j = 0; j < sol.inst.Gamma+1; j++){
			for(int i = 0; i <= j; i++){
				if(pi_subopt_bool[t][j] and j <= i+sol.inst.deltat[t-1] and (t != 1 or i == 0)){ // last and is specific for first layer of the graph
					if(abs(pi_value[t][j] - (pi_value[t-1][i]+costs[t][i][j][0])) < eps){
						arcbool[t][i][j][0] = 1;
						pi_subopt_bool[t-1][i] = true;
						scenario[t-1] = sol.inst.Dt[t-1] - (j-i);		
						break;
					}

					if(abs(pi_value[t][j] - (pi_value[t-1][i]+costs[t][i][j][1])) < eps){
						arcbool[t][i][j][1] = 1;
						pi_subopt_bool[t-1][i] = true;
						scenario[t-1] = sol.inst.Dt[t-1] + (j-i);
						break;
					}
				}
			}
		}
	}

	// cout << "================= BEGIN TEST ==================" << endl;

	// Solution test = KC_benders_Master(sol.inst, arcbool);

	// cout<<"new sol (with Graph LP): ";
	// display_vector_float(test.Xt);
	// cout<<"new sol value: "<<test.obj_val<<endl;

	// cout << "================= END TEST ==================" << endl;


	// Display the subgraph
	// cout<<"subgraph:"<<endl;
	// stringstream bufft;
	// for(int t = 0; t<sol.inst.T+2; t++){
	// 		bufft<<t;
	// 		bufft<<" ";
	// }

	// cout<<bufft.str()<<endl;
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 		string buff = "";
	// 		for(int t = 0; t<sol.inst.T+2; t++){
	// 			buff += BoolToString(pi_subopt_bool[t][i])+" " ;
	// 		}
	// 	
	//		cout<<buff<<endl;
	// }

	// cout<<endl;
	// Display the subgraph
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	stringstream buff;
	// 	for(int t = 0; t<sol.inst.T+2; t++){
	// 		buff<< " ";
	// 		buff<<pi_value[t][i];
	// 	}

	// 	cout<<buff.str()<<endl;
	// }

	sol_adv.Dt = scenario;

	return sol_adv;
}

Benders_Result BA_benders_Main(Instance inst, float eps, int max_iter, int max_time_s){
	auto start = high_resolution_clock::now();
	Solution sol;
	Solution new_sol;
	Solution_ADV sol_adv;
	float total_time_master = 0.0;
	float total_time_subproblem = 0.0;
	bool stopCriterion = false;
	float proc_time;
	int i = 1;
	vector<vector<float> > scenarios;
	scenarios.resize(0);

	scenarios.push_back(inst.Dt);						// Initialization with the nominal scenario

	auto start_m = high_resolution_clock::now();
	sol = BA_benders_Master(inst, scenarios);			// Proposal for an initial plan X (initial lower bound)
	auto stop_m = high_resolution_clock::now();
	total_time_master += duration_cast<microseconds>(stop_m - start_m).count() * 1e-6;

	while(!stopCriterion){
		auto start_s = high_resolution_clock::now();
		sol_adv = BA_benders_Subproblem_DP(sol, eps);	// The subproblem seeks the worst-case scenario D' against plan X (solution).
		auto stop_s = high_resolution_clock::now();
		total_time_subproblem += duration_cast<microseconds>(stop_s - start_s).count() * 1e-6;

		float upper_bound_cost = objective_value(sol, sol_adv.Dt);

		if(upper_bound_cost <= sol.obj_val + eps){
			stopCriterion = true;
			break;
		}

		scenarios.push_back(sol_adv.Dt);

		//display_vector_float(sol_adv.Dt);
		//display_vector_float(sol_adv.Dt);

		auto start_m2 = high_resolution_clock::now();
		sol = BA_benders_Master(inst, scenarios);
		auto stop_m2 = high_resolution_clock::now();
		total_time_master += duration_cast<microseconds>(stop_m2 - start_m2).count() * 1e-6;
		auto current_time = high_resolution_clock::now();
        long long elapsed_s = duration_cast<seconds>(current_time - start).count();

        if(elapsed_s >= max_time_s){
            cout << "TIMEOUT: benders stop." << endl;
            break;
        }

        if(i >= max_iter){
            cout << "MAXITER: benders stop." << endl;
            break;
        }

		i++;
			
		cout << "Iteration " << i << " - Master Obj (LB): " << sol.obj_val << " | Subproblem Cost (UB): " << upper_bound_cost << endl;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	proc_time = duration.count() * 1e-6;
	
	return {i, proc_time, total_time_master, total_time_subproblem, sol.obj_val, sol};
}

// ========================================================================================================================================================================================================
// ========================================================================================================================================================================================================
// ================================================================================================= Main =================================================================================================
// ========================================================================================================================================================================================================
// ========================================================================================================================================================================================================

vector<string> list_dir(const char *path) {
vector<string> allfile;
   struct dirent *entry;
   DIR *dir = opendir(path);
   
   if (dir == NULL) {
      return allfile;
   }
   
   while((entry = readdir(dir)) != NULL){
		allfile.push_back(entry -> d_name);
   }

   closedir(dir);
   
   return allfile;
}

int main(int argc, const char* argv[]){
	Benders_Result benders_sol_BA;
	Benders_Result benders_sol_BAO;
	Benders_Result benders_sol_KCA;
	Benders_Result benders_sol_KCF;
	Benders_Result benders_sol_KCRDK;
	Benders_Result benders_sol_KCRDKL;
	Benders_Result benders_sol_KCU;
	Benders_Result benders_sol_KCUD;
	Benders_Result benders_sol_KCOG;
	Benders_Result benders_sol_KCOGL;
	// Output parameters
	float approx_coeff;
	// Monte Carlo parameter
	bool use_monte_carlo = false;
	// Simulation parameters
	int limit_number_paths = 600;		// Used for the DFS algo
	int nb_path_to_select = 2;			// Number of paths the subprobleme give to the master at each iteration (it's the upperbound like the max of the parameter and the number found)
	float eps = 1e-1;
	bool use_graph_export = false;		// To be corrected before use
	bool use_result_export = true;		// True if you want to export the results
	float max_iter = 100;				// Security
	float max_time_s = 3600;
	float p_few = 1/52;					// Probabilitie for SelectFew heuristique
	string validation_status;
	// Read instances randomized parameters
	float read_instance_rd_lb = 0.4;	// The production plan will be between lb% and ub% of the cumulative demand
	float read_instance_rd_ub = 0.8;
	float adv_margin          = 2.0;	// Margin allowed to the opponent in the calculation of deltats: adv_margin*Gamma (=1: no marge, =2: a lot of)
	// Création of the folder architecture
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	ostringstream oss_exp;
	oss_exp << "/result_" << put_time(&tm, "%Y-%m-%d_%H%M") << "_n=" << nb_path_to_select << "_l=" << limit_number_paths << "_pfew=" << p_few << "_advmargin=" << adv_margin;
	std::string experience_name = oss_exp.str();

	ostringstream oss_folder;
	oss_folder << "./results" << experience_name;
	std::string folder_path = oss_folder.str();

	if(use_result_export){
		if(fs::create_directories(folder_path)){
			std::cout << "File '" << folder_path << "' has been successfully created" << std::endl;
		} else{
			std::cout << "File '" << folder_path << "' already exists" << std::endl;
		}

		cout << "Recording of results in: " << folder_path << endl;
	}

	ostringstream oss_BA_BAO;
	ostringstream oss_BA_KCF;
	ostringstream oss_BA_KCRDK;
	ostringstream oss_BA_KCRDKL;
	ostringstream oss_BA_KCU;
	ostringstream oss_KCU_KCUD;
	ostringstream oss_BA_OG;
	ostringstream oss_BA_OGL;
	ostringstream oss_KCA_KCF;
	ostringstream oss_validation;
	ostringstream oss_time_m_s;
	ostringstream oss_stats;
	oss_BA_BAO		<< folder_path << experience_name << "_BA_BAO.csv";
	oss_BA_KCF 	   	<< folder_path << experience_name << "_BA_KCF.csv";
	oss_BA_KCRDK    << folder_path << experience_name << "_BA_KCRDK.csv";
	oss_BA_KCRDKL   << folder_path << experience_name << "_BA_KCRDKL.csv";
	oss_BA_KCU		<< folder_path << experience_name << "_BA_KCU.csv";
	oss_KCU_KCUD	<< folder_path << experience_name << "_KCU_KCUD.csv";
	oss_BA_OG 	  	<< folder_path << experience_name << "_BA_OG.csv";
	oss_BA_OGL		<< folder_path << experience_name << "_BA_OGL.csv";
	oss_KCA_KCF  	<< folder_path << experience_name << "_KCA_KCF.csv";
	oss_validation  << folder_path << experience_name << "_validation.txt";
	oss_time_m_s    << folder_path << experience_name << "_time_m_s.csv";
	oss_stats		<< folder_path << experience_name << "_stats.csv";
	ofstream output_BA_BAO(oss_BA_BAO.str());
	ofstream output_BA_KCF(oss_BA_KCF.str());
	ofstream output_BA_KCRDK(oss_BA_KCRDK.str());
	ofstream output_BA_KCRDKL(oss_BA_KCRDKL.str());
	ofstream output_BA_KCU(oss_BA_KCU.str());
	ofstream output_BA_KCUD(oss_KCU_KCUD.str());
	ofstream output_BA_OG(oss_BA_OG.str());
	ofstream output_BA_OGL(oss_BA_OGL.str());
	ofstream output_KCA_KCF(oss_KCA_KCF.str());
	ofstream output_validation;
	ofstream output_time_m_s(oss_time_m_s.str());
	ofstream output_stats(oss_stats.str());

	if(use_result_export){
		output_validation.open(oss_validation.str());
		output_validation << "Fichier 			| Gamma 	| tau 	| 	Validation\n";
		output_validation << "------------------------------------------------------------\n";
		output_time_m_s << "Gamma, tau, nb_path_to_select, limit_number_paths, time_master_BA, time_subproblem_BA, time_master_BAO, time_subproblem_BAO, time_master_KCA, time_subproblem_KCA, time_master_KCF, time_subproblem_KCF, time_master_KCRDK, time_subproblem_KCRDK, time_master_KCRDKL, time_subproblem_KCRDKL, time_master_KCU, time_subproblem_KCU, time_master_KCUD, time_subproblem_KCUD, time_master_KCOG, time_subproblem_KCOG, time_master_KCOGL, time_subproblem_KCOGL\n";
	}

	vector<string> file_list;
	bool toy_instances = true;
	if(toy_instances){
		file_list = list_dir("./toy_instances/");
	} else{
		file_list = list_dir("./hand_benders_instances/parsed_instances/");
	}
	
	int total_files = file_list.size();
	int nbInst = 0;
	for(int i = 0; i < total_files; i++){
		if(file_list[i] != "." && file_list[i] != ".."){
			nbInst++;
		}
	}

	// Security
	if(nbInst == 0){
		cerr << "Error: No valid instance file found. Check the path." << endl;
		return -1;
	} else{
		cout << "Success: " << nbInst << " files found in the instances folder." << endl;
	}

	Instance inst;
	string filename;

	for(int i = 0; i < total_files; i++ ){
		if(file_list[i] == "." || file_list[i] == "..") continue;

		if(toy_instances){
			filename = "toy_instances/" + file_list[i];
		} else{
			filename = "hand_benders_instances/parsed_instances/" + file_list[i];
		}

		for(int Gamma = 1; Gamma < 110; Gamma += 20){
			//int Gamma = 101;
			//for(int tau = 80; tau < 101; tau += 10){
				int tau = 90;
				approx_coeff = float(tau)/100;
				//for(int nb_path_to_select = 1; nb_path_to_select < 12; nb_path_to_select += 2){
					//for(int limit_number_paths = 100; limit_number_paths < 1001; limit_number_paths += 100){
						if(toy_instances){
							inst = read_instance_py(filename, Gamma, adv_margin);
							// inst = read_instance_randomized(filename, Gamma, read_instance_rd_lb, read_instance_rd_ub, adv_margin);	// If you don't want to use py script
						} else {
							inst = read_hand_instance(filename, 2);
						}
											
						cout << "\n" << filename << " " << Gamma << " " << tau << " " << endl;

						// BA		
						benders_sol_BA = BA_benders_Main(inst, eps, max_iter, max_time_s);
						cout << "\nBA-----done (Obj:" << benders_sol_BA.obj_value << ")" << endl;
					
						// BAO
						benders_sol_BAO = BAO_benders_Main(inst, approx_coeff, nb_path_to_select, eps, max_iter, max_time_s);
						cout << "\nBAO----done (Obj:" << benders_sol_BAO.obj_value << ")" << endl;
						
						// KC All (BAEA)
						benders_sol_KCA = KC_benders_Main(inst, approx_coeff, KC_Method::KC, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, 1); // First bool is to use HOG, the other is to use the KC random K method
						cout << "\nKCA----done (Obj:" << benders_sol_KCA.obj_value << ")" <<endl;

						// KC Few (BAEF)
						benders_sol_KCF = KC_benders_Main(inst, approx_coeff, KC_Method::KC, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, p_few);
						cout << "\nKCF----done (Obj:" << benders_sol_KCF.obj_value << ")" <<endl;

						// KCRDK (RanDom K)
						benders_sol_KCRDK = KC_benders_Main(inst, approx_coeff, KC_Method::RDK, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, p_few);
						cout << "\nKCRDK--done (Obj:" << benders_sol_KCRDK.obj_value << ")" <<endl;

						// KCRDKL (KCRDK Lexicographical)
						benders_sol_KCRDKL = KC_benders_Main(inst, approx_coeff, KC_Method::RDKL, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, p_few);
						cout << "\nKCRDKL-done (Obj:" << benders_sol_KCRDKL.obj_value << ")" <<endl;
						
						// KCU	(KC Unique)
						benders_sol_KCU = KC_benders_Main(inst, approx_coeff, KC_Method::Unique, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, p_few); //KCU_benders_Main(inst, eps, max_iter, max_time_s);
						cout << "\nKCU----done (Obj:" << benders_sol_KCU.obj_value << ")" << endl;

						// KCUD (KC with two optimal paths)
						benders_sol_KCUD = KC_benders_Main(inst, approx_coeff, KC_Method::UniqueDual, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, p_few);
						cout << "\nKCUD---done (Obj:" << benders_sol_KCUD.obj_value << ")" << endl;

						// KCOG
						benders_sol_KCOG = KC_benders_Main(inst, approx_coeff, KC_Method::HOG, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, p_few);
						cout << "\nKCOG---done (Obj:" << benders_sol_KCOG.obj_value << ")"<< endl;

						// KCOGL (OG Lexicographical)
						benders_sol_KCOGL = KC_benders_Main(inst, approx_coeff, KC_Method::HOGL, use_graph_export, limit_number_paths, nb_path_to_select, eps, max_iter, max_time_s, p_few);
						cout << "\nKCOGL--done (Obj:" << benders_sol_KCOGL.obj_value << ")"<< endl;

						// Quality control of the solution
						if(abs(benders_sol_BA.obj_value - benders_sol_BAO.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCA.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCF.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCRDK.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCRDKL.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCU.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCUD.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCOG.obj_value) > eps || abs(benders_sol_BA.obj_value - benders_sol_KCOGL.obj_value) > eps){
							cout << "\nDegraded quality" << endl;
							validation_status = "Degraded quality";
						} else{
							cout << "\nValid quality" << endl;
							validation_status = "Valid quality";
						}

						// We export the exact time for each instance and each parameters
						if(use_result_export){
							// Validation file
							output_validation << filename << " | " << Gamma << " | " << tau << " | " << nb_path_to_select << " | " << limit_number_paths << " | " << validation_status << endl;
						
							// BA with BAO
							output_BA_BAO	<< "BA," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BA.iter << "," 
											<< benders_sol_BA.time << ","
											<< benders_sol_BA.time_master << ","
											<< benders_sol_BA.time_subproblem << endl;
							
							output_BA_BAO 	<< "BAO," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BAO.iter << "," 
											<< benders_sol_BAO.time << ","
											<< benders_sol_BAO.time_master << ","
											<< benders_sol_BAO.time_subproblem << endl;

							// BA with KC standard 
							output_BA_KCF 	<< "BA," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BA.iter << "," 
											<< benders_sol_BA.time << ","
											<< benders_sol_BA.time_master << ","
											<< benders_sol_BA.time_subproblem << endl;

							output_BA_KCF 	<< "KCF," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCF.iter << ","
											<< benders_sol_KCF.time << ","
											<< benders_sol_KCF.time_master << ","
											<< benders_sol_KCF.time_subproblem << endl;
							
							// BA with KCRDK
							output_BA_KCRDK << "BA,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BA.iter << ","
											<< benders_sol_BA.time << ","
											<< benders_sol_BA.time_master << ","
											<< benders_sol_BA.time_subproblem << endl;

							output_BA_KCRDK << "KCRDK,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCRDK.iter << ","
											<< benders_sol_KCRDK.time << ","
											<< benders_sol_KCRDK.time_master << ","
											<< benders_sol_KCRDK.time_subproblem << endl;

							// BA with KCRDKL
							output_BA_KCRDKL	<< "BA,"
												<< Gamma << ","
												<< tau << ","
												<< nb_path_to_select << ","
												<< limit_number_paths << ","
												<< benders_sol_BA.iter << ","
												<< benders_sol_BA.time << ","
												<< benders_sol_BA.time_master << ","
												<< benders_sol_BA.time_subproblem << endl;
										
							output_BA_KCRDKL	<< "KCRDKL,"
												<< Gamma << ","
												<< tau << ","
												<< nb_path_to_select << ","
												<< limit_number_paths << ","
												<< benders_sol_KCRDKL.iter << ","
												<< benders_sol_KCRDKL.time << ","
												<< benders_sol_KCRDKL.time_master << ","
												<< benders_sol_KCRDKL.time_subproblem << endl;

							// BA with KCU
							output_BA_KCU	<< "BA," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BA.iter << "," 
											<< benders_sol_BA.time << ","
											<< benders_sol_BA.time_master << ","
											<< benders_sol_BA.time_subproblem << endl;

							output_BA_KCU 	<< "KCU," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCU.iter << "," 
											<< benders_sol_KCU.time << ","
											<< benders_sol_KCU.time_master << ","
											<< benders_sol_KCU.time_subproblem << endl;

							// KCU with KCUD
							output_BA_KCUD 	<< "KCU," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCU.iter << "," 
											<< benders_sol_KCU.time << ","
											<< benders_sol_KCU.time_master << ","
											<< benders_sol_KCU.time_subproblem << endl;

							output_BA_KCUD	<< "KCUD,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCUD.iter << ","
											<< benders_sol_KCUD.time << ","
											<< benders_sol_KCUD.time_master << ","
											<< benders_sol_KCUD.time_subproblem << endl;

							// BA with HOG
							output_BA_OG	<< "BA," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BA.iter << "," 
											<< benders_sol_BA.time << ","
											<< benders_sol_BA.time_master << ","
											<< benders_sol_BA.time_subproblem << endl;
											
							output_BA_OG 	<< "KCOG," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCOG.iter << "," 
											<< benders_sol_KCOG.time << ","
											<< benders_sol_KCOG.time_master << ","
											<< benders_sol_KCOG.time_subproblem << endl;

							// BA with OGL
							output_BA_OGL	<< "BA," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BA.iter << "," 
											<< benders_sol_BA.time << ","
											<< benders_sol_BA.time_master << ","
											<< benders_sol_BA.time_subproblem << endl;
											
							output_BA_OGL 	<< "KCOGL," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCOGL.iter << "," 
											<< benders_sol_KCOGL.time << ","
											<< benders_sol_KCOGL.time_master << ","
											<< benders_sol_KCOGL.time_subproblem << endl;

							// KCA with KCF
							output_KCA_KCF	<< "KCA," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCA.iter << "," 
											<< benders_sol_KCA.time << ","
											<< benders_sol_KCA.time_master << ","
											<< benders_sol_KCA.time_subproblem << endl;
							
							output_KCA_KCF	<< "KCF," 
											<< Gamma << "," 
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCF.iter << "," 
											<< benders_sol_KCF.time << ","
											<< benders_sol_KCF.time_master << ","
											<< benders_sol_KCF.time_subproblem << endl;

							// Time master subproblem
							output_time_m_s << Gamma 								<< ","
											<< tau 									<< ","
											<< nb_path_to_select 					<< ","
											<< limit_number_paths 					<< ","
											<< benders_sol_BA.time_master 			<< ","
											<< benders_sol_BA.time_subproblem 		<< ","
											<< benders_sol_BAO.time_master 			<< ","
											<< benders_sol_BAO.time_subproblem 		<< ","
											<< benders_sol_KCA.time_master 			<< ","
											<< benders_sol_KCA.time_subproblem 		<< ","
											<< benders_sol_KCF.time_master 			<< ","
											<< benders_sol_KCF.time_subproblem 		<< ","
											<< benders_sol_KCRDK.time_master 		<< ","
											<< benders_sol_KCRDK.time_subproblem 	<< ","
											<< benders_sol_KCRDKL.time_master 		<< ","
											<< benders_sol_KCRDKL.time_subproblem 	<< ","
											<< benders_sol_KCU.time_master 			<< ","
											<< benders_sol_KCU.time_subproblem 		<< ","
											<< benders_sol_KCUD.time_master 		<< ","
											<< benders_sol_KCUD.time_subproblem 	<< ","
											<< benders_sol_KCOG.time_master 		<< ","
											<< benders_sol_KCOG.time_subproblem 	<< ","
											<< benders_sol_KCOGL.time_master 		<< ","
											<< benders_sol_KCOGL.time_subproblem 	<< endl;

							// Stats
							output_stats	<< "BA,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BA.time << ","
											<< benders_sol_BA.iter << endl;

							output_stats	<< "BAO,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_BAO.time << ","
											<< benders_sol_BAO.iter << endl;

							output_stats 	<< "KCA,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCA.time << ","
											<< benders_sol_KCA.iter << endl;

							output_stats 	<< "KCF,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCF.time << ","
											<< benders_sol_KCF.iter << endl;

							output_stats 	<< "KCRDK,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCRDK.time << ","
											<< benders_sol_KCRDK.iter << endl;

							output_stats	<< "KCRDKL,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCRDKL.time << ","
											<< benders_sol_KCRDKL.iter << endl;

							output_stats 	<< "KCU,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCU.time << ","
											<< benders_sol_KCU.iter << endl;

							output_stats 	<< "KCUD,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCUD.time << ","
											<< benders_sol_KCUD.iter << endl;

							output_stats 	<< "KCOG,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCOG.time << ","
											<< benders_sol_KCOG.iter << endl;
							
							output_stats 	<< "KCOGL,"
											<< Gamma << ","
											<< tau << ","
											<< nb_path_to_select << ","
											<< limit_number_paths << ","
											<< benders_sol_KCOGL.time << ","
											<< benders_sol_KCOGL.iter << endl;
						}
					//}	// limit_dfs
				//}		// number_path_to_take
			//}			// tau
		}
	}

	if(use_monte_carlo){
		MonteCarlo_Result mc_robuste;
		int num_scenarios = 10000;

		// ==================== BA ====================
		mc_robuste = run_monte_carlo(benders_sol_BA.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the BA plan on " << num_scenarios 				  << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev 									<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_BA.obj_value 							<< endl;

		// ==================== BAO====================
		mc_robuste = run_monte_carlo(benders_sol_BAO.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the BAO plan on " << num_scenarios 		      << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev 									<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_BAO.obj_value 							<< endl;
		
		// ==================== BAEA ====================
		mc_robuste = run_monte_carlo(benders_sol_KCA.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the BAEA plan on " << num_scenarios 			  << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev 									<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_KCA.obj_value 							<< endl;

		// ==================== BAEF ====================
		mc_robuste = run_monte_carlo(benders_sol_KCF.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the BAEF plan on " << num_scenarios			  << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev 									<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_KCF.obj_value 							<< endl;
	
		// ==================== KCU ====================
		mc_robuste = run_monte_carlo(benders_sol_KCU.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the KCU plan on " << num_scenarios 			  << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev	 								<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_KCU.obj_value 							<< endl;

		// ==================== KCUD ====================
		mc_robuste = run_monte_carlo(benders_sol_KCUD.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the KCUD plan on " << num_scenarios 			  << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev 									<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_KCUD.obj_value 							<< endl;

		// ==================== KCRDKL ====================
		mc_robuste = run_monte_carlo(benders_sol_KCRDKL.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the KCRDKL plan on " << num_scenarios 			  << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev 									<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_KCRDKL.obj_value 						<< endl;

		// ==================== KCOGL ====================
		mc_robuste = run_monte_carlo(benders_sol_KCOGL.final_solution, num_scenarios);
		cout << "\n--- Monte Carlo evaluation of the KCOGL plan on " << num_scenarios 			  << " demands ---" << endl;
		cout << "Expected average cost                : " 	<< mc_robuste.mean_cost 								<< endl;
		cout << "Standard deviation                   : " 	<< mc_robuste.std_dev 									<< endl;
		cout << "CI (95%)                             : [" << mc_robuste.ci_lower << "; " << mc_robuste.ci_upper 	<< "]" << endl;
		cout << "Worst case simulated                 : " 	<< mc_robuste.worst_case_simulated 						<< endl;
		cout << "Theoretical worst-case cost (Benders): " 	<< benders_sol_KCOGL.obj_value	 						<< endl;
	}


	if(use_result_export){
		output_BA_BAO.close();
		output_BA_KCF.close();
		output_BA_KCRDK.close();
		output_BA_KCRDKL.close();
		output_BA_KCU.close();
		output_BA_KCUD.close();
		output_BA_OG.close();
		output_BA_KCRDKL.close();
		output_KCA_KCF.close();
		output_validation.close();
		output_time_m_s.close();
		output_stats.close();
	}
	
	cout << "========== END OF THE PROGRAM ==========" << endl;
}