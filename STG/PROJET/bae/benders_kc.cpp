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

#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

//=========================================== Structures

struct Instance
{
	int T, cI, cB, bP;
	float Gamma;	// Uncertainty budget
	vector<float> deltat;
	vector<float> Dt;
	vector<float> dt;
	vector<float> X;
};

struct Solution
{
	Instance inst;
	vector<float> Xt;
	vector<float> xt;
	float obj_val;
};

struct Instance_ADV
{
	int T, cI, cB, bP;
	float Gamma;
	vector<float> deltat;
	vector<float> Xt;
	vector<float> xt;
};

struct Solution_ADV
{
	Instance_ADV inst;
	vector<float> Dt;
	vector<float> dt;

};

// ====================================================================== IN PROGRESS ==============================================================================================================
// ====================================================================================================================================================================================
struct Benders_Result{
	int iter;
	float time;
	float obj_value;
};

struct Arc_Decision
{
	int t;	// Time
	int i;	// Budget at the start
	int j;	// Budget at the end
	int type;	// 0 or 1 (overstock / stockout)

	// Definition of two identical arcs
	bool operator == (const Arc_Decision& other) const{
		return (t == other.t && i == other.i && j == other.j && type == other.type);
	}
};

typedef vector<Arc_Decision> Path;

// =========================================== Calculate the Jaccard distance for the orthogonality heuristic

float calculate_jaccard_distance(const Path& pathA, const Path& pathB){
	int intersection_size = 0;

	for(size_t k = 0; k < pathA.size(); k++){
		if(pathA[k] == pathB[k]){
			intersection_size++;
		}
	}

	int union_size = pathA.size() + pathB.size() - intersection_size;

	// Calculate the Jaccard index
	float jaccard_similarity = (float)intersection_size / (float)union_size;

	// Return the dissimilarity
	return 1.0f - jaccard_similarity;
}

// =========================================== Calculate the Brays-Curtis distance for the orthogonality heuristic

float calculate_BC_distance(const Path& pathA, const Path& pathB){
	if(pathA.size() != pathB.size()){
		return 0.0f;
	}

	int sum_diff = 0;
	int sum_total = 0;

	for(size_t k = 0; k < pathA.size(); k++){
		int delta_A = pathA[k].j - pathB[k].i;
		int delta_B = pathB[k].j - pathA[k].i;

		sum_diff += abs(delta_A - delta_B);	// Manhattan local distance
		sum_total += (delta_A + delta_B);	// Total budget consumed by both routes
	}

	if(sum_total == 0){
		return 0.0f;
	}

	// Return the dissimilarity
	return static_cast<float>(sum_diff) / static_cast<float>(sum_total);
}

// =========================================== Recursive extraction of worst-case scenarios following a Depth-First Search

void extract_paths_dfs(
			int t,													
			int j,													// Current node in the backtrack
			const vector<vector<float> >& pi_value,					// Dynamic programming matrix
			const vector<vector<vector<vector<float> > > >& costs,	// Original costs
			const Solution& sol,
			Path& current_path,
			vector<Path>& all_paths,
			int limit_number_paths
			){
	
	if(all_paths.size() >= limit_number_paths){	// This solution can be a problem if we come in a choke point, maybe we can add a seed
		return;
	}

	// If we went back to the beginning we stop
	if(t == 0){
		Path reversed_path = current_path;	// Because we start at the end, we have to reverse the path of this branch
		reverse(reversed_path.begin(), reversed_path.end());
		all_paths.push_back(reversed_path);
		return;
	}

	// We are trying to figure out where we came from, its like : what was the budget i? to arrive at j at step t
	for(int i = 0; i <= j; i++){
		if(j <= i + sol.inst.deltat[t-1] && (t != 1 || i == 0)){
			// Check for arc of type 0
			if(abs(pi_value[t][j] - (pi_value[t-1][i] + costs[t][i][j][0])) < 1e-4){
				Arc_Decision arc = {t, i, j, 0};
				current_path.push_back(arc);
				extract_paths_dfs(t-1, i, pi_value, costs, sol, current_path, all_paths, limit_number_paths);
				current_path.pop_back();
			}

			// Check for arc of type 1
			if(abs(pi_value[t][j] - (pi_value[t-1][i] + costs[t][i][j][1]) < 1e-4)){
				Arc_Decision arc = {t, i, j, 1};
				current_path.push_back(arc);
				extract_paths_dfs(t-1, i, pi_value, costs, sol, current_path, all_paths, limit_number_paths);
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
		const vector<vector<vector<vector<int>>>>& arcbool
		){

	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	ostringstream oss;
	oss << "graph_visualization/benders_graph_output_" << put_time(&tm, "%Y-%m-%d_%H%M%S") << ".json";

	ofstream output(oss.str());

	if (!output.is_open()) {
		cerr << "ERREUR CRITIQUE : Le dossier 'graph_visualization' n'existe peut-être pas !" << endl;
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
                
                for(int type = 0; type < 2; type++) {
                    if (t == sol.inst.T + 1 && j != 0) continue;	// Do not export the arcs of the last layer if it is not the well
                    
                    if(!first_arc) output << ",\n";
                    float cost = (t == sol.inst.T+1) ? 0 : costs[t][i][j][type];
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

// ====================================================================================================================================================================================
// ====================================================================================================================================================================================

//=========================================== Misc.

// Displan an float vector, for debugging
void display_vector_float(vector<float> v){
	for (int i = 0; i < v.size(); ++i){
		cout<<v[i]<<" ";
	}

	cout<<endl;
}

// Displan an int vector, for debugging
void display_vector_int(vector<int> v){
	for (int i = 0; i < v.size(); ++i){
		cout<<v[i]<<" ";
	}

	cout<<endl;
}

string BoolToString(bool b){
  return b ? "1" : "0";
}

//=========================================== Instance related code

vector<float> standardToCumul(vector<float> data){
	vector<float> cumul;
	float tmp = 0;
	for(int i = 0; i<data.size(); i++){
		tmp += data[i];
		cumul.push_back(tmp);
	}

	return cumul;
}

vector<float> cumulToStandard(vector<float> data){
	vector<float> stand;
	stand.push_back(data[0]);
	for(int i = 1; i<data.size(); i++){
		stand.push_back(data[i]-data[i-1]);
	}

	return stand;
}

// Read an instance from psplib problem 58 pspInstance. Products are agregated.
Instance read_instance(string filename, int budget){
	Instance inst;

	int nbProd;
	int tmp;

	ifstream file(filename.c_str());
	if (!file){
		cout << "problem with file" << endl;
		exit(-1);
	} 

	file >> inst.T;
	file >> nbProd;
	inst.dt.resize(inst.T);
	for(int i = 0; i<nbProd; i++){
		for(int j = 0; j<inst.T; j++){
			file >> tmp;
			inst.dt[j] += tmp;
		}
	}
	
	inst.Dt = standardToCumul(inst.dt);

	inst.cI = 3; 	// Stock cost
	inst.cB = 6; 	// Backorder cost
	inst.bP = 10; 	// Selling price
	inst.Gamma = budget;
	inst.deltat.resize(inst.T);
	
	for(int t = 0; t<inst.T;t++){
		inst.deltat[t] = int(inst.Dt[t]/float(2));
	}

	file.close();

	return inst;
}

Instance read_instance_randomized(string filename, int budget){
	Instance inst;

	int nbProd;
	int tmp;

	ifstream file(filename.c_str());
	if (!file){
		cout << "problem with file" << endl;
		exit(-1);
	} 

	file >> inst.T;
	file >> nbProd;
	inst.dt.resize(inst.T);

	for(int i = 0; i<nbProd; i++){
		for(int j = 0; j<inst.T; j++){
			file >> tmp;
			inst.dt[j] += tmp + int(rand() % 2);	// Arbitrary to have nice instances
		}
	}

	inst.Dt = standardToCumul(inst.dt);

	
	inst.cB = rand() % 10 + 10;	// Backorder cost
	inst.cI = rand() % 10 + 10; // Stock cost
	inst.bP = rand() % 10 + 10; // Selling price
	inst.Gamma = budget;
	inst.Dt = standardToCumul(inst.dt);
	inst.deltat.resize(inst.T);

	// Display_vector_float(inst.Dt);
	for(int t = 0; t<inst.T;t++){
		// Recours temporaire (pas propre)
		if(inst.Dt[t] == 0){
			inst.deltat[t] = 0;
		}

		else if(t==0){
			inst.deltat[t] = rand()%(int(inst.Dt[t]));
			while(inst.deltat[t] > inst.dt[t+1]){
				inst.deltat[t] = rand()%(int(inst.Dt[t]));
			}
		}

		else if(t<inst.T-1){
			inst.deltat[t] = rand()%(int(inst.Dt[t]));
			while(inst.Dt[t-1] + inst.deltat[t-1] > inst.Dt[t] - inst.deltat[t] or inst.deltat[t] > inst.dt[t+1]){
				inst.deltat[t] = rand()%(int(inst.Dt[t]));
			}
		}

		else{
			inst.deltat[t] = rand()%(int(inst.Dt[t]));
			while(inst.Dt[t-1] + inst.deltat[t-1] > inst.Dt[t] - inst.deltat[t]){
				inst.deltat[t] = rand()%(int(inst.Dt[t]));
			}
		}
	}

	inst.X.resize(inst.T);
	
	// Recours temporaire (pas propre)
	for(int t = 0; t<inst.T;t++){
		if(inst.Dt[t]==0){
			inst.X[t] = 0;
		}
	
		else{
			inst.X[t] = int(rand()%(int(0.4*inst.Dt[t])+1) + 0.8*inst.Dt[t]);
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
	// costs[t][i][j][k]
	// t: time
	// i: uncertainty budget already consumed at the beginning of the period
	// j: uncertainty budget consumed at the end of the period
	// (i-j): delta_t deflection chosen by the opponent
	// k=2: to store the 2 sides of the maximum of the cost function
	vector<vector<vector<vector<float> > > > costs;

	costs.resize(sol.inst.T+2);

	for(int t = 1; t<sol.inst.T+2;t++){
		costs[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			costs[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j<sol.inst.Gamma+1; j++){
				costs[t][i][j].resize(2);
				if(j<=i+sol.inst.deltat[t-1] and j>=i){
					if(t<sol.inst.T){
																	// (sol.inst.Dt[t-1] - (j-i)): reduced demand
																	// (sol.inst.Dt[t-1] + (j-i)): increased demand
						costs[t][i][j][0] = sol.inst.cI*(sol.Xt[t-1]- (sol.inst.Dt[t-1] - (j-i)));	// Cost of Inventory
						costs[t][i][j][1] = sol.inst.cB*(sol.inst.Dt[t-1] + (j-i) - sol.Xt[t-1]);	// Cost of backorders
					}

					else if(t==sol.inst.T){
						costs[t][i][j][0] = sol.inst.cI*(sol.Xt[t-1]- (sol.inst.Dt[t-1] - (j-i))) - sol.inst.bP*(sol.inst.Dt[t-1] - (j-i));
						costs[t][i][j][1] = sol.inst.cB*(sol.inst.Dt[t-1] + (j-i) - sol.Xt[t-1]) - sol.inst.bP*sol.Xt[t-1];
					}
				}
			}

			if(t==sol.inst.T+1){
				costs[t][i][0][0] = 0;
				costs[t][i][0][1] = 0;
			}
			
		}
	}

	return costs;
}

//=========================================== Solution related code

float objective_value(Solution sol, vector<float> Dt){
	float obj = 0;
	
	for(int i = 0; i<sol.inst.T; i++){
		obj += max(sol.inst.cI*(sol.Xt[i]-Dt[i]), 
			sol.inst.cB*(Dt[i]-sol.Xt[i]));
	}

	obj -=  sol.inst.bP*min(Dt[sol.inst.T-1],sol.Xt[sol.inst.T-1]); 
	
	return obj;
}

//=========================================== Knowledge Compilation Benders Decomposition code

// It works in two phases: 
//		- forward pass : it calculates the worst-case scenario by traversing the graph
//		- backpropagation : it retrieves the most interesting sub-graph (containing the worst-case scenario)
Solution KC_benders_Master(Instance inst, vector<vector<vector<vector<int> > > > arcsol){
	Solution sol;

	IloEnv env;
	IloModel model(env);

	// Used to retrieve relevant pi after the solve
	vector<vector<int> > pibool;
	pibool.resize(inst.T+2);

	// Vars
	IloArray<IloNumVarArray> pi(env,inst.T+2);
	for(int t = 0; t<inst.T+2;t++){
		// for t = 0 and t = T+1 only one variable in the array
		pi[t] = IloNumVarArray(env, inst.Gamma+1);
		pibool[t].resize(inst.Gamma+1);
		for(int i = 0; i<inst.Gamma+1; i++){
			char name[80];
			pi[t][i] = IloNumVar(env, -IloInfinity, IloInfinity);
			sprintf(name,"pi_%d_%d",t,i);
			pi[t][i].setName(name);
		}
	}

	IloNumVarArray X(env, inst.T);
	for(int t = 0; t<inst.T;t++){
		X[t] = IloNumVar(env);
		char name[80];
		sprintf(name,"X_%d",t);
		X[t].setName(name);
	}

	// Bounbds for X variables
	for(int t = 1; t<inst.T+1;t++){
		model.add(X[t-1]<=inst.X[t-1]);
	}

	// for t in 1...T-1
	for(int t = 1; t<inst.T;t++){
		for(int i = 0; i<inst.Gamma+1; i++){
			for(int j = i; j<inst.Gamma+1; j++){
				// The arc is in the graph =>
				if(j<=i+inst.deltat[t-1] and (arcsol[t][i][j][0] or arcsol[t][i][j][1])){
					IloExpr expr(env);
					
					if(t==1) {
						expr = pi[0][0];
						pibool[0][0] = 1;
					}

					else{ expr = pi[t-1][i];}
					model.add(pi[t][j] - expr >= inst.cI*(X[t-1]- (inst.Dt[t-1] - (j-i))));
					model.add(pi[t][j] - expr >= inst.cB*(inst.Dt[t-1] + (j-i) - X[t-1]));
					pibool[t][j] = 1;
				}
			}
		}
	}

	// for t = T (T-1 -> T)
	int t = inst.T;
	for(int i = 0; i<inst.Gamma+1; i++){
		for(int j = i; j<inst.Gamma+1; j++){
			// cout<<t<<" "<<i<<" "<<j<<endl;
			// cout<<pisol[t][j]
			// cout<<t-1<<" "<<i<<" -> "<<t<<" "<<j<<" "<<arcsol[t][i][j][0]<<" "<<inst.deltat[t]<<endl;
			// cout<<t-1<<" "<<i<<" -> "<<t<<" "<<j<<" "<<arcsol[t][i][j][1]<<" "<<inst.deltat[t]<<endl;
			if(j<=i+inst.deltat[t-1] and (arcsol[t][i][j][0] or arcsol[t][i][j][1])){
				// cout<<t-1<<" "<<i<<" -> "<<t<<" "<<j<<" "<<endl;
				model.add(pi[t][j] - pi[t-1][i] >= inst.cI*(X[t-1]- (inst.Dt[t-1] - (j-i))) - inst.bP*(inst.Dt[t-1] - (j-i)));
				model.add(pi[t][j] - pi[t-1][i] >= inst.cB*(inst.Dt[t-1] + (j-i) - X[t-1]) - inst.bP*X[t-1]);
				pibool[t][j] = 1;
			}
		}
	}

	// Modeling the arcs of the last (T -> T+1) layer
	for(int i = 0; i<inst.Gamma+1; i++){
		if(pibool[t][i]==1){
			model.add(pi[t+1][0] - pi[t][i] >= 0);
		}
	}
	
	model.add(pi[0][0]==0);

	pibool[0][0] = 1;
	pibool[inst.T+1][0] = 1;

	// Obj
	model.add(IloMinimize(env, pi[inst.T+1][0]));

	// IloExpr expr(env);
	// for(int t = 0; t<inst.T+2;t++){
	// 	// for t = 0 and t = T+1 only one variable in the array
	// 	for(int i = 0; i<inst.Gamma+1; i++){
	// 		if(pibool[t][i]){
	// 			expr += pi[t][i];
	// 		}
	// 	}
	// }

	IloCplex cplex(model);
	// cplex.exportModel ("lpex1.lp");
	// cplex.setParam(IloCplex::Param::MIP::Display, 1); //<- displays a bit of info
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
	cplex.setOut(env.getNullStream());

    if ( !cplex.solve() ) {
    	env.error() << "Failed to optimize LP." << endl;
    	throw(-1);
	}

	sol.obj_val = cplex.getValue(pi[inst.T+1][0]);
	sol.Xt.resize(inst.T);
	
	for(int t = 0; t<inst.T; t++){
		sol.Xt[t] = cplex.getValue(X[t]);
	}

	sol.inst = inst;

	env.end();
	
	return sol;
}

// It works in two phases: 
//		- forward pass : it calculates the worst-case scenario by traversing the graph
//		- backpropagation : it retrieves the most interesting sub-graph (containing the worst-case scenario)
vector<vector<vector<vector<int> > > > KC_benders_Subproblem(Solution sol, float approx_coeff, bool use_export){
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios (if a specific decision by the opponent is part of the subgraph)
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j] (It stores the "maximum cumulative cost" to reach period t having consumed j budget units)
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t<sol.inst.T+2;t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j<sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path
	float tmp;
	pi_value[0][0] = 0;	// Start at period 0 cost 0
	for(int t=1; t<sol.inst.T+1;t++){
		for(int j = 0; j<sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j]+costs[t][j][j][0]; 	// Init of pi_value[t][j]
														// It's the value of the dual problem that will store the value of the longest path from the start to t, having consumed j units of budget
			for(int i = 0; i<=j; i++){
				if(j<=i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i]+costs[t][i][j][0] > tmp){
						// if(j== 0){
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][0]<<endl;
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][1]<<endl;
						// }
						tmp = pi_value[t-1][i]+costs[t][i][j][0];	// Maximum cost to reach i + the cost of arc(i, j) in overstock/ stockout and keeps the worst of the two
					}

					if(pi_value[t-1][i]+costs[t][i][j][1] > tmp){
						// if(j== 0){
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][0]<<endl;
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][1]<<endl;
						// }
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					} 
				}
			}
			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		if(pi_value[sol.inst.T][i]>tmp){
			tmp = pi_value[sol.inst.T][i];
		}
	}
	
	pi_value[sol.inst.T+1][0] = tmp;	// Longuest path
	
	//========================== Now the backtrack

	// The variable approx 
	float sub_OPT;
	if(pi_value[sol.inst.T+1][0]>=0){
		sub_OPT = approx_coeff*pi_value[sol.inst.T+1][0];
	}

	else{
		sub_OPT = (1-approx_coeff)*pi_value[sol.inst.T+1][0]+pi_value[sol.inst.T+1][0];
	}
	
	// t = T+1
	pi_subopt_bool[sol.inst.T+1][0] = true;
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		if(pi_value[sol.inst.T][i]>=sub_OPT){
			arcbool[sol.inst.T+1][i][0][0] = 1;
			arcbool[sol.inst.T+1][i][0][1] = 1;
			pi_subopt_bool[sol.inst.T][i] = true;
		} 
	}

	for(int t=sol.inst.T; t>0; t--){
		for(int j = 0; j<sol.inst.Gamma+1; j++){
			for(int i = 0; i<=j; i++){
				// cout<<t<<" "<<j<<" -> "<<t-1<<" "<<i<<endl;
				// cout<<" "<<BoolToString(pi_subopt_bool[t][j])<<" "<<pi_value[t][j]<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][0]<<endl;
				// cout<<" "<<BoolToString(pi_subopt_bool[t][j])<<" "<<pi_value[t][j]<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][1]<<endl;
				if(pi_subopt_bool[t][j] and j<=i+sol.inst.deltat[t-1] and (t!=1 or i==0)){ //last and is specific for first layer of the graph
					// if ==, we go through (t, j) in the matrix and we start again form the node (t-1, i)
					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][0]){
						arcbool[t][i][j][0] = 1;
						pi_subopt_bool[t-1][i] = true;
					}
					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][1]){
						arcbool[t][i][j][1] = 1;
						pi_subopt_bool[t-1][i] = true;
					}
				}
			}
		}
	}

	if(use_export){
		export_budget_graph_json(sol, pi_value, costs, arcbool);
	}
	
	return arcbool;
}

// =================================================================== IN PROGRESS =================================================================================================================
// ====================================================================================================================================================================================

// It works in two phases: 
//		- forward pass : it calculates the worst-case scenario by traversing the graph
//		- backpropagation : it retrieves the most interesting sub-graph (containing the worst-case scenario)
// In this alternative, we had initialized the orthogonality heuristic with greedy & Brays-Curtis
vector<vector<vector<vector<int> > > > KC_benders_Subproblem_HOG(Solution sol, float approx_coeff, int nb_path_to_select, bool use_export, int limit_number_paths){
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios (if a specific decision by the opponent is part of the subgraph)
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j] (It stores the "maximum cumulative cost" to reach period t having consumed j budget units)
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t<sol.inst.T+2;t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j<sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path
	float tmp;
	pi_value[0][0] = 0;	// Start at period 0 cost 0
	for(int t=1; t<sol.inst.T+1;t++){
		for(int j = 0; j<sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j]+costs[t][j][j][0]; 	// Init of pi_value[t][j]
														// It's the value of the dual problem that will store the value of the longest path from the start to t, having consumed j units of budget
			for(int i = 0; i<=j; i++){
				if(j<=i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i]+costs[t][i][j][0] > tmp){
						tmp = pi_value[t-1][i]+costs[t][i][j][0];	// Maximum cost to reach i + the cost of arc(i, j) in overstock/ stockout and keeps the worst of the two
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
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		if(pi_value[sol.inst.T][i]>tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
	}

	pi_value[sol.inst.T+1][0] = tmp;	// Longuest path
	
	//========================== Now the backtrack with HOG

	// Tolerance threshold in relation to the worst possible cost
	float sub_OPT;
	if(pi_value[sol.inst.T+1][0]>=0){
		sub_OPT = approx_coeff*pi_value[sol.inst.T+1][0];
	}

	else{
		sub_OPT = (1-approx_coeff)*pi_value[sol.inst.T+1][0]+pi_value[sol.inst.T+1][0];
	}
	
	// ================================================== IN PROGRESS ==================================================

	vector<Path> candidates_path;	// Store the worst-case scenarios
	Path current_path_buffer;		// Use for recursion

	// After the recursion, candidates_path contains all the worst-case paths that exceed the sub_OPT budget
	//for(int i = 0; i < sol.inst.Gamma+1; i++){
	//	if(pi_value[sol.inst.T][i] >= sub_OPT){
	//		extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, candidates_path, limit_number_paths);
	//	}
	//}

	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(abs(pi_value[sol.inst.T][i] - pi_value[sol.inst.T+1][0]) < 1e-4){	// pi_value[sol.inst.T+1][0] is the worst path, so we will prioritize taking paths with minimal difference in price.
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, candidates_path, limit_number_paths);
		}
	}

	for(int i = 0; i < sol.inst.Gamma+1; i++){
		if(pi_value[sol.inst.T][i] >= sub_OPT && abs(pi_value[sol.inst.T][i] - pi_value[sol.inst.T+1][0]) >= 1e-4){	
			extract_paths_dfs(sol.inst.T, i, pi_value, costs, sol, current_path_buffer, candidates_path, limit_number_paths);
		}
	}

	vector<Path> selected_paths;

	if(!candidates_path.empty()){
		// We take the first scénario we have for reference
		selected_paths.push_back(candidates_path[0]);
		candidates_path.erase(candidates_path.begin());	
		
		// MaxMin loop
		while(selected_paths.size() < nb_path_to_select && !candidates_path.empty()){
			float best_max_min_distance = -1.0;
			int best_candidate_index = -1;

			for(size_t c = 0; c < candidates_path.size(); c++){
				float min_distance_selected = 2.0;	// Brays-Curtis distance is in [0, 1]
				
				for(size_t s = 0; s < selected_paths.size(); s++){
					//float dist = calculate_jaccard_distance(candidates_path[c], selected_paths[s]);
					float dist = calculate_BC_distance(candidates_path[c], selected_paths[s]);
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

		// Reconnected the end of the path to node T+1 like in the previuos code
		int final_budget = selected_paths[p].back().j;
		arcbool[sol.inst.T+1][final_budget][0][0] = 1;
		arcbool[sol.inst.T+1][final_budget][0][1] = 1;
	}

	if(use_export){
		export_budget_graph_json(sol, pi_value, costs, arcbool);
	}
	
	// =================================================================================================================

	return arcbool;
}

// ====================================================================================================================================================================================
// ====================================================================================================================================================================================

// Initialize budget graph  with nominal scenario
vector<vector<vector<vector<int> > > > init_graph(Instance inst){
	vector<vector<vector<vector<int> > > > arcbool;
	arcbool.resize(inst.T+2);

	for(int t = 0; t<inst.T+2;t++){
		// for t = 0 and t = T+1 only one variable in the array
		arcbool[t].resize(inst.Gamma+1);
		for(int i = 0; i<inst.Gamma+1; i++){
			arcbool[t][i].resize(inst.Gamma+1);
			for(int j = 0; j<inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
			if(t>=1) arcbool[t][0][0][0] = 1;
			if(t>=1) arcbool[t][0][0][1] = 1;
		}
	}
	return arcbool;
}

// Initialize budgetgraph with ALL scenarios
vector<vector<vector<vector<int> > > > init_graph_full(Instance inst){
	vector<vector<vector<vector<int> > > > arcbool;
	arcbool.resize(inst.T+2);

	for(int t = 0; t<inst.T+2;t++){
		// for t = 0 and t = T+1 only one variable in the array
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

// Return the union of budget graphs (technically is larger)
vector<vector<vector<vector<int> > > > merge_budget_graph(vector<vector<vector<vector<int> > > > arcbool1, vector<vector<vector<vector<int> > > > arcbool2){
	vector<vector<vector<vector<int> > > > merge_arcbool;

	merge_arcbool.resize(arcbool1.size());

	for(int t = 0; t<arcbool1.size();t++){
		merge_arcbool[t].resize(arcbool1[t].size());
		for(int i = 0; i<arcbool1[t].size(); i++){
			merge_arcbool[t][i].resize(arcbool1[t][i].size());
			for(int j = 0; j<arcbool1[t][i].size(); j++){
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


Benders_Result KC_benders_Main(Instance inst, float approx_coeff, bool use_HOG, bool use_export, int limit_number_paths, int number_orthogonal_axes){
	Solution sol;
	Solution new_sol;
	Solution_ADV sol_adv;
	auto start = high_resolution_clock::now();
	int iter = 0;
	bool stopCriterion = false;
	float proc_time;
	vector<vector<vector<vector<int> > > > arcsol = init_graph(inst);
	vector<vector<vector<vector<int> > > > arcsol_new;

	sol = KC_benders_Master(inst, arcsol);	// Proposes a first (unsatisfactory) solution compared to the nominal scenario
	
	while(!stopCriterion){	// While the solution is not satisfactory we do the merge loop		
		if(use_HOG){
			arcsol_new = KC_benders_Subproblem_HOG(sol, approx_coeff, number_orthogonal_axes, use_export, limit_number_paths);	// Proposes a new worst solution according to the master solution
		} else{
			arcsol_new = KC_benders_Subproblem(sol, approx_coeff, use_export);	// Proposes a new worst solution according to the master solution
		}

		arcsol = merge_budget_graph(arcsol, arcsol_new);		// Merge the worst solution with the current solution
		new_sol = KC_benders_Master(inst, arcsol);				// Proposes a new solution according to the merge
		iter++;
			
		// If the cost increased, the opponent has found a computer breach and we continue in the loop
		if(abs(new_sol.obj_val - sol.obj_val) < 1e-5){	
			stopCriterion = true;
			break;
		}

		sol = new_sol;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	float accuracy = (1./100000);
	proc_time = accuracy*float(duration.count());

	//return make_pair(iter, proc_time);
	return {iter, proc_time, sol.obj_val};
}


//=========================================== Standard Benders Decomposition code
Solution benders_Master(Instance inst, vector<vector<float> > scenarios){
	Solution sol;
	sol.inst = inst;

	// cpo model creation and solving

	IloEnv env;
	IloModel model(env);

	// Vars
	IloNumVar z(env, -IloInfinity, IloInfinity);
	// IloNumVar z(env, -1000, 70);
	z.setName("z");
	IloNumVarArray  X(env, inst.T);
	IloArray<IloNumVarArray> s(env,scenarios.size());
	IloArray<IloNumVarArray> B(env,scenarios.size());
	IloArray<IloNumVarArray> I(env,scenarios.size());
	for(int o = 0; o<scenarios.size();o++){
		s[o] = IloNumVarArray(env, inst.T);
		B[o] = IloNumVarArray(env, inst.T);
		I[o] = IloNumVarArray(env, inst.T);
		for(int t = 0; t<inst.T; t++){
			char name[80];
			s[o][t] = IloNumVar(env, -IloInfinity, IloInfinity);	// Sales
			sprintf(name,"s_%d_%d",o,t);
			s[o][t].setName(name);

			B[o][t] = IloNumVar(env);	// Backorders
			sprintf(name,"B_%d_%d",o,t);
			B[o][t].setName(name);
	
			I[o][t] = IloNumVar(env);	// Inventory
			sprintf(name,"I_%d_%d",o,t);
			I[o][t].setName(name);

			// X[t] = IloNumVar(env, 0, IloInfinity,  IloNumVar::Int);
			X[t] = IloNumVar(env, 0, IloInfinity,  IloNumVar::Float);
			sprintf(name,"X_%d",t);
			X[t].setName(name);

		}
	}

	// Consts
	for(int t = 1; t<inst.T+1;t++){
		model.add(X[t-1]<=inst.X[t-1]);
	}

	for(int t = 0; t<inst.T; t++){
		for(int o = 0; o<scenarios.size();o++){
			model.add(B[o][t] - I[o][t] == scenarios[o][t] - X[t]);	// (2): Flow balance
			IloExpr expr(env);
			for(int i = 0; i<=t; i++){
				expr += s[o][i];
			}
			model.add(expr == scenarios[o][t]-B[o][t]);				// (3)
		}
	}

	// model.add(X[inst.T-1]==14);
	for(int o = 0; o<scenarios.size();o++){
		IloExpr expr(env);
		for(int t = 0; t<inst.T; t++){
			expr += (inst.cI*I[o][t] + inst.cB*B[o][t] - inst.bP*s[o][t]);	// Total cost incurred by the previous scenario
		}

		model.add(z >= expr);
	}

	// Obj
	model.add(IloMinimize(env, z));

	// Solve
	IloCplex cplex(model);
	
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
	cplex.setOut(env.getNullStream());
    if ( !cplex.solve() ) {
    	env.error() << "Failed to optimize LP." << endl;
    	throw(-1);
	}

	vector<float> Xt;
	Xt.resize(inst.T);
	for(int t = 0; t<inst.T; t++){
		Xt[t] = cplex.getValue(X[t]);
	}

	sol.Xt = Xt;
	sol.xt = cumulToStandard(sol.Xt);
	sol.obj_val = cplex.getObjValue();

	env.end();
	return sol;
}

Solution benders_Master_integer(Instance inst, vector<vector<float> > scenarios){
	Solution sol;
	sol.inst = inst;

	// Cpo model creation and solving
	IloEnv env;
	IloModel model(env);

	// Setup cost
	int cP = 2;

	// Vars
	IloNumVar z(env, -IloInfinity, IloInfinity);
	// IloNumVar z(env, -1000, 70);
	z.setName("z");
	IloNumVarArray  X(env, inst.T);
	IloNumVarArray y(env, inst.T);
	IloArray<IloNumVarArray> s(env,scenarios.size());
	IloArray<IloNumVarArray> B(env,scenarios.size());
	IloArray<IloNumVarArray> I(env,scenarios.size());
	for(int o = 0; o<scenarios.size();o++){
		s[o] = IloNumVarArray(env, inst.T);
		B[o] = IloNumVarArray(env, inst.T);
		I[o] = IloNumVarArray(env, inst.T);
		for(int t = 0; t<inst.T; t++){
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

			// X[t] = IloNumVar(env, 0, IloInfinity,  IloNumVar::Int);
			X[t] = IloNumVar(env, 0, IloInfinity,  IloNumVar::Float);
			sprintf(name,"X_%d",t);
			X[t].setName(name);

			y[t] = IloNumVar(env, 0, 1,  IloNumVar::Int);
			sprintf(name,"y_%d",t);
			y[t].setName(name);

		}
	}

	// Consts
	for(int t = 1; t<inst.T+1;t++){
		model.add(X[t-1]<inst.X[t-1]);
	}

	for(int t = 0; t<inst.T; t++){
		for(int o = 0; o<scenarios.size();o++){
			model.add(B[o][t] - I[o][t] == scenarios[o][t] - X[t]); // (2)
			IloExpr expr(env);
			for(int i = 0; i<=t; i++){
				expr += s[o][i];
			}
			model.add(expr == scenarios[o][t]-B[o][t]);				// (3)
		}
	}

	// model.add(X[inst.T-1]==14);
	for(int o = 0; o<scenarios.size();o++){
		IloExpr expr(env);
		for(int t = 0; t<inst.T; t++){
			expr += (inst.cI*I[o][t] + inst.cB*B[o][t] - inst.bP*s[o][t]+cP*y[t]);
		}
		model.add(z>= expr);
	}


	// Constraints for y
	float M = inst.Dt[inst.T-1] + inst.Gamma;
	model.add(X[1]<=y[1]*M);

	for(int t = 1; t<inst.T; t++){
		model.add(X[t]-X[t-1]<=y[t]*M);
	}

	// Obj
	model.add(IloMinimize(env, z));

	// Solve
	IloCplex cplex(model);
	
	// cplex.setParam(IloCplex::Param::MIP::Display, 1);
	// cplex.exportModel ("ben_main.lp");
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
	cplex.setOut(env.getNullStream());
    if ( !cplex.solve() ) {
    	env.error() << "Failed to optimize LP." << endl;
    	throw(-1);
	}

	// cout<<"solved"<<endl;
	vector<float> Xt;
	Xt.resize(inst.T);
	for(int t = 0; t<inst.T; t++){
		Xt[t] = cplex.getValue(X[t]);
	}

	sol.Xt = Xt;
	sol.xt = cumulToStandard(sol.Xt);
	sol.obj_val = cplex.getObjValue();

	env.end();
	return sol;
}


// Solve subproblem using CPLEX: NOT USED, probably not working
Solution_ADV benders_Subproblem(Solution sol){
	Solution_ADV sol_adv;

	IloEnv env;
	IloModel model(env);

	// Used to retrieve relevant pi after the solve
	vector<vector<bool> > pibool;
	pibool.resize(sol.inst.T+2);

	// Vars
	IloArray<IloNumVarArray> pi(env,sol.inst.T+2);
	for(int t = 0; t<sol.inst.T+2;t++){
		// for t = 0 and t = T+1 only one variable in the array
		pi[t] = IloNumVarArray(env, sol.inst.Gamma+1);
		pibool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			char name[80];
			pi[t][i] = IloNumVar(env, -IloInfinity, IloInfinity);
			sprintf(name,"pi_%d_%d",t,i);
			pi[t][i].setName(name);
		}
	}

	// Consts
	// for t in 1...T-1
	for(int t = 1; t<sol.inst.T;t++){
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			for(int j = i; j<sol.inst.Gamma+1; j++){
				// The arc is in the graph =>
				if(j<=i+sol.inst.deltat[t]){
					IloExpr expr(env);
					if(t==1) {
						expr = pi[0][0];
						pibool[0][0] = true;
					}

					else{ expr = pi[t-1][i];}
					model.add(pi[t][j] - expr >= sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (j-i))));
					model.add(pi[t][j] - expr >= sol.inst.cB*(sol.inst.Dt[t] + (j-i) - sol.Xt[t]));
					pibool[t][j] = true;
				}
			}
		}
	}

	// for t = T
	int t = sol.inst.T;
	for(int i = 0; i<sol.inst.Gamma+1; i++){
		for(int j = i; j<sol.inst.Gamma+1; j++){
			if(j<=i+sol.inst.deltat[t]){
				model.add(pi[t][j] - pi[t-1][i] >= sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (j-i))) - sol.inst.bP*(sol.inst.Dt[t] - (j-i)));
				model.add(pi[t][j] - pi[t-1][i] >= sol.inst.cB*(sol.inst.Dt[t] + (j-i) - sol.Xt[t]) - sol.inst.bP*sol.Xt[t]);
				pibool[t][j] = true;
			}
		}
	}


	for(int i = 0; i<sol.inst.Gamma+1; i++){
		if(pibool[t][i]) model.add(pi[sol.inst.T+1][0] - pi[sol.inst.T][i] >= 0);
	}

	// model.add(pi[sol.inst.T][0]-pi[0][0]>= 0);
	model.add(pi[0][0]==0);
	pibool[0][0] = true;
	pibool[sol.inst.T+1][0] = true;

	// Obj
	// model.add(IloMinimize(env, pi[sol.inst.T][0]));
	IloExpr expr(env);
	for(int t = 0; t<sol.inst.T+2;t++){
		// for t = 0 and t = T+1 only one variable in the array
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			if(pibool[t][i]){
				expr += pi[t][i];
			}
		}
	}

	model.add(IloMinimize(env, expr)); // Bon objectif pour etre sur que les potentiels collent leur borne

	IloCplex cplex(model);
	// cplex.exportModel ("lpex1.lp");
	cplex.setParam(IloCplex::Param::MIP::Display, 1);
	cplex.setOut(env.getNullStream());
    if(!cplex.solve()) {
    	env.error() << "Failed to optimize LP." << endl;
    	throw(-1);
	}

	// Retrieve value for pi var
	vector<vector<float> > pisol;
	pisol.resize(sol.inst.T+2);
	for(int t = 0; t<sol.inst.T+2;t++){
		// for t = 0 and t = T+1 only one variable in the array
		pisol[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			if(pibool[t][i]){
				pisol[t][i] = cplex.getValue(pi[t][i]);
				// cout<<t<<" "<<i<<" "<<pisol[t][i]<<endl;
			}
		}
	}

	cout<<"======================LONGEST PATH : "<<pisol[sol.inst.T+1][0]<<endl;
	float previous_val = pisol[sol.inst.T+1][0];
	//================= DANGEROUS, WORKS ONLY IF GAMMA IS FULLY USED, SHOULD CHANGE THIS ====================
	int previous_budget = sol.inst.Gamma;

	vector<int> offset;
	offset.resize(sol.inst.T);

	for(int t = sol.inst.T; t>=1; t--){
		for(int i=previous_budget; i>=0; i--){
			if(t == sol.inst.T){
				if( pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (previous_budget-i))) - sol.inst.bP*(sol.inst.Dt[t] - (previous_budget-i))) { // worse for deltat negative
					offset[t-1] = -(previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}

				else if(pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cB*(sol.inst.Dt[t] + (previous_budget-i) - sol.Xt[t]) - sol.inst.bP*sol.Xt[t]) { // worse for deltat negative
					offset[t-1] = (previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
				}
			}

			else if (t == 1) {
				if( pisol[t][previous_budget] - pisol[0][0] == sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (previous_budget-i)))){ // worse for deltat negative
					offset[t-1] = -(previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}
			
				else if(pisol[t][previous_budget] - pisol[0][0] == sol.inst.cB*(sol.inst.Dt[t] + (previous_budget-i) - sol.Xt[t])) { // worse for deltat negative
					offset[t-1] = (previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}
			}
			
			// Generic case
			else{
				if( pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (previous_budget-i)))){ // worse for deltat negative
					offset[t-1] = -(previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
				}
			
				else if(pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cB*(sol.inst.Dt[t] + (previous_budget-i) - sol.Xt[t])) { // worse for deltat negative
					offset[t-1] = (previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}
			
			}
		}
	}

	vector<float> Dt;
	Dt.resize(sol.inst.T);
	for(int t = 0; t<sol.inst.T; t++){
		Dt[t] = sol.inst.Dt[t] + offset[t];
	}

	env.end();

	sol_adv.Dt = Dt;
	return sol_adv;
}

// Solve subproblem with dynamic prog
Solution_ADV benders_Subproblem_DP(Solution sol){
	Solution_ADV sol_adv;
	vector<vector<vector<vector<int> > > > arcbool; // Bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; 				// Value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); // Costs of all arcs of the budget graph

	pi_value.resize(sol.inst.T+2);
	pi_subopt_bool.resize(sol.inst.T+2);
	arcbool.resize(sol.inst.T+2);
	for(int t = 0; t<sol.inst.T+2;t++){
		pi_value[t].resize(sol.inst.Gamma+1);
		pi_subopt_bool[t].resize(sol.inst.Gamma+1);
		arcbool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			arcbool[t][i].resize(sol.inst.Gamma+1);
			for(int j = 0; j<sol.inst.Gamma+1; j++){
				arcbool[t][i][j].resize(2);
			}
		}
	}

	// Dynamic prog. for longest path
	float tmp;
	pi_value[0][0] = 0;

	for(int t=1; t<sol.inst.T+1;t++){
		for(int j = 0; j<sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j]+costs[t][j][j][0];	// Init of pi_value[t][j]
			for(int i = 0; i<=j; i++){
				if(j<=i+sol.inst.deltat[t-1]){
					if(pi_value[t-1][i]+costs[t][i][j][0] > tmp){
						// if(j== 0){
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][0]<<endl;
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][1]<<endl;
						// }
						tmp = pi_value[t-1][i]+costs[t][i][j][0];
					} 
	
					if(pi_value[t-1][i]+costs[t][i][j][1] > tmp){
						// if(j== 0){
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][0]<<endl;
						// 	cout<<"=============="<<t<<" "<<pi_value[t-1][i]<<" "<<costs[t][i][j][1]<<endl;
						// }
						tmp = pi_value[t-1][i]+costs[t][i][j][1];
					} 
				}
			}
			pi_value[t][j] = tmp;
		}
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		if(pi_value[sol.inst.T][i]>tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
	}

	pi_value[sol.inst.T+1][0] = tmp;
	
	//========================== Now the backtrack

	vector<float> scenario;
	scenario.resize(sol.inst.T);

	// t = T+1
	pi_subopt_bool[sol.inst.T+1][0] = true;
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		if(pi_value[sol.inst.T][i]==pi_value[sol.inst.T+1][0]){
			arcbool[sol.inst.T+1][i][0][0] = 1;
			arcbool[sol.inst.T+1][i][0][1] = 1;
			pi_subopt_bool[sol.inst.T][i] = true;
		} 
	}

	for(int t=sol.inst.T; t>0; t--){
		for(int j = 0; j<sol.inst.Gamma+1; j++){
			for(int i = 0; i<=j; i++){
				if(pi_subopt_bool[t][j] and j<=i+sol.inst.deltat[t-1] and (t!=1 or i==0)){ // Last and is specific for first layer of the graph
					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][0]){
						arcbool[t][i][j][0] = 1;
						pi_subopt_bool[t-1][i] = true;
						scenario[t-1] = sol.inst.Dt[t-1] - (j-i);
						break;
					}

					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][1]){
						arcbool[t][i][j][1] = 1;
						pi_subopt_bool[t-1][i] = true;
						scenario[t-1] = sol.inst.Dt[t-1] + (j-i);
						break;
					}
				}
			}
		}
	}

	// cout<<"================= BEGIN TEST =================="<<endl;

	// Solution test = KC_benders_Master(sol.inst, arcbool);

	// cout<<"new sol (with Graph LP): ";
	// display_vector_float(test.Xt);
	// cout<<"new sol value : "<<test.obj_val<<endl;

	// cout<<"================= END TEST =================="<<endl;


	// // Display the subgraph
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
	// // Display the subgraph
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


Benders_Result benders_Main(Instance inst){
	auto start = high_resolution_clock::now();

	Solution sol;
	Solution new_sol;
	Solution_ADV sol_adv;
	bool stopCriterion = false;
	float proc_time;
	vector<vector<float> > scenarios;
	scenarios.resize(0);

	// Only nominal scenario
	// cout<<"nominal scenario : ";
	// vector<float> debug = {0,0,0,0,1,1,2,4,5,5,6,9,9,11,11};
	// scenarios.push_back(debug);
	// display_vector_float(inst.Dt);
	scenarios.push_back(inst.Dt);

	sol = benders_Master(inst, scenarios);
	
	int i = 2;
	while(!stopCriterion){
		// cout<<"============= ITERATION"<<i<<endl;
		sol_adv = benders_Subproblem_DP(sol);
		scenarios.push_back(sol_adv.Dt);
		new_sol = benders_Master(inst, scenarios);
		
		i++;
		
		if(objective_value(new_sol, sol_adv.Dt) == objective_value(sol, sol_adv.Dt)){
			stopCriterion = true;
			break;
		}
		sol = new_sol;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	float accuracy = (1./100000);
	proc_time = accuracy*float(duration.count());
	
	//return make_pair(i, proc_time);
	return {i, proc_time, sol.obj_val};
}

// =========================================== Main

vector<string> list_dir(const char *path) {
vector<string> allfile;
   struct dirent *entry;
   DIR *dir = opendir(path);
   
   if (dir == NULL) {
      return allfile;
   }

   while ((entry = readdir(dir)) != NULL) {
   allfile.push_back(entry->d_name);
   }

   closedir(dir);
   return allfile;
}

int main(int argc, const char* argv[]){
	Benders_Result benders_sol;
	Benders_Result benders_sol_augmented;
	Benders_Result benders_sol_augmented_HOG;
	float approx_coeff;
	float time, timeKC, timeKCHOG;
	int iter, iterKC, iterKCHOG;
	// Simulation parameters
	int limit_number_paths;		// Used for the DFS algo
	limit_number_paths = 2000;	// 2000
	int number_orthogonal_axes;	// Number of orthogonal axes we want for the heuristics
	number_orthogonal_axes = 3;	// 5
	bool use_export = false;	// To be corrected before use

	//====================================================================== IN PROGRESS ======================================================================

	// Création of the folder architecture
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	ostringstream oss_exp;
    oss_exp << "/result_benders_" << put_time(&tm, "%Y-%m-%d_%H%M") << "_n=" << number_orthogonal_axes << "_l=" << limit_number_paths;
    std::string experience_name = oss_exp.str();

    ostringstream oss_folder;
    oss_folder << "results" << experience_name;
    std::string folder_path = oss_folder.str();

    if(fs::create_directories(folder_path)){
        std::cout << "Le dossier '" << folder_path << "' a été créé avec succès." << std::endl;
    } else{
        std::cout << "Le dossier '" << folder_path << "' existe déjà." << std::endl;
    }

	ostringstream oss_classic; 
	ostringstream oss_augmented; 
	ostringstream oss_augmented_HOG;
	oss_classic 	  << folder_path << experience_name << "_classic.csv";
	oss_augmented 	  << folder_path << experience_name << "_HOG.csv";
	oss_augmented_HOG << folder_path << experience_name << "_augmented_HOG.csv";

	ofstream output_classic(oss_classic.str());
	ofstream output_augmented(oss_augmented.str());
	ofstream output_augmented_HOG(oss_augmented_HOG.str());

	cout << "Enregistrement des résultats dans : " << folder_path << endl;

	//================================================================= TEMPORAIRE ========================================================================================
	vector<string> file_list;
	int choice_instances;
	choice_instances = 1;
	
	if(choice_instances == 1){
		file_list = list_dir("/home/mfrancineh/Documents/REPO/STG_1RO_LAASCNRS/STG/PROJET/bae/parsed_large_instances/");
	} else if(choice_instances == 2){
		file_list = list_dir("/home/mfrancineh/Documents/REPO/STG_1RO_LAASCNRS/STG/PROJET/bae/other_instances/");
	} else{
		file_list = list_dir("/home/mfrancineh/Documents/REPO/STG_1RO_LAASCNRS/STG/PROJET/bae/test/");
	}
  	
	//=========================================================================================================================================================

  	int total_files = file_list.size();

	if (total_files <= 2) {
    	cerr << "ERREUR : Aucun fichier d'instance trouvé. Vérifiez le chemin du dossier." << endl;
    	return -1;
	}

	int nbInst = total_files-2;
	cout << "Succès : " << nbInst << " fichiers trouvés dans le dossier." << endl;

	//=========================================================================================================================================================

  	Instance inst;
  	vector<int> debug;
  	vector<int> debug2;
	string filename;

  	int seed = 31415;
  	srand (seed);

	for(int Gamma=1; Gamma<100; Gamma+=20){
		for(int tau=0; tau<11; tau+=2){		
			iter      = 0;
			iterKC    = 0;
			iterKCHOG = 0;
			time      = 0;
			timeKC    = 0;
			timeKCHOG = 0;
			debug.resize(0);
			debug2.resize(0);
			for(int i = 2; i<total_files; i++ ){
				
				//=========================================================================================================================================================
					
				if(choice_instances == 1){
					filename = "parsed_large_instances/" + file_list[i];
				} else if(choice_instances == 2){
					filename = "other_instances/" + file_list[i];
				} else{
					filename = "test/" + file_list[i];
				}
				
				//=========================================================================================================================================================

				cout << "\n" << filename << " " << Gamma << " " << tau << " " << endl;
				
				inst = read_instance_randomized(filename, Gamma);

				// Classique
				benders_sol = benders_Main(inst);
				iter += benders_sol.iter;
				time += benders_sol.time;
				cout << "STANDARD-done (Obj :" << benders_sol.obj_value << ")" << endl;

				approx_coeff = float(tau)/10;

				// KC
				benders_sol_augmented = KC_benders_Main(inst, approx_coeff, false, use_export, limit_number_paths, number_orthogonal_axes);	// First bool is to use KC with HOG
				iterKC += benders_sol_augmented.iter;
				timeKC += benders_sol_augmented.time;
				cout << "KC-------done (Obj :" << benders_sol_augmented.obj_value << ")" <<endl;

				// KC with HOG
				benders_sol_augmented_HOG = KC_benders_Main(inst, approx_coeff, true, use_export, limit_number_paths, number_orthogonal_axes);
				iterKCHOG += benders_sol_augmented_HOG.iter;
				timeKCHOG += benders_sol_augmented_HOG.time;
				cout << "KC_HOG---done (Obj :" << benders_sol_augmented_HOG.obj_value << ")"<< endl;

				// Quality control of the solution
				float eps = 1e-4;
				if(abs(benders_sol.obj_value - benders_sol_augmented.obj_value) > eps || abs(benders_sol.obj_value - benders_sol_augmented_HOG.obj_value) > eps){
					cout << "ALERTE DEGRADATION" << endl;
				} else{
					cout << "Qualité valide" << endl;
				}
			}

			// Standard with KC standard 
			output_classic << "STANDARD," << Gamma << "," << tau << ","
                           << float(iter)/nbInst << "," << float(time)/nbInst << endl;
            output_classic << "KC," << Gamma << "," << tau << ","
                           << float(iterKC)/nbInst << "," << float(timeKC)/nbInst << endl;

			// Standard with KC augmented (HOG)
            output_augmented << "STANDARD," << Gamma << "," << tau << ","
                             << float(iter)/nbInst << "," << float(time)/nbInst << endl;
            output_augmented << "KC_HOG," << Gamma << "," << tau << ","
                             << float(iterKCHOG)/nbInst << "," << float(timeKCHOG)/nbInst << endl;

			// KC standard with KC augmented
			output_augmented_HOG << "KC," << Gamma << "," << tau << ","
                             << float(iterKC)/nbInst << "," << float(timeKC)/nbInst << endl;
            output_augmented_HOG << "KC_HOG," << Gamma << "," << tau << ","
                             << float(iterKCHOG)/nbInst << "," << float(timeKCHOG)/nbInst << endl;
		}
	}

	output_classic.close();
	output_augmented.close();
	output_augmented_HOG.close();
	cout<<"========== END OF THE PROGRAM =========="<<endl;
}