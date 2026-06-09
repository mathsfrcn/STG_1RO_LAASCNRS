#include <ilcplex/ilocplex.h>
#include<vector>
#include<string>
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

using namespace std;
using namespace std::chrono;


//=========================================== structures

struct Instance
{
	int T, cI, cB, bP;
	float Gamma;
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

//===========================================misc.

//displan an float vector, for debugging
void display_vector_float(vector<float> v){
	for (int i = 0; i < v.size(); ++i)
	{
		cout<<v[i]<<" ";
	}
	cout<<endl;
}

//displan an int vector, for debugging
void display_vector_int(vector<int> v){
	for (int i = 0; i < v.size(); ++i)
	{
		cout<<v[i]<<" ";
	}
	cout<<endl;
}

string BoolToString(bool b)
{
  return b ? "1" : "0";
}

//===========================================instance related code

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

//read an instance from psplib problem 58 pspInstance. Products are agregated.
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

	inst.cI = 3; //stock cost
	inst.cB = 6; //backorder cost
	inst.bP = 10; //selling price
	//inst.Gamma = int(inst.Dt[inst.Dt.size()-1]);
	inst.Gamma = budget;
	inst.deltat.resize(inst.T);
	for(int t = 0; t<inst.T;t++){
		inst.deltat[t] = int(inst.Dt[t]/float(2));
		// cout<<inst.deltat[t]<<endl;
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
			//arbitrary to have nice instances
			inst.dt[j] += tmp + int(rand() % 2);
		}
	}
	inst.Dt = standardToCumul(inst.dt);

	
	inst.cB = rand() % 10 + 10; //backorder cost
	inst.cI = rand() % 10 + 10; //stock cost
	inst.bP = rand() % 10 + 10; //selling price
	//inst.Gamma = int(inst.Dt[inst.Dt.size()-1]);
	inst.Gamma = budget;
	inst.Dt = standardToCumul(inst.dt);
	inst.deltat.resize(inst.T);

	//display_vector_float(inst.Dt);
	for(int t = 0; t<inst.T;t++){
		//recours temporaire(pas propre)
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
				//cout<<inst.Dt[t-1] <<" "<< inst.deltat[t-1] <<" " <<inst.Dt[t] <<" "<< inst.deltat[t]<<" "<<t<<endl;
				inst.deltat[t] = rand()%(int(inst.Dt[t]));
				//cout<<"fin"<<endl;
			}
		}
		else{
			inst.deltat[t] = rand()%(int(inst.Dt[t]));
			while(inst.Dt[t-1] + inst.deltat[t-1] > inst.Dt[t] - inst.deltat[t]){
				//cout<<inst.Dt[t-1] <<" "<< inst.deltat[t-1] <<" " <<inst.Dt[t] <<" "<< inst.deltat[t]<<" "<<t<<endl;
				inst.deltat[t] = rand()%(int(inst.Dt[t]));
			}
		}
		// cout<<inst.deltat[t]<<endl;
	}
	//cout<<"poueeeeeeeeeeeeeeeet000"<<endl;
	inst.X.resize(inst.T);
	//recours temporaire(pas propre)
	for(int t = 0; t<inst.T;t++){
		if(inst.Dt[t]==0){
			inst.X[t] = 0;
		}
		else{
			inst.X[t] = int(rand()%(int(0.4*inst.Dt[t])+1) + 0.8*inst.Dt[t]);
		}
	}

	file.close();
	//cout<<"poueeeeeeeeeeeeeeeet"<<endl;
	return inst;
}

Instance duplicate_instance(Instance inst){
	return inst;
}


//returns the graph with all the costs for KC_subproblem
vector<vector<vector<vector<float> > > > budget_graph_cost(Solution sol){
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
						costs[t][i][j][0] = sol.inst.cI*(sol.Xt[t-1]- (sol.inst.Dt[t-1] - (j-i)));
						costs[t][i][j][1] = sol.inst.cB*(sol.inst.Dt[t-1] + (j-i) - sol.Xt[t-1]);
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

	// cout<<"nominal scenario cost :"<<endl;
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

//===========================================solution related code

float objective_value(Solution sol, vector<float> Dt){
	//cout<<"=============obj"<<endl;
	float obj = 0;
	for(int i = 0; i<sol.inst.T; i++){
		//cout<<sol.inst.cI*(sol.Xt[i]-Dt[i])<<" "<<sol.inst.cB*(Dt[i]-sol.Xt[i]) - sol.inst.bP*min(Dt[i],sol.Xt[i])<<endl;
		obj += max(sol.inst.cI*(sol.Xt[i]-Dt[i]), 
			sol.inst.cB*(Dt[i]-sol.Xt[i]));
	}
	obj -=  sol.inst.bP*min(Dt[sol.inst.T-1],sol.Xt[sol.inst.T-1]); 
	return obj;
}

//=========================================== Knowledge Compilation Benders Decomposition code

Solution KC_benders_Master(Instance inst, vector<vector<vector<vector<int> > > > arcsol){
	Solution sol;

	IloEnv env;
	IloModel model(env);

	//used to retrieve relevant pi after the solve
	vector<vector<int> > pibool;
	pibool.resize(inst.T+2);

	//vars
	IloArray<IloNumVarArray> pi(env,inst.T+2);
	for(int t = 0; t<inst.T+2;t++){
		//for t = 0 and t = T+1 only one variable in the array
		pi[t] = IloNumVarArray(env, inst.Gamma+1);
		pibool[t].resize(inst.Gamma+1);
		for(int i = 0; i<inst.Gamma+1; i++){
			char name[80];
			pi[t][i] = IloNumVar(env, -IloInfinity, IloInfinity);
			// cout<<pi[t][i].getLB()<<endl;
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
	//csts

	//bounbds for X variables
	for(int t = 1; t<inst.T+1;t++){
		model.add(X[t-1]<=inst.X[t-1]);
	}

	//for t in 1...T-1
	for(int t = 1; t<inst.T;t++){
		for(int i = 0; i<inst.Gamma+1; i++){
			for(int j = i; j<inst.Gamma+1; j++){
				//the arc is in the graph =>
				if(j<=i+inst.deltat[t-1] and (arcsol[t][i][j][0] or arcsol[t][i][j][1])){
					// cout<<t-1<<" "<<i<<" -> "<<t<<" "<<j<<endl;
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
	//for t = T (T-1 -> T)
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

	//modeling the arcs of the last (T -> T+1) layer
	for(int i = 0; i<inst.Gamma+1; i++){
		if(pibool[t][i]==1){
			model.add(pi[t+1][0] - pi[t][i] >= 0);
		}
	}
	
	// model.add(X[inst.T-1]==12.8947);
	// model.add(X[inst.T-1]==14);

	model.add(pi[0][0]==0);

	pibool[0][0] = 1;
	pibool[inst.T+1][0] = 1;

	// model.add(X[inst.T-1]==11.9474);

	//obj

	model.add(IloMinimize(env, pi[inst.T+1][0]));

	// IloExpr expr(env);
	// for(int t = 0; t<inst.T+2;t++){
	// 	//for t = 0 and t = T+1 only one variable in the array
	// 	for(int i = 0; i<inst.Gamma+1; i++){
	// 		if(pibool[t][i]){
	// 			expr += pi[t][i];
	// 		}
	// 	}
	// }

	// model.add(IloMinimize(env, expr)); //objectif pour etre sur que les potentiels collent leur borne, et ne laisse plus de marge à X (pas sur de ça)

	IloCplex cplex(model);
	// cplex.exportModel ("lpex1.lp");
	// cplex.setParam(IloCplex::Param::MIP::Display, 1); //<- displays a bit of info
	cplex.setParam(IloCplex::Param::MIP::Display, 0);
	cplex.setOut(env.getNullStream());
    if ( !cplex.solve() ) {
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

	// cout<<inst.T<<endl;
	// cout<<"LONGEST PATH IN GRAPH:"<<cplex.getValue(pi[inst.T+1][0])<<endl;
	sol.obj_val = cplex.getValue(pi[inst.T+1][0]);
	sol.Xt.resize(inst.T);
	for(int t = 0; t<inst.T; t++){
		// cout<<t<<" "<<cplex.getValue(X[t])<<" "<<inst.Dt[t]<<" "<<cplex.getValue(pi[t][0])<<endl;
		sol.Xt[t] = cplex.getValue(X[t]);
		// cout<<sol.Xt[t]<<endl;
	}
	sol.inst = inst;

	// cout<<objective_value(sol, inst.Dt)<<endl;;


	env.end();
	return sol;
}

vector<vector<vector<vector<int> > > > KC_benders_Subproblem(Solution sol, float approx_coeff){
	vector<vector<vector<vector<int> > > > arcbool; //bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; //value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); //costs of all arcs

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

	//dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;
	for(int t=1; t<sol.inst.T+1;t++){
		for(int j = 0; j<sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j]+costs[t][j][j][0];//init of pi_value[t][j]
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
		// cout<<"=============="<<t<<" "<<pi_value[t][0]<<endl;
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		if(pi_value[sol.inst.T][i]>tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
		// cout<<i<<", "<<pi_value[sol.inst.T][i]<<endl;
	}
	pi_value[sol.inst.T+1][0] = tmp;
	// cout<<"longest path : "<<pi_value[sol.inst.T+1][0]<<endl;
	

	//==========================now the backtrack
	//float sub_OPT = -114;
	float sub_OPT;
	if(pi_value[sol.inst.T+1][0]>=0){
		sub_OPT = approx_coeff*pi_value[sol.inst.T+1][0];
	}
	else{
		sub_OPT = (1-approx_coeff)*pi_value[sol.inst.T+1][0]+pi_value[sol.inst.T+1][0];
	}
	

	//t = T+1
	pi_subopt_bool[sol.inst.T+1][0] = true;
	//cout<<"poeut"<<endl;
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		// cout<<i<<" "<<pi_value[sol.inst.T][i]<<" "<<sub_OPT<<endl;
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
					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][0]){
						arcbool[t][i][j][0] = 1;
						pi_subopt_bool[t-1][i] = true;
						// cout<<"test passed : "<<t<<" "<<j<<" -> "<<t-1<<" "<<i<<" "<<"0"<<endl;
					}
					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][1]){
						arcbool[t][i][j][1] = 1;
						pi_subopt_bool[t-1][i] = true;
						// cout<<"test passed : "<<t<<" "<<j<<" -> "<<t-1<<" "<<i<<" "<<"1"<<endl;
					}
				}
			}
		}
	}

	//display the subgraph
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

//initialize budget graph  with nominal scenario
vector<vector<vector<vector<int> > > > init_graph(Instance inst){

	vector<vector<vector<vector<int> > > > arcbool;
	arcbool.resize(inst.T+2);

	for(int t = 0; t<inst.T+2;t++){
		//for t = 0 and t = T+1 only one variable in the array
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

//initialize budgetgraph with ALL scenarios
vector<vector<vector<vector<int> > > > init_graph_full(Instance inst){

	vector<vector<vector<vector<int> > > > arcbool;
	arcbool.resize(inst.T+2);

	for(int t = 0; t<inst.T+2;t++){
		//for t = 0 and t = T+1 only one variable in the array
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

//return the union of budget graphs (technically is larger)
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
					// cout<<t<<" "<<i<<" "<<j<<" "<<arcbool1[t][i][j][0]<<" "<<arcbool2[t][i][j][0]<<endl;
					merge_arcbool[t][i][j][0] = 1;
				}
				if(arcbool1[t][i][j][1] == 1 or arcbool2[t][i][j][1] == 1){
					// cout<<t<<" "<<i<<" "<<j<<" "<<arcbool1[t][i][j][1]<<" "<<arcbool2[t][i][j][1]<<endl;
					merge_arcbool[t][i][j][1] = 1;
				}
			}
		}
	}


	return merge_arcbool;
}


pair<int, float> KC_benders_Main(Instance inst, float approx_coeff){

	auto start = high_resolution_clock::now();

	Solution sol;
	Solution new_sol;
	Solution_ADV sol_adv;
	int iter = 0;
	bool stopCriterion = false;
	float proc_time;
	vector<vector<vector<vector<int> > > > arcsol = init_graph(inst);
	// vector<vector<vector<vector<int> > > > arcsol = init_graph_full(inst);
	vector<vector<vector<vector<int> > > > arcsol_new;
	//only nominal scenario
	// cout<<"nominal scenario : ";
	// display_vector_float(inst.Dt);

	sol = KC_benders_Master(inst, arcsol);

	// cout<<"first sol cost : "<<sol.obj_val<<endl;;
	// display_vector_float(sol.Xt);

	// vector<vector<vector<vector<int> > > > tmp =  KC_benders_Subproblem(sol,approx_coeff);

	while(!stopCriterion){
		arcsol_new = KC_benders_Subproblem(sol, approx_coeff);
		// cout<<"subproblem solved"<<endl;
		arcsol = merge_budget_graph(arcsol, arcsol_new);
		// cout<<"budget graphs merged"<<endl;
		new_sol = KC_benders_Master(inst, arcsol);
		iter++;
		// cout<<"master problem solved"<<endl;
		// cout<<"new sol : ";
		// display_vector_float(new_sol.Xt);
		// cout<<"new sol value : "<<new_sol.obj_val<<endl;

		// cout<<"============ "<< new_sol.obj_val << " " << sol.obj_val<<endl;
		if(new_sol.obj_val == sol.obj_val){
			stopCriterion = true;
			break;
		}

		
		// cout<<"worst case : ";
		// display_vector_float(sol_adv.Dt);
		sol = new_sol;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	float accuracy = (1./100000);
	proc_time = accuracy*float(duration.count());

	return make_pair(iter, proc_time);
}


//===========================================Standard Benders Decomposition code

Solution benders_Master(Instance inst, vector<vector<float> > scenarios){
	Solution sol;
	sol.inst = inst;

	//cpo model creation and solving

	IloEnv env;
	IloModel model(env);


	//vars
	IloNumVar z(env, -IloInfinity, IloInfinity);
	//IloNumVar z(env, -1000, 70);
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

		}
	}

	//cout<<"pouet1"<<endl;

	//consts

	for(int t = 1; t<inst.T+1;t++){
		model.add(X[t-1]<=inst.X[t-1]);
	}

	for(int t = 0; t<inst.T; t++){
		for(int o = 0; o<scenarios.size();o++){
			model.add(B[o][t] - I[o][t] == scenarios[o][t] - X[t]); //(2)
			IloExpr expr(env);
			for(int i = 0; i<=t; i++){
				expr += s[o][i];
			}
			model.add(expr == scenarios[o][t]-B[o][t]);			//(3)
		}
	}

	// model.add(X[inst.T-1]==14);

	//cout<<"pouet2"<<endl;
	for(int o = 0; o<scenarios.size();o++){
		IloExpr expr(env);
		for(int t = 0; t<inst.T; t++){
			expr += (inst.cI*I[o][t] + inst.cB*B[o][t] - inst.bP*s[o][t]);
		}
		model.add(z>= expr);
	}

	//obj
	model.add(IloMinimize(env, z));

	//cout<<"pouet3"<<endl;

	//solve
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

	// for(int o = 0; o<scenarios.size();o++){
	// 	for(int t = 0; t<inst.T; t++){
	// 		cout<<cplex.getValue(B[o][t])<<" "<<cplex.getValue(I[o][t])<<" "<<scenarios[o][t]<<" "<<cplex.getValue(X[t])<<" "<<cplex.getValue(s[o][t])<<endl;
	// 	}
	// }

	// float debug = 0;
	// for(int o = 0; o<scenarios.size();o++){
	// 	for(int t = 0; t<inst.T; t++){
	// 		//debug += inst.cI*cplex.getValue(I[o][t]) + inst.cB*cplex.getValue(B[o][t]) - inst.bP*cplex.getValue(s[o][t]);
	// 		debug += - inst.bP*cplex.getValue(s[o][t]);
	// 	}
	// }
	// cout<<"solution cost (cplex issue ?): "<<debug<<endl;
	env.end();
	return sol;
}

Solution benders_Master_integer(Instance inst, vector<vector<float> > scenarios){
	Solution sol;
	sol.inst = inst;

	//cpo model creation and solving

	IloEnv env;
	IloModel model(env);

	//setup cost
	int cP = 2;


	//vars
	IloNumVar z(env, -IloInfinity, IloInfinity);
	//IloNumVar z(env, -1000, 70);
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

	//cout<<"pouet1"<<endl;

	//consts

	for(int t = 1; t<inst.T+1;t++){
		model.add(X[t-1]<inst.X[t-1]);
	}

	for(int t = 0; t<inst.T; t++){
		for(int o = 0; o<scenarios.size();o++){
			model.add(B[o][t] - I[o][t] == scenarios[o][t] - X[t]); //(2)
			IloExpr expr(env);
			for(int i = 0; i<=t; i++){
				expr += s[o][i];
			}
			model.add(expr == scenarios[o][t]-B[o][t]);			//(3)
		}
	}

	// model.add(X[inst.T-1]==14);

	//cout<<"pouet2"<<endl;
	for(int o = 0; o<scenarios.size();o++){
		IloExpr expr(env);
		for(int t = 0; t<inst.T; t++){
			expr += (inst.cI*I[o][t] + inst.cB*B[o][t] - inst.bP*s[o][t]+cP*y[t]);
		}
		model.add(z>= expr);
	}


	//constraints for y............
	float M = inst.Dt[inst.T-1] + inst.Gamma;
	model.add(X[1]<=y[1]*M);

	for(int t = 1; t<inst.T; t++){
		model.add(X[t]-X[t-1]<=y[t]*M);
	}

	//obj
	model.add(IloMinimize(env, z));

	//cout<<"pouet3"<<endl;

	//solve
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

	// for(int o = 0; o<scenarios.size();o++){
	// 	for(int t = 0; t<inst.T; t++){
	// 		cout<<cplex.getValue(B[o][t])<<" "<<cplex.getValue(I[o][t])<<" "<<scenarios[o][t]<<" "<<cplex.getValue(X[t])<<" "<<cplex.getValue(s[o][t])<<endl;
	// 	}
	// }

	// float debug = 0;
	// for(int o = 0; o<scenarios.size();o++){
	// 	for(int t = 0; t<inst.T; t++){
	// 		//debug += inst.cI*cplex.getValue(I[o][t]) + inst.cB*cplex.getValue(B[o][t]) - inst.bP*cplex.getValue(s[o][t]);
	// 		debug += - inst.bP*cplex.getValue(s[o][t]);
	// 	}
	// }
	// cout<<"solution cost (cplex issue ?): "<<debug<<endl;
	env.end();
	return sol;
}


//solve subproblem using CPLEX: NOT USED, probably not working
Solution_ADV benders_Subproblem(Solution sol){
	Solution_ADV sol_adv;

	IloEnv env;
	IloModel model(env);

	//used to retrieve relevant pi after the solve
	vector<vector<bool> > pibool;
	pibool.resize(sol.inst.T+2);

	//vars
	IloArray<IloNumVarArray> pi(env,sol.inst.T+2);
	for(int t = 0; t<sol.inst.T+2;t++){
		//for t = 0 and t = T+1 only one variable in the array
		pi[t] = IloNumVarArray(env, sol.inst.Gamma+1);
		pibool[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			char name[80];
			pi[t][i] = IloNumVar(env, -IloInfinity, IloInfinity);
			sprintf(name,"pi_%d_%d",t,i);
			pi[t][i].setName(name);
		}
	}

	//csts
	//for t in 1...T-1
	for(int t = 1; t<sol.inst.T;t++){
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			for(int j = i; j<sol.inst.Gamma+1; j++){
				//the arc is in the graph =>
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
	//for t = T
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

	//obj

	//model.add(IloMinimize(env, pi[sol.inst.T][0]));

	IloExpr expr(env);
	for(int t = 0; t<sol.inst.T+2;t++){
		//for t = 0 and t = T+1 only one variable in the array
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			if(pibool[t][i]){
				expr += pi[t][i];
			}
		}
	}

	model.add(IloMinimize(env, expr)); //bon objectif pour etre sur que les potentiels collent leur borne

	IloCplex cplex(model);
	//cplex.exportModel ("lpex1.lp");
	cplex.setParam(IloCplex::Param::MIP::Display, 1);
	cplex.setOut(env.getNullStream());
    if ( !cplex.solve() ) {
    	env.error() << "Failed to optimize LP." << endl;
    	throw(-1);
	}

	//retrieve value for pi var
	vector<vector<float> > pisol;
	pisol.resize(sol.inst.T+2);
	for(int t = 0; t<sol.inst.T+2;t++){
		//for t = 0 and t = T+1 only one variable in the array
		pisol[t].resize(sol.inst.Gamma+1);
		for(int i = 0; i<sol.inst.Gamma+1; i++){
			if(pibool[t][i]){
				pisol[t][i] = cplex.getValue(pi[t][i]);
				// cout<<t<<" "<<i<<" "<<pisol[t][i]<<endl;
			}
		}
	}

	// string buff;
	// for(int i = sol.inst.Gamma; i>=0; i--){
	// 	//for t = 0 and t = T+1 only one variable in the array
	// 	buff = ' ';
	// 	for(int t = 0; t<sol.inst.T+2;t++){
	// 		if(pibool[t][i]){
	// 			buff+= ' ' + to_string(pisol[t][i]);
	// 		}
	// 	}
	// 	cout<<buff<<endl;
	// }



	// pisol[0][0] = cplex.getValue(pi[0][0]);
	// pisol[sol.inst.T][0] = cplex.getValue(pi[sol.inst.T][0]);

	//retrieve uncertainty from pi var (where constraints are tights)

	cout<<"======================LONGEST PATH : "<<pisol[sol.inst.T+1][0]<<endl;
	float previous_val = pisol[sol.inst.T+1][0];
	//=================DANGEROUS, WORKS ONLY IF GAMMA IS FULLY USED, SHOULD CHANGE THIS====================
	int previous_budget = sol.inst.Gamma;

	vector<int> offset;
	offset.resize(sol.inst.T);

	for(int t = sol.inst.T; t>=1; t--){
		for(int i=previous_budget; i>=0; i--){

			if(t == sol.inst.T){
				if( pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (previous_budget-i))) - sol.inst.bP*(sol.inst.Dt[t] - (previous_budget-i))) { //worse for deltat negative
					offset[t-1] = -(previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}
				else if(pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cB*(sol.inst.Dt[t] + (previous_budget-i) - sol.Xt[t]) - sol.inst.bP*sol.Xt[t]) { //worse for deltat negative
					offset[t-1] = (previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
				}
			}
			else if (t == 1) {
				if( pisol[t][previous_budget] - pisol[0][0] == sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (previous_budget-i)))){ //worse for deltat negative
					offset[t-1] = -(previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}
				else if(pisol[t][previous_budget] - pisol[0][0] == sol.inst.cB*(sol.inst.Dt[t] + (previous_budget-i) - sol.Xt[t])) { //worse for deltat negative
					offset[t-1] = (previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}
			}
			//generic case
			else{
				if( pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cI*(sol.Xt[t]- (sol.inst.Dt[t] - (previous_budget-i)))){ //worse for deltat negative
					offset[t-1] = -(previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
				}
				else if(pisol[t][previous_budget] - pisol[t-1][i] == sol.inst.cB*(sol.inst.Dt[t] + (previous_budget-i) - sol.Xt[t])) { //worse for deltat negative
					offset[t-1] = (previous_budget - i);
					// cout<<"tight constraint found for t = "<<t<<" from " <<previous_budget<< " to "<<i<<endl;
					previous_budget = i;
					break;
					
				}
			
			}
		}

	}

	// cout<<"maxindex :";
	// display_vector_int(maxindex);

	vector<float> Dt;
	Dt.resize(sol.inst.T);
	for(int t = 0; t<sol.inst.T; t++){
		Dt[t] = sol.inst.Dt[t] + offset[t];
	}

	env.end();

	// cout<<"hope this work : "<<endl;
	// display_vector_int(offset);
	// display_vector_float(Dt);

	sol_adv.Dt = Dt;
	return sol_adv;
}

//solve subproblem with dynamic prog
Solution_ADV benders_Subproblem_DP(Solution sol){
	Solution_ADV sol_adv;
	vector<vector<vector<vector<int> > > > arcbool; //bool flag to arcs within the worsts scenarios
	vector<vector<float> > pi_value; //value of the longest path to pi[t][j]
	vector<vector<bool> > pi_subopt_bool;
	vector<vector<vector<vector<float> > > > costs = budget_graph_cost(sol); //costs of all arcs of the budget graph

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

	//dynamic prog. for longest path

	float tmp;
	pi_value[0][0] = 0;
	for(int t=1; t<sol.inst.T+1;t++){
		for(int j = 0; j<sol.inst.Gamma+1; j++){
			tmp = pi_value[t-1][j]+costs[t][j][j][0];//init of pi_value[t][j]
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
		// cout<<"=============="<<t<<" "<<pi_value[t][0]<<endl;
	}

	tmp = pi_value[sol.inst.T][0];
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		if(pi_value[sol.inst.T][i]>tmp){
			tmp = pi_value[sol.inst.T][i];
		} 
		// cout<<i<<", "<<pi_value[sol.inst.T][i]<<endl;
	}
	pi_value[sol.inst.T+1][0] = tmp;
	// cout<<"longest path : "<<pi_value[sol.inst.T+1][0]<<endl;
	

	//==========================now the backtrack

	// cout<<"LONGEST PATH VALUE : "<<pi_value[sol.inst.T+1][0]<<endl;;

	vector<float> scenario;
	scenario.resize(sol.inst.T);

	//t = T+1
	pi_subopt_bool[sol.inst.T+1][0] = true;
	//cout<<"poeut"<<endl;
	for(int i = 0; i<sol.inst.Gamma+1;i++){
		// cout<<i<<" "<<pi_value[sol.inst.T][i]<<" "<<sub_OPT<<endl;
		if(pi_value[sol.inst.T][i]==pi_value[sol.inst.T+1][0]){
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
					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][0]){
						arcbool[t][i][j][0] = 1;
						pi_subopt_bool[t-1][i] = true;
						scenario[t-1] = sol.inst.Dt[t-1] - (j-i);
						// cout<<"local cost from : "<<t-1<<","<<i<<" to "<<t<<","<<j<<" : "<<costs[t][i][j][0]<<" : "<<endl;
						
						break;
						// cout<<"test passed : "<<t<<" "<<j<<" -> "<<t-1<<" "<<i<<" "<<"0"<<endl;
					}
					if(pi_value[t][j] == pi_value[t-1][i]+costs[t][i][j][1]){
						arcbool[t][i][j][1] = 1;
						pi_subopt_bool[t-1][i] = true;
						scenario[t-1] = sol.inst.Dt[t-1] + (j-i);
						// cout<<"local cost from : "<<t-1<<","<<i<<" to "<<t<<","<<j<<" : "<<costs[t][i][j][1]<<endl;
						// cout<<sol.inst.cB<<" "<<sol.inst.Dt[t]<<" "<<(j-i)<<" "<<-sol.Xt[t]<<" : "<<sol.inst.cB*(sol.inst.Dt[t] + (j-i) - sol.Xt[t])<<endl;
						break;
						// cout<<"test passed : "<<t<<" "<<j<<" -> "<<t-1<<" "<<i<<" "<<"1"<<endl;
					}
				}
			}
		}
	}

	// cout<<"=================BEGIN TEST=================="<<endl;

	// Solution test = KC_benders_Master(sol.inst, arcbool);

	// cout<<"new sol (with Graph LP): ";
	// display_vector_float(test.Xt);
	// cout<<"new sol value : "<<test.obj_val<<endl;

	// cout<<"=================END TEST=================="<<endl;


	//display the subgraph
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
	sol_adv.Dt = scenario;
	return sol_adv;
}


pair<int, float> benders_Main(Instance inst){

	auto start = high_resolution_clock::now();

	Solution sol;
	Solution new_sol;
	Solution_ADV sol_adv;
	bool stopCriterion = false;
	float proc_time;
	vector<vector<float> > scenarios;
	scenarios.resize(0);
	//only nominal scenario
	// cout<<"nominal scenario : ";
	// vector<float> debug = {0,0,0,0,1,1,2,4,5,5,6,9,9,11,11};
	// scenarios.push_back(debug);
	// display_vector_float(inst.Dt);
	scenarios.push_back(inst.Dt);

	sol = benders_Master(inst, scenarios);
	// sol = benders_Master_integer(inst, scenarios);

	// cout<<"first sol cost : "<<sol.obj_val<<endl;
	// cout<<"first sol : ";
	// display_vector_float(sol.Xt);

	int i = 2;
	while(!stopCriterion){
		// cout<<"=============ITERATION "<<i<<endl;
		sol_adv = benders_Subproblem_DP(sol);
		// sol_adv = benders_Subproblem(sol);
		scenarios.push_back(sol_adv.Dt);
		// cout<<"NEW SCENARIO :"<<endl;
		// display_vector_float(sol_adv.Dt);
		new_sol = benders_Master(inst, scenarios);
		// new_sol = benders_Master_integer(inst, scenarios);

		// cout<<"new sol : ";
		// display_vector_float(new_sol.Xt);
		i++;
		// cout<<"Obj Values of new sol and prec sol on new scenario :  "<< objective_value(new_sol, sol_adv.Dt) << " <= " << objective_value(sol, sol_adv.Dt)<<endl;
		if(objective_value(new_sol, sol_adv.Dt) == objective_value(sol, sol_adv.Dt)){
			stopCriterion = true;
			break;
		}

		
		// cout<<"new sol value : "<<new_sol.obj_val<<endl;
		// cout<<"worst case : ";
		// display_vector_float(sol_adv.Dt);
		sol = new_sol;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	float accuracy = (1./100000);
	proc_time = accuracy*float(duration.count());
	return make_pair(i, proc_time);
}

//===========================================main

vector<string> list_dir(const char *path) {
vector<string> allfile;
   struct dirent *entry;
   DIR *dir = opendir(path);
   
   if (dir == NULL) {
      return allfile;
   }
   while ((entry = readdir(dir)) != NULL) {
   // cout << entry->d_name << endl;
   allfile.push_back(entry->d_name);
   }
   closedir(dir);
   return allfile;
}

int main(int argc, const char* argv[]){
	//cout<<"test"<<endl;
	// pair<Solution, Solution_ADV> benders_sol;
	pair<int, float> benders_sol;
	float approx_coeff;
	float time, timeKC;
	int iter, iterKC;
	ofstream result;
  	result.open ("result.txt");
  	string filename;
  	vector<string> filelist = list_dir("/home/tportole/Documents/These/benders_kc/Parsed_Large_Instances");
  	int nbInst = filelist.size() - 2;
  	Instance inst;
  	vector<int> debug;
  	vector<int> debug2;

  	int seed = 31415;
  	srand (seed);

	for(int Gamma = 1; Gamma<100; Gamma+=20){
		for(int tau=0; tau <11; tau +=2){
			
			iter = 0;
			iterKC = 0;
			time = 0;
			timeKC=0;
			debug.resize(0);
			debug2.resize(0);
			for(int i = 2; i<nbInst; i++ ){

				filename = "Parsed_Large_Instances/" + filelist[i];
				cout<<filename<<" "<<Gamma<<" "<<tau<<" "<<endl;
				//cout<<"pouet1"<<endl;
				inst = read_instance_randomized(filename, Gamma);
				//cout<<"pouet2"<<endl;
				//display_vector_float(inst.Dt);
				//display_vector_float(inst.deltat);
				benders_sol = benders_Main(inst);
				cout<<"STANDARD done"<<endl;
				// debug2.push_back(benders_sol.first);
				// iter += benders_sol.first;
				iter += benders_sol.first;
				time += benders_sol.second;

				approx_coeff = float(tau)/10;
				benders_sol = KC_benders_Main(inst, approx_coeff);
				cout<<"KC done"<<endl;
				// debug.push_back(benders_sol.first);
				// iterKC += benders_sol.first;
				iterKC+= benders_sol.first;
				timeKC+= benders_sol.second;

				// timeKC += benders_sol.second;

				// std::this_thread::sleep_for(std::chrono::milliseconds(1000));
			}
			result<<"STANDARD "<<Gamma<<" "<<tau<<" ";
			result<<float(iter)/nbInst<<" "<<float(time)/nbInst<<endl;
			// display_vector_int(debug2);
			result<<"KC "<<Gamma<<" "<<tau<<" ";
			result<<float(iterKC)/nbInst<<" "<<float(timeKC)/nbInst<<endl;
			// display_vector_int(debug);
			
		}	
	}

	result.close();
	//cout<<"test"<<endl;
}