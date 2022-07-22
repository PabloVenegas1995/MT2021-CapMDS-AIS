/***************************************************************************
                         brkga.cpp  -  description
                             -------------------
    begin                : Wed Sept 4 2019
    copyright            : (C) 2019 by Christian Blum
    email                : christian.blum@iiia.csic.es
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "callbacks.hpp"

#ifndef CALLBACKS_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Timer.h"
#include "Random.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <list>
#include <set>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <limits>
#include <filesystem>
#include <ilcplex/ilocplex.h>
namespace fs = std::filesystem;

#endif


struct Option {
	int vertex;
	double value;
};

struct Solution {
	int score;
	set<int> vertices;
};

struct Individual {
	vector<double> vec;
	double ofv;
	set<int> selected_vertex;
};

struct ValPos {
	double val;
	int pos;
};

time_t t;
Random* rnd;
int verbose =0;
// instance data
int n_of_vertices;
vector<int> capacity;
//vector<vector<int> > matrix;
vector<set<int> > adj;
bool uniform_capacities = false;
bool tuning =false;
string file_directory;
string inputFile="none";
//BRKGA parameters
int population_size = 30;
double elite_proportion = 0.15; // normally between 0.1 and 0.25
double mutant_proportion = 0.20; // normally between 0.1 and 0.3
double elite_inheritance_probability = 0.7; // normally greater than 0.5 and <= 0.8

//general parameters
double computation_time_limit = 1000.0;

//Barrakuda-ILP parameters
int n_barrakuda_sols = 1; //Debe ser estrictamente menor al número de elites
int mip_emphasis = 1; //cambio de comportamiento de solver
double CPLEX_t_limit = 3600.0;

int m_limit = 10000000;
double  best_solution = std::numeric_limits<int>::max();
double best_time      = std::numeric_limits<int>::max();
double  integrated_best_solution = std::numeric_limits<int>::max();
double integrated_best_time;
bool ilp_feedback_disable = false;
int feeder_mode = 0; //0 by default: toma los "n_barrakuda_sols" mejores //1: toma el best + n_barrakuda_sols -1 random cromosomas.

// Parametros nuevos
int n_chunks = -1;
double Chunck_CPLEX_t_limit = CPLEX_t_limit/n_chunks;
int brkga_it = 1;

// #0 Cambiar feeder_mode si la instancia no cambia
int flag_0 = 0; 
int num_iguales = 0;
int max_iguales = -1;
int input_feedermode = feeder_mode;

// #1 Aumentar n_barrakuda_sols si encuentra el optimo muy rapido
int flag_1 = 0;

// #2  Disminuir n_barrakuda_sols si se demora mucho en mejorar el incumbent
int flag_2 = 0;

// #3 Abortar ejecucion de CPLEX si pasan muchos chunks sin acercarse a la mejor solucion conocida
int flag_3 = 0;
double max_percent = -1;

// #4 Disminuir n_barrakuda_sols si se demora mucho en mejorar el gap
int flag_4 = 0;

// Random seed
int seed = 8;

ILOMIPINFOCALLBACK2(ramLimitCallback,
	IloBool,  aborted,
	int , limit)
{

	int vmrss_kb;
	ifstream infile("/proc/self/status");
	string line;
	bool found = false;
	while (not found and getline(infile, line)) {
		istringstream iss(line);
		string s;
		iss >> s;
		if (s == "VmRSS:") {
			found = true;
			iss >> vmrss_kb;
			if (vmrss_kb >= limit) abort();
		}
	}
}

//EXPERIMENTAL PROCEDURES CHECK IT!!


//Usado antes por errores de aproximación...¡será necesario???
bool verNum(double num)
{
	bool isInt;

	if(round(num)-num==0.0) isInt=true; else isInt=false;
	if (verbose && ! isInt) cout <<"Error de CPLEX"; 

	return isInt;
}

ILOMIPINFOCALLBACK6(loggingCallback,
	Timer&, timer,
	vector<double>&, results,
	vector<double>&, times,
	vector<double>&, gaps,
	IloNum,         lastIncumbent,
	int, iter)
{

	IloNum nv = getIncumbentObjValue();
	double newTime = timer.elapsed_time(Timer::VIRTUAL);
	double newGap = 100.0*getMIPRelativeGap();
	if (not (lastIncumbent == nv)) {
        //cout << "best " << int(nv) << "\ttime " << newTime << "\tgap " << newGap << endl;
        if(int(nv)<best_solution && int(nv)>0 && verNum(nv)) {  //ORIGINAL
        //if(int(round(nv))<best_solution && int(nv)>0) {  //Modificado 
           //newline
            best_solution=int(nv); //ORIGINAL
            best_time=newTime;

            if(best_solution<integrated_best_solution)
            {
            	if(!tuning) cout << "Mejor solucion en CPLEX: " << best_solution 
            	<< "| tiempo = " << newTime << "\n";
            	integrated_best_solution=best_solution;
            	integrated_best_time=newTime;
            	if(verbose)cout<<"best: "<<integrated_best_solution<<" "<<integrated_best_time<<endl;
            }
            //best_solution=int(round(nv));  //modificado


            //********************************************************************************************
            //COMENTADO PARA TUNING
            //********************************************************************************************
            if(!tuning) cout<<"best "<< int(best_solution)<<"\ttime " <<newTime<<"\tgap "<<newGap<<"\tnv "<<nv<<endl;
            //********************************************************************************************

        }
        results[iter] = double(nv);  //ORIGINAL
       // results[iter] = double(round(nv));  //modificado
        
        times[iter] = newTime;
        gaps[iter] = newGap;
    }
    lastIncumbent = nv;
}

bool vertex_set_comparison(const set<int> &A, const set <int> &B)
{
    //bool hipo=true;

	if(A.size()!=B.size()) { return false;}
	set<int>::iterator itB=B.begin();
	for(set<int>::iterator itA=A.begin();itA!=A.end();itA++)
	{ 
		if (*itA!=*itB) {return false;}
		itB++;
	}

	return true;
}

bool option_compare(const Option& o1, const Option& o2) {

	return (o1.value < o2.value);
}

bool individual_compare(const Individual& i1, const Individual& i2) {

	return i1.ofv < i2.ofv;
}

inline int stoi(string &s) {

	return atoi(s.c_str());
}

inline double stof(string &s) {

	return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

	int iarg = 1;;


	while (iarg < argc) {
		if (strcmp(argv[iarg],"-i")==0) inputFile=(argv[++iarg]);
		if (strcmp(argv[iarg],"-dir")==0) file_directory = argv[++iarg];
		else if (strcmp(argv[iarg],"-p")==0) population_size = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-unicap")==0) uniform_capacities = true;
		else if (strcmp(argv[iarg],"-pe")==0) elite_proportion = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-pm")==0) mutant_proportion = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-rhoe")==0) elite_inheritance_probability = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-t")==0) computation_time_limit = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-verbose")==0) verbose = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-n_barrakuda_sols")==0) n_barrakuda_sols = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-mipemph")==0) mip_emphasis = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-cpl_time")==0) CPLEX_t_limit = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-mlim")==0) m_limit = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-ilp_feedback_disable")==0) ilp_feedback_disable=true;
		else if (strcmp(argv[iarg],"-feeder_mode")==0) feeder_mode= atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-tuning")==0) tuning=true;
		else if (strcmp(argv[iarg],"-n_chunks")==0) n_chunks = atoi(argv[++iarg]); //Nuevos
		else if (strcmp(argv[iarg],"-brkga_it")==0) brkga_it = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-pol_0")==0) max_iguales = atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-pol_1")==0) flag_1= atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-pol_2")==0) flag_2= atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-pol_3")==0) max_percent = atof(argv[++iarg]);
		else if (strcmp(argv[iarg],"-pol_4")==0) flag_4= atoi(argv[++iarg]);
		else if (strcmp(argv[iarg],"-seed")==0) seed = atoi(argv[++iarg]);

		iarg++;
	}
}

int produce_random_integer(int max) {

	int num = int(double(max) * rnd->next());
	if (num == max) num = num - 1;
	return num;
}

int get_random_element(const set<int>& s) {

	double r = produce_random_integer(int(s.size()));
	set<int>::iterator it = s.begin();
	advance(it, r);
	return *it; 
}

void check_solution(vector<int>& covered_by, vector<bool>& chosen, Individual& ind) {

    /*
    cout << "Solution:";
    for (int i = 0; i < n_of_vertices; ++i) if (chosen[i]) cout << " " << i;
    cout << endl;
    for (int i = 0; i < n_of_vertices; ++i) {
        if (not chosen[i]) cout << "Node " << i << " is covered by: " << covered_by[i] << endl;
    }
    */

	int score = 0;
	vector<int> number_of_covered(n_of_vertices, 0);
	for (int i = 0; i < n_of_vertices; ++i) {
		if (chosen[i]) ++score;
		else {
			if (covered_by[i] == -1) cout << "Node " << i << " is not covered." << endl;
			else {
				++number_of_covered[covered_by[i]];
				if (not chosen[covered_by[i]])
					cout << "Node " << i << " is covered by node " << covered_by[i] << " which is NOT chosen for the solution." << endl;
			}
		}
	}

	for (int i = 0; i < n_of_vertices; ++i) {
		if (chosen[i]) {
			if (capacity[i] < number_of_covered[i])
				cout << "Node " << i << " has capacity " << capacity[i] << " and covers " << number_of_covered[i] << " nodes." << endl;
		}
	}
}

void show_individual(Individual & ind)
{
	for(int i=0;i<n_of_vertices;i++)
		cout<<ind.vec[i]<<",";

	cout<<endl;
	for(int i=n_of_vertices;i<2*n_of_vertices;i++)
		cout<<ind.vec[i]<<",";

	cout<<endl;
	cout<<"cost: "<<ind.ofv<<endl;
	cout<<"solution components: "<<ind.selected_vertex.size()<<endl;
	for(set<int>::iterator it=ind.selected_vertex.begin();it!=ind.selected_vertex.end();it++)
		cout<<*it<<",";

	cout<<endl;
}

void integrate(Individual &ind, set<int> &vertex_in_solution)
{
	for(int i=0;i<n_of_vertices;i++){
		ind.vec[i]=0.0;
		ind.vec[i+n_of_vertices]=0.5;
	}

	for(set<int>::iterator it=vertex_in_solution.begin();it!=vertex_in_solution.end();it++)
		ind.vec[*it]=1.0;
}

void evaluate(Individual& ind) {

	int score = 0;
	set<int> chosen_vertices;
	vector<bool> allready_covered(n_of_vertices, false);
	vector<bool> chosen(n_of_vertices, false);
	vector<int> covered_by(n_of_vertices, -1);
	int num_nodes_uncovered = n_of_vertices;

	vector<set<int> > uncovered_neighbors = adj;

	set<int> candidates;
	for (int i = 0; i < n_of_vertices; ++i) candidates.insert(i);

		while (num_nodes_uncovered > 0)
		{
			int max_val = -1.0;
			int chosen_vertex = -1;
			for (set<int>::iterator cit = candidates.begin(); cit != candidates.end(); ++cit) {
				int h_func = int(uncovered_neighbors[*cit].size());

				if (capacity[*cit] < h_func) h_func = capacity[*cit];

	            double cit_value = (double(h_func) + 1.0)*(ind.vec)[*cit];  //Unifica valor heurístico con valor biased
	            if (cit_value > max_val) {
	            	max_val = cit_value;
	            	chosen_vertex = *cit;
	            }
        	}

	        chosen_vertices.insert(chosen_vertex);
	        chosen[chosen_vertex] = true;
	        if (not allready_covered[chosen_vertex]) {
	        	allready_covered[chosen_vertex] = true;
	        	--num_nodes_uncovered;
	        }
        	score += 1;

	        if (capacity[chosen_vertex] >= int(uncovered_neighbors[chosen_vertex].size())) {
	        	num_nodes_uncovered -= int(uncovered_neighbors[chosen_vertex].size());
	        	set<int> to_cover = uncovered_neighbors[chosen_vertex];
	        	for (set<int>::iterator sit = to_cover.begin(); sit != to_cover.end(); ++sit) {
	        		allready_covered[*sit] = true;
	        		covered_by[*sit] = chosen_vertex;

	        		for (set<int>::iterator ssit = adj[*sit].begin(); ssit != adj[*sit].end(); ssit++)
	        			uncovered_neighbors[*ssit].erase(*sit);
	        	}

	        	uncovered_neighbors[chosen_vertex].clear();
	    	}
	    	else
	    	{
	    		num_nodes_uncovered -= capacity[chosen_vertex];
	    		vector<Option> to_cover;
		    	for (set<int>::iterator sit = uncovered_neighbors[chosen_vertex].begin(); sit != uncovered_neighbors[chosen_vertex].end(); ++sit) {
		    		Option opt;
		    		opt.vertex = *sit;
		    		opt.value = double(uncovered_neighbors[*sit].size())*(ind.vec)[*sit + n_of_vertices];
		    		to_cover.push_back(opt);
		    	}

	    		sort(to_cover.begin(), to_cover.end(), option_compare);
		    	for (int cc = 0; cc < capacity[chosen_vertex]; ++cc) {
		    		int el = to_cover[cc].vertex;
		    		uncovered_neighbors[chosen_vertex].erase(el);
		    		allready_covered[el] = true;
		    		covered_by[el] = chosen_vertex;

		    		for (set<int>::iterator ssit = adj[el].begin(); ssit != adj[el].end(); ssit++)
		    			uncovered_neighbors[*ssit].erase(el);
		    	}
			}

			for (set<int>::iterator sit = adj[chosen_vertex].begin(); sit != adj[chosen_vertex].end(); sit++)
				uncovered_neighbors[*sit].erase(chosen_vertex);

			candidates.erase(chosen_vertex);
			set<int> to_delete;
			for (set<int>::iterator cit = candidates.begin(); cit != candidates.end(); ++cit) {
				if (int(uncovered_neighbors[*cit].size()) == 0 and allready_covered[*cit])
					to_delete.insert(*cit);
			}

			for (set<int>::iterator sit = to_delete.begin(); sit != to_delete.end(); ++sit)
				candidates.erase(*sit);
		}

	ind.selected_vertex=chosen_vertices;
	ind.ofv = score;
    //check_solution(covered_by, chosen, ind);
}

void generate_random_solution(Individual& ind) {
	ind.vec = vector<double>(2*n_of_vertices);
	for (int i = 0; i < 2*n_of_vertices; ++i) (ind.vec)[i] = rnd->next();
		evaluate(ind);
}


void show_population(vector<Individual> &pop, int n_e, int n_m, int n_o)
{   
	for(int j=0;j<population_size;j++)
	{
		cout<<j<<"--";
		for(int i=0; i<n_of_vertices*2;i++)
		{
			cout<<pop[j].vec[i]<<", ";
			if(i==n_of_vertices-1) cout<<" - ";
		}

		cout<<"("<<pop[j].ofv<<")--Sol: ";

		for(set<int>::iterator it=pop[j].selected_vertex.begin();it!=pop[j].selected_vertex.end();it++)
			cout<<*it<<"-";

		if(j<n_e) cout <<" E";
		else if(j>=n_e && j<n_m+n_e) cout <<" M";
		else if(j>=n_m+n_e && j<population_size-1) cout <<" O";
		else cout <<" ILP";

		cout<<endl;
	}
	//exit(0);
}

void show_vertex_in_c(set<int> & v){
	for(set<int>::iterator it=v.begin();it!=v.end();it++)
		cout <<*it<<"-";

	cout<<endl;
}

/**********
Main function
**********/

int main( int argc, char **argv ) {

	read_parameters(argc,argv);

	input_feedermode = feeder_mode;
	if(max_iguales > 0) flag_0 = 1;
	if(max_percent > 0) flag_3 = 1;
	if(n_chunks < 0) n_chunks = 1;
	Chunck_CPLEX_t_limit = CPLEX_t_limit/n_chunks;

	//Check parameters	
	/*
	if(flag_0) cout << "Politica 0 activada, max_iguales = " << max_iguales << "\n";
	if(flag_1) cout << "Politica 1 activada\n";
	if(flag_2) cout << "Politica 2 activada\n";
	if(flag_3) cout << "Politica 3 activada, max_percent = " << max_percent << "\n";
	if(flag_4) cout << "Politica 4 activada\n";
	cout << "n_chunks = " << n_chunks << "\n";
	cout << "brkga_it = " << brkga_it << "\n";
	*/

	std::cout << std::setprecision(3) << std::fixed;

	//rnd = new Random((unsigned) time(&t));
	rnd = new Random(seed);
	rnd->next();

    vector<double> results;// (5, 0.0); //5 es arbitrario
    vector<double> times; //(5, 0.0);
    vector<double> gaps; // (((5, 0.0);

    int iter = 0;

  	//for (const auto & entry : fs::directory_iterator(file_directory)) {

    results.push_back(0.0);
    times.push_back(0.0);
    gaps.push_back(0.0);

	// string file_name = entry.path();
    ifstream indata;

    if(inputFile!="none"){
    	indata.open(inputFile.c_str());
    }
    int cap_captured;
    indata >> n_of_vertices >> cap_captured;    
    if(cap_captured==0) uniform_capacities=false; else uniform_capacities=true;
    indata.close();



    if(inputFile!="none"){
    	indata.open(inputFile.c_str());
    }

   	// else
   	// {     
   	// indata.open(file_name.c_str());
   	// }
    //cout << "start file " << file_name << endl; //--silenciado tuning
    if(!indata) { // file couldn't be opened
    	cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_vertices;
    capacity = vector<int>(n_of_vertices);
    if (uniform_capacities) {
    	int uniform_capacity;
    	indata >> uniform_capacity;
    	for (int i = 0; i < n_of_vertices; ++i) capacity[i] = uniform_capacity;
    }
	else {
		for (int i = 0; i < n_of_vertices; ++i) {
			int v, cap;
			indata >> v >> cap;
			capacity[v] = cap;
		}
	}

    //matrix = vector< vector<int> >(n_of_vertices, vector<int>(n_of_vertices, 0));
	adj = vector< set<int> >(n_of_vertices);

	int v1, v2;
	while (indata >> v1 >> v2) {
		adj[v1].insert(v2);
		adj[v2].insert(v1);
	}
	indata.close();

	// el tiempo empieza
	Timer timer;


	int best_sol_value = std::numeric_limits<int>::max();
	best_solution =std::numeric_limits<int>::max();
	best_time      =std::numeric_limits<int>::max();

	int n_elites = int(double(population_size)*elite_proportion);
	if (n_elites < 1) n_elites = 1;

	int n_mutants = int(double(population_size)*mutant_proportion);
	if (n_mutants < 1) n_mutants = 1;

	int n_offspring = population_size - n_elites - n_mutants;
	if (n_offspring < 1) {
		cout << "OHOHOH: wrong parameter settings" << endl;
		exit(0);
	}

	if(verbose){//configuration
    	cout<<"Tamaño de población "<<population_size<<endl;
    	cout<<"n_elites: "<<n_elites<<endl;
    	cout<<"n_mutants: "<<n_mutants<<endl;
    	cout<<"n_offspring: "<<n_offspring<<endl;
    }

    double ctime;
    vector<Individual> population(population_size);
    for (int pi = 0; pi < population_size; ++pi)
    {
    	generate_random_solution(population[pi]);
    	if (population[pi].ofv < best_sol_value) {
    		best_sol_value = population[pi].ofv;
    		ctime = timer.elapsed_time(Timer::VIRTUAL);
    		results[iter] = best_sol_value;
    		times[iter] = ctime;
         	if(!tuning) cout << "best " << best_sol_value << "\ttime " << ctime << endl; //silneciado tuning

	        if(best_sol_value<integrated_best_solution){
	         	integrated_best_solution=best_sol_value;
	        	integrated_best_time=ctime;
	        }
    	}
	}

    //Población inicial creada y evaluada
    //if (verbose) show_population(population,n_elites,n_mutants,n_offspring);

	ctime = timer.elapsed_time(Timer::VIRTUAL);
	set <int> vertex_in_c_old;

	//CICLO PRINCIPAL DE ALGORITMO
	while (ctime < computation_time_limit)
	{
		for(int iciclo = 0; iciclo < brkga_it; iciclo++){
			//Copia de elites en nueva población
			sort(population.begin(), population.end(), individual_compare); //ORDENAMIENTO DE NUEVA POBLACION
			if(population[0].ofv<best_sol_value){
				best_sol_value=population[0].ofv;
				ctime = timer.elapsed_time(Timer::VIRTUAL);
				results[iter] = best_sol_value;
				times[iter] = ctime;
			    if(!tuning) cout << "best " << best_sol_value << "\ttime " << ctime << endl; //silenciado tuning

			    if(best_sol_value<integrated_best_solution)
			    {
			    	if(!tuning) cout << "Mejor solucion en BRKGA (ORDENAMIENTO): " << best_sol_value 
			    	<< "| tiempo = " << ctime << "\n";
			    	integrated_best_solution=best_sol_value;
			    	integrated_best_time=ctime;
			    }                

			}

			vector<Individual> new_population(population_size);
			for (int ic = 0; ic < n_elites; ++ic) {
				new_population[ic].vec = population[ic].vec;
				new_population[ic].ofv = population[ic].ofv;
				new_population[ic].selected_vertex =population[ic].selected_vertex;
			}

			//if (verbose) show_individual(population[0]); //mostrando un elemento de elite
			//creación de mutantes
			for (int ic = 0; ic < n_mutants; ++ic) {
				generate_random_solution(new_population[n_elites + ic]);
				if (new_population[n_elites + ic].ofv < best_sol_value) {
					best_sol_value = new_population[n_elites + ic].ofv;
					ctime = timer.elapsed_time(Timer::VIRTUAL);
					results[iter] = best_sol_value;
					times[iter] = ctime;
			        if(!tuning) cout << "best " << best_sol_value << "\ttime " << ctime << endl; //silenciado tuning

			        if(best_sol_value<integrated_best_solution)
			        {
			        	if(!tuning) cout << "Mejor solucion en BRKGA (MUTAR): " << best_sol_value 
			        	<< "| tiempo = " << ctime << "\n";
			        	integrated_best_solution=best_sol_value;
			        	integrated_best_time=ctime;
			        }
			    }
			}

			//creación de hijos por cruza
			for (int ic = 0; ic < n_offspring; ++ic) {
				int first_parent = produce_random_integer(n_elites);
				int second_parent = n_elites + produce_random_integer(population_size - n_elites);
				new_population[n_elites + n_mutants + ic].vec = vector<double>(2*n_of_vertices);
				for (int i = 0; i < 2*n_of_vertices; ++i) {
					double rnum = rnd->next();
					if (rnum <= elite_inheritance_probability) (new_population[n_elites + n_mutants + ic].vec)[i] = (population[first_parent].vec)[i];
					else (new_population[n_elites + n_mutants + ic].vec)[i] = (population[second_parent].vec)[i];
				}

				evaluate(new_population[n_elites + n_mutants + ic]);
				if (new_population[n_elites + n_mutants + ic].ofv < best_sol_value) {
					best_sol_value = new_population[n_elites + n_mutants + ic].ofv;
					ctime = timer.elapsed_time(Timer::VIRTUAL);
					results[iter] = best_sol_value;
					times[iter] = ctime;
			        if(!tuning) cout << "best " << best_sol_value << "\ttime " << ctime << endl; //silenciado tuning
			        if(best_sol_value<integrated_best_solution)
			        {
			        	if(!tuning) cout << "Mejor solucion en BRKGA (CRUZA): " << best_sol_value 
			        	<< "| tiempo = " << ctime << "\n";
			        	integrated_best_solution=best_sol_value;
			        	integrated_best_time=ctime;
			        }
			    }
			}

			population.clear();
			population = new_population;
			ctime = timer.elapsed_time(Timer::VIRTUAL);

			//Fin de iteración
			//if (verbose) show_population(population,n_elites,n_mutants,n_offspring);
			if (verbose) cout <<"end of iteration"<<endl;
		}

		set <int> vertex_in_c;
		set <int> vertex_in_c_complement;
		if(n_barrakuda_sols>n_elites) n_barrakuda_sols=n_elites;

		// creación de set para cplex
		if(feeder_mode==0)
		{    
		//agregando componentes a subinstancia
			for (int i=0;i<n_barrakuda_sols;i++)
			{
				for(set<int>::iterator it=population[i].selected_vertex.begin();it!=population[i].selected_vertex.end();it++)
					{vertex_in_c.insert(*it);}
			}
		}
		else if(feeder_mode==1)
		{
		    //agregando al mejor
			for(set<int>::iterator it=population[0].selected_vertex.begin();it!=population[0].selected_vertex.end();it++)
				{vertex_in_c.insert(*it);}

			for(int i=0;i<n_barrakuda_sols-1;i++)
			{
				int selected=produce_random_integer(population_size);
				for(set<int>::iterator it=population[selected].selected_vertex.begin();it!=population[selected].selected_vertex.end();it++)
					{vertex_in_c.insert(*it);}

			}
		}

		if(verbose) cout <<"Tamaño conjunto de componentes de solución para ILP: "<<vertex_in_c.size()<<endl;

		//Comparación para verificar si el conjunto ha cambiado
		if(vertex_set_comparison(vertex_in_c,vertex_in_c_old) && flag_0) {
		    //cout<<"IGUALES!!!!!!!!!!!!!!!!!: "<<vertex_in_c.size()<<" "<<vertex_in_c_old.size()<<endl;
        	num_iguales++;
        	if(num_iguales >= max_iguales){
            	//cout << "CAMBIO FEEDER MODE!!!!!!!!!!!!!!!\n";
            	feeder_mode = 1; // No elite
            }
		}
		else{
			//cout<<"DISTINTOS!!!!"<<endl;

			if(flag_0){
				feeder_mode = input_feedermode;
				num_iguales = 0;
			}

			//if(verbose) show_vertex_in_c(vertex_in_c);
			//Creating complement set of C'
			vertex_in_c_complement.clear();
			for(int i=0;i<n_of_vertices;i++)
				vertex_in_c_complement.insert(i);

			for(set<int>::iterator it=vertex_in_c.begin();it!=vertex_in_c.end();it++)
				vertex_in_c_complement.erase(*it);


			//*****************************************
			//LOADING STATIC PART OF THE MODEL
			//*****************************************

			IloEnv env;

			try{
				env.setOut(env.getNullStream());
				IloModel model(env);
				IloNumVarArray x(env, n_of_vertices, 0, 1, ILOINT);
				typedef IloArray<IloIntArray> IntArray;
				typedef IloArray<IloNumArray> NumMatrix;
				typedef IloArray<IloNumVarArray> NumVarMatrix;

				int ind=0;
				NumVarMatrix y(env,n_of_vertices);
				for(int i=0;i<n_of_vertices;i++)
					y[i]=IloNumVarArray(env,n_of_vertices,0,1,ILOINT);

				IloExpr obj(env);
				for (int i=0;i<n_of_vertices;i++) obj += x[i];
					model.add(IloMinimize(env, obj));

				obj.end();

				//***************************************
				//RESTRICCIONES ORIGINALES (basadas del ILP puro) + Restriccion 4 para CMSA
				//***************************************


	            // Restricción (1) "Cobertura a todos los nodos"
	            // Sum Yji + Xi >=1
	            // Sum Yji + Xi >=0.5
				for (int i=0; i<n_of_vertices;i++)
				{
					IloExpr expr(env);
					expr += -1*(x[i]);
					for(set<int>::iterator sit =adj[i].begin(); sit != adj[i].end(); ++sit)
					{
						expr += -1*y[*sit][i];
						//expr += -1*x[i];
					}
					model.add(expr == -1);
					expr.end();
				}


	            // Restriccion (2) "Restricción de capacidad"
	            // Sum Yij <=Xi*Cap(Vi), para todo Vi en V.
				for (int i=0; i<n_of_vertices;i++) 
				{

					IloExpr expr(env);
					expr += -1*capacity[i]*x[i];
					for(set<int>::iterator sit =adj[i].begin(); sit != adj[i].end(); ++sit)
					{
						expr += y[i][*sit];
					}
					model.add(expr <= 0);
					expr.end();
				}

	            //Restricción (3) "Relación entre variables"
	            // Yij<=Xi
	            //Yij-Xi<=0
	            //Yij-Xi<=0.4
				for(int i=0;i<n_of_vertices;i++)
				{      
					for (int j=0;j<n_of_vertices;j++)
					{           
						IloExpr expr(env);
						expr += -1*x[i];
						expr += y[i][j];
						model.add(expr <= 0);
						expr.end();              
					}
				}

	            //Ingresar nueva restricción
	            //IloExtractableArray new_constraints(env);
				for(set<int>::iterator it=vertex_in_c_complement.begin();it!=vertex_in_c_complement.end();it++)
				{ 
					IloExpr expr(env);
					expr+=x[*it];
					//IloRange r=(expr == 0);
					model.add(expr==0);
					expr.end();
					//new_constraints.add(r);
				}


	            //model.add(new_constraints);
	            //Resolver
				IloCplex cpl(model);

				//double time_aux = timer.elapsed_time(Timer::VIRTUAL);
				//CPLEX_t_limit = CPLEX_t_limit - (time_aux - ctime);

				cpl.setParam(IloCplex::TiLim, CPLEX_t_limit);   
				cpl.setParam(IloCplex::MIPEmphasis,mip_emphasis);
				cpl.setParam(IloCplex::NodeFileInd, 2);
				cpl.setParam(IloCplex::Threads, 1);
				cpl.use(ramLimitCallback(env, IloFalse, m_limit));

				IloNum lastObjVal = 100000000;
				cpl.use(loggingCallback(env, timer, results, times, gaps, lastObjVal, iter));


				if(verbose) cout << "Tiempo antes: " << timer.elapsed_time(Timer::REAL) << "\n";
				if(verbose) cout<<"S(init Solver)";

	            /*********************** CICLO DINAMICO ****************************/

	            /*cpl.solve();
	            double lastVal_aux = double(cpl.getObjValue());
	            double lastGap_aux = 100.0*cpl.getMIPRelativeGap();
	            if(lastGap_aux == 0){
	            	n_barrakuda_sols++;
	            	cout << "Aumenta n_barrakuda_sols: " << n_barrakuda_sols<< "\n";
		        }*/

	            //bool* abort = new bool;
	            //*abort = false;
	            //double start_time = timer.elapsed_time(Timer::VIRTUAL);
	            // CAMBIO
	            //cpl.use(simplexCallback(env, timer, start_time, abort));
	            // CAMBIO

				if(verbose){
					cout << "\n--------------------ENTRA-----------------\n";
					cout << "Tiempo chunck: " << Chunck_CPLEX_t_limit << "\n";
					cout << "Numero de chuncks: " << n_chunks << "\n";
					cout << "n_barrakuda_sols: " << n_barrakuda_sols << "\n";
					cout << "n_elites: " << n_elites << "\n";
					cout << "tamanio: " << vertex_in_c.size() << "\n\n";
				}

				double gap_anterior = std::numeric_limits<int>::max();
				double solucion_anterior = std::numeric_limits<int>::max();
				int cnt_2 = 0;
				int cnt_3 = 0;
				int cnt_4 = 0;
				for(int cnt_it = 0; cnt_it < n_chunks; cnt_it++)
				{
					cpl.setParam(IloCplex::TiLim, Chunck_CPLEX_t_limit);
					cpl.solve();
					double lastVal_aux = double(cpl.getObjValue()); //cambiar a incumbent
					double lastGap_aux = 100.0*cpl.getMIPRelativeGap();
					if(verbose) cout << "d " << cnt_it << "\nd " << "obj_value = " << lastVal_aux;
					if(verbose) cout <<"; " << "gap = " << lastGap_aux << "\n";

					if(verbose) cout << "d Tiempo entre: " << timer.elapsed_time(Timer::REAL) << "\n";

					// #1
					if(lastGap_aux == 0 && cnt_it <= n_chunks/4 && flag_1){
						n_barrakuda_sols++;
						if(verbose) cout << "d Aumenta n_barrakuda_sols (p#1): " << n_barrakuda_sols<< "\n";
						break;
					}

					// #2
					if(lastVal_aux < solucion_anterior){
						solucion_anterior = lastVal_aux;
						cnt_2 = 0;
					}else{
						cnt_2++;
						if(cnt_2 >= n_chunks/4 && flag_2){
							
							if(n_barrakuda_sols > 1)
								n_barrakuda_sols--;

							if(verbose) cout << "d Disminuye n_barrakuda_sols (p#2): " << n_barrakuda_sols<< "\n";
							break;
						}
					}

					// #3
					if(lastVal_aux >= best_solution*max_percent + best_solution){
						cnt_3++;
						if(cnt_3 >= n_chunks/2 && flag_3){
							if(verbose) cout << "d Aborta CPLEX (p#3)\n";
							break;
						}
					}else{
						cnt_3 = 0;
					}

					// #4
					if(lastGap_aux < gap_anterior){
						gap_anterior = lastGap_aux;
						cnt_4 = 0;
					}else{
						cnt_4++;
						if(cnt_4 >= n_chunks/4 && flag_4){

							if(n_barrakuda_sols > 1)
								n_barrakuda_sols--;
							
							if(verbose) cout << "d Disminuye n_barrakuda_sols (p#4): " << n_barrakuda_sols<< "\n";
							break;
						}
					}
					if(verbose) cout << "\n";
				}

				if(verbose){
					cout << "Tiempo despues: " << timer.elapsed_time(Timer::REAL) << "\n";
					cout << "--------------------SALE-----------------\n";
				}

				double lastVal_aux = double(cpl.getObjValue());
				double lastGap_aux = 100.0*cpl.getMIPRelativeGap();
				//double node_count = ;
				if(verbose) cout << "d obj_value = " << lastVal_aux << "\n";
				if(verbose) cout << "d gap = " << lastGap_aux << "\n";
				//cout << "node_count = " << node_count << "\n";*/

				double finish_time = timer.elapsed_time(Timer::VIRTUAL);

	            //cout<<"saliendo"<<endl;
	        	if(cpl.getStatus() ==IloAlgorithm::Infeasible)
	        	{
	        		cout<<"ERROR EN CPLEX"<<endl;
	        		exit(0);
	        	}

				if (cpl.getStatus() == IloAlgorithm::Optimal or cpl.getStatus() == IloAlgorithm::Feasible)
				{
					double newTime = timer.elapsed_time(Timer::VIRTUAL);
					double lastVal = double(cpl.getObjValue());  //ORIGINAL
					//double lastVal = double(round(cpl.getObjValue())); //modificado
					double lastGap = 100.0*cpl.getMIPRelativeGap();
					if (lastGap < 0.0) lastGap *= -1.0;

					if (lastVal < results[iter] or lastGap < gaps[iter])
					{
						results[iter] = lastVal;
						times[iter] = newTime;
						gaps[iter] = lastGap;
				    	//cout << "best " << int(lastVal) << "\ttime " << newTime << "\tgap " << lastGap << endl;
					}

					if(int(lastVal)<best_solution && lastVal>0 && verNum(lastVal))
					{
				        //********************************************************************************************
				        //COMENTADO PARA TUNING
				        //********************************************************************************************

						if(!tuning) cout << "best " << int(lastVal) << "\ttime " << newTime << "\tgap " << lastGap << "\tnv "<<lastVal<<endl;
				    	//********************************************************************************************

						best_solution=int(lastVal); 
						best_time=newTime;
						results[iter]=best_solution;

				    	if(best_solution<integrated_best_solution)
				    	{
				    		if(!tuning) cout << "Mejor solucion en CPLEX2: " << best_solution 
				    		<< "| tiempo = " << newTime << "\n";
				    		integrated_best_solution=best_solution;
				    		integrated_best_time=newTime;
				    	}
					}


				    //POSIBILIDAD DE CREACION DE NUEVO CROMOSOMA AQUI EN BASE  SOLUCION DE SOLVER
				    //solution reconstitution

				    //cout<<"ANTES"<<endl;
				    //show_population(population, n_elites, n_mutants, n_offspring);  
				    //show_individual(population[n_of_vertices-1]);

				    set <int> vertex_in_solution;
				    vector < set<int> > vertex_dominated (n_of_vertices);

					//cout<<"A"<<endl;
					if(!ilp_feedback_disable)
					{    
						for(int i=0;i<n_of_vertices;i++)
							if(cpl.getValue(x[i]>0))
							{ 
								vertex_in_solution.insert(i);
					    		// cout<<"vertex: "<<i<<" in solution"<<endl;
					    		//population[population_size-1].vec[i]=1.0;  //opción inicial --mal resultado--
							}
					    	//else population[population_size-1].vec[i]=0.0; //opción inicial --mal resultado--

							for(int i=0;i<n_of_vertices;i++)
							{
								for(int j=0;j<n_of_vertices;j++)
								{
									if(cpl.getValue(y[i][j]>0))
									{
										vertex_dominated[i].insert(j);
										// cout<<i<<" domina a: "<<j<<endl;
									}
								}
							}

							integrate(population[population_size-1],vertex_in_solution);
							evaluate(population[population_size-1]);
							if(verbose) cout<<"Calidad inyectado luego de integración: "<<population[population_size-1].ofv<<endl;
							//evaluate(population[population_size-1]);
					}

					// cout<<"C"<<endl;
					//cout<<"DESPUES"<<endl;
					//show_individual(population[population_size-1]);
					//cout<<"D"<<endl;
					//show_population(population, n_elites, n_mutants, n_offspring);    
					//exit(0);
				}

				//Eliminar restricción 
				//cpl.clearModel(); //EXPERIMEntAL
				//model.remove(new_constraints);

				vertex_in_c_old=vertex_in_c;
	        	
	        	//cout<<"Momento de copia: "<<vertex_in_c.size()<<" "<<vertex_in_c_old.size()<<endl;


				cpl.end();

			}catch(IloException &e){
				cerr << "Cplex Exception " << e.getMessage() << endl;
				env.end();
			}catch(...){
				cerr << "Non Cplex Exception" << endl;
				env.end();
			}

			env.end();

			//cout<<"fin cplex"<<endl;
		}

    } //fin while   

    //cout << "end file " << file_name << endl;

    ++iter;

    //}   //cambio de cargador 


    double r_mean = 0.0;
    double g_mean = 0.0;
    for (int i = 0; i < results.size(); i++){
		//cout<<results[i]<<","<<endl; //--silenciado para tuning
    	r_mean = r_mean + results[i];
    	g_mean = g_mean + times[i];
    }

	r_mean = r_mean / ((double)results.size());
	g_mean = g_mean / ((double)times.size());
	double rsd = 0.0;
	double gsd = 0.0;

	for (int i = 0; i < results.size(); i++) {
		rsd = rsd + pow(results[i]-r_mean,2.0);
		gsd = gsd + pow(times[i]-g_mean,2.0);
	}

	rsd = rsd / ((double)(results.size()-1.0));
	if (rsd > 0.0) {
		rsd = sqrt(rsd);
	}

	gsd = gsd / ((double)(times.size()-1.0));
	if (gsd > 0.0) {
		gsd = sqrt(gsd);
	}

    //cout << r_mean << "\t" << rsd << "\t" << g_mean << "\t" << gsd<< endl; //--silenciado para tuning
    //cout<<"ibest: "<<integrated_best_solution<<" itime: "<<integrated_best_time<<endl;
    //cout<<"ibest: "<<integrated_best_solution<<" itime: "<<integrated_best_time<<endl;

	double composed_result;
	double timeComponent = integrated_best_time/1000;
	if (timeComponent >1) timeComponent = 0.99;
	composed_result=integrated_best_solution+timeComponent;
	if(tuning) cout<<integrated_best_solution<<endl;
	if(!tuning) cout<<"ibest:\t"<<integrated_best_solution<<"\titime:\t"<<integrated_best_time<<"\tComposed:\t"<<composed_result<<endl;
}