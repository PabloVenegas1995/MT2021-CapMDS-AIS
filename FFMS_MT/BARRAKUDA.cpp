/***************************************************************************
                          greedy.cpp  -  description
                             -------------------
    begin                : Fri Dec 17 2021
    copyright            : (C) 2021 by Christian Blum
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Timer.h"
#include <array>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <limits>
#include <iomanip>
#include <chrono>

#include <random>
#include <utility>
#include <iterator>
#include <algorithm>
#include <random>

#include "createHeuristicSolution.h"

#include <ilcplex/ilocplex.h>

using namespace std;

ILOSTLBEGIN

// Data structures for the problem data
vector<string> input_sequence;
int n_of_sequences;
int sequence_length;
int alphabet_size = 4;
map<int,char> mapping;
map<char,int> rev_mapping;
double threshold;
int t_value;

double t_limit = 90.0;

int flag = 0;

//Genetic Algorithm Hyper-Parameters
int n_population = 50; //Cambiar por un limite de tiempo
int population_Creation = 0; // 0 creates the population randomly, 1 creates the population with Prof. Pinacho's Heuristics
int parents_selection_type = 0; // 0 tournement selection, 1 Fitness proportionate selection
int k_contestants = 5; 
int crossover_type = 0; // 0 uniform; 1 single point crossover; 2 double point crossover.

int r_population = 0;   // 0 new + elite; 1 new + elite + random_previous;
double p_newGen = 0.7;       // parameter of percentage of new offspring permited for r_population 0 & 1
int nng_value;
double p_elite     = 0.2;    // parameter of percentage of elites
int elite_value;
double mutation_rate = 0.05; //parameter of probability to change a individual

float determinism = 1;

//Barrakuda specific 
int barrakuda_strategy = 0; // 0 best + n_elites; 1 best + random
int barrakuda_n = 0; // number of elites or randoms to add to S' barrakuda_n < n_population 1-100
vector<pair<int,pair<string,int>>> barrakuda_subinstance;
string barrakuda_solution;

// Genetic Algorithm Data Structures
vector<pair<int,pair<string,int>>> population;
vector<pair<int,pair<string,int>>> generation_new;
vector<pair<int,pair<string,int>>> gen_elites;

vector<pair<int,pair<string,int>>> parents_FPS;

tuple <int, double, string> best_so_far = make_tuple(-1,0.0, "");



// vector for keeping all the names of the input files
vector<string> inputFiles;

double cplex_time_limit = 5.0;


// inline int stoi(string &s) {

//   return atoi(s.c_str());
// }

// inline double stof(string &s) {

//   return atof(s.c_str());
// }

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFiles.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tlim")==0) t_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-np")==0) n_population = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pc")==0) population_Creation = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-dt")==0) determinism = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pst")==0) parents_selection_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-ct")==0) crossover_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-png")==0) p_newGen = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pe")==0) p_elite = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-rp")==0) r_population = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-mr")==0) mutation_rate = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-cpl_t")==0) cplex_time_limit = atoi(argv[++iarg]);// reading t he computation time limit of CPLEX for the application in the hybrid metaheuristic from the command line (if provided)
        else if (strcmp(argv[iarg],"-bks")==0) barrakuda_strategy = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-bkn")==0) barrakuda_n = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-flag")==0) flag = atoi(argv[++iarg]);
        iarg++;
    }
}

//  CPLEX Callback for logging the solution

ILOSOLVECALLBACK4(loggingCallback,
    Timer&, timer,
    double&, time_stamp,
    double&, result,
    double&, gap) {

    if (hasIncumbent()) {
        IloNum nv = getIncumbentObjValue();
        //clock_t current = clock();
        //double newTime = double(current - start) / CLOCKS_PER_SEC;
        double newTime = timer.elapsed_time(Timer::VIRTUAL);
        double newGap = 100.0 * getMIPRelativeGap();
        if (result < double(nv)) {
            if(flag) cout << "value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
            result = double(nv);
            time_stamp = newTime;
            gap = newGap;
        }
    }
}

void show_subinstance(vector<array<int,4>> &subinstance) {

    for (int i = 0; i < subinstance.size(); i++) {
        for (int j = 0; j < 4; j++) {
            cout << subinstance[i][j] << " ";
        }
        cout << endl;
    }
}

void run_cplex(vector<array<int,4>> &subinstance, Timer & timer, vector<double>& results, vector<double>& times, vector<double>& gaps, int& na) {


    IloEnv env;
    // this instruction redirects the standard output of CPLEX to the null stream in order to avoid printing it to the screen
    env.setOut(env.getNullStream());
    try{

        IloModel model(env);
        //cout<<"Definiendo variables"<<endl;

        // Here you have to implement the ILP model for its resolution with CPLEX
        //***********************
        //      Variables:
        //***********************
        //Variable 2D para representar letra en posicion del string 'x'
        typedef IloArray<IloNumVarArray> IntVar2DMatrix;
        //typedef IloArray<IntVar2DMatrix> IntVar3DMatrix;
        //typedef IloArray<IntVar3DMatrix> IntVar4DMatrix;
        //typedef IloArray<IntVar4DMatrix> IntVar5DMatrix;
        
        //x[posicion][BASE-GEN]
        IntVar2DMatrix x(env,sequence_length);
        for(int i=0; i<sequence_length; i++){
            x[i]=IloNumVarArray(env,alphabet_size);
            for(int j=0; j<alphabet_size;j++)
                x[i][j]=IloNumVar(env,0,1, ILOINT);
        }

        //cout<<"x ya definido"<<endl;
        //Variable binaria y que identifica que secuencia está en omega.
        IloNumVarArray y(env,n_of_sequences,0,1, ILOINT);
        //cout<<"y ya definido"<<endl;
        //Minimizando la función objetivo
        IloExpr obj(env);
        for(int i=0; i<n_of_sequences;i++) obj += y[i];
        model.add(IloMaximize(env,obj));
        obj.end();
        //cout<<"Variables definidas"<<endl;
        //********************************
        //        Restricciones:
        //********************************
        //cout<<"ingresando Const. 1"<<endl;
        //Restriccion 1: construccion correcta de x
        for(int i=0; i<sequence_length;i++){
            IloExpr expr(env);
            for (int j=0; j<alphabet_size;j++)
                expr +=x[i][j];
            model.add(expr == 1);
            expr.end();
        }
        //cout<<"saliendo Const. 1"<<endl;
        //cout<<"ingresando Const. 2"<<endl;

        //restriccion 2: condiciones para que yr sea 1
        for(int i=0; i<n_of_sequences;i++){
            IloExpr expr(env);
            for (int j=0; j< sequence_length;j++){
                expr +=x[j][rev_mapping[input_sequence[i][j]]];
            }
            model.add(expr <= sequence_length - t_value * y[i] );
            expr.end();


        }
         //cout<<"saliendo Const. 2"<<endl;

        //Add constraints of subinstance
        ////////////////////////////////////         CHANGE FROM HEREEE ON OUT              //////////////////////////////////////////////
        for(int i=0; i<sequence_length;i++){
            for(int j=0; j<alphabet_size;j++){
                if(subinstance[i][j]==0){
                    IloExpr expr(env);
                    expr +=x[i][j];
                    model.add(expr == 0);
                    expr.end();
                }
            }
        }

       IloCplex cpl(model);

        cpl.setParam(IloCplex::TiLim, cplex_time_limit);
        cpl.setParam(IloCplex::EpGap, 0.0);
        cpl.setParam(IloCplex::EpAGap, 0.0);
        cpl.setParam(IloCplex::Threads, 1);
        cpl.setWarning(env.getNullStream());
        cpl.use(loggingCallback(env, timer, times[na], results[na], gaps[na]));

        cpl.solve();
    
        if (cpl.getStatus() == IloAlgorithm::Optimal || cpl.getStatus() == IloAlgorithm::Feasible) {
            //if (cpl.getStatus() == IloAlgorithm::Optimal) cout << "optimality proven" << endl;
            //clock_t current = clock();
            //double newTime = double(current - start) / CLOCKS_PER_SEC;
            double newTime = timer.elapsed_time(Timer::VIRTUAL);
            double lastVal = double(cpl.getObjValue());
            double lastGap = 100.0 * cpl.getMIPRelativeGap();
            if (lastGap < 0.0) lastGap *= -1.0;
            if (lastVal > results[na]) {
                if(flag) cout << "value " << lastVal << "\ttime " << newTime << "\tgap " << lastGap << endl;
                results[na] = lastVal;
                times[na] = newTime;
                gaps[na] = lastGap;
            }

            //cout<<"Adapt: zeroing..."<<endl;
            for(int i=0;i<sequence_length;i++){
                for(int j=0;j<alphabet_size;j++){
                    IloNum xval =cpl.getValue(x[i][j]);
                    if(xval>0){ 
                        //cout<<mapping[j];  //solución aquí!
                        barrakuda_solution.push_back(mapping[j]);
                        break;
                    }
                }
            }
            for(int i=0;i<sequence_length;i++)
                for(int j=0;j<alphabet_size;j++)
                    subinstance[i][j]=0;

            
        }
        env.end();
        //Adapt: removing older elements in subinstance
    
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
    }
    env.end();
}


//Fitness Calculation
// Hamming Distance
int hammingDist(string str1, string str2)
{
    int count = 0;
    for( int it = 0; it < sequence_length ; it++){
        if (str1.at(it) != str2.at(it)) count++;
    }
    return count;
}

void fitness_calculation(){

    //cout << t_value << " tvalue " << endl;
    for (int i = 0; i < population.size(); i++)
    {
        int finesse = 0;
        int missmatchs = 0;
        for (int j = 0; j < n_of_sequences; j++)
        {
                missmatchs = hammingDist(population[i].second.first, input_sequence[j]);
                if (flag > 2) cout << "missmatch for; [" << i << "] " << missmatchs << endl;
                if (missmatchs >= t_value) finesse += 1;
        }
        population[i].second.second = finesse;
    }
    if (flag > 1) cout << population[0].second.first << endl;        
}

int fitness_single_debug(string str){
    int finesse = 0;
    int missmatchs = 0;
    for (int j = 0; j < n_of_sequences; j++)
    {
            missmatchs = hammingDist(str, input_sequence[j]);
            if (missmatchs >= t_value) finesse += 1;
    }
    return finesse;
}

int fitness_sum(vector<pair<int,pair<string,int>>> parents){
    int sum = 0;
    for (int i; i < parents.size(); i++)
        sum += parents[i].second.second;
    return sum;
}

void parents_selection(vector<pair<int,pair<string,int>>> parents){
    //Shuffle
    double offset = 0.0;
    int sum = fitness_sum(parents);
    //if(flag) cout << sum <<endl;
    double r = ((double) rand() / (RAND_MAX));
    if(sum < 1){
        int i = rand() % parents.size();
        parents_FPS.push_back(parents[i]);
        parents.erase(parents.begin() + i);
    }
    else{
        auto rng = default_random_engine {};
        shuffle(parents.begin(), parents.end(), rng);
        for(int i = 0; i < parents.size(); i++){
            //if(flag) cout << r << " adas " << offset<<endl;
            //if(flag) cout << ((double)parents[i].second.second/ (double) sum) << " adas " <<endl;    
            offset += ((double)parents[i].second.second/ (double) sum);
            //if(flag) cout << r << " adas " << offset<<endl;
            if(r < offset){
                parents_FPS.push_back(parents[i]);
                //if(flag >0) cout << parents_FPS.size() <<endl;
                parents.erase(parents.begin() + i);
                break;
            } 
        }
    }
    //if(flag) cout << "First selected" <<endl;
    offset = 0.0;
    sum = fitness_sum(parents);
    r = ((double) rand() / (RAND_MAX));
    if(sum < 1){
        int i = rand() % parents.size();
        parents_FPS.push_back(parents[i]);
        parents.erase(parents.begin() + i);
    }
    else{
        for(int i = 0; i < parents.size(); i++){    
            offset += ((double)parents[i].second.second/ (double) sum);
            if(r < offset){
                parents_FPS.push_back(parents[i]);
                parents.erase(parents.begin() + i);
                break;
            } 
        }
    }
    //if(flag) cout << "second selected" <<endl;
}


void mating(pair<int,pair<string,int>> parent1, pair<int,pair<string,int>> parent2){

    string child1 = "", child2 = "";        
    if (crossover_type == 0){ //Uniform Crossover
        for (int i = 0 ; i < sequence_length; i++){
            if (rand() % 2 == 0){    
                child1 += parent1.second.first.at(i);
                child2 += parent2.second.first.at(i);
            }
            else{
                child1 += parent2.second.first.at(i);
                child2 += parent1.second.first.at(i);
            }
        }
    }
    else if (crossover_type == 1){ // Crossover Point, pivote al azar
        int crossover_point = rand () % (sequence_length - 1);
        for (int i = 0; i < sequence_length; i++)
        {
            if ( i >= crossover_point){
                child1.push_back(parent2.second.first.at(i));
                child2.push_back(parent1.second.first.at(i));    
            }
            else{
                child1.push_back(parent1.second.first.at(i));
                child2.push_back(parent2.second.first.at(i));
            }
        }
    }
    else if (crossover_type == 2){ //Dual Crossover point, pivotes al azar
        int crossover_point1 = rand () % (sequence_length -1), crossover_point2 = rand () % (sequence_length -1);
        while (crossover_point1 == crossover_point2){
            crossover_point2 = rand () % (sequence_length -1);}
        if (crossover_point1 > crossover_point2) swap(crossover_point1,crossover_point2);
        for (int i = 0; i < sequence_length; i++){
            if (i < crossover_point1 || i > crossover_point2){
                child1.push_back(parent1.second.first.at(i));
                child2.push_back(parent2.second.first.at(i));
            }
            else{
                child1.push_back(parent2.second.first.at(i));
                child2.push_back(parent1.second.first.at(i));
            }
        }
    }
    
    generation_new.push_back(make_pair(parent1.first + 100,make_pair(child1,0)));
    generation_new.push_back(make_pair(parent2.first + 100,make_pair(child2,0)));
}

bool sortbysec( const pair<int,pair<string,int>> &a, const pair<int,pair<string,int>> &b)
{
    return (a.second.second > b.second.second);
}

void replace_population(){
    sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
    switch (r_population)
    {
    case 0:
        for (int l = 0; generation_new.size() < n_population; l++) 
            generation_new.push_back(population[l]);
        break;
    case 1:
        elite_value = (n_population - generation_new.size()) * p_elite;
        int rndm_pop_value = n_population - generation_new.size() - elite_value;
        for(int l = 0; l < elite_value; l++){
            gen_elites.push_back(population[0]);
            population.erase(population.begin());
        }
        for (int i = 0; i < rndm_pop_value; i++)
        {
            string str = "";
            for (int j = 0; j < sequence_length; j++)
            {
                char ch = mapping[rand () % 4];
                str.push_back(ch);
            }
            generation_new.push_back(make_pair(i + 100,make_pair(str,0)));
        }
        break;
    }    
    population.clear();
    population.insert(population.begin(), generation_new.begin(), generation_new.end());
    generation_new.clear();        
}

void mutation(){
    int mutated = 0;
    for (int i = 0; i < population.size(); i++){
        for (int j = 0; j < sequence_length; j++){
            if(rand() % 100 <= mutation_rate * 100){
                char ch = mapping[rand() % 4];
                while (ch == population[i].second.first.at(j)) ch = mapping[rand() % 4];
                population[i].second.first.at(j) = ch;
                mutated++;
            }
        }
    }
    if (flag > 1) cout << mutated << " individals mutated" << endl;
    population.insert(population.end(), gen_elites.begin(), gen_elites.end());
    gen_elites.clear();
}

/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // number of input files
    int n_files = int(inputFiles.size());

    // vectors for storing the result and the computation time obtained by the applications of the greedy heuristic
    vector<double> results(n_files, std::numeric_limits<int>::min());
    vector<double> times(n_files, 0.0);
    vector<double> gaps(n_files, 0.0);

    // letter to index mapping (note: this only works for instances on alphabet Sigma= {A, C, T, G}
    mapping[0] = 'A';
    mapping[1] = 'C';
    mapping[2] = 'T';
    mapping[3] = 'G';
    rev_mapping['A'] = 0;
    rev_mapping['C'] = 1;
    rev_mapping['T'] = 2;
    rev_mapping['G'] = 3;

    // main loop over all input files (problem instances)
    for (int na = 0; na < n_files; ++na) {

        // opening the corresponding input file and reading the problem data
        ifstream indata;
        indata.open(inputFiles[na].c_str());
        if(!indata) { // file couldn't be opened
            if (flag) cout << "Error: file could not be opened" << endl;
        }

        input_sequence.clear();
        string seq;
        while (indata >> seq) input_sequence.push_back(seq);
        n_of_sequences = input_sequence.size();
        sequence_length = input_sequence[0].size();
        indata.close();

        if (flag) cout << "n of sequences " << n_of_sequences << endl;
        if (flag) cout << "sequence lenght " << sequence_length << endl;

        // minimum required Hamming distance
        t_value = int(threshold * sequence_length);

        //Data Structures for BARRAKUDA
        vector < array<int,4> >  subinstance;
        for (int i=0; i<sequence_length; ++i) {
            array <int, 4> aux = {0,0,0,0};
            subinstance.push_back(aux);
        }


        if(flag) cout << threshold << " threshold " << sequence_length << " sequence lenght " << t_value << " t_value" << endl;

        // Example for requesting the elapsed computation time (in seconds) at any moment: 
        // clock_t end = clock();
        // double elapsed = double(end - start)/CLOCKS_PER_SEC;

        if (flag) cout << "start file " << inputFiles[na] << endl;

        // generar dt = t * porcentago (0.50) * iteraciones (generaciones)
        srand(time(0));

        Timer timer;

        
        //Population Creation
        switch (population_Creation)
        {
        case 0:
            for (int i = 0; i < n_population; i++)
            {
                string str = "";
                for (int j = 0; j < sequence_length; j++)
                {
                    char ch = mapping[rand() % 4];
                    str.push_back(ch);
                }
                //cout << str <<endl<<endl; 
                population.push_back(make_pair(i,make_pair(str,0)));
            }
            break;
        case 1:
            for (int i = 0; i < n_population; i++)
                population.push_back(make_pair(i,make_pair(createHeuristicSolution(inputFiles[na], threshold, determinism),0)));
            break;
        default:
            break;
        }

        nng_value = population.size() * p_newGen;
        clock_t start, end; 
        start = clock();
        int terminate_tournements = 1;
        int barrakuda_gen = 0;
        int pst_original = parents_selection_type;
        while (timer.elapsed_time(Timer::VIRTUAL)<t_limit){
            //Population Finesse calculation
            fitness_calculation();

            if(flag > 1)for (int i = 0; i < population.size(); i++) cout << "[" << population[i].second.first <<"] " << population[i].second.second << endl;
            sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
            if(flag) cout << population[0].second.second << " best    ";//<< endl; 

            if (population[0].second.second > get<0>(best_so_far)){  
                end = clock();
                double elapsed = double(end - start)/CLOCKS_PER_SEC;
                best_so_far = make_tuple(population[0].second.second, elapsed, population[0].second.first);
                //if(flag)  
                    cout << get<0>(best_so_far) << '\t' << get<1>(best_so_far) <<endl;
            }
            if (flag > 1) 
                for (int a = 0; a < population.size(); a++)
                    cout <<population[a].second.second << endl;

            //Tournament ARC
            //select k candidates
            //vector<pair<int,pair<string,int>>> group3 = population;
            // if (pst_original == 1 ){
            //     //cout << "entre if";
            //     if (fitness_sum(population) < 1)
            //         parents_selection_type = 0;
            //     else parents_selection_type = 1;
            // }
            switch (parents_selection_type)
            {
            case 0:
                terminate_tournements = 1;
                while (terminate_tournements){            
                    vector<pair<int,pair<string,int>>> group1, group2;
                    auto rng = default_random_engine {};
                    shuffle(population.begin(), population.end(), rng);
                    
                    for (int l = 0; l < k_contestants * 2; l++){
                        if (l <k_contestants) group1.push_back(population[l]);
                        else group2.push_back(population[l]);
                    } 
                    
                    //pick 2 best
                    pair<int,pair<string,int>> better_parent = group1[0], best_parent = group2[0];
                    for (int j = 0; j < group1.size(); j++)
                    {
                        if(group1[j].second.second > better_parent.second.second)
                            better_parent = group1[j];
                        if(group2[j].second.second > best_parent.second.second)
                            best_parent = group2[j];
                    }
                    
                    //Make them procreate
                    mating(better_parent, best_parent);
                    if (generation_new.size() >= nng_value){
                        while(generation_new.size() > nng_value) generation_new.pop_back();
                        terminate_tournements = 0;
                    } 
                }
                break;
            case 1:
                terminate_tournements = 1;
                while (terminate_tournements){
                    vector<pair<int,pair<string,int>>> group3 = population;
                    parents_selection(group3);
                    mating(parents_FPS[0], parents_FPS[1]);
                    //if(flag) cout << "Parents selection time" <<endl;
                    parents_FPS.clear();
                    if (generation_new.size() >= nng_value){
                            while(generation_new.size() > nng_value) generation_new.pop_back();
                            terminate_tournements = 0;
                        }
                }
                break;
            
            default:
                break;
            }

            if (flag > 2)for (int ii = 0; ii < population.size(); ii++) cout << generation_new[ii].second.first<<endl<<endl;
            replace_population();
            mutation();

            fitness_calculation();
            
            // BARRAKUDA GOES HERE
            vector<pair<int, pair<string, int>>> barrakuda_population = population;
            sort(barrakuda_population.begin(), barrakuda_population.end(), sortbysec);
            barrakuda_subinstance.push_back(barrakuda_population[0]);
                        
            barrakuda_population.erase(barrakuda_population.begin());
            switch (barrakuda_strategy)
            {
            case 0:
                for(int i = 0; barrakuda_subinstance.size() < barrakuda_n; i++)
                    barrakuda_subinstance.push_back(barrakuda_population[i]);
                break;
            case 1:
                {
                    auto rng = default_random_engine {};
                    shuffle(barrakuda_population.begin(), barrakuda_population.end(), rng);
                    for (int i = 0; i < barrakuda_n; i++)
                        barrakuda_subinstance.push_back(barrakuda_population[i]);                
                    break;
                }
            default:
                break;
            }
            barrakuda_population.clear();
            // Pass barrakuda subinstance to solver and get solution
            for(int m =0; m < barrakuda_subinstance.size(); m++){
                for(int idx=0;idx<sequence_length;idx++){
                    if(subinstance[idx][rev_mapping[barrakuda_subinstance[m].second.first.at(idx)]]==0){
                        //case: first time seeing this character
                        subinstance[idx][rev_mapping[barrakuda_subinstance[m].second.first.at(idx)]]++;
                    }
                    else{
                        subinstance[idx][rev_mapping[barrakuda_subinstance[m].second.first.at(idx)]]++;
                        //case: already seen this character
                        //nothing to do.
                    }
                }
            }
            barrakuda_subinstance.clear();

            //show_subinstance(subinstance);

            run_cplex(subinstance, timer, results, times, gaps, na);


//            cout << barrakuda_solution << endl;
            
            sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
            population.pop_back();
            population.push_back(make_pair(1000,make_pair(barrakuda_solution, 0)));

//            cout << fitness_single_debug(barrakuda_solution) << endl;

            barrakuda_solution.clear();

            fitness_calculation();

            //mutacion sobre nP menos elite
            end = clock();
            double elapsed = double(end - start)/CLOCKS_PER_SEC;
            if(flag)cout << elapsed << " elapsed"<<endl;
            
        }
        if (flag) cout << "end file " << inputFiles[na] << endl;
    }

    //cout << "best so far is " << get<0>(best_so_far) << " at " << get<1>(best_so_far) << " seconds" <<endl;
    //if(flag) 
        cout << get<0>(best_so_far) << '\t' << get<1>(best_so_far) <<endl;
    // calculating the average of the results and computation times and write them to the screen
    double r_mean = 0.0;
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    //if (flag) cout << r_mean << "\t" << t_mean << endl; else cout << get<0>(best_so_far) * -1 <<endl;
    
}