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
#include <numeric>
#include <climits>
#include <limits.h>

#include "createHeuristicSolution.h"

using namespace std;

// Data structures for the problem data
vector<string> input_sequence;
int n_of_sequences;
int sequence_length;
int alphabet_size = 4;
map<int,char> mapping;
map<char,int> rev_mapping;
double threshold = 0.70;
int t_value;

double t_limit;

int flag = 0;

//Genetic Algorithm Hyper-Parameters
int n_population = 50; //Cambiar por un limite de tiempo
int population_Creation = 0;
int p_pop_heuristic = 0;
int n_tournement = 1;
double p_tournement_contestans = 0.1;
int k_contestants = 5; 
int crossover_type = 0; // 0 uniform; 1 single point crossover.

int r_population = 0;   // 0 new + random_previous; 1 new + elite; 2 new + elite + random_previous;
double p_newGen = 0.7;       // parameter of percentage of new offspring permited for r_population 1 & 2
int nng_value;
double p_elite = 0.2;    // parameter of percentage of elites
int elite_value;
int rndm_replace_type = 0; //0 select random from this gen; 1 create completely random new individual
double mutation_rate = 0.05; //parameter of probability to change a individual

float determinism = 1;

// Genetic Algorithm Data Structures
vector<pair<int,pair<string,int>>> population;
vector<pair<int,pair<string,int>>> generation_new;
vector<pair<int,pair<string,int>>> gen_elites;


// vector for keeping all the names of the input files
vector<string> inputFiles;

vector<string> testingPop;

// Annealing variables
int nAscending = 0;
int nDescending = 0;
int steady = 0;
int nToChange = 3; // how many steadys are needed to change the treshold;
double steps = 1;
int prevVal = INT_MIN;
enum Direction { Ascending, Descending}; 

Direction direction = Ascending;

vector<tuple<double, int, double, string>> bests_per_threshold;
tuple <double, int, double, string> best_so_far = make_tuple(0.0,0,0.0, "");



// dummy parameter as an example for creating command line parameters -> see function read_parameters(...)
int dummy_integer_parameter = 0;


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
        else if (strcmp(argv[iarg],"-np")==0) n_population = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pc")==0) population_Creation = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pph")==0) p_pop_heuristic = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-dt")==0) determinism = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-ntc")==0) nToChange = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-stp")==0) steps = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-ct")==0) crossover_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-mr")==0) mutation_rate = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-png")==0) p_newGen = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pe")==0) p_elite = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-rp")==0) r_population = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-rrt")==0) rndm_replace_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-flag")==0) flag = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-param1")==0) dummy_integer_parameter = atoi(argv[++iarg]); // example for creating a command line parameter param1 -> integer value is stored in dummy_integer_parameter
        iarg++;
    }
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
    auto rng = default_random_engine {};        
    switch (r_population)
    {
    case 0:
        shuffle(population.begin(), population.end(), rng);
        for (int aa = 0; generation_new.size() < n_population; aa++)
            generation_new.push_back(population[aa]);
        break;
    case 1:
        sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
        for (int l = 0; generation_new.size() < n_population; l++) 
            generation_new.push_back(population[l]);
        break;
    case 2:
        sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
        elite_value = (n_population - generation_new.size()) * p_elite;
        int rndm_pop_value = n_population - generation_new.size() - elite_value;
        for(int l = 0; l < elite_value; l++){
            gen_elites.push_back(population[0]);
            population.erase(population.begin());
        }
        if (rndm_replace_type == 0){
            shuffle(population.begin(), population.end(), rng);
            for(int aa = 0; aa < rndm_pop_value; aa++)
                generation_new.push_back(population[aa]);
        }                
        else if (rndm_replace_type == 1){
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
            if(rand() % 100 < mutation_rate * 100){
                char ch = mapping[rand() % 4];
                while (ch == population[i].second.first.at(j)) ch = mapping[rand() % 4];
                population[i].second.first.at(j) = ch;
                mutated++;
            }
        }
    }
    if (flag > 1) cout << mutated << " individuals mutated" << endl;
    population.insert(population.end(), gen_elites.begin(), gen_elites.end());
    gen_elites.clear();
}


void annealing(int best){
    
    bool change = false;
    t_value = int(threshold * sequence_length);

    int curVal = best;

    if (prevVal < curVal) {  // (still) ascending?
        if (direction != Ascending) {
            direction = Ascending;
            nAscending++;
        }
        else{ nAscending++; steady++;}
        direction = Ascending;
    }
    else if (prevVal > curVal) { // (still) descending?
        if (direction != Descending) { // starts descending?
            direction = Descending;
            nDescending++;
        }
        else {nDescending++; steady++;}
        direction = Descending;
    }
    if (curVal == prevVal) steady++;
    if ( steady + nToChange < nAscending + nDescending  ||  nAscending + nDescending + nToChange < steady) change = true;

    //if(flag)cout << nAscending <<endl << nDescending << endl<<steady <<endl <<endl;


    prevVal = curVal;

    if (change){
        threshold+= (15.0/steps)/100.0;
        steady = 0;
        nAscending = 0;
        nDescending = 0;
        t_value = int(threshold * sequence_length);
        best_so_far = make_tuple(0.0,-1,0.0, "");
        if(flag)cout << threshold <<endl;
    }

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
            cout << "Error: file could not be opened" << endl;
        }

        input_sequence.clear();
        string seq;
        while (indata >> seq) input_sequence.push_back(seq);
        n_of_sequences = input_sequence.size();
        sequence_length = input_sequence[0].size();
        indata.close();

        if(flag) cout << "n of sequences " << n_of_sequences << endl;
        if(flag) cout << "sequence lenght " << sequence_length << endl;

        // minimum required Hamming distance
        t_value = int(threshold * sequence_length);

        if(flag) cout << threshold << " threshold " << sequence_length << " sequence lenght " << t_value << " t_value" << endl;
        // the computation time starts now
        //clock_t start = clock();

        // Example for requesting the elapsed computation time (in seconds) at any moment: 
        // clock_t end = clock();
        // double elapsed = double(end - start)/CLOCKS_PER_SEC;

        if(flag) cout << "start file " << inputFiles[na] << endl;    

        // generar dt = t * porcentago (0.50) * iteraciones (generaciones)
        srand(time(0));
        
        //Population Creation
        int heuristic_population  = 0;
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
        case 2:
            heuristic_population = n_population * p_pop_heuristic; 
            for (int i = 0; i < heuristic_population; i++)
                population.push_back(make_pair(i,make_pair(createHeuristicSolution(inputFiles[na], threshold, determinism),0)));
            for (int i = 0; i < n_population - heuristic_population; i++)
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
        default:
            break;
        }
        
        nng_value = population.size() * p_newGen;
        clock_t start, end; 
        start = clock();
        int terminate_algorithm = 1, terminate_tournements = 1;
        end = clock();
        double elapsed = double(end - start)/CLOCKS_PER_SEC;                
        
        double endPoint = 0.85;

        while (terminate_algorithm){
            //Population Finesse calculation
            fitness_calculation();

            if(flag > 1)for (int i = 0; i < population.size(); i++) cout << "[" << population[i].second.first <<"] " << population[i].second.second << endl;
            sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
            //if(flag) cout << population[0].second.second << " best    ";//<< endl; 
            

            if (population[0].second.second > get<1>(best_so_far)){  
                end = clock();
                elapsed = double(end - start)/CLOCKS_PER_SEC;
                //cout << threshold << '\t' << population[0].second.second << '\t' << elapsed << endl;
                best_so_far = make_tuple(threshold, population[0].second.second, elapsed, population[0].second.first);
                bests_per_threshold.push_back(best_so_far);
            }

            if (threshold < 0.85)
                annealing(population[0].second.second);
            
            if (flag > 1) for (int a = 0; a < population.size(); a++){
                cout <<population[a].second.second << endl;
            }

            //Tournament ARC
            //select k candidates//
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

            if (flag > 2)for (int ii = 0; ii < population.size(); ii++) cout << generation_new[ii].second.first<<endl<<endl;

            replace_population();
            mutation();
            //mutacion sobre nP menos elite
            end = clock();
            elapsed = double(end - start)/CLOCKS_PER_SEC;
            //if(flag)cout << elapsed << " elapsed"<<endl;
            if (elapsed >= t_limit)
                terminate_algorithm = 0;

        }

        // HERE GOES YOUR GREEDY HEURISTIC
        // When finished with generating a solution, first take the computation time as explained above. Say you store it in variable 'time'.
        // Then write the following to the screen: cout << "value " << <value of your solution> << "\ttime " << time << endl;

        if(flag) cout << "end file " << inputFiles[na] << endl;
    }

   // if (flag){ 
    for (int k = 0; k < bests_per_threshold.size(); k++)
        cout << /* get<0>(bests_per_threshold[k]) << '\t' << */ get<1>(bests_per_threshold[k]) << '\t' << get<2>(bests_per_threshold[k]) << '\t' << endl;// << get<3>(bests_per_threshold[k]) << endl;
    //}
    //cout << get<1>(best_so_far) * -1 <<endl;
    
    
    // calculating the average of the results and computation times and write them to the screen
    double r_mean = 0.0;
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    if(flag) cout << r_mean << "\t" << t_mean << endl;
    
}

