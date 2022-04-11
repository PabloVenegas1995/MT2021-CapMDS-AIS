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

using namespace std;

// Data structures for the problem data
vector<string> input_sequence;
int n_of_sequences;
int sequence_length;
int alphabet_size = 4;
map<int,char> mapping;
map<char,int> rev_mapping;
double threshold;
int t_value;

int flag = 0;

//Genetic Algorithm Hyper-Parameters
int n_population = 50;
int n_generation = 1;
int n_tournement = 1;
int k_tournement_contestant = n_population / 10;
 
int crossover_type = 0; // 0 uniform; 1 single point crossover.
int crossover_point = 0;

double mutation_rate = 0.05;


// Genetic Algorithm Data Structures
vector<string> population;
vector<pair<string, int>> population_finesse;
vector<string> generation_new;



// vector for keeping all the names of the input files
vector<string> inputFiles;

// dummy parameter as an example for creating command line parameters -> see function read_parameters(...)
int dummy_integer_parameter = 0;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFiles.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-np")==0) n_population = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-ng")==0) n_generation = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-nt")==0) n_tournement = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-ct")==0) crossover_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-cp")==0) crossover_point = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-flag")==0) flag = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-param1")==0) dummy_integer_parameter = atoi(argv[++iarg]); // example for creating a command line parameter param1 -> integer value is stored in dummy_integer_parameter
        iarg++;
    }
}

//Fitness Calculation

// Hamming Distance
int hammingDist(string str1, string str2)
{
    int i = 0, count = 0;
    for ( std::string::iterator it=str1.begin(); it!=str1.end(); ++it) 
    {
        if (str1.at(i) != str2.at(i))
            count++;
        i++;
    }
    return count;
}

void fitness_calculation(vector<string> population_tested, vector<string> initial_input){

    for (int i = 0; i < population_tested.size(); i++)
    {
        int finesse = 0;
        int missmatchs = 0;
        for (int j = 0; j < n_of_sequences; j++)
        {
                missmatchs = hammingDist(population_tested[i], initial_input[j]);
                //cout << "missmatchs found: "<< missmatchs << endl;
                if (missmatchs >= t_value) finesse += 1;
        }
        //if(flag) cout << "finesse for [" << i << "]: " << finesse << endl;
        population_finesse.push_back( make_pair(population_tested[i], finesse) );        
    }
    
}

vector<string> mating(string parent1, string parent2){

    vector<string> offspring;
    string child1, child2;
    if (crossover_type == 0){ //Uniform Crossover
        for (int i = 0 ; i < sequence_length; i++){
            if ( rand() % 2 == 0){
                child1.push_back(parent1.at(i));
                child2.push_back(parent2.at(i));
            }
            else{
                child1.push_back(parent2.at(i));
                child2.push_back(parent1.at(i));
            }

        }
    }
    else if (crossover_type == 1){
        if (crossover_point == 0) crossover_point = sequence_length / 2;
        for (int i = 0; i < sequence_length; i++)
        {
            if ( i >= crossover_point){
                child1.push_back(parent2.at(i));
                child2.push_back(parent1.at(i));    
            }
            else{
                child1.push_back(parent1.at(i));
                child2.push_back(parent2.at(i));
            }
        }
    }
    offspring.push_back(child1);
    offspring.push_back(child2);
    return offspring;
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

        cout << "n of sequences " << n_of_sequences << endl;
        cout << "sequence lenght " << sequence_length << endl;

        // minimum required Hamming distance
        t_value = int(threshold * sequence_length);

        // the computation time starts now
        clock_t start = clock();

        // Example for requesting the elapsed computation time (in seconds) at any moment: 
        // clock_t end = clock();
        // double elapsed = double(end - start)/CLOCKS_PER_SEC;

        cout << "start file " << inputFiles[na] << endl;

        random_device dev;
        mt19937 rng(dev());
        uniform_int_distribution<mt19937::result_type> dist6(0,3); // distribution in range [1, 6]
        
        //Population Creation
        for (int i = 0; i < n_population; i++)
        {
            string str = "";
            for (int j = 0; j < sequence_length; j++)
            {
                char ch = mapping[dist6(rng)];
                str.push_back(ch);
            }
            //cout << str << endl;
            population.push_back(str);
        }
        

        for (int generation = 0; generation < n_generation; generation++)
        {
            //Population Finesse calculation
            fitness_calculation(population, input_sequence);

            //Tournament ARC
            //select k candidates

            for (int i = 0; i < n_tournement; i++)
            {    
                shuffle(population_finesse.begin(), population_finesse.end(), mt19937{std::random_device{}()});
                vector<pair<string,int>> contestant_group1(population_finesse.begin(), population_finesse.begin() + k_tournement_contestant);
                vector<pair<string,int>> contestant_group2(population_finesse.end() - k_tournement_contestant, population_finesse.end());

                //pick 2 best

                int better = contestant_group1[0].second, best = contestant_group2[0].second; 
                string parent1 = contestant_group1[0].first, parent2 = contestant_group1[0].first;
                for (int j = 0; j < contestant_group1.size(); j++)
                {
                    if(contestant_group1[j].second > better){
                        parent1 = contestant_group1[j].first;
                        better = contestant_group1[j].second;  
                    }
                    if(contestant_group2[j].second > best){
                        parent2 = contestant_group2[j].first;
                        best = contestant_group2[j].second;
                    }
                }
                if(flag) 
                cout << better << ", " << best << endl;
            

                //for ( int i = 0; i < k_tournement_contestant ; i++)cout << contestant_group1[i].first << ", " << contestant_group1[i].second << endl;
                //cout << endl<<endl;
                //for ( int i = 0; i < k_tournement_contestant ; i++)cout << contestant_group2[i].first << ", " << contestant_group2[i].second << endl;
            
                //Make them procreate
                vector<string> tmp = mating(parent1,parent2);
                generation_new.insert(generation_new.end(), tmp.begin(), tmp.end());
            
            }
            //for (string x : generation_new) cout << x << endl<<endl;

            if (generation_new.size() < n_population){
                shuffle(population.begin(), population.end(), mt19937{std::random_device{}()});
                for (int i = generation_new.size() - 1 ; i < n_population; i++) generation_new.push_back(population[i]);
            }
            else if ( generation_new.size() > n_population){
                for (int i = generation_new.size() ; i > n_population; i--) generation_new.pop_back();
            }
            population.clear();
            ////////////////////////////////////////////////////////////////

            //Choose population changer criteria????????????????????????????
            
            ////////////////////////////////////////////////////////////////
            population.insert(population.begin(), generation_new.begin(), generation_new.end());
            generation_new.clear();
        }

        // HERE GOES YOUR GREEDY HEURISTIC
        // When finished with generating a solution, first take the computation time as explained above. Say you store it in variable 'time'.
        // Then write the following to the screen: cout << "value " << <value of your solution> << "\ttime " << time << endl;

        cout << "end file " << inputFiles[na] << endl;
    }

    // calculating the average of the results and computation times and write them to the screen
    double r_mean = 0.0;
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    cout << r_mean << "\t" << t_mean << endl;
    
}

