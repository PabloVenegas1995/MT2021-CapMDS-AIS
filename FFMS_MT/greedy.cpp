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

double t_limit;

int flag = 0;

//Genetic Algorithm Hyper-Parameters
int n_population = 50; //Cambiar por un limite de tiempo
int n_tournement = 1;
double p_tournement_contestans = 0.1; 
int crossover_type = 0; // 0 uniform; 1 single point crossover.

int r_population = 0;   // 0 new + random_previous; 1 new + elite; 2 new + elite + random_previous;
double p_newGen = 0.7;       // parameter of percentage of new offspring permited for r_population 1 & 2
int nng_value;
double p_elite     = 0.2;    // parameter of percentage of elites
int elite_value;
int rndm_replace_type = 0; //0 select random from this gen; 1 create completely random new individual
double mutation_rate = 0.05; //parameter of probability to change a individual
int mutation_type = 0;  // 0 change sequence of individual for a new random; 1 change allel(char) of individual under a percentage
double p_allele_mutation = 0.1; //parameter of probability to change an allel
int forced_mutation = 0; // 0 change but can be equal, 1 force change on the allele


// Genetic Algorithm Data Structures
vector<pair<int,pair<string,int>>> population;
vector<pair<int,pair<string,int>>> generation_new;
vector<pair<int,pair<string,int>>> gen_elites;




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
        else if (strcmp(argv[iarg],"-tlim")==0) t_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-np")==0) n_population = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-ptc")==0) p_tournement_contestans = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-ct")==0) crossover_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-png")==0) p_newGen = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pe")==0) p_elite = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-rp")==0) r_population = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-rrt")==0) rndm_replace_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-mr")==0) mutation_rate = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-mt")==0) mutation_type = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-pam")==0) p_allele_mutation = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-fm")==0) forced_mutation = atoi(argv[++iarg]);        
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

void fitness_calculation(){

    //cout << t_value << " tvalue " << endl;
    for (int i = 0; i < population.size(); i++)
    {
        int finesse = 0;
        int missmatchs = 0;
        for (int j = 0; j < n_of_sequences; j++)
        {
                missmatchs = hammingDist(population[i].second.first, input_sequence[j]);
                //cout << "missmatch for; [" << i << "] " << missmatchs << endl;
                if (missmatchs >= t_value) finesse += 1;
        }
        population[i].second.second = finesse;
    }
    //cout << population[0].second.first << endl;        
}

int fitness_offspring(string parent){
    int finesse = 0;
    int missmatchs = 0;
    for (int j = 0; j < n_of_sequences; j++){
        missmatchs = hammingDist(parent, input_sequence[j]);
        if (missmatchs >= t_value) finesse += 1;
    }
    return finesse;
}

void mating(pair<int,pair<string,int>> parent1, pair<int,pair<string,int>> parent2){

    string child1 = "", child2 = "";
    if (crossover_type == 0){ //Uniform Crossover
        for (int i = 0 ; i < sequence_length; i++){
            random_device rd;
            mt19937 rng(rd());
            uniform_int_distribution<mt19937::result_type> dist(0, 1);    
            if ( dist(rng) == 0){
                child1.push_back(parent1.second.first.at(i));
                child2.push_back(parent2.second.first.at(i));
            }
            else{
                child1.push_back(parent2.second.first.at(i));
                child2.push_back(parent1.second.first.at(i));
            }
        }
    }
    else if (crossover_type == 1){ // Crossover Point, pivote al azar
        random_device rd;
        mt19937 rng(rd());
        uniform_int_distribution<int> dist(0, sequence_length - 1); 
        int crossover_point = dist(rng);
        //if(flag)cout << crossover_point <<endl;
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
        random_device rd;
        mt19937 rng(rd());
        uniform_int_distribution<int> dist(0, sequence_length - 1); 
        int crossover_point1 = dist(rng), crossover_point2 = dist(rng);
        while (crossover_point1 == crossover_point2){
            //if(flag) cout<<"SON IGUALES"<<endl; 
            crossover_point2 = dist(rng);}
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
    
    generation_new.push_back(make_pair(parent1.first + 100,make_pair(child1,fitness_offspring(child1))));
    generation_new.push_back(make_pair(parent2.first + 100,make_pair(child2,fitness_offspring(child2))));
}

bool sortbysec( const pair<int,pair<string,int>> &a, const pair<int,pair<string,int>> &b)
{
    return (a.second.second > b.second.second);
}

void replace_population(){
            
    if (r_population == 0)
    {
        sample(population.begin(), population.end(), back_inserter(generation_new), n_population - generation_new.size(), mt19937{std::random_device{}()});
    }
    else if(r_population == 1 || r_population == 2){        
        sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
        if (r_population == 1) for (int l = 0; generation_new.size() < n_population; l++) generation_new.push_back(population[l]);
        else if (r_population == 2){        
            elite_value = n_population * p_elite;
            int rndm_pop_value = n_population - elite_value - generation_new.size();
            for(int l = 0; l < elite_value; l++){
                gen_elites.push_back(population[0]);
                population.erase(population.begin());
            }
            if (rndm_replace_type == 0)
                sample(population.begin(), population.end(), back_inserter(generation_new), rndm_pop_value, mt19937{std::random_device{}()});                
            else if (rndm_replace_type == 1){
                random_device dev;
                mt19937 rng(dev());
                uniform_int_distribution<mt19937::result_type> dist6(0,3);
                
                for (int i = 0; i < rndm_pop_value; i++)
                {
                    string str = "";
                    for (int j = 0; j < sequence_length; j++)
                    {
                        char ch = mapping[dist6(rng)];
                        str.push_back(ch);
                    }
                    generation_new.push_back(make_pair(i + 100,make_pair(str,fitness_offspring(str))));
                }
            }
        } 
    }
    population.clear();
    population.insert(population.begin(), generation_new.begin(), generation_new.end());
    generation_new.clear();        
}

void mutation(){
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist(0,100); // distribution in range [1, 6]
    uniform_int_distribution<mt19937::result_type> dist3(0, 3); 
    int mutated = 0;                
    for(int i = 0; i < population.size(); i++){
        int mutation = dist(rng);
        if(mutation <= mutation_rate * 100){
            int identity = population[i].first;
            if (mutation_type == 0){
                string str = "";
                for (int j = 0; j < sequence_length; j++){
                    char ch = mapping[dist3(rng)];
                    str.push_back(ch);
                }
                //population[i].first = identity + 22000;
                population[i].second.first = str;
            }
            else if(mutation_type == 1){
                for (int l = 0; l < sequence_length; l++){
                    int allele_mutation = dist(rng);
                    char ch = mapping[dist3(rng)];
                    if (allele_mutation <= p_allele_mutation * 100){
                        //cout << "mutando " << population[i].first;
                        if (forced_mutation == 0)
                            population[i].second.first.at(l) = ch;
                        else if (forced_mutation == 1){
                            if (ch == population[i].second.first.at(l)){
                                //cout << "explota" << endl;
                                while (ch == population[i].second.first.at(l)) ch = mapping[dist3(rng)];
                                population[i].second.first.at(i) = ch;}
                        }
                    }
                }
                //population[i].first = identity + 33000;
            }
        mutated++;
        }
    }
    //cout << mutated << " individals mutated" << endl;
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

        if(flag) cout << threshold << " threshold " << sequence_length << " sequence lenght " << t_value << " t_value" << endl;
        // the computation time starts now
        //clock_t start = clock();

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
            //cout << str <<endl<<endl; 
            population.push_back(make_pair(i,make_pair(str,0)));
        }
        //exit(0);

        // Heuristic????????
        // for (int itt = 0; itt < 10; itt++){
        //     for (int l = 0; l < sequence_length; l++){
        //         char ch = mapping[rand () % 4];
        //         while (ch == input_sequence[3].at(l)) ch = mapping[rand () % 4];
        //             population[itt].second.first.at(l) = ch;
        //     }        
        // }
        nng_value = population.size() * p_newGen;
        clock_t start, end; 
        start = clock();
        int terminate_algorithm = 1, terminate_tournements = 1;
        while (terminate_algorithm){

            //Population Finesse calculation
            fitness_calculation();

            //if(flag)for (int i = 0; i < population.size(); i++) cout << "[" << population[i].first <<"] " << population[i].second.second << endl;
            //exit(0);
            sort(population.begin(), population.end(), sortbysec);  // Nlog(N)  on worst case
            if(flag) cout << population[0].second.second << " best fitnes so far"<< endl; 

            //exit(0);

            //Tournament ARC
            //select k candidates
            int k_tournement_contestant = n_population * p_tournement_contestans; //change for a new parameter variable percentage???
            //cout << k_tournement_contestant << " number of contestants" <<endl; 
            terminate_tournements = 1;
            while (terminate_tournements){            
                vector<pair<int,pair<string,int>>> get_sample, group1, group2;
                sample(population.begin(), population.end(), back_inserter(group1), k_tournement_contestant, mt19937{std::random_device{}()});
                sample(population.begin(), population.end(), back_inserter(group2), k_tournement_contestant, mt19937{std::random_device{}()});

                // for (int l = 0; l < k_tournement_contestant * 2; l++){
                //     if (l <k_tournement_contestant) group1.push_back(get_sample[l]);
                //     else group2.push_back(get_sample[l]);
                // } 

                //pick 2 best
                pair<int,pair<string,int>> better_parent = group1[0], best_parent = group2[0];
                for (int j = 0; j < group1.size(); j++)
                {
                    if(group1[j].second.second > better_parent.second.second)
                        better_parent = group1[j];
                    if(group2[j].second.second > best_parent.second.second)
                        best_parent = group2[j];
                }
                //if(flag) cout << better_parent.second.second << ", " << best_parent.second.second << endl;
                
                //Make them procreate
                mating(better_parent, best_parent);
                if (generation_new.size() >= nng_value){
                    while(generation_new.size() > nng_value) generation_new.pop_back();
                    terminate_tournements = 0;
                } 
            }
            replace_population();
            mutation();
            //cout <<endl;
            //
            //mutacion sobre nP menos elite
            end = clock();
            double elapsed = double(end - start)/CLOCKS_PER_SEC;
            if(flag)cout << elapsed << ", elapsed"<<endl;
            if (elapsed >= t_limit)
                terminate_algorithm = 0;
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

