/***************************************************************************
                          cplex.cpp  -  description
                             -------------------
    begin                : Thu Dec 16 2021
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

#include <string>
#include <stdio.h>
#include <time.h>
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
#include <random>
#include <chrono>
  // the following "include" is necessary for the correct working/compilation of CPLEX.
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

// vector for keeping all the names of the input files
vector<string> inputFiles;

// time limit for CPLEX (can be supplied to the algorithm via the -t comand line parameter)
double time_limit = 3200.0;


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
        else if (strcmp(argv[iarg],"-tlim")==0) time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        iarg++;
    }
}

ILOSOLVECALLBACK4(loggingCallback,
    clock_t&, start,
    double&, time_stamp,
    double&, result,
    double&, gap) {

    if (hasIncumbent()) {
        IloNum nv = getIncumbentObjValue();
        clock_t current = clock();
        double newTime = double(current - start) / CLOCKS_PER_SEC;
        double newGap = 100.0 * getMIPRelativeGap();
        if (result < double(nv)) {
            cout << "value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
            result = double(nv);
            time_stamp = newTime;
            gap = newGap;
        }
    }
}


void run_cplex(clock_t& start, vector<double>& results, vector<double>& times, vector<double>& gaps, int& na) {


    IloEnv env;
    // this instruction redirects the standard output of CPLEX to the null stream in order to avoid printing it to the screen
    env.setOut(env.getNullStream());
    try{

        IloModel model(env);

        // Here you have to implement the ILP model for its resolution with CPLEX


        IloCplex cpl(model);

        cpl.setParam(IloCplex::TiLim, time_limit);
        cpl.setParam(IloCplex::EpGap, 0.0);
        cpl.setParam(IloCplex::EpAGap, 0.0);
        cpl.setParam(IloCplex::Threads, 1);
        cpl.setWarning(env.getNullStream());
        cpl.use(loggingCallback(env, start, times[na], results[na], gaps[na]));

        cpl.solve();
    
        if (cpl.getStatus() == IloAlgorithm::Optimal || cpl.getStatus() == IloAlgorithm::Feasible) {
            clock_t current = clock();
            double newTime = double(current - start) / CLOCKS_PER_SEC;
            double lastVal = double(cpl.getObjValue());
            double lastGap = 100.0 * cpl.getMIPRelativeGap();
            if (lastGap < 0.0) lastGap *= -1.0;
            if (lastVal > results[na]) {
                cout << "value " << lastVal << "\ttime " << newTime << "\tgap " << lastGap << endl;
                results[na] = lastVal;
                times[na] = newTime;
                gaps[na] = lastGap;
            }
            if (cpl.getStatus() == IloAlgorithm::Optimal) cout << "optimality proven" << endl;
        }
        env.end();
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
    }
    env.end();
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

    // vectors for storing the result and the computation time obtained by the applications of CPLEX
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
            cout << "Error: file could not be opened" << endl;
        }

        input_sequence.clear();
        string seq;
        while (indata >> seq) input_sequence.push_back(seq);
        n_of_sequences = input_sequence.size();
        sequence_length = input_sequence[0].size();
        indata.close();

        // minimum required Hamming distance
        t_value = int(threshold * sequence_length);

        // the computation time starts now
        clock_t start = clock();

        cout << "start file " << inputFiles[na] << endl;

        run_cplex(start, results, times, gaps, na);

        cout << "end file " << inputFiles[na] << endl;
    }

    // calculating the average of the results, computation times, and optimality gaps and write them to the screen
    double r_mean = 0.0;
    double t_mean = 0.0;
    double g_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        t_mean = t_mean + times[i];
        g_mean = g_mean + gaps[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    g_mean = g_mean/double(gaps.size());
    cout << r_mean << "\t" << t_mean << "\t" << g_mean << endl;
}

