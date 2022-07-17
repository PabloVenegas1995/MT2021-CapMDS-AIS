/***************************************************************************
                    hybrid-metaheuristic.cpp  -  description
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


//Data Structures for Greedy
vector<int> accumulated_difference;
float determinism = 1.0;

//Data Structures for CMSA
int sols = 1;
int max_age = 1;


//Tuning mode
int tuning = 0;



// vector for keeping all the names of the input files
vector<string> inputFiles;

// computing time limit for each application of the hybrid metaheuristic
double time_limit = 90.0;

// computing time limit for each application of CPLEX within the hybrid metaheuristic
double cplex_time_limit = 5.0;


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
        else if (strcmp(argv[iarg],"-tlim")==0) time_limit = atoi(argv[++iarg]); // reading the computation time limit from the command line (if provided)
        else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-d")==0) determinism = stof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-sols")==0) sols = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-cpl_t")==0) cplex_time_limit = atoi(argv[++iarg]);// reading t he computation time limit of CPLEX for the application in the hybrid metaheuristic from the command line (if provided)
        else if (strcmp(argv[iarg],"-max_age")==0) max_age = atoi(argv[++iarg]);
        else if (strcmp(argv[iarg],"-tuning")==0) tuning = atoi(argv[++iarg]);
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
            if(!tuning) cout << "value " << nv << "\ttime " << newTime << "\tgap " << newGap << endl;
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
        for(int i=0; i<sequence_length;i++){
            for(int j=0; j<alphabet_size;j++){
                if(subinstance[i][j]==-1){
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
                if(!tuning) cout << "value " << lastVal << "\ttime " << newTime << "\tgap " << lastGap << endl;
                results[na] = lastVal;
                times[na] = newTime;
                gaps[na] = lastGap;
            }

            //Adapt: increasing ages in subinstance
            //cout<<"Adapt: increasing..."<<endl;
            for(int i=0; i<sequence_length;i++){
                for(int j=0; j<alphabet_size;j++){
                    if(subinstance[i][j]>=0){
                        subinstance[i][j]++;
                    }
                }
            }
            //cout<<"Adapt: zeroing..."<<endl;
            for(int i=0;i<sequence_length;i++){
                for(int j=0;j<alphabet_size;j++){
                    IloNum xval =cpl.getValue(x[i][j]);
                    if(xval>0) cout<<mapping[j];  //solución aquí!

                    //Adapt: zeroing subinstance
                    if(xval>0){
                        subinstance[i][j]=0;
                    }
                    
                    



                }
            }
        }
        env.end();
        //Adapt: removing older elements in subinstance
        //cout<<"Adapt: removing..."<<endl;
        for(int i=0; i<sequence_length;i++){
            for(int j=0; j<alphabet_size;j++){
                if(subinstance[i][j]>=max_age){
                    subinstance[i][j]=-1;
                }
            }
        }
    }
    catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
    }
    env.end();
}

//Greedy functions for FFMSP

int calculate_best_heuristic_option(int pos, int t_value) {
   
    int previous_value = 0;
    //calculado previous-values
    for(int i=0; i<n_of_sequences;i++)
        if(accumulated_difference[i]>=t_value) previous_value++;
    //calculando after-values
    //calcular el valor con cada una de las 4 opciones
    int valores[alphabet_size];
    int best_option = -1;
    int best_value = 0;
    int add;
    for(int i=0; i<alphabet_size; i++) {
        int actual_value=0;
        for(int j=0;j<n_of_sequences;j++){
            if(input_sequence[j][pos]!=mapping[i]) add = 1; else add = 0;
            if(accumulated_difference[j]+add >= t_value) 
                actual_value++;    
        }
        valores[i] = actual_value;
    }
    int firstS,secondS;
    //resolver empates...
    if(valores[0]>valores[1]) firstS = 0;
        else if(valores[0]<valores[1]) firstS = 1; else firstS = rand()%2;
    if(valores[2]>valores[3]) secondS = 2;
        else if(valores[2]<valores[3]) secondS = 3; else secondS = rand()%2+2;
    if(valores[firstS]>valores[secondS]) best_option = firstS;
        else if(valores[firstS]<valores[secondS]) best_option = secondS;
            else  if(rand()%2) best_option = firstS;
                else best_option = secondS;

    if(previous_value < valores[best_option]) return best_option; else return -1;
}



//complementary greedy function

void show_accumulated_difference() {
    for(int i=0;i<accumulated_difference.size();i++) {
        cout <<"seq: "<<i<<" "<< accumulated_difference[i] << " ";
        cout<<endl;
    }
}

int calculate_solution_quality() {
    int solution_quality = 0;
    for(int i=0;i<accumulated_difference.size();i++) {
        if(accumulated_difference[i]>=t_value) solution_quality++;
    }
    return solution_quality;

}

/**********
Main function
**********/

int main( int argc, char **argv ) {


    // random seed initilization
    srand (time(NULL));
    read_parameters(argc,argv);
    
    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // number of input files
    int n_files = int(inputFiles.size());

    // vectors for storing the result and the computation time obtained by the applications of the hybrid metaheuristic
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


        //Data Structures for CMSA
        //int subinstance[4][sequence_length]={-1};
        vector < array<int,4> >  subinstance;
        for (int i=0; i<sequence_length; ++i) {
            array <int, 4> aux = {-1,-1,-1,-1};
            subinstance.push_back(aux);
        }



        
        // the computation time starts now
        //clock_t start = clock();
        Timer timer;

        // Example for requesting the elapsed computation time (in seconds) at any moment: 
        //clock_t end = clock();
        //double elapsed = double(end - start)/CLOCKS_PER_SEC;

        if(!tuning) cout << "start file " << inputFiles[na] << endl;

        // HERE GOES YOUR HYBRID METAHEURISTIC

        // For implementing the hybrid metaheuristic you probably want to take profit from the greedy heuristic and/or the local search method that you already developed.
        // Whenever the best found solution is improved, first take the computation time as explained above. Say you store it in variable ct.
        // Then, write the following to the screen: cout << "value " << <value of the new best found solution> << "\ttime " << ct << endl;
        // Moreover, store the value of the new best found solution in vector results: results[na] = <value of the new best found solution>;
        // And store the current computation time (that is, the time measured at that moment and stored in variable "ct") in vector times: times[na] = ct;

        // Stop the execution of the hybrid metaheuristic once the time limit "time_limit" is reached.

        //Calculing support numbers for Greedy heuristics
        int n_of_X[sequence_length][4]={0};        

        for (int i=0;i<sequence_length;i++)
            for(int j=0;j<n_of_sequences;j++){
                n_of_X[i][rev_mapping[input_sequence[j][i]]]++;
                
            }

        
        //Identificando el menor en 0 por caracter   
        // ***Nota***: aquí se toma el primer menor, se puede mejorar para que tome de todos los posibles menores ( implementar tiebreak) (tie break implementado)
        // encontrar el numero menor por posición...luego
        //Transformar n_of_X en un vector de conjuntos, y agregar ahí cada vez que se encuentre menor o igual
        //construir en gredy seleccionando al azar del conjunto cuando sea necesario.
        int menor_n_of_X[sequence_length]={0};
        int numero_menor_n_of_X[sequence_length]={0};

        for(int i=0; i<sequence_length;i++){
            int lower = 10000;
            for(int j=0;j<4;j++){
                if( (n_of_X[i][j] < lower) || (n_of_X[i][j] == lower && rand()%100 <0.5)) {
                    lower = n_of_X[i][j];
                    menor_n_of_X[i] = j;
                    numero_menor_n_of_X[i]=lower;
                    }
                
                  
                }
        }
        //Generating vector of maximum difference and stablish initial point
        int lower =10000;
        int lower_pos = -1;
        for(int i=0; i<sequence_length;i++){
            //cout<<mapping[menor_n_of_X[i]];
            if (lower > numero_menor_n_of_X[i]) {
                lower = numero_menor_n_of_X[i];
                lower_pos = i;
            }
        }



        //WHILE {} ciclo principal de metaheuristica
        while(timer.elapsed_time(Timer::VIRTUAL)<time_limit){

            //Construct phase
            for(int cs=0; cs<sols;cs++){
                //Structure initialization for Greedy
                accumulated_difference= vector<int>(n_of_sequences,0); //carefull con la ram
                string solution;
                solution.resize(sequence_length);
                int pos = lower_pos;
                int i = 0;
                int number_of_sequences_over_treshold = 0;
                if(rand() % 100 < determinism * 100) {
                    solution[pos] = mapping[menor_n_of_X[pos]];
                    //cout<<solution<<endl;
                    //update the vector of accumulated difference
                    for(int j=0;j<n_of_sequences;j++){
                        if(input_sequence[j][pos] != solution[pos]){
                            accumulated_difference[j]++;
                        }
                    }
                    //show_accumulated_difference();
                    //calculate_best_heuristic_option(pos + 1,t_value);
                }
                else{
                    solution[pos] = mapping[rand()%4];
                }
            
                i++;
                pos++;
                pos=pos%sequence_length; 
                int h_option;
                
                while(i<sequence_length){
                //verify whats is better the lower recurrent or other in terms of number of sequences over threshold
                if(rand() % 100 < determinism * 100) {
                    h_option=calculate_best_heuristic_option(pos,t_value);
                    if(h_option==-1) solution[pos] += mapping[menor_n_of_X[pos]]; //importante: considerar que hay realmente varias opciones quizás..
                    else{
                        solution[pos] += mapping[h_option];
                    }
                    //update the vector of accumulated difference
                }
                else{solution[pos]=mapping[rand()%4];}

                // update the vector of accumulated difference
                for(int j=0;j<n_of_sequences;j++){
                    if(input_sequence[j][pos] != solution[pos]){
                    accumulated_difference[j]++;
                    }
                }

                i++;
                pos++;
                pos=pos%sequence_length; 
            
                }
            
            //solution constructed
            //cout<<solution<<endl;
            //cout<<"calidad: "<<calculate_solution_quality()<<endl;


            //Merge
            for(int idx=0;idx<sequence_length;idx++){
                if(subinstance[idx][rev_mapping[solution[idx]]]==-1){
                    //case: first time seeing this character
                    subinstance[idx][rev_mapping[solution[idx]]]=0;
                }
                else{
                    //case: already seen this character
                    //nothing to do.
                }

            }
            
            
            
            }

        //solve
        run_cplex(subinstance, timer, results, times, gaps, na);
        }

        //show_subinstance(subinstance);
        if (!tuning) cout << "end file " << inputFiles[na] << endl;
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
    if(!tuning) cout << r_mean << "\t" << t_mean << endl; else {cout << -1 * r_mean <<endl;}
}

