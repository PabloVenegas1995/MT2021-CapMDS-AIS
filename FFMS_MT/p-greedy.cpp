/***************************************************************************
                          greedy.cpp  -  description
                             -------------------
    begin                : Fri Dec 17 2021
    copyright            : (C) 2021 by Christian Blum
    email                : christian.blum@iiia.csic.es
                         : (C) 2022 by Pedro Pinacho-Davidson
                         : ppinacho@udec.cl

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
float determinism = 1.0;

// Data Structures for the Greedy

vector<int> accumulated_difference;




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
        else if (strcmp(argv[iarg],"-d")==0) determinism = stof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-param1")==0) dummy_integer_parameter = atoi(argv[++iarg]); // example for creating a command line parameter param1 -> integer value is stored in dummy_integer_parameter
        iarg++;
    }
}


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
    //cout<<"t_value "<<t_value<<endl;
    //cout<<"Previous value "<<previous_value<<endl;
    //cout<<valores[0]<<"-"<<valores[1]<<"-"<<valores[2]<<"-"<<valores[3]<<endl;

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

    srand (time(NULL));
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

        // minimum required Hamming distance
        t_value = int(threshold * sequence_length);
        //estructura para el greedy
        accumulated_difference= vector<int>(n_of_sequences,0);

        // the computation time starts now
        clock_t start = clock();

        // Example for requesting the elapsed computation time (in seconds) at any moment: 
        // clock_t end = clock();
        // double elapsed = double(end - start)/CLOCKS_PER_SEC;

        cout << "start file " << inputFiles[na] << endl;
        //Calcular números de soporte a greedy
       
        int n_of_X[sequence_length][4]={0};        

        for (int i=0;i<sequence_length;i++)
            for(int j=0;j<n_of_sequences;j++){
                n_of_X[i][rev_mapping[input_sequence[j][i]]]++;
                
            }
        
        //Identificando el menor en 0 por caracter   
        // Nota: aquí se toma el primer menor, se puede mejorar para que tome de todos los posibles menores ( implementar tiebreak) 
        int menor_n_of_X[sequence_length]={0};
        int numero_menor_n_of_X[sequence_length]={0};

        for(int i=0; i<sequence_length;i++){
            int lower = 10000;
            for(int j=0;j<4;j++){
                if( n_of_X[i][j] < lower ) {
                    lower = n_of_X[i][j];
                    menor_n_of_X[i] = j;
                    numero_menor_n_of_X[i]=lower;
                }
                  
                }
        }

        //Interpretando vector de maxima diferencia
        int lower =10000;
        int lower_pos = -1;
        for(int i=0; i<sequence_length;i++){
            cout<<mapping[menor_n_of_X[i]];
            if (lower > numero_menor_n_of_X[i]) {
                lower = numero_menor_n_of_X[i];
                lower_pos = i;
            }
        }
        cout<<endl;
        cout<<"Posición de inicio de construcción: "<<lower_pos<<endl;


        // HERE GOES YOUR GREEDY HEURISTIC
        // When finished with generating a solution, first take the computation time as explained above. Say you store it in variable 'time'.
        // Then write the following to the screen: cout << "value " << <value of your solution> << "\ttime " << time << endl;
        

        // comencemos a construir la solución
        
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
                if(h_option==-1) solution[pos] += mapping[menor_n_of_X[pos]];
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
        cout<<solution<<endl;
        cout<<"calidad: "<<calculate_solution_quality()<<endl;




        

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

