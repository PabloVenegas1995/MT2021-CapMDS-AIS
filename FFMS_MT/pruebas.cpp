
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


#include <omp.h>

using namespace std;


vector<string> input_sequence;
int n_of_sequences;
int sequence_length;
int alphabet_size = 4;
map<int,char> mapping;
map<char,int> rev_mapping;
double threshold;

int t_value = 225;
vector<string> inputFiles;


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

vector<pair<string, int>> population_finesse;

void fitness_calculation(vector<string> population_tested, vector<string> initial_input){
    cout << t_value << " tvalue"<< endl;
    for (int i = 0; i < population_tested.size(); i++)
    {
        int finesse = 0;
        int missmatchs = 0;
        for (int j = 0; j < n_of_sequences; j++)
        {
                missmatchs = hammingDist(population_tested[i], initial_input[j]);
                cout << "missmatchs found for: ["<< i <<"] " << missmatchs << endl;
                if (missmatchs >= t_value) finesse += 1;
        }
        cout << "finesse for [" << i << "]: " << finesse << endl;
        population_finesse.push_back( make_pair(population_tested[i], finesse) );        
    }   
}

int dummy_integer_parameter = 0;

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFiles.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-th")==0) threshold = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-param1")==0) dummy_integer_parameter = atoi(argv[++iarg]); // example for creating a command line parameter param1 -> integer value is stored in dummy_integer_parameter
        iarg++;
    }
}



vector<string> population;


int main(int argc, char **argv){
    
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

        t_value = int(threshold * sequence_length);

       cout << threshold << " threshold " << sequence_length << " sequence lenght " << endl;


        random_device dev;
        mt19937 rng(dev());
        uniform_int_distribution<mt19937::result_type> dist6(0,3); // distribution in range [1, 6]
        
        //Population Creation
        for (int i = 0; i < 3; i++)
        {
            string str = "";
            for (int j = 0; j < sequence_length; j++)
            {
                char ch = mapping[dist6(rng)];
                str.push_back(ch);
            }
            cout << str << endl;
            population.push_back(str);
        }
        
        clock_t start = clock();

        for (int i = 0; i < 1; i++)
        {
        
        fitness_calculation(population, input_sequence);

        // Example for requesting the elapsed computation time (in seconds) at any moment: 
        // clock_t end = clock();
        // double elapsed = double(end - start)/CLOCKS_PER_SEC;

        }

        clock_t end = clock();
        double elapsed = double(end - start)/CLOCKS_PER_SEC;
        cout << elapsed << endl;

        
    }

//    Prueba x = individuos[2];
    

}


