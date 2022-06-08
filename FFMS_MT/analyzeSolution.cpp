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

vector<string> inputFiles;
vector<string> solution;
int n_of_sequences;
int sequence_length;
int alphabet_size = 4;
map<int,char> mapping;
map<char,int> rev_mapping;
double threshold;
int t_value;
float determinism = 1.0;

int flag = 0;

vector<int> accumulated_difference;

// vector for keeping all the names of the input files
vector<string> input_sequence;

// dummy parameter as an example for creating command line parameters -> see function read_parameters(...)
int dummy_integer_parameter = 0;




void analyze(){
    for (int i = 0; i < n_of_sequences; i++)
    {
        int diff = 0;
        for (int j = 0; j < solution[0].length(); j++)
        {
            if (solution[0].at(j) != input_sequence[i].at(j))
                diff++;
        }
        cout << input_sequence[i] << " " << diff <<endl;
    }
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFiles.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-sol")==0) solution.push_back(argv[++iarg]);
        iarg++;
    }
}


int main(int argc, char **argv){
    

    srand (time(NULL));
    read_parameters(argc,argv);
    
    int n_files = int(inputFiles.size());

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
    }

    analyze();

}