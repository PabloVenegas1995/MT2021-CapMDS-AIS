#include <string>
#include <cstring>
#include <vector>
#include <set>
#include <utility> // PAIR
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


// Data structures for the problem data
vector<string> input_sequence;
int n_of_sequences;
int sequence_length_m;
int alphabet_size = 4;
map<int,char> mapping;
map<char,int> rev_mapping;
//double threshold;
int threshold;
int t_value;


// vector for keeping all the names of the input files
vector<string> inputFiles;

// dummy parameter as an example for creating command line parameters -> see function read_parameters(...)
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


int main( int argc, char **argv ) {

    read_parameters(argc,argv);

    // number of input files
    int n_files = int(inputFiles.size());

    // vectors for storing the result and the computation time obtained by the applications of the greedy heuristic
    vector<double> results(n_files, std::numeric_limits<int>::min());
    vector<double> times(n_files, 0.0);

    // letter to index mapping (note: this only works for instances on alphabet Sigma= {A, C, T, G}
    mapping[0] = '0';
    mapping[1] = '1';
    rev_mapping['0'] = 0;
    rev_mapping['1'] = 1;

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
        sequence_length_m = input_sequence[0].size();
        indata.close();

        cout << "N of sequences "<< n_of_sequences << endl;
        for (int i = 0; i < n_of_sequences; i++) cout << input_sequence[i] << endl;

        // minimum required Hamming distance
        //t_value = int(threshold * sequence_length_m);
        cout << "Threshold given " << threshold << endl;
        t_value = threshold;

        // the computation time starts now
        //clock_t start = clock();

        // Example for requesting the elapsed computation time (in seconds) at any moment: 
        // clock_t end = clock();
        // double elapsed = double(end - start)/CLOCKS_PER_SEC;

        cout << "start file " << inputFiles[na] << endl;


////////////////////// CODE GOES HERE /////////////////////


//////////////////////////////////////////////////////////
        cout << "end file " << inputFiles[na] << endl;

    }

    return 0;
}