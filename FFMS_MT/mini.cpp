#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <utility>
#include <string.h>
#include <algorithm>

using namespace std;

vector<string> inputFiles;

int dummy_integer_parameter = 0;

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFiles.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-param1")==0) dummy_integer_parameter = atoi(argv[++iarg]); // example for creating a command line parameter param1 -> integer value is stored in dummy_integer_parameter
        iarg++;
    }
}


vector<string> input_sequence;
int n_of_sequences;
int sequence_length;
int alphabet_size = 4;

//pair<char, vector<int,int>> aa;


int main(int argc, char **argv ) {

    read_parameters(argc,argv);
    

    // number of input_sequence files
    int n_files = int(inputFiles.size());

    // vectors for storing the result and the computation time obtained by the applications of the greedy heuristic
    vector<double> times(n_files, 0.0);

    for (int na = 0; na < n_files; ++na) {

        // opening the corresponding input_sequence file and reading the problem data
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

        int a[alphabet_size][sequence_length][1];




        for (int ii = 0; ii < alphabet_size; ii++){
            for (int jj = 0; jj < sequence_length; jj++){
                a[ii][jj][0] = 0;
            }
        }

        for (int j = 0; j < n_of_sequences; j++){
            for (int k = 0; k < sequence_length; k++){
                if (input_sequence[j].at(k) == 'A')
                    a[0][k][0]++;
                else if (input_sequence[j].at(k) == 'C')
                    a[1][k][0]++;
                else if (input_sequence[j].at(k) == 'T')
                    a[2][k][0]++;
                else if (input_sequence[j].at(k) == 'G')
                    a[3][k][0]++;
            }
        }
        // for (int i = 0; i < sequence_length; i++){
        //     cout << "in position " << i << endl;
        //     cout << "As: " << a[0][i][0] <<endl;
        //     cout << "Cs: " << a[1][i][0] <<endl;
        //     cout << "Ts: " << a[2][i][0] <<endl;
        //     cout << "Gs: " << a[3][i][0] <<endl;
        // }

   
        for (int l = 0; l < sequence_length; l++){
            string out;
            if ( a[0][l][0] < a[1][l][0] && a[0][l][0] < a[2][l][0] && a[0][l][0] < a[3][l][0]) out+= 'A';
            if ( a[1][l][0] < a[0][l][0] && a[1][l][0] < a[2][l][0] && a[1][l][0] < a[3][l][0]) out+= 'C';
            if ( a[2][l][0] < a[1][l][0] && a[2][l][0] < a[0][l][0] && a[2][l][0] < a[3][l][0]) out+= 'T';
            if ( a[3][l][0] < a[1][l][0] && a[3][l][0] < a[2][l][0] && a[3][l][0] < a[0][l][0]) out+= 'G';
            if (a[0][l][0] == a[1][l][0])      out+= 'A';
            else if (a[0][l][0] == a[2][l][0]) out+= 'A';
            else if (a[0][l][0] == a[3][l][0]) out+= 'A';
            

        }


    }
}