#include <iostream>
#include <map>

int n_sequences = 10;
int sequence_length = 20;
std::map<int,char> mapping;

int main(){
    mapping[0] = 'A';
    mapping[1] = 'C';
    mapping[2] = 'T';
    mapping[3] = 'G';
    srand(time(0));
    for (int i = 0; i < n_sequences; i++){
        std::string str = "";
            for (int j = 0; j < sequence_length; j++){
                char ch = mapping[rand() % 4];
                str.push_back(ch);
            }
            std::cout << str << std::endl; 
        }
}