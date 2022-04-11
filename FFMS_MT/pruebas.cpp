#include <string>
#include <vector>

using namespace std;

class Prueba{
    public:
        string testeo;
};

int hammingDist(Prueba str1){
            int i = 0, count = 0;
            for ( std::string::iterator it=str1.testeo.begin(); it!=str1.testeo.end(); ++it) 
            {
                if (str1.testeo.at(i) != str1.testeo.at(i))
                    count++;
                i++;
            }
            return count;
}

vector<Prueba> individuos;

int main(){
    for (int i = 0; i < 10; i++)
    {
        Prueba nuevo;
        nuevo.testeo = to_string(i);
        individuos.push_back(nuevo);
    }

    Prueba x = individuos[2];
    

}


