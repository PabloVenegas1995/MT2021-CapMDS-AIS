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


void fitness_calculation(vector<string> population_tested, vector<string> initial_input){

    for (int i = 0; i < population_tested.size(); i++)
    {
        int finesse = 0;
        int missmatchs = 0;
        for (int j = 0; j < n_of_sequences; j++)
        {
                missmatchs = hammingDist(population_tested[i], initial_input[j]);
                //cout << "missmatchs found: "<< missmatchs << endl;
                if (missmatchs >= t_value) finesse += 1;
        }
        //if(flag) cout << "finesse for [" << i << "]: " << finesse << endl;
        population_finesse.push_back( make_pair(population_tested[i], finesse) );        
    }    
}
