#include <iostream>
#include <fstream>
#include "fress.h"
using namespace std;

int main(int argc, char *argv[]) { 
    cout.precision(3);

  string filename, hp_sequence;

    filename=argv[1];
    std::ifstream in(filename);
    vector <int> sequence;
    getline(in, hp_sequence);
    cout << hp_sequence<< endl;
    for (char c: hp_sequence){
        if(c=='H'){
            sequence.push_back(1);
        }
        if(c=='P'){
            sequence.push_back(0);
        }

    }
    cout << "length is " << sequence.size() << endl;



Protein p(sequence);


    int e = p.count_contacts();
    cout << "Start energy : " << e << endl;

    time_t start, end;
    time (&start);


    p.find_minimum();
    time (&end);

    double dif = difftime(end, start);
    printf ("Time = %lf \n", dif);
    cout <<  " Number of moves : " <<p.iterator << std:: endl;
    cout << "Founden minimum : " << p.min_E << endl;

    cout << endl;
    return 0;
}