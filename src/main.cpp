//
//  main.cpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
//


#include <iostream>
#include <vector>
using namespace std;

int main()
{
    //Create a vector of particles
    vector<int> vect;
     
        vect.push_back(1);
        vect.push_back(2);
        vect.push_back(3);
        vect.push_back(4);
        vect.push_back(5);
        vect.push_back(6);
        vect.push_back(7);
        vect.push_back(8);
        vect.push_back(9);
        vect.push_back(10);
     
        for (int x : vect)
            cout << x << " ";
     
        return 0;
}
