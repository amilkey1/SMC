//
//  particle.cpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
// particle class contains a forest + weight

#pragma once
#include <vector>
using namespace std;

class Particle {
    friend class forest;
    double weight;
    //Forest forest;
};

vector<int> assignWeight (vector<int> vec) {
    //multiply each particle by a random weight
    for (int particle : vec) {
        particle = particle*(((double) rand() / (RAND_MAX))); //chooses random value btw 0 and 1
        return vec;
    };
    
vector<int> assignForest (vector<int> vec) {
    for (int particle : vec) {
        particle = particle*forest; //???
        return vec;
    }
};
    

//return particle population containing complete states
