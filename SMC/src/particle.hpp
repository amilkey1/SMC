//
//  particle.cpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
// particle class contains a forest + weight

#pragma once
#include <vector>
#include "conditionals.hpp"
#include "forest.hpp"
using namespace std;

#if defined(POLSUGGESTION)
#include "lot.hpp"
using namespace strom;
extern Lot rng;
#endif

class Particle {
#if defined(POLSUGGESTION)
    // You probably don't need to declare Forest to be a friend of Particle
    // because Forest objects will never need to mess with private members of Particle
    // We will need to make Particle a friend of Forest, however, so that Particle 
    // objects can access private members of their Forest data member  
    //friend class forest;

    public:
    
        // Although a default constructor will be provided for you, it is not
        // a bad idea to create your own constructor. That would allow you to 
        // assign a random weight as it is being constructed
        Particle(); 
    
        // member functions of Particle class
        
    private:
        // data members of Particle class
        // I suggest using the underscore convention for data members
        Forest _forest;
        double _weight;
        
#else
    friend class forest;
    double weight;
    //Forest forest;
#endif
};

#if defined(POLSUGGESTION)
// Constructor assigns a random weight
inline Particle::Particle() {
    // You can just let the Forest constructor automatically initialize the _forest object
    // I recommend using the Lot class for choosing the random weight.
    // Lot is a global variable created in main.cpp, so it can be used anywhere.
    _weight = rng.uniform();
};
#else
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
#endif
    

//return particle population containing complete states
