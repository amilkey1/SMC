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
#include "boost/format.hpp"
using namespace std;

#include "lot.hpp"

extern Lot rng;


class Particle {
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
        void showParticle();
        void advance();
    private:
        // data members of Particle class
        // I suggest using the underscore convention for data members
        Forest _forest;
        double _weight;
};

// Constructor assigns a random weight
inline Particle::Particle() {
    // You can just let the Forest constructor automatically initialize the _forest object
    // I recommend using the Lot class for choosing the random weight.
    // Lot is a global variable created in main.cpp, so it can be used anywhere.
    _weight = rng.uniform();
};

inline void Particle::showParticle() {
    //print out weight of each particle
    cout << "Particle weight: " << _weight << "\n" ;
    cout << "Forest: " << "\n";
    _forest.showForest();
//    _forest.chooseTrees();
}

inline void Particle::advance() {
    _forest.nextStep();
}
