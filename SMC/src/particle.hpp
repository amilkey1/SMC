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
        ~Particle();
    private:
        // data members of Particle class
        // I suggest using the underscore convention for data members
        Forest _forest;
        double _weight;
        Particle(const Particle & other);
        void pruneParticles();
        void resampleParticles();
};

// Constructor assigns a random weight
inline Particle::Particle() {
    // You can just let the Forest constructor automatically initialize the _forest object
    // I recommend using the Lot class for choosing the random weight.
    // Lot is a global variable created in main.cpp, so it can be used anywhere.
    _weight = rng.uniform();
    resampleParticles();
};

inline void Particle::showParticle() {
    //print out weight of each particle
    cout << "Particle weight: " << _weight << "\n" ;
    cout << "Forest: " << "\n";
    _forest.showForest();
}

inline void Particle::advance() {
    _forest.nextStep();
}

inline Particle::Particle(const Particle & other) {
    assert(false);
}

inline void Particle::resampleParticles(){
    //prune particles with low weight
    //expected number of times each particle is resampled is proportional to particle weight
    
    //so let's say # times resampled = particle weight
    
    //prune lower half?
    //not sure when this step happens
    if (_weight < 0.5) {
//        _weight = 0;
        Particle::~Particle();
    }
}

inline Particle::~Particle() {
    std::cout << "Destroying a Particle" << std::endl;
    
//    delete &_forest;
//    delete &_weight;
}
