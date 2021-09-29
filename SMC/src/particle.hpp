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
//using namespace strom;
extern strom::Lot rng;

namespace strom {
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
}
}



//importance sampling:
//0) create proposal distribution q_imp (use prior as proposal?)
//1) sample from proposal distribution with density q_imp to obtain tree t_k
//2) correct discrepancy between proposal and posterior distribution by assigning weight w_k to proposed tree
    //weight = (prior*likelihood) / (proposal distribution) ?
    //to form particle population

//sequential importance sampling:
//compute particles in stages
//proposal is a conditional density ?
//weight proposal is more complex
//to form particle population

//SMC:
//0) create particle pop
//1) resampling
//  prune particles using randomized scheme
//  throw K darts on one-dimensional target, then modify weight of each particle
//2) proposal
//3) weighting
