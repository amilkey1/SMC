#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
using namespace std;

#include "lot.hpp"

extern proj::Lot rng;

namespace proj {

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
        void                                    showParticle();
        void                                    proposal();
        void                                    reweightParticles();
        void                                  resampleParticles(int n);
//        ~Particle();
        void                                    setData(Data::SharedPtr d) {
            _data = d;
            _forest.setData(d);
            }
        double                                  sumParticleWeights(double particle_sum);
        void                normalizeParticleWeights(double particle_sum);
    private:
        // data members of Particle class
        // I suggest using the underscore convention for data members
        Forest                              _forest;
        double                              _weight;
//        double                              _sum_particles;
        Particle(const Particle & other);
        Data::SharedPtr                     _data;
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

inline void Particle::proposal() {
    _forest.proposeParticles(); //proposal step
}

inline Particle::Particle(const Particle & other) {
    assert(false);
}

//reweight a particle after proposal
inline void Particle::reweightParticles() {
    //this function should assign weight proportional to likelihood of each forest in each particle
    //for now, assign weights randomly
    _weight=rng.uniform();
}

inline double Particle::sumParticleWeights(double particle_sum){
    particle_sum = particle_sum + _weight;
    return particle_sum;
}

inline void Particle::normalizeParticleWeights(double particle_sum){
    _weight = _weight / particle_sum;
}
inline void Particle::resampleParticles(int n){
    // prune out particles with low weights
    // draw K uniform random variables between 0 and 1
    // save particle if n falls within particle interval
}
}
