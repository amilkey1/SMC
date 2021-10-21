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
        void                                   resampleParticles(double running_sum);
//        ~Particle();
        void                                    setData(Data::SharedPtr d) {
            _data = d;
            _forest.setData(d);
            }
        double                                  sumParticleWeights(double particle_sum);
        void                                    normalizeParticleWeights(double particle_sum);
        void    saveForest(std::string treefilename) const;
    void savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
    double calcLogLikelihood();

    private:
        // data members of Particle class
        // I suggest using the underscore convention for data members
        Forest                              _forest;
        double                              _weight;
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

inline double Particle::calcLogLikelihood() {
    double log_likelihood = _forest.calcLogLikelihood();
    cout << "log likelihood equals " << log_likelihood << endl;
//    savePaupFile("paup.nex", "green4.nex", "forest.tre", log_likelihood);
//    saveForest("forest.tre");
    
    return log_likelihood;
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
//    _weight=rng.uniform();
}

inline double Particle::sumParticleWeights(double particle_sum){
    particle_sum = particle_sum + _weight;
    return particle_sum;
}

inline void Particle::normalizeParticleWeights(double particle_sum){
    _weight = _weight / particle_sum;
}
inline void Particle::resampleParticles(double running_sum){
    // save particle if n falls below running_sum
}

inline void Particle::saveForest(std::string treefilename) const {
    ofstream treef(treefilename);
    treef << "#nexus\n\n";
    treef << "begin trees;\n";
    treef << "  tree test = [&R] " << _forest.makeNewick(8, true)  << ";\n";
    treef << "end;\n";
    treef.close();
}

inline void Particle::savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const {
    ofstream paupf(paupfilename);
    paupf << "#nexus\n\n";
    paupf << "begin paup;\n";
    paupf << "exe " << datafilename << ";\n";
    paupf << "set crit=like forcepoly;\n";
    paupf << "lset nst=1 basefreq=equal rates=equal pinvar=0;\n";
    paupf << "gettrees file=" << treefilename << " storebrlen;\n";
    paupf << "lscores 1 / userbrlen;\n";
    paupf << "[!expected lnL = " << expected_lnL << "]\n";
    paupf << "end;\n";
    paupf.close();
}
}
