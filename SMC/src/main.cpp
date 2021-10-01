//
//  main.cpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
//

#include <iostream>
#include <vector>

#include "conditionals.hpp"
#include "particle.hpp"

// Initialize our random number generator here so it will be a global variable
#include "lot.hpp"
Lot rng;

// Must initialize _nspecies here because it is a static data member of Forest
#include "forest.hpp"
unsigned Forest::_nspecies = 4;
const double Node::_smallest_edge_length=1.0e-12;

using namespace std;

int main(int argc, const char * argv[]) {
    rng.setSeed(12345);
    unsigned nspecies=4;
    Forest::setNumSpecies(nspecies);
    //initialize vector of particles
    unsigned nparticles = 1;
    vector<Particle> my_vec(nparticles);
    
    cout << "\n Particles at the start: " << endl;
    
    for (auto p:my_vec ) {
        p.showParticle();
    }

    for (unsigned g=0; g<nspecies-2; g++){
        cout << "\n Particles after generation " << g << endl;
        for (auto p:my_vec ) {
            p.advance();
            p.showParticle();
        }
    }
    return 0;
}
