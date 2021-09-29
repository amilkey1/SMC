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
strom::Lot rng;

// Must initialize _nspecies here because it is a static data member of Forest
#include "forest.hpp"
unsigned strom::Forest::_nspecies = 4;
const double strom::Node::_smallest_edge_length=1.0e-12;

using namespace std;

int main(int argc, const char * argv[]) {
    //initialize vector of particles
    unsigned nparticles = 10;
    vector<strom::Particle> my_vec(nparticles);
    
    for (auto p:my_vec ) {
        p.showParticle();
    }
    return 0;
}
