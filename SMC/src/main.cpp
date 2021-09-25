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
unsigned Forest::_nspecies = 5;

using namespace std;

int main(int argc, const char * argv[]) {
    //initialize vector of particles
    unsigned nparticles = 10;
    vector<Particle> my_vec(nparticles);
    
    return 0;
}
