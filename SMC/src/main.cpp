//
//  main.cpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
//

#include <iostream>
#include <vector>

//#include "conditionals.hpp"
#include "proj.hpp"
#include "particle.hpp"


// Initialize our random number generator here so it will be a global variable
#include "lot.hpp"
proj::Lot rng;

using namespace proj;

// Must initialize _nspecies here because it is a static data member of Forest
#include "forest.hpp"
std::string  Proj::_program_name        = "proj";
unsigned     Proj::_major_version       = 1;
unsigned     Proj::_minor_version       = 0;

unsigned Forest::_nspecies = 4;
const double Node::_smallest_edge_length=1.0e-12;

using namespace std;

int main(int argc, const char * argv[]) {
    Proj proj;
    try {
        proj.processCommandLineOptions(argc, argv);
        proj.run();
    }
    catch(std::exception & x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    
    rng.setSeed(1234);
    unsigned nspecies=26;
    Forest::setNumSpecies(nspecies);
    //initialize vector of particles
    unsigned nparticles = 1;
    vector<Particle> my_vec(nparticles);
    
    cout << "\n Particles at the start: " << endl;
    
    for (auto & p:my_vec ) {
        p.showParticle();
    }

    for (unsigned g=0; g<nspecies-2; g++){
        cout << "\n Particles after generation " << g << endl;
        for (auto & p:my_vec ) {
            p.advance();
            p.showParticle();
        }
    }
    return 0;
}
