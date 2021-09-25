//
//  forest.hpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
// creates the forest class
//forest class contains a set of clades

#pragma once

#include <stack>
#include <memory>
#include <iostream>
#include "conditionals.hpp"

#if defined(POLSUGGESTION)
#include <list>
#include "clade.hpp"
#endif
using namespace std;

class Forest {
    //set of clades
    //construct a vector of clades
#if defined(POLSUGGESTION)
    public:
        Forest();
        
    private:
        static unsigned _nspecies;
        std::list<Clade> _clades;
#endif
    
};

#if defined(POLSUGGESTION)
inline Forest::Forest() {
    for (unsigned i = 0; i < _nspecies; i++) {
        Clade clade;
        clade._split.resize(_nspecies);
        clade._split.setBitAt(i);
        _clades.push_back(clade);
    }
}
#endif
