//
//  clade.cpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
//clade class contains set of splits + edge lengths

#include "conditionals.hpp"

#if defined(POLSUGGESTION)
#include "lot.hpp"
#include "split.hpp"
using namespace strom;
extern Lot rng;

class Clade {
    friend class Forest;    // lets a Forest object access private data members
    
    public:
        Clade();
        
    private:
        Split  _split;
        double _edge_length;
};

inline Clade::Clade() {
    // For now, just assign the edge length an arbitrary exponential random value
    _edge_length = rng.gamma(1,1);
}
#else
class Clade {
    friend class split;
};

inline Clade::assignSplit() {
    //split is type const
}

inline Clade::assignEdgeLength() {
    //
    double EdgeLength;
}
#endif
