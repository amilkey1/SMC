//
//  clade.cpp
//  SMC
//
//  Created by Analisa Milkey on 9/24/21.
//clade class contains set of splits + edge lengths

#include "conditionals.hpp"


#include "lot.hpp"
#include "split.hpp"
using namespace strom;
using namespace std;
extern Lot rng;

class Clade {
    friend class Forest;    // lets a Forest object access private data members
    
    public:
        Clade();
        void printSplits();
        
    private:
        Split  _split;
        double _edge_length;
        void printEdgeLengths();
};

inline Clade::Clade() {
    // For now, just assign the edge length an arbitrary exponential random value
    _edge_length = rng.gamma(1,1);
    printEdgeLengths();
}

inline void Clade::printEdgeLengths() {
    //print out edge lengths associated with all the clades in each forest in each particle
    cout << "Clade edge length: " << _edge_length << "\n";
}
