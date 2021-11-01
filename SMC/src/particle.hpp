#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
using namespace std;

#include "lot.hpp"

extern proj::Lot rng;

namespace proj {

class Particle {
    public:
        Particle();
        Particle(const Particle & other);

        void                                    showParticle();
        double                                  proposal();
        void                                    setData(Data::SharedPtr d) {
            _data = d;
            _forest.setData(d);
            }
        void                                    saveForest(std::string treefilename) const;
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    operator=(const Particle & other);
        const Forest &                          getForest() const {return _forest;}
        std::string                             saveForestNewick() const {
            return _forest.makeNewick(8, true);
        }
    

    private:
        Forest                                  _forest;
        double                                  _log_weight;
        Data::SharedPtr                          _data;
        double                                  _log_likelihood;
        unsigned                                _n;
};

// Constructor assigns a random weight
inline Particle::Particle() {
    _log_weight = rng.uniform();
    _n = 0;
    _log_likelihood = 0.0; //initialize log likelihood so generation 0 is able to calculate log weights
};

inline void Particle::showParticle() {
    //print out weight of each particle
    cout << "Particle weight: " << _log_weight << "\n" ;
    cout << "Forest: " << "\n";
    _forest.showForest();
}

inline double Particle::calcLogLikelihood() {
    double log_likelihood = _forest.calcLogLikelihood();
//    cout << "log likelihood equals " << log_likelihood << endl;
//    savePaupFile("paup.nex", "green4.nex", "forest.tre", log_likelihood);
//    saveForest("forest.tre");
    
    return log_likelihood;
}

inline double Particle::proposal() {
    _forest.proposeParticles(); //proposal step
    double prev_log_likelihood = _log_likelihood;
    _log_likelihood = calcLogLikelihood();
    unsigned n = Forest::_nspecies -_n;
    double log_q = log(_forest._speciation_rate) + log(n-1) - _forest._speciation_rate*(n-1)*_forest._last_edge_length - (boost::math::lgamma(n+1) - log(2) - boost::math::lgamma(n-1));
    _log_weight = _log_likelihood - prev_log_likelihood - log_q;
    _n++;
    return _log_weight;
}

inline Particle::Particle(const Particle & other) {
    *this = other;
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

inline double Particle::calcHeight() {
    double sum_height = 0.0;
    for (auto nd : _forest._preorder) {
        if (nd->getParent()!=_forest._root) {
            sum_height += nd->getEdgeLength();
        }
        if (!nd->getLeftChild()) {
            break;
        }
    }
    return sum_height;
}

inline void Particle::operator=(const Particle & other) {
    _log_weight     = other._log_weight;
    _log_likelihood = other._log_likelihood;
    _forest         = other._forest;
    _data           = other._data;
};
}
