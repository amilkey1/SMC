#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
using namespace std;
using namespace boost;

#include "lot.hpp"

extern proj::Lot rng;

namespace proj {

class Particle {
    public:
        Particle();
        Particle(const Particle & other);

        void                                    debugParticle(std::string name);
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
        void                                    firstPair(pair<unsigned, unsigned>);
    bool operator<(const Particle & other) const {
        return _log_weight<other._log_weight;
    }


    private:
        Forest                                  _forest;
        double                                  _log_weight;
        Data::SharedPtr                          _data;
        double                                  _log_likelihood;
        unsigned                                _n;
};

inline Particle::Particle() {
    //log weight and log likelihood are 0 for first generation
    _log_weight = 0.0;
    _n = 0;
    _log_likelihood = 0.0;
};

inline void Particle::showParticle() {
    //print out weight of each particle
    cout << "\nParticle:\n";
    cout << "  _log_weight: " << _log_weight << "\n" ;
    cout << "  _forest: " << "\n";
//    _forest.showForest();
//    cout << "  _log_likelihood: " << _log_likelihood << endl;
}

//more detailed version of showParticle
inline void Particle::debugParticle(std::string name) {
    //print out weight of each particle
    cout << "\nParticle " << name << ":\n";
    cout << "  _log_weight:               " << _log_weight                 << "\n" ;
    cout << "  _log_likelihood:           " << _log_likelihood             << "\n";
    cout << "  _forest._nleaves:          " << _forest._nleaves            << "\n";
    cout << "  _forest._ninternals:       " << _forest._ninternals         << "\n";
    cout << "  _forest._npatterns:        " << _forest._npatterns          << "\n";
    cout << "  _forest._nstates:          " << _forest._nstates            << "\n";
    cout << "  _forest._nsubtrees:        " << _forest._nsubtrees          << "\n";
//    cout << "  _forest._speciation_rate:  " << _forest._speciation_rate    << "\n";
    cout << "  _forest._last_edge_length: " << _forest._last_edge_length   << "\n";
    cout << "  _forest._new_basal_height: " << _forest._new_basal_height.first << ", " << _forest._new_basal_height.second << "\n";
    cout << "  _forest._old_basal_height: " << _forest._old_basal_height.first << ", " << _forest._old_basal_height.second << "\n";
    cout << "  newick description:        " << _forest.makeNewick(5,false) << "\n";
}

inline double Particle::calcLogLikelihood() {
    double log_likelihood = _forest.calcLogLikelihood();
    assert(!isnan (log_likelihood));
//    cout << "log likelihood equals " << log_likelihood << endl;

    return log_likelihood;
}

inline double Particle::proposal() {
    auto taxon_pair = _forest.chooseTaxaToJoin(); //proposal step
    _forest.createNewSubtree(taxon_pair.first, taxon_pair.second);
    double prev_log_likelihood = _log_likelihood;
    _log_likelihood = calcLogLikelihood();

    //cout << format("\nold weight: %.5f + %.5f\n") % _forest._old_basal_height.second % prev_log_likelihood;
    //cout << format("new weight: %.5f + %.5f\n") % _forest._new_basal_height.second % _log_likelihood;

    _log_weight = _forest._new_basal_height.second + _log_likelihood - _forest._old_basal_height.second - prev_log_likelihood;
    _n++;
//    firstPair(taxon_pair);
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


inline void Particle::firstPair(pair<unsigned, unsigned> p) {
    std::string s;
    unsigned pair_1 = p.first;
    unsigned pair_2 = p.second;

    for (auto nd:_forest._preorder) {
        if (nd->_number == pair_1) {
            s += boost::str(boost::format("%s (%d) || ")%nd->_name %nd->_number);
        }
        else if (nd->_number == pair_2) {
            s += boost::str(boost::format("%s (%d) | log likelihood is %.5f | branch length is %d ||")%nd->_name %nd->_number %_log_likelihood %nd->_edge_length);
        }
    }
    cout << s << endl;
    showParticle();

}
}
