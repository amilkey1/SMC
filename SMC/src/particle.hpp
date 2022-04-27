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
                                                    _nsubsets = d->getNumSubsets();
                                                    _data = d;
                                                    int index = 0;
                                                    _forests.resize(_nsubsets+1);
                                                    for (auto &_forest:_forests) {
                                                        _forest.setData(d, index);
                                                        index++;
                                                    }
                                                }
        void                                    mapSpecies(map<string, string> &taxon_map, vector<string> &species_names);
        void                                    saveForest(std::string treefilename) const;
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    operator=(const Particle & other);
        const vector<Forest> &                  getForest() const {return _forests;}
        std::string                             saveForestNewick() const {
                return _forests[0].makeNewick(8, true);
        }
        bool operator<(const Particle & other) const {
            return _log_weight<other._log_weight;
        }
        
        bool operator>(const Particle & other) const {
            return _log_weight>other._log_weight;
        }
    
        static void                                   setNumSubsets(unsigned n);
    
        vector<Forest> &                             getForests() {return _forests;}

    private:
    
        static unsigned                         _nsubsets;
        vector<Forest>                         _forests;
        double                                  _log_weight;
        Data::SharedPtr                          _data;
        double                                  _log_likelihood;
};

    inline Particle::Particle() {
        //log weight and log likelihood are 0 for first generation
        _log_weight = 0.0;
        _log_likelihood = 0.0;
    };

    inline void Particle::showParticle() {
        //print out weight of each particle
        cout << "\nParticle:\n";
        cout << "  _log_weight: " << _log_weight << "\n" ;
        cout << " _log_likelihood: " << _log_likelihood << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        for (auto &_forest:_forests) {
            _forest.showForest();
        }
    }

    //more detailed version of showParticle
    inline void Particle::debugParticle(std::string name) {
        //print out weight of each particle
        cout << "\nParticle " << name << ":\n";
        for (auto &_forest:_forests) {
            cout << "  _log_weight:               " << _log_weight                 << "\n" ;
            cout << "  _log_likelihood:           " << _log_likelihood             << "\n";
            cout << "  _forest._nleaves:          " << _forest._nleaves            << "\n";
            cout << "  _forest._ninternals:       " << _forest._ninternals         << "\n";
            cout << "  _forest._npatterns:        " << _forest._npatterns          << "\n";
            cout << "  _forest._nstates:          " << _forest._nstates            << "\n";
            cout << "  _forest._nsubtrees:        " << _forest._nsubtrees          << "\n";
            cout << "  _forest._last_edge_length: " << _forest._last_edge_length   << "\n";
            cout << "  newick description:        " << _forest.makeNewick(5,false) << "\n";
        }
    }

    inline double Particle::calcLogLikelihood() {
        //calculate likelihood for each gene tree
        double log_likelihood = 0.0;
        for (unsigned i=1; i<_forests.size(); i++) {
            double gene_tree_log_likelihood = _forests[i].calcLogLikelihood();
            assert(!isnan (log_likelihood));
//            cout << "gene tree log like: " << gene_tree_log_likelihood << endl;

            //total log likelihood is sum of gene tree log likelihoods?
            log_likelihood += gene_tree_log_likelihood;
        }
//        cout << "total log like: " << log_likelihood << endl;
        return log_likelihood;
    }

    inline double Particle::proposal() {
        //species tree
        tuple<string, string, string> t = _forests[0].speciesTreeProposal();
        
        //gene trees
        _forests[0].showForest();
        for (unsigned i=1; i<_forests.size(); i++){
            cout << "gene " << i << endl;
            _forests[i].geneTreeProposal(t, _forests[0]._last_edge_length);
        }
        
        double prev_log_likelihood = _log_likelihood;
        _log_likelihood = calcLogLikelihood();
        _log_weight = _log_likelihood - prev_log_likelihood;
        return _log_weight;
    }

    inline Particle::Particle(const Particle & other) {
        *this = other;
    }

    inline void Particle::saveForest(std::string treefilename) const {
        for (auto &_forest:_forests) {
            ofstream treef(treefilename);
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            treef << "  tree test = [&R] " << _forest.makeNewick(8, true)  << ";\n";
            treef << "end;\n";
            treef.close();
        }
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
        //species tree
            double sum_height = 0.0;
            for (auto nd : _forests[0]._preorder) {
                if (nd->getParent()!=_forests[0]._root) {
                    sum_height += nd->getEdgeLength();
                }
                if (!nd->getLeftChild()) {
                    break;
                }
            }
        return sum_height;
    }

    inline void Particle::setNumSubsets(unsigned n) {
        _nsubsets = n;
    }

    inline void Particle::mapSpecies(map<string, string> &taxon_map, vector<string> &species_names) {
        //species tree
        _forests[0].setUpSpeciesForest(species_names);
        
        //gene trees
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].setUpGeneForest(taxon_map);
        }
    }

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        
//        _forests.resize(other._forests.size());
//
////        _forests.clear();
//        unsigned i = 0;
//        for (auto forest : other._forests) {
////            _forests.push_back(forest);
//            _forests[i] = forest;
//            i++;
//        }

        _forests         = other._forests;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
    };
}
