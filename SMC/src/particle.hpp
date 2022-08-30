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
        void                                    setData(Data::SharedPtr d, map<string, string> &taxon_map) {
                                                    _nsubsets = d->getNumSubsets();
                                                    _data = d;
                                                    int index = 0;
                                                    _forests.resize(_nsubsets+1);
                                                    for (auto &_forest:_forests) {
                                                        _forest.setData(d, index, taxon_map);
                                                        index++;
                                                    }
                                                }
        void                                    mapSpecies(map<string, string> &taxon_map, vector<string> &species_names);
        void                                    saveForest(std::string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    operator=(const Particle & other);
        const vector<Forest> &                  getForest() const {return _forests;}
        std::string                             saveForestNewick() {
            return _forests[0].makeNewick(8, true);
        }
        bool operator<(const Particle & other) const {
            return _log_weight<other._log_weight;
        }

        bool operator>(const Particle & other) const {
            return _log_weight>other._log_weight;
        }

        static void                                     setNumSubsets(unsigned n);
        vector<Forest> &                                getForests() {return _forests;}
        void                                            showSpeciesIncrement();
        void                                            showSpeciesJoined();
        void                                            showSpeciesTree();
        void                                            showHybridNodes();
        string                                          saveHybridNodes();
        void                                            showGamma();
        string                                          saveGamma();
        void                                            calculateGamma();

    private:

        static unsigned                         _nsubsets;
        vector<Forest>                          _forests;
        double                                  _log_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        int                                     _generation = 0;
    vector<tuple<string, string, string>> _triple;
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

    inline void Particle::showSpeciesTree() {
        //print out weight of each particle
        cout << "\nParticle:\n";
        cout << "  _log_weight: " << _log_weight << "\n" ;
        cout << " _log_likelihood: " << _log_likelihood << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _forests[0].showForest();
    }

    //more detailed version of showParticle
    inline void Particle::debugParticle(std::string name) {
        cout << "debugging particle" << endl;
        //print out weight of each particle
        cout << "\nParticle " << name << ":\n";
        for (auto &_forest:_forests) {
            cout << "  _log_weight:               " << _log_weight                 << "\n" ;
            cout << "  _log_likelihood:           " << _log_likelihood             << "\n";
            cout << "  _forest._nleaves:          " << _forest._nleaves            << "\n";
            cout << "  _forest._ninternals:       " << _forest._ninternals         << "\n";
            cout << "  _forest._npatterns:        " << _forest._npatterns          << "\n";
            cout << "  _forest._nstates:          " << _forest._nstates            << "\n";
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
            //total log likelihood is sum of gene tree log likelihoods
            log_likelihood += gene_tree_log_likelihood;
        }
        _generation++;
        
        // set _generation for each forest
        for (int i=0; i < (int) _forests.size(); i++ ){
            _forests[i].setGeneration(_generation);
        }
        return log_likelihood;
    }

    inline double Particle::proposal() {
        string event;
        if (_generation == 0) {
            _forests[0].chooseSpeciesIncrement();
            for (unsigned i=1; i<_forests.size(); i++){
                _forests[i].firstGeneTreeProposal(_forests[0]._last_edge_length);
            }
        }
        else if (_forests[0]._lineages.size()==1) {
            for (unsigned i=1; i<_forests.size(); i++) {
                list<Node*> lineages_list(_forests[i]._lineages.begin(), _forests[i]._lineages.end());
                _forests[i].fullyCoalesceGeneTree(lineages_list);
            }
        }
        else {
            event = _forests[0].chooseEvent();
            if (event == "hybridization") {
                vector<string> hybridized_nodes = _forests[0].hybridizeSpecies();
                if (_forests[0]._lineages.size()>1) {
                    _forests[0].addSpeciesIncrement();
                }
                for (unsigned i=1; i<_forests.size(); i++) {
                    _forests[i].hybridizeGene(hybridized_nodes, _forests[0]._last_edge_length);
                }
                calculateGamma();
            }
            else if (event == "speciation") {
                tuple<string, string, string> t = _forests[0].speciesTreeProposal();
                if (_forests[0]._lineages.size()>1) {
                    _forests[0].addSpeciesIncrement();
                }
                for (unsigned i=1; i<_forests.size(); i++){
                    _forests[i].geneTreeProposal(t, _forests[0]._last_edge_length);
                }
            }
        }
        double prev_log_likelihood = _log_likelihood;
        _log_likelihood = calcLogLikelihood();
        _log_weight = _log_likelihood - prev_log_likelihood;
        
        return _log_weight;
    }

    inline Particle::Particle(const Particle & other) {
        *this = other;
    }

    inline void Particle::calculateGamma() {
        double major = 0.0;
        double total = _forests.size()-1;
        for (int i=1; i < (int) _forests.size(); i++) {
            if (_forests[i]._last_direction == "major") {
                major++;
            }
        }
        double gamma = major / total;
        _forests[0]._gamma.push_back(gamma);
    }

    inline void Particle::saveForest(std::string treefilename)  {
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

        // calculate height of lineage
        Node* base_node = _forests[0]._lineages[0];
        sum_height += base_node->getEdgeLength();
        for (Node* child=base_node->_left_child; child; child=child->_left_child) {
            sum_height += child->getEdgeLength();
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

    inline void Particle::showSpeciesIncrement(){
        cout << "species tree increment: " << "     " << _forests[0]._last_edge_length << endl;
    }

    inline void Particle::showSpeciesJoined(){
        _forests[0].showSpeciesJoined();
    }
    
    inline void Particle::showHybridNodes() {
        cout << "particle" << endl;
        showGamma();
        for (auto &nd:_forests[0]._nodes) {
            if (nd._major_parent) {
                cout << "       " << "hybridized node is: " << nd._name << " with minor parent " << nd._minor_parent->_name << " and major parent " << nd._major_parent->_name << endl;
            }
        }
    }

    inline string Particle::saveHybridNodes() {
        string nodes = "";
        int i = 0;
        for (auto &nd:_forests[0]._nodes) {
//        for (Node* nd = _forests[0]._lineages[0]; nd; nd=nd->_left_child->_right_sib) {
//            for (Node * child=new_nd->_left_child; child; child=child->_right_sib) {
            if (nd._major_parent) {
                string gammastr = to_string(_forests[0]._gamma[i]);
                nodes +=  "hybridized node is: " + nd._name + " with minor parent " + nd._minor_parent->_name + " and major parent " + nd._major_parent->_name + "\n" + "gamma is: " + gammastr + "\n";
                i++;
            }
        }
        return nodes;
    }

    inline void Particle::showGamma() {
        if (_forests[0]._gamma.size() > 0) {
            cout << "   " << "gamma is: " << endl;
            for (auto &g:_forests[0]._gamma) {
                cout << g << "   ";
            }
            cout << "\n";
        }
    }

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _triple         = other._triple;
    };
}
