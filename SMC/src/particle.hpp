#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
using namespace std;
using namespace boost;

#include "lot.hpp"
#include "conditionals.hpp"

extern proj::Lot rng;

namespace proj {

class Particle {
    public:

        Particle();
        Particle(const Particle & other);
        typedef std::shared_ptr<Particle>               SharedPtr;

        void                                    debugParticle(std::string name);
        void                                    showParticle();
        void                                    proposal();
        void                                    altProposal();
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
        double                                  getLogLikelihood();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    operator=(const Particle & other);
        const vector<Forest> &                  getForest() const {return _forests;}
        string                                  saveForestNewick() {
            return _forests[0].makeNewick(8, true);}
            
            string                             saveGeneNewick(unsigned i) {
            return _forests[i].makeNewick(8, true);}
    
        bool operator<(const Particle::SharedPtr & other) const {
            return _log_weight<other->_log_weight;
        }

        bool operator>(const Particle::SharedPtr & other) const {
            return _log_weight>other->_log_weight;
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
        int                                             selectPair(vector<double> weight_vec);
        double                                          getTopologyPrior(unsigned i);
        vector<pair<double, double>>                    getIncrementPriors(unsigned i);
        double                                          getCoalescentLikelihood(unsigned g);
        bool                                            speciesJoinProposed();

    private:

        static unsigned                         _nsubsets;
        vector<Forest>                          _forests;
        double                                  _log_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        int                                     _generation = 0;
        unsigned                                _prev_forest_number;
        bool                                    _species_join_proposed;
        double                                  _prev_log_coal_like;
        double                                  _prev_increment;
};

    inline Particle::Particle() {
        //log weight and log likelihood are 0 for first generation
        _log_weight = 0.0;
        _log_likelihood = 0.0;
        _prev_forest_number = -1;
        _species_join_proposed = false;
        _prev_log_coal_like = 0.0;
        _prev_increment = 0.0;
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
            if (_generation == 0) {
                _forests[i]._combined_likelihood = gene_tree_log_likelihood; // coalescent likelihood starts at 0
            }
        }
        if (_generation == 0) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline double Particle::getLogLikelihood() {
        //retrieve likelihood for each gene tree
        double log_likelihood = 0.0;
        for (unsigned i=1; i<_forests.size(); i++) {
            double gene_tree_log_likelihood = _forests[i]._gene_tree_log_likelihood;
            assert(!isnan (log_likelihood));
            //total log likelihood is sum of gene tree log likelihoods
            log_likelihood += gene_tree_log_likelihood;
        }
        if (_generation == 0) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline void Particle::proposal() {
#if defined (JOIN_FIRST)
        altProposal();
#else
        _species_join_proposed = false;
        bool done = false;
        
        while (!done) {
#if !defined(TWO_JOINS)
            done = true;
#endif
    
        vector<double> forest_rates; // this vector contains total rate of species tree, gene 1, etc.
        vector<vector<double>> gene_forest_rates; // this vector contains rates by species for each gene forest
        gene_forest_rates.resize(_forests.size()-1);
        vector<unsigned> event_choice_index;
        vector<string> event_choice_name;
            
        for (int i=0; i<_forests.size(); i++) {
            if (i > 0) {
#if defined (SNAKE)
                if (i == 1) {
                    Forest::_theta *= 3.25926;
                }
                else if (i == 2) {
                    Forest::_theta *= 0.64160;
                }
                else if (i == 3) {
                    Forest::_theta *= 0.75517;
                }
                else if (i == 4) {
                    Forest::_theta *= 0.98977;
                }
                else if (i == 5) {
                    Forest::_theta *= 0.73621;
                }
                else if (i == 6) {
                    Forest::_theta *= 1.45336;
                }
                else if (i == 7) {
                    Forest::_theta *= 0.51891;
                }
                else if (i == 8) {
                    Forest::_theta *= 1.77001;
                }
                else if (i == 9) {
                    Forest::_theta *= 0.57665;
                }
                else if (i == 10) {
                    Forest::_theta *= 0.54812;
                }
                else if (i == 11) {
                    Forest::_theta *= 0.60184;
                }
                else if (i == 12) {
                    Forest::_theta *= 0.54980;
                }
                else if (i == 13) {
                    Forest::_theta *= 0.39104;
                }
                else if (i == 14) {
                    Forest::_theta *= 0.45743;
                }
                else if (i == 15) {
                    Forest::_theta *= 0.52157;
                }
#endif
                vector<pair<double, string>> rates_by_species = _forests[i].calcForestRate();
                double total_gene_rate = 0.0;
                for (auto &r:rates_by_species) {
                    gene_forest_rates[i-1].push_back(r.first);
                    event_choice_name.push_back(r.second);
                    total_gene_rate += r.first;
                    event_choice_index.push_back(i);
                }
                if (total_gene_rate > 0.0) {
                    forest_rates.push_back(total_gene_rate);
                }
#if defined (SNAKE)
                Forest::_theta = 0.001;
#endif
            }
            else {
                if (_forests[0]._lineages.size() > 1) {
                    forest_rates.push_back(Forest::_speciation_rate * _forests[0]._lineages.size());
                    event_choice_index.push_back(0);
                    event_choice_name.push_back("species");
                }
            }
        }
        
        double total_rate = 0.0;
        for (auto &r:forest_rates) {
            total_rate += r;
        }
        
        // draw an increment
        double increment = rng.gamma(1.0, 1.0/(total_rate));
            _prev_increment = increment;
            
        // add increment to all nodes in all forests
        for (int i=0; i<_forests.size(); i++) {
            if (_forests[i]._lineages.size() > 1) {
                _forests[i].addIncrement(increment); // if forest is finished, don't add another increment
            }
            else {
                _forests[i]._done = true;
            }
        }
            
        vector<double> event_choice_rates;
        if (_forests[0]._lineages.size() > 1) {
            event_choice_rates.push_back(forest_rates[0]); // push back species tree rate
        }
        for (int i=0; i<gene_forest_rates.size(); i++) {
            for (auto &r:gene_forest_rates[i]) {
                event_choice_rates.push_back(r);
            }
        }
            
        // choose an event
        for (auto &p:event_choice_rates) {
             p = log(p/total_rate);
         }
        unsigned index = selectPair(event_choice_rates);
        unsigned forest_number = event_choice_index[index];
        string species_name = event_choice_name[index];
            
            // need to calculate coalescent likelihood before joining anything or updating species partition
        for (int f=0; f<_forests.size(); f++) {
            bool new_increment = false;
            bool coalescence = false;
            bool gene_tree = false;
            
            if (f == _prev_forest_number) {
                // add to existing increment + prior
                new_increment = true;
            }
            if (f == forest_number) {
                coalescence = true;
            }
            if (_generation == 0) {
                new_increment = true;
            }
            if (f > 0) {
                gene_tree = true;
            }
            
            _forests[f].calcIncrementPrior(increment, species_name, new_increment, coalescence, gene_tree);
        }
            
        if (species_name == "species") {
            // species tree proposal, need to update species partition in all gene forests
            assert (index == 0);
            assert (forest_number == 0);
            tuple <string, string, string> species_joined = _forests[0].speciesTreeProposal(increment);
            for (int i=1; i<_forests.size(); i++) {
                // reset species partitions for all gene forests
                _forests[i].updateSpeciesPartition(species_joined);
            }
        }
        
        else {
            _forests[forest_number].allowCoalescence(species_name, increment);
            done = true;
        }
                    
        if (species_name != "species") {
            _log_weight = _forests[forest_number]._log_weight;
//            cout << "log weight is " << _log_weight << endl;
        }
        
        else {
            // species log weight is always 0
            _species_join_proposed = true;
            _log_weight = 0.0;
        }
            
        _prev_forest_number = forest_number;
        }
        _generation++;
    #endif
    }

    inline void Particle::altProposal() {
//        showParticle();
        _species_join_proposed = false;
        
        vector<double> forest_rates; // this vector contains total rate of species tree, gene 1, etc.
        vector<vector<double>> gene_forest_rates; // this vector contains rates by species for each gene forest
        gene_forest_rates.resize(_forests.size()-1);
        vector<unsigned> event_choice_index;
        vector<string> event_choice_name;
        
        if (_generation == 0) {
            // choose an increment
            for (int i=0; i<_forests.size(); i++) {
                if (i > 0) {
                    vector<pair<double, string>> rates_by_species = _forests[i].calcForestRate();
                    double total_gene_rate = 0.0;
                    for (auto &r:rates_by_species) {
                        gene_forest_rates[i-1].push_back(r.first);
                        event_choice_name.push_back(r.second);
                        total_gene_rate += r.first;
                        event_choice_index.push_back(i);
                    }
                    if (total_gene_rate > 0.0) {
                        forest_rates.push_back(total_gene_rate);
                    }
                }
                else {
                    if (_forests[0]._lineages.size() > 1) {
                        forest_rates.push_back(Forest::_speciation_rate * _forests[0]._lineages.size());
                        event_choice_index.push_back(0);
                        event_choice_name.push_back("species");
                    }
                }
            }

            double total_rate = 0.0;
            for (auto &r:forest_rates) {
                total_rate += r;
            }
            
            // draw an increment
            double increment = rng.gamma(1.0, 1.0/(total_rate));
            _prev_increment = increment;
            
            // add increment to all nodes in all forests
            for (int i=0; i<_forests.size(); i++) {
                if (_forests[i]._lineages.size() > 1) {
                    _forests[i].addIncrement(increment); // if forest is finished, don't add another increment
                }
                else {
                    _forests[i]._done = true;
                }
            }
//            _log_weight = 0.0;
            for (int f=0; f<_forests.size(); f++) {
                bool gene_tree = false;
                if (f > 0) {
                    gene_tree = true;
                }
                _forests[f].calcIncrementPrior(increment, "null", true, false, gene_tree);
            }
        }
            
        else {
            // join taxa
            for (int i=0; i<_forests.size(); i++) {
                if (i > 0) {
                    vector<pair<double, string>> rates_by_species = _forests[i].calcForestRate();
                    double total_gene_rate = 0.0;
                    for (auto &r:rates_by_species) {
                        gene_forest_rates[i-1].push_back(r.first);
                        event_choice_name.push_back(r.second);
                        total_gene_rate += r.first;
                        event_choice_index.push_back(i);
                    }
                    if (total_gene_rate > 0.0) {
                        forest_rates.push_back(total_gene_rate);
                    }
                }
                else {
                    if (_forests[0]._lineages.size() > 1) {
                        forest_rates.push_back(Forest::_speciation_rate * _forests[0]._lineages.size());
                        event_choice_index.push_back(0);
                        event_choice_name.push_back("species");
                    }
                }
            }

            double total_rate = 0.0;
            for (auto &r:forest_rates) {
                total_rate += r;
            }
            
            vector<double> event_choice_rates;
            if (_forests[0]._lineages.size() > 1) {
                event_choice_rates.push_back(forest_rates[0]); // push back species tree rate
            }
            for (int i=0; i<gene_forest_rates.size(); i++) {
                for (auto &r:gene_forest_rates[i]) {
                    event_choice_rates.push_back(r);
                }
            }
            
            // choose an event
            for (auto &p:event_choice_rates) {
                 p = log(p/total_rate);
             }
            unsigned index = selectPair(event_choice_rates);
            unsigned forest_number = event_choice_index[index];
            string species_name = event_choice_name[index];
            
            // calculate coalescent likelihood before actually joining
            for (int f=0; f<_forests.size(); f++) {
                bool gene_tree = false;
                
                if (f == forest_number) { // TODO: double check new_increment for lorad output - doesn't matter for coalescent likelihood
                    if (f > 0) {
                        gene_tree = true;
                    }
                    
                    _forests[f].calcIncrementPriorForJoinOnly(_prev_increment, species_name, gene_tree);
                    break;
                }
            }
            
            if (species_name == "species") {
                // species tree proposal, need to update species partition in all gene forests
                assert (index == 0);
                assert (forest_number == 0);
                tuple <string, string, string> species_joined = _forests[0].speciesTreeProposal(_prev_increment);
                for (int i=1; i<_forests.size(); i++) {
                    // reset species partitions for all gene forests
                    _forests[i].updateSpeciesPartition(species_joined);
                }
            }
            else {
                _forests[forest_number].allowCoalescence(species_name, _prev_increment);
            }
            
            forest_rates.clear();
            gene_forest_rates.clear();
            gene_forest_rates.resize(_forests.size()-1);
            
            // draw a new increment and add to all nodes
            for (int i=0; i<_forests.size(); i++) {
                if (i > 0) {
                    vector<pair<double, string>> rates_by_species = _forests[i].calcForestRate();
                    double total_gene_rate = 0.0;
                    for (auto &r:rates_by_species) {
                        gene_forest_rates[i-1].push_back(r.first);
                        event_choice_name.push_back(r.second);
                        total_gene_rate += r.first;
                        event_choice_index.push_back(i);
                    }
                    if (total_gene_rate > 0.0) {
                        forest_rates.push_back(total_gene_rate);
                    }
                }
                else {
                    if (_forests[0]._lineages.size() > 1) {
                        forest_rates.push_back(Forest::_speciation_rate * _forests[0]._lineages.size());
                        event_choice_index.push_back(0);
                        event_choice_name.push_back("species");
                    }
                }
            }

            double new_total_rate = 0.0;
            for (auto &r:forest_rates) {
                new_total_rate += r;
            }
            
            // draw an increment
            double increment = rng.gamma(1.0, 1.0/(new_total_rate));
            _prev_increment = increment;
            
            // add increment to all nodes in all forests
            for (int i=0; i<_forests.size(); i++) {
                if (_forests[i]._lineages.size() > 1) {
                    _forests[i].addIncrement(increment); // if forest is finished, don't add another increment
                }
                else {
                    _forests[i]._done = true;
                }
            }
                        
            // calculate coalescent likelihood after extending increment - no join
            
            for (int f=0; f<_forests.size(); f++) {
                bool new_increment = false;
                bool coalescence = false;
                bool gene_tree = false;
                
                if (f == _prev_forest_number) {
                    // add to existing increment + prior
                    new_increment = true;
                }
                
                // bool coalescence is always false now
                
                if (_generation == 1) {
                    new_increment = false;
                }
                if (f > 0) {
                    gene_tree = true;
                }
                
                _forests[f].calcIncrementPrior(increment, species_name, new_increment, coalescence, gene_tree);
            }
            
            if (species_name != "species") {
                _log_weight = _forests[forest_number]._log_weight;
            }
            else {
                
                // species log weight is always 0
                _species_join_proposed = true;
                _log_weight = 0.0;
                
            }
            _prev_forest_number = forest_number;
        }
        _generation++;
//        showParticle();
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

    inline int Particle::selectPair(vector<double> weight_vec) {
        // choose a random number [0,1]
        double u = rng.uniform();
        double cum_prob = 0.0;
        unsigned index = 0;
        for (int i=0; i < (int) weight_vec.size(); i++) {
            cum_prob += exp(weight_vec[i]);
            if (u <= cum_prob) {
                index = i;
                break;
            }
        }
        // return index of choice
        return index;
    }

    inline double Particle::getTopologyPrior(unsigned i) {
        return _forests[i]._log_joining_prob;
    }

    inline vector<pair<double, double>> Particle::getIncrementPriors(unsigned i) {
        return _forests[i]._increments_and_priors;
    }

    inline double Particle::getCoalescentLikelihood(unsigned g) {
        assert (g>0); // no coalescent likelihood for species tree
        return _forests[g]._log_coalescent_likelihood;
    }

    inline bool Particle::speciesJoinProposed() {
        if (_species_join_proposed) {
            return true;
        }
        else {
            return false;
        }
    }

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _prev_forest_number = other._prev_forest_number;
        _species_join_proposed = other._species_join_proposed;
        _prev_log_coal_like = other._prev_log_coal_like;
        _prev_increment = other._prev_increment;
    };
}

