#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace boost;

#include "lot.hpp"

extern proj::Lot rng;

namespace proj {

class Particle {
    public:

        Particle();
        Particle(const Particle & other);
        ~Particle();
        typedef std::shared_ptr<Particle>               SharedPtr;


        void                                    debugParticle(std::string name);
        void                                    showParticle();
        double                                  proposal(bool gene_trees_only, bool unconstrained, bool deconstruct);
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
        void                                    resetGeneTreePartials(Data::SharedPtr d, map<string, string> taxon_map);
        void                                    mapSpecies(map<string, string> &taxon_map, vector<string> &species_names);
        void                                    mapGeneTrees(map<string, string> &taxon_map, vector<string> &species_names);
        void                                    saveForest(std::string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        void                                    setLogLikelihood(double log_likelihood);
        void                                    setLogCoalescentLikelihood(double log_coal_likelihood);
        void                                    setParticleGeneration(int n);
        double                                  calcHeight();
        void                                    buildEntireSpeciesTree();
        string                                  getSpeciesNewick() {return _forests[0].makeNewick(9, true);}
        double                                  getLogWeight(string a);

        void                                    setLogWeight(double w, string a);
        void                                    operator=(const Particle & other);
        const vector<Forest> &                  getForest() const {return _forests;}
        std::string                             saveForestNewick() {
            return _forests[0].makeNewick(8, true);
        }
        bool operator<(const Particle::SharedPtr & other) const {
            return _log_weight<other->_log_weight;
        }

        bool operator>(const Particle::SharedPtr & other) const {
            return _log_weight>other->_log_weight;
        }

        static void                                     setNumSubsets(unsigned n);
        void                                            showSpeciesIncrement();
        void                                            showSpeciesJoined();
        void                                            showSpeciesTree();
        void                                            showHybridNodes();
        string                                          saveHybridNodes();
        double                                          saveParticleWeights();
        double                                          saveParticleLikelihoods();
        void                                            showGamma();
        string                                          saveGamma();
        void                                            calculateGamma();
        void                                            setRunOnEmpty(bool a) {_running_on_empty = a;}
        vector<double>                                  getBranchLengths();
        vector<double>                                  getBranchLengthPriors();
        vector<double>                                  getGeneTreeLogLikelihoods();
        vector<double>                                  getTopologyPriors();
        void                                            hybridizationProposal();
        void                                            calcParticleWeight();
        vector<string>                                  getGeneTreeNames(int j);
        vector<string>                                  getGeneTreeNewicks();
        double                                          getNumSpecies(){return _forests[0]._lineages.size();}
        vector<pair<string, string>>                    getSpeciesJoined(){return _forests[1]._names_of_species_joined;}
        double                                          getCoalescentLikelihood(){return _log_coalescent_likelihood;}
        void                                            geneTreeProposal(bool unconstrained, bool deconstruct);
        void                                            speciesProposal();
        void                                            resetSpeciesInfo(){_t.clear();}
        void                                            resetSpeciesTreeHeight(){ _species_tree_height = 0.0;}
        void                                            resetSpecies();
        void                                            resetGeneIncrements();
        int                                             getNGenes(){return (int) _forests.size() - 1;}
        void                                            processSpeciesNewick(vector<string> newicks);
        void                                            processGeneNewicks(vector<string> newicks);
        void                                            remakeGeneTrees(map<string, string> &taxon_map);
        void                                            drawHeightsFromPrior();
    
    private:

        static unsigned                         _nsubsets;
        vector<Forest>                          _forests;
        double                                  _log_weight;
        double                                  _species_log_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        double                                  _log_coalescent_likelihood;
        int                                     _generation = 0;
        bool                                    _running_on_empty = false;
        bool                                    _no_species_joined = true;
        double                                  _species_tree_height;
        vector<pair<tuple<string, string, string>, double>> _t;
        void                                    resetVariables();
        bool                                    _inf = false;
};

    inline Particle::Particle() {
        //log weight and log likelihood are 0 for first generation
        _log_weight = 0.0;
        _log_likelihood = 0.0;
        _log_coalescent_likelihood = 0.0;
        _species_tree_height = 0.0;
        _species_log_weight = 0.0;
    };

    inline Particle::~Particle() {
//        cout << "destroying a particle" << endl;
    }
    inline void Particle::showParticle() {
//        cout << "log coalescent likelihood: " << _log_coalescent_likelihood << endl;
//        _forests[0].showForest();
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
        double gene_tree_log_likelihood = 0.0;
        for (unsigned i=1; i<_forests.size(); i++) {
            gene_tree_log_likelihood = _forests[i].calcLogLikelihood();
            assert(!isnan (log_likelihood));
            assert(!isnan (gene_tree_log_likelihood));
            //total log likelihood is sum of gene tree log likelihoods
            log_likelihood += gene_tree_log_likelihood;
            log_likelihood += _forests[i]._gene_tree_log_coalescent_likelihood;
        }

        // set _generation for each forest
        for (int i=0; i < (int) _forests.size(); i++ ){
            _forests[i].setGeneration(_generation);
        }
        _generation++;
        return log_likelihood;
    }

    inline void Particle::remakeGeneTrees(map<string, string> &taxon_map) {
        for (int i=1; i<_forests.size(); i++) {
            _forests[i].clear();
            _forests[i].remakeGeneTree(taxon_map);
        }
    }

    inline double Particle::proposal(bool gene_trees_only, bool unconstrained, bool deconstruct) {
        string event;
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._theta = _forests[i]._starting_theta;
        }

        if (_generation == 0) {
            for (int i=1; i<_forests.size(); i++) {
                _forests[i]._nincrements = 0;
                _forests[i]._gene_tree_log_coalescent_likelihood = 0.0;
            }
        }
        if (gene_trees_only) {
            geneTreeProposal(unconstrained, deconstruct);
        }
        else if (!gene_trees_only) {
            if (_generation == 0 || _generation == Forest::_ntaxa - 1) {
                for (int i=1; i<_forests.size(); i++) {
                    _forests[i].calcMinDepth();
                    _forests[i]._nincrements = 0;
                }
            }
            speciesProposal();
            _generation++;
        }
        
        if (_running_on_empty) {
            _generation++;
            _log_weight = 0.0;
        }
        resetVariables();
        return _log_weight;
    }

    inline void Particle::geneTreeProposal(bool unconstrained, bool deconstruct) {
        if (_generation == 0 && !unconstrained && deconstruct) {
            for (int i=1; i<_forests.size(); i++) {
                _forests[i].deconstructGeneTree();
            }
        }
        
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._theta = _forests[i]._starting_theta;
            assert (_t.size() == _forests[0]._nspecies);
            pair<double, string> species_info = _forests[i].chooseDelta(_t, false);
            _forests[i].geneTreeProposal(species_info, _t);
            if (_forests[i]._lineages.size() == 1) {
                _forests[i].refreshPreorder();
            }
        }

        if (!_running_on_empty) {
            double prev_log_likelihood = _log_likelihood;
            _log_likelihood = calcLogLikelihood();
            if (Forest::_proposal == "prior-prior") {
                _log_weight = _log_likelihood - prev_log_likelihood;
            }
            else {
                calcParticleWeight();
            }
        }
    }

    inline double Particle::getLogWeight(string a) {
        if (a == "s") {
            return _species_log_weight;
        }
        else {
            return _log_weight;
        }
        
    }

    inline void Particle::setLogWeight(double w, string a) {
        if (a == "s") {
            _species_log_weight = w;
        }
        else {
            _log_weight = w;
        }
        
    }


    inline void Particle::speciesProposal() {
        assert (!_inf);
        vector<double> species_tree_proposals;
        if (!_inf) {
            tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
            double prev_log_coalescent_likelihood = _log_coalescent_likelihood;
            vector<pair<int, int>> species_choices;
            
            if (_forests[0]._last_edge_length > 0.0) {
            // choose species to join if past the first species generation for each forest vector
                species_joined = _forests[0].speciesTreeProposal();
            }
                vector<double> max_depth_vector;
                double max_depth = 0.0;

                for (int i=1; i<_forests.size(); i++) {
                    
                    _forests[i].updateSpeciesPartitionTwo(species_joined);
                    string species1 = get<0>(species_joined);
                    string species2 = get<1>(species_joined);
                        
                    if (_forests[0]._lineages.size() > 1) {
                        _forests[i].resetDepthVector(species_joined);
                        max_depth = (_forests[i].getMinDepths())[0].first;
                        max_depth_vector.push_back(max_depth);
                    }
                }
                if (_forests[0]._lineages.size() > 1) {
                    max_depth = *min_element(max_depth_vector.begin(), max_depth_vector.end());
                    max_depth -= _forests[0].getTreeHeight();
                    // choose a species tree increment
                }
                
                if (_forests[0]._lineages.size() > 1) {
                    assert (max_depth > 0.0);
                    _forests[0].chooseSpeciesIncrement(max_depth);
                    species_tree_proposals.push_back(_forests[0]._last_edge_length);
                    _species_tree_height += _forests[0]._last_edge_length;
                }
                if (_forests[0]._lineages.size() == 1) {
                    _forests[0]._last_edge_length = 0.0;
                }
                        
                _t.push_back(make_pair(species_joined, _forests[0]._last_edge_length));
                
                for (int i = 1; i<_forests.size(); i++) {
                    double coal_like_increment = _forests[i].calcCoalescentLikelihood(_forests[0]._last_edge_length, species_joined, _species_tree_height, true);
                    _log_coalescent_likelihood += coal_like_increment;
                }
            
            if (!_running_on_empty) {
                _species_log_weight = _log_coalescent_likelihood - prev_log_coalescent_likelihood;
                double test = 1/_species_log_weight;
                if (test == -0) {
                    _inf = true;
                }
            }
        }
        }


    inline void Particle::calcParticleWeight() {
        // use the gene tree weights to calculate the particle weight
//        double prev_log_weight = _log_weight;
        _log_weight = 0.0;
        for (int i = 1; i<_forests.size(); i++) {
            _log_weight += _forests[i]._gene_tree_log_weight;
//            _log_weight += _forests[i]._gene_tree_log_likelihood;
        }
//        _log_weight = prev_log_weight - _log_weight;
//        _generation++;
    }


    inline void Particle::resetVariables() {
        for (int i = 1; i<_forests.size(); i++) {
            _forests[i]._num_coalescent_events_in_generation = 0;
            _forests[i]._searchable_branch_lengths.clear();
            _forests[i]._new_nodes.clear();
        }
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
        ofstream treef(treefilename);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        treef << "  tree test = [&R] " << _forests[0].makeNewick(8, true)  << ";\n";
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

    inline void Particle::processSpeciesNewick(vector<string> newicks) {
        if (newicks.size() > 1) {
            throw XProj(boost::str(boost::format("only one species tree newick may be specified")));
        }
        
        if (newicks.size() > 0) {
            string newick = newicks[0];
            _t = _forests[0].buildFromNewickTopology(newick, true, false);
            drawHeightsFromPrior();
        }
        
        else {
            buildEntireSpeciesTree();
        }
    }

    inline void Particle::drawHeightsFromPrior() {
        int nlineages = Forest::_nspecies;
        _forests[0].resetIncrements();
        
       vector<string> existing_lineages = _forests[0].setUpExistingLineagesVector();
        
        for (int a=0; a<Forest::_nspecies-1; a++) {
            existing_lineages = _forests[0].updateExistingLineagesVector(existing_lineages, _t[a].first);
            _forests[0].chooseSpeciesIncrementFromNewick(nlineages, existing_lineages);
            _t[a].second = _forests[0]._last_edge_length;
            
            nlineages--;
        }
        
        assert (_t[Forest::_nspecies-1].second == 0); // last element of _t should not draw a branch length
    }

    inline void Particle::processGeneNewicks(vector<string> newicks) {
        assert (newicks.size() == _forests.size() - 1);
        for (int i=1; i<_forests.size(); i++) {
            _forests[i].buildFromNewick(newicks[i-1], true, false);
            _forests[i].refreshPreorder();
        }
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

    inline void Particle::hybridizationProposal() {
        vector<string> hybridized_nodes = _forests[0].hybridizeSpecies();
        if (_forests[0]._lineages.size()>1) {
            _forests[0].addSpeciesIncrement();
        }
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].hybridizeGene(hybridized_nodes, _forests[0]._last_edge_length);
        }
        calculateGamma();
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

    inline void Particle::mapGeneTrees(map<string, string> &taxon_map, vector<string> &species_names) {
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
            if (nd._major_parent) {
                string gammastr = to_string(_forests[0]._gamma[i]);
                nodes +=  "hybridized node is: " + nd._name + " with minor parent " + nd._minor_parent->_name + " and major parent " + nd._major_parent->_name + "\n" + "gamma is: " + gammastr + "\n";
                i++;
            }
        }
        return nodes;
    }

    inline double Particle::saveParticleWeights() {
        return _log_weight;
    }

    inline double Particle::saveParticleLikelihoods() {
        return _log_likelihood;
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

    inline void Particle::setLogLikelihood(double log_likelihood) {
        _log_likelihood = log_likelihood;
    }

    inline void Particle::setLogCoalescentLikelihood(double log_coal_likelihood) {
        _log_coalescent_likelihood = log_coal_likelihood;
    }

    inline void Particle::setParticleGeneration(int n) {
        _generation = n;
    }

    inline vector<double> Particle::getBranchLengths() {
        vector<double> divergence_times;
        for (auto &f:_forests) {
            for (auto &b:f._increments) {
                divergence_times.push_back(b.first);
            }
        }
        return divergence_times;
    }

    inline vector<double> Particle::getBranchLengthPriors() {
        vector<double> priors;
        for (auto &f:_forests) {
            for (auto &b:f._increments) {
                priors.push_back(b.second);
            }
        }
        return priors;
    }


    inline vector<double> Particle::getGeneTreeLogLikelihoods() {
        vector<double> gene_tree_log_likelihoods;
        for (int i=1; i<_forests.size(); i++) {
            gene_tree_log_likelihoods.push_back(_forests[i]._gene_tree_log_likelihood);
        }
        return gene_tree_log_likelihoods;
    }

    inline vector<double> Particle::getTopologyPriors() {
        // calculate species tree topology probability
        vector<double> topology_priors;
        for (auto &f:_forests) {
            topology_priors.push_back(f._log_joining_prob);
        }
        return topology_priors;
    }

    inline vector<string> Particle::getGeneTreeNames(int j) {
        vector<string> names;
        string particle_number = to_string(j);
        for (int i=1; i<_forests.size(); i++) {
            string forest_number = to_string(i-1);
            string name = "\"gene-" + forest_number + "-" + particle_number + "\"";
            names.push_back(name);
        }
        return names;
    }

    inline vector<string> Particle::getGeneTreeNewicks() {
        vector<string> newicks;
        for (int i=1; i<_forests.size(); i++) {
//            string newick = "newick:\"" + _forests[i].makeNewick(9, true) + "\"";
            string newick = _forests[i].makeNewick(9, true);
            newicks.push_back(newick);
        }
        return newicks;
    }

    inline void Particle::resetSpecies() {
        setLogLikelihood(0.0);
        setLogWeight(0.0, "s");
        resetSpeciesInfo();
        resetSpeciesTreeHeight();
        _forests[0]._increments.clear();
        _log_coalescent_likelihood = 0.0;
    }

    inline void Particle::resetGeneIncrements() {
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._increments.clear();
        }
    }

    inline void Particle::buildEntireSpeciesTree() {
        double max_depth = 0.0;
        _forests[0].chooseSpeciesIncrement(max_depth);
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        _t.push_back(make_pair(species_joined, 0.0));

//        cout << "new particle" << endl;
        for (int i=0; i < _forests[0]._nspecies-1; i++) {
            if (_forests[0]._lineages.size() > 1) {
                species_joined = _forests[0].speciesTreeProposal();
                
                // if the species tree is not finished, add another species increment
                if (_forests[0]._lineages.size()>1) {
                    _forests[0].addSpeciesIncrement();
                }
                
                double edge_len = 0.0;
                if (_forests[0]._lineages.size() > 1) {
                    edge_len = _forests[0]._last_edge_length;
                }
                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline void Particle::resetGeneTreePartials(Data::SharedPtr d, map<string, string> taxon_map) {
        _nsubsets = d->getNumSubsets();
        _data = d;
        for (int i=1; i<_forests.size(); i++) {
            _forests[i].setData(d, i, taxon_map);
        }
    }


    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _t = other._t;
        _running_on_empty = other._running_on_empty;
        _no_species_joined = other._no_species_joined;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _species_tree_height = other._species_tree_height;
        _inf = other._inf;
        _species_log_weight = other._species_log_weight;
    };
}
