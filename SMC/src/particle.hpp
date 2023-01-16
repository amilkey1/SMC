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
        void                                    saveMSCTrees(string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        void                                    setLogLikelihood(double log_likelihood);
        void                                    setParticleGeneration(int n);
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
        map<int, vector<double>>                     getRandomSeeds() {return _random_seeds;}
        void                                    setLogWeight(double w){_log_weight = w;}
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
        vector<Forest> &                                getForests() {return _forests;}
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

    private:

        static unsigned                         _nsubsets;
        vector<Forest>                          _forests;
        double                                  _log_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        int                                     _generation = 0;
        double                                  _gene_tree_marg_like;
        double                                  _prev_gene_tree_marg_like;
        map<int, vector<double>>                     _random_seeds;
        bool                                    _running_on_empty = false;
        tuple<string, string, string>           _t;
        unsigned                                _gene_tree_proposal_attempts;
        unsigned                                _coalescent_attempts_within_species_generation;
        unsigned                                _num_coalescent_attempts_needed;
        bool                                    _ready_to_join_species;
    
        bool                                    firstProposal();
        void                                    priorPostIshChoice(int i);
        void                                    resetVariables();
        bool                                    checkIfReadyToJoinSpecies();
};

    inline Particle::Particle() {
        //log weight and log likelihood are 0 for first generation
        _log_weight = 0.0;
        _log_likelihood = 0.0;
        _gene_tree_proposal_attempts = 0;
        _coalescent_attempts_within_species_generation = 0;
    };

    inline Particle::~Particle() {
//        cout << "destroying a particle" << endl;
    }
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
        double gene_tree_log_likelihood = 0.0;
        for (unsigned i=1; i<_forests.size(); i++) {
            gene_tree_log_likelihood = _forests[i].calcLogLikelihood();
            assert(!isnan (log_likelihood));
            assert(!isnan (gene_tree_log_likelihood));
            //total log likelihood is sum of gene tree log likelihoods
            log_likelihood += gene_tree_log_likelihood;
        }

        // set _generation for each forest
        for (int i=0; i < (int) _forests.size(); i++ ){
            _forests[i].setGeneration(_generation);
        }
        _generation++;
        return log_likelihood;
    }

    inline double Particle::proposal() {
        // TODO: check _generation = 0 for generation 0 on run on empty
        string event;
        showParticle();
        bool coalescence = false;
        
        // set starting variables
        if (_generation == 0 && _gene_tree_proposal_attempts == 0) {
            coalescence = firstProposal();
        }
//        showParticle();

        // attempt coalescent events until there has been at least 1 coalescence
        while (!coalescence) {
            showParticle();
            for (int i = 1; i < _forests.size(); i++) {
                if (_forests[0]._lineages.size() > 1) {
                checkIfReadyToJoinSpecies();
                    if (_coalescent_attempts_within_species_generation == _num_coalescent_attempts_needed || _ready_to_join_species) {
    //                if (_coalescent_attempts_within_species_generation == _num_coalescent_attempts_needed && _forests[i]._lineages.size() == 1) {
                        if (_forests[i]._lineages.size() > 1) {
                        _forests[i].extendGeneTreeLineages(_forests[0].getTreeHeight());
                            if (_forests[0]._lineages.size() > 1) {
                                _ready_to_join_species = true;
                                break;
                            }
                            else {
                                _ready_to_join_species = false;
                            }
                        }
                    }
                }
                else if (_forests[0]._lineages.size() == 1) { // if species tree is done
                    _forests[i].finishGeneTree();
                }
                }
//                if (!_ready_to_join_species && Forest::_proposal == "prior-prior") {
//                if (!_ready_to_join_species) {
//                    checkIfReadyToJoinSpecies();
//                    if (_ready_to_join_species) {
//                        _forests[i].extendGeneTreeLineages(_forests[0].getTreeHeight());
//                    }
//                }
//            }
            showParticle();
            
    //            event = _forests[0].chooseEvent();
                event = "speciation"; // for now, assume event is speciation
                if (event == "hybridization") {
                    hybridizationProposal();
                }
                else if (event == "speciation") {
                    // ready to join species
//                    if (_ready_to_join_species && _forests[0]._lineages.size() > 1) {
                    if (_ready_to_join_species) {
                        showParticle();
                        // reset _ready_to_join_species
                        _ready_to_join_species = false;
                        assert (_forests[0]._lineages.size() > 1); // TODO: figure out why it's going in here
//                        showParticle();
                        _t = _forests[0].speciesTreeProposal();
                        _coalescent_attempts_within_species_generation = 0; // reset coalescent attempts for this new species generation
                        _num_coalescent_attempts_needed = (unsigned) _forests[1]._lineages.size() - (unsigned) _forests[0]._lineages.size(); // TODO: check if this is always true
//                        showParticle();
                        
                        // if the species tree is not finished, add another species increment
                        if (_forests[0]._lineages.size()>1) {
                            _forests[0].addSpeciesIncrement();
                        }
                    }
                    
                    // don't join any more species, now must finish the gene trees first
                    for (unsigned i=1; i<_forests.size(); i++){
                        // if species have been joined previously
                        if (_forests[0]._lineages.size() != Forest::_nspecies) {
                               if (_forests[i]._lineages.size() > 1) {
                                   _gene_tree_proposal_attempts++;
                                   _coalescent_attempts_within_species_generation++;
//                                   showParticle();
                                    _forests[i].geneTreeProposal(_t, _forests[0]._last_edge_length, _forests[0].getTreeHeight());
                                   if (Forest::_proposal == "prior-post-ish") {
                                       priorPostIshChoice(i);
                                   }
//                                   showParticle();
                               }
                        }
                        // if no species have been joined at all
                        else {
                            for (unsigned i=1; i<_forests.size(); i++){ // TODO: fixed this
                                assert (_forests[i]._lineages.size() > 1);
                                _gene_tree_proposal_attempts++;
                                _coalescent_attempts_within_species_generation++;
                                _forests[i].firstGeneTreeProposal(_forests[0]._last_edge_length, _forests[0].getTreeHeight());
//                                showParticle();
                                if (Forest::_proposal == "prior-post-ish") {
                                    priorPostIshChoice(i);
                                }
                            }
                        }
                    }
                }
            for (int i=1; i<_forests.size(); i++) { // TODO: fixed this
                // if one forest has a deep coalescent event, need to continue in this generation
                coalescence = true;
                if (_forests[i]._num_coalescent_events_in_generation == 0) {
                    coalescence = false;
                    break;
                }
            }
        }
        
        if (Forest::_proposal != "prior-post-ish") {
        for (int i = 1; i<_forests.size(); i++) {
//            showParticle();
            if (_forests[i]._num_coalescent_events_in_generation > 1) {
                double smallest_branch = _forests[i].findShallowestCoalescence();
                _forests[i].revertToShallowest(smallest_branch);
//                showParticle();
                }
            }
        }
        
        if (_running_on_empty == false) {
            if (Forest::_proposal != "prior-post-ish") {
                double prev_log_likelihood = _log_likelihood;
                _log_likelihood = calcLogLikelihood();
                _log_weight = _log_likelihood - prev_log_likelihood;
            }
            else {
                // use the gene tree weights to calculate the particle weight
                _log_weight = 0.0;
                for (int i = 1; i<_forests.size(); i++) {
                    _log_weight += _forests[i]._gene_tree_log_weight;
                }
                _generation++;
            }
        }
        else {
            _generation++;
        }
        resetVariables();
        showParticle();
//        assert(_forests[0].getTreeHeight() >= _forests[1].getTreeHeight());
        return _log_weight;
    }

    inline bool Particle::firstProposal() {
        bool coalescence = true;
        // set starting variables
        assert (_generation == 0);
        assert (_gene_tree_proposal_attempts == 0);
        
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._theta = _forests[i]._starting_theta;
        }
        
        _ready_to_join_species = false;
//        _num_coalescent_attempts_needed = (unsigned) _forests[1]._lineages.size() - (unsigned) _forests[1]._species_partition.size();
        _num_coalescent_attempts_needed = Forest::_ntaxa - Forest::_nspecies;
            
        // choose a species tree height but don't join species yet (TODO: try both ways?)
        _forests[0].chooseSpeciesIncrement();
        
        for (unsigned i=1; i<_forests.size(); i++){
            assert (_forests[i]._lineages.size() > 1);
            _gene_tree_proposal_attempts++;
//            if (_proposal == "prior-post-ish") {
                _coalescent_attempts_within_species_generation++;
//            }
            _forests[i].firstGeneTreeProposal(_forests[0]._last_edge_length, _forests[0].getTreeHeight());
            if (Forest::_proposal == "prior-post-ish") {
                priorPostIshChoice(i);
            }
            showParticle();
        }
        
        // check if there has been a deep coalescent event
        for (unsigned i=1; i<_forests.size(); i++) {
            if (_forests[i]._num_coalescent_events_in_generation == 0) {
                coalescence = false;
                break;
            }
        }
        return coalescence;
    }

    inline bool Particle::checkIfReadyToJoinSpecies() {
        _ready_to_join_species = false;
        double species_tree_height = _forests[0].getTreeHeight();
        
        for (int i=1; i<_forests.size(); i++) {
            double gene_tree_height = _forests[i].getTreeHeight();
            if (gene_tree_height >= species_tree_height && _forests[0]._lineages.size() > 1) {
                _ready_to_join_species = true;
            }
            break;
        }
            
        return _ready_to_join_species;
    }

    inline void Particle::priorPostIshChoice(int i) {
        // attempt first gene tree proposal for all lineages, then select the one to keep
        _forests[i].chooseCoalescentEvent();
        _forests[i].mergeChosenPair(_forests[0]._last_edge_length);
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

    inline void Particle::saveMSCTrees(string treefilename) {
        // this function writes the species tree all particles to the same file
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

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _gene_tree_marg_like = other._gene_tree_marg_like;
        _prev_gene_tree_marg_like = other._prev_gene_tree_marg_like;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _t = other._t;
        _gene_tree_proposal_attempts = other._gene_tree_proposal_attempts;
        _coalescent_attempts_within_species_generation = other._coalescent_attempts_within_species_generation;
        _num_coalescent_attempts_needed = other._num_coalescent_attempts_needed;
        _ready_to_join_species = other._ready_to_join_species;
    };
}
