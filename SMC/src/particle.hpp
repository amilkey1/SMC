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
        double                                  proposal(bool gene_trees_only, bool unconstrained);
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
        void                                    mapGeneTrees(map<string, string> &taxon_map, vector<string> &species_names);
        void                                    saveForest(std::string treefilename);
        void                                    saveMSCTrees(string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        void                                    setLogLikelihood(double log_likelihood);
        void                                    setLogCoalescentLikelihood(double log_coal_likelihood);
        void                                    setParticleGeneration(int n);
        double                                  calcHeight();
        void                                    buildEntireSpeciesTree();
        string                                  getSpeciesNewick() {return _forests[0].makeNewick(9, true);}
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
        void                                            setBuildEntireSpeciesTree(bool a) {_species_first = a;}
        vector<double>                                  getBranchLengths();
        vector<double>                                  getBranchLengthPriors();
        vector<double>                                  getGeneTreeLogLikelihoods();
        vector<double>                                  getTopologyPriors();
        void                                            hybridizationProposal();
        bool                                            checkForDeepCoalescence();
        void                                            priorPriorProposal();
        void                                            calcParticleWeight();
        void                                            speciesJoinedProposal();
        vector<string>                                  getGeneTreeNames(int j);
        vector<string>                                  getGeneTreeNewicks();
        double                                          getNumSpecies(){return _forests[0]._lineages.size();}
        vector<pair<string, string>>                    getSpeciesJoined(){return _forests[1]._names_of_species_joined;}
        void                                            buildFakeSpeciesTree();
        void                                            panmicticSpeciesPartition();
        double                                          getCoalescentLikelihood(){return _log_coalescent_likelihood;}
        void                                            firstGen();
        void                                            geneTreeProposal(bool unconstrained);
        void                                            speciesProposal();
        void                                            filterSpeciesTrees();
        void                                            normalizeSpeciesTreeWeights();
        void                                            resampleParticles();
        void                                            resetSpeciesInfo(){_t.clear();}
        void                                            resetSpeciesTreeHeight(){ for (auto &s:_species_tree_height) {s = 0.0;};}
        void                                            resetSpecies();
        void                                            resetGeneIncrements();
        void                                            increaseNumberOfForests();
        double                                          getRunningSum();

    private:

        static unsigned                         _nsubsets;
        vector<Forest>                          _forests;
        vector<vector<Forest>>                  _alt_forests;
        double                                  _log_weight;
        vector<double>                          _log_weight_options;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        double                                  _log_coalescent_likelihood;
        vector<double>                          _log_coalescent_likelihood_options;
        int                                     _generation = 0;
        map<int, vector<double>>                     _random_seeds;
        bool                                    _running_on_empty = false;
        bool                                    _no_species_joined = true;
        vector<double>                                  _species_tree_height;
        vector<vector<pair<tuple<string, string, string>, double>>> _t;
        unsigned                                _gene_tree_proposal_attempts;
        bool                                    _ready_to_join_species;
        bool                                    _species_first;
        bool                                    firstProposal();
        void                                    priorPostIshChoice(int i, vector<pair<tuple<string, string, string>, double>> _t);
        void                                    resetVariables();
        bool                                    checkIfReadyToJoinSpecies();
        int                                     _nspecies_forests = 1000;
        bool                                    _inf = false;
//        bool                                    _reset_min = false;
        vector<bool>                             _reset_min;
};

    inline Particle::Particle() {
        //log weight and log likelihood are 0 for first generation
        _log_weight = 0.0;
        _log_likelihood = 0.0;
        _gene_tree_proposal_attempts = 0;
        _log_coalescent_likelihood = 0.0;
        _species_tree_height.clear();
//        _reset_min = false;
    };

    inline Particle::~Particle() {
//        cout << "destroying a particle" << endl;
//	cout << "test" << endl;
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

    inline double Particle::proposal(bool gene_trees_only, bool unconstrained) {
        string event;
        
        if (unconstrained && _generation == 0) {
            firstGen();
        }
//        _forests[0].showForest();
        if (gene_trees_only) {
            geneTreeProposal(unconstrained);
//            _forests[1].showForest();
        }
        else if (!gene_trees_only) {
            if (_generation == Forest::_ntaxa - 1) {
                increaseNumberOfForests();
                for (int a = 0; a<_alt_forests.size(); a++) {
                    for (int i=1; i<_alt_forests[a].size(); i++) {
                        _alt_forests[a][i].calcMinDepth();
                    }
                }
            }
            speciesProposal();
            filterSpeciesTrees();
            
            if (_generation == Forest::_ntaxa+Forest::_nspecies-2) {
                // once species trees are fully resolved, choose a species tree based on their weights
                int choice = _forests[0].selectPair(_log_weight_options);
                // reset _forests to the chosen species tree in _alt_forests
                _forests[0] = _alt_forests[choice][0];
                _log_weight = _log_weight_options[choice];
                _alt_forests.clear();
                vector<pair<tuple<string, string, string>, double>> chosen_species = _t[choice];
                _t.clear();
                _t.push_back(chosen_species);
                
                _species_tree_height.clear();
                _log_coalescent_likelihood_options.clear();
                _log_weight_options.clear();
                if (_inf) {
                    double neg_inf = -1*numeric_limits<double>::infinity();
                    _log_weight = neg_inf;
                }
            }
            
            _generation++;
        }
        
        if (_running_on_empty) {
            _generation++;
            _log_weight = 0.0;
        }
        resetVariables();
        return _log_weight;
    }

    inline void Particle::firstGen() {
        assert (_generation == 0); // species tree should only be rebuilt in the first round of the first iteration of gene trees
//        buildFakeSpeciesTree();
        buildEntireSpeciesTree();
//        _forests[0].showForest();
//        panmicticSpeciesPartition();
    }

    inline void Particle::geneTreeProposal(bool unconstrained) {
        if (_generation == 0 && !unconstrained) {
            for (int i=1; i<_forests.size(); i++) {
                _forests[i].deconstructGeneTree();
            }
        }
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._theta = _forests[i]._starting_theta;
            assert (_t.size() == 1);
            pair<double, string> species_info = _forests[i].chooseDelta(_t[0], false);
//            pair<double, string> species_info = _forests[i].chooseDelta(_t[0], unconstrained);
            _forests[i].geneTreeProposal(species_info, _t[0]);
//            _forests[i].showForest();
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

    inline double Particle::getRunningSum() {
        double running_sum = 0.0;
        double log_particle_sum = 0.0;

        double log_max_weight = *max_element(_log_weight_options.begin(), _log_weight_options.end());
        
        for (auto & i:_log_weight_options) {
            running_sum += exp(i - log_max_weight);
        }
        log_particle_sum = log(running_sum) + log_max_weight;

        return log_particle_sum;
    }

    inline void Particle::normalizeSpeciesTreeWeights() {
        double log_particle_sum = getRunningSum();
        
        for (int a=0; a<_log_weight_options.size(); a++) {
            _log_weight_options[a] = _log_weight_options[a] - log_particle_sum;
        }
    }

    inline void Particle::resampleParticles() {
        vector<pair<double, double>> cum_probs;
        vector<double> alt_log_weights;
        vector<vector<Forest>> to_forests;
        vector<double> heights;
        heights.resize(_nspecies_forests);
        
        vector<double> log_coal_likelihoods;
        log_coal_likelihoods.resize(_nspecies_forests);
        
        vector<vector<pair<tuple<string, string, string>, double>>> t;
        t.resize(_nspecies_forests); // TODO: check this one
        
        alt_log_weights.resize(_nspecies_forests);
        // Create vector of pairs p, with p.first = log weight and p.second = particle index
        cum_probs.resize(_nspecies_forests);
        unsigned i = 0;
//        for (int i=0; i<_alt_forests.size(); i++) {
//            _alt_forests[i][0].showForest();
//        }
        for(int a=0; a<_log_weight_options.size(); a++) {
            cum_probs[i].first = _log_weight_options[a];
            cum_probs[i].second = i;
            ++i;
        }
            
        // Sort cum_probs vector so that largest log weights come first
        sort(cum_probs.begin(), cum_probs.end(), greater< pair<double,unsigned> >());

        // Convert vector from storing log weights to storing cumulative weights
        double cumpr = 0.0;
        for (auto & w : cum_probs) {
            cumpr += exp(w.first);
            w.first = cumpr;
        }
        
        // Last element in cum_probs should hold 1.0 if weights were indeed normalized coming in
        if (!isnan(cumpr)) {
        assert( fabs( 1.0 - cum_probs[_nspecies_forests-1].first ) < 0.0000001);
        }

        
        // Draw new set of particles by sampling with replacement according to cum_probs
        to_forests.resize(_nspecies_forests);
        for (unsigned i = 0; i < _nspecies_forests; i++) {
        
            // Select a particle to copy to the ith slot in to_particles
            int sel_index = -1;
            double u = rng.uniform();
            for(unsigned j = 0; j < _nspecies_forests; j++) {
                if (u < cum_probs[j].first) {
                    sel_index = cum_probs[j].second;
                    break;
                }
            }
            assert(sel_index > -1);
            vector<Forest> forests = _alt_forests[sel_index];
            to_forests[i] = forests;
            heights[i] = _species_tree_height[sel_index];
            alt_log_weights[i] = _log_weight_options[sel_index];
            log_coal_likelihoods[i] = _log_coalescent_likelihood_options[sel_index];
            assert(_nspecies_forests == to_forests.size());
            t[i] = _t[sel_index];
        }
        _log_weight_options.clear();
        _log_weight_options = alt_log_weights;
        _log_coalescent_likelihood_options.clear();
        _log_coalescent_likelihood_options = log_coal_likelihoods;
        _species_tree_height.clear();
        _species_tree_height = heights;
        _t.clear();
        _t = t;
        _alt_forests = to_forests;
    }

    inline void Particle::filterSpeciesTrees() {
//        for (auto &forest:_alt_forests) {
//            forest[0].showForest();
//        }
        normalizeSpeciesTreeWeights();
        
        double ess_inverse = 0.0;
        for (int b = 0; b<_log_weight_options.size(); b++) {
            ess_inverse += exp(_log_weight_options[b]);
        }
        
        double ess = 1.0/ess_inverse;
//        cout << "   " << "ESS = " << ess << endl;
        
        if (isnan(ess)) {
            _inf = true;
            throw XProj(format("no species trees can accommodate the gene trees; increase the number of particles and try again"));
        }
     
        if (!_inf) {
        resampleParticles();
        }
//        for (auto &forest:_alt_forests) {
//            forest[0].showForest();
//        }
        //if use_first is true, my_vec = my_vec_2
        //if use_first is false, my_vec = my_vec_1
        
//        _alt_forests = _use_first ? _alt_forests2:_alt_forests;

        //change use_first from true to false or false to true
//        _use_first = !_use_first;
    }

    inline void Particle::speciesProposal() {
//        assert (!_inf);
        if (!_inf) {
        for (int a = 0; a<_alt_forests.size(); a++) {
//            _alt_forests[a][0].showForest();
            tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
            double prev_log_coalescent_likelihood = _log_coalescent_likelihood_options[a];
            if (_alt_forests[a][0]._last_edge_length > 0.0) {
            // choose species to join if past the first species generation for each forest vector
                species_joined = _alt_forests[a][0].speciesTreeProposal();
            }
            
            vector<double> max_depth_vector;
            double max_depth = 0.0;

            for (int i=1; i<_alt_forests[a].size(); i++) {
                _alt_forests[a][i].updateSpeciesPartitionTwo(species_joined);
//                    _alt_forests[a][i].resetDepthVector(species_joined);
                string species1 = get<0>(species_joined);
                string species2 = get<1>(species_joined);
                
                bool match1 = false;
                bool match2 = false;
                
                if (_alt_forests[a][i].getMinDepths()[0].second.first == species1 || _alt_forests[a][i].getMinDepths()[0].second.first == species2) {
                    match1 = true;
                }
                if (_alt_forests[a][i].getMinDepths()[0].second.second == species2 || _alt_forests[a][i].getMinDepths()[0].second.second == species1) {
                    match2 = true;
                }
                    
                if (_alt_forests[a][0]._lineages.size() > 1) {
                    max_depth = (_alt_forests[a][i].getMinDepths())[0].first;
                    max_depth_vector.push_back(max_depth);
                    _alt_forests[a][i].resetDepthVector(species_joined);
                }
                
//                if (match1 && match2) { // TODO: I don't think is necessary - resetdepthvector will take care of it
////                    if (_alt_forests[a][i].getMinDepths()[0].second == species_joined) {
//                    _reset_min.push_back(true);
//                    _alt_forests[a][i].trimDepthVector();
//                }
                else {
                    _reset_min.push_back(false);
                }
                
//                if (_alt_forests[a][0]._lineages.size() > 1) {
//                    max_depth = (_alt_forests[a][i].getMinDepths())[0].first;
//                    max_depth_vector.push_back(max_depth);
//                    _alt_forests[a][i].resetDepthVector(species_joined);
//                }
            }
            if (_alt_forests[a][0]._lineages.size() > 1) {
                max_depth = *min_element(max_depth_vector.begin(), max_depth_vector.end());
//                cout << _alt_forests[a][0].getTreeHeight() << endl;
                max_depth -= _alt_forests[a][0].getTreeHeight();
                // choose a species tree increment
            }
            
            if (_alt_forests[a][0]._lineages.size() > 1) {
                _alt_forests[a][0].chooseSpeciesIncrement(0.0);
//                assert (max_depth > 0.0);
//                _alt_forests[a][0].chooseSpeciesIncrement(max_depth);
                _species_tree_height[a] += _alt_forests[a][0]._last_edge_length;
            }
            if (_alt_forests[a][0]._lineages.size() == 1) {
                _alt_forests[a][0]._last_edge_length = 0.0;
            }
                    
            _t[a].push_back(make_pair(species_joined, _alt_forests[a][0]._last_edge_length));
            
            for (int i = 1; i<_alt_forests[a].size(); i++) {
                double coal_like_increment = _alt_forests[a][i].calcCoalescentLikelihood(_alt_forests[a][0]._last_edge_length, species_joined, _species_tree_height[a], true);
                _log_coalescent_likelihood_options[a] += coal_like_increment;
            }
            if (!_running_on_empty) {
                double nlineages = _alt_forests[a][0]._lineages.size();
                double phi = 0.0;
                if (nlineages > 1) {
                    phi = log(1-exp(-1*max_depth*nlineages*Forest::_speciation_rate));
                }
                _log_weight_options[a] = _log_coalescent_likelihood_options[a] - prev_log_coalescent_likelihood + phi;
            }
        }
//            _alt_forests[a][0].showForest();
        }
        
//         recalculate log weight now
//        _log_weight = 0.0;
//        for (int b=0; b<_log_weight_options.size(); b++) {
//            bool is_infinity = isinf(_log_weight_options[b]);
//            if (!is_infinity) {
//                _log_weight += _log_weight_options[b]; // TODO: need to include prior in weight
//            }
//        }
//        if (_log_weight == 0.0) {
//            double neg_inf = -1*numeric_limits<double>::infinity();
//            _log_weight = neg_inf;
//        }
//        assert (_log_weight !=0.0);
//    cout << "stop";
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

    inline void Particle::buildFakeSpeciesTree() {
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        _forests[0]._last_edge_length = 0.0;

        double edge_len = _forests[0]._last_edge_length;
        _t.resize(1);

//        _t.push_back(make_pair(species_joined, edge_len));

        for (int i=0; i < _forests[0]._nspecies-1; i++) {
            if (_forests[0]._lineages.size() > 1) {
//                tuple<string, string, string> species_joined = _forests[0].speciesTreeProposal();
                tuple<string, string, string> species_joined = make_tuple("null", "null", "null");

                double edge_len = 0.0;
                if (_forests[0]._lineages.size() > 1) {
                    edge_len = _forests[0]._last_edge_length;
                }
                _t[0].push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline void Particle::panmicticSpeciesPartition() {
        for (int i=1; i<_forests.size(); i++) {
            _forests[i].combineSpeciesPartition();
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
            string newick = "newick:\"" + _forests[i].makeNewick(9, true) + "\"";
            newicks.push_back(newick);
        }
        return newicks;
    }

    inline void Particle::resetSpecies() {
        setLogLikelihood(0.0);
        setLogWeight(0.0);
        resetSpeciesInfo();
        resetSpeciesTreeHeight();
        _forests[0]._increments.clear();
    }

    inline void Particle::resetGeneIncrements() {
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._increments.clear();
        }
    }

    inline void Particle::increaseNumberOfForests() {
        _t.resize(_nspecies_forests);
        _species_tree_height.resize(_nspecies_forests);
        _log_coalescent_likelihood_options.resize(_nspecies_forests);
        _log_weight_options.resize(_nspecies_forests);
        
        for (int i=0; i<_nspecies_forests; i++) {
            vector<Forest> f = _forests;
            _alt_forests.push_back(f);
        }
    }

    inline void Particle::buildEntireSpeciesTree() {
        double max_depth = 0.0;
        _forests[0].chooseSpeciesIncrement(max_depth);
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        double edge_len = _forests[0]._last_edge_length;
        vector<pair<tuple<string, string, string>, double>> test;
        test.push_back(make_pair(species_joined, edge_len));

        for (int i=0; i < _forests[0]._nspecies-1; i++) {
            if (_forests[0]._lineages.size() > 1) {
//                species_joined = _forests[0].speciesTreeProposal();
                pair<unsigned, unsigned> example;
                if (i == 0) {
//                    species_joined = make_tuple("s1", "s4", "node-5");
                    example = make_pair(1, 4);
                }
                else if (i == 1) {
//                    species_joined = make_tuple("s2", "s3", "node-6");
                    example = make_pair(1,2 );
                }
                else if (i == 2) {
//                    species_joined = make_tuple("node-5", "node-6", "node-7");
                    example = make_pair(1, 2);
                }
                else if (i == 3) {
//                    species_joined = make_tuple("s0", "node-7", "node-8");
                    example = make_pair(0, 1);
                }
                species_joined = _forests[0].preDeterminedSpeciesTreeProposal(example);

                // if the species tree is not finished, add another species increment
                if (_forests[0]._lineages.size()>1) {
                    _forests[0].addSpeciesIncrement();
                }

                double edge_len = 0.0;
                if (_forests[0]._lineages.size() > 1) {
                    edge_len = _forests[0]._last_edge_length;
                }
                test.push_back(make_pair(species_joined, edge_len));
//                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
//        _forests[0].showForest();
        _t.push_back(test);
    }

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _alt_forests = other._alt_forests;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _t = other._t;
        _gene_tree_proposal_attempts = other._gene_tree_proposal_attempts;
        _ready_to_join_species = other._ready_to_join_species;
        _running_on_empty = other._running_on_empty;
        _random_seeds = other._random_seeds;
        _species_first = other._species_first;
        _no_species_joined = other._no_species_joined;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _species_tree_height = other._species_tree_height;
        _nspecies_forests = other._nspecies_forests;
        _log_coalescent_likelihood_options = other._log_coalescent_likelihood_options;
        _log_weight_options = other._log_weight_options;
        _inf = other._inf;
        _reset_min = other._reset_min;
    };
}
