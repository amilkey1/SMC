#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
#include <mutex>
#include <unordered_map>

using namespace std;
using namespace boost;

#include "lot.hpp"
#include "conditionals.hpp"
#include "stopwatch.hpp"

extern proj::Lot rng;
std::mutex mutx;

namespace proj {

class Particle {
    public:

        Particle();
        Particle(const Particle & other);
        typedef std::shared_ptr<Particle>               SharedPtr;

        void                                    debugParticle(std::string name);
        void                                    showParticle();
        void                                    showSpeciesParticle();
        void                                    proposal();
        void                                    setData(Data::SharedPtr d, map<string, string> &taxon_map, bool partials) {
                                                    int index = 0;
                                                    _forests.resize(G::_nloci+1);
                                                    for (auto &_forest:_forests) {
                                                        if (index > 0) {
                                                            _forest.setData(d, index, taxon_map, partials);
                                                        }
                                                        index++;
                                                    }
                                                }
        void                                    setSimData(Data::SharedPtr d, map<string, string> &taxon_map, unsigned ntaxa) {
                                                    int index = 0;
                                                    _forests.resize(G::_nloci+1);
                                                    for (auto &_forest:_forests) {
                                                        if (index > 0) {
                                                            _forest.setSimData(d, index, taxon_map, ntaxa);
                                                        }
                                                        index++;
                                                    }
                                                }
        void                                    mapSpecies(map<string, string> &taxon_map);
        void                                    saveForest(std::string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        double                                  getLogLikelihood();
        vector<double>                          calcGeneTreeLogLikelihoods();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
#if defined (DEBUG_MODE)
        double                                  getSpeciesIncrement () {return _forests[0]._last_edge_length;}
#endif
        double                                  getSpeciesLogWeight() const {return _log_weight;}
        unsigned                                getPartialCount();
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    setLogSpeciesWeight(double w){_log_weight = w;}
#if defined (UNUSED_FUNCTIONS)
        void                                    setLogLikelihood(vector<double> forest_likelihoods);
#endif
        void                                    setLogCoalescentLikelihood(double coalescent_like);
        void                                    operator=(const Particle & other);
        const vector<Forest> &                  getForest() const {return _forests;}
        vector<double>                          getThetaMap();
        double                                  getThetaMean(){return _theta_mean;}
        string                                  saveForestNewick() {
            return _forests[0].makeNewick(8, true);}
        string                                  saveForestNewickAlt() {return _forests[0].makeAltNewick(8, false);}
            
        string                                  saveGeneNewick(unsigned i) {
            return _forests[i].makeNewick(8, true);}
        string                                  saveChangedForest() {return _forests[_gene_order[G::_generation-1]].makePartialNewick(8, true);}
#if defined (USING_MPI)
        unsigned                                getNextGene(){return _gene_order[G::_generation];}
        void                                    initSpeciesForest(string newick);
        void                                    initGeneForest(string newick);
        void                                    checkPartition();
        void                                    initGeneForestSpeciesPartition(string species_partition);
        string                                  saveForestSpeciesPartition();
#endif
        void                                    setGeneUPGMAMatrices();
        void                                    createThetaMap();
        void                                    createThetaMapFixedTheta();
        void                                    updateThetaMap(string new_species_name);
        void                                    updateThetaMapFixedTheta(string new_species_name);
#if !defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
        void                                    resetThetaMap(Lot::SharedPtr lot, unordered_map<string, double> &theta_map);
#endif
        void                                    drawNewTheta(string new_species);
    
        bool operator<(const Particle::SharedPtr & other) const {
            return _log_weight<other->_log_weight;
        }

        bool operator>(const Particle::SharedPtr & other) const {
            return _log_weight>other->_log_weight;
        }

        void                                            setGroupNumber(unsigned n) {_group_number = n;} // group number for parallelization
        unsigned                                        getGroupNumber() {return _group_number;}// group number for parallelization
        vector<Forest> &                                getForests() {return _forests;}
#if defined (DEBUG_MODE)
        void                                            showSpeciesJoined();
        void                                            showSpeciesIncrement();
        pair<string, string>                            getSpeciesJoined(){return make_pair(_forests[0]._species_joined.first->_name, _forests[0]._species_joined.second->_name);}
#endif
        void                                            showSpeciesTree();
        int                                             selectEventLinearScale(vector<double> weight_vec);
        double                                          getTopologyPrior(unsigned i);
        vector<pair<double, double>>                    getIncrementPriors(unsigned i);
        vector<pair<double, double>>                    getSpeciesTreeIncrementPriors();
        double                                          getCoalescentLikelihood(unsigned g);
        void                                            clear();
        void                                            speciesOnlyProposal();
#if !defined (FASTER_SECOND_LEVEL)
        void                                            speciesOnlyProposalIntegratingOutTheta();
#endif
        void                                            drawTheta();
        void                                            fixTheta();
        void                                            clearPartials();
        Lot::SharedPtr getLot() const {return _lot;}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}
        double                                          getSpeciesTreeHeight();
        double                                          getSpeciesTreeLength();
        vector<double>                                  getGeneTreeHeights();
        vector<double>                                  getGeneTreeLengths();
        vector<double>                                  getGeneTreeLogLikelihoods();
        vector<double>                                  getGeneTreePriors();
        inline vector<double>                           getGeneTreeCoalescentLikelihoods();
        double                                          getSpeciesTreePrior();
        double                                          getAllPriors();
        double                                          getAllPriorsFirstRound();
#if !defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
        vector<double>                                  getVectorPrior();
#endif
        void                                            simulateData(vector<unsigned> sites_vector);
        unsigned                                        getNumDeepCoalescences() {return _num_deep_coalescences;}
        unsigned                                        getMaxDeepCoalescences(){return _max_deep_coal;}
        void                                            resetSpecies();
        void                                            setForest(Forest f, unsigned forest_number);
        Forest                                          getForest(unsigned i) {return _forests[i];} // TODO: should return a pointer?
        void                                            setNewTheta(bool fix_theta);
        vector<double>                                  getThetaVector();
        double                                          getPopMean();
        void                                            setPsuffix(unsigned psuffix) {_psuffix = psuffix;}
        double                                          calcInitialCoalescentLikelihood();
        void                                            processGeneNewicks(vector<string> newicks);
        void                                            processSpeciesNewick(string newick_string);
        void                                            setNextSpeciesNumber() {_next_species_number = G::_nspecies;}
        string                                          getTranslateBlock();
        void                                            buildEntireSpeciesTree();
        void                                            rebuildSpeciesTree();
        void                                            setGeneOrder(vector<unsigned> gene_order) {_gene_order = gene_order;}
        void                                            resetGeneOrder(unsigned step, vector<unsigned> gene_order);
        void                                            trimSpeciesTree();
#if defined (OLD_UPGMA)
        void                                            calcStartingRowCount();
        void                                            calcStartingUPGMAMatrix();
        vector<vector<double>>                          getStartingUPGMAMatrix();
        void                                            setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene);
        void                                            createSpeciesIndices();
#endif
        void                                            setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
#if defined (FASTER_SECOND_LEVEL)
        void                                            saveCoalInfoInitial();
        unsigned                                        proposeSpeciationEvent();
        double                                          findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect);
        double                                          calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose);
        void                                            resetPrevLogCoalLike();
        void                                            clearGeneForests();
#endif
        void                                            setNodeHeights();
    
    private:

#if defined (FASTER_SECOND_LEVEL)
        double                                  _log_coal_like;
        double                                  _prev_log_coal_like;
#endif
        vector<Forest>                          _forests;
        double                                  _log_weight;
        double                                  _log_likelihood;
        double                                  _log_coalescent_likelihood;
        mutable                                 Lot::SharedPtr _lot;
        unsigned                                _num_deep_coalescences;
        unsigned                                _max_deep_coal;
        unsigned                                _psuffix;
        unsigned                                _next_species_number;
        vector<unsigned>                        _next_species_number_by_gene;
        vector<pair<tuple<string, string, string>, double>> _t;
        vector<vector<pair<tuple<string, string, string>, double>>> _t_by_gene; // TODO: should these be maps instead?
        vector<unsigned>                        _gene_order;
    
        unordered_map<string, double>           _theta_map;
        double                                  _theta_mean;

#if defined (UNUSED_FUNCTIONS)
        vector<double>                          _starting_log_likelihoods;
#endif
#if !defined (FASTER_SECOND_LEVEL)
        unsigned                                _species_branches;
        double                                  _species_tree_height;
#endif
        unsigned                                _group_number;
};

    inline Particle::Particle() {
        _lot.reset(new Lot());
        clear();
    };

    inline void Particle::showParticle() {
        //print out weight of each particle
        cout << "\nParticle:\n";
        cout << "  _log_weight: " << _log_weight << "\n" ;
        cout << " _log_likelihood: " << _log_likelihood << "\n";
        cout << " _log_coalescent_likelihood: " << _log_coalescent_likelihood << "\n";
        cout << "  _forest: " << "\n";
#if defined (DRAW_NEW_THETA)
        unsigned a = 1;
        for (auto &t:_theta_map) {
            cout << "theta " << a << " = " << t.second << endl;
            a++;
        }
#endif
        cout << "\n";
        for (auto &_forest:_forests) {
            _forest.showForest();
        }
    }

    inline void Particle::showSpeciesParticle() {
        //print out weight of each species forest part of the particle
        cout << "\nParticle:\n";
        cout << " _log_species_weight: " << _log_weight << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _forests[0].showForest();
    }

    inline void Particle::clear() {
        _log_weight     = 0.0;
        _log_likelihood = 0.0;
        _forests.clear();
        _log_coalescent_likelihood = 0.0;
        _num_deep_coalescences = 0.0;
        _max_deep_coal = 0.0;
        _t.clear();
        _psuffix = 0;
        _next_species_number = G::_nspecies;
        _t_by_gene.clear();
        _next_species_number_by_gene.clear();
        _gene_order.clear();
        _group_number = 0;
#if !defined (FASTER_SECOND_LEVEL)
        _species_branches = 0;
        _species_tree_height = 0.0;
#endif
#if defined (UNUSED_FUNCTIONS)
        _starting_log_likelihoods.clear();
#endif
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
        if (G::_generation == 0 && !G::_run_on_empty) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline vector<double> Particle::calcGeneTreeLogLikelihoods() {
        vector<double> gene_forest_likelihoods;
        gene_forest_likelihoods.resize(_forests.size()-1);
        
        //calculate likelihood for each gene tree
        for (unsigned i=1; i<_forests.size(); i++) {
            double gene_tree_log_likelihood  = 0.0;
            if (!G::_run_on_empty) {
                gene_tree_log_likelihood = _forests[i].calcLogLikelihood();
                assert(!isnan (gene_tree_log_likelihood));
            }
            gene_forest_likelihoods[i-1] = gene_tree_log_likelihood;
        }

        return gene_forest_likelihoods;
    }

    inline void Particle::clearPartials() {
        for (int i=1; i<_forests.size(); i++) {
            _forests[i].clearPartials();
        }
    }

    inline void Particle::setLogCoalescentLikelihood(double log_coal_like) {
        _log_coalescent_likelihood = log_coal_like;
        _log_weight = _log_coalescent_likelihood;
    }

#if defined (UNUSED_FUNCTIONS)
    inline void Particle::setLogLikelihood(vector<double> forest_log_likelihoods) {
        double total_log_likelihood = 0.0;
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._gene_tree_log_likelihood = forest_log_likelihoods[i-1];
            total_log_likelihood += forest_log_likelihoods[i-1];
            _forests[i]._log_weight = forest_log_likelihoods[i-1];
            _starting_log_likelihoods.push_back(forest_log_likelihoods[i-1]);
        }
        _log_likelihood = total_log_likelihood;
        _log_weight = total_log_likelihood;
    }
#endif

    inline double Particle::getSpeciesTreeHeight() {
        return _forests[0]._forest_height;
    }

    inline double Particle::getSpeciesTreeLength() {
        return _forests[0].getTreeLength();
    }

    inline vector<double> Particle::getGeneTreeHeights() {
        vector<double> gene_tree_heights;
        for (int i=1; i<_forests.size(); i++) {
            gene_tree_heights.push_back(_forests[i]._forest_height);
        }
        return gene_tree_heights;
    }

    inline vector<double> Particle::getGeneTreeLengths() {
        vector<double> gene_tree_heights;
        for (int i=1; i<_forests.size(); i++) {
            gene_tree_heights.push_back(_forests[i].getTreeLength());
        }
        return gene_tree_heights;
    }

    inline vector<double> Particle::getGeneTreeLogLikelihoods() {
        vector<double> gene_tree_log_likelihoods;
        for (int i=1; i<_forests.size(); i++) {
            gene_tree_log_likelihoods.push_back(_forests[i]._gene_tree_log_likelihood);
            assert (_forests[i]._gene_tree_log_likelihood <= 0.0);
        }
        return gene_tree_log_likelihoods;
    }

    inline vector<double> Particle::getGeneTreePriors() {
        vector<double> gene_tree_priors;
        for (int i=1; i<_forests.size(); i++) {
            double prior = 0.0;
            for (auto &p:_forests[i]._increments_and_priors) {
                prior += p.second;
            }
            
            gene_tree_priors.push_back(prior);
        }
        return gene_tree_priors;
    }

    inline vector<double> Particle::getGeneTreeCoalescentLikelihoods() {
#if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
        // do not have separate coalescent likelihood for each gene tree
        vector<double> gene_tree_priors;
        for (int i=1; i<_forests.size(); i++) {
            gene_tree_priors.push_back(0);
        }
        return gene_tree_priors;
#else
        vector<double> gene_tree_priors;
        for (int i=1; i<_forests.size(); i++) {
            gene_tree_priors.push_back(_forests[i]._log_coalescent_likelihood);
        }
        return gene_tree_priors;
#endif
    }

    inline double Particle::getSpeciesTreePrior() {
        double prior = 0.0;
        for (auto &p:_forests[0]._increments_and_priors) {
            prior += p.second;
        }
        prior += _forests[0]._log_joining_prob;
        
        return prior;
    }

    inline double Particle::getAllPriorsFirstRound() {
        double species_tree_prior = getSpeciesTreePrior();
        double total_prior = 0.0;
        // starbeast3 does not include gene tree priors separately - this is already accounted for in coalescent likelihood
        total_prior += species_tree_prior;
    #if defined (DRAW_NEW_THETA)
        if (G::_theta_prior_mean > 0.0) {
            total_prior += log(G::_theta_prior_mean) - (G::_theta * G::_theta_prior_mean);
        }
        // TODO: what if theta mean is set?
    #endif
        return total_prior;
    }

    inline double Particle::getAllPriors() {
        double species_tree_prior = getSpeciesTreePrior();
        double total_prior = 0.0;
        // starbeast3 does not include gene tree priors separately - this is already accounted for in coalescent likelihood
        total_prior += species_tree_prior;
        // no prior on theta for second round
        return total_prior;
    }

#if !defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
inline vector<double> Particle::getVectorPrior() {
// this is the InverseGamma(2, psi) prior on the population sizes
    return _forests[1]._vector_prior;
}
#endif

    inline double Particle::getLogLikelihood() {
        //retrieve likelihood for each gene tree
        double log_likelihood = 0.0;
        for (unsigned i=1; i<_forests.size(); i++) {
            double gene_tree_log_likelihood = _forests[i]._gene_tree_log_likelihood;
            assert(!isnan (log_likelihood));
            //total log likelihood is sum of gene tree log likelihoods
            log_likelihood += gene_tree_log_likelihood;
        }
        if (G::_generation == 0 && !G::_run_on_empty) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline void Particle::proposal() {
        double inv_gamma_modifier = 0.0;
        
        unsigned next_gene = _gene_order[G::_generation];
        bool calc_weight = false;
        
        if (G::_species_newick_name == "null") { // if specifying species newick, keep the species tree the same for all particles and never reset it in first level
            if (G::_generation == 0) {
                buildEntireSpeciesTree();
                // make a separate species tree information vector for each gene
                for (unsigned i=1; i<_forests.size(); i++) {
                    _t_by_gene.push_back(_t);
                    _next_species_number_by_gene.push_back(0);
                }
            }
#if !defined (USING_MPI)
            else if (G::_generation % G::_nloci == 0 && G::_start_mode_type == G::StartModeType::START_MODE_SMC) { // after every locus has been filtered once, trim back the species tree as far as possible & rebuild it
                    trimSpeciesTree();
                if (_forests[0]._lineages.size() > 1) {
                    rebuildSpeciesTree();
                }
            }
#endif
        }
        else { // TODO: for now, don't rebuild tree for MPI
            if (G::_generation == 0) {
                for (unsigned i=1; i<_forests.size(); i++) {
                    _t_by_gene.push_back(_t);
                    _next_species_number_by_gene.push_back(0);
                }
            }
        }
        
        bool done = false;
                
        while (!done) {
            vector<pair<double, string>> rates_by_species = _forests[next_gene].calcForestRate(_lot, _theta_map);
            double total_rate = 0.0;
            double gene_increment = -1.0;
            if (rates_by_species.size() > 0) {
                for (auto &r:rates_by_species) {
                    total_rate += r.first;
                }
                assert (total_rate > 0.0);
                gene_increment = -log(1.0 - _lot->uniform())/total_rate;
                assert (gene_increment > 0.0);
            }
            
            unsigned next_species_index = _next_species_number_by_gene[next_gene-1];
            double species_increment = _t_by_gene[next_gene-1][next_species_index].second;
            
            string species_name = "species";
            
          // if total rate is 0, gene increment will be -1.0, which will be taken care of

            if ((gene_increment < species_increment || species_increment == 0.0) && gene_increment != -1.0) { // if species increment is 0.0, choose a coalescent event because the species tree is finished
                
                assert (gene_increment > 0.0);
                _forests[next_gene].addIncrement(gene_increment);
                
                vector<double> event_choice_rates;
                for (auto &r:rates_by_species) {
                    event_choice_rates.push_back(r.first / total_rate);
                }
                
                unsigned index = selectEventLinearScale(event_choice_rates);
                string species_name = rates_by_species[index].second;
                _forests[next_gene].allowCoalescence(species_name, gene_increment, _lot);
                
#if defined (OLD_UPGMA)
                if (G::_start_mode_type == G::StartModeType::START_MODE_SMC) {
                    if (G::_upgma) {
                        if (!G::_run_on_empty) {
                            _forests[next_gene].buildRestOfTreeFaster();
                        }
                    }
                }
#endif
                    
                if (species_increment > 0.0) { // otherwise, species tree is done and there is nothing left to update
                    _t_by_gene[next_gene-1][next_species_index].second -= gene_increment; // update species tree increments
                }
                    calc_weight = true;
                }
                else {
                    // carry out speciation event
                    
                    assert (species_increment > 0.0);
                    assert (_forests[next_gene]._species_partition.size() > 1);
                    _forests[next_gene].addIncrement(species_increment);
                    
                    if (G::_start_mode_type == G::StartModeType::START_MODE_SIM) {
                        _num_deep_coalescences += _forests[next_gene].getDeepCoal(_t_by_gene[next_gene - 1][next_species_index + 1].first);
                        _max_deep_coal += _forests[next_gene].getMaxDeepCoal(_t_by_gene[next_gene - 1][next_species_index + 1].first);
                    }
                    
                    _forests[next_gene].updateSpeciesPartition(_t_by_gene[next_gene-1][next_species_index+1].first);
                    assert (next_species_index < _t_by_gene[next_gene-1].size());
                    _t_by_gene[next_gene-1][next_species_index].second -= species_increment; // update species tree increments
                    assert (_t_by_gene[next_gene-1][next_species_index].second == 0.0);
                    if (_forests[next_gene]._species_partition.size() > 1) {
                        _next_species_number_by_gene[next_gene-1]++;
                }
            }
                
        
            if (calc_weight) { // calc weight just means coalescent event has been proposed
                done = true;
            }
        }
         
         done = true;
        
#if defined (INV_GAMMA_PRIOR_TWO)
        // turn this off when using 2.0, not 2.01, for the inverse gamma - don't include a correction in this case
            // include inverse gamma prior correction for every species population for every locus at every step
        double theta_mean = _forests[1]._theta_mean;
        double eps = 0.01;
        double a = 2.0;
        
        inv_gamma_modifier = lgamma(a + eps) - lgamma(a) + a * log(a-1.0) - (a + eps) * log(a + eps - 1.0);
        
        for (auto &t:_forests[next_gene]._theta_map) {

            double y = t.second; // theta
            inv_gamma_modifier += eps * (log(y) - log(theta_mean)) + (theta_mean * eps / y);
        }
#endif
        
#if defined (WEIGHT_MODIFIER)
        // modifier only happens on first round
        if (G::_generation == 0 && G::_theta_prior_mean > 0.0 && G::_theta_proposal_mean > 0.0) {
            if (G::_theta_prior_mean != G::_theta_proposal_mean) {
                // else, log weight modifier is 0
                double prior_rate = 1.0/G::_theta_prior_mean;
                double proposal_rate = 1.0/G::_theta_proposal_mean;
                double log_weight_modifier = log(prior_rate) - log(proposal_rate) - (prior_rate - proposal_rate)*_theta_mean;

                _log_weight += log_weight_modifier;
            }
        }
#endif
        
        if (G::_start_mode_type == G::StartModeType::START_MODE_SMC && !G::_run_on_empty) {
            _log_weight = _forests[next_gene]._log_weight + inv_gamma_modifier;
        }
        else {
            _log_weight = 0.0;
        }
        
    }

#if defined (FASTER_SECOND_LEVEL)
    struct bitless {
        bool operator()(const G::species_t a, const G::species_t b) const {
            bool returned_value = ((a & b) > 0 ? false : a < b);
            return returned_value;
        }
    };
#endif

#if !defined (FASTER_SECOND_LEVEL)
    inline void Particle::speciesOnlyProposalIntegratingOutTheta() {
        if (_generation == 0) {
            _species_branches = G::_nspecies;
            for (int i=1; i<_forests.size(); i++) {
                _forests[i].refreshPreorder();
                _forests[i].calcMinDepth();
                _forests[i]._nincrements = 0;

                if (i > 1) {
                    _forests[i]._theta_mean = _forests[1]._theta_mean;
                }
            }
        }
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        double prev_log_coalescent_likelihood = _log_coalescent_likelihood;
        
#if !defined (COAL_LIKE_TEST)
        
            if (_forests[0]._last_edge_length > 0.0) {
            // choose species to join if past the first species generation for each forest vector
                species_joined = _forests[0].speciesTreeProposal(_lot);
            }
                vector<double> max_depth_vector;
                double max_depth = 0.0;

                for (int i=1; i<_forests.size(); i++) {
                    string species1 = get<0>(species_joined);
                    string species2 = get<1>(species_joined);
                        
                    if (_forests[0]._lineages.size() > 1 && species1 != "null") {
                        // if using Jones formula, species partition update will happen in coalescent likelihood calculation
                        _forests[i].resetDepthVector(species_joined);
                    }

                    max_depth = (_forests[i].getMinDepths())[0].first;
                    max_depth_vector.push_back(max_depth);
                }
                if (_forests[0]._lineages.size() > 1) {
                    max_depth = *min_element(max_depth_vector.begin(), max_depth_vector.end());
                    max_depth -= _forests[0]._forest_height;
                    // choose a species tree increment
                }
                
                if (_forests[0]._lineages.size() > 1) {
#if !defined (UNCONSTRAINED_PROPOSAL)
                    
                    assert (max_depth > 0.0);
                    
                    _forests[0].chooseSpeciesIncrementOnly(_lot, max_depth);
#else
                    _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
#endif
                    _species_tree_height += _forests[0]._last_edge_length;
                }
                if (_forests[0]._lineages.size() == 1) {
                    _forests[0]._last_edge_length = 0.0;
                }
        
                _log_coalescent_likelihood = 0.0;
        
            assert (_log_coalescent_likelihood == 0.0);
#endif

        if (!G::_run_on_empty || G::_run_on_empty_first_level_only) {
            _species_branches += 1;
            _log_coalescent_likelihood = 2 * _species_branches * log(_forests[1]._theta_mean) - _species_branches * boost::math::lgamma(2);
            
            vector<double> gamma_jb;
            vector<unsigned> q_jb;
            
            double neg_inf = -1*numeric_limits<double>::infinity();
            
            for (int i = 1; i<_forests.size(); i++) {
                if (_log_coalescent_likelihood != neg_inf) {
                    pair<vector<double>, vector<unsigned>> params;
#if defined (COAL_LIKE_TEST)
                    showParticle();
                    Forest::_theta = 0.07;
//                    calcInitialCoalescentLikelihood();
                    params = _forests[i].calcInitialCoalescentLikelihoodIntegratingOutTheta();
#else
                    sw.start();
                    
                    if (_forests[0]._lineages.size() > 1) {
                        params = _forests[i].calcCoalescentLikelihoodIntegratingOutTheta(_forests[0]._species_build);
                    }
                    else {
                        params = _forests[i].calcCoalescentLikelihoodIntegratingOutThetaLastStep(_forests[0]._species_build); // not using this right now
                    }
#endif
                    if (i == 1) {
                        for (auto &g:params.first) {
                            if (g == neg_inf) {
                                _log_coalescent_likelihood = neg_inf;
                                break;
                            }
                            gamma_jb.push_back(g);
                        }
                        for (auto &q:params.second) {
                            q_jb.push_back(q);
                        }
                    }
                    else {
                        for (unsigned p=0; p<params.first.size(); p++) {
                            gamma_jb[p] += params.first[p];
                            if (params.first[p] == neg_inf) {
                                _log_coalescent_likelihood = neg_inf;
                                break;
                            }
                            q_jb[p] += params.second[p];
                        }
                    }
                }
            }
                   
            sw.start();
            if (_forests[0]._lineages.size() > 2) {
                for (unsigned p=0; p<gamma_jb.size(); p++) {
                    double log_rb = q_jb[p] * log((4 / _forests[1]._ploidy));
                    double q_b = q_jb[p];
                    double gamma_b = gamma_jb[p];
                    
                    if (gamma_b == neg_inf) {
                        _log_coalescent_likelihood = neg_inf;
                    }
                    
                    if (_log_coalescent_likelihood != neg_inf) {
                        double test = log_rb - (2+q_b)*log(_forests[1]._theta_mean + gamma_b) + boost::math::lgamma(2 + q_b);
                        _log_coalescent_likelihood += log_rb - (2+q_b)*log(_forests[1]._theta_mean + gamma_b) + boost::math::lgamma(2 + q_b);
                    }
                }
            }
            
            else {
                species_joined = _forests[0].speciesTreeProposal(_lot);
                for (unsigned p=0; p<gamma_jb.size(); p++) {
                    double log_rb = q_jb[p] * log((4 / _forests[1]._ploidy));
                    double q_b = q_jb[p];
                    double gamma_b = gamma_jb[p];

                    if (gamma_b == neg_inf) {
                        _log_coalescent_likelihood = neg_inf;
                    }

                    if (_log_coalescent_likelihood != neg_inf) {
                        _log_coalescent_likelihood += log_rb - (2+q_b)*log(_forests[1]._theta_mean + gamma_b) + boost::math::lgamma(2 + q_b);
                    }
                }
            }
        }

#if defined (COAL_LIKE_TEST)
        double max_depth = 0.0;
        double constrained_factor = 0.0;
#endif
        
#if !defined (COAL_LIKE_TEST)
        double constrained_factor = 0.0;
#if !defined (UNCONSTRAINED_PROPOSAL)
        if (!G::_run_on_empty) {
            assert (max_depth > 0.0);
            double nlineages = _forests[0]._lineages.size();
            if (nlineages == 1) {
                nlineages = 2; // for last step, constraint was before final two species were joined
            }
//            constrained_factor = log(1 - exp(-1*nlineages*Forest::_lambda*max_depth));
            constrained_factor = log(1 - exp(-1*nlineages*G::_lambda*max_depth));
        }
#endif
        
#endif
            _log_weight = _log_coalescent_likelihood - prev_log_coalescent_likelihood + constrained_factor;
#if !defined (UNCONSTRAINED_PROPOSAL)
            double test = 1/_log_weight;
            assert(test != -0); // assert coalescent likelihood is not -inf
#endif
        
        if (G::_run_on_empty && !G::_run_on_empty_first_level_only) {
            assert (_log_coalescent_likelihood == 0.0);
            assert (_log_weight == 0.0);
        }
        _generation++;
    }
#endif

    inline void Particle::speciesOnlyProposal() {
#if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
#if !defined (FASTER_SECOND_LEVEL)
        speciesOnlyProposalIntegratingOutTheta();
#endif
#else
        
        if (_generation == 0) {
#if !defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
            _forests[1]._vector_prior.clear();
#endif
            for (int i=1; i<_forests.size(); i++) {
                _forests[i].refreshPreorder();
                _forests[i].calcMinDepth();
                _forests[i]._nincrements = 0;
#if defined (DRAW_NEW_THETA)
                _forests[i]._theta_map.clear(); // clear old thetas
#endif
                if (i > 1) {
                    _forests[i]._theta_mean = _forests[1]._theta_mean;
                    _forests[i]._species_names = _forests[1]._species_names;
                }
            }
#if defined (DRAW_NEW_THETA)
            for (int i=1; i<_forests.size(); i++) {
                assert (_forests[i]._theta_mean > 0.0);
            }
            _forests[1].resetThetaMap(_lot); // reset tip thetas and ancestral pop theta
            if (_forests.size() > 2) {
                for (int i=2; i<_forests.size(); i++) {
                    _forests[i]._theta_map = _forests[1]._theta_map;
                }
            }
#endif
        }
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        double prev_log_coalescent_likelihood = _log_coalescent_likelihood;
        
            if (_forests[0]._last_edge_length > 0.0) {
            // choose species to join if past the first species generation for each forest vector
                species_joined = _forests[0].speciesTreeProposal(_lot);
            }
                vector<double> max_depth_vector;
                double max_depth = 0.0;

                for (int i=1; i<_forests.size(); i++) {
                    string species1 = get<0>(species_joined);
                    string species2 = get<1>(species_joined);
                    
                    if (species1 != "null") {
                        _forests[i].updateSpeciesPartition(species_joined); // if using Jones formula, this will happen in coalescent likelihood calculation
                    }
                    
                    if (_forests[0]._lineages.size() > 1 && species1 != "null") {
                        _forests[i].resetDepthVector(species_joined);
                    }

                    max_depth = (_forests[i].getMinDepths())[0].first;
                    max_depth_vector.push_back(max_depth);
                }
                if (_forests[0]._lineages.size() > 1) {
                    max_depth = *min_element(max_depth_vector.begin(), max_depth_vector.end());
                    max_depth -= _forests[0]._forest_height;
                    // choose a species tree increment
                }
                
                if (_forests[0]._lineages.size() > 1) {
#if !defined (UNCONSTRAINED_PROPOSAL)
                    assert (max_depth > 0.0);
                    _forests[0].chooseSpeciesIncrementOnly(_lot, max_depth);
#else
                    _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
#endif
                    _species_tree_height += _forests[0]._last_edge_length;
                }
                if (_forests[0]._lineages.size() == 1) {
                    _forests[0]._last_edge_length = 0.0;
                }
        
                _t.push_back(make_pair(species_joined, _forests[0]._last_edge_length));
                
                _log_coalescent_likelihood = 0.0;
        
#if defined (DRAW_NEW_THETA)
        if (_generation > 0) {
            _forests[1].drawNewTheta(get<2>(species_joined), _lot); // each time species are joined, draw a new theta for the new population and ancestral pop
            if (_forests.size() > 2) {
                for (int i=2; i<_forests.size(); i++) {
                    _forests[i]._theta_map = _forests[1]._theta_map;
                }
            }
        }
#endif
            assert (_log_coalescent_likelihood == 0.0);
        
        if (_forests[0]._last_edge_length > max_depth) {
            _log_coalescent_likelihood = -1*numeric_limits<double>::infinity();
        }
        else {
            for (int i = 1; i<_forests.size(); i++) {
                _forests[i].calcCoalescentLikelihood(_forests[0]._last_edge_length, species_joined, _species_tree_height);
                _log_coalescent_likelihood += _forests[i]._log_coalescent_likelihood + _forests[i]._panmictic_coalescent_likelihood;
            }
        }
        
        if (_forests[0]._lineages.size() == 2) {
            // join remaining species lineages, no change in coalescent likelihood, just need to add panmictic coalescent for each gene tree (to avoid total recalculation)
            // no need to draw a new theta because we are at the ancestral population now
            species_joined = _forests[0].speciesTreeProposal(_lot);
            for (int i=1; i<_forests.size(); i++) {
                _forests[i]._log_coalescent_likelihood += _forests[i]._panmictic_coalescent_likelihood;
                _forests[i]._panmictic_coalescent_likelihood = 0.0; // for clarity, reset to 0
            }
        }
        
        double constrained_factor = 0.0;
#if !defined (UNCONSTRAINED_PROPOSAL)
        assert (max_depth > 0.0);
        double nlineages = _forests[0]._lineages.size();
        constrained_factor = log(1 - exp(-1*nlineages*Forest::_lambda*max_depth));
#endif
        _log_weight = _log_coalescent_likelihood - prev_log_coalescent_likelihood + constrained_factor;
#if !defined (UNCONSTRAINED_PROPOSAL)
        double test = 1/_log_weight;
        assert(test != -0); // assert coalescent likelihood is not -inf
#endif
        double neg_inf = -1*numeric_limits<double>::infinity();
        if (_log_coalescent_likelihood == neg_inf) {
            assert (_forests[0]._last_edge_length > max_depth);
        }
        
        _generation++;
#endif
    }

#if defined (OLD_UPGMA)
    inline void Particle::calcStartingUPGMAMatrix() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].buildStartingUPGMAMatrix();
        }
    }
#endif

#if defined (OLD_UPGMA)
    inline vector<vector<double>> Particle::getStartingUPGMAMatrix() {
        vector<vector<double>> starting_upgma_matrices_by_gene;
        for (unsigned i=1; i<_forests.size(); i++) {
            starting_upgma_matrices_by_gene.push_back(_forests[i]._starting_dij);
        }
        return starting_upgma_matrices_by_gene;
    }
#endif

#if defined (OLD_UPGMA)
    inline void Particle::setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene) {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._starting_dij = starting_upgma_matrices_by_gene[i-1];
        }
    }
#endif

#if defined (OLD_UPGMA)
    inline void Particle::calcStartingRowCount() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].buildStartingRow();
        }
    }
#endif

#if !defined (FASTER_SECOND_LEVEL)
    inline double Particle::calcInitialCoalescentLikelihood() {
#if defined GRAHAM_JONES_COALESCENT_LIKELIHOOD
        unsigned nbranches = G::_nspecies*2 - 1;
        _log_coalescent_likelihood = 2 * nbranches * log(_forests[1]._theta_mean) - nbranches * boost::math::lgamma(2);
        
        vector<double> gamma_jb;
        vector<unsigned> q_jb;
        
        double neg_inf = -1*numeric_limits<double>::infinity();
        for (int i = 1; i<_forests.size(); i++) {
            
            _forests[i].refreshPreorder();
            
            if (_log_coalescent_likelihood != neg_inf) {
                pair<vector<double>, vector<unsigned>> params;
                params = _forests[i].calcInitialCoalescentLikelihoodIntegratingOutTheta();
                
                if (i == 1) {
                    for (auto &g:params.first) {
                        if (g == neg_inf) {
                            _log_coalescent_likelihood = neg_inf;
                            break;
                        }
                        gamma_jb.push_back(g);
                    }
                    for (auto &q:params.second) {
                        q_jb.push_back(q);
                    }
                }
                else {
                    for (unsigned p=0; p<params.first.size(); p++) {
                        gamma_jb[p] += params.first[p];
                        if (params.first[p] == neg_inf) {
                            _log_coalescent_likelihood = neg_inf;
                            break;
                        }
                        q_jb[p] += params.second[p];
                    }
                }
            }
        }
        return _log_coalescent_likelihood;
#else
        cerr << "calling calc initial coalescent likelihood when Graham Jones coalescent likelihood integrating out theta is not defined" << endl;
        assert (1 == 2); // only call this function when integrating out theta
        return 0;
#endif
    }
#endif

    inline void Particle::setNewTheta(bool fix_theta) {
        // gamma mean = shape * scale
        // draw mean from lognormal distribution
        // shape = 2.0 to be consistent with starbeast3
        // scale = 1 / mean;
        
        if (G::_theta_proposal_mean > 0.0) {
            assert (_theta_mean == 0.0);
            _theta_mean = _lot->gamma(1, G::_theta_proposal_mean); // equivalent to exponential(exponential_rate)
        }
        else {
            _theta_mean = G::_theta; // if no proposal distribution specified, use one theta mean for all particles
        }
        
        if (fix_theta) { // fix theta for all populations
            // map should be 2*nspecies - 1 size
            unsigned number = 0;
            vector<string> species_names;
            
            for (auto &s:_forests[1]._species_partition) {
                species_names.push_back(s.first);
                number++;
            }
            for (int i=0; i<G::_nspecies-1; i++) {
                string name = boost::str(boost::format("node-%d")%number);
                number++;
                species_names.push_back(name);
            }
            
            assert (species_names.size() == 2*G::_nspecies - 1);
            
            for (auto &name:species_names) {
                assert (G::_theta > 0.0);
                _theta_map[name] = G::_theta;      // create a theta map with all the same theta for simulations, set theta_mean to theta
            }
            
        }
        else {
            // create theta map
            createThetaMap();
        }
        
    }
    
    inline vector<double> Particle::getThetaVector() {
        vector<double> theta_vec;
        for (auto &t:_theta_map) {
            theta_vec.push_back(t.second);
        }
        return theta_vec;
    }

    inline double Particle::getPopMean() {
#if defined (DRAW_NEW_THETA)
        return _theta_mean;
#else
        return G::_theta;
#endif
    }

    inline void Particle::fixTheta() {
        createThetaMapFixedTheta();
        
#if defined (OLD_UPGMA)
        map<string, unsigned> species_indices = _forests[1]._species_indices;
#endif
        if (_forests.size() > 2) {
            for (int i=2; i<_forests.size(); i++) {
#if defined (OLD_UPGMA)
                _forests[i]._species_indices = species_indices;
#endif
            }
        }
    }

    inline void Particle::drawTheta() {
        // set seed first
//        assert (_psuffix > 0);
//        setSeed(rng.randint(1,9999) + _psuffix);
        
        // gamma mean = shape * scale
        // draw mean from lognormal distribution
        // shape = 2.0 to be consistent with starbeast3
        // scale = 1 / mean;
        
        if (G::_theta_proposal_mean > 0.0) {
            assert (_theta_mean == 0.0);
            _theta_mean = _lot->gamma(1, G::_theta_proposal_mean); // equivalent to exponential(exponential_rate)
        }
        else {
            _theta_mean = G::_theta; // if no proposal distribution specified, use one theta mean for all particles
        }
        
        assert (_theta_mean > 0.0);
            
        if (G::_fix_theta) {
            createThetaMapFixedTheta();
        }
        else {
            createThetaMap(); // create map for particle
        }
        
#if defined (OLD_UPGMA)
        map<string, unsigned> species_indices = _forests[1]._species_indices;
#endif
        if (_forests.size() > 2) {
            for (unsigned i=2; i<_forests.size(); i++) {
//                _forests[i]._theta_map = theta_map;
#if defined (OLD_UPGMA)
                _forests[i]._species_indices = species_indices;
#endif
//                _forests[i]._theta_mean = theta_mean;
            }
        }
    }

    inline void Particle::createThetaMap() {
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        for (auto &name:G::_species_names) {
            double new_theta = 0.0;
            if (new_theta < G::_small_enough) {
                new_theta = 1 / (_lot->gamma(2.0, scale));
                assert (new_theta > 0.0);
                _theta_map[name] = new_theta;
            }

        }
    }

    inline void Particle::createThetaMapFixedTheta() {
        for (auto &name:G::_species_names) {
            _theta_map[name] = G::_theta;
        }
    }

    inline void Particle::updateThetaMap(string new_species_name) {
        // add a new theta for the most recently drawn species
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        double new_theta = 0.0;
        if (new_theta < G::_small_enough) {
    //            new_theta = 1 / (lot->gamma(2.01, scale));
            new_theta = 1 / (_lot->gamma(2.0, scale));
            assert (new_theta > 0.0);
            _theta_map[new_species_name] = new_theta;
        }
    }

    inline void Particle::updateThetaMapFixedTheta(string new_species_name) {
        _theta_map[new_species_name] = G::_theta;
    }

#if !defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
    inline void Particle::resetThetaMap(Lot::SharedPtr lot, unordered_map<string, double> &theta_map) {
        assert (_theta_map.size() == 0);
        // map should be 2*nspecies - 1 size
        unsigned number = 0;
        vector<string> species_names;
        
        for (auto &s:_forests[1]._species_partition) {
            species_names.push_back(s.first);
            number++;
        }
        
        for (int i=0; i<G::_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            species_names.push_back(name);
        }
        
        assert (species_names.size() == 2*G::_nspecies - 1);
        
        // draw thetas for tips of species trees and ancestral population
        // for all other populations, theta = -1
        
        if (G::_theta_proposal_mean == 0.0) {
            assert (G::_theta > 0.0);
            G::_theta_proposal_mean = G::_theta;
        }
        double scale = 1 / G::_theta_proposal_mean;
        
        unsigned count = 0;
        for (auto &name:species_names) {
            if (count < G::_nspecies || count == 2*G::_nspecies-2) {
                double new_theta = 0.0;
                if (new_theta < G::_small_enough) {
                    new_theta = 1 / (lot->gamma(2.0, scale));
                    assert (new_theta > 0.0);
                    _theta_map[name] = new_theta;
                }
            }
            else {
                _theta_map[name] = -1;
            }
            count++;
        }
    }
#endif

    inline void Particle::drawNewTheta(string new_species) {
        // draw a new theta for the newest species population
        double scale = 1 / _theta_mean;
        double new_theta = 0.0;
        if (new_theta < G::_small_enough) {
            new_theta = 1 / _lot->gamma(2.0, scale);
            _theta_map[new_species] = new_theta;
        }
    }

    inline Particle::Particle(const Particle & other) {
        _lot.reset(new Lot());
        *this = other;
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

    inline void Particle::mapSpecies(map<string, string> &taxon_map) {
        //species tree
        _forests[0].setUpSpeciesForest();

        if (_forests[1]._lineages.size() > 0) { // don't redo this for faster second level when nodes have been cleared from gene trees
            //gene trees
            for (unsigned i=1; i<_forests.size(); i++) {
                _forests[i].setUpGeneForest(taxon_map);
            }
        }
    }

#if defined (DEBUG_MODE)
    inline void Particle::showSpeciesIncrement(){
        cout << "species tree increment: " << "     " << _forests[0]._last_edge_length << endl;
    }
#endif

#if defined (DEBUG_MODE)
    inline void Particle::showSpeciesJoined(){
        _forests[0].showSpeciesJoined();
    }
#endif
        
    inline int Particle::selectEventLinearScale(vector<double> weight_vec) {
        // choose a random number [0,1]
        double u =  _lot->uniform();
        assert (u > 0.0);
        assert (u < 1.0);
        double cum_prob = 0.0;
        unsigned index = 0;
        for (int i=0; i < (int) weight_vec.size(); i++) {
            cum_prob += weight_vec[i];
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

    inline vector<pair<double, double>> Particle::getSpeciesTreeIncrementPriors() {
        return _forests[0]._increments_and_priors;
    }

    inline double Particle::getCoalescentLikelihood(unsigned g) {
#if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
#if defined (FASTER_SECOND_LEVEL)
        return _log_coal_like;
#else
        return _log_coalescent_likelihood; // can't get coalescent likelihood separately for each gene tree
#endif
#else
        assert (g>0); // no coalescent likelihood for species tree
        return _forests[g]._log_coalescent_likelihood;
#endif
    }

    inline void Particle::simulateData(vector<unsigned> sites_vector) {
    // Simulate sequence data
        unsigned starting_site = 0;
        for (int i=1; i<_forests.size(); i++) {
            unsigned nsites = sites_vector[i-1];
            _forests[i].simulateData(_lot, starting_site, nsites);
            starting_site += sites_vector[i-1];
        }
    }

    inline vector<double> Particle::getThetaMap() {
        vector<double> thetas;
        for (auto &t:_theta_map) {
            thetas.push_back(t.second);
        }
        return thetas;
    }

    inline void Particle::processSpeciesNewick(string newick_string) {
        assert (newick_string != "");
       _t =  _forests[0].buildFromNewickMPI(newick_string, true, false, _lot);
//        _t = _forests[0].resetT();
        // TODO: need to set _t and _t_by_gene too
//        vector<tuple<string, string, string>> species_order = _forests[0].buildFromNewickTopology(newick_string);
//        _forests[0].resetLineages();
//        cout << "test";
//        _forests[0]._lineages.clear();
//        _species_order.erase(_species_order.begin()); // don't need "null", "null", "null"
    }

    inline string Particle::getTranslateBlock() {
        string block = "";
        block += "  Translate\n";
        unsigned count = 1;
        for (auto &nd:_forests[0]._nodes) {
            if (count < G::_nspecies + 1) {
                string name = nd._name;
                block += to_string(count) + " ";
                block += name;
                if (count != G::_nspecies) {
                    block += ",";
                    block += "\n";
                }
                else {
                    block += "\n";
                }
                count ++;
            }
            else {
                break;
            }
        }
        block += ";\n";
        return block;
    }

    inline void Particle::processGeneNewicks(vector<string> newicks) {
        for (int i=1; i<_forests.size(); i++) {
//            _forests[i]._nodes.clear();
            _forests[i].buildFromNewick(newicks[i-1], true, false); // newicks starts at 0
//            _forests[i].showForest();
            _forests[i].refreshPreorder();
//            _forests[i]._theta_mean = G::_theta; // for now, set theta mean equal to whatever theta user specifies
        }
        _theta_mean = G::_theta; // for now, set theta mean equal to whatever theta user specifies
    }

    inline void Particle::resetSpecies() {
        if (!G::_gene_newicks_specified) {
            _forests[0].clear();
        } // otherwise, starting from complete gene trees and species tree is already set
        setLogWeight(0.0);
        _log_coalescent_likelihood = 0.0;
#if !defined (FASTER_SECOND_LEVEL)
        _species_tree_height = 0.0;
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._log_coalescent_likelihood = 0.0;
            _forests[i]._data = nullptr;
            _forests[i]._log_weight = 0.0;
            // do not clear species indices - save this for use in jones coalescent likelihood calculation
        }
        _gene_order.clear();
        _next_species_number_by_gene.clear();
#if defined (UNUSED_FUNCTIONS)
        _starting_log_likelihoods.clear();
#endif
#endif
        _forests[0]._last_edge_length = 0.0;
        _forests[0]._increments_and_priors.clear();
        _t.clear();
        _t_by_gene.clear();
#if defined (FASTER_SECOND_LEVEL)
        _forests[0].refreshAllPreorders();
#endif
    }

    inline void Particle::setForest(Forest f, unsigned forest_number) {
        _forests[forest_number] = f;
    }

    inline void Particle::trimSpeciesTree() {
//        unordered_map<string, double> theta_map = _forests[1]._theta_map;

        unsigned spp_count = G::_nspecies*2 - 1;

        bool trim = true;
        vector<double> gene_tree_heights;
        for (unsigned i=1; i<_forests.size(); i++) {
            if (_forests[i]._species_partition.size() > 1) {
                gene_tree_heights.push_back(_forests[i]._forest_height);
            }
            else {
                trim = false;
                break;
            }
        }

        if (trim) {
            double max_gene_tree_height = *max_element(gene_tree_heights.begin(), gene_tree_heights.end());

            bool done = false;
            unsigned count = (unsigned) _t.size();

            while (!done) {
                double species_tree_height = _forests[0]._forest_height;
                double amount_to_trim = 0.0;

                Node* nd = _forests[0]._lineages.back();
                if (_forests[0]._lineages.size() < G::_nspecies) {
                    Node* subtree1 = nd->_left_child;
                    Node* subtree2 = nd->_left_child->_right_sib;

                    _forests[0].revertNodeVector(_forests[0]._lineages, subtree1, subtree2, nd);

                    // reset siblings and parents of original nodes back to 0
                    subtree1->resetNode(); //subtree1
                    subtree2->resetNode(); //subtree2

                    // clear new node from _nodes
                    //clear new node that was just created
                    nd->clear(); //new_nd

                    amount_to_trim = _t[count - 2].second;

                    _t.pop_back();
                    
                    for (auto &g:_t_by_gene) {
                        g.pop_back();
                    }

                    _forests[0]._ninternals--;

#if defined (DRAW_NEW_THETA)
                    _theta_map[G::_species_names[spp_count-1]] = -1.0;
#endif

                    spp_count--;
                }
                assert (amount_to_trim > 0.0);

                if (species_tree_height - amount_to_trim > max_gene_tree_height) {
                    for (auto &nd:_forests[0]._lineages) {
                        nd->_edge_length -= amount_to_trim;
                    }
                    _forests[0]._forest_height -= amount_to_trim;
                }
                else {
                    amount_to_trim = species_tree_height - max_gene_tree_height;
                    assert (amount_to_trim > 0.0);
                    for (auto &nd:_forests[0]._lineages) {
                        nd->_edge_length -= amount_to_trim;
                    }
                    _forests[0]._forest_height -= amount_to_trim;

                    _t[count-2].second -= amount_to_trim;

                    for (auto &g:_t_by_gene) {
                        g[count-2].second -= amount_to_trim;
                    }
                    done = true;
                }
                count--;
            }
        }

#if defined (DRAW_NEW_THETA)
//        for (unsigned i=1; i<_forests.size(); i++) {
//            _forests[i]._theta_map = theta_map; // TODO: make a theta_map for particle and pass to all forests?
//        }
#endif
    }

//    inline void Particle::trimSpeciesTree() {
//        map<string, double> theta_map = _forests[1]._theta_map;
//        vector<string> species_names = _forests[1]._species_names;
//
//        unsigned spp_count = (unsigned) species_names.size();
//
//        bool trim = true;
//        vector<double> gene_tree_heights;
//        for (unsigned i=1; i<_forests.size(); i++) {
//            if (_forests[i]._species_partition.size() > 1) {
//                gene_tree_heights.push_back(_forests[i]._forest_height);
//            }
//            else {
//                trim = false;
//                break;
//            }
//        }
//
//        if (trim) {
//            double max_gene_tree_height = *max_element(gene_tree_heights.begin(), gene_tree_heights.end());
//
//            bool done = false;
//            unsigned count = (unsigned) _t.size();
//
//            while (!done) {
//                double species_tree_height = _forests[0]._forest_height;
//                double amount_to_trim = 0.0;
//
//                Node* nd = _forests[0]._lineages.back();
//                if (_forests[0]._lineages.size() < G::_nspecies) {
//                    Node* subtree1 = nd->_left_child;
//                    Node* subtree2 = nd->_left_child->_right_sib;
//
//                    _forests[0].revertNodeVector(_forests[0]._lineages, subtree1, subtree2, nd);
//
//                    // reset siblings and parents of original nodes back to 0
//                    subtree1->resetNode(); //subtree1
//                    subtree2->resetNode(); //subtree2
//
//                    // clear new node from _nodes
//                    //clear new node that was just created
//                    nd->clear(); //new_nd
//
//                    _forests[0]._nodes.pop_back();
//
//                    amount_to_trim = _t[count - 2].second;
//
//                    _t.pop_back();
//
//                    for (auto &g:_t_by_gene) {
//                        g.pop_back();
//                    }
//                    _forests[0]._ninternals--;
//
//    //                theta_map[G::_species_names[spp_count-1]] = -1.0;
//
//                    spp_count--;
//                }
//
//                if (species_tree_height - amount_to_trim > max_gene_tree_height) {
//                    for (auto &nd:_forests[0]._lineages) {
//                        nd->_edge_length -= amount_to_trim;
//                    }
//                    _forests[0]._forest_height -= amount_to_trim;
//                }
//                else {
//                    amount_to_trim = species_tree_height - max_gene_tree_height;
//                    assert (amount_to_trim > 0.0);
//                    for (auto &nd:_forests[0]._lineages) {
//                        nd->_edge_length -= amount_to_trim;
//                    }
//                    _forests[0]._forest_height -= amount_to_trim;
//
//                    _t[count-2].second -= amount_to_trim;
//
//                    for (auto &g:_t_by_gene) {
//                        g[count-2].second -= amount_to_trim;
//                    }
//                    done = true;
//                }
//                count--;
//            }
//
//        }
//    //    for (unsigned i=1; i<_forests.size(); i++) {
//    //        _forests[i]._theta_map = theta_map;
//    //    }
//    }


#if defined (OLD_UPGMA)
    inline void Particle::createSpeciesIndices() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].createSpeciesIndices();
        }
    }
#endif

    inline void Particle::rebuildSpeciesTree() {
        bool trim_to_previous_join = false;
        
        if (trim_to_previous_join) {
            double min_branch_length = _forests[0]._lineages.back()->_edge_length; // species tree must remain at least as tall as it was after initial trimming
            for (auto &nd:_forests[0]._lineages) {
                nd->_edge_length -= min_branch_length;
            }

            _t.back().second -= min_branch_length;
            
            for (auto &g:_t_by_gene) {
                g.back().second -= min_branch_length;
            }
            
            tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
            
            double edge_increment = 0.0;
            while (edge_increment < min_branch_length) {
                // draw an increment and add to existing species lineages, don't join anything else at this stage
                edge_increment = _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0).first;
                if (edge_increment < min_branch_length) {
                    for (auto &nd:_forests[0]._lineages) {
                        nd->_edge_length -= edge_increment;
                    }
                }
            }
            
            assert (edge_increment >= min_branch_length);
            _t.back().second += edge_increment;
            
            for (auto &g:_t_by_gene) {
                g.back().second += edge_increment;
            }
            
            // now walk through loop,
            while (_forests[0]._lineages.size() > 1) {
                if (_forests[0]._lineages.size() > 1) {
                    species_joined = _forests[0].speciesTreeProposal(_lot);
                    double edge_len = 0.0;
                    if (_forests[0]._lineages.size() > 1) {
                        edge_len = _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0).first;
                    }
                    _t.push_back(make_pair(species_joined, edge_len));
                    
                    for (auto &g:_t_by_gene) {
                        g.push_back(make_pair(species_joined, edge_len));
                    }
                }
            }
            
            // update theta map
            for (auto &s:G::_species_names) {
                if (_theta_map[s] == -1.0) {
                    if (G::_fix_theta) {
                      updateThetaMapFixedTheta(s);
                    }
                    else {
                       updateThetaMap(s);
                    }
                }
            }
            
        }
        else {
            tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
            
            // draw an increment and add to existing species lineages, don't join anything else at this stage
            double edge_increment = _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0).first;
            _t.back().second += edge_increment;
            
            for (auto &g:_t_by_gene) {
                g.back().second += edge_increment;
            }
            
            // now walk through loop,
            while (_forests[0]._lineages.size() > 1) {
                if (_forests[0]._lineages.size() > 1) {
                    species_joined = _forests[0].speciesTreeProposal(_lot);
                    
                    double edge_len = 0.0;
                    if (_forests[0]._lineages.size() > 1) {
                        edge_len = _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0).first;
                    }
                    _t.push_back(make_pair(species_joined, edge_len));
                    
                    for (auto &g:_t_by_gene) {
                        g.push_back(make_pair(species_joined, edge_len));
                    }
                }
            }
            
#if defined (DRAW_NEW_THETA)
            // update theta map
            for (auto &s:G::_species_names) {
                if (_theta_map[s] == -1.0) {
                    updateThetaMap(s);
                }
            }
#endif
        }
 
    }

    inline void Particle::buildEntireSpeciesTree() {
        double edge_len = _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0).first;
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        _t.push_back(make_pair(species_joined, edge_len));

        for (unsigned i=0; i < G::_nspecies-1; i++) {
            if (_forests[0]._lineages.size() > 1) {
                species_joined = _forests[0].speciesTreeProposal(_lot);
                
                double edge_len = 0.0;
                if (_forests[0]._lineages.size() > 1) {
                    edge_len = _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0).first;
                }
                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline void Particle::setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species) {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].setNTaxaPerSpecies(ntaxa_per_species);
        }
    }

    inline unsigned Particle::getPartialCount() {
        unsigned partial_count = 0;
        for (auto &f:_forests) {
            partial_count += f._partials_calculated_count;
        }
        return partial_count;
    }

#if defined (FASTER_SECOND_LEVEL)
    inline void Particle::saveCoalInfoInitial() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].saveCoalInfoInitial();
        }
    }
#endif

#if defined (FASTER_SECOND_LEVEL)
    inline unsigned Particle::proposeSpeciationEvent() {
        // This function is only used for proposing speciation events when there are
        // complete gene trees available. It thus draws increments from a truncated
        // exponential distribution where the trunction point is the next height at
        // which at least one coalescent event combines lineages from two different
        // species.
        unsigned num_species_tree_lineages = 0;
        
        // Stores tuple (height, 0, vector of species) for each join in the current species forest.
        // Do not cap with ancestral species at this point.
        vector<Forest::coalinfo_t> sppinfo_vect;
        _forests[0].saveCoalInfoSpeciesTree(sppinfo_vect, false); // don't save ancestral pop

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
        // coalinfo_vect stores a tuple (height, gene + 1, vector of species)
        // for each join in any gene tree. To get 0-offset gene index,
        // subtract 1 from 2nd element of tuple (if 2nd element is 0, then
        // tuple represents a species tree join)
        vector<Forest::coalinfo_t> coalinfo_vect;
        
        // Add gene tree joins to coalinfo_vect
        // Just need coalescent events at this point in order to choose
        // limits for species tree increments
        
        for (unsigned f=1; f<_forests.size(); f++) {
            _forests[f].saveCoalInfoGeneForest(coalinfo_vect);
        }
        
        // Sort coalinfo_vect from smallest to largest height
        sort(coalinfo_vect.begin(), coalinfo_vect.end());
        
        // Get maximum height of any gene tree
        double max_height = get<0>((*coalinfo_vect.rbegin()));
        
//        if (G::_generation > 0) {
        if (_forests[0]._last_edge_length != 0.0) {

            // Create speciation event
            _forests[0].speciesTreeProposal(_lot);

            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _forests[0].buildCoalInfoVect();
            _forests[0].saveCoalInfoSpeciesTree(sppinfo_vect, false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());
                
            // Adjust elements of coalinfo_vect affected by species tree joins
            _forests[0].fixupCoalInfo(coalinfo_vect, sppinfo_vect);
        }
        
        // Draw a speciation increment dt.
        double forest_height = _forests[0]._forest_height;
        //speclog_element._height = forest_height;
        double h = findHeightNextCoalescentEvent(forest_height, coalinfo_vect);
        assert(h <= max_height);
        
        pair<double,double> tmp = _forests[0].chooseSpeciesIncrementOnlySecondLevel(_lot, h);
        
        double log_weight_factor = tmp.second;
        assert (log_weight_factor == log_weight_factor); // check for NaN
        
        num_species_tree_lineages = (unsigned) _forests[0]._lineages.size();
        
        if (num_species_tree_lineages == 2) {
            // Create final speciation event
            _forests[0].speciesTreeProposal(_lot);
        }
        
        // Add species tree joins to sppinfo_vect. Cap with ancestral species
        // in order to compute complete coalescent likelihood.
        sppinfo_vect.clear();
        _forests[0].buildCoalInfoVect();
        _forests[0].saveCoalInfoSpeciesTree(sppinfo_vect, /*cap*/true);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
                
        // Adjust elements of coalinfo_vect affected by species tree joins
        _forests[0].fixupCoalInfo(coalinfo_vect, sppinfo_vect);

        // Add speciations into coalinfo_vect
        coalinfo_vect.insert(coalinfo_vect.begin(), sppinfo_vect.begin(), sppinfo_vect.end());
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Compute coalescent likelihood and log weight
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
        _log_weight = _log_coal_like - _prev_log_coal_like + log_weight_factor;

        resetPrevLogCoalLike();
        return num_species_tree_lineages;
    }
#endif

#if defined (FASTER_SECOND_LEVEL)
    double Particle::findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect) {
        double min_height = G::_infinity;
        for (auto & ci : coalinfo_vect) {
            double h = get<0>(ci);
            if (h >= hstart) {
                vector<G::species_t> & v = get<2>(ci);
                assert(v.size() == 2);
                if (v[0] != v[1]) {
                    if (h < min_height)
                        min_height = h;
                    break;
                }
            }
        }
        
        assert(min_height < G::_infinity);
        return min_height;
    }
#endif

#if defined (FASTER_SECOND_LEVEL)
    inline double Particle::calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose) {
        assert (G::_taxon_names.size() > 0);
        // This function assumes gene forests are complete gene trees (not partial states) and
        // that preorders and heights have been precalculated.
        
        // The symbol b is used for a branch (i.e. segment) of the species tree to match the use of b
        // by Graham Jones (2017). Thus, a branch b is synonymous with a a leaf or ancestral species.
        
        // Uses eq. 4, p. 454, in G. Jones. 2017. Algorithmic improvements
        // to species delimitation and phylogeny estimation under the
        // multispecies coalescent. J. Math. Biol. 74:447-467.
        // doi: 10.1007/s00285-016-1034-0
        
        // eq. 4: log p(G | alpha, beta, sigma) = sum_b { log(r_b)
        //          + alpha log(sigma beta) - (alpha + q_b) log(sigma beta + gamma_b)
        //          + log Gamma(alpha + q_b) - log Gamma(alpha) }
        // where (eq. 2)
        //    p_j      = ploidy (2 = diploid, 1 = haploid) for gene j
        //    k_jb     = no. coal. events for gene j, species tree edge b
        //    c_jbi    = sojourn i to i+1 coalescence for gene j, edge b
        //    q_b      = sum_j k_jb
        //    log(r_b) = sum_j -k_jb log(p_j)
        //    gamma_b  = sum_j (1/p_j) sum_i^{k_jb} binom{n_jb - i}{2} c_jbi
        //
        // If p_j = 2 for all genes:
        //    log(r_b) = -log(2) q_b
        //    gamma_b  = 0.5 sum_j sum_i^{k_jb} binom{n_jb - i}{2} c_jbi
        //
        // For compatibility with starbeast3, we assume:
        //    alpha = 2
        //    beta  = (mean of theta)/4
        //    sigma = 1
        
    #if defined(DEBUG_COALLIKE)
        if (!coalinfo_vect.empty()) {
            for (auto & cinfo : coalinfo_vect) {
                double      height = get<0>(cinfo);
                unsigned     locus = get<1>(cinfo);
                auto & sppvect = get<2>(cinfo);
                vector<string> s;
                for (auto x : sppvect) {
                    s.push_back(to_string(x));
                }
            }
        }
    #endif
                
        // Create maps to hold quantities needed by Graham Jones' (2017) formula
        // for coalescent likelihood integrating out theta.
        typedef map<G::species_t, double> dmap;
        
        // n_jb[g][b] holds starting (leaf-end) number of lineages for locus g, branch b
        vector<dmap> n_jb(G::_nloci);
        
        // log_r_b[b] holds sum of log(r_b) for branch b
        dmap log_r_b;

        // q_b[b] holds number of coalescent events for branch b
        dmap q_b;

        // gamma_b[b] holds cumulative gamma_b for branch b
        dmap gamma_b;
        
        // Initialize n_jb for locus 0
        for (auto & nm : G::_taxon_names) {
            // Find index of species corresponding to taxon name nm
            unsigned i = G::_taxon_to_species.at(nm);
            
            // Get species from index
            G::species_t b = (G::species_t)1 << i;
            
            // Increment count of species b in locus 0
            if (n_jb[0].count(b) == 0)
                n_jb[0][b] = 1;
            else
                n_jb[0][b]++;
        }
        
        // Initialize n_jb for other loci by copying n_jb for locus 0
        // Assumes the same number of individuals have been sampled from each species for all loci
        for (unsigned g = 1; g < G::_nloci; g++) {
            n_jb[g] = n_jb[0];
        }
        
        // prev_height[g][b] holds previous height for locus g and branch b
        vector<map<G::species_t, double> > prev_height(G::_nloci);
        
        // Ploidy for each locus (currently all loci assumed to be diploid)
        vector<double> p_j(G::_nloci, 2.0);
        
        // Vector branches stores branches (including ancestral ones)
        vector<G::species_t> branches;
        for (unsigned i = 0; i < G::_nspecies; i++) {
            G::species_t b = (G::species_t)1 << i;
            branches.push_back(b);
        }
        
        // Vector bavail stores branches that are in the species tree at the current height
        set<G::species_t, bitless> bavail(branches.begin(), branches.end());
        
        // Walk down coalinfo_vect, accumulating Graham Jones r_b, q_b, and gamma_b
        for (auto & cinfo : coalinfo_vect) {
            double               height = get<0>(cinfo);
            unsigned             locus_plus_one   = get<1>(cinfo);
            int                  locus   = locus_plus_one - 1;
            
            vector<G::species_t> & b_vect = get<2>(cinfo);
            
            G::species_t banc = 0;
            for (auto b : b_vect)
                banc |= b;
                
            if (locus_plus_one == 0) {
                // Speciation event
                
                // If banc represents a new branch in the species tree,
                // add it to the vector of branches
                if (find(branches.begin(), branches.end(), banc) == branches.end()) {
                    branches.push_back(banc);
                }
                    
                for (unsigned g = 0; g < G::_nloci; g++) {
                    // Account for non-coalescence since the last coalescent
                    unsigned nsum = 0;
                    for (auto b : b_vect) {
                        if (prev_height[g].count(b) == 0) {
                            prev_height[g][b] = 0.0;
                        }
                        double nb = n_jb[g][b];
                        nsum += nb;
                        double cb = height - prev_height[g].at(b);
                        double gammab = 2.0*nb*(nb-1)*cb/p_j[g];
                        gamma_b[b] += gammab;
                    }
                    
                    // Record the beginning height of the new ancestral branch banc
                    if (prev_height[g].count(banc) == 0) {
                        prev_height[g][banc] = height;
                    }
                    
                    // Pool lineages from descendant species to create
                    // counts for the new ancestral species
                    n_jb[g][banc] = 0;
                    for (auto b : b_vect) {
                        n_jb[g].erase(b);
                    }
                    n_jb[g][banc] = nsum; //nleft + nright;
                    
                }
                      
                for (auto b : b_vect) {
                    bavail.erase(b);
                }
                bavail.insert(banc);
            }
            else {
                // Coalescent event

                // Determine b, the edge of the species tree we're dealing with
                // Note that b may be an edge that is not yet in the species tree (i.e. the root edge)
                auto it = find_if(bavail.begin(), bavail.end(), [banc](G::species_t v){return (banc & v) > 0;});
                assert(it != bavail.end());
                G::species_t b = *it;
                
                if (prev_height[locus].count(b) == 0) {
                    prev_height[locus][b] = 0.0;
                }
                double height0 = prev_height[locus].at(b);
                
                // Coalescent event: update Jones quantities
                log_r_b[b] += log(4/p_j[locus]);
                q_b[b]++;
                double n_jbk = n_jb[locus][b];
                double c_jbk = height - height0;
                gamma_b[b] += 2.0*n_jbk*(n_jbk-1)*c_jbk/p_j[locus];

                prev_height[locus][b] = height;
                n_jb[locus][b]--;
            }
        }

        double log_likelihood = 0.0;
        
        if (integrate_out_thetas) {
            double alpha = 2.0; // G::_invgamma_shape;
            double beta = G::_theta;
            assert(beta > 0.0);
            double B = (double)branches.size();
            log_likelihood  = B*alpha*log(beta);
            log_likelihood -= B*boost::math::lgamma(alpha);
            for (auto b : branches) {
                log_likelihood += log_r_b[b];
                log_likelihood -= (alpha + q_b[b])*log(beta + gamma_b[b]);
                log_likelihood += boost::math::lgamma(alpha + q_b[b]);
            }
        }
        else {
            double theta_b = G::_theta;
            assert(theta_b > 0.0);
            double sum_log_rb = 0.0;
            double sum_gamma_b = 0.0;
            unsigned sum_qb = 0;
            for (auto b : branches) {
                sum_qb += q_b[b];
                sum_log_rb += log_r_b[b];
                sum_gamma_b += gamma_b[b];
                log_likelihood += log_r_b[b] - gamma_b[b]/theta_b - q_b[b]*log(theta_b);
            }
        }

        _log_coal_like = log_likelihood;
        return log_likelihood;
    }
#endif

    inline void Particle::setNodeHeights() {
        for (unsigned f=1; f<_forests.size(); f++) {
            _forests[f].setNodeHeights();
        }
    }

#if defined (FASTER_SECOND_LEVEL)
    inline void Particle::resetPrevLogCoalLike() {
        _prev_log_coal_like = _log_coal_like;
    }
#endif

#if defined (FASTER_SECOND_LEVEL)
    inline void Particle::clearGeneForests() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].saveCoalInfoInitial();
            _forests[i].setTreeHeight();
            _forests[i]._data = nullptr;
            _forests[i]._nodes.clear();
            _forests[i]._lineages.clear();
            _forests[i]._preorder.clear();
#if !defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
            _forests[i]._vector_prior.clear();
#endif
        }
    }
#endif

#if defined (USING_MPI)
    inline void Particle::initSpeciesForest(string newick) {
        // TODO: need to rebuild _t and _t_by_gene first?
//        cout << "INITIALIZING SPECIES NEWICK generation " << _generation << "   " << newick << endl;
//        _forests[0].clear();
        _forests[0].buildFromNewickMPI(newick, true, false, _lot);
//        cout << "after initializing species forest, lineages are: " << endl;
//        cout << "\t";
//        for (auto &nd:_forests[0]._lineages) {
//            cout << nd->_name << " position " << nd->_position_in_lineages << endl;
//        }
//        _forests[0].showForest();
    }
#endif

#if defined (USING_MPI)
    inline void Particle::initGeneForest(string newick) {
        unsigned gene_number = _gene_order[_generation-1]; // updating for the previous step
//        _forests[gene_number].clear();
//        cout << "BUILDING NEWICK  " << newick << endl;
//        cout << "gene number is " << gene_number << endl;
        _forests[gene_number].buildFromNewickMPI(newick, true, true, _lot);
        
//        for (auto &nd:_forests[gene_number]._nodes) {
//            cout << nd._position_in_lineages << endl;
//        }
//        cout << "FOREST BUILT IS ";
//        _forests[gene_number].showForest();
        for (auto &s:_forests[gene_number]._species_partition) {
            cout << "x";
        }
    }
#endif

    inline void Particle::resetGeneOrder(unsigned step, vector<unsigned> gene_order) {
        if (_gene_order.size() == 0) {
            _gene_order.resize((G::_ntaxa-1)*G::_nloci);
        }
        unsigned count=0;
        if (step == 0) { // step 0 is different because this happens before any proposals have occurred
            for (unsigned s=step; s<G::_nloci; s++) {
                _gene_order[s] = gene_order[count];
                count++;
            }
        }
        else {
            for (unsigned s=step+1; s<step+G::_nloci+1; s++) {
                _gene_order[s] = gene_order[count];
                count++;
            }
        }
    }

#if defined (USING_MPI)
    inline void Particle::checkPartition() {
        for (auto &s:_forests[1]._species_partition) {
            cout << "x";
        }
    }
#endif

#if defined (USING_MPI)
    inline void Particle::initGeneForestSpeciesPartition(string species_partition) {
        unsigned gene_number = _gene_order[G::_generation-1];
        _forests[gene_number].resetSpeciesPartition(species_partition); // TODO: fill this in
    }
#endif

#if defined (USING_MPI)
    inline string Particle::saveForestSpeciesPartition() { // TODO: can't save by number because the order of nodes might change
        unsigned gene_number = _gene_order[_generation-1];
        map<string, vector<string>> partition = _forests[gene_number].saveSpeciesPartition();
        string partition_for_message = "";
        for (auto &m:partition) {
            partition_for_message += " spp " + m.first + "";
            for (auto &nd:m.second) {
                partition_for_message += " nd " + nd;
            }
        }
        return partition_for_message;
    }
#endif

    inline void Particle::setGeneUPGMAMatrices() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].setGeneUPGMAMatrices();
        }
    }

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _t = other._t;
        _psuffix = other._psuffix;
        _next_species_number = other._next_species_number;
        _t_by_gene = other._t_by_gene;
        _next_species_number_by_gene = other._next_species_number_by_gene;
        _gene_order = other._gene_order;
        _group_number = other._group_number;
        _theta_map = other._theta_map;
        _theta_mean = other._theta_mean;

        // the following data members apply only when simulating and do not need to be copied because simulating data only deals with one particle at a time
//            _num_deep_coalescences = other._num_deep_coalescences;
//            _max_deep_coal = other._max_deep_coal;

#if !defined (FASTER_SECOND_LEVEL)
        _species_branches = other._species_branches;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _species_tree_height = other._species_tree_height;
#endif
#if defined (FASTER_SECOND_LEVEL)
        _log_coal_like = other._log_coal_like;
        _prev_log_coal_like = other._prev_log_coal_like;
#endif
        
#if defined (UNUSED_FUNCTIONS)
        _starting_log_likelihoods = other._starting_log_likelihoods;
#endif
    };
}

