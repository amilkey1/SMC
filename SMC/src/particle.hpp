#pragma once
#include <vector>
#include "forest.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
#include <mutex>

using namespace std;
using namespace boost;

#include "lot.hpp"
#include "conditionals.hpp"

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
                                                    _data = d;
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
        void                                    mapSpecies(map<string, string> &taxon_map, vector<string> &species_names);
        void                                    saveForest(std::string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        double                                  getLogLikelihood();
        vector<double>                          calcGeneTreeLogLikelihoods();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
        double                                  getSpeciesIncrement () {return _forests[0]._last_edge_length;}
        double                                  getSpeciesLogWeight() const {return _log_species_weight;}
        unsigned                                getPartialCount();
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    setLogSpeciesWeight(double w){_log_species_weight = w;}
        void                                    setLogLikelihood(vector<double> forest_likelihoods);
        void                                    setLogCoalescentLikelihood(double coalescent_like);
        void                                    operator=(const Particle & other);
        const vector<Forest> &                  getForest() const {return _forests;}
        vector<double>                          getThetaMap();
        double                                  getThetaMean(){return _forests[1]._theta_mean;}
        string                                  saveForestNewick() {
            return _forests[0].makeNewick(8, true);}
        string                                  saveForestNewickAlt() {return _forests[0].makeAltNewick(8, false);}
            
        string                                  saveGeneNewick(unsigned i) {
            return _forests[i].makeNewick(8, true);}
    
        bool operator<(const Particle::SharedPtr & other) const {
            return _log_weight<other->_log_weight;
        }

        bool operator>(const Particle::SharedPtr & other) const {
            return _log_weight>other->_log_weight;
        }

        void                                            setGroupNumber(unsigned n) {_group_number = n;} // group number for parallelization
        unsigned                                        getGroupNumber() {return _group_number;}// group number for parallelization
        vector<Forest> &                                getForests() {return _forests;}
        void                                            showSpeciesIncrement();
        void                                            showSpeciesJoined();
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
        vector<double>                                  getVectorPrior();
        void                                            simulateData(vector<unsigned> sites_vector);
        unsigned                                        getNumDeepCoalescences() {return _num_deep_coalescences;}
        unsigned                                        getMaxDeepCoalescences(){return _max_deep_coal;}
        void                                            resetSpecies();
        void                                            setForest(Forest f, unsigned forest_number);
        Forest                                          getForest(unsigned i) {return _forests[i];} // TODO: should return a pointer?
        void                                            setNewTheta(bool fix_theta);
        vector<double>                                  getThetaVector();
        double                                          getPopMean();
        pair<string, string>                            getSpeciesJoined(){return make_pair(_forests[0]._species_joined.first->_name, _forests[0]._species_joined.second->_name);}
        void                                            setPsuffix(unsigned psuffix) {_psuffix = psuffix;}
        double                                          calcInitialCoalescentLikelihood();
        void                                            processGeneNewicks(vector<string> newicks);
        void                                            processSpeciesNewick(string newick_string);
        void                                            setNextSpeciesNumber() {_next_species_number = G::_nspecies;}
        string                                          getTranslateBlock();
        void                                            buildEntireSpeciesTree();
        void                                            rebuildSpeciesTree();
        void                                            setGeneOrder(vector<unsigned> gene_order) {_gene_order = gene_order;}
        void                                            trimSpeciesTree();
        void                                            setFixTheta(bool fix) {_fix_theta = fix;}
        void                                            setRelativeRatesByGene(vector<double> rel_rates);
        void                                            calcStartingUPGMAMatrix();
        vector<vector<double>>                          getStartingUPGMAMatrix();
        void                                            setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene);
        void                                            calcStartingRowCount();
        void                                            createSpeciesIndices();
        void                                            setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
#if defined (FASTER_SECOND_LEVEL)
        void                                            saveCoalInfoInitial();
        unsigned                                        proposeSpeciationEvent();
        double                                          findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect);
        double                                          calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose);
        void                                            resetPrevLogCoalLike();
        void                                            clearGeneForests();
#endif

    private:

#if defined (FASTER_SECOND_LEVEL)
        double                                  _log_coal_like;
        double                                  _prev_log_coal_like;
#endif
        vector<Forest>                          _forests;
        double                                  _log_weight;
        double                                  _log_species_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        int                                     _generation = 0;
        double                                  _log_coalescent_likelihood;
        mutable                                 Lot::SharedPtr _lot;
        unsigned                                _num_deep_coalescences;
        unsigned                                _max_deep_coal;
        double                                  _species_tree_height;
        unsigned                                _psuffix;
        unsigned                                _next_species_number;
        vector<unsigned>                        _next_species_number_by_gene;
        vector<tuple<string, string, string>>   _species_order;
        vector<pair<tuple<string, string, string>, double>> _t;
        vector<vector<pair<tuple<string, string, string>, double>>> _t_by_gene;
        vector<double>                          _starting_log_likelihoods;
        vector<unsigned>                        _gene_order;
        bool                                    _fix_theta;
        vector<double>                          _relative_rates_by_gene;
#if !defined (FASTER_SECOND_LEVEL)
        unsigned                                _species_branches;
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
        for (auto &t:_forests[1]._theta_map) {
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
        cout << " _log_species_weight: " << _log_species_weight << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _forests[0].showForest();
    }

    inline void Particle::clear() {
        _log_weight     = 0.0;
        _log_species_weight = 0.0;
        _log_likelihood = 0.0;
        _forests.clear();
        _data           = nullptr;
        _generation     = 0;
        _log_coalescent_likelihood = 0.0;
        _num_deep_coalescences = 0.0;
        _max_deep_coal = 0.0;
        _species_tree_height = 0.0;
        _t.clear();
        _psuffix = 0;
        _next_species_number = G::_nspecies;
        _species_order.clear();
        _starting_log_likelihoods.clear();
        _t_by_gene.clear();
        _next_species_number_by_gene.clear();
        _gene_order.clear();
        _fix_theta = false;
        _relative_rates_by_gene.clear();
#if !defined (FASTER_SECOND_LEVEL)
        _species_branches = 0;
#endif
        _group_number = 0;
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
        if (_generation == 0 && !G::_run_on_empty) {
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
        _log_species_weight = _log_coalescent_likelihood;
    }

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

    inline double Particle::getSpeciesTreeHeight() {
        return _forests[0].getTreeHeight();
    }

    inline double Particle::getSpeciesTreeLength() {
        return _forests[0].getTreeLength();
    }

    inline vector<double> Particle::getGeneTreeHeights() {
        vector<double> gene_tree_heights;
        for (int i=1; i<_forests.size(); i++) {
            gene_tree_heights.push_back(_forests[i].getTreeHeight());
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

inline vector<double> Particle::getVectorPrior() {
// this is the InverseGamma(2, psi) prior on the population sizes
    return _forests[1]._vector_prior;
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
        if (_generation == 0 && !G::_run_on_empty) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline void Particle::proposal() {
        double inv_gamma_modifier = 0.0;
//        double log_weight_modifier = 0.0;
        
        unsigned next_gene = _gene_order[_generation];
        bool calc_weight = false;
        
        if (_generation == 0) {
            buildEntireSpeciesTree();
            // make a separate species tree information vector for each gene
            for (unsigned i=1; i<_forests.size(); i++) {
                _t_by_gene.push_back(_t);
                _next_species_number_by_gene.push_back(0);
            }
        }
        else if (_generation % G::_nloci == 0 && G::_start_mode == "smc") { // after every locus has been filtered once, trim back the species tree as far as possible & rebuild it
            // don't rebuild the species tree at every step when simulating data
            trimSpeciesTree();
            if (_forests[0]._lineages.size() > 1) {
                rebuildSpeciesTree();
            }
        }
        
        bool done = false;
                
        while (!done) {
            vector<pair<double, string>> rates_by_species = _forests[next_gene].calcForestRate(_lot);
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
                           
                if (G::_start_mode == "smc") {
                    if (G::_upgma) {
                        if (!G::_run_on_empty) {
                            _forests[next_gene].buildRestOfTreeFaster();
                        }
                    }
                }
                    
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
                    
                    if (G::_start_mode == "sim") {
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
            

         if (_forests[1]._theta_mean == 0.0) {
             assert (G::_theta > 0.0);
             for (int i=1; i<_forests.size(); i++) {
                 _forests[i]._theta_mean = G::_theta;
             }
         }
         assert (_forests[1]._theta_mean > 0.0);
         
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
        if (_generation == 0 && G::_theta_prior_mean > 0.0 && G::_theta_proposal_mean > 0.0) {
            if (G::_theta_prior_mean != G::_theta_proposal_mean) {
                // else, log weight modifier is 0
                double prior_rate = 1.0/G::_theta_prior_mean;
                double proposal_rate = 1.0/G::_theta_proposal_mean;
                double log_weight_modifier = log(prior_rate) - log(proposal_rate) - (prior_rate - proposal_rate)*_forests[1]._theta_mean;

                _log_weight += log_weight_modifier;
            }
        }
#endif
        
        _generation++;
        
        if (G::_start_mode == "smc" && !G::_run_on_empty) {
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
                    max_depth -= _forests[0].getTreeHeight();
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
//            unsigned nbranches = Forest::_nspecies*2 - 1;
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
//        _forests[0].showForest();
//        _forests[1].showForest();
            _log_species_weight = _log_coalescent_likelihood - prev_log_coalescent_likelihood + constrained_factor;
#if !defined (UNCONSTRAINED_PROPOSAL)
            double test = 1/_log_species_weight;
            assert(test != -0); // assert coalescent likelihood is not -inf
#endif
        
        if (G::_run_on_empty && !G::_run_on_empty_first_level_only) {
            assert (_log_coalescent_likelihood == 0.0);
            assert (_log_species_weight == 0.0);
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
            _forests[1]._vector_prior.clear();
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
                    _forests[i]._ancestral_species_name = _forests[1]._ancestral_species_name;
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
                    max_depth -= _forests[0].getTreeHeight();
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
            _log_species_weight = _log_coalescent_likelihood - prev_log_coalescent_likelihood + constrained_factor;
#if !defined (UNCONSTRAINED_PROPOSAL)
            double test = 1/_log_species_weight;
            assert(test != -0); // assert coalescent likelihood is not -inf
#endif
        double neg_inf = -1*numeric_limits<double>::infinity();
        if (_log_coalescent_likelihood == neg_inf) {
            assert (_forests[0]._last_edge_length > max_depth);
        }
        
        _generation++;
#endif
    }

    inline void Particle::calcStartingUPGMAMatrix() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].buildStartingUPGMAMatrix(); // TODO: can do this once and copy to all particles
        }
    }

    inline vector<vector<double>> Particle::getStartingUPGMAMatrix() {
        vector<vector<double>> starting_upgma_matrices_by_gene;
        for (unsigned i=1; i<_forests.size(); i++) {
            starting_upgma_matrices_by_gene.push_back(_forests[i]._starting_dij);
        }
        return starting_upgma_matrices_by_gene;
    }

    inline void Particle::setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene) {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._starting_dij = starting_upgma_matrices_by_gene[i-1];
        }
    }

    inline void Particle::calcStartingRowCount() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].buildStartingRow();
        }
    }

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

    inline void Particle::setNewTheta(bool fix_theta) {
        if (fix_theta) { // fix theta for all populations
            // map should be 2*nspecies - 1 size
            unsigned number = 0;
            vector<string> species_names;
            map<string, double> theta_map;
            
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
                theta_map[name] = G::_theta;      // create a theta map with all the same theta for simulations, set theta_mean to theta
            }
            
            for (int i=1; i<_forests.size(); i++) {
                _forests[i]._theta_map = theta_map;
                _forests[i]._theta_mean = G::_theta;
            }
            
        }
        else {
            // create theta map
            _forests[1].createThetaMap(_lot);
            if (_forests.size() > 2) {
                for (int i=2; i<_forests.size(); i++) {
                    _forests[i]._theta_map = _forests[1]._theta_map;
                }
            }
        }
    }
    
    inline vector<double> Particle::getThetaVector() {
        vector<double> theta_vec;
        for (auto &t:_forests[1]._theta_map) {
            theta_vec.push_back(t.second);
        }
        return theta_vec;
    }

    inline double Particle::getPopMean() {
#if defined (DRAW_NEW_THETA)
        return _forests[1]._theta_mean;
#else
        return Forest::_theta;
#endif
    }

    inline void Particle::fixTheta() {
        _forests[1].createThetaMapFixedTheta(_lot); // create map for one forest, then copy it to all forests
        double theta_mean = _forests[1]._theta_mean;
        map<string, double> theta_map = _forests[1]._theta_map;
        map<string, unsigned> species_indices = _forests[1]._species_indices;
        if (_forests.size() > 2) {
            for (int i=2; i<_forests.size(); i++) {
                _forests[i]._theta_map = theta_map;
                _forests[i]._species_indices = species_indices;
                _forests[i]._theta_mean = theta_mean;
            }
        }
    }

    inline void Particle::drawTheta() {
        // set seed first
//        assert (_psuffix > 0);
//        setSeed(rng.randint(1,9999) + _psuffix);
            
        _forests[1].createThetaMap(_lot); // create map for one forest, then copy it to all forests
        double theta_mean = _forests[1]._theta_mean;
        map<string, double> theta_map = _forests[1]._theta_map;
        map<string, unsigned> species_indices = _forests[1]._species_indices;
        if (_forests.size() > 2) {
            for (int i=2; i<_forests.size(); i++) {
                _forests[i]._theta_map = theta_map;
                _forests[i]._species_indices = species_indices;
                _forests[i]._theta_mean = theta_mean;
            }
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

    inline void Particle::mapSpecies(map<string, string> &taxon_map, vector<string> &species_names) {
        //species tree
        _forests[0].setUpSpeciesForest(species_names);

        if (_forests[1]._lineages.size() > 0) {
            //gene trees
            for (unsigned i=1; i<_forests.size(); i++) {
                _forests[i].setUpGeneForest(taxon_map);
            }
        }
    }

    inline void Particle::showSpeciesIncrement(){
        cout << "species tree increment: " << "     " << _forests[0]._last_edge_length << endl;
    }

    inline void Particle::showSpeciesJoined(){
        _forests[0].showSpeciesJoined();
    }
        
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
        for (auto &t:_forests[1]._theta_map) {
            thetas.push_back(t.second);
        }
        return thetas;
    }

    inline void Particle::processSpeciesNewick(string newick_string) {
        assert (newick_string != "");
        _species_order = _forests[0].buildFromNewickTopology(newick_string);
        _forests[0]._lineages.clear();
        _species_order.erase(_species_order.begin()); // don't need "null", "null", "null"
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
            _forests[i].refreshPreorder();
            _forests[i]._theta_mean = G::_theta; // for now, set theta mean equal to whatever theta user specifies
        }
    }

    inline void Particle::resetSpecies() {
        _data = nullptr;
        _forests[0].clear();
        _generation = 0;
        setLogWeight(0.0);
        _species_tree_height = 0.0;
        _log_coalescent_likelihood = 0.0;
#if !defined (FASTER_SECOND_LEVEL)
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._log_coalescent_likelihood = 0.0;
            _forests[i]._data = nullptr;
            _forests[i]._log_weight = 0.0;
            // do not clear species indices - save this for use in jones coalescent likelihood calculation
        }
        _gene_order.clear();
        _next_species_number_by_gene.clear();
        _starting_log_likelihoods.clear();
        _t.clear();
        _t_by_gene.clear();
#endif
        _forests[0]._last_edge_length = 0.0;
        _forests[0]._increments_and_priors.clear();
#if defined (FASTER_SECOND_LEVEL)
        _forests[0].refreshAllPreorders();
#endif
    }

    inline void Particle::setForest(Forest f, unsigned forest_number) {
        _forests[forest_number] = f;
    }

    inline void Particle::trimSpeciesTree() {
        map<string, double> theta_map = _forests[1]._theta_map;
        vector<string> species_names = _forests[1]._species_names;
        
        unsigned spp_count = (unsigned) species_names.size();
        
        bool trim = true;
        vector<double> gene_tree_heights;
        for (unsigned i=1; i<_forests.size(); i++) {
            if (_forests[i]._species_partition.size() > 1) {
                gene_tree_heights.push_back(_forests[i].getTreeHeight());
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
                double species_tree_height = _forests[0].getTreeHeight();
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
                    
                    _forests[0]._nodes.pop_back();
                    
                    amount_to_trim = _t[count - 2].second;
                    
                    _t.pop_back();
                    
                    for (auto &g:_t_by_gene) {
                        g.pop_back();
                    }
                    _forests[0]._ninternals--;
                    
                    theta_map[species_names[spp_count-1]] = -1.0;
                    
                    spp_count--;
                }
                
                if (species_tree_height - amount_to_trim > max_gene_tree_height) {
                    for (auto &nd:_forests[0]._lineages) {
                        nd->_edge_length -= amount_to_trim;
                    }

                }
                else {
                    amount_to_trim = species_tree_height - max_gene_tree_height;
                    assert (amount_to_trim > 0.0);
                    for (auto &nd:_forests[0]._lineages) {
                        nd->_edge_length -= amount_to_trim;
                    }
                    _t[count-2].second -= amount_to_trim;
                    
                    for (auto &g:_t_by_gene) {
                        g[count-2].second -= amount_to_trim;
                    }
                    done = true;
                }
                count--;
            }

        }
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._theta_map = theta_map;
        }
    }


    inline void Particle::createSpeciesIndices() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].createSpeciesIndices();
        }
    }

    inline void Particle::rebuildSpeciesTree() {
        bool trim_to_previous_join = false;
        
        if (trim_to_previous_join) {
            map<string, double> theta_map = _forests[1]._theta_map;
                        
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
                _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
                edge_increment = _forests[0]._last_edge_length;
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
                        _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
                        edge_len = _forests[0]._last_edge_length;
                    }
                    _t.push_back(make_pair(species_joined, edge_len));
                    
                    for (auto &g:_t_by_gene) {
                        g.push_back(make_pair(species_joined, edge_len));
                    }
                }
            }
            
            // update theta map
            for (auto t:theta_map) {
                if (t.second == -1.0)  {
                    if (_fix_theta) {
                        _forests[1].updateThetaMapFixedTheta(_lot, t.first);
                    }
                    else {
                        _forests[1].updateThetaMap(_lot, t.first);
                    }
                }
            }
            
            if (_forests.size() > 2) {
                for (unsigned i=2; i<_forests.size(); i++) {
                    _forests[i]._theta_map = _forests[1]._theta_map;
                }
            }
            
        }
        else {
            map<string, double> theta_map = _forests[1]._theta_map;
            
            tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
            
            // draw an increment and add to existing species lineages, don't join anything else at this stage
            _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
            double edge_increment = _forests[0]._last_edge_length;
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
                        _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
                        edge_len = _forests[0]._last_edge_length;
                    }
                    _t.push_back(make_pair(species_joined, edge_len));
                    
                    for (auto &g:_t_by_gene) {
                        g.push_back(make_pair(species_joined, edge_len));
                    }
                }
            }
            
            // update theta map
            for (auto t:theta_map) {
                if (t.second == -1.0)  {
                    _forests[1].updateThetaMap(_lot, t.first);
                }
            }
            
            if (_forests.size() > 2) {
                for (unsigned i=2; i<_forests.size(); i++) {
                    _forests[i]._theta_map = _forests[1]._theta_map;
                }
            }
        }
 
    }

    inline void Particle::buildEntireSpeciesTree() {
        _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
        double edge_len = _forests[0]._last_edge_length;
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        _t.push_back(make_pair(species_joined, edge_len));

        for (unsigned i=0; i < G::_nspecies-1; i++) {
            if (_forests[0]._lineages.size() > 1) {
                species_joined = _forests[0].speciesTreeProposal(_lot);
                
                double edge_len = 0.0;
                if (_forests[0]._lineages.size() > 1) {
                    _forests[0].chooseSpeciesIncrementOnly(_lot, 0.0);
                    edge_len = _forests[0]._last_edge_length;
                }
                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline void Particle::setRelativeRatesByGene(vector<double> rel_rates) {
        _relative_rates_by_gene = rel_rates;
        
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._relative_rate = rel_rates[i-1];
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
//        cout << "x";
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

        // //temporary!
        // G::SpecLog speclog_element;
        // speclog_element._seed = seed;
        
        // Get maximum height of any gene tree
        double max_height = get<0>((*coalinfo_vect.rbegin()));
        
        if (_generation > 0) {
#if defined(DEBUG_COALLIKE)
//            output("\nSpecies tree before creating speciation:\n", 0);
//            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), 0);
#endif

            // Create speciation event
//            G::species_t left_spp, right_spp, anc_spp;
            _forests[0].speciesTreeProposal(_lot);
//            _forests[0].speciationEvent(_lot, left_spp, right_spp, anc_spp);
            
#if defined(DEBUG_COALLIKE)
//            output("\nSpecies tree after creating speciation:\n", 0);
//            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), 0);
#endif

            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _forests[0].buildCoalInfoVect();
            _forests[0].saveCoalInfoSpeciesTree(sppinfo_vect, false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            _species_forest.debugShowCoalInfo("sppinfo_vect after speciation event", sppinfo_vect, /*fn*/"");
            _species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
#endif
                
            // Adjust elements of coalinfo_vect affected by species tree joins
            _forests[0].fixupCoalInfo(coalinfo_vect, sppinfo_vect);

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            _species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");
#endif
                
            // //temporary!
            // speclog_element._left = left_spp;
            // speclog_element._right = right_spp;
            // speclog_element._anc = anc_spp;
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
#if defined(DEBUG_COALLIKE)
//            output("\nSpecies tree before creating FINAL speciation:\n", 0);
//            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), 0);
#endif
            // Create final speciation event
//            G::species_t left_spp, right_spp, anc_spp;
            _forests[0].speciesTreeProposal(_lot);
//            _forests[0].speciationEvent(_lot, left_spp, right_spp, anc_spp);
            
#if defined(DEBUG_COALLIKE)
//            output("\nSpecies tree after creating FINAL speciation:\n", 0);
//            output(format("  %s\n") % _species_forest.makeNewick(9, true, false), 0);
#endif
            // Let sppinfo_vect reflect current state of species forest
            //sppinfo_vect.clear();
            //_species_forest.buildCoalInfoVect();
            //_species_forest.saveCoalInfo(sppinfo_vect, /*cap*/false);
            
            // Sort sppinfo_vect from smallest height to largest height
            //sort(sppinfo_vect.begin(), sppinfo_vect.end());

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            //_species_forest.debugShowCoalInfo("sppinfo_vect after FINAL speciation event", sppinfo_vect, /*fn*/"");
            //_species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
#endif
            // Adjust elements of coalinfo_vect affected by species tree joins
            //_species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

#if defined(DEBUG_COALLIKE)
            // Show coalinfo_vect before fixing up
            //_species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");
#endif
        }
        
        // Add species tree joins to sppinfo_vect. Cap with ancestral species
        // in order to compute complete coalescent likelihood.
        sppinfo_vect.clear();
        _forests[0].buildCoalInfoVect();
        _forests[0].saveCoalInfoSpeciesTree(sppinfo_vect, /*cap*/true);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
        
#if defined(DEBUG_COALLIKE)
        // Show coalinfo_vect before fixing up
        _species_forest.debugShowCoalInfo("sppinfo_vect", sppinfo_vect, /*fn*/"");
        _species_forest.debugShowCoalInfo("coalinfo_vect before fixup", coalinfo_vect, /*fn*/"");
#endif
                
        // Adjust elements of coalinfo_vect affected by species tree joins
        _forests[0].fixupCoalInfo(coalinfo_vect, sppinfo_vect);

#if defined(DEBUG_COALLIKE)
        // Show coalinfo_vect after fixing up
        _species_forest.debugShowCoalInfo("coalinfo_vect after fixup", coalinfo_vect, /*fn*/"");
#endif

        // Add speciations into coalinfo_vect
        //BUG lines below fix 2nd level bug 2024-06-19
        // see also Particle::recordAllForests
        coalinfo_vect.insert(coalinfo_vect.begin(), sppinfo_vect.begin(), sppinfo_vect.end());
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Compute coalescent likelihood and log weight
        //double prev_log_coallike = _prev_log_coallike;
#if defined(DEBUG_COALLIKE)
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/true);
#else
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
#endif
        _log_species_weight = _log_coal_like - _prev_log_coal_like + log_weight_factor;

        resetPrevLogCoalLike();
        
        // //temporary!
        // speclog_element._logw = log_weight;
        // speclog_element._filtered = true;
        // G::_speclog[step].push_back(speclog_element);

        _generation++;
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
//            output("\nContents of coalinfo_vect:\n", 0);
//            output(format("%12s %6s %s\n") % "height" % "locus" % "species", 0);
            for (auto & cinfo : coalinfo_vect) {
                double      height = get<0>(cinfo);
                unsigned     locus = get<1>(cinfo);
                auto & sppvect = get<2>(cinfo);
                vector<string> s;
                for (auto x : sppvect) {
                    s.push_back(to_string(x));
                }
//                output(format("%12.9f %6d %s\n") % height % locus % boost::algorithm::join(s, ","), 0);
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
        
//        if (verbose) {
//            output("\nParticle::calcLogCoalescentLikelihood:\n", 1);
//            output("Key:\n", 1);
//            output("  height:   height above leaf-level of this coalescence or speciation event\n", 1);
//            output("  locus:    0 if speciation event, locus index for coalescence events\n", 1);
//            output("  spp:      species of the two lineages involved in a coalescence/speciation event\n", 1);
//            output("  b:        the species tree branch (i.e. species) to which this event contributes\n", 1);
//            output("  n_jbk:    the number of coalescences for branch b\n", 1);
//            output("  c_jbk:    the sojourn time prior to coalescence k\n", 1);
//            output("  log(r_b): the log of 4/pj, where pj is the ploidy of locus j\n", 1);
//            output("  gamma_b:  cumulative 2*n*(n-1)*c_jbk/pj over all loci for branch b\n", 1);
//            output("  q_b:      the cumulative number of coalescences over all loci for branch b\n", 1);
//            output(format("%12s %12s %25s %12s %12s %12s %12s %12s\n") % "height" % "locus" % "spp" % "b" % "n_jbk" % "c_jbk" % "gamma_b" % "q_b", 1);
//        }

        // Walk down coalinfo_vect, accumulating Graham Jones r_b, q_b, and gamma_b
        for (auto & cinfo : coalinfo_vect) {
            double               height = get<0>(cinfo);
            unsigned             locus_plus_one   = get<1>(cinfo);
            int                  locus   = locus_plus_one - 1;
            
            ostringstream oss;
            vector<G::species_t> & b_vect = get<2>(cinfo);
            copy(b_vect.begin(), b_vect.end(), ostream_iterator<G::species_t>(oss, " "));
            string spp_joined = oss.str();
            boost::algorithm::trim_right(spp_joined);
            
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

//                if (verbose) {
//                    output(format("%12.9f %12d %25s %12s %12s %12s %12s %12s\n") % height % 0 % spp_joined % banc % "-" % "-" % "-" % "-", 1);
//                }
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

//                if (verbose) {
//                    output(format("%12.9f %12d %25s %12d %12d %12.9f %12.9f %12d\n") % height % locus_plus_one % spp_joined % b % n_jbk % c_jbk % gamma_b[b] % q_b[b], 1);
//                }
            }
        }

//        if (verbose) {
//            output(format("\n%12s %12s %12s %12s %12s %12s\n") % "b" % "q_b" % "log(r_b)" % "gamma_b" % "theta_b" % "logL", 1);
//        }

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
                
//                if (verbose) {
//                    output(format("%12d %12d %12.9f %12.9f %12.9f %12.9f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % beta % log_likelihood, 1);
//                }
            }

//            if (verbose) {
//                output(format("\nalpha = %.9f\n") % alpha, 1);
//                output(format("beta  = %.9f\n") % beta, 1);
//                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood, 1);
//            }
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
                
//                if (verbose) {
//                    output(format("%12d %12d %12.9f %12.9f %12.9f %12.9f\n") % b % q_b[b] % log_r_b[b] % gamma_b[b] % theta_b % log_likelihood, 1);
//                }
            }
            
//            if (verbose) {
//                output(format("\nsum log(r_b) = %.9f\n") % sum_log_rb, 1);
//                output(format("sum q_b      = %d\n") % sum_qb, 1);
//                output(format("sum gamma_b  = %d\n") % sum_gamma_b, 1);
//                output(format("log(coalescent likelihood) = %.5f\n") % log_likelihood, 1);
//            }
        }

        _log_coal_like = log_likelihood;
        return log_likelihood;
    }
#endif

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
            _forests[i]._vector_prior.clear();
        }
    }
#endif

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_species_weight = other._log_species_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _data           = other._data;
        _generation     = other._generation;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _num_deep_coalescences = other._num_deep_coalescences;
        _max_deep_coal = other._max_deep_coal;
        _species_tree_height = other._species_tree_height;
        _t = other._t;
        _psuffix = other._psuffix;
        _next_species_number = other._next_species_number;
        _species_order = other._species_order;
        _starting_log_likelihoods = other._starting_log_likelihoods;
        _t_by_gene = other._t_by_gene;
        _next_species_number_by_gene = other._next_species_number_by_gene;
        _gene_order = other._gene_order;
        _fix_theta = other._fix_theta;
        _relative_rates_by_gene = other._relative_rates_by_gene;
#if !defined (FASTER_SECOND_LEVEL)
        _species_branches = other._species_branches;
#endif
        _group_number = other._group_number;
#if defined (FASTER_SECOND_LEVEL)
        _log_coal_like = other._log_coal_like;
        _prev_log_coal_like = other._prev_log_coal_like;
#endif
    };
}

