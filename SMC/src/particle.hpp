#pragma once
#include <vector>
#include "forest.hpp"
#include "species-forest.hpp"
#include "forest-extension.hpp"
#include "forest-func.hpp"
#include "boost/format.hpp"
#include "boost/math/special_functions/gamma.hpp"
#include <mutex>
#include <unordered_map>
#include "opvertex.hpp"

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
        void                                    proposalSim();
        void                                    setData(Data::SharedPtr d, map<string, string> &taxon_map, bool partials) {
                                                    int index = 1;
                                                    _prev_log_likelihoods.resize(G::_nloci);
                                                    _gene_forest_ptrs.resize(G::_nloci);
                                                    _gene_forest_extensions.resize(G::_nloci);
                                                    for (unsigned g = 0; g < G::_nloci; g++) {
                                                        _prev_log_likelihoods[g] = 0.0;
                                                        _gene_forest_ptrs[g] = Forest::SharedPtr(new Forest());
                                                        Forest::SharedPtr gfp = _gene_forest_ptrs[g];
                                                        gfp->setData(d, index, taxon_map, partials);
                                                        index++;
                                                    }
                                                }
        void                                    setSimData(Data::SharedPtr d, map<string, string> &taxon_map, unsigned ntaxa) {
                                                    int index = 1;
                                                    _gene_forests.resize(G::_nloci);
                                                    for (auto &_gene_forest:_gene_forests) {
                                                        _gene_forest.setSimData(d, index, taxon_map, ntaxa);
                                                        index++;
                                                    }
                                                }
        void                                    mapSpecies(map<string, string> &taxon_map);
        void                                    saveForest(std::string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        double                                  calcLogLikelihoodLocus(unsigned locus, bool join_finalized);
        double                                  getLogLikelihood();
        vector<double>                          calcGeneTreeLogLikelihoods();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
#if defined (DEBUG_MODE)
        double                                  getSpeciesIncrement () {return _species_forest._last_edge_length;}
#endif
        double                                  getSpeciesLogWeight() const {return _log_weight;}
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    setLogSpeciesWeight(double w){_log_weight = w;}
        void                                    setLogCoalescentLikelihood(double coalescent_like);
        void                                    operator=(const Particle & other);
        const vector<Forest> &                  getGeneForests() const {return _gene_forests;}
        const SpeciesForest &                   getSpeciesForest() const {return _species_forest;}
        vector<double>                          getThetaMap();
        double                                  getThetaMean(){return _theta_mean;}
        string                                  saveForestNewick() {
            return _species_forest.makeNewick(8, true);}
        string                                  saveForestNewickAlt() {return _species_forest.makeAltNewick(8, false);}
            
        string                                  saveGeneNewick(unsigned i);
        string                                  saveChangedForest() {return _gene_forests[_gene_order[G::_generation-1]-1].makePartialNewick(8, true);}
        unsigned                                getNextGene(){return _gene_order[G::_generation];}
        void                                    createThetaMap();
        void                                    createThetaMapFixedTheta();
        void                                    updateThetaMap(G::species_t new_species_name);
        void                                    updateThetaMapFixedTheta(G::species_t new_species_name);
    
        void                                    drawNewTheta(G::species_t new_species);
        void                                    resetSubgroupPointers(vector<Forest::SharedPtr> gene_forest_copies);
        void                                    resetSpeciesForestPointers(SpeciesForest species_forest_copy);
    
        void                                    proposeMCMCMove(bool last_round);
    
        bool operator<(const Particle::SharedPtr & other) const {
            return _log_weight<other->_log_weight;
        }

        bool operator>(const Particle::SharedPtr & other) const {
            return _log_weight>other->_log_weight;
        }

        void                                            setGroupNumber(unsigned n) {_group_number = n;} // group number for parallelization
        unsigned                                        getGroupNumber() {return _group_number;}// group number for parallelization
#if defined (DEBUG_MODE)
        void                                            showSpeciesJoined();
        void                                            showSpeciesIncrement();
        pair<string, string>                            getSpeciesJoined(){return make_pair(_species_forest._species_joined.first->_name, _species_forest._species_joined.second->_name);}
#endif
        void                                            showSpeciesTree();
        int                                             selectEventLinearScale(vector<double> weight_vec);
        double                                          getTopologyPrior(unsigned i);
        vector<pair<double, double>>                    getIncrementPriors(unsigned i);
        vector<pair<double, double>>                    getSpeciesTreeIncrementPriors();
        double                                          getCoalescentLikelihood(unsigned g);
        double                                          getCoalescentLikelihoodSecondLevel();
        void                                            clear();
        void                                            drawTheta();
        void                                            fixTheta();
        void                                            clearPartials();
        Lot::SharedPtr getLot() const {return _lot;}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}
        double                                          getSpeciesTreeHeight();
        double                                          getSpeciesTreeLength();
        double                                          getSpeciesTreeHeightAfterFirstRound(){return _species_forest._forest_height;}
        vector<double>                                  getGeneTreeHeights();
        vector<double>                                  getGeneTreeLengths();
        void                                            calcGeneTreeLengths();
        void                                            calcSpeciesTreeLength();
        vector<double>                                  getGeneTreeLogLikelihoods();
        vector<double>                                  getGeneTreePriors();
        inline vector<double>                           getGeneTreeCoalescentLikelihoods();
        double                                          getSpeciesTreePrior();
        double                                          getAllPriors();
        double                                          getAllPriorsFirstRound();
        void                                            simulateData(vector<unsigned> sites_vector);
        unsigned                                        getNumDeepCoalescences() {return _num_deep_coalescences;}
        unsigned                                        getMaxDeepCoalescences(){return _max_deep_coal;}
        void                                            resetSpecies();
        void                                            setNewTheta(bool fix_theta);
        vector<double>                                  getThetaVector();
        double                                          getPopMean();
        void                                            setPsuffix(unsigned psuffix) {_psuffix = psuffix;}
        void                                            processGeneNewicks(vector<string> newicks);
        void                                            processSpeciesNewick(string newick_string);
        string                                          getTranslateBlock();
        void                                            buildEntireSpeciesTree();
        void                                            buildEntireSpeciesTreeSim();
        void                                            rebuildSpeciesTree();
        void                                            setGeneOrder(vector<unsigned> gene_order) {_gene_order = gene_order;}
        void                                            resetGeneOrder(unsigned step, vector<unsigned> gene_order);
        void                                            trimSpeciesTree();
        void                                            setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
        void                                            saveCoalInfoInitial();
        unsigned                                        proposeSpeciationEvent();
        double                                          findHeightNextCoalescentEvent(double hstart, vector<Forest::coalinfo_t> & coalinfo_vect);
        double                                          calcLogCoalescentLikelihood(vector<Forest::coalinfo_t> & coalinfo_vect, bool integrate_out_thetas, bool verbose);
        void                                            clearGeneForests();
    
        void                                            setNodeHeights();
        void                                            setSortedThetaVector();
    
            vector<Forest::SharedPtr> &                 getGeneForestPtrs();
            const vector<Forest::SharedPtr> &           getGeneForestPtrsConst() const;
            Forest::SharedPtr                           getGeneForestPtr(unsigned locus);
            const Forest::SharedPtr                     getGeneForestPtrConst(unsigned locus) const;
            void                                        finalizeLatestJoin(int locus, unsigned index, map<const void *, list<unsigned> > & nonzero_map);
            void                                        finalizeLatestJoinMCMC(int locus, unsigned index);
            void                                        resetGeneForests(bool compute_partials);
            void                                        buildEnsembleCoalInfo();
            unsigned                                    getPartialCount();
            double                                      calcBHVDistance();
            vector<double>                              calcBHVDistanceGeneTrees();
            double                                      calcBHVDistanceTrueTree(string true_newick);
            double                                      calcBHVDistanceTrueTreeGene(string true_newick);
            double                                      opCalcLeafContribution(const Split::treeid_t & Alvs, const Split::treeid_t & Blvs) const;
            double                                      opFindCommonEdges(const Split::treeid_t & A, const Split::treeid_t & B, vector<Split> & common_edges) const;
            double                                      opCalcGeodesicDist(vector<Split::treeid_pair_t> & ABpairs) const;
            bool                                        opRefineSupport(const Split::treeid_pair_t & AB, Split::treeid_pair_t & AB1, Split::treeid_pair_t & AB2) const;
            static double                               opCalcTreeIDLength(const Split::treeid_t & splits);
            static void                                 opEdmondsKarp(op::OPVertex & source, op::OPVertex & sink, vector<op::OPVertex> & avect, vector<op::OPVertex> & bvect, Split::treeid_t & C1, Split::treeid_t & C2, Split::treeid_t & D1, Split::treeid_t & D2);
    
    private:
        SpeciesForest                           _species_forest;
        vector<Forest>                          _gene_forests;
        double                                  _log_weight;
        double                                  _log_likelihood;
        double                                  _log_coalescent_likelihood;
        mutable Lot::SharedPtr                  _lot;
        unsigned                                _num_deep_coalescences;
        unsigned                                _max_deep_coal;
        unsigned                                _psuffix;
        vector<unsigned>                        _next_species_number_by_gene;
        double                                  _prev_next_species_number_by_gene;
        double                                  _prev_increment_prior;
        unsigned                                _total_particle_partials = 0;
    
        vector<pair<tuple<G::species_t, G::species_t, G::species_t>, double>> _t;
        vector<vector<pair<tuple<G::species_t, G::species_t, G::species_t>, double>>> _t_by_gene;
        vector<pair<tuple<G::species_t, G::species_t, G::species_t>, double>> _prev_t_by_gene;
        vector<G::species_t>                    _sorted_species_names; // only use for simulations
        vector<Forest::coalinfo_t>              _ensemble_coalinfo;  // second level only
    
        vector<pair<tuple<string, string, string>, double>> _t_sim;
        vector<vector<pair<tuple<string, string, string>, double>>> _t_by_gene_sim;
    
        vector<unsigned>                        _gene_order;
    
        unordered_map<G::species_t, double>     _theta_map;
        double                                  _theta_mean;
        vector<double>                          _theta_vector;

        unsigned                                _group_number;
            mutable vector<ForestExtension>     _gene_forest_extensions;
            vector<Forest::SharedPtr>           _gene_forest_ptrs;
    
        vector<double>                          _prev_log_likelihoods;
        double                                  _prev_gene_increment;
        double                                  _prev_total_to_add;
        vector<G::species_t>                    _prev_species_assignments_before_coalescence;
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
        _species_forest.showForest();
        for (auto &_forest:_gene_forest_ptrs) {
            _forest->showForest();
        }
    }

    inline void Particle::showSpeciesParticle() {
        //print out weight of each species forest part of the particle
        cout << "\nParticle:\n";
        cout << " _log_species_weight: " << _log_weight << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _species_forest.showForest();
    }

    inline void Particle::clear() {
        _log_weight     = 0.0;
        _log_likelihood = 0.0;
        _species_forest.clear();
        _gene_forest_ptrs.clear();
        _prev_log_likelihoods.clear();
        _log_coalescent_likelihood = 0.0;
        _num_deep_coalescences = 0.0;
        _max_deep_coal = 0.0;
        _t.clear();
        _psuffix = 0;
        _t_by_gene.clear();
        _prev_t_by_gene.clear();
        _next_species_number_by_gene.clear();
        _gene_order.clear();
        _group_number = 0;
        _prev_gene_increment = 0;
        _prev_total_to_add = 0;
        _prev_next_species_number_by_gene = 0;
        _prev_increment_prior = 0;
        _prev_species_assignments_before_coalescence.clear();
    }

    inline void Particle::showSpeciesTree() {
        //print out weight of each particle
        cout << "\nParticle:\n";
        cout << "  _log_weight: " << _log_weight << "\n" ;
        cout << " _log_likelihood: " << _log_likelihood << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _species_forest.showForest();
    }

    //more detailed version of showParticle
    inline void Particle::debugParticle(std::string name) {
        cout << "debugging particle" << endl;
        //print out weight of each particle
        cout << "\nParticle " << name << ":\n";
        cout << "  _log_weight:               " << _log_weight                 << "\n" ;
        // species forest
        cout << "  _species_forest._ninternals:       " << _species_forest._ninternals         << "\n";
        cout << "  _species_forest._npatterns:        " << _species_forest._npatterns          << "\n";
        cout << "  _species_forest._last_edge_length: " << _species_forest._last_edge_length   << "\n";
        cout << "  newick description:        " << _species_forest.makeNewick(5,false) << "\n";
        
        // gene forests
        for (auto &_forest:_gene_forests) {
            cout << "  _log_likelihood:           " << _forest._gene_tree_log_likelihood << "\n";
            cout << "  _forest._ninternals:       " << _forest._ninternals         << "\n";
            cout << "  _forest._npatterns:        " << _forest._npatterns          << "\n";
            cout << "  _forest._last_edge_length: " << _forest._last_edge_length   << "\n";
            cout << "  newick description:        " << _forest.makeNewick(5,false) << "\n";
        }
    }

    inline double Particle::calcLogLikelihoodLocus(unsigned locus, bool join_finalized) {
        if (join_finalized) {
            return _gene_forest_ptrs[locus]->calcLogLikelihood();
        }
        else {
            return _gene_forest_ptrs[locus]->_gene_tree_log_likelihood;
//            return _prev_log_likelihoods[locus] + _gene_forest_extensions[locus].getLogWeight();
        }
    }

    inline double Particle::calcLogLikelihood() {
        //calculate likelihood for each gene tree
        double log_likelihood = 0.0;
        for (unsigned g = 0; g < G::_nloci; g++) {
            log_likelihood += _gene_forest_ptrs[g]->calcLogLikelihood();
        }
        if (G::_generation == 0 && !G::_run_on_empty) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline vector<double> Particle::calcGeneTreeLogLikelihoods() {
        vector<double> gene_forest_likelihoods;
        
        //calculate likelihood for each gene tree
        gene_forest_likelihoods.resize(_gene_forest_ptrs.size());
        for (unsigned i=0; i<_gene_forest_ptrs.size(); i++) {
            double gene_tree_log_likelihood  = 0.0;
            if (!G::_run_on_empty) {
                gene_tree_log_likelihood = _gene_forest_ptrs[i]->calcLogLikelihood();
                assert(!isnan (gene_tree_log_likelihood));
            }
            gene_forest_likelihoods[i] = gene_tree_log_likelihood;
        }

        return gene_forest_likelihoods;
    }

    inline void Particle::clearPartials() {
        for (auto &forest:_gene_forests) {
            forest.clearPartials();
        }
    }

    inline void Particle::setLogCoalescentLikelihood(double log_coal_like) {
        _log_coalescent_likelihood = log_coal_like;
        _log_weight = _log_coalescent_likelihood;
    }

    inline double Particle::getSpeciesTreeHeight() {
        return _species_forest._forest_height;
    }

    inline double Particle::getSpeciesTreeLength() {
        return _species_forest.getTreeLength();
    }

    inline vector<double> Particle::getGeneTreeHeights() {
        vector<double> gene_tree_heights;
        for (int i=0; i<_gene_forest_ptrs.size(); i++) {
            gene_tree_heights.push_back(_gene_forest_ptrs[i]->_forest_height);
        }
        if (G::_start_mode == "sim") {
            for (int i=0; i<_gene_forests.size(); i++) {
                gene_tree_heights.push_back(_gene_forests[i]._forest_height);
            }
        }
        return gene_tree_heights;
    }

    inline void Particle::calcGeneTreeLengths() {
        for (unsigned f=0; f<_gene_forest_ptrs.size(); f++) {
            _gene_forest_ptrs[f]->calcTreeLength();
        }
    }

    inline void Particle::calcSpeciesTreeLength() {
        _species_forest.calcTreeLength();
    }

    inline vector<double> Particle::getGeneTreeLengths() {
        vector<double> gene_tree_heights;
        for (int i=0; i<_gene_forest_ptrs.size(); i++) {
            if (_gene_forest_ptrs[i]->_forest_length == 0) {
                _gene_forest_ptrs[i]->calcTreeLength();
            }
            gene_tree_heights.push_back(_gene_forest_ptrs[i]->getTreeLength());
        }
        return gene_tree_heights;
    }

    inline vector<double> Particle::getGeneTreeLogLikelihoods() {
        vector<double> gene_tree_log_likelihoods;
        for (int i=0; i<_gene_forest_ptrs.size(); i++) {
            gene_tree_log_likelihoods.push_back(_gene_forest_ptrs[i]->_gene_tree_log_likelihood);
//            if (!G::_run_on_empty) {
//                assert (_gene_forest_ptrs[i]->_gene_tree_log_likelihood <= 0.0); // if there is all missing data for a gene, likelihood could be 0
//            }
        }
        return gene_tree_log_likelihoods;
    }

    inline vector<double> Particle::getGeneTreePriors() {
        vector<double> gene_tree_priors;
        for (int i=0; i<_gene_forest_ptrs.size(); i++) {
            double prior = 0.0;
            for (auto &p:_gene_forest_ptrs[i]->_increments_and_priors) {
                prior += p.second;
            }
            
            gene_tree_priors.push_back(prior);
        }
        return gene_tree_priors;
    }

    inline vector<double> Particle::getGeneTreeCoalescentLikelihoods() {
        // do not have separate coalescent likelihood for each gene tree
        vector<double> gene_tree_priors;
        for (int i=0; i<_gene_forest_ptrs.size(); i++) {
            gene_tree_priors.push_back(-1);
        }
        return gene_tree_priors;
    }

    inline double Particle::getSpeciesTreePrior() {
        double prior = 0.0;
        for (auto &p:_species_forest._increments_and_priors) {
            prior += p.second;
        }
        prior += _species_forest._log_joining_prob;
        
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

    inline double Particle::getLogLikelihood() {
        //retrieve likelihood for each gene tree
        double log_likelihood = 0.0;
        for (unsigned i=0; i<_gene_forest_ptrs.size(); i++) {
            double gene_tree_log_likelihood = _gene_forest_ptrs[i]->_gene_tree_log_likelihood;
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
        
//        _prev_log_likelihoods[next_gene-1] = _gene_forest_ptrs[next_gene-1]->getLogLikelihood();

        bool calc_weight = false;
        
        if (!G::_species_newick_specified) { // if specifying species newick, keep the species tree the same for all particles and never reset it in first level
            if (G::_generation == 0) {
                buildEntireSpeciesTree();
                // make a separate species tree information vector for each gene
                for (unsigned i=0; i<_gene_forest_extensions.size(); i++) {
                    _t_by_gene.push_back(_t);
                    _next_species_number_by_gene.push_back(0);
                }
            }
            else {
                trimSpeciesTree();
                if (_species_forest._lineages.size() > 1) {
                    rebuildSpeciesTree();
                }
            }
        }
        else {
            if (G::_generation == 0) {
                for (unsigned i=0; i<_gene_forest_extensions.size(); i++) {
                    _t_by_gene.push_back(_t);
                    _next_species_number_by_gene.push_back(0);
                }
            }
        }
        
        if (G::_mcmc) {
            _prev_t_by_gene = _t_by_gene[next_gene-1]; // preserve previous _t_by_gene
            
            _prev_species_assignments_before_coalescence = _gene_forest_ptrs[next_gene-1]->getAllPrevSpecies();
            
            _prev_next_species_number_by_gene = _next_species_number_by_gene[next_gene-1];
        }
        
        bool done = false;
        
        // Create temporary gene forest extending existing forest
        // without touching existing forest (which may be used
        // by many particles)
        assert(_gene_forest_extensions.size() == G::_nloci);
        _gene_forest_extensions[next_gene-1].dock(_gene_forest_ptrs[next_gene-1], _gene_forest_ptrs[next_gene-1]->pullPartial(), _lot);
                
        _prev_total_to_add = 0.0;
        _prev_increment_prior = 0.0;
        
        while (!done) {
            vector<pair<double, unsigned long>> rates_by_species = _gene_forest_extensions[next_gene-1].calcForestRate(_lot, _theta_map);
            double total_rate = 0.0;
            double gene_increment = -1.0;
            if (rates_by_species.size() > 0) {
                for (auto &r:rates_by_species) {
                    total_rate += r.first;
                }
                assert (total_rate > 0.0);
                
                // Draw coalescence increment delta ~ Exponential(total_rate)
                gene_increment = -log(1.0 - _lot->uniform())/total_rate;
                assert (gene_increment > 0.0);
            }
            
            unsigned next_species_index = _next_species_number_by_gene[next_gene-1];
            double species_increment = _t_by_gene[next_gene-1][next_species_index].second;
            
          // if total rate is 0, gene increment will be -1.0, which will be taken care of

            if ((gene_increment < species_increment || species_increment == 0.0) && gene_increment != -1.0) { // if species increment is 0.0, choose a coalescent event because the species tree is finished
                
                assert (gene_increment > 0.0);
                _prev_gene_increment = gene_increment;
                _prev_total_to_add += gene_increment;
                
                // tell gene forest extension about gene increment
                _gene_forest_extensions[next_gene-1].addIncrement(gene_increment);
                
                vector<double> event_choice_rates;
                for (auto &r:rates_by_species) {
                    event_choice_rates.push_back(r.first / total_rate);
                }
                
                unsigned index = selectEventLinearScale(event_choice_rates);
                G::species_t species_name = rates_by_species[index].second;
                
                _gene_forest_extensions[next_gene-1].coalesce(total_rate, species_name, _theta_map, true);
                _total_particle_partials++;
                
                    
                    if (species_increment > 0.0) { // otherwise, species tree is done and there is nothing left to update
                        _t_by_gene[next_gene-1][next_species_index].second -= gene_increment; // update species tree increments
                    }
                
                    _prev_increment_prior += log(total_rate) - (total_rate * gene_increment);
                
                    calc_weight = true;
                }
                else {
                    // carry out speciation event
                    
                    assert (species_increment > 0.0);
                
                    _prev_increment_prior -= (total_rate * species_increment);
                    
                    // tell gene forest extension about the species increment
                    _prev_total_to_add += species_increment;
                    _gene_forest_extensions[next_gene-1].addIncrement(species_increment);
                    
                    // tell gene forest extension about species merged
                    G::species_t left_spp = get<0>(_t_by_gene[next_gene-1][next_species_index+1].first);
                    G::species_t right_spp = get<1>(_t_by_gene[next_gene-1][next_species_index+1].first);
                    
                    if (left_spp == 0) {
                        assert (right_spp == 0);
                    }
                    else if (right_spp == 0) {
                        assert (left_spp == 0);
                    }
                    _gene_forest_extensions[next_gene-1].mergeSpecies(left_spp, right_spp);
                    
#if defined (DRAW_NEW_THETA)
                    if (!_theta_map.count(left_spp + right_spp)) { // don't overwrite theta map if it already got updated in the mcmc proposal
                        updateThetaMap(left_spp + right_spp);
                    }
#endif

                    assert (next_species_index < _t_by_gene[next_gene-1].size());
                    _t_by_gene[next_gene-1][next_species_index].second -= species_increment; // update species tree increments
                    assert (_t_by_gene[next_gene-1][next_species_index].second == 0.0);
                    
                    if (_gene_forest_extensions[next_gene-1].getSpeciesPartitionSize() > 1) {
                        _next_species_number_by_gene[next_gene-1]++;
                }
            }
                
        
            if (calc_weight) { // calc weight just means coalescent event has been proposed
                done = true;
            }
        }
        
//        _species_forest.showForest();
         
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
        
        if (!G::_run_on_empty) {
            double log_weight = _gene_forest_extensions[next_gene-1].getLogWeight();
            _log_weight = log_weight + inv_gamma_modifier;
        }
        else {
            _log_weight = 0.0;
        }
        
    }

    inline void Particle::proposalSim() {
        unsigned next_gene = _gene_order[G::_generation];
        
        bool calc_weight = false;
        
        if (G::_generation == 0) {
            _theta_map.clear();
            for (auto &t:G::_species_names_typed) {
                updateThetaMap(t);
            }
            buildEntireSpeciesTreeSim();
            // make a separate species tree information vector for each gene
            for (unsigned i=0; i<_gene_forests.size(); i++) {
                _t_by_gene_sim.push_back(_t_sim);
                _next_species_number_by_gene.push_back(0);
            }
        }
        
        bool done = false;
                
        while (!done) {
            vector<pair<double, string>> rates_by_species = _gene_forests[next_gene-1].calcForestRateSim(_lot, _theta_map);
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
            double species_increment = _t_by_gene_sim[next_gene-1][next_species_index].second;
            
          // if total rate is 0, gene increment will be -1.0, which will be taken care of

            if ((gene_increment < species_increment || species_increment == 0.0) && gene_increment != -1.0) { // if species increment is 0.0, choose a coalescent event because the species tree is finished
                
                assert (gene_increment > 0.0);
                _gene_forests[next_gene-1].addIncrement(gene_increment);
                
                vector<double> event_choice_rates;
                for (auto &r:rates_by_species) {
                    event_choice_rates.push_back(r.first / total_rate);
                }
                
                unsigned index = selectEventLinearScale(event_choice_rates);
                string species_name = rates_by_species[index].second;
                _gene_forests[next_gene-1].allowCoalescence(species_name, gene_increment, _lot);
                    
                if (species_increment > 0.0) { // otherwise, species tree is done and there is nothing left to update
                    _t_by_gene_sim[next_gene-1][next_species_index].second -= gene_increment; // update species tree increments
                }
                    calc_weight = true;
                }
                else {
                    // carry out speciation event
                    
                    assert (species_increment > 0.0);
                    assert (_gene_forests[next_gene-1]._species_partition.size() > 1);
                    _gene_forests[next_gene-1].addIncrement(species_increment);
                    
                    // need to tally up number of deep coalescences for simulations
                    _num_deep_coalescences += _gene_forests[next_gene-1].getDeepCoal(_t_by_gene_sim[next_gene - 1][next_species_index + 1].first);
                    _max_deep_coal += _gene_forests[next_gene-1].getMaxDeepCoal(_t_by_gene_sim[next_gene - 1][next_species_index + 1].first);
                    
                    _gene_forests[next_gene-1].updateSpeciesPartitionSim(_t_by_gene_sim[next_gene-1][next_species_index+1].first, _sorted_species_names[_next_species_number_by_gene[next_gene-1]]);
                    assert (next_species_index < _t_by_gene_sim[next_gene-1].size());
                    _t_by_gene_sim[next_gene-1][next_species_index].second -= species_increment; // update species tree increments
                    assert (_t_by_gene_sim[next_gene-1][next_species_index].second == 0.0);
                    if (_gene_forests[next_gene-1]._species_partition.size() > 1) {
                        _next_species_number_by_gene[next_gene-1]++;
                }
            }
                
        
            if (calc_weight) { // calc weight just means coalescent event has been proposed
                done = true;
            }
        }
         
         done = true;
        
        _log_weight = 0.0; // log weight is always 0 for simulations
    }

    struct bitless {
        bool operator()(const G::species_t a, const G::species_t b) const {
            bool returned_value = ((a & b) > 0 ? false : a < b);
            return returned_value;
        }
    };

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
            for (auto &name:G::_species_names_typed) {
                assert (G::_theta > 0.0);
                _theta_map[name] = G::_theta;      // create a theta map with all the same theta for simulations, set theta_mean to theta
            }
            
        }
        else {
            // create theta map
            createThetaMap();
        }
        
    }

    inline void Particle::setSortedThetaVector() {
        // _theta_map is unsorted, so get species names in the correct order
        for (auto &nd:_species_forest._nodes) {
            _theta_vector.push_back(_theta_map[nd._species]);
        }
    }
    
    inline vector<double> Particle::getThetaVector() {
        return _theta_vector;
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
    }

    inline void Particle::createThetaMap() {
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        for (auto &name:G::_species_names_typed) {
            double new_theta = 0.0;
            if (new_theta < G::_small_enough) {
                new_theta = 1 / (_lot->gamma(2.0, scale));
                assert (new_theta > 0.0);
                _theta_map[name] = new_theta;
            }
        }
        
    }

    inline void Particle::createThetaMapFixedTheta() {
        for (auto &name:G::_species_names_typed) {
            _theta_map[name] = G::_theta;
        }
    }

    inline void Particle::updateThetaMap(G::species_t new_species_name) {
        
        double new_theta = G::_theta;
        if (G::_start_mode_type == G::StartModeType::START_MODE_SIM && G::_fix_theta_for_simulations) {
            if (_theta_map.count(new_species_name) == 0) {
                _theta_map[new_species_name] = G::_theta; // only update theta map if the species does not already exist in the map
            }
        }
        else {
#if defined (DRAW_NEW_THETA)
        // add a new theta for the most recently drawn species
        double scale = (2.0 - 1.0) / _theta_mean;
        assert (scale > 0.0);
        new_theta = 0.0;
#endif
        if (new_theta < G::_small_enough) {
    //            new_theta = 1 / (lot->gamma(2.01, scale));
#if defined (DRAW_NEW_THETA)
            new_theta = 1 / (_lot->gamma(2.0, scale));
            assert (new_theta > 0.0);
#endif
            if (G::_fix_theta) {
                new_theta = G::_theta;
            }
            if (_theta_map.count(new_species_name) == 0) {
                _theta_map[new_species_name] = new_theta; // only update theta map if the species does not already exist in the map
            }
        }
        }
    }

        inline void Particle::updateThetaMapFixedTheta(G::species_t new_species_name) {
        _theta_map[new_species_name] = G::_theta;
    }

        inline void Particle::drawNewTheta(G::species_t new_species) {
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
        ofstream treef1(treefilename);
        treef1 << "#nexus\n\n";
        treef1 << "begin trees;\n";
        treef1 << "  tree test = [&R] " << _species_forest.makeNewick(8, true)  << ";\n";
        treef1 << "end;\n";
        treef1.close();
        
        for (auto &_forest:_gene_forests) {
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
        Node* base_node = _species_forest._lineages[0];
        sum_height += base_node->getEdgeLength();
        for (Node* child=base_node->_left_child; child; child=child->_left_child) {
            sum_height += child->getEdgeLength();
        }
        return sum_height;
    }

    inline void Particle::mapSpecies(map<string, string> &taxon_map) {
        //species tree
        _species_forest.setUpSpeciesForest();

        if (!G::_in_second_level) {
            //gene trees
            for (unsigned i=0; i<_gene_forests.size(); i++) {
                _gene_forests[i].setUpGeneForest(taxon_map);
            }
        }
    }

#if defined (DEBUG_MODE)
    inline void Particle::showSpeciesIncrement(){
        cout << "species tree increment: " << "     " << _species_forest._last_edge_length << endl;
    }
#endif

#if defined (DEBUG_MODE)
    inline void Particle::showSpeciesJoined(){
        _species_forest.showSpeciesJoined();
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
        // assuming i = gene number (starting at 1)
        assert (i > 0);
        return _gene_forests[i-1]._log_joining_prob;
    }

    inline vector<pair<double, double>> Particle::getIncrementPriors(unsigned i) {
        if (i == 0) {
            return _species_forest._increments_and_priors;
        }
        else {
            return _gene_forests[i-1]._increments_and_priors;
        }
    }

    inline vector<pair<double, double>> Particle::getSpeciesTreeIncrementPriors() {
        return _species_forest._increments_and_priors;
    }

    inline double Particle::getCoalescentLikelihood(unsigned g) {
        // first level coalescent likelihood
        return _gene_forest_ptrs[g]->_log_coalescent_likelihood;
//        return _log_coalescent_likelihood; // can't get coalescent likelihood separately for each gene tree // TODO: need to get coal like for first level
    }

    inline double Particle::getCoalescentLikelihoodSecondLevel() {
        return _log_coalescent_likelihood;
    }

    inline void Particle::simulateData(vector<unsigned> sites_vector) {
    // Simulate sequence data
        unsigned starting_site = 0;
        for (int i=0; i<_gene_forests.size(); i++) {
            unsigned nsites = sites_vector[i];
            _gene_forests[i].simulateData(_lot, starting_site, nsites);
            starting_site += sites_vector[i];
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
//        _t = _species_forest.resetT();
        // TODO: need to set _t and _t_by_gene too
//        vector<tuple<string, string, string>> species_order = _species_forest.buildFromNewickTopology(newick_string);
//        _species_forest.resetLineages();
//        cout << "test";
//        _species_forest._lineages.clear();
//        _species_order.erase(_species_order.begin()); // don't need "null", "null", "null"
    }

    inline string Particle::getTranslateBlock() {
        string block = "";
        block += "  Translate\n";
        unsigned count = 1;
        for (auto &nd:_species_forest._nodes) {
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
        for (int i=0; i<_gene_forest_ptrs.size(); i++) {
            _gene_forest_ptrs[i]->buildFromNewick(newicks[i], true, false); // newicks starts at 0
            _gene_forest_ptrs[i]->refreshPreorder();
            _gene_forest_ptrs[i]->setNodeHeights();
        }
        _theta_mean = G::_theta; // for now, set theta mean equal to whatever theta user specifies
        for (int i=0; i<_gene_forests.size(); i++) {
            _gene_forests[i].buildFromNewick(newicks[i], true, false); // newicks starts at 0
            _gene_forests[i].refreshPreorder();
            _gene_forests[i].setNodeHeights();
        }
        _theta_mean = G::_theta; // for now, set theta mean equal to whatever theta user specifies
    }

    inline void Particle::resetSpecies() {
        if (!G::_gene_newicks_specified) {
            _species_forest.clear();
        } // otherwise, starting from complete gene trees and species tree is already set
        setLogWeight(0.0);
        _log_coalescent_likelihood = 0.0;
        _species_forest._last_edge_length = 0.0;
        _species_forest._increments_and_priors.clear();
//        _t.clear();  // don't reset these because they will not get copied in the copier
//        _t_by_gene.clear();
        _species_forest.refreshAllPreorders();
        _species_forest._log_joining_prob = 0.0;
        for (auto &p:_gene_forest_ptrs) {
            p->_log_coalescent_likelihood = 0.0;
        }
    }

    inline void Particle::trimSpeciesTree() {
        unsigned spp_count = G::_nspecies*2 - 1;

        bool trim = true;
        vector<double> gene_tree_heights;
        
        for (unsigned i=0; i<_gene_forest_ptrs.size(); i++) {
            if (_gene_forest_ptrs[i]->checkNumberOfUniqueSpeciesInExistence() > 1) {
                gene_tree_heights.push_back(_gene_forest_ptrs[i]->_forest_height);
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
                double species_tree_height = _species_forest._forest_height;
                double amount_to_trim = 0.0;

                Node* nd = _species_forest._lineages.back();
                if (_species_forest._lineages.size() < G::_nspecies) {
                    Node* subtree1 = nd->_left_child;
                    Node* subtree2 = nd->_left_child->_right_sib;

                    _species_forest.revertNodeVector(_species_forest._lineages, subtree1, subtree2, nd);

                    // reset siblings and parents of original nodes back to 0
                    subtree1->resetSpeciesNode(); //subtree1
                    subtree2->resetSpeciesNode(); //subtree2

                    // clear new node from _nodes
                    //clear new node that was just created
                    nd->clear(); //new_nd
                    nd->resetSpeciesNode();

                    amount_to_trim = _t[count - 2].second;

                    _t.pop_back();
                    
                    for (auto &g:_t_by_gene) {
                        g.pop_back();
                    }

                    _species_forest._ninternals--;

                    spp_count--;
                }
                assert (amount_to_trim > 0.0);

                if (species_tree_height - amount_to_trim > max_gene_tree_height) {
                    for (auto &nd:_species_forest._lineages) {
                        nd->_edge_length -= amount_to_trim;
                    }
                    _species_forest._forest_height -= amount_to_trim;
                }
                else {
                    amount_to_trim = species_tree_height - max_gene_tree_height;
                    assert (amount_to_trim > 0.0);
                    for (auto &nd:_species_forest._lineages) {
                        nd->_edge_length -= amount_to_trim;
                    }
                    _species_forest._forest_height -= amount_to_trim;

                    _t[count-2].second -= amount_to_trim;

                    for (auto &g:_t_by_gene) {
                        g[count-2].second -= amount_to_trim;
                    }
                    done = true;
                }
                count--;
            }
        }
    }

    inline void Particle::rebuildSpeciesTree() {
        bool trim_to_previous_join = false;
        
        if (trim_to_previous_join) {
            double min_branch_length = _species_forest._lineages.back()->_edge_length; // species tree must remain at least as tall as it was after initial trimming
            for (auto &nd:_species_forest._lineages) {
                nd->_edge_length -= min_branch_length;
            }

            _t.back().second -= min_branch_length;
            
            for (auto &g:_t_by_gene) {
                g.back().second -= min_branch_length;
            }
            
            tuple<G::species_t, G::species_t, G::species_t> species_joined = make_tuple(0,0,0);
            
            double edge_increment = 0.0;
            while (edge_increment < min_branch_length) {
                // draw an increment and add to existing species lineages, don't join anything else at this stage
                edge_increment = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;
                if (edge_increment < min_branch_length) {
                    for (auto &nd:_species_forest._lineages) {
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
            while (_species_forest._lineages.size() > 1) {
                if (_species_forest._lineages.size() > 1) {
                    species_joined = _species_forest.speciesTreeProposal(_lot);
                    if (get<0>(species_joined) == 0) {
                        assert (get<1>(species_joined) == 0);
                    }
                    if (get<1>(species_joined) == 0) {
                        assert (get<0>(species_joined) == 0);
                    }
                    double edge_len = 0.0;
                    if (_species_forest._lineages.size() > 1) {
                        edge_len = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;
                    }
                    _t.push_back(make_pair(species_joined, edge_len));
                    
                    for (auto &g:_t_by_gene) {
                        g.push_back(make_pair(species_joined, edge_len));
                    }
                }
            }
            
            // update theta map
            for (auto &s:G::_species_names_typed) {
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
            tuple<G::species_t, G::species_t, G::species_t> species_joined = make_tuple(0,0,0);
            
            // draw an increment and add to existing species lineages, don't join anything else at this stage
            double edge_increment = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;
            _t.back().second += edge_increment;
            
            for (auto &g:_t_by_gene) {
                g.back().second += edge_increment;
            }
            
            // now walk through loop,
            while (_species_forest._lineages.size() > 1) {
                if (_species_forest._lineages.size() > 1) {
                    species_joined = _species_forest.speciesTreeProposal(_lot);
                    
                    if (get<0>(species_joined) == 0) {
                        assert (get<1>(species_joined) == 0);
                    }
                    else if (get<1>(species_joined) == 0) {
                        assert (get<0>(species_joined) == 0);
                    }
                    
                    double edge_len = 0.0;
                    if (_species_forest._lineages.size() > 1) {
                        edge_len = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;
                    }
                    _t.push_back(make_pair(species_joined, edge_len));
                    for (auto &g:_t_by_gene) {
                        g.push_back(make_pair(species_joined, edge_len));
                    }
                }
            }
            
#if defined (DRAW_NEW_THETA)
            // update theta map
            for (auto &s:G::_species_names_typed) {
                if (_theta_map[s] == -1.0) {
                    updateThetaMap(s);
                }
            }
#endif
        }
 
    }

    inline void Particle::buildEntireSpeciesTree() {
        double edge_len = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;
        
        tuple<G::species_t, G::species_t, G::species_t> species_joined = make_tuple(0,0,0);
        _t.push_back(make_pair(species_joined, edge_len));

        for (unsigned i=0; i < G::_nspecies-1; i++) {
            if (_species_forest._lineages.size() > 1) {
                species_joined = _species_forest.speciesTreeProposal(_lot);
                
                double edge_len = 0.0;
                if (_species_forest._lineages.size() > 1) {
                    edge_len = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;
                }
                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline void Particle::buildEntireSpeciesTreeSim() {
        double edge_len = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;

        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        _t_sim.push_back(make_pair(species_joined, edge_len));

        for (unsigned i=0; i < G::_nspecies-1; i++) {
            if (_species_forest._lineages.size() > 1) {
                species_joined = _species_forest.speciesTreeProposalSim(_lot);
                G::species_t new_species_name = _species_forest._lineages.back()->_left_child->_species + _species_forest._lineages.back()->_left_child->_right_sib->_species;
                updateThetaMap(new_species_name);
                _sorted_species_names.push_back(new_species_name);

                double edge_len = 0.0;
                if (_species_forest._lineages.size() > 1) {
                    edge_len = _species_forest.chooseSpeciesIncrementOnly(_lot, 0.0).first;
                }
                _t_sim.push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline void Particle::setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species) {
        for (unsigned i=0; i<_gene_forests.size(); i++) {
            _gene_forests[i].setNTaxaPerSpecies(ntaxa_per_species);
        }
    }

    inline void Particle::saveCoalInfoInitial() {
        for (unsigned i=0; i<_gene_forests.size(); i++) {
            _gene_forests[i].saveCoalInfoInitial();
        }
    }

    inline unsigned Particle::proposeSpeciationEvent() {
        double prev_log_coal_like = _log_coalescent_likelihood;
        
        // This function is only used for proposing speciation events when there are
        // complete gene trees available. It thus draws increments from a truncated
        // exponential distribution where the trunction point is the next height at
        // which at least one coalescent event combines lineages from two different
        // species.
        unsigned num_species_tree_lineages = 0;
        
        // Stores tuple (height, 0, vector of species) for each join in the current species forest.
        // Do not cap with ancestral species at this point.
        vector<Forest::coalinfo_t> sppinfo_vect;
        _species_forest.saveCoalInfoSpeciesTree(sppinfo_vect, false); // don't save ancestral pop

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
        
        assert (_ensemble_coalinfo.size() > 0);
        coalinfo_vect = _ensemble_coalinfo;
        
        // coalinfo_vect should already be sorted
        
        // Get maximum height of any gene tree
        double max_height = get<0>((*coalinfo_vect.rbegin()));
        
        if (_species_forest._last_edge_length != 0.0) { // don't draw a speciation event for first step

            // Create speciation event
            _species_forest.speciesTreeProposal(_lot);

            // Let sppinfo_vect reflect current state of species forest
            sppinfo_vect.clear();
            _species_forest.buildCoalInfoVect();
            _species_forest.saveCoalInfoSpeciesTree(sppinfo_vect, false);
            
            // Sort sppinfo_vect from smallest height to largest height
            sort(sppinfo_vect.begin(), sppinfo_vect.end());
                
            // Adjust elements of coalinfo_vect affected by species tree joins
            _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);
        }
        
        // Draw a speciation increment dt.
        double forest_height = _species_forest._forest_height;
        //speclog_element._height = forest_height;
        double h = findHeightNextCoalescentEvent(forest_height, coalinfo_vect);
        assert(h <= max_height);
        
        pair<double,double> tmp = _species_forest.chooseSpeciesIncrementOnlySecondLevel(_lot, h);
        
        double log_weight_factor = tmp.second;
        assert (log_weight_factor == log_weight_factor); // check for NaN
        
        num_species_tree_lineages = (unsigned) _species_forest._lineages.size();
        
        if (num_species_tree_lineages == 2) {
            // Create final speciation event
            _species_forest.speciesTreeProposal(_lot);
        }
        
        // Add species tree joins to sppinfo_vect. Cap with ancestral species
        // in order to compute complete coalescent likelihood.
        sppinfo_vect.clear(); // TODO: try overwriting this and not clearing?
        _species_forest.buildCoalInfoVect();
        _species_forest.saveCoalInfoSpeciesTree(sppinfo_vect, /*cap*/true);

        // Sort sppinfo_vect from smallest height to largest height
        sort(sppinfo_vect.begin(), sppinfo_vect.end());
                
        // Adjust elements of coalinfo_vect affected by species tree joins
        _species_forest.fixupCoalInfo(coalinfo_vect, sppinfo_vect);

        // Add speciations into coalinfo_vect
        coalinfo_vect.insert(coalinfo_vect.begin(), sppinfo_vect.begin(), sppinfo_vect.end());
        sort(coalinfo_vect.begin(), coalinfo_vect.end());

        // Compute coalescent likelihood and log weight
        calcLogCoalescentLikelihood(coalinfo_vect, /*integrate_out_thetas*/true, /*verbose*/false);
        _log_weight = _log_coalescent_likelihood - prev_log_coal_like + log_weight_factor;
        
        if (G::_run_on_empty) {
            _log_weight = 0.0;
        }

        return num_species_tree_lineages;
    }
           
    inline void Particle::buildEnsembleCoalInfo() {
        // Store a tuple (height, gene + 1, vector of species)
        // for each join in any gene tree in _ensemble_coalinfo.
        // To get 0-offset gene index,
        // subtract 1 from 2nd element of tuple (if 2nd element
        // is 0, then tuple represents a species tree join)
        _ensemble_coalinfo.clear();

        // Add gene tree joins to _ensemble_coalinfo_vect
        for (unsigned g = 0; g < G::_nloci; g++) {
            Forest::SharedPtr gfp = _gene_forest_ptrs[g];
            gfp->saveCoalInfoGeneForest(_ensemble_coalinfo);
        }
        
        sort(_ensemble_coalinfo.begin(), _ensemble_coalinfo.end());
    }

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
                for (auto & x : sppvect) {
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
            for (auto & b : b_vect)
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
                    for (auto & b : b_vect) {
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
                    for (auto & b : b_vect) {
                        n_jb[g].erase(b);
                    }
                    n_jb[g][banc] = nsum; //nleft + nright;
                    
                }
                      
                for (auto & b : b_vect) {
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
        
#if !defined (DRAW_NEW_THETA)
        integrate_out_thetas = false;
#endif
        
        if (G::_fix_theta) {
            integrate_out_thetas = false;
        }
        
        if (integrate_out_thetas) {
            double alpha = 2.0; // G::_invgamma_shape;
#if defined (DRAW_NEW_THETA)
            double beta = _theta_mean; // TODO: double check - G::_theta or _theta_mean?
#else
            double beta = G::_theta;
#endif
            assert(beta > 0.0);
            double B = (double)branches.size();
            log_likelihood  = B*alpha*log(beta);
            log_likelihood -= B*boost::math::lgamma(alpha);
            for (auto & b : branches) {
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
            for (auto & b : branches) {
                sum_qb += q_b[b];
                sum_log_rb += log_r_b[b];
                sum_gamma_b += gamma_b[b];
                log_likelihood += log_r_b[b] - gamma_b[b]/theta_b - q_b[b]*log(theta_b);
            }
        }

        _log_coalescent_likelihood = log_likelihood;
        return log_likelihood;
    }

    inline void Particle::setNodeHeights() {
        for (unsigned f=0; f<_gene_forests.size(); f++) {
            _gene_forests[f].setNodeHeights();
        }
    }

    inline void Particle::clearGeneForests() {
        for (unsigned i=0; i<_gene_forest_ptrs.size(); i++) {
            _gene_forest_ptrs[i]->saveCoalInfoInitial();
            _gene_forest_ptrs[i]->setTreeHeight(); // TODO: only need to do these things once per pointer
//            _gene_forest_ptrs[i]->_data = nullptr; // don't need to clear these because they will not get copied in the copier
//            _gene_forest_ptrs[i]->_nodes.clear();
//            _gene_forest_ptrs[i]->_lineages.clear(); // don't need clear these because they will not get copied in the copier
            _gene_forest_ptrs[i]->_preorder.clear();
        }
            
            buildEnsembleCoalInfo();
        }

    inline void Particle::resetGeneOrder(unsigned step, vector<unsigned> gene_order) {
        unsigned count=0;
        
        if (_gene_order.size() == 0) {
            _gene_order.resize((G::_ntaxa-1)*G::_nloci);
            
            assert (step == 0);
            // step 0 is different because this happens before any proposals have occurred
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
        inline void Particle::proposeMCMCMove(bool last_round) {
            // _prev_t_by_gene represents species / gene forest before the coalescent event we are attemping to change
            // save a copy of it to reset _prev_t_by_gene for the next mcmc round
            vector<pair<tuple<G::species_t, G::species_t, G::species_t>, double>> unchanged_prev_t_by_gene = _prev_t_by_gene;
            vector<G::species_t> unchanged_prev_species_assignments_before_coalescence = _prev_species_assignments_before_coalescence;
                    
            double unchanged_prev_next_species_number_by_gene = _prev_next_species_number_by_gene;
            
//            double hastings_ratio = 0.0;
            
            unsigned locus_number = _gene_order[G::_generation];
            
            // create a new gene forest extension to propose on
            ForestExtension new_gfx;
            new_gfx.dock(_gene_forest_ptrs[locus_number-1], _gene_forest_ptrs[locus_number-1]->pullPartial(), _lot);

            double u = _lot->uniform();
//            double prev_total_to_add = _prev_total_to_add;
            
            double prev_gene_tree_height = _gene_forest_ptrs[locus_number-1]->_forest_height;
            double prev_proposed_gene_tree_height = prev_gene_tree_height + _prev_total_to_add;
            
            double prev_species_tree_height = 0.0;
            
            for (unsigned i=0; i< _t.size(); i++) {
                if (_prev_t_by_gene[i].second == 0) {
                    prev_species_tree_height += _t[i].second;
                }
                else if (_prev_t_by_gene[i].second != _t[i].second) {
//                    prev_species_tree_height += _prev_t_by_gene[i].second; // include part of species tree that's been eaten into but not speciated
                    prev_species_tree_height += _t[i].second; // TODO: not sure if this is right - this is true if there is only 1 lineage per species, but if there are multiple lineages per species, you don't necessarily have to go beyond the species barrier? in that case, the boundary would be the previous coalescent event height?
                }
                else {
                    break;
                }
            }
            

            double proposed_height = (prev_proposed_gene_tree_height - G::_sliding_window / 2.0) + (u*G::_sliding_window);
            
            if (proposed_height < 0.0) {
                proposed_height *= -1; // reflect proposed height into the positive zone, then check if it violates the species tree
            }
            
            if ((proposed_height < prev_species_tree_height) || (proposed_height < prev_gene_tree_height)) { // if the chosen height eats into the fixed part of the species tree or gene trees, reflect it back into the legal region
                if (proposed_height < prev_species_tree_height) { // TODO: not sure if the reflection is correct
                    double amount_to_add = 2 * (prev_species_tree_height - proposed_height);
                    proposed_height += amount_to_add;
                }
                if (proposed_height < prev_gene_tree_height) {
                    double amount_to_add = 2 * (prev_gene_tree_height - proposed_height);
                    proposed_height += amount_to_add;
                }
//                proposed_height += G::_sliding_window / 2.0;
            }
            
            assert (proposed_height > prev_species_tree_height);
            
//            double proposed_increment = 0.0;
//            if (prev_species_tree_height > prev_gene_tree_height) {
//                proposed_increment = proposed_height - prev_species_tree_height; // proposed increment is the remaining height afer the specie tree
//            }
//            else {
                // there is some untouched gene tree after the species barrier
            double proposed_increment = proposed_height - prev_gene_tree_height; // the total increment to be added is the remainder of the height after the previous gene tree height
            assert (proposed_increment > 0);
//                if (proposed_increment < 0) {
//                    proposed_increment *= -1;
//                }
//            }
            
//            double proposed_increment = proposed_height - prev_species_tree_height; // proposed increment is the remaining height after the species tree
            
//            bool recalculate_increment = false;
            
//            double a = prev_gene_tree_height - G::_sliding_window / 2.0;
//            double new_sliding_window = 0.0;
//
//            if (a < prev_species_tree_height) {
//                a = prev_species_tree_height;
//                recalculate_increment = true;
//                // TODO: sliding window bounds change
//            }
//            double b = prev_gene_tree_height + G::_sliding_window / 2.0;
//            assert (b > 0.0);

            // b = proposed_increment + previous height of gene tree
            // a = previous height of gene tree - proposed increment or species tree height, whichever is smaller
            // a' to b' = window of new proposal
            
//            double proposed_increment = (prev_total_to_add - G::_sliding_window / 2.0) + (u*G::_sliding_window);
//
//            if (proposed_increment < 0) { // TODO: not sure if it's 0 anymore
//                proposed_increment *= -1;
//            }
            
            double total_proposed_increment = proposed_increment;
//            if (recalculate_increment) {
//                // make a new sliding window centered on proposed increment to calculate hastings ratio
//
//                double a_prime = prev_species_tree_height + proposed_increment - G::_sliding_window / 2.0;
//                assert (a_prime > 0.0);
//
//                double b_prime = prev_gene_tree_height + G::_sliding_window / 2.0;
//                assert (b_prime > 0.0);
//
//                hastings_ratio = (b - a) / (b_prime - a_prime);
//            }
            
//            double total_proposed_increment = proposed_increment;
//            double total_proposed_increment = new_proposed_increment;
            
            
            double new_increment_prior = 0.0;

            // new proposal
            bool done = false;
            bool calc_weight = false;
            bool impossible_increment = false;

            while (!done) {
                vector<pair<double, unsigned long>> rates_by_species = new_gfx.calcForestRate(_lot, _theta_map);

                double total_rate = 0.0;
                if (rates_by_species.size() > 0) {
                    for (auto &r:rates_by_species) {
                        total_rate += r.first;
                    }
                }

                unsigned next_species_index = _prev_next_species_number_by_gene;

                double species_increment = _prev_t_by_gene[next_species_index].second;
//                a += species_increment;

                if ((proposed_increment < species_increment || species_increment == 0.0) && ( rates_by_species.size() != 0)) {
                    // propose coalescence
                    new_gfx.addIncrement(proposed_increment);

                    vector<double> event_choice_rates;
                    for (auto &r:rates_by_species) {
                        event_choice_rates.push_back(r.first / total_rate);
                    }

                    unsigned index = selectEventLinearScale(event_choice_rates);
                    
                    assert (rates_by_species.size() > 0);
                    G::species_t species_name = rates_by_species[index].second;

                    new_gfx.coalesce(total_rate, species_name, _theta_map, false); // don't recalculate coal like
                    
                    _total_particle_partials++;

                    if (species_increment > 0.0) { // otherwise, species tree is done and there is nothing left to update
                        _prev_t_by_gene[next_species_index].second -= proposed_increment; // update species tree increments
                    }
                    
                    new_increment_prior += log(total_rate) - (total_rate * proposed_increment);
                    
                    calc_weight = true;
                }

                else {
                    // carry out speciation event
                    assert (species_increment > 0.0);

                    // add species increment to the gene tree
                    new_gfx.addIncrement(species_increment);

                    // merge the species in the gene tree
                    G::species_t left_spp = get<0>(_prev_t_by_gene[next_species_index+1].first);
                    G::species_t right_spp = get<1>(_prev_t_by_gene[next_species_index+1].first);

                    if (left_spp == 0) {
                        assert (right_spp == 0);
                    }
                    else if (right_spp == 0) {
                        assert (left_spp == 0);
                    }

                    new_gfx.mergeSpecies(left_spp, right_spp);

    #if defined (DRAW_NEW_THETA)
                    if (!_theta_map.count(left_spp + right_spp)) {
                        updateThetaMap(left_spp + right_spp);
                    }
    #endif
                    assert (next_species_index < _prev_t_by_gene.size());
                    _prev_t_by_gene[next_species_index].second -= species_increment; // update species tree increments
                    assert (_prev_t_by_gene[next_species_index].second == 0.0);

                    unsigned nremaining_species = new_gfx.getNRemainingSpecies();
                    if (nremaining_species > 1) {
                        _prev_next_species_number_by_gene++;
                    }
                    proposed_increment -= species_increment;
                    
                    if (rates_by_species.size() == 0) {
                        assert (total_rate == 0);
                    }
                    
                    new_increment_prior -= (total_rate * species_increment);
                    
                    if (proposed_increment <= 0.0) {
                        done = true;
                        impossible_increment = true;
                    }
                }

                if (calc_weight) { // calc weight just means coalescent event has been proposed
                    done = true;
                }
            }
            
            assert (!impossible_increment);
            
//            double hastings_ratio = 0;
//            if (impossible_increment) {
//                // TODO: propose a new increment
//                // a is species increment total, accumulated above
//                double b = total_proposed_increment;
//                double a_prime = 0;
//                double b_prime = (b - G::_sliding_window / 2) + (u*G::_sliding_window); // new proposed increment centered on old increment
//                if (b_prime < 0) {
//                    b_prime *= -1;
//                }
//                hastings_ratio = (b - a) / (b_prime - a_prime);
//                cout << "x";
//            }

            // calculate weight of new proposal
            double inv_gamma_modifier = 0.0;

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

            double new_log_weight = new_gfx.getLogWeight();
            new_log_weight += inv_gamma_modifier;
            
            double new_likelihood = _prev_log_likelihoods[locus_number-1]+new_gfx.getLogWeight();
            double prev_log_likelihood = _prev_log_likelihoods[locus_number-1]+_gene_forest_extensions[locus_number-1].getLogWeight();

            // accept or reject new move
            double prev_likelihood_x_prior = prev_log_likelihood + _prev_increment_prior;
            double proposed_likelihood_x_prior = new_likelihood + new_increment_prior;
            double log_ratio = proposed_likelihood_x_prior - prev_likelihood_x_prior;
            
            double u2 = log(_lot->uniform());

            bool accept = false;
            if (u2 < log_ratio && proposed_increment >= 0.0) {
                accept = true;
            }
            
            if (impossible_increment) {
                accept = false;
            }
            
//            cout << endl;
//            cout << prev_log_likelihood << "\t" << _prev_increment_prior << "\t" << prev_likelihood_x_prior << endl;
//            cout << new_likelihood << "\t" << new_increment_prior << "\t" << proposed_likelihood_x_prior << endl;
//            cout << "accept = " << accept << endl;
//            cout << endl;

            if (accept) {
                G::_nmcmc_moves_accepted++;
                _t_by_gene[locus_number-1] = _prev_t_by_gene;

                _next_species_number_by_gene[locus_number-1] = _prev_next_species_number_by_gene;
                _log_weight = new_log_weight; // only matters for keeping track of ESS / unique particles

                // _prev_log_likelihoods gets reset in the next proposal and is not used in this function
                
                _gene_forest_extensions[locus_number - 1] = new_gfx;
                
                if (!last_round) {
                    _prev_total_to_add = total_proposed_increment;
                    _prev_increment_prior = new_increment_prior;
                }
                // TODO: need to finalize latest join
            }
            
            // regardless of whether move is accepted, some variables need to be reset for the next round of mcmc
            if (!last_round) {
                _prev_species_assignments_before_coalescence = unchanged_prev_species_assignments_before_coalescence;
                _prev_t_by_gene = unchanged_prev_t_by_gene;
                _prev_next_species_number_by_gene = unchanged_prev_next_species_number_by_gene;
            }
            
            // if move is not accepted, nothing changes to the existing gene forest pointer
            
            // don't reset theta map - keep thetas the same until species tree is rebuilt
        }

    inline void Particle::finalizeLatestJoin(int locus, unsigned index, map<const void *, list<unsigned> > & nonzero_map) {
        // Makes join closest to leaf-level in _gene_forest_extensions[locus]
        // permanent, then undocks _gene_forest_extensions[locus]
        
        // Get reference to gene forest extension for this locus
        ForestExtension & gfx = _gene_forest_extensions[locus];
        
        // Get pointer to gene forest for this locus
        Forest::SharedPtr gfp = _gene_forest_ptrs[locus];
        
        // If we are not finalizing the last particle for this
        // gene forest object, make a copy that can be modified
        // without affecting other surviving particles
        unsigned nz = (unsigned)nonzero_map[gfp.get()].size();
        if (nz > 1) {
            // Remove the element corresponding to index
            list<unsigned> & v = nonzero_map[gfp.get()];
            auto it = find(v.begin(), v.end(), index);
            assert(it != v.end());
            v.erase(it);
            
            // Make a copy of the object pointed to by gfp
            Forest::SharedPtr gfcpy = Forest::SharedPtr(new Forest());
            *gfcpy = *gfp;
            _gene_forest_ptrs[locus] = gfcpy;
            
            // Let gpf point to the copy
            gfp = gfcpy;
        }
        
        // Copy log likelihood
//        gfp->setLogLikelihood(_prev_log_likelihoods[locus] + gfx.getLogWeight());
                        
        // Get splits for children of _proposed_anc
        const Node * anc = gfx.getProposedAnc();
        assert(anc);
        const Node * lchild = gfx.getProposedLChild();
        assert(lchild);
        const Node * rchild = gfx.getProposedRChild();
        assert(rchild);
        Split lsplit = lchild->_split;
        Split rsplit = rchild->_split;
        
        assert(anc->_split.isEquivalent(lsplit + rsplit));
        
        // Recreate extension's join in the actual gene forest
        double incr = gfx.getProposedDelta();
        assert(incr > 0.0);
        
        gfp->addIncrAndJoin(incr, lsplit, rsplit, gfx);
        
        _gene_forest_ptrs[locus]->_log_coalescent_likelihood = gfx.getLogCoalescentLikelihood();

        // Can now get rid of extension
        _gene_forest_extensions[locus].undock();
    }
        
    inline void Particle::finalizeLatestJoinMCMC(int locus, unsigned index) {
        // Makes join closest to leaf-level in _gene_forest_extensions[locus]
        // permanent, then undocks _gene_forest_extensions[locus]
        
        // Get reference to gene forest extension for this locus
        ForestExtension & gfx = _gene_forest_extensions[locus];
        
        // Get pointer to gene forest for this locus
        Forest::SharedPtr gfp = _gene_forest_ptrs[locus];
        
        // If we are not finalizing the last particle for this
        // gene forest object, make a copy that can be modified
        // without affecting other surviving particles
        
        // TODO: assuming every particle is unique here
        // Make a copy of the object pointed to by gfp
        Forest::SharedPtr gfcpy = Forest::SharedPtr(new Forest());
        *gfcpy = *gfp;
        _gene_forest_ptrs[locus] = gfcpy;
        
        // Let gpf point to the copy
        gfp = gfcpy;
        
        // Copy log likelihood
        gfp->setLogLikelihood(_prev_log_likelihoods[locus] + gfx.getLogWeight());
                        
        // Get splits for children of _proposed_anc
        const Node * anc = gfx.getProposedAnc();
        assert(anc);
        const Node * lchild = gfx.getProposedLChild();
        assert(lchild);
        const Node * rchild = gfx.getProposedRChild();
        assert(rchild);
        Split lsplit = lchild->_split;
        Split rsplit = rchild->_split;
        
        assert(anc->_split.isEquivalent(lsplit + rsplit));
        
        // Recreate extension's join in the actual gene forest
        double incr = gfx.getProposedDelta();
        assert(incr > 0.0);
        
        gfp->addIncrAndJoin(incr, lsplit, rsplit, gfx);

        // Can now get rid of extension
        _gene_forest_extensions[locus].undock();
    }

    inline Forest::SharedPtr Particle::getGeneForestPtr(unsigned locus) {
        assert(locus < G::_nloci);
        assert(_gene_forest_ptrs[locus]);
        return _gene_forest_ptrs[locus];
    }
        
    inline unsigned Particle::getPartialCount() {
        if (G::_plus_G) {
            _total_particle_partials *= G::_gamma_rate_cat.size();
        }
        return _total_particle_partials;
    }
        
    inline string Particle::saveGeneNewick(unsigned i) {
        if (G::_start_mode_type == G::StartModeType::START_MODE_SIM) {
            return _gene_forests[i-1].makeNewick(8, true);
        }
        else {
            return _gene_forest_ptrs[i-1]->makeNewick(8, true);
        }
    }
        
    inline void Particle::resetSubgroupPointers(vector<Forest::SharedPtr> gene_forest_copies) {
        for (unsigned l=0; l<G::_nloci; l++) {
            Forest::SharedPtr gfp = _gene_forest_ptrs[l];
            Forest::SharedPtr gfcpy = gene_forest_copies[l];
            
            *gfcpy = *gfp;
            _gene_forest_ptrs[l] = gfcpy;
            
            // Let gpf point to the copy
            gfp = gfcpy;
        }
    }

    inline void Particle::resetSpeciesForestPointers(SpeciesForest species_forest_copy) {
        SpeciesForest sfp = _species_forest;
        SpeciesForest sfcpy = species_forest_copy;
        
        sfcpy = sfp;
        _species_forest = species_forest_copy;
    }
        
    inline double Particle::calcBHVDistanceTrueTreeGene(string true_newick) {
        // Build the reference tree
        Forest reftree;
        
        reftree.buildFromNewick(G::_bhv_reference, true, false);
        
        Forest truetree;
        truetree.buildFromNewick(true_newick, true, false);

        // Store splits from the reference tree
        Split::treeid_t A0;
        Split::treeid_t Alvs;
        reftree.storeSplits(A0, Alvs);

        // Only save splits with non-zero edge length
        Split::treeid_t A;
        for (auto & a : A0) {
            if (a.getEdgeLen() > 0.0) {
                A.insert(a);
            }
        }

        // test tree is the species tree
        // store splits from the test tree
        Split::treeid_t B0;
        Split::treeid_t Blvs;
        truetree.storeSplits(B0, Blvs);

        // Only save splits with non-zero edge length
        Split::treeid_t B;
        for (auto & b : B0) {
            if (b.getEdgeLen() > 0.0) {
                B.insert(b);
            }
        }

        // Calculate the contribution of leaf edges to the geodesic
        double leaf_contribution_squared = opCalcLeafContribution(Alvs, Blvs);

        // Find common edges and calculate the contribution of common edge lengths to the geodesic
        vector<Split> common_edges;
        double common_edge_contribution_squared = opFindCommonEdges(A, B, common_edges);

        if (!common_edges.empty()) {
            // Remove splits representing common edge from both A and B
            for (auto & c : common_edges) {
                A.erase(c);
                B.erase(c);
            }
        }
        vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
        in_pairs.emplace_back(A,B);

        unsigned pair_index = 1;
        vector<double> geodesic_distances;
        for (auto & inpair : in_pairs) {

            vector<Split::treeid_pair_t> ABpairs;
            ABpairs.push_back(inpair);

            double L = opCalcGeodesicDist(ABpairs);

            geodesic_distances.push_back(L);
            ++pair_index;
        }

        // Calculate total geodesic distance
        double total_geodesic_distance = 0.0;
        for (double x : geodesic_distances) {
            total_geodesic_distance += pow(x, 2);
        }
        total_geodesic_distance += leaf_contribution_squared;
        total_geodesic_distance += common_edge_contribution_squared;
        total_geodesic_distance = sqrt(total_geodesic_distance);

        return total_geodesic_distance;
    }
        
    inline double Particle::calcBHVDistanceTrueTree(string true_newick) {
        // Build the reference tree
        SpeciesForest reftree;
        
        reftree.buildFromNewick(G::_bhv_reference, true, false);
        
        SpeciesForest truetree;
        truetree.buildFromNewick(true_newick, true, false);

        // Store splits from the reference tree
        Split::treeid_t A0;
        Split::treeid_t Alvs;
        reftree.storeSplits(A0, Alvs);

        // Only save splits with non-zero edge length
        Split::treeid_t A;
        for (auto & a : A0) {
            if (a.getEdgeLen() > 0.0) {
                A.insert(a);
            }
        }

        // test tree is the species tree
        // store splits from the test tree
        Split::treeid_t B0;
        Split::treeid_t Blvs;
        truetree.storeSplits(B0, Blvs);

        // Only save splits with non-zero edge length
        Split::treeid_t B;
        for (auto & b : B0) {
            if (b.getEdgeLen() > 0.0) {
                B.insert(b);
            }
        }

        // Calculate the contribution of leaf edges to the geodesic
        double leaf_contribution_squared = opCalcLeafContribution(Alvs, Blvs);

        // Find common edges and calculate the contribution of common edge lengths to the geodesic
        vector<Split> common_edges;
        double common_edge_contribution_squared = opFindCommonEdges(A, B, common_edges);

        if (!common_edges.empty()) {
            // Remove splits representing common edge from both A and B
            for (auto & c : common_edges) {
                A.erase(c);
                B.erase(c);
            }
        }
        vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
        in_pairs.emplace_back(A,B);

        unsigned pair_index = 1;
        vector<double> geodesic_distances;
        for (auto & inpair : in_pairs) {

            vector<Split::treeid_pair_t> ABpairs;
            ABpairs.push_back(inpair);

            double L = opCalcGeodesicDist(ABpairs);

            geodesic_distances.push_back(L);
            ++pair_index;
        }

        // Calculate total geodesic distance
        double total_geodesic_distance = 0.0;
        for (double x : geodesic_distances) {
            total_geodesic_distance += pow(x, 2);
        }
        total_geodesic_distance += leaf_contribution_squared;
        total_geodesic_distance += common_edge_contribution_squared;
        total_geodesic_distance = sqrt(total_geodesic_distance);

        return total_geodesic_distance;
    }
        
    inline double Particle::calcBHVDistance() {
        // Build the reference tree
        SpeciesForest reftree;
        
        reftree.buildFromNewick(G::_bhv_reference, true, false);

        // Store splits from the reference tree
        Split::treeid_t A0;
        Split::treeid_t Alvs;
        reftree.storeSplits(A0, Alvs);

        // Only save splits with non-zero edge length
        Split::treeid_t A;
        for (auto & a : A0) {
            if (a.getEdgeLen() > 0.0) {
                A.insert(a);
            }
        }

        // test tree is the species tree
        // store splits from the test tree
        Split::treeid_t B0;
        Split::treeid_t Blvs;
        _species_forest.storeSplits(B0, Blvs);

        // Only save splits with non-zero edge length
        Split::treeid_t B;
        for (auto & b : B0) {
            if (b.getEdgeLen() > 0.0) {
                B.insert(b);
            }
        }

        // Calculate the contribution of leaf edges to the geodesic
        double leaf_contribution_squared = opCalcLeafContribution(Alvs, Blvs);

        // Find common edges and calculate the contribution of common edge lengths to the geodesic
        vector<Split> common_edges;
        double common_edge_contribution_squared = opFindCommonEdges(A, B, common_edges);

        if (!common_edges.empty()) {
            // Remove splits representing common edge from both A and B
            for (auto & c : common_edges) {
                A.erase(c);
                B.erase(c);
            }
        }
        vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
        in_pairs.emplace_back(A,B);

        unsigned pair_index = 1;
        vector<double> geodesic_distances;
        for (auto & inpair : in_pairs) {

            vector<Split::treeid_pair_t> ABpairs;
            ABpairs.push_back(inpair);

            double L = opCalcGeodesicDist(ABpairs);

            geodesic_distances.push_back(L);
            ++pair_index;
        }

        // Calculate total geodesic distance
        double total_geodesic_distance = 0.0;
        for (double x : geodesic_distances) {
            total_geodesic_distance += pow(x, 2);
        }
        total_geodesic_distance += leaf_contribution_squared;
        total_geodesic_distance += common_edge_contribution_squared;
        total_geodesic_distance = sqrt(total_geodesic_distance);

        return total_geodesic_distance;
    }
        
    inline vector<double> Particle::calcBHVDistanceGeneTrees() {
        vector<double> gene_tree_distances;
        
        for (unsigned l=0; l<G::_nloci; l++) {
            // build the reference tree - rebuild each time to avoid issues with modifying A
            
            Forest reftree;
            
            reftree.buildFromNewick(G::_bhv_reference, true, false);
            
            // Store splits from the reference tree
            Split::treeid_t A0;
            Split::treeid_t Alvs;
            reftree.storeSplits(A0, Alvs);

            // Only save splits with non-zero edge length
            Split::treeid_t A;
            for (auto & a : A0) {
                if (a.getEdgeLen() > G::_small_enough) {
                    A.insert(a);
                }
            }
            
            // test tree is the gene tree
            // store splits from the test tree
            Split::treeid_t B0;
            Split::treeid_t Blvs;
            _gene_forest_ptrs[l]->storeSplits(B0, Blvs);

            // Only save splits with non-zero edge length
            Split::treeid_t B;
            for (auto & b : B0) {
                if (b.getEdgeLen() > 0.0) {
                    B.insert(b);
                }
            }

            // Calculate the contribution of leaf edges to the geodesic
            double leaf_contribution_squared = opCalcLeafContribution(Alvs, Blvs);

            // Find common edges and calculate the contribution of common edge lengths to the geodesic
            vector<Split> common_edges;
            double common_edge_contribution_squared = opFindCommonEdges(A, B, common_edges);

            if (!common_edges.empty()) {
                // Remove splits representing common edge from both A and B
                for (auto & c : common_edges) {
                    A.erase(c);
                    B.erase(c);
                }
            }
            vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
            in_pairs.emplace_back(A,B);

            unsigned pair_index = 1;
            vector<double> geodesic_distances;
            for (auto & inpair : in_pairs) {

                vector<Split::treeid_pair_t> ABpairs;
                ABpairs.push_back(inpair);

                double L = opCalcGeodesicDist(ABpairs);

                geodesic_distances.push_back(L);
                ++pair_index;
            }

            // Calculate total geodesic distance
            double total_geodesic_distance = 0.0;
            for (double x : geodesic_distances) {
                total_geodesic_distance += pow(x, 2);
            }
            total_geodesic_distance += leaf_contribution_squared;
            total_geodesic_distance += common_edge_contribution_squared;
            total_geodesic_distance = sqrt(total_geodesic_distance);

            gene_tree_distances.push_back(total_geodesic_distance);
//            return total_geodesic_distance;
        }
        return gene_tree_distances;
    }
        
    inline double Particle::opCalcLeafContribution(const Split::treeid_t & Alvs, const Split::treeid_t & Blvs) const {
        // Compute leaf contribution
        double leaf_contribution_squared = 0.0;
        for (auto & b : Blvs) {
            auto it = find(Alvs.begin(), Alvs.end(), b);
            assert(it != Alvs.end());
            double leafa = it->getEdgeLen();
            double leafb = b.getEdgeLen();
            leaf_contribution_squared += pow(leafa-leafb, 2);
        }

        return leaf_contribution_squared;
    }
        
    inline double Particle::opFindCommonEdges(const Split::treeid_t & A, const Split::treeid_t & B, vector<Split> & common_edges) const {
        // Find splits in the intersection of A and B
        set_intersection(
            A.begin(), A.end(),
            B.begin(), B.end(),
            back_inserter(common_edges)
        );

        double common_edge_contribution_squared = 0.0;
        for (auto & s : common_edges) {
            auto itA = find(A.begin(), A.end(), s);
            double edgeA = itA->getEdgeLen();
            auto itB = find(B.begin(), B.end(), s);
            double edgeB = itB->getEdgeLen();
            common_edge_contribution_squared += pow(edgeA-edgeB, 2);
        }

        // Count the number of splits in B compatible with each split in A (and vice versa)
        map<const Split *, unsigned> acompatibilities;
        map<const Split *, unsigned> bcompatibilities;
        for (auto & a : A) {
            for (auto & b : B) {
                if (a.compatibleWith(b)) {
                    acompatibilities[&a]++;
                    bcompatibilities[&b]++;
                }
            }
        }

        // Add splits in A that are compatible with all splits in B to common_edges
        for (auto & apair : acompatibilities) {
            if (apair.second == B.size()) {
                // This split is compatible with every split in the other tree, so add it to the vector
                // of common edges if it is not already in the vector of common edges
                if (find(common_edges.begin(), common_edges.end(), *apair.first) == common_edges.end()) {
                    common_edges.push_back(*apair.first);
                    common_edge_contribution_squared += pow(apair.first->getEdgeLen(), 2);
                }
            }
        }

        // Add splits in B that are compatible with all splits in A to common_edges
        for (auto & bpair : bcompatibilities) {
            if (bpair.second == A.size()) {
                // This split is compatible with every split in the other tree, so add it to the vector
                // of common edges if it is not already in the vector of common edges
                if (find(common_edges.begin(), common_edges.end(), *bpair.first) == common_edges.end()) {
                    common_edges.push_back(*bpair.first);
                    common_edge_contribution_squared += pow(bpair.first->getEdgeLen(), 2);
                }
            }
        }

        return common_edge_contribution_squared;
    }
        
    inline double Particle::opCalcGeodesicDist(vector<Split::treeid_pair_t> & ABpairs) const {
        // Assumes a_splits and b_splits have no common edges
        vector<Split::treeid_pair_t> support;
        bool done = false;
        while (!done) {
            unsigned nrefinements = 0;
            for (auto & ABpair : ABpairs) {
                Split::treeid_pair_t AB1;
                Split::treeid_pair_t AB2;
                bool success = opRefineSupport(ABpair, AB1, AB2);
                if (success) {
                    // ABpair was successfully refined, so add AB1 and AB2 to support
                    support.push_back(AB1);
                    support.push_back(AB2);
                    nrefinements++;
                }
                else {
                    // ABpair was not successfully refined, so add ABpair to support
                    support.push_back(ABpair);
                }
            }
            done = (nrefinements == 0);
            if (!done) {
                ABpairs = support;
                support.clear();
            }
        }

        // Calculate geodesic distance
        unsigned ratio_index = 1;
        double geodesic_distance = 0.0;
        for (auto & AB : support) {
            double dropped_length = opCalcTreeIDLength(AB.first);
            double added_length   = opCalcTreeIDLength(AB.second);
            double ratio = dropped_length/added_length;
            geodesic_distance += pow(dropped_length + added_length, 2);

            ++ratio_index;
        }
        geodesic_distance = sqrt(geodesic_distance);
        return geodesic_distance;
    }
        
    inline bool Particle::opRefineSupport(const Split::treeid_pair_t & AB, Split::treeid_pair_t & AB1, Split::treeid_pair_t & AB2) const {
        // Create a vector of incompatibility graph vertices
        vector<op::OPVertex> avect(AB.first.size());
        vector<op::OPVertex> bvect(AB.second.size());

        //                                      23  123  456  123456
        //vector<unsigned> a_splits_in_order = { 6,   7,  56,     63};

        // Calculate weights for the "A" vertices.
        // For example, if A = {a1,a2,a3}, then the weight of a1 is
        // weight a1 = a1^2 / (a1^2 + a2^2 + a3^2)
        unsigned aindex = 0;
        double asum = 0.0;
        for (auto & a : AB.first) {
            asum += pow(a.getEdgeLen(),2);
        }
        for (auto & a : AB.first) {
            avect[aindex]._split = &a;
            avect[aindex]._weight = pow(a.getEdgeLen(),2)/asum;
            aindex++;
        }

        // Calculate weights for the "B" vertices.
        // For example, if B = {b1,b2,b3}, then the weight of b1 is
        // weight b1 = b1^2 / (b1^2 + b2^2 + b3^2)
        unsigned bindex = 0;
        double bsum = 0.0;
        for (auto & b : AB.second) {
            bsum += pow(b.getEdgeLen(),2);
        }
        
        for (auto & b : AB.second) {
            bvect[bindex]._split = &b;
            bvect[bindex]._weight = pow(b.getEdgeLen(),2)/bsum;
            bindex++;
        }

        // Create the incompatibility graph
        // A vertices go on the left, B vertices go on the right, and edges connect an A vertex
        // to a B vertex only if the two vertices are incompatible.
        vector<op::OPEdge> all_edges;
        all_edges.reserve(avect.size() * bvect.size() + avect.size() + bvect.size());
        unsigned nincompatibilities = 0;
        auto asize = static_cast<unsigned>(avect.size());
        auto bsize = static_cast<unsigned>(bvect.size());

        // Create edges from source to the A-vertices
        op::OPVertex source;
        source._name = "source";
        for (unsigned i = 0; i < asize; i++) {
            // Assign a name to avect[i]
            avect[i]._name = str(format("a%d") % i);

            // //temporary!
            // cerr << "avect[" << i << "] split: \n";
            // for (auto v : avect[i]._split->getBits()) {
            //     cerr << " " << v << endl;
            // }

            // Create the forward edge
            all_edges.emplace_back();
            op::OPEdge & source_forward_edge = all_edges.back();
            source_forward_edge._from = &source;
            source_forward_edge._to = &avect[i];
            source_forward_edge._capacity = avect[i]._weight;
            source_forward_edge._flow = 0.0;
            source_forward_edge._reverse_flow = 0.0;
            source_forward_edge._open = true;
            source._edges.push_back(&source_forward_edge);
        }

        // Create edges from A-vertices to B-vertices
        for (unsigned i = 0; i < asize; i++) {
            // Get split associated with avect[i]
            const Split * a = avect[i]._split;
            assert(a);
            for (unsigned j = 0; j < bsize; j++) {
                // Get split associated with bvect[j]
                const Split * b = bvect[j]._split;
                assert(b);
                if (!a->compatibleWith(*b)) {
                    nincompatibilities++;

                    // Create a forward edge in the incompatibility graph
                    all_edges.emplace_back();
                    op::OPEdge & forward_edge = all_edges.back();
                    forward_edge._from = &avect[i];
                    forward_edge._to = &bvect[j];
                    forward_edge._capacity = 1.0;
                    forward_edge._flow = 0.0;
                    forward_edge._reverse_flow = 0.0;
                    forward_edge._open = true;
                    avect[i]._edges.emplace_back(&forward_edge);
                }
            }
        }

        // Create edges from the B-vertices to the sink
        op::OPVertex sink;
        sink._name = "sink";
        for (unsigned j = 0; j < bsize; j++) {
            // Assign a name to bvect[j]
            bvect[j]._name = str(format("b%d") % j);

            // //temporary!
            // cerr << "bvect[" << j << "] split: \n";
            // for (auto v : bvect[j]._split->getBits()) {
            //     cerr << " " << v << endl;
            // }

            // Create the forward edge
            all_edges.emplace_back();
            op::OPEdge & sink_forward_edge = all_edges.back();
            sink_forward_edge._from = &bvect[j];
            sink_forward_edge._to = &sink;
            sink_forward_edge._capacity = bvect[j]._weight;
            sink_forward_edge._flow = 0.0;
            sink_forward_edge._reverse_flow = 0.0;
            sink_forward_edge._open = true;
            bvect[j]._edges.emplace_back(&sink_forward_edge);
        }

        bool success = false;
        if (nincompatibilities < asize*bsize) {
            // At least one independent pair of vertices exists
            // Carry out Edmonds-Karp algorithm to find min-weight vertex cover (identifies max weight independent set)
            // In Owens-Provan terminology,
            //   C1 = A's contribution to vertex cover     (equals AB1.first)
            //   C2 = A's contribution to independent set  (equals AB2.first)
            //   D1 = B's contribution to independent set  (equals AB1.second)
            //   D2 = B's contribution to vertex cover     (equals AB2.second)
            // A1 and A2 must represent a non-trivial partition of A
            // B1 and B2 must represent a non-trivial partition of B
            // C2 and D1 (AB2.first and AB1.second) are compatible sets of splits
            // C1 and D2 (AB1.first and AB2.second) compose the minimum weight vertex cover
            // ||C1||/||D1|| < ||C2||/||D2|| must be true
            // If all of the above conditions hold, then the support can be refined:
            // A = (A1,A2) and B = (B1,B2)
            // where A1 = C1, A2 = C2, B1 = D1, and B2 = D2
            // Length of this segment of the geodesic is
            //   L = sqrt{ (||A1|| + ||B1||)^2 +  (||A2|| + ||B2||)^2 }
            // Orthants crossed:
            //   start:        A1, A2
            //      edges in A1 decreasing, edges in B1 increasing
            //   intermediate: B1, A2
            //      edges in A2 decreasing, edges in B2 increasing
            //   finish:       B1, B2
            // The point at which A1 edges become 0 is
            //   ||A1||/(||A1|| + ||B1||)
            // The point at which A2 edges become 0 is
            //   ||A2||/(||A2|| + ||B2||)
            opEdmondsKarp(source, sink, avect, bvect, AB1.first, AB2.first, AB1.second, AB2.second);

            // Check conditions for successful refinement
            bool A_is_trivial = (AB1.first.empty()) || (AB2.first.empty());
            bool B_is_trivial = (AB1.second.empty())|| (AB2.second.empty());
            double C1len = opCalcTreeIDLength(AB1.first);
            double C2len = opCalcTreeIDLength(AB2.first);
            double D1len = opCalcTreeIDLength(AB1.second);
            double D2len = opCalcTreeIDLength(AB2.second);
            bool P2 = C1len/D1len < C2len/D2len;
            success = !A_is_trivial && !B_is_trivial && P2;
        }
        return success;
    }
        
    inline double Particle::opCalcTreeIDLength(const Split::treeid_t & splits) {
        double length = 0.0;
        for (auto & split : splits) {
            length += pow(split.getEdgeLen(),2);
        }
        return sqrt(length);
    }
        
    inline void Particle::opEdmondsKarp(
                                        op::OPVertex & source,
                                        op::OPVertex & sink,
                                        vector<op::OPVertex> & avect,
                                        vector<op::OPVertex> & bvect,
                                        Split::treeid_t & C1,
                                        Split::treeid_t & C2,
                                        Split::treeid_t & D1,
                                        Split::treeid_t & D2) {

                double cumulative_flow = 0.0;
                bool done_augmenting_path = false;
                while (!done_augmenting_path) {
                    vector<op::OPVertex *> route;

                    // All vertices begin as unvisited, have residual capacity 0, and all "A" to "B" edges begin as non-reversed
                    source._parent_edge = nullptr;
                    for (auto & a : avect) {
                        a._parent_edge = nullptr;
                        a._residual_capacity = 0.0;
                        for (auto & a_edge : a._edges) {
                            a_edge->_edge_is_reversed = false;
                        }
                    }
                    for (auto & b : bvect) {
                        b._parent_edge = nullptr;
                    }
                    sink._parent_edge = nullptr;

                    // Add the source to the route
                    route.push_back(&source);
                    unsigned route_cursor = 0;

                    // Find the next augmenting route
                    bool sink_found = false;
                    while (!sink_found) {
                        op::OPVertex * current = route[route_cursor];
                        for (auto & edge : current->_edges) {
                            // //temporary!
                            // cout << "checking " << edge->_from->_name << " -> " << edge->_to->_name << endl;
                            if (edge->_capacity > 0.0 && edge->_to->_parent_edge == nullptr) {
                                edge->_to->_parent_edge = edge;

                                // //temporary!
                                // cout << "  adding edge to route" << endl;

                                route.push_back(edge->_to);
                                if (edge->_to == &sink) {
                                    sink_found = true;
                                    break;
                                }
                            }
                            else {
                                bool already_visited = edge->_to->_parent_edge != nullptr;
                                bool capacity_zero = edge->_capacity == 0.0;
                                if (already_visited && capacity_zero) {
                                    // //temporary!
                                    // cout << "  rejected edge because capacity is zero and to-vertex already visited" << endl;
                                }
                                else if (already_visited) {
                                    // //temporary!
                                    // cout << "  rejected edge because to-vertex already visited" << endl;
                                }
                                else {
                                    // //temporary!
                                    // cout << "  rejected edge because capacity is zero" << endl;

                                    // If it is a "B" vertex that has zero capacity, then the route cannot
                                    // go from thix "B" vertex to the sink but it may be able to go back to
                                    // an "A" vertex if that "A" vertex was not accessible from the source
                                    // and there is residual flow on the edge
                                    if (edge->_to == &sink) {
                                        op::OPVertex * B_vertex = edge->_from;
                                        for (unsigned j = 0; j < source._edges.size(); j++) {
                                            op::OPEdge * s_to_A_edge = source._edges[j];  // edge from source to an "A" vertex
                                            if (s_to_A_edge->_capacity == 0) {
                                                // See if any of the "A" vertex edges lead to the "B" vertex in question
                                                op::OPVertex * A_vertex = s_to_A_edge->_to;
                                                for (unsigned k = 0; k < A_vertex->_edges.size(); k++) {
                                                    op::OPEdge * A_to_B_edge = A_vertex->_edges[k];  // edge from "A" to a "B" vertex
                                                    if (A_vertex->_parent_edge == nullptr && A_to_B_edge->_to == B_vertex && A_to_B_edge->_reverse_flow > 0.0) {
                                                        // If we are here, we know that:
                                                        // 1. the "A" vertex has not been visited
                                                        // 2. the "A" vertex is joined to the "B" vertex in question
                                                        // 3. the "A" vertex is not accessible from the source, and
                                                        // 4. there is residual flow on the A_to_B_edge

                                                        // //temporary!
                                                        // cout << "  adding REVERSE edge to route (" << B_vertex->_name << " -> " << A_vertex->_name << ")" << endl;

                                                        // //temporary!
                                                        // if (A_to_B_edge->_from->_name == "a2" && A_to_B_edge->_to->_name == "b4") {
                                                        //     cerr << "*** debug breakpoint ***" << endl;
                                                        // }

                                                        A_vertex->_parent_edge = A_to_B_edge;
                                                        A_to_B_edge->_edge_is_reversed = true;
                                                        A_vertex->_residual_capacity = A_to_B_edge->_reverse_flow;
                                                        route.push_back(A_vertex);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        route_cursor++;
                        if (route_cursor >= route.size()) {
                            // Not able to reach sink, so all augmenting routes have been found
                            done_augmenting_path = true;
                            break;
                        }
                    }

                    if (!done_augmenting_path) {
                        // Follow _from pointers from sink back to source to identify
                        // the bottleneck and determine flow through this route
                        double min_capacity = 1.0;
                        op::OPVertex * current = &sink;
                        while (current != &source) {
                            // //temporary!
                            //cerr << "current = " << current->_name << endl;
                            op::OPEdge * edge = current->_parent_edge;
                            if (edge->_edge_is_reversed) {
                                if (edge->_reverse_flow < min_capacity) {
                                    min_capacity = edge->_reverse_flow;
                                }
                                current = edge->_to;
                            }
                            else {
                                if (edge->_capacity < min_capacity) {
                                    min_capacity = edge->_capacity;
                                }
                                current = edge->_from;
                            }
                        }

                        // Adjust capacities and flows along the route
                        cumulative_flow += min_capacity;
                        current = &sink;
                        while (current != &source) {
                            op::OPEdge * edge = current->_parent_edge;
                            if (edge->_edge_is_reversed) {
                                edge->_reverse_flow -= min_capacity;
                                if (fabs(edge->_reverse_flow) < 1e-10)
                                    edge->_reverse_flow = 0.0;
                                edge->_capacity += min_capacity;
                                current = edge->_to;
                            }
                            else {
                                edge->_capacity -= min_capacity;
                                if (fabs(edge->_capacity) < 1e-10)
                                    edge->_capacity = 0.0;
                                edge->_reverse_flow += min_capacity;
                                current = edge->_from;
                            }
                        }
                    }
                }   // while (!done_augmenting_path)

            // Identify C1, C2, D1, and D2
            // C1 and D2 compose the min weight vertex cover
            // C2 and D` compose the independent set
            for (auto & source_edge : source._edges) {
                op::OPVertex * avertex = source_edge->_to;
                if (source_edge->_capacity > 0.0 || avertex->_residual_capacity > 0.0) {
                    // Because this source edge has remaining capacity, its distal vertex is part of the independent set
                    C2.insert(*(avertex->_split));

                    // This A vertex allows access to the B side, so any connected B vertices with
                    // zero capacity are part of the vertex cover
                    for (const auto central_edge : avertex->_edges) {
                        op::OPVertex * bvertex = central_edge->_to;
                        assert(bvertex->_edges.size() == 1);
                        op::OPEdge * sink_edge = bvertex->_edges[0];
                        if (sink_edge->_capacity == 0.0) {
                            D2.insert(*(bvertex->_split));
                        }
                    }
                }
                else {
                    // Because this A-vertex has zero capacity, it is part of the vertex cover
                    C1.insert(*(avertex->_split));

                    // No need to consider connected B vertices because this A-vertex already
                    // covers all connected edges
                }
            }

            // D1 includes every split not in D2
            D1.clear();
            for (auto & b : bvect) {
                if (D2.count(*(b._split)) == 0) {
                    D1.insert(*(b._split));
                }
            }
        }

    inline void Particle::operator=(const Particle & other) {
        if (!G::_in_second_level) {
            _log_likelihood = other._log_likelihood;
            _t_by_gene = other._t_by_gene;
            _t = other._t;
            _next_species_number_by_gene = other._next_species_number_by_gene;
            _gene_order = other._gene_order;
            _theta_map = other._theta_map;
            _prev_gene_increment = other._prev_gene_increment;
            _prev_total_to_add = other._prev_total_to_add;
            _prev_next_species_number_by_gene = other._prev_next_species_number_by_gene;
            _prev_t_by_gene = other._prev_t_by_gene; // TODO: not sure if these need to be copied
            _prev_increment_prior = other._prev_increment_prior;
            _prev_species_assignments_before_coalescence = other._prev_species_assignments_before_coalescence;
        }

        _prev_log_likelihoods = other._prev_log_likelihoods;
        
        if (!G::_mcmc) {
            // Ensure that _gene_forest_extensions is allocated
            if (_gene_forest_extensions.size() == 0) {
                _gene_forest_extensions.resize(G::_nloci);
            }
            else {
                assert(_gene_forest_extensions.size() == G::_nloci);
                
                // Undock all gene forest extensions
                for_each(_gene_forest_extensions.begin(), _gene_forest_extensions.end(),
                [](ForestExtension & f){f.undock();});
            }
        }
        else {
            _gene_forest_extensions = other._gene_forest_extensions; // for mcmc proposals, forest extensions have not been docked yet
        }
        
        _gene_forest_ptrs = other._gene_forest_ptrs;
        _ensemble_coalinfo = other._ensemble_coalinfo;
        
        _log_weight     = other._log_weight;
        _species_forest = other._species_forest;
        _psuffix = other._psuffix;
        _theta_mean = other._theta_mean;
        _theta_vector = other._theta_vector;

        // the following data members apply only when simulating and do not need to be copied because simulating data only deals with one particle at a time
//            _num_deep_coalescences = other._num_deep_coalescences;
//            _max_deep_coal = other._max_deep_coal;

        if (G::_in_second_level) {
            _log_coalescent_likelihood = other._log_coalescent_likelihood;
            _group_number = other._group_number;
        }
    };
}

