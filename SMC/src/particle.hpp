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

#if defined (FOSSILS)
    #include "fossil.hpp"
    #include "taxset.hpp"
#endif

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
                                                    _nsubsets = d->getNumSubsets();
                                                    _data = d;
                                                    int index = 0;
                                                    _forests.resize(_nsubsets+1);
                                                    for (auto &_forest:_forests) {
                                                        if (index > 0) {
                                                            _forest.setData(d, index, taxon_map, partials);
                                                        }
                                                        index++;
                                                    }
                                                }
        void                                    setSimData(Data::SharedPtr d, map<string, string> &taxon_map, unsigned nsubsets, unsigned ntaxa) {
                                                    _nsubsets = nsubsets;
                                                    int index = 0;
                                                    _forests.resize(_nsubsets+1);
                                                    for (auto &_forest:_forests) {
                                                        _forest.setStartMode("sim");
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
            
        string                                  saveGeneNewick(unsigned i) {
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
        int                                             selectEventLinearScale(vector<double> weight_vec);
        double                                          getTopologyPrior(unsigned i);
        vector<pair<double, double>>                    getIncrementPriors(unsigned i);
        vector<pair<double, double>>                    getSpeciesTreeIncrementPriors();
        double                                          getCoalescentLikelihood(unsigned g);
        bool                                            speciesJoinProposed() {return _species_join_proposed;}
        void                                            clear();
        vector<double>                                  chooseIncrements(vector<double> event_choice_rates);
        void                                            speciesOnlyProposal();
        void                                            speciesOnlyProposalIntegratingOutTheta();
        void                                            calculateIncrementPriors(double increment, string species_name, unsigned forest_number, bool speciation, bool first_step);
        void                                            changeTheta(unsigned i);
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
        double                                          getLambda(){return _forests[0]._lambda;}
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
        void                                            setNextSpeciesNumber() {_next_species_number = Forest::_nspecies;}
        unsigned                                        showPrevForestNumber(){return _prev_forest_number;}
        string                                          getTranslateBlock();
        void                                            buildEntireSpeciesTree();
        void                                            rebuildSpeciesTree();
        void                                            setGeneOrder(vector<unsigned> gene_order) {_gene_order = gene_order;}
        void                                            trimSpeciesTree();
        void                                            setFixTheta(bool fix) {_fix_theta = fix;}
        void                                            setRelativeRatesByGene(vector<double> rel_rates);
#if defined (FASTER_UPGMA_TREE)
        void                                            calcStartingUPGMAMatrix();
        vector<vector<double>>                          getStartingUPGMAMatrix();
        void                                            setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene);
        vector<map<Node*,  unsigned>>                   getStartingRowCount();
        void                                            setStartingRowCount(vector<map<Node*,  unsigned>> starting_row_count_by_gene);
        void                                            calcStartingRowCount();
#endif
        void                                            createSpeciesIndices();
        void                                            drawParticleLambda();
        double                                          getParticleLambda(){return _forests[0]._lambda;}
        void                                            setParticleLambda(double lambda){_forests[0]._lambda = lambda;}
        void                                            setParticleExtinctionRate(double extinction_rate){_forests[0]._extinction_rate = extinction_rate;}
        void                                            setNTaxaPerSpecies(vector<unsigned> ntaxa_per_species);
//
        static double                                   _lambda_prior_mean;
//        vector<string>                                  _fossil;
#if defined(FOSSILS)
        static vector<Fossil>           _fossils;
        static vector<TaxSet>           _taxsets;
#endif

    private:

        static unsigned                         _nsubsets;
        vector<Forest>                          _forests;
        double                                  _log_weight;
        double                                  _log_species_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        int                                     _generation = 0;
        unsigned                                _prev_forest_number;
        bool                                    _species_join_proposed;
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
        unsigned                                _species_branches;
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
        _nsubsets       = 0;
        _generation     = 0;
        _prev_forest_number = 0;
        _species_join_proposed = false;
        _log_coalescent_likelihood = 0.0;
        _num_deep_coalescences = 0.0;
        _max_deep_coal = 0.0;
        _species_tree_height = 0.0;
        _t.clear();
        _psuffix = 0;
        _next_species_number = Forest::_nspecies;
        _species_order.clear();
        _starting_log_likelihoods.clear();
        _t_by_gene.clear();
        _next_species_number_by_gene.clear();
        _gene_order.clear();
        _fix_theta = false;
        _relative_rates_by_gene.clear();
        _species_branches = 0;
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
        if (_generation == 0 && !Forest::_run_on_empty) {
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
            if (!Forest::_run_on_empty) {
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
        total_prior += log(Forest::_theta_prior_mean) - (_forests[1]._theta * Forest::_theta_prior_mean);
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
        if (_generation == 0 && !Forest::_run_on_empty) {
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
//        else {
        else if (_generation % _nsubsets == 0 && Forest::_start_mode == "smc") { // after every locus has been filtered once, trim back the species tree as far as possible & rebuild it
            // don't rebuild the species tree at every step when simulating data
            // TODO: do this at every step, not just after filtering all loci once? - did not work as well on simulation grid
            trimSpeciesTree();
            if (_forests[0]._lineages.size() > 1) {
                rebuildSpeciesTree();
            }
        }
        
        _species_join_proposed = false;
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
                        
                        unsigned index = selectEventLinearScale(event_choice_rates); // TODO: do the rates have to be in order?
                        string species_name = rates_by_species[index].second;
                        _forests[next_gene].allowCoalescence(species_name, gene_increment, _lot);
                                   
                        if (Forest::_start_mode == "smc") {
# if defined (BUILD_UPGMA_TREE)
# if defined (BUILD_UPGMA_TREE_CONSTRAINED)
                        _forests[next_gene].buildRestOfTree(_lot, _t);
#else
                        _forests[next_gene].buildRestOfTree(_lot);
#endif
#endif
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
                            
                            if (Forest::_start_mode == "sim") {
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
                 assert (_forests[1]._theta > 0.0);
                 for (int i=1; i<_forests.size(); i++) {
                     _forests[i]._theta_mean = _forests[1]._theta;
                 }
             }
             assert (_forests[1]._theta_mean > 0.0);
             
             done = true;
        
#if defined (INV_GAMMA_PRIOR_TWO)
        // TODO: use 2.0, not 2.01, for the inv gamma, then don't include a correction
            // include inverse gamma prior correction for every species population for every locus at every step
            // TODO: can make this faster by calculating these as the thetas are drawn, but for now, calculate them all each time
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
        // modifier only happens on first round // TODO: unsure - does gene order matter?
        if (_generation == 0 && _forests[1]._theta_prior_mean > 0.0 && _forests[1]._theta_proposal_mean > 0.0) {
            if (_forests[1]._theta_prior_mean != _forests[1]._theta_proposal_mean) {
                // else, log weight modifier is 0
                double prior_rate = 1.0/_forests[1]._theta_prior_mean;
                double proposal_rate = 1.0/_forests[1]._theta_proposal_mean;
                double log_weight_modifier = log(prior_rate) - log(proposal_rate) - (prior_rate - proposal_rate)*_forests[1]._theta_mean;

                _log_weight += log_weight_modifier;
            }
        }
#endif
        
        _generation++;
        
        if (Forest::_start_mode == "smc") {
            _log_weight = _forests[next_gene]._log_weight + inv_gamma_modifier;
        }
        else {
            _log_weight = 0.0;
        }
        
            }

    vector<double> Particle::chooseIncrements(vector<double> event_choice_rates) {
        vector<double> increments;
        increments.resize(event_choice_rates.size());
        
        for (int p=0; p<event_choice_rates.size(); p++) {
            increments[p] = -log(1.0 - _lot->uniform())/event_choice_rates[p];
        }
        return increments;
    }

    inline void Particle::speciesOnlyProposalIntegratingOutTheta() {
//        for (unsigned i=0; i<10; i++) {
//            double u = _lot->uniform();
//            cout << u << endl;
//        }
        
        if (_generation == 0) {
            _species_branches = Forest::_nspecies;
            if (_lambda_prior_mean > 0.0) {
                _forests[0]._lambda = _lot->gamma(1, _lambda_prior_mean);
            }
            for (int i=1; i<_forests.size(); i++) {
                _forests[i].refreshPreorder();
                _forests[i].calcMinDepth();
                _forests[i]._nincrements = 0;

                if (i > 1) {
                    _forests[i]._theta_mean = _forests[1]._theta_mean;
                }
//                else if (i == 1) {
//                    _forests[1]._theta_map.clear();
//                    _forests[1].resetThetaMap(_lot); // reset tip thetas and ancestral pop theta
//                }
                // TODO: careful - trying so random seed order matches non jones coalescent like version
            }
        }
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        double prev_log_coalescent_likelihood = _log_coalescent_likelihood;
        
//        _forests[0].showForest();
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
                        
//                    _forests[i].showForest();
                    if (_forests[0]._lineages.size() > 1 && species1 != "null") {
                        // if using Jones formula, species partition update will happen in coalescent likelihood calculation
                        _forests[i].resetDepthVector(species_joined);
                    }

//                    _forests[i].showForest();
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

        if (!Forest::_run_on_empty) {
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
        if (!Forest::_run_on_empty) {
            assert (max_depth > 0.0);
            double nlineages = _forests[0]._lineages.size();
            if (nlineages == 1) {
                nlineages = 2; // for last step, constraint was before final two species were joined
            }
//            constrained_factor = log(1 - exp(-1*nlineages*Forest::_lambda*max_depth));
            constrained_factor = log(1 - exp(-1*nlineages*_forests[0]._lambda*max_depth));
        }
#endif
        
#endif
            _log_species_weight = _log_coalescent_likelihood - prev_log_coalescent_likelihood + constrained_factor;
#if !defined (UNCONSTRAINED_PROPOSAL)
            double test = 1/_log_species_weight;
            assert(test != -0); // assert coalescent likelihood is not -inf
#endif
        
        if (Forest::_run_on_empty) {
            assert (_log_coalescent_likelihood == 0.0);
            assert (_log_species_weight == 0.0);
        }
        _generation++;
    }

    inline void Particle::speciesOnlyProposal() {
#if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
        speciesOnlyProposalIntegratingOutTheta();
#else
        
        if (_generation == 0) {
            if (_lambda_prior_mean > 0.0) {
                _forests[0]._lambda = _lot->gamma(1, _lambda_prior_mean);
            }
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

#if defined (FASTER_UPGMA_TREE)
    inline void Particle::calcStartingUPGMAMatrix() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].buildStartingUPGMAMatrix(); // TODO: can do this once and copy to all particles
        }
    }
#endif

#if defined (FASTER_UPGMA_TREE)
    inline vector<vector<double>> Particle::getStartingUPGMAMatrix() {
        vector<vector<double>> starting_upgma_matrices_by_gene;
        for (unsigned i=1; i<_forests.size(); i++) {
            starting_upgma_matrices_by_gene.push_back(_forests[i]._starting_dij);
        }
        return starting_upgma_matrices_by_gene;
    }
#endif

#if defined (FASTER_UPGMA_TREE)
    inline void Particle::setStartingUPGMAMatrix(vector<vector<double>> starting_upgma_matrices_by_gene) {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._starting_dij = starting_upgma_matrices_by_gene[i-1];
        }
    }
#endif

#if defined (FASTER_UPGMA_TREE)
    inline vector<map<Node*,  unsigned>> Particle::getStartingRowCount() {
        vector<map<Node*,  unsigned>> starting_row_count_by_gene;
        for (unsigned i=1; i<_forests.size(); i++) {
            starting_row_count_by_gene.push_back(_forests[i]._starting_row);
        }
        return starting_row_count_by_gene;
    }
#endif

#if defined (FASTER_UPGMA_TREE)
    inline void Particle::setStartingRowCount(vector<map<Node*,  unsigned>> starting_row_count_by_gene) {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i]._starting_row = starting_row_count_by_gene[i-1];
        }
    }
#endif

#if defined (FASTER_UPGMA_TREE)
    inline void Particle::calcStartingRowCount() {
        for (unsigned i=1; i<_forests.size(); i++) {
            _forests[i].buildStartingRow();
        }
    }
#endif

    inline double Particle::calcInitialCoalescentLikelihood() {
#if defined GRAHAM_JONES_COALESCENT_LIKELIHOOD
        unsigned nbranches = Forest::_nspecies*2 - 1;
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
            for (int i=0; i<Forest::_nspecies-1; i++) {
                string name = boost::str(boost::format("node-%d")%number);
                number++;
                species_names.push_back(name);
            }
            
            assert (species_names.size() == 2*Forest::_nspecies - 1);
            
            for (auto &name:species_names) {
                assert (Forest::_theta > 0.0);
                theta_map[name] = Forest::_theta;         // create a theta map with all the same theta for simulations, set theta_mean to theta
            }
            
            for (int i=1; i<_forests.size(); i++) {
                _forests[i]._theta_map = theta_map;
                _forests[i]._theta_mean = Forest::_theta;
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
        double theta_proposal_mean = _forests[1]._theta_proposal_mean;
        double theta_prior_mean = _forests[1]._theta_prior_mean;
        map<string, double> theta_map = _forests[1]._theta_map;
        map<string, unsigned> species_indices = _forests[1]._species_indices;
        if (_forests.size() > 2) {
            for (int i=2; i<_forests.size(); i++) {
                _forests[i]._theta_map = theta_map;
                _forests[i]._species_indices = species_indices;
                _forests[i]._theta_mean = theta_mean;
                _forests[i]._theta_proposal_mean = theta_proposal_mean;
                _forests[i]._theta_prior_mean = theta_prior_mean;
            }
        }
    }

    inline void Particle::drawTheta() {
        // set seed first
//        assert (_psuffix > 0);
//        setSeed(rng.randint(1,9999) + _psuffix);
            
        _forests[1].createThetaMap(_lot); // create map for one forest, then copy it to all forests
        double theta_mean = _forests[1]._theta_mean;
        double theta_proposal_mean = _forests[1]._theta_proposal_mean;
        double theta_prior_mean = _forests[1]._theta_prior_mean;
        map<string, double> theta_map = _forests[1]._theta_map;
        map<string, unsigned> species_indices = _forests[1]._species_indices;
        if (_forests.size() > 2) {
            for (int i=2; i<_forests.size(); i++) {
                _forests[i]._theta_map = theta_map;
                _forests[i]._species_indices = species_indices;
                _forests[i]._theta_mean = theta_mean;
                _forests[i]._theta_proposal_mean = theta_proposal_mean;
                _forests[i]._theta_prior_mean = theta_prior_mean;
            }
        }
    }

    inline void Particle::changeTheta(unsigned i) {
        // for snake data set?
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
    }

    inline void Particle::calculateIncrementPriors(double increment, string species_name, unsigned forest_number, bool speciation, bool first_step) {
        // need to calculate coalescent likelihood before joining anything or updating species partition
        for (int f=0; f<_forests.size(); f++) {
            bool new_increment = false;
            bool coalescence = false;
            bool gene_tree = false;
            
            if (first_step) {
                new_increment = true;
            }
//            if ((f == _prev_forest_number || _generation == 0) && !speciation) { // if previous join was a species join, new increment is false
            if (f == _prev_forest_number) {
                // add to existing increment + prior
                new_increment = true;
            }
            if (f == forest_number) {
                coalescence = true;
            }
            if (f > 0) {
                gene_tree = true;
            }
            _forests[f].calcIncrementPrior(increment, species_name, new_increment, coalescence, gene_tree);
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
        return _log_coalescent_likelihood; // can't get coalescent likelihood separately for each gene tree
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
            if (count < Forest::_nspecies + 1) {
                string name = nd._name;
                block += to_string(count) + " ";
                block += name;
                if (count != Forest::_nspecies) {
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
            _forests[i]._theta_mean = Forest::_theta; // for now, set theta mean equal to whatever theta user specifies
        }
    }

    inline void Particle::resetSpecies() {
        _data = nullptr;
        _forests[0].clear();
        _generation = 0;
        setLogWeight(0.0);
        _species_tree_height = 0.0;
        _log_coalescent_likelihood = 0.0;
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._log_coalescent_likelihood = 0.0;
            _forests[i]._data = nullptr;
            _forests[i]._log_weight = 0.0;
//            _forests[i]._species_indices.clear(); // TODO: need to save this for use in jones coalescent likelihood calculation
        }
        _gene_order.clear();
        _next_species_number_by_gene.clear();
        _starting_log_likelihoods.clear();
        _t.clear();
        _t_by_gene.clear();
        _forests[0]._last_edge_length = 0.0;
        _forests[0]._increments_and_priors.clear();
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
                if (_forests[0]._lineages.size() < Forest::_nspecies) {
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

    inline void Particle::drawParticleLambda() {
        if (_lambda_prior_mean > 0.0) {
            // otherwise, lambda is fixed and not estimated
            _forests[0]._lambda = _lot->gamma(1, _lambda_prior_mean);
        }
    }

    inline void Particle::rebuildSpeciesTree() {
        bool trim_to_previous_join = false;
        
        if (trim_to_previous_join) {
            map<string, double> theta_map = _forests[1]._theta_map;
            
            // TODO: testing trimming species tree all the way back to the previous join
            
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
                if (edge_increment < min_branch_length) { // TODO: there must be a better way to do this
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
                        _forests[1].updateThetaMap(_lot, t.first); // TODO: make sure this works regardless of gene order
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
                    _forests[1].updateThetaMap(_lot, t.first); // TODO: make sure this works regardless of gene order
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

        for (unsigned i=0; i < _forests[0]._nspecies-1; i++) {
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

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_species_weight = other._log_species_weight;
        _log_likelihood = other._log_likelihood;
        _forests         = other._forests;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _prev_forest_number = other._prev_forest_number;
        _species_join_proposed = other._species_join_proposed;
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
        _species_branches = other._species_branches;
        _lambda_prior_mean = other._lambda_prior_mean;
    };
}

