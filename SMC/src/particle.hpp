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
        int                                             selectEvent(vector<double> weight_vec);
        double                                          getTopologyPrior(unsigned i);
        vector<pair<double, double>>                    getIncrementPriors(unsigned i);
        vector<pair<double, double>>                    getSpeciesTreeIncrementPriors();
        double                                          getCoalescentLikelihood(unsigned g);
        bool                                            speciesJoinProposed();
        void                                            clear();
        vector<double>                                  chooseIncrements(vector<double> event_choice_rates);
        void                                            speciesProposal();
        void                                            speciesOnlyProposal();
        void                                            speciesOnlyProposalIntegratingOutTheta();
        void                                            geneProposal(vector<unsigned> event_choice_index, unsigned forest_number, vector<string> event_choice_name, double increment, string species_name);
        void                                            calculateIncrementPriors(double increment, string species_name, unsigned forest_number, bool speciation, bool first_step);
        void                                            changeTheta(unsigned i);
        void                                            drawTheta();
        void                                            clearPartials();
        Lot::SharedPtr getLot() const {return _lot;}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}
        double                                          getRunningSumChoices(vector<double> choices);
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
        void                                            resetSpecies();
        void                                            setForest(Forest f, unsigned forest_number);
        Forest                                          getForest(unsigned i) {return _forests[i];} // TODO: should return a pointer?
        void                                            setNewTheta();
        vector<double>                                  getThetaVector();
        double                                          getPopMean();
        pair<string, string>                            getSpeciesJoined(){return make_pair(_forests[0]._species_joined.first->_name, _forests[0]._species_joined.second->_name);}
        void                                            setPsuffix(unsigned psuffix) {_psuffix = psuffix;}
        double                                          calcInitialCoalescentLikelihood();
        void                                            processGeneNewicks(vector<string> newicks);
    
//        static bool                                     _run_on_empty;

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
        bool                                    _deep_coal;
        double                                  _species_tree_height;
        unsigned                                _psuffix;
        vector<pair<tuple<string, string, string>, double>> _t;
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
        _species_tree_height = 0.0;
        _t.clear();
        _psuffix = 0;
        _deep_coal = false;
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
        }
        _log_likelihood = total_log_likelihood;
        _log_weight = total_log_likelihood;
    }

    inline double Particle::getRunningSumChoices(vector<double> choices) {
        double running_sum = 0.0;
        double log_weight_choices_sum = 0.0;
        double log_max_weight = *max_element(choices.begin(), choices.end());
        for (auto & i:choices) {
            running_sum += exp(i - log_max_weight);
        }
        log_weight_choices_sum = log(running_sum) + log_max_weight;
        return log_weight_choices_sum;
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
        _species_join_proposed = false;
        bool done = false;
                
        while (!done) {
    
            bool speciation = false;
            
        vector<double> forest_rates; // this vector contains total rate of species tree, gene 1, etc.
        vector<vector<double>> gene_forest_rates; // this vector contains rates by species for each gene forest
        gene_forest_rates.resize(_forests.size()-1);
        vector<unsigned> event_choice_index;
        vector<string> event_choice_name;
            
        for (int i=0; i<_forests.size(); i++) {
            if (i > 0) {
                
#if defined (SNAKE)
                changeTheta(i);
#endif
                vector<pair<double, string>> rates_by_species = _forests[i].calcForestRate(_lot);
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
                    forest_rates.push_back(Forest::_lambda * _forests[0]._lineages.size());
                    event_choice_index.push_back(0);
                    event_choice_name.push_back("species");
                }
                else {
                    _forests[0]._done = true;
                }
            }
        }
            
            vector<double> event_choice_rates;
            if (_forests[0]._lineages.size() > 1) {
                event_choice_rates.push_back(forest_rates[0]); // push back species tree rate
            }
            for (int i=0; i<gene_forest_rates.size(); i++) {
                for (auto &r:gene_forest_rates[i]) {
                    assert (r > 0.0);
                    event_choice_rates.push_back(r);
                }
            }
            
            bool no_speciation = false;
            double speciation_time = -1;
            unsigned index = 0;
            double increment = 0.0;
            
#if defined (USE_TOTAL_RATE)
            unsigned forest_number = 0;
            string species_name = "species";
            double total_rate = 0.0;
            
            for (auto &r:event_choice_rates) {
                total_rate += r;
            }
            
            if (event_choice_name[0] == "species") {
                total_rate -= event_choice_rates[0]; // remove speciation rate from total rate used for choosing increment
                speciation_time = -log(1.0 - _lot->uniform()) / event_choice_rates[0];
            }
            
            if (total_rate > 0.0) {
                double gene_increment = -log(1.0 - _lot->uniform())/total_rate;
                
                if (gene_increment > speciation_time && speciation_time != -1) {
                    _deep_coal = true;
                    // num deep coalescences += (n-1), where n is number of lineages in each affected species lineage for all genes
                    // need to know which species joined to calculate this
                }
                
                if (gene_increment < speciation_time || speciation_time == -1) {
                    // choose a coalescent event if coalescent increment < speciation increment or if species tree is finished
                    increment = gene_increment;
                    if (event_choice_name[0] == "species") { // remove the speciation event from consideration
                        event_choice_index.erase(event_choice_index.begin() + 0);
                        event_choice_name.erase(event_choice_name.begin() + 0);
                        event_choice_rates.erase(event_choice_rates.begin() + 0);
                    }
                    
                    for (auto &p:event_choice_rates) {
                         p = log(p/total_rate);
                     }
                    
                    index = selectEvent(event_choice_rates);
                                
                    forest_number = event_choice_index[index];
                    species_name = event_choice_name[index];
                    assert (species_name != "species");
                    assert (forest_number != 0);
                    
                }
                else {
                    forest_number = 0;
                    increment = speciation_time;
                    assert (increment != -1.0);
                    no_speciation = false;
                }
            }
            else {
                forest_number = 0;
                increment = speciation_time;
                assert (increment != -1.0);
                no_speciation = false;
            }
            
#else
            vector<double> increments = chooseIncrements(event_choice_rates);
        
//            double speciation_time = -1;
            if (!_forests[0]._done) {
                speciation_time = increments[0];
            }
            
            if (speciation_time > -1) { // otherwise, species tree is done, and there is no constraint on gene tree increments
                for (int i = (int) increments.size()-1; i>0; i--) {
                    if (increments [i] > speciation_time) {
                        event_choice_index.erase(event_choice_index.begin() + i);
                        event_choice_name.erase(event_choice_name.begin() + i);
                        event_choice_rates.erase(event_choice_rates.begin() + i);
                        increments.erase(increments.begin() + i);
                    }
                }
            }
        
            // if a gene forest coalescence is possible, do not pick a speciation event
//            bool no_speciation = false;
            if (event_choice_name[0] == "species" && event_choice_name.size() > 1) {
//                 erase speciation event possibility
                event_choice_index.erase(event_choice_index.begin() + 0);
                event_choice_name.erase(event_choice_name.begin() + 0);
                event_choice_rates.erase(event_choice_rates.begin() + 0);
                increments.erase(increments.begin() + 0);
                no_speciation = true;
            }
            
            double total_rate = 0.0; // normalize rates before selecting an event
            for (auto &r:event_choice_rates) {
                assert (r > 0.0);
                total_rate += r;
            }
        
            // choose the minimum coalescence time
            double min_coalescence_time = 0.0;
            index = 0;
            unsigned forest_number = 0;
            
            if (event_choice_name.size() == 1 && event_choice_name[0] == "species") {
                // choose the speciation event
                increment = speciation_time;
                assert (speciation_time == increments[0]);
            }
            else {
                // choose the minimum event
                min_coalescence_time = *min_element(std::begin(increments), std::end(increments));
                increment = min_coalescence_time;
                bool entered = false;
                for (int i=0; i<increments.size(); i++) {
                    if (increments[i] == min_coalescence_time) {
                        index = i;
                        entered = true;
                        break;
                    }
                }
                assert (entered);
                forest_number = event_choice_index[index];
                if (no_speciation) {
                    assert (forest_number != 0);
                }
            }
        
            string species_name = event_choice_name[index];
#endif
                
            // add increment to all nodes in all forests
            for (int i=0; i<_forests.size(); i++) {
                if (_forests[i]._lineages.size() > 1) {
                    _forests[i].addIncrement(increment); // if forest is finished, don't add another increment
                }
                else {
                    _forests[i]._done = true;
                }
            }
            
            if (forest_number == 0) {
                speciation = true;
            }
            
            bool first_step = true;
            if (_forests[0]._lineages.size() != Forest::_nspecies) {
                first_step = false;
            }
            
            if (first_step) {
                for (int i=1; i<_forests.size(); i++) {
                    if (_forests[i]._lineages.size() != Forest::_ntaxa) {
                        first_step = false;
                    }
                }
            }
            calculateIncrementPriors(increment, species_name, forest_number, speciation, first_step);
            
            if (species_name == "species") {
                unsigned n = (unsigned) _forests[0]._lineages.size();
                assert (n > 1);
                assert (index == 0);
                assert (forest_number == 0);
                
                speciesProposal();
                _species_join_proposed = true;
                assert (increment > 0.0);
            }
        
            else {
                assert (increment > 0.0);
                double log_speciation_term = 0.0;
                geneProposal(event_choice_index, forest_number, event_choice_name, increment, species_name);
                double log_likelihood_term = _forests[forest_number]._log_weight;

                _log_weight = log_speciation_term + log_likelihood_term;
                
//                assert (_forests[1]._theta_prior_mean > 0.0);
//                assert (_forests[1]._theta_proposal_mean > 0.0);
                assert (_forests[1]._theta_mean > 0.0);
                
//                assert (_forests[1]._theta_mean == _forests[10]._theta_mean);
                
                // modifier only happens on first round
                if (_generation == 0 && _forests[1]._theta_prior_mean > 0.0 && _forests[1]._theta_proposal_mean > 0.0) {
                    double prior_rate = 1.0/_forests[1]._theta_prior_mean;
                    double proposal_rate = 1.0/_forests[1]._theta_proposal_mean;
                    double log_weight_modifier = log(prior_rate) - log(proposal_rate) - (prior_rate - proposal_rate)*_forests[1]._theta_mean;
                    
                    _log_weight += log_weight_modifier;
                }
       
                
                done = true;
            }
                  
            if (Forest::_run_on_empty) {
                _log_weight = 0.0;
            }
            
            
            _prev_forest_number = forest_number;
            
            }
        _generation++;
    }

    vector<double> Particle::chooseIncrements(vector<double> event_choice_rates) {
        vector<double> increments;
        increments.resize(event_choice_rates.size());
        
        for (int p=0; p<event_choice_rates.size(); p++) {
            increments[p] = -log(1.0 - _lot->uniform())/event_choice_rates[p];
        }
        return increments;
    }

    inline void Particle::speciesProposal() {
        // species tree proposal, need to update species partition in all gene forests
        assert (_lot != nullptr);
        tuple <string, string, string> species_joined = _forests[0].speciesTreeProposal(_lot);
        for (int i=1; i<_forests.size(); i++) {
            if (_deep_coal) {
                // number of deep coalescences is (n-1) where n is the total number of lineages in the joining species
                string new_spp = get<2>(species_joined);
                _num_deep_coalescences += _forests[i].getDeepCoal(species_joined);
            }
            _forests[i].updateSpeciesPartition(species_joined);
        }
        _deep_coal = false;
    }

    inline void Particle::speciesOnlyProposalIntegratingOutTheta() {
        if (_generation == 0) {
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
        
            if (_forests[0]._last_edge_length > 0.0) {
            // choose species to join if past the first species generation for each forest vector
                species_joined = _forests[0].speciesTreeProposal(_lot);
            }
                vector<double> max_depth_vector;
                double max_depth = 0.0;

//#if !defined (UNCONSTRAINED_PROPOSAL)
                for (int i=1; i<_forests.size(); i++) {
                    string species1 = get<0>(species_joined);
                    string species2 = get<1>(species_joined);
                    
                    if (species1 != "null") {
                        // if using Jones formula, species partition update will happen in coalescent likelihood calculation
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
//#endif
                
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
        
            assert (_log_coalescent_likelihood == 0.0);

        if (!Forest::_run_on_empty) {
            unsigned nbranches = Forest::_nspecies*2 - 1;
            _log_coalescent_likelihood = 2 * nbranches * log(_forests[1]._theta_mean) - nbranches * boost::math::lgamma(2);
            
            vector<double> gamma_jb;
            vector<unsigned> q_jb;
            
            double neg_inf = -1*numeric_limits<double>::infinity();
            for (int i = 1; i<_forests.size(); i++) {
                if (_log_coalescent_likelihood != neg_inf) {
                    pair<vector<double>, vector<unsigned>> params;
                    if (_forests[0]._lineages.size() > 1) {
                        params = _forests[i].calcCoalescentLikelihoodIntegratingOutTheta(_forests[0]._species_build);
                    }
                    else {
                        params = _forests[i].calcCoalescentLikelihoodIntegratingOutThetaLastStep(_forests[0]._species_build); // not using this right now
                    }
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

        double constrained_factor = 0.0;
#if !defined (UNCONSTRAINED_PROPOSAL)
        if (!Forest::_run_on_empty) {
            assert (max_depth > 0.0);
            double nlineages = _forests[0]._lineages.size();
            if (nlineages == 1) {
                nlineages = 2; // for last step, constraint was before final two species were joined
            }
            constrained_factor = log(1 - exp(-1*nlineages*Forest::_lambda*max_depth));
        }
        
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
        
#if defined (UNCONSTRAINED_PROPOSAL)
        if (_forests[0]._last_edge_length > max_depth) {
            _log_coalescent_likelihood = -1*numeric_limits<double>::infinity();
        }
        else {
            for (int i = 1; i<_forests.size(); i++) {
                _forests[i].calcCoalescentLikelihood(_forests[0]._last_edge_length, species_joined, _species_tree_height);
                _log_coalescent_likelihood += _forests[i]._log_coalescent_likelihood + _forests[i]._panmictic_coalescent_likelihood;
            }
        }
#endif
        
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

    inline void Particle::geneProposal(vector<unsigned> event_choice_index, unsigned forest_number, vector<string> event_choice_name, double increment, string species_name) {
        _forests[forest_number].allowCoalescence(species_name, increment, _lot);
    }

    inline void Particle::setNewTheta() {
        // map should be 2*nspecies - 1 size
        // create a theta map with all the same theta for simulations, set theta_mean to theta
        unsigned number = 0;
        vector<string> species_names;
        map<string, double> theta_map;
        
        for (auto &s:_forests[1]._species_partition) {
            species_names.push_back(s.first);
            number++;
        }
        for (int i=0; i<Forest::Forest::_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            number++;
            species_names.push_back(name);
        }
        
        assert (species_names.size() == 2*Forest::_nspecies - 1);
        
        for (auto &name:species_names) {
            theta_map[name] = Forest::_theta;
        }
        
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._theta_map = theta_map;
            _forests[i]._theta_mean = Forest::_theta;
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

#if defined (SIM_TEST)
            if (_species_join_proposed) {
                new_increment = false;
            }
#endif
#if defined (SIM_TES3T)
            if (_species_join_proposed) {
                new_increment = false;
            }
            if (f == 0) {
                new_increment = true;
            }
#endif
//            new_increment = true;
            _forests[f].calcIncrementPrior(increment, species_name, new_increment, coalescence, gene_tree);
#if defined (SIM_TEST3)
            _species_join_proposed = false;
#endif
        }
    }

    inline Particle::Particle(const Particle & other) {
        _lot.reset(new Lot());
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

    inline int Particle::selectEvent(vector<double> weight_vec) {
        // choose a random number [0,1]
        double u =  _lot->uniform();
        assert (u > 0.0);
        assert (u < 1.0);
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

    inline bool Particle::speciesJoinProposed() {
        if (_species_join_proposed) {
            return true;
        }
        else {
            return false;
        }
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

    inline void Particle::processGeneNewicks(vector<string> newicks) {
        for (int i=1; i<_forests.size(); i++) {
//            _forests[i]._nodes.clear();
            _forests[i].buildFromNewick(newicks[i-1], true, false); // newicks starts at 0
            _forests[i].refreshPreorder();
            _forests[i]._theta_mean = Forest::_theta; // for now, set theta mean equal to whatever theta user specifies
        }
    }

    inline void Particle::resetSpecies() {
        _forests[0].clear();
        _generation = 0;
        setLogWeight(0.0);
        _species_tree_height = 0.0;
        _log_coalescent_likelihood = 0.0;
        for (int i=1; i<_forests.size(); i++) {
            _forests[i]._log_coalescent_likelihood = 0.0;
            _forests[i]._data = nullptr;
            _forests[i]._log_weight = 0.0;
        }
        _t.clear();
        _forests[0]._last_edge_length = 0.0;
        _forests[0]._increments_and_priors.clear();
    }

    inline void Particle::setForest(Forest f, unsigned forest_number) {
        _forests[forest_number] = f;
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
        _deep_coal = other._deep_coal;
        _species_tree_height = other._species_tree_height;
        _t = other._t;
        _psuffix = other._psuffix;
    };
}

