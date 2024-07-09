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
    friend class Bundle;
    
    public:

        Particle();
        Particle(const Particle & other);
        typedef std::shared_ptr<Particle>               SharedPtr;

        void                                    debugParticle(std::string name);
        void                                    showParticle();
        void                                    showSpeciesParticle();
        void                                    proposal();
        void                                    setData(Data::SharedPtr d, map<string, string> &taxon_map, bool partials, string type, unsigned index) {
                                                    _type = type;
                                                    _data = d;
                                                    if (_type == "gene") {
                                                        _forest.setData(d, index, taxon_map, partials); // TODO: go to one forest, not a vector
                                                    }
                                                }
        void                                    mapSpecies(map<string, string> &taxon_map, vector<string> &species_names, string type);
        double                                  addNewIncrement();
        void                                    saveForest(std::string treefilename);
        void                                    savePaupFile(std::string paupfilename, std::string datafilename, std::string treefilename, double expected_lnL) const;
        double                                  calcLogLikelihood();
        double                                  calcConstrainedWeightFactor(double max_depth);
        double                                  getLogLikelihood();
        double                                  calcGeneTreeLogLikelihood();
        double                                  calcHeight();
        double                                  getLogWeight() const {return _log_weight;}
        void                                    setLogWeight(double w){_log_weight = w;}
        void                                    setLogLikelihood(double forest_likelihood);
        void                                    setLogCoalescentLikelihood(double coalescent_like);
        void                                    operator=(const Particle & other);
        const Forest &                          getForest() const {return _forest;}
        string                                  saveForestNewick() {
            return _forest.makeNewick(8, true);}
            
        string                             saveGeneNewick(unsigned i) {
            return _forest.makeNewick(8, true);}
    
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
        int                                             selectEvent(vector<double> weight_vec);
        double                                          getTopologyPrior(unsigned i);
        vector<pair<double, double>>                    getIncrementPriors(unsigned i);
        vector<pair<double, double>>                    getSpeciesTreeIncrementPriors();
        double                                          getCoalescentLikelihood(unsigned g);
        void                                            clear();
        vector<double>                                  chooseIncrements(vector<double> event_choice_rates);
        double                                          speciesOnlyProposal(vector<double> max_depth_vector);
        double                                          calcLogCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height);
        tuple<string, string, string>                   speciesJoinProposal();
        double                                          getMaxDepth(tuple<string, string, string> species_joined, bool reset);
        void                                            setUpSpeciesOnlyProposal();
        void                                            clearPartials();
        Lot::SharedPtr getLot() const {return _lot;}
        void setSeed(unsigned seed) const {_lot->setSeed(seed);}
        double                                          getSpeciesTreeHeight();
        double                                          getSpeciesTreeLength();
        double                                          getGeneTreeHeight();
        double                                          getGeneTreeLength();
        double                                          getGeneTreeLogLikelihood();
        vector<double>                                  getGeneTreePriors();
        inline vector<double>                           getGeneTreeCoalescentLikelihoods();
        double                                          getSpeciesTreePrior();
        void                                            resetSpecies();
        Forest                                          getForest() {return _forest;} // TODO: should return a pointer?
        pair<string, string>                            getSpeciesJoined(){return make_pair(_forest._species_joined.first->_name, _forest._species_joined.second->_name);}
        void                                            setPsuffix(unsigned psuffix) {_psuffix = psuffix;}
        void                                            updateSpeciesPartition(pair<tuple<string, string, string>, double> species_info);
        void                                            buildEntireSpeciesTree();
        void                                            rebuildEntireSpeciesTree();
        void                                            setLot(Lot::SharedPtr lot) {_lot = lot;}
        void                                            updateSpeciesTree();
        void                                            unJoinSpecies(unsigned max_current_species);
        void                                            trimSpeciesTree(double amount_to_trim);
        void                                            drawTheta();
        vector<string>                                  getExistingSpeciesNames();
        void                                            rebuildThetaMap();
    
    private:

        static unsigned                         _nsubsets;
        Forest                                  _forest;
        double                                  _log_weight;
        Data::SharedPtr                         _data;
        double                                  _log_likelihood;
        int                                     _generation = 0;
        double                                  _log_coalescent_likelihood;
        mutable                                 Lot::SharedPtr _lot;
        unsigned                                _num_deep_coalescences;
        unsigned                                _psuffix;
        vector<pair<tuple<string, string, string>, double>> _t;
        string                                  _type;
        unsigned                                _current_species;
        double                                  _constrained_factor;
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
        cout << "\n";
        
#if defined (DRAW_NEW_THETA)
        unsigned a = 1;
        for (auto &t:_forest._theta_map) {
            cout << "theta " << a << " = " << t.second << endl;
            a++;
        }
#endif
        cout << "\n";
        _forest.showForest();
    }

    inline void Particle::showSpeciesParticle() {
        assert (_type == "species");
        //print out weight of each species forest part of the particle
        cout << "\nParticle:\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _forest.showForest();
    }

    inline void Particle::clear() {
        _log_weight     = 0.0;
        _log_likelihood = 0.0;
        _forest.clear();
        _data           = nullptr;
        _nsubsets       = 0;
        _generation     = 0;
        _log_coalescent_likelihood = 0.0;
        _num_deep_coalescences = 0.0;
        _t.clear();
        _psuffix = 0;
        _type = "null";
        _current_species = 0;
        _constrained_factor = 0.0;
    }

    inline void Particle::showSpeciesTree() {
        //print out weight of each particle
        cout << "\nParticle:\n";
        cout << "  _log_weight: " << _log_weight << "\n" ;
        cout << " _log_likelihood: " << _log_likelihood << "\n";
        cout << "  _forest: " << "\n";
        cout << "\n";
        _forest.showForest();
    }

    //more detailed version of showParticle
    inline void Particle::debugParticle(std::string name) {
        cout << "debugging particle" << endl;
        //print out weight of each particle
        cout << "\nParticle " << name << ":\n";
            cout << "  _log_weight:               " << _log_weight                 << "\n" ;
            cout << "  _log_likelihood:           " << _log_likelihood             << "\n";
            cout << "  _forest._nleaves:          " << _forest._nleaves            << "\n";
            cout << "  _forest._ninternals:       " << _forest._ninternals         << "\n";
            cout << "  _forest._npatterns:        " << _forest._npatterns          << "\n";
            cout << "  _forest._nstates:          " << _forest._nstates            << "\n";
            cout << "  _forest._last_edge_length: " << _forest._last_edge_length   << "\n";
            cout << "  newick description:        " << _forest.makeNewick(5,false) << "\n";
    }

    inline double Particle::calcLogLikelihood() {
        //calculate likelihood for each gene tree
        double log_likelihood = _forest.calcLogLikelihood();
        assert(!isnan (log_likelihood));
        if (_generation == 0 && !Forest::_run_on_empty) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline double Particle::calcGeneTreeLogLikelihood() {
        return _forest.calcLogLikelihood();
    }

    inline void Particle::clearPartials() {
        _forest.clearPartials();
    }

    inline void Particle::setUpSpeciesOnlyProposal() {
        assert (_generation == 0);
        assert (_type == "gene");
        _forest.refreshPreorder();
        _forest.calcMinDepth();
        _forest._nincrements = 0;
    }

    inline tuple<string, string, string> Particle::speciesJoinProposal() {
        assert (_type == "species");
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        
        if (_forest._last_edge_length > 0.0) {
        // choose species to join if past the first species generation for each forest vector
            species_joined = _forest.speciesTreeProposal(_lot);
        }
        
        return species_joined;
    }

    inline double Particle::getMaxDepth(tuple<string, string, string> species_joined, bool reset) {
        assert (_type == "gene");
        
        string species1 = get<0>(species_joined);
        string species2 = get<1>(species_joined);
        
        if (species1 != "null") {
            _forest.updateSpeciesPartition(species_joined); // if using Jones formula, this will happen in coalescent likelihood calculation
        }
        
        if (reset && species1 != "null") {
            _forest.resetDepthVector(species_joined);
        }
        double max_depth = (_forest.getMinDepths())[0].first;
        
        return max_depth;
    }

    inline double Particle::speciesOnlyProposal(vector<double> max_depth_vector) {
        assert (_type == "species");
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        double prev_log_coalescent_likelihood = _log_coalescent_likelihood;
        
                double max_depth = 0.0;
        
                if (_forest._lineages.size() > 1) {
                    max_depth = *min_element(max_depth_vector.begin(), max_depth_vector.end());
                    max_depth -= _forest.getTreeHeight();
                    // choose a species tree increment
                }
                
                if (_forest._lineages.size() > 1) {
    #if !defined (UNCONSTRAINED_PROPOSAL)
                    assert (max_depth > 0.0);
                    _forest.chooseSpeciesIncrementOnly(_lot, max_depth);
    #else
                    _forest.chooseSpeciesIncrementOnly(_lot, 0.0);
    #endif
                }
                if (_forest._lineages.size() == 1) {
                    _forest._last_edge_length = 0.0;
                }
        
                _t.push_back(make_pair(species_joined, _forest._last_edge_length));
        
        _constrained_factor = 0.0;
#if !defined (UNCONSTRAINED_PROPOSAL)
        if (!Forest::_run_on_empty) {
            assert (max_depth > 0.0);
            double nlineages = _forest._lineages.size();
            if (nlineages == 1) {
                nlineages = 2; // for last step, constraint was before final two species were joined
            }
            _constrained_factor = log(1 - exp(-1*nlineages*Forest::_lambda*max_depth));
        }
        
#endif
//        _log_species_weight = _log_coalescent_likelihood - prev_log_coalescent_likelihood + constrained_factor;
        
        _generation++;
        return max_depth;
    }

    inline double Particle::calcConstrainedWeightFactor(double max_depth) {
        assert (max_depth > 0.0);
        double nlineages = _forest._lineages.size();
        double constrained_factor = log(1 - exp(-1*nlineages*Forest::_lambda*max_depth));
        return constrained_factor;
    }

    inline double Particle::calcLogCoalescentLikelihood(double species_increment, tuple<string, string, string> species_joined, double species_tree_height) {
        _forest.calcCoalescentLikelihood(species_increment, species_joined, species_tree_height);
        _log_coalescent_likelihood = _forest._log_coalescent_likelihood + _forest._panmictic_coalescent_likelihood;
        return _log_coalescent_likelihood;

    }

    inline void Particle::setLogCoalescentLikelihood(double log_coal_like) {
        _log_coalescent_likelihood = log_coal_like;
    }

    inline void Particle::setLogLikelihood(double forest_log_likelihood) {
        _log_likelihood = forest_log_likelihood;
        _log_weight = forest_log_likelihood;
    }

    inline double Particle::getSpeciesTreeHeight() {
        assert (_type == "species");
        return _forest.getTreeHeight();
    }

    inline double Particle::getSpeciesTreeLength() {
        assert (_type == "species");
        return _forest.getTreeLength();
    }

    inline double Particle::getGeneTreeHeight() {
        return _forest.getTreeHeight();
    }

    inline double Particle::getGeneTreeLength() {
        return _forest.getTreeLength();
    }

    inline double Particle::getGeneTreeLogLikelihood() {
        return _forest._gene_tree_log_likelihood;
    }

    inline vector<double> Particle::getGeneTreePriors() {
        vector<double> gene_tree_priors;
        double prior = 0.0;
        for (auto &p:_forest._increments_and_priors) {
            prior += p.second;
        }
        
        gene_tree_priors.push_back(prior);
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
        gene_tree_priors.push_back(_forest._log_coalescent_likelihood);
        return gene_tree_priors;
#endif
    }

    inline double Particle::getSpeciesTreePrior() {
        double prior = 0.0;
        for (auto &p:_forest._increments_and_priors) {
            prior += p.second;
        }
        prior += _forest._log_joining_prob;
        
        return prior;
    }

    inline double Particle::getLogLikelihood() {
        //retrieve likelihood for each gene tree
        double log_likelihood = 0.0;
        double gene_tree_log_likelihood = _forest._gene_tree_log_likelihood;
        assert(!isnan (log_likelihood));
        //total log likelihood is sum of gene tree log likelihoods
        log_likelihood += gene_tree_log_likelihood;
        
        if (_generation == 0 && !Forest::_run_on_empty) {
            _log_weight = log_likelihood;
        }

        return log_likelihood;
    }

    inline void Particle::proposal() {
        double inv_gamma_modifier = 0.0;
        bool done = false;
                
        while (!done) {
            
            vector<pair<double, string>> rates_by_species = _forest.calcForestRate(_lot);
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
            
            double species_increment = _t[_current_species].second;

            if ((gene_increment < species_increment || _forest._species_partition.size() == 1) && gene_increment != -1.0) {
                // if the species forest is done, don't choose a speciation event
                assert (gene_increment > 0.0);
                
                _forest.addIncrement(gene_increment);
                
                vector<double> event_choice_rates;
                for (auto &r:rates_by_species) {
                    event_choice_rates.push_back(r.first);
                }
                for (auto &p:event_choice_rates) {
                     p = log(p/total_rate);
                 }

                unsigned index = selectEvent(event_choice_rates);
                string species_name = rates_by_species[index].second;
                _forest.allowCoalescence(species_name, gene_increment, _lot);
                
                if (_forest._species_partition.size() > 1) { // otherwise, species tree is done and there is no need to update it further
                    _t[_current_species].second -= gene_increment;
                }
                done = true;
            }
            else {
                // carry out speciation event
                
                assert (_forest._species_partition.size() > 1);
                _forest.addIncrement(species_increment);
                _forest.updateSpeciesPartition(_t[_current_species+1].first);
                assert (_current_species < _t.size());
                _t[_current_species].second -= species_increment;
                assert (_t[_current_species].second == 0.0);

                if (_forest._species_partition.size() > 1) {
                    _current_species++;
                }
                
#if defined (INV_GAMMA_PRIOR_TWO) // modifier each time a new theta is drawn
                double theta_mean = _forest._theta_mean;
                double eps = 0.01;
                double a = 2.0;
                
                inv_gamma_modifier = lgamma(a + eps) - lgamma(a) + a * log(a-1.0) - (a + eps) * log(a + eps - 1.0);
                
                // new theta drawn for the new population
                string new_spp = get<2> (_t[_current_species].first);
                double y = _forest._theta_map[new_spp];
                assert (y > 0.0);
            
                inv_gamma_modifier += eps * (log(y) - log(theta_mean)) + (theta_mean * eps / y);
#endif
                }
            }
        
#if defined (INV_GAMMA_PRIOR_TWO)
            if (_generation == 0) { // modifier for each theta drawn for each tip // TODO: is this for every gene or just once for the whole particle? if just once, how does it affect the weights?
                inv_gamma_modifier = 0.0;
                double theta_mean = _forest._theta_mean;
                double eps = 0.01;
                double a = 2.0;
                
                inv_gamma_modifier = lgamma(a + eps) - lgamma(a) + a * log(a-1.0) - (a + eps) * log(a + eps - 1.0);
                
                // for first step, new theta drawn for each tip population
                for (int i=0; i<Forest::_nspecies; i++) {
                    auto it = _forest._theta_map.begin();
                    std::advance(it, i);
                    double y = it->second;
                    
                    inv_gamma_modifier += eps * (log(y) - log(theta_mean)) + (theta_mean * eps / y);
                }
            }
#endif
        
        _log_weight = _forest._log_weight;
        
#if defined (DRAW_NEW_THETA)
        // modifier only happens on first round
        if (_generation == 0 && _forest._theta_prior_mean > 0.0 && _forest._theta_proposal_mean > 0.0) {
            double prior_rate = 1.0/_forest._theta_prior_mean;
            double proposal_rate = 1.0/_forest._theta_proposal_mean;
            double log_weight_modifier = log(prior_rate) - log(proposal_rate) - (prior_rate - proposal_rate)*_forest._theta_mean;
            
            _log_weight += log_weight_modifier;
        }
#endif
        
        _log_weight += inv_gamma_modifier;
        
        _generation++;
    }

    inline void Particle::updateSpeciesTree() {
        assert (_forest._lineages.size() > 1);
        tuple<string, string, string> species_joined = _forest.speciesTreeProposal(_lot);
        _forest.chooseSpeciesIncrementOnly(_lot, 0.0);
        double edge_len = _forest._last_edge_length;
        _t.push_back(make_pair(species_joined, edge_len));
    }

    inline void Particle::unJoinSpecies(unsigned max_species_count) {
        unsigned nsteps_to_revert = Forest::_nspecies - max_species_count - 1; // -1 because gene trees use current_species + 1 when updating species partitions
        for (unsigned i=0; i<nsteps_to_revert; i++) {
            _forest.revertSpeciesTreeOneJoin(_t.back().second);
            _t.pop_back(); // edge length to use is one back from species join (going through _t backwards)?
            _forest._ninternals--;
        }
    }

    inline vector<string> Particle::getExistingSpeciesNames() {
        assert (_type == "species");
        return _forest.getExistingSpeciesNames();
    }

    inline void Particle::rebuildThetaMap() {
//        assert (_type == "gene");
        _forest.rebuildThetaMap(_lot);
    }

    inline void Particle::drawTheta() {
        // create map for one forest, then copy it to all forests
        assert (_type == "gene");
        _forest.createThetaMap(_lot);
    }

    inline void Particle::trimSpeciesTree(double amount_to_trim) {
        assert (_type == "species");
//        double current_height = _forest.getTreeHeight();
//        double amount_to_trim = current_height - deepest_coalescent_event;
//        if (amount_to_trim > 0.0) {
            _forest.trimTree(amount_to_trim);
            _t.back().second -= amount_to_trim;
        
//            assert (fabs(_forest.getTreeHeight() - deepest_coalescent_event) < 0.0000001);
//        }
        // else, the gene trees are below the species tree
//        return amount_to_trim;
    }

    vector<double> Particle::chooseIncrements(vector<double> event_choice_rates) {
        vector<double> increments;
        increments.resize(event_choice_rates.size());
        
        for (int p=0; p<event_choice_rates.size(); p++) {
            increments[p] = -log(1.0 - _lot->uniform())/event_choice_rates[p];
        }
        return increments;
    }

    inline double Particle::addNewIncrement() {
        _forest.chooseSpeciesIncrementOnly(_lot, 0.0);
        _t.back().second += _forest._last_edge_length;
        
        return _forest._last_edge_length;
    }

    inline void Particle::updateSpeciesPartition(pair<tuple<string, string, string>, double> species_info) {
        _forest.updateSpeciesPartition(species_info.first);
    }

    inline void Particle::rebuildEntireSpeciesTree() {
        unsigned current_size = (unsigned) _t.size();
        double edge_len = 0.0;
        
        for (unsigned i=current_size; i < _forest._nspecies; i++) {
            if (_forest._lineages.size() > 2) {
                tuple<string, string, string> species_joined = _forest.speciesTreeProposal(_lot);
                
                _forest.chooseSpeciesIncrementOnly(_lot, 0.0);
                edge_len = _forest._last_edge_length;
                _t.push_back(make_pair(species_joined, edge_len));
            }
            else {
                // join remaining species
                tuple<string, string, string> species_joined = _forest.speciesTreeProposal(_lot);
                edge_len = 0.0;
                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
        assert (_forest._lineages.size() == 1);
    }
    inline void Particle::buildEntireSpeciesTree() {
        _forest.chooseSpeciesIncrementOnly(_lot, 0.0);
        double edge_len = _forest._last_edge_length;
        
        tuple<string, string, string> species_joined = make_tuple("null", "null", "null");
        _t.push_back(make_pair(species_joined, edge_len));

        for (unsigned i=0; i < _forest._nspecies-1; i++) {
            if (_forest._lineages.size() > 1) {
                species_joined = _forest.speciesTreeProposal(_lot);
                
                double edge_len = 0.0;
                if (_forest._lineages.size() > 1) {
                    _forest.chooseSpeciesIncrementOnly(_lot, 0.0);
                    edge_len = _forest._last_edge_length;
                }
                _t.push_back(make_pair(species_joined, edge_len));
            }
        }
    }

    inline Particle::Particle(const Particle & other) {
        _lot.reset(new Lot());
        *this = other;
    }

    inline void Particle::saveForest(std::string treefilename)  {
        ofstream treef(treefilename);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        treef << "  tree test = [&R] " << _forest.makeNewick(8, true)  << ";\n";
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
        Node* base_node = _forest._lineages[0];
        sum_height += base_node->getEdgeLength();
        for (Node* child=base_node->_left_child; child; child=child->_left_child) {
            sum_height += child->getEdgeLength();
        }
        return sum_height;
    }

    inline void Particle::setNumSubsets(unsigned n) {
        _nsubsets = n;
    }

    inline void Particle::mapSpecies(map<string, string> &taxon_map, vector<string> &species_names, string type) {
        if (type == "species" ) {
            //species tree
            _forest.setUpSpeciesForest(species_names);
        }
        else {
            // gene trees
            _forest.setUpGeneForest(taxon_map);
        }
    }

    inline void Particle::showSpeciesIncrement(){
        assert (_type == "species");
        cout << "species tree increment: " << "     " << _forest._last_edge_length << endl;
    }

    inline void Particle::showSpeciesJoined(){
        assert (_type == "species");
        _forest.showSpeciesJoined();
    }

    inline int Particle::selectEvent(vector<double> weight_vec) {
        // choose a random number [0,1]
        double u =  _lot->uniform();
        assert (u > 0.0);
        assert (u < 1.0);
        double cum_prob = 0.0;
        unsigned index = 0;
        for (int i=0; i < (int) weight_vec.size(); i++) {
            cum_prob += exp(weight_vec[i]); // TODO: can leave weights on linear scale
            if (u <= cum_prob) {
                index = i;
                break;
            }
        }
        // return index of choice
        return index;
    }

    inline double Particle::getTopologyPrior(unsigned i) {
        return _forest._log_joining_prob;
    }

    inline vector<pair<double, double>> Particle::getIncrementPriors(unsigned i) {
        return _forest._increments_and_priors;
    }

    inline vector<pair<double, double>> Particle::getSpeciesTreeIncrementPriors() {
        return _forest._increments_and_priors;
    }

    inline double Particle::getCoalescentLikelihood(unsigned g) {
        // TODO: don't need g
#if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
        return _log_coalescent_likelihood; // can't get coalescent likelihood separately for each gene tree
#else
        assert (g>0); // no coalescent likelihood for species tree
        return _forest._log_coalescent_likelihood;
#endif
    }

    inline void Particle::resetSpecies() {
        if (_type == "species" ) {
            _forest.clear();
            setLogWeight(0.0);
            _log_coalescent_likelihood = 0.0;
            _forest._last_edge_length = 0.0;
            _forest._increments_and_priors.clear();
            _forest._species_tree_height = 0.0;
#if defined (DRAW_NEW_THETA)
            _forest._theta_map.clear(); // clear old thetas
            
#endif
        }
        
        else {
            _forest._log_coalescent_likelihood = 0.0;
            _forest._data = nullptr;
            _forest._log_weight = 0.0;
        }
        
        _t.clear();
        _generation = 0;
    }

    inline void Particle::operator=(const Particle & other) {
        _log_weight     = other._log_weight;
        _log_likelihood = other._log_likelihood;
        _forest         = other._forest;
        _data           = other._data;
        _nsubsets       = other._nsubsets;
        _generation     = other._generation;
        _log_coalescent_likelihood = other._log_coalescent_likelihood;
        _num_deep_coalescences = other._num_deep_coalescences;
        _t = other._t;
        _psuffix = other._psuffix;
        _type = other._type;
        _current_species = other._current_species;
        _constrained_factor = other._constrained_factor;
    };
}

