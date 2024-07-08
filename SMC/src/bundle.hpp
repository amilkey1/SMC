#pragma once
#include <vector>
#include "boost/format.hpp"

using namespace std;
using namespace boost;

#include "lot.hpp"
#include "conditionals.hpp"
#include "particle.hpp"

extern proj::Lot rng;

    namespace proj {

    class Bundle {
        public:
        
            Bundle();
            Bundle(const Bundle & other);
            typedef std::shared_ptr<Bundle>               SharedPtr;
            Lot::SharedPtr getLot() const {return _lot;}
        
            void setSeed(unsigned seed) const {_lot->setSeed(seed);}
            void                                            clear();
            static void                                     setNumGenes(unsigned nsubsets){_ngenes = nsubsets;}
            void                                            operator=(const Bundle & other);
            void                                            setData(Data::SharedPtr d, map<string, string> &taxon_map, bool partials);
            void                                            mapSpecies(map<string, string> &taxon_map, vector<string> &species_names);
            void                                            setNGeneParticles(unsigned nparticles) {_ngene_particles = nparticles;}
            void                                            runBundle();
            void                                            filterLoci();
            void                                            filterSpecies();
            double                                          getRunningSum(vector<double> log_weight_vec);
            void                                            resetWeights(vector<Particle> & particles);
            vector<pair<double, double>>                    getSpeciesTreeIncrementPriors();
            string                                          saveSpeciesNewick(){return _species_particle.saveForestNewick();}
            string                                          saveGeneNewick(unsigned i){return _gene_particles[i][0].saveForestNewick();} // TODO: for now, just saving first gene tree - later, save all of them?
            void                                            showGeneNewick(unsigned i) {_gene_particles[i][0].showParticle();} // TODO: for now, just saving first gene tree - later, save all of them?
            double                                          getBundleLogWeight(){return _bundle_log_weight;}
            void                                            setBundleLogWeight(double w) {_bundle_log_weight = w;}
            void                                            clearPartials();
            void                                            resetSpecies();
            void                                            proposeSpeciesParticles();
            void                                            deleteExtraGeneParticles();
            vector<double>                                  getLogLikelihoods();
            void                                            drawTheta();
            vector<double>                                  getThetas();
            double                                          getThetaMean();
        
        private:
        
            vector<vector<Particle>>                                        _gene_particles;
            Particle                                                        _species_particle;
            mutable                                                         Lot::SharedPtr _lot;
            static unsigned                                                 _ngenes;
            static unsigned                                                 _ngene_particles;
            vector<double>                                                  _log_marginal_likelihood_by_gene;
            bool                                                            _use_first;
            double                                                          _small_enough;
            unsigned                                                        _generation;
            double                                                          _bundle_log_weight;
            vector<double>                                                  _prev_log_marginal_likelihood_by_gene;
            double                                                          _prev_log_coalescent_likelihood;
    };

    inline Bundle::Bundle() {
        _lot.reset(new Lot());
        clear();
    };

    inline void Bundle::clear() {
        _ngenes = 1;
        _ngene_particles = 1;
        _gene_particles.clear();
        _use_first = true;
        _small_enough = 0.0000001;
        _generation = 0;
        _bundle_log_weight = 0.0;
        _prev_log_marginal_likelihood_by_gene.clear();
        _log_marginal_likelihood_by_gene.clear();
        _prev_log_coalescent_likelihood = 0.0;
    }

    inline Bundle::Bundle(const Bundle & other) {
        _lot.reset(new Lot()); // TODO: reset lot or not?
        *this = other;
    }

    inline void Bundle::setData(Data::SharedPtr d, map<string, string> &taxon_map, bool partials) {
        _gene_particles.resize(_ngenes);
        for (auto &g:_gene_particles) {
            g.resize(_ngene_particles);
        }
        
        unsigned index = 1;
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                _gene_particles[g][i] = Particle(*new Particle);
                _gene_particles[g][i].setData(d, taxon_map, partials, "gene", index);
            }
            index++;
        }
        _species_particle = Particle(*new Particle);
        _species_particle.setData(d, taxon_map, partials, "species", index);
    }

    inline void Bundle::mapSpecies(map<string, string> &taxon_map, vector<string> &species_names) {
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                _gene_particles[g][i].mapSpecies(taxon_map, species_names, "gene");
            }
        }
        _species_particle.mapSpecies(taxon_map, species_names, "species");
    }

    inline void Bundle::runBundle() {
        // TODO: create full theta map, then erase and rebuild at each step as needed
//        unsigned nsteps = Forest::_ntaxa - 1;
        if (_generation == 0) {
            
            _species_particle.setLot(_lot);
            _species_particle.buildEntireSpeciesTree();

            // set a starting likelihood for each gene particle  - // TODO: can copy some of this later
            for (unsigned g=0; g<_ngenes; g++) {
                for (unsigned i=0; i<_ngene_particles; i++) {
                    _gene_particles[g][i].setLot(_lot); // TODO: set lot or rnseed?
                    _gene_particles[g][i]._t = _species_particle._t;
                    _gene_particles[g][i].calcLogLikelihood();
                }
            }
            
            _log_marginal_likelihood_by_gene.resize(_ngenes);
            _prev_log_marginal_likelihood_by_gene.resize(_ngenes);
        }
        
        else {
            // don't do this for first gen
//            unsigned species_tree_size = (unsigned) _species_particle._forest._lineages.size();
            // TODO: erase members of theta map that aren't needed
            // TODO: rebuild theta map as needed
            vector<string> existing_species = _species_particle.getExistingSpeciesNames();
            
            _species_particle.rebuildThetaMap();
//            _gene_particles[0][0].rebuildThetaMap(existing_species); // TODO: do this in species particle, then copy to gene particles
            
            map<string, double> theta_map = _gene_particles[0][0]._forest._theta_map;
            for (unsigned g=0; g<_ngenes; g++) {
                for (unsigned i=0; i<_ngene_particles; i++) {
                    _gene_particles[g][i]._forest._theta_map = theta_map;
                }
            }
        }
        
        vector<unsigned> current_species;
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                current_species.push_back(_gene_particles[g][i]._current_species);
            }
        }
        
        unsigned max_current_species = *max_element(current_species.begin(), current_species.end());
//        _species_particle.unJoinSpecies(max_current_species); // this function will unjoin species back to the max current species
        // i.e. if max_current_species = 1, this function will erase so only _t[0] through _t[1] remain
        // i.e. if max_current_species = 3, this function will erase so _t[0] through _t[3] remain
        if (_generation > 0) {
            _species_particle.rebuildEntireSpeciesTree();
        }
                    
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
//                unsigned t_size = (unsigned) _gene_particles[g][i]._t.size();
//                unsigned n_to_remove = t_size - max_current_species - 1; // -1 because gene trees use current_species + 1 when updating species partitions
                
//                for (unsigned a=0; a<n_to_remove; a++) {
//                    _gene_particles[g][i]._t.pop_back();
//                }
//                assert (_gene_particles[g][i]._t.size() > 0); // _t should never be completely removed because first increments always count
                
                for (unsigned s=max_current_species + 1; s < Forest::_nspecies; s++) {
                    _gene_particles[g][i]._t.push_back(_species_particle._t[s]);
                }
            }
        }
            
        unsigned psuffix = 1;
            
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                
                // set particle random number seeds
                _gene_particles[g][i].setSeed(rng.randint(1,9999) + psuffix);
                psuffix += 2;
                
                _gene_particles[g][i].proposal(); // TODO: make this random, not sequential

            }
        }
            
        filterLoci();
        
        
        for (unsigned g=0; g<_ngenes; g++) {
            resetWeights(_gene_particles[g]);
        }
        
        double log_weight = 0.0;
        for (auto &m:_log_marginal_likelihood_by_gene) {
            assert (m != 0.0);
            log_weight += m;
        }
        _bundle_log_weight = log_weight;
        
        
        if (_generation < Forest::_ntaxa - 2) {
            current_species.clear();
            for (unsigned g=0; g<_ngenes; g++) {
                for (unsigned i=0; i<_ngene_particles; i++) {
                    current_species.push_back(_gene_particles[g][i]._current_species);
                }
            }
            
            max_current_species = *max_element(current_species.begin(), current_species.end());
            
            vector<double> gene_tree_heights;
            for (unsigned g=0; g<_ngenes; g++) {
                for (unsigned i=0; i<_ngene_particles; i++) {
                    gene_tree_heights.push_back(_gene_particles[g][i].getGeneTreeHeight());
                }
            }
            
            double current_height = 0.0;
            for (int t=0; t<max_current_species+1; t++) {
                current_height += _species_particle._t[t].second;
            }
//           current_height = _species_particle._forest.getTreeHeight();
            double deepest_coalescent_event = *max_element(gene_tree_heights.begin(), gene_tree_heights.end());
            double amount_to_trim = current_height - deepest_coalescent_event;
            
            if (amount_to_trim > 0.0) {
            
                _species_particle.unJoinSpecies(max_current_species); // this function will unjoin species back to the max current species
                
//                vector<double> gene_tree_heights;
//                for (unsigned g=0; g<_ngenes; g++) {
//                    for (unsigned i=0; i<_ngene_particles; i++) {
//                        gene_tree_heights.push_back(_gene_particles[g][i].getGeneTreeHeight());
//                    }
//                }
                
    //            double max_height = *max_element(gene_tree_heights.begin(), gene_tree_heights.end());
    //            double amount_to_trim = _species_particle.trimSpeciesTree(max_height);
                _species_particle.trimSpeciesTree(amount_to_trim);
                double new_increment = _species_particle.addNewIncrement();
                
                double increment_difference = -amount_to_trim + new_increment;
                unsigned species_tree_size = (unsigned) _species_particle._t.size();
                
                // TODO: trim gene tree _t back to species_tree_size
                // TODO: then add new_increment to _t.back() in all gene trees
                
                for (unsigned g=0; g<_ngenes; g++) {
                    for (unsigned i=0; i<_ngene_particles; i++) {
                        unsigned t_size = (unsigned) _gene_particles[g][i]._t.size();
                        unsigned n_to_remove = t_size - species_tree_size; // -1 because gene trees use current_species + 1 when updating species partitions ???
                        
                        for (unsigned a=0; a<n_to_remove; a++) {
                            _gene_particles[g][i]._t.pop_back();
                        }
                        assert (_gene_particles[g][i]._t.size() == species_tree_size);
                        _gene_particles[g][i]._t.back().second += increment_difference;
                    }
                }
            }
        }
        _generation++;
    }

    inline void Bundle::filterLoci() {
        for (unsigned l = 0; l<_ngenes; l++) {
            // Copy log weights for all particles in this locus to prob vector
            vector<double> probs(_ngene_particles, 0.0);
            for (unsigned p=0; p<_ngene_particles; p++) {
                probs[p] = _gene_particles[l][p].getLogWeight();
            }

            // Normalize log_weights to create discrete probability distribution
            double log_sum_weights = getRunningSum(probs);
            transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});

            // Compute component of the log marginal likelihood due to this step
            _prev_log_marginal_likelihood_by_gene[l] = _log_marginal_likelihood_by_gene[l];
            _log_marginal_likelihood_by_gene[l] += log_sum_weights - log(_ngene_particles);

            // Compute cumulative probabilities
            partial_sum(probs.begin(), probs.end(), probs.begin());
            
            // Initialize vector of counts storing number of darts hitting each particle
            vector<unsigned> counts(_ngene_particles, 0);
                    
            // Throw _nparticles darts
            for (unsigned i = 0; i < _ngene_particles; ++i) {
                double u = _lot->uniform();
                auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
                assert(it != probs.end());
                unsigned which = (unsigned)std::distance(probs.begin(), it);
                counts[which]++;
            }
                
            // Copy particles
            // Locate first donor
            unsigned donor = 0;
            while (counts[donor] < 2) {
                donor++;
            }
            
            // Locate first recipient
            unsigned recipient = 0;
            while (counts[recipient] != 0) {
                recipient++;
            }
            
            // Count number of cells that need copying to
            unsigned nzeros = 0;
            for (unsigned i = 0; i < _ngene_particles; i++) {
                if (counts[i] == 0)
                    nzeros++;
            }
            
            while (nzeros > 0) {
                assert(donor < _ngene_particles);
                assert(recipient < _ngene_particles);
                _gene_particles[l][recipient] = _gene_particles[l][donor];
                counts[donor]--;
                counts[recipient]++;
                nzeros--;
                if (counts[donor] == 1) {
                    // Move donor to next slot with count > 1
                    donor++;
                    while (donor < _ngene_particles && counts[donor] < 2) {
                        donor++;
                    }
                }
                
                // Move recipient to next slot with count equal to 0
                recipient++;
                while (recipient < _ngene_particles && counts[recipient] > 0) {
                    recipient++;
                }
            }
        }
    }

    inline double Bundle::getRunningSum(vector<double> log_weight_vec) {
        double running_sum = 0.0;
        double log_particle_sum = 0.0;

        double log_max_weight = *max_element(log_weight_vec.begin(), log_weight_vec.end());
        for (auto & i:log_weight_vec) {
            running_sum += exp(i - log_max_weight);
        }
        log_particle_sum = log(running_sum) + log_max_weight;

        return log_particle_sum;
    }

    inline void Bundle::resetWeights(vector<Particle> & particles) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p.setLogWeight(logw);
        }
    }

    inline vector<pair<double, double>> Bundle::getSpeciesTreeIncrementPriors() {
        return _species_particle.getSpeciesTreeIncrementPriors();
    }

    inline void Bundle::clearPartials() {
        // TODO: for now, just condition on first particle of each gene
        for (unsigned g=0; g<_ngenes; g++) {
            _gene_particles[g][0].clearPartials();
        }
    }

    inline void Bundle::resetSpecies() {
        _generation = 0;
        // TODO: for now, just condition on first particle of each gene
        for (unsigned g=0; g<_ngenes; g++) {
            _gene_particles[g][0].resetSpecies();
        }
        _species_particle.resetSpecies();
    }
    
    inline void Bundle::proposeSpeciesParticles() {
        if (_generation == 0) {
            for (unsigned g=0; g<_ngenes; g++) {
                _gene_particles[g][0].setUpSpeciesOnlyProposal();
            }
        }
        
        unsigned psuffix = 1;
        
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                
                // set particle random number seeds
                _gene_particles[g][i].setSeed(rng.randint(1,9999) + psuffix);
                psuffix += 2;
            }
        }
        
        _species_particle.setSeed(rng.randint(1,9999) + psuffix);
                
        vector<double> max_depths;
        
        tuple<string, string, string> species_joined = _species_particle.speciesJoinProposal();
        
        bool reset = true;
        
        if (_species_particle._forest._lineages.size() == 1) {
            reset = false;
        }
        for (unsigned g=0; g<_ngenes; g++) {
            max_depths.push_back(_gene_particles[g][0].getMaxDepth(species_joined, reset));
        }
        
        // TODO: switch to jones coalescent likelihood or draw new thetas for second round of species only filtering
        double max_depth = _species_particle.speciesOnlyProposal(max_depths);
        
        double log_coalescent_likelihood = 0.0;
        for (unsigned g=0; g<_ngenes; g++) {
            log_coalescent_likelihood += _gene_particles[g][0].calcLogCoalescentLikelihood(_species_particle._forest._last_edge_length, species_joined, _species_particle._forest._species_tree_height); // TODO: double check species tree height - don't need it in both forest and particle
        }
        
        _bundle_log_weight = log_coalescent_likelihood - _prev_log_coalescent_likelihood;
        
        _prev_log_coalescent_likelihood = log_coalescent_likelihood;
        
#if !defined (UNCONSTRAINED_PROPOSAL)
            double test = 1/_bundle_log_weight;
            assert(test != -0); // assert coalescent likelihood is not -inf
        if (max_depth > 0.0) {
            _bundle_log_weight += _species_particle.calcConstrainedWeightFactor(max_depth);
        }
#endif
        
        if (_species_particle._forest._lineages.size() == 2) {
            _species_particle.speciesJoinProposal(); // join remaining species - no change in coalescent likelihood
        }
        _generation++;
        
    }
    
    inline void Bundle::deleteExtraGeneParticles() {
        for (unsigned g=0; g<_ngenes; g++) {
            _gene_particles[g].erase(_gene_particles[g].begin(), _gene_particles[g].end() - 1);
        }
        _ngene_particles = 1;
    }
    
    inline vector<double> Bundle::getLogLikelihoods() {
        vector<double> log_likelihoods;
            for (unsigned i=0; i<_ngene_particles; i++) { // TODO: this function return likelihoods by gene, not by particle
                for (unsigned g=0; g<_ngenes; g++) {
                    log_likelihoods.push_back(_gene_particles[g][i].getLogLikelihood());
            }
        }
        return log_likelihoods;
    }
    
    inline void Bundle::drawTheta() {
        // set seeds first
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                _gene_particles[g][i].setLot(_lot); // TODO: set lot or rnseed?
            }
        }
        _species_particle.setLot(_lot); // TODO: set lot or rnseed?
        
        // TODO: should probably make some functions to do this rather than setting _forest directly from bundle
        _gene_particles[0][0].drawTheta(); // create theta map, then copy to all particles
        
        double theta_mean = _gene_particles[0][0]._forest._theta_mean;
        double theta_proposal_mean = _gene_particles[0][0]._forest._theta_proposal_mean;
        double theta_prior_mean = _gene_particles[0][0]._forest._theta_prior_mean; // TODO: do this in speices particle, not first gene particle
        
        map<string, double> theta_map = _gene_particles[0][0]._forest._theta_map;
        map<string, unsigned> species_indices = _gene_particles[0][0]._forest._species_indices;
        string ancestral_spp_name = _gene_particles[0][0]._forest._ancestral_species_name;
        
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
//                if (g != 0 && i != 0) {
                if (g != 0 || i != 0) {
                        // don't need to do this for first particle
                        _gene_particles[g][i]._forest._theta_map = theta_map;
                        _gene_particles[g][i]._forest._species_indices = species_indices;
                        _gene_particles[g][i]._forest._theta_mean = theta_mean;
                        _gene_particles[g][i]._forest._theta_proposal_mean = theta_proposal_mean;
                        _gene_particles[g][i]._forest._theta_prior_mean = theta_prior_mean;
                        _gene_particles[g][i]._forest._ancestral_species_name = ancestral_spp_name;
                }
            }
        }
        
        _species_particle._forest._theta_map = theta_map;
        _species_particle._forest._species_indices = species_indices;
        _species_particle._forest._theta_mean = theta_mean;
        _species_particle._forest._theta_proposal_mean = theta_proposal_mean;
        _species_particle._forest._theta_prior_mean = theta_prior_mean;
        _species_particle._forest._ancestral_species_name = ancestral_spp_name;
    }
    
    inline vector<double> Bundle::getThetas() {
        vector<double> thetas;
        for (auto &t:_species_particle._forest._theta_map) {
            thetas.push_back(t.second);
        }
        return thetas;
    }
    
    inline double Bundle::getThetaMean() {
#if defined (DRAW_NEW_THETA)
        return _species_particle._forest._theta_mean;
#else
        return Forest::_theta;
#endif
    }

    inline void Bundle::operator=(const Bundle & other) {
//        _gene_particles = other._gene_particles;
        _species_particle = other._species_particle;
        _ngenes = other._ngenes;
        _ngene_particles = other._ngene_particles;
        _log_marginal_likelihood_by_gene = other._log_marginal_likelihood_by_gene;
        _use_first = other._use_first;
        _generation = other._generation;
        _bundle_log_weight = other._bundle_log_weight;
        _small_enough = other._small_enough;
        _prev_log_marginal_likelihood_by_gene = other._prev_log_marginal_likelihood_by_gene;
        _prev_log_coalescent_likelihood = other._prev_log_coalescent_likelihood;
        
        // copy particles // TODO: is this necessary?
        _gene_particles.resize(_ngenes);
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i <_ngene_particles; i++) {
                _gene_particles[g].resize(_ngene_particles);
                _gene_particles[g][i] = other._gene_particles[g][i];
            }
        }
    };

}
