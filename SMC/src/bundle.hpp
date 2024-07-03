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
            double                                          getBundleLogWeight(){return _bundle_log_weight;}
            void                                            setBundleLogWeight(double w) {_bundle_log_weight = w;}
            void                                            clearPartials();
            void                                            resetSpecies();
            void                                            proposeSpeciesParticles();
            void                                            deleteExtraGeneParticles();
        
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
                _gene_particles[g][i] = Particle(*new Particle); // TODO: why * ?
                _gene_particles[g][i].setData(d, taxon_map, partials, "gene", index);
            }
            index++;
        }
        _species_particle = Particle(*new Particle); // TODO: why * ?
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
//        unsigned nsteps = Forest::_ntaxa - 1;
        
        if (_generation == 0) {

            // TODO: draw a species increment to start with - for now, build entire species tree from the start
//            _species_particle.drawFirstSpeciesIncrement();
            
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
        
        // TODO: erase / break down species tree as far as possible
//        cout << "stop";
        vector<unsigned> current_species;
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                current_species.push_back(_gene_particles[g][i]._current_species);
            }
        }
        
        unsigned max_current_species = *max_element(current_species.begin(), current_species.end());
        if (max_current_species > 0) {
//            cout << "stop";
            // TODO: un-join species tree back to this point
            // TODO: erase last elements of _t in species tree
            // TODO: rebuild the species tree after that point
            // TODO: update _t in all the gene particles after the rebuilding point
        }
            
            unsigned psuffix = 1;
            
            for (unsigned g=0; g<_ngenes; g++) {
                for (unsigned i=0; i<_ngene_particles; i++) {
                    
                    // set particle random number seeds
                    _gene_particles[g][i].setSeed(rng.randint(1,9999) + psuffix);
                    psuffix += 2;
                    
                    _gene_particles[g][i].proposal(); // TODO: make this random, not sequential
//                    if (_gene_particles[g][i]._speciation_event_next) {
//                        _species_particle.updateSpeciesTree();
//                    }

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
//            for (unsigned i=1; i<_ngene_particles; i++) {
//                _gene_particles[g][i].clear();
//            }
        }
        _ngene_particles = 1;
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
