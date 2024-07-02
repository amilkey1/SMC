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
        void                                            filterLocus(unsigned g);
        void                                            normalizeWeights(vector<Particle::SharedPtr> particles);
        double                                          getRunningSum(vector<double> log_weight_vec);
        void                                            resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
        void                                            resetWeights(vector<Particle::SharedPtr> & particles);
        vector<pair<double, double>>                    getSpeciesTreeIncrementPriors();
        string                                          saveSpeciesNewick(){return _species_particle->saveForestNewick();}
        string                                          saveGeneNewick(unsigned i){return _gene_particles[i][0]->saveForestNewick();} // TODO: for now, just saving first gene tree - later, save all of them?
    
    private:
    
        vector<vector<Particle::SharedPtr>>                             _gene_particles;
        vector<vector<Particle::SharedPtr>>                             _my_vec_1;
        vector<vector<Particle::SharedPtr>>                             _my_vec_2;
        Particle::SharedPtr                                             _species_particle;
        mutable                                                         Lot::SharedPtr _lot;
        static unsigned                                                 _ngenes;
        static unsigned                                                 _ngene_particles;
        vector<double>                                                  _log_marginal_likelihood_by_gene;
        bool                                                            _use_first;
        double                                                          _small_enough;
};

    inline Bundle::Bundle() {
        _lot.reset(new Lot());
        clear();
    };

    inline void Bundle::clear() {
        _ngenes = 1;
        _ngene_particles = 1;
        _gene_particles.clear();
        _species_particle = nullptr;
        _use_first = true;
        _my_vec_1.clear();
        _my_vec_2.clear();
        _small_enough = 0.0000001;
    }

    inline void Bundle::setData(Data::SharedPtr d, map<string, string> &taxon_map, bool partials) {
        _gene_particles.resize(_ngenes);
        for (auto &g:_gene_particles) {
            g.resize(_ngene_particles);
        }
        
        _my_vec_1 = _gene_particles;
        _my_vec_2 = _gene_particles;
        
        unsigned index = 1;
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                _gene_particles[g][i] = Particle::SharedPtr(new Particle);
                _gene_particles[g][i]->setData(d, taxon_map, partials, "gene", index);
            }
            index++;
        }
        _species_particle = Particle::SharedPtr(new Particle);
        _species_particle->setData(d, taxon_map, partials, "species", index);
    }

    inline void Bundle::mapSpecies(map<string, string> &taxon_map, vector<string> &species_names) {
        for (unsigned g=0; g<_ngenes; g++) {
            for (unsigned i=0; i<_ngene_particles; i++) {
                _gene_particles[g][i]->mapSpecies(taxon_map, species_names, "gene");
            }
        }
        _species_particle->mapSpecies(taxon_map, species_names, "species");
    }

    inline void Bundle::runBundle() {
        unsigned nsteps = Forest::_ntaxa - 1;

            // TODO: draw a species increment to start with
            // TODO: for now, build entire species tree from the start
    //        _species_particle->drawFirstSpeciesIncrement();
            _species_particle->setLot(_lot);
            _species_particle->buildEntireSpeciesTree();
        // TODO: check gene forests are consistent with species increments
            // each gene particle has its own _t since particles may be at different places in the species tree
            // also set a starting likelihood for each gene particle  - // TODO: can copy some of this later
            for (unsigned g=0; g<_ngenes; g++) {
                for (unsigned i=0; i<_ngene_particles; i++) {
                    _gene_particles[g][i]->setLot(_lot); // TODO: set lot or rnseed?
                    _gene_particles[g][i]->_t = _species_particle->_t;
                    _gene_particles[g][i]->calcLogLikelihood();
                }
            }
            
        
        for (unsigned n=0; n<nsteps; n++) {
//            cout << "beginning step " << n << endl;
            _log_marginal_likelihood_by_gene.clear();
            
            unsigned psuffix = 1;
            
            for (unsigned g=0; g<_ngenes; g++) {
                for (unsigned i=0; i<_ngene_particles; i++) {
                    
                    // set particle random number seeds
                    _gene_particles[g][i]->setSeed(rng.randint(1,9999) + psuffix);
                    psuffix += 2;
                    
                    _gene_particles[g][i]->proposal(); // TODO: make this random, not sequential

                }
            }
            
            for (unsigned g=0; g<_ngenes; g++) {
                filterLocus(g);
            }
//            _gene_particles[0][0]->showParticle();
        }
    }

    inline void Bundle::filterLocus(unsigned g) {
//        vector<Particle::SharedPtr> particles = _gene_particles[g];
        normalizeWeights(_gene_particles[g]);
        
        resampleParticles(_gene_particles[g], _use_first ? _my_vec_2[g]:_my_vec_1[g]);
        //if use_first is true, my_vec = my_vec_2
        //if use_first is false, my_vec = my_vec_1
        
        _gene_particles[g] = _use_first ? _my_vec_2[g]:_my_vec_1[g];
        
        //change use_first from true to false or false to true
        _use_first = !_use_first;
        
        resetWeights(_gene_particles[g]);
    }

    inline void Bundle::normalizeWeights(vector<Particle::SharedPtr> particles) {
        unsigned i = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            log_weight_vec[i++] = p->getLogWeight();
        }

        double log_particle_sum = getRunningSum(log_weight_vec);

        for (auto & p : particles) {
            p->setLogWeight(p->getLogWeight() - log_particle_sum);
        }

    //    _log_marginal_likelihood_by_gene += log_particle_sum - log(_nparticles);
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

    inline void Bundle::resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles) {
         assert (from_particles.size() == _ngene_particles);
         assert (to_particles.size() == _ngene_particles);

         vector<pair<double, double>> cum_probs;
             // Create vector of pairs p, with p.first = log weight and p.second = particle index
         cum_probs.resize(_ngene_particles);
         unsigned i = 0;

        for (unsigned p=0; p < _ngene_particles; p++) {
             cum_probs[i].first = from_particles[p]->getLogWeight();
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
             assert( fabs( 1.0 - cum_probs[_ngene_particles-1].first ) < _small_enough);


         // Draw new set of particles by sampling with replacement according to cum_probs
         to_particles.resize(_ngene_particles);
         for(unsigned i = 0; i < _ngene_particles; i++) {

             // Select a particle to copy to the ith slot in to_particles
             int sel_index = -1;
             double u = rng.uniform();
             for(unsigned j = 0; j < _ngene_particles; j++) {
                 if (u < cum_probs[j].first) {
                     sel_index = cum_probs[j].second;
                     break;
                 }
             }
             assert(sel_index > -1);

             Particle::SharedPtr p0 = from_particles[sel_index];
             to_particles[i]=Particle::SharedPtr(new Particle(*p0));

             assert(_ngene_particles == to_particles.size());
             }
        }

    inline void Bundle::resetWeights(vector<Particle::SharedPtr> & particles) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p->setLogWeight(logw);
        }
    }

    inline vector<pair<double, double>> Bundle::getSpeciesTreeIncrementPriors() {
        return _species_particle->getSpeciesTreeIncrementPriors();
    }

    inline void Bundle::operator=(const Bundle & other) {
        _gene_particles = other._gene_particles;
        _species_particle = other._species_particle;
        _ngenes = other._ngenes;
        _ngene_particles = other._ngene_particles;
        _log_marginal_likelihood_by_gene = other._log_marginal_likelihood_by_gene;
        _my_vec_1 = other._my_vec_1;
        _my_vec_2 = other._my_vec_2;
        _use_first = other._use_first;
    };

}
