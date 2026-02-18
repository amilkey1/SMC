#pragma once

#include <iostream>
#include "data.hpp"
#include "partition.hpp"
#include "stopwatch.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/distributions/gamma.hpp>
#include "xproj.hpp"
#include "particle.hpp"
#include <vector>
#include <thread>
#include <boost/algorithm/string/split.hpp>
#include "conditionals.hpp"
#include <algorithm>
#include <random>
#include <cstdio>
#include <mutex>
#include <random>

extern int my_rank;

extern void output(string msg);

using namespace std;
using namespace boost;
using namespace boost::algorithm;

#include "partial_store.hpp"
#include "g.hpp"
extern proj::PartialStore ps;
extern proj::Lot rng;
extern proj::Lot rng_mcmc;
extern proj::StopWatch stopwatch;

namespace proj {

    class Proj {
        public:

                                Proj();
                                ~Proj();

            void                clear();
            void                processCommandLineOptions(int argc, const char * argv[]);
            void                run();
            void                saveAllForests(vector<Particle> &v) const ;
            void                saveSpeciesTrees(vector<Particle> &v) const;
            void                saveAllSpeciesTrees(vector<Particle> &v) const;
            void                saveSpeciesTreesAfterFirstRound(vector<Particle> &v) const;
            void                saveSpeciesTreesHierarchical(vector<Particle> &v, string filename1, string filename2, unsigned group_number);
            void                saveSpeciesTreesAltHierarchical(vector<Particle> &v, unsigned group_number);
            void                saveGeneTrees(vector<Particle> &v) const;
            void                writeLoradFile(vector<Particle> &v) const;
            void                writeLoradFileAfterSpeciesFiltering(vector<Particle> &v) const;
            void                writeDeepCoalescenceFile(vector<Particle> &v);
            void                writeThetaFile(vector<Particle> &v);
            void                writeParamsFileForBeastComparison(vector<Particle> &v);
            void                writeParamsFileForBeastComparisonTestA(vector<Particle> &v, string filename) const;
            void                writeParamsFileForBeastComparisonTestB(vector<Particle> &v, string filename) const;
            void                writeParamsFileForBeastComparisonAfterSpeciesFiltering(vector<Particle> &v, string filename, unsigned group_number);
            void                writeParamsFileForBeastComparisonAfterSpeciesFilteringSpeciesOnly(vector<Particle> &v, string filename, unsigned group_number);
            void                writePartialCountFile(vector<Particle> &particles);
            void                createSpeciesMap(bool data);
            void                createSpeciesMapTyped();
            void                simSpeciesMap();
            string              inventName(unsigned k, bool lower_case);
            void                showFinal(vector<Particle> my_vec);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle> &particles);
            void                proposeParticleGroupRange(unsigned first, unsigned last, vector<vector<Particle>> &particles);
            void                proposeParticles(vector<Particle> &particles);
            void                mcmcMoves(vector<Particle> &particles, bool last_round);
            void                proposeMCMCMoveRange(unsigned first, unsigned last, vector<Particle> &particles, bool last_round);
            void                proposeParticlesSim(vector<Particle> &particles);
            void                proposeParticlesParallelizeByGroup(vector<vector<Particle>> &particles);
            void                simulateData();
            void                drawFromPrior();
            void                writePaupFile(vector<Particle> particles, vector<string> taxpartition);
            void                initializeParticle(Particle &particle);
            void                handleGeneNewicks();
            string              handleSpeciesNewick();
            string              readNewickFromFile(string file_name);
            void                buildNonzeroMap(vector<Particle> & particles, unsigned locus, map<const void *, list<unsigned> > & nonzero_map, const vector<unsigned> & nonzeros, vector<unsigned> particle_indices, unsigned start, unsigned end);
            double              filterParticles(unsigned step, vector<Particle> & particles, vector<unsigned> &particle_indices, unsigned start, unsigned end);
            void                filterParticlesThreading(vector<Particle> &particles, unsigned g, vector<unsigned> particle_indices);
            void                filterParticlesRange(unsigned first, unsigned last, vector<Particle> &particles, unsigned g, vector<unsigned> particle_indices);
            unsigned            multinomialDraw(const vector<double> & probs);
            double              filterSpeciesParticles(unsigned step, vector<Particle> & particles, unsigned id_number);
            double              computeEffectiveSampleSize(const vector<double> & probs) const;
        
#if defined (DRAW_NEW_THETA)
            void                updateSpeciesNames();
#endif
        
#if defined(SPECIES_IN_CONF)
        static void     parseSpeciesDefinition(string s);
#endif


        private:

            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood = 0.0;
            double                      _log_species_tree_marginal_likelihood = 0.0;
            unsigned                    _random_seed;
            unsigned long               _partials_needed;
            map<string, string>         _taxon_map;
        
            void                        summarizeData(Data::SharedPtr);
            double                      getRunningSum(const vector<double> &) const;
            void                        handleBaseFrequencies();
            void                        handleRelativeRates();
            void                        handleNTaxaPerSpecies();
            void                        checkOutgroupName();
        
#if defined (DEBUG_MODE)
            void                        debugSpeciesTree(vector<Particle> &particles);
#endif
        
            void                        secondLevel(vector<Particle> &particles);
            void                        buildSpeciesMap(bool taxa_from_data);
            void                        proposeSpeciesGroupsFaster(vector<Particle> &particle, unsigned ngroups, string filename1, string filename2, string filename3);
            void                        proposeSpeciesGroupRangeFaster(unsigned first, unsigned last, vector<Particle> &particles, unsigned ngroups, string filename1, string filename2, string filename3);
        
            double                      _small_enough;
            bool                        _first_line;
            unsigned                    _count; // counter for params output file
            vector<double>              _starting_log_likelihoods;
            vector<Lot::SharedPtr>      _group_rng;
            vector<vector<unsigned>>    _second_level_indices_to_keep; // particles to write to output files after second level
            vector<unsigned>            _particle_indices_to_thin;
    
            vector<pair<double, bool>>  _ranks;
            vector<vector<pair<double, bool>>> _gene_tree_ranks;
            vector<pair<double, bool>> _species_tree_ranks_after_first_round;
            vector<double>              _species_tree_heights_after_first_round;
            vector<pair<double, double>> _hpd_values;
            vector<vector<pair<double, double>>> _hpd_values_genes;
            vector<pair<double, double>>              _hpd_first_level_values;
            vector<double>              _species_tree_heights;
            vector<vector<double>>      _gene_tree_heights;
            vector<double>              _bhv_distances;
            vector<vector<double>>      _bhv_distances_genes;

    };

    inline Proj::Proj() {
//        std::cout << "Constructing a Proj" << std::endl;
        clear();
    }

    inline Proj::~Proj() {
//        std::cout << "Destroying a Proj" << std::endl;
    }

    inline void Proj::clear() {
        _partition.reset(new Partition());
        _data = nullptr;
        _small_enough = 0.0000001;
    }

    inline void Proj::saveAllForests(vector<Particle> &v) const {
        ofstream treef("forest.trees");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:v) {
            treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
        }
        treef << "end;\n";
        treef.close();
    }

    inline void Proj::writePartialCountFile(vector<Particle> &particles) {
        ofstream partialf("partial_count.txt");
        partialf << "total partials needed: ";
        
//        unsigned num_partials_needed = ((G::_ntaxa - 1) * G::_nparticles * G::_nloci * G::_ngroups);
        unsigned total_partials = 0;
        for (auto &p:particles) {
            total_partials += p.getPartialCount();
        }
        partialf << total_partials << "\n";
        
        partialf.close();
    }

    inline void Proj::writeParamsFileForBeastComparisonTestA(vector<Particle> &v, string filename) const {
        // this function creates a params file that is comparable to output from starbeast3
        ofstream logf(filename);
        logf << "iter ";
        logf << "\t" << "posterior ";
        logf << "\t" << "likelihood ";
        logf << "\t" << "prior ";
        logf << "\t" << "vectorPrior "; // log joint prior on population sizes (vectorPrior)
        logf << "\t" << "speciescoalescent ";
        logf << "\t" << "Tree.t:Species.height ";
        logf << "\t" << "Tree.t:Species.treeLength ";

        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "Tree.t:gene" + to_string(i) + "height";
            logf << "\t" << "Tree.t:gene" + to_string(i) + "treeLength";
        }

        logf << "\t" << "YuleModel.t:Species "; // this is the log probability of the species tree (multiply by log(3!) to get increment log prob)
        logf << "\t" << "popMean "; // this is psi in the InverseGamma(2,psi) distribution of popSize

        for (int i=0; i<(G::_nspecies*2-1); i++) {
            logf << "\t" << "popSize." + to_string(i+1);
        }

        logf << "\t" << "speciationRate.t:Species ";

        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "treeLikelihood:gene" + to_string(i);
        }
        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "treePrior:gene" + to_string(i);
        }
        logf << endl;

        int iter = 0;
        
        for (auto &p:v) {
                logf << iter;
                iter++;

                double vector_prior = 0.0; // no vector prior with jones coalescent likelihood

                double log_coalescent_likelihood = 0.0;
                for (unsigned g=1; g<G::_nloci+1; g++) {
                    log_coalescent_likelihood += p.getCoalescentLikelihood(g);
                }

                unsigned locus = p.getNextGene() - 1;
                double log_likelihood = p.calcLogLikelihoodLocus(locus, false);
                double log_prior = p.getAllPriorsFirstRound();

                double log_posterior = log_likelihood + log_prior + log_coalescent_likelihood + vector_prior;
                // no vector prior under Jones method

                logf << "\t" << log_posterior;

                logf << "\t" << log_likelihood;

                logf << "\t" << log_prior; // starbeast3 does not include coalescent likelihood in this prior


                logf << "\t" << vector_prior;

                logf << "\t" << log_coalescent_likelihood;

                double species_tree_height = p.getSpeciesTreeHeight();
                logf << "\t" << species_tree_height;

                double species_tree_length = p.getSpeciesTreeLength();
                logf << "\t" << species_tree_length;

                vector<double> gene_tree_heights = p.getGeneTreeHeights();
                vector<double> gene_tree_lengths = p.getGeneTreeLengths();
                assert (gene_tree_heights.size() == gene_tree_lengths.size());

                for (int i=0; i<gene_tree_heights.size(); i++) {
                    logf << "\t" << gene_tree_heights[i];
                    logf << "\t" << gene_tree_lengths[i];
                }

                double yule_model = p.getSpeciesTreePrior();
                logf << "\t" << yule_model;

                logf << "\t" << p.getPopMean() / 4.0; // beast uses Ne * u = theta / 4

                for (int i=0; i<(G::_nspecies*2-1); i++) {
    #if defined (DRAW_NEW_THETA)
                    vector<double> theta_vec = p.getThetaVector();
                    
                    if (i < theta_vec.size()) {
                        logf << "\t" << theta_vec[i] / 4.0;
                    }
                    else {
                        logf << "\t" << -1;
                    }
    #else
                    logf << "\t" << G::_theta / 4.0; // all pop sizes are the same under this model, Ne*u = theta / 4?
    #endif
                }

                logf << "\t" << G::_lambda; // TODO: not estimating lambda for now
            
                vector<double> gene_tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
                vector<double> gene_tree_priors = p.getGeneTreePriors();
                assert (gene_tree_log_likelihoods.size() == gene_tree_priors.size());

                for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                    logf << "\t" << gene_tree_log_likelihoods[i];
                }

                for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                    logf << "\t" << gene_tree_priors[i];
                }

                logf << endl;
            }

        logf.close();
    }

    inline void Proj::writeParamsFileForBeastComparisonTestB(vector<Particle> &v, string filename) const {
        // this function creates a params file that is comparable to output from starbeast3
        ofstream logf(filename);
        logf << "iter ";
        logf << "\t" << "posterior ";
        logf << "\t" << "likelihood ";
        logf << "\t" << "prior ";
        logf << "\t" << "vectorPrior "; // log joint prior on population sizes (vectorPrior)
        logf << "\t" << "speciescoalescent ";
        logf << "\t" << "Tree.t:Species.height ";
        logf << "\t" << "Tree.t:Species.treeLength ";

        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "Tree.t:gene" + to_string(i) + "height";
            logf << "\t" << "Tree.t:gene" + to_string(i) + "treeLength";
        }

        logf << "\t" << "YuleModel.t:Species "; // this is the log probability of the species tree (multiply by log(3!) to get increment log prob)
        logf << "\t" << "popMean "; // this is psi in the InverseGamma(2,psi) distribution of popSize

        for (int i=0; i<(G::_nspecies*2-1); i++) {
            logf << "\t" << "popSize." + to_string(i+1);
        }

        logf << "\t" << "speciationRate.t:Species ";

        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "treeLikelihood:gene" + to_string(i);
        }
        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "treePrior:gene" + to_string(i);
        }
        logf << endl;

        int iter = 0;
        
        for (auto &p:v) {
                logf << iter;
                iter++;

                double vector_prior = 0.0; // no vector prior with jones coalescent likelihood

                double log_coalescent_likelihood = 0.0;
                for (unsigned g=1; g<G::_nloci+1; g++) {
                    log_coalescent_likelihood += p.getCoalescentLikelihood(g);
                }

                unsigned locus = p.getNextGene() - 1;
                double log_likelihood = p.calcLogLikelihoodLocus(locus, true);
                double log_prior = p.getAllPriorsFirstRound();

                double log_posterior = log_likelihood + log_prior + log_coalescent_likelihood + vector_prior;
                // no vector prior under Jones method

                logf << "\t" << log_posterior;

                logf << "\t" << log_likelihood;

                logf << "\t" << log_prior; // starbeast3 does not include coalescent likelihood in this prior


                logf << "\t" << vector_prior;

                logf << "\t" << log_coalescent_likelihood;

                double species_tree_height = p.getSpeciesTreeHeight();
                logf << "\t" << species_tree_height;

                double species_tree_length = p.getSpeciesTreeLength();
                logf << "\t" << species_tree_length;

                vector<double> gene_tree_heights = p.getGeneTreeHeights();
                vector<double> gene_tree_lengths = p.getGeneTreeLengths();
                assert (gene_tree_heights.size() == gene_tree_lengths.size());

                for (int i=0; i<gene_tree_heights.size(); i++) {
                    logf << "\t" << gene_tree_heights[i];
                    logf << "\t" << gene_tree_lengths[i];
                }

                double yule_model = p.getSpeciesTreePrior();
                logf << "\t" << yule_model;

                logf << "\t" << p.getPopMean() / 4.0; // beast uses Ne * u = theta / 4

                for (int i=0; i<(G::_nspecies*2-1); i++) {
    #if defined (DRAW_NEW_THETA)
                    vector<double> theta_vec = p.getThetaVector();
                    
                    if (i < theta_vec.size()) {
                        logf << "\t" << theta_vec[i] / 4.0;
                    }
                    else {
                        logf << "\t" << -1;
                    }
    #else
                    logf << "\t" << G::_theta / 4.0; // all pop sizes are the same under this model, Ne*u = theta / 4?
    #endif
                }

                logf << "\t" << G::_lambda; // TODO: not estimating lambda for now
            
                vector<double> gene_tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
                vector<double> gene_tree_priors = p.getGeneTreePriors();
                assert (gene_tree_log_likelihoods.size() == gene_tree_priors.size());

                for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                    logf << "\t" << gene_tree_log_likelihoods[i];
                }

                for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                    logf << "\t" << gene_tree_priors[i];
                }

                logf << endl;
            }

        logf.close();
    }

    inline void Proj::writeParamsFileForBeastComparison(vector<Particle> &v) {
        // this function creates a params file that is comparable to output from starbeast3
        ofstream logf("params-beast-comparison.log");
        logf << "iter ";
        logf << "\t" << "posterior ";
        logf << "\t" << "likelihood ";
        logf << "\t" << "prior ";
        logf << "\t" << "vectorPrior "; // log joint prior on population sizes (vectorPrior)
        logf << "\t" << "speciescoalescent ";
        logf << "\t" << "Tree.t:Species.height ";
        logf << "\t" << "Tree.t:Species.treeLength ";

        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "Tree.t:gene" + to_string(i) + "height";
            logf << "\t" << "Tree.t:gene" + to_string(i) + "treeLength";
        }

        logf << "\t" << "YuleModel.t:Species "; // this is the log probability of the species tree (multiply by log(3!) to get increment log prob)
        logf << "\t" << "popMean "; // this is psi in the InverseGamma(2,psi) distribution of popSize

        for (int i=0; i<(G::_nspecies*2-1); i++) {
            logf << "\t" << "popSize." + to_string(i+1);
        }

        logf << "\t" << "speciationRate.t:Species ";

        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "treeLikelihood:gene" + to_string(i);
        }
        for (int i=1; i<G::_nloci+1; i++) {
            logf << "\t" << "treePrior:gene" + to_string(i);
        }
        logf << endl;
        
        if (_gene_tree_ranks.size() == 0) {
            _gene_tree_ranks.resize(G::_nloci);
        }
        
        if (_bhv_distances_genes.size() == 0) {
            _bhv_distances_genes.resize(G::_nloci);
        }
        
        if (_gene_tree_heights.size() == 0) {
            _gene_tree_heights.resize(G::_nloci);
        }
        
        int iter = 0;
        for (auto &p:v) {
            vector<double> gene_tree_heights;
            gene_tree_heights = p.getGeneTreeHeights();
            
            double species_tree_height_after_first_round = p.getSpeciesTreeHeightAfterFirstRound();
            
//            if (G::_ruv || G::_hpd) {
//                for (unsigned i=0; i<G::_nloci; i++) {
//                    ofstream heightf;
//                    string filename = "gene_tree_heights" + to_string(i+1) + ".txt";
//                    heightf.open(filename, std::ios::app);
//                    heightf << gene_tree_heights[i] << endl;
//                }
//            }
            
            if (G::_ruv) {
                for (unsigned l=0; l<G::_nloci; l++) {
                    _gene_tree_ranks[l].push_back(make_pair(gene_tree_heights[l], false));
                    _gene_tree_heights[l].push_back(gene_tree_heights[l]);
                }
            }
            
            if (G::_ruv_first_level_species) {
                _species_tree_ranks_after_first_round.push_back(make_pair(species_tree_height_after_first_round, false));
                _species_tree_heights_after_first_round.push_back(species_tree_height_after_first_round);
            }
            
//            if (G::_hpd_first_level_species) {
//                ofstream heightf;
//                string filename = "species_tree_heights_after_first_round.txt";
//                heightf.open(filename, std::ios::app);
//                heightf << species_tree_height_after_first_round << endl;
//            }
            
            if (G::_bhv_reference != "" || G::_bhv_reference_path != ".") {
                if (G::_bhv_reference_path != ".") {
                    G::_bhv_reference = readNewickFromFile(G::_bhv_reference_path);
                }
                vector<double> bhv_dists = p.calcBHVDistanceGeneTrees();
                for (unsigned l=0; l<G::_nloci; l++) {
                    _bhv_distances_genes[l].push_back(bhv_dists[l]);
                }
            }
            
            logf << iter;
            iter++;

            double vector_prior = 0.0; // no vector prior with jones coalescent likelihood

            double log_coalescent_likelihood = 0.0;
            for (unsigned g=0; g<G::_nloci; g++) {
                log_coalescent_likelihood += p.getCoalescentLikelihood(g);
            }

            double log_likelihood = p.getLogLikelihood();
            double log_prior = p.getAllPriorsFirstRound();

            double log_posterior = log_likelihood + log_prior + log_coalescent_likelihood + vector_prior;
            
            if (G::_hpd) {
                _hpd_values_genes.resize(G::_nloci);
                for (unsigned l=0; l<G::_nloci; l++) {
                    _hpd_values_genes[l].push_back(make_pair(log_posterior, gene_tree_heights[l]));
                }
            }
            
            if (G::_hpd_first_level_species) {
                _hpd_first_level_values.push_back(make_pair(log_posterior, species_tree_height_after_first_round));
            }
            // no vector prior under Jones method

            logf << "\t" << log_posterior;

            logf << "\t" << log_likelihood;

            logf << "\t" << log_prior; // starbeast3 does not include coalescent likelihood in this prior


            logf << "\t" << vector_prior;

            logf << "\t" << log_coalescent_likelihood;

            double species_tree_height = p.getSpeciesTreeHeight();
            logf << "\t" << species_tree_height;

            double species_tree_length = p.getSpeciesTreeLength();
            logf << "\t" << species_tree_length;

//            vector<double> gene_tree_heights = p.getGeneTreeHeights();
            vector<double> gene_tree_lengths = p.getGeneTreeLengths();
            assert (gene_tree_heights.size() == gene_tree_lengths.size());

            for (int i=0; i<gene_tree_heights.size(); i++) {
                logf << "\t" << gene_tree_heights[i];
                logf << "\t" << gene_tree_lengths[i];
            }

            double yule_model = p.getSpeciesTreePrior();
            logf << "\t" << yule_model;

            logf << "\t" << p.getPopMean() / 4.0; // beast uses Ne * u = theta / 4

            for (int i=0; i<(G::_nspecies*2-1); i++) {
#if defined (DRAW_NEW_THETA)
                vector<double> theta_vec = p.getThetaVector();
                logf << "\t" << theta_vec[i] / 4.0;
#else
                logf << "\t" << G::_theta / 4.0; // all pop sizes are the same under this model, Ne*u = theta / 4?
#endif
            }

            logf << "\t" << G::_lambda; // TODO: not estimating lambda for now
        
            vector<double> gene_tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
            vector<double> gene_tree_priors = p.getGeneTreePriors();
            assert (gene_tree_log_likelihoods.size() == gene_tree_priors.size());

            for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                logf << "\t" << gene_tree_log_likelihoods[i];
            }

            for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                logf << "\t" << gene_tree_priors[i];
            }

            logf << endl;
        }

        logf.close();
    }

    inline void Proj::writeParamsFileForBeastComparisonAfterSpeciesFiltering(vector<Particle> &v, string filename, unsigned group_number) {
        // this function creates a params file that is comparable to output from starbeast3
        std::ofstream logf;

        logf.open(filename, std::ios_base::app);

        if (_first_line) {
            _first_line = false;
            logf << "\t" << "posterior ";
            logf << "\t" << "likelihood ";
            logf << "\t" << "prior ";
            logf << "\t" << "vectorPrior "; // log joint prior on population sizes (vectorPrior)
            logf << "\t" << "speciescoalescent ";
            logf << "\t" << "Tree.t:Species.height ";
            logf << "\t" << "Tree.t:Species.treeLength ";

            for (int i=1; i<G::_nloci+1; i++) {
                logf << "\t" << "Tree.t:gene" + to_string(i) + "height";
                logf << "\t" << "Tree.t:gene" + to_string(i) + "treeLength";
            }

            logf << "\t" << "YuleModel.t:Species "; // this is the log probability of the species tree (multiply by log(3!) to get increment log prob)
            logf << "\t" << "popMean "; // this is psi in the InverseGamma(2,psi) distribution of popSize

            for (int i=0; i<(G::_nspecies*2-1); i++) {
                logf << "\t" << "popSize." + to_string(i+1);
            }

            logf << "\t" << "speciationRate.t:Species ";

            for (int i=1; i<G::_nloci+1; i++) {
                logf << "\t" << "treeLikelihood:gene" + to_string(i);
            }
            for (int i=1; i<G::_nloci+1; i++) {
                logf << "\t" << "treePrior:gene" + to_string(i);
            }
            logf << endl;
        }

        unsigned sample_size = round(double (G::_particle_increase) / double(G::_save_every) );
        if (sample_size == 0) {
            sample_size = G::_particle_increase;
        }
        
        for (unsigned i=0; i<_second_level_indices_to_keep[group_number].size(); i++) {
            Particle p = v[_second_level_indices_to_keep[group_number][i]];
//        for (auto &p:v) {
            double log_coalescent_likelihood = 0.0;
            log_coalescent_likelihood += p.getCoalescentLikelihoodSecondLevel();

            double vector_prior = 0.0;

//            double log_likelihood = p.getLogLikelihood();
            double log_likelihood = 0.0;
            vector<double> gene_tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
            for (auto &g:gene_tree_log_likelihoods) {
                log_likelihood += g;
            }
            double log_prior = p.getAllPriors();

            double log_posterior = log_likelihood + log_prior + log_coalescent_likelihood + vector_prior;

            logf << "\t" << log_posterior;

            logf << "\t" << log_likelihood;

            logf << "\t" << log_prior; // starbeast3 does not include coalescent likelihood in the prior

            logf << "\t" << vector_prior;

            logf << "\t" << log_coalescent_likelihood;

            double species_tree_height = p.getSpeciesTreeHeight();
            logf << "\t" << species_tree_height;

            double species_tree_length = p.getSpeciesTreeLength();
            logf << "\t" << species_tree_length;

            vector<double> gene_tree_heights = p.getGeneTreeHeights();
            vector<double> gene_tree_lengths = p.getGeneTreeLengths();
            assert (gene_tree_heights.size() == gene_tree_lengths.size());
            
            _hpd_values.push_back(make_pair(log_posterior, species_tree_height));

            for (int i=0; i<gene_tree_heights.size(); i++) {
                logf << "\t" << gene_tree_heights[i];
                logf << "\t" << gene_tree_lengths[i];
            }

            double yule_model = p.getSpeciesTreePrior();
            logf << "\t" << yule_model;

            logf << "\t" << p.getPopMean() / 4.0; // beast uses Ne * u = theta / 4

            for (int i=0; i<(G::_nspecies*2-1); i++) {
#if defined (DRAW_NEW_THETA)
                if (!G::_gene_newicks_specified) {
                    vector<double> theta_vec = p.getThetaVector();
                    logf << "\t" << theta_vec[i] / 4.0;
                }
                else {
                    logf << "\t" << G::_theta / 4.0;
                }
#else
                logf << "\t" << G::_theta / 4.0; // all pop sizes are the same under this model, Ne*u = theta / 4?
#endif
            }

            logf << "\t" << G::_lambda; // TODO: not estimating lambda for now

//            vector<double> gene_tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
            vector<double> gene_tree_priors = p.getGeneTreeCoalescentLikelihoods();


            for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                logf << "\t" << gene_tree_log_likelihoods[i];
            }

            for (int i=0; i<gene_tree_log_likelihoods.size(); i++) {
                logf << "\t" << gene_tree_priors[i];
            }
            
            if (G::_calc_bhv_distances_to_true_tree) {
                string fname = "bhvdists.txt";
                ofstream speciestrbhvf(fname, std::ios::app);
                double bhvdist = p.calcBHVDistance();
                speciestrbhvf << bhvdist << endl;
            }

            logf << endl;
            
            ofstream heightf;
            heightf.open("species_tree_heights.txt", std::ios::app);
            heightf << p.getSpeciesTreeHeight() << endl;
        }

        logf.close();
        
//        // write species tree heights to a file
//        string fname = "species_tree_heights.txt";
//        std::ofstream heightf;
//
//        heightf.open(fname, std::ios::app);
//
//        for (auto &p:v) {
//            double species_tree_height = p.getSpeciesTreeHeight();
//            heightf << species_tree_height << "\n";
//        }
    }

    inline void Proj::writeParamsFileForBeastComparisonAfterSpeciesFilteringSpeciesOnly(vector<Particle> &v, string filename, unsigned group_number) {
        // this function creates a params file that is comparable to output from starbeast3
        std::ofstream logf;

        logf.open(filename, std::ios_base::app);

        // no gene tree parameters now
        if (_first_line) {
            _first_line = false;
            logf << "\t" << "posterior ";
            logf << "\t" << "prior ";
            logf << "\t" << "vectorPrior "; // log joint prior on population sizes (vectorPrior)
            logf << "\t" << "speciescoalescent ";
            logf << "\t" << "Tree.t:Species.height ";
            logf << "\t" << "Tree.t:Species.treeLength ";

            logf << "\t" << "YuleModel.t:Species "; // this is the log probability of the species tree (multiply by log(3!) to get increment log prob)
            logf << "\t" << "popMean "; // here, this will be user specified theta

            logf << "\t" << "speciationRate.t:Species ";

            logf << endl;
        }

        unsigned sample_size = round(double (G::_particle_increase) / double(G::_save_every) );
        if (sample_size == 0) {
            sample_size = G::_particle_increase;
        }
        
//        for (auto &p:v) {
        for (unsigned i=0; i<_second_level_indices_to_keep[group_number].size(); i++) {
            Particle p = v[_second_level_indices_to_keep[group_number][i]];

            double log_coalescent_likelihood = 0.0;
            log_coalescent_likelihood += p.getCoalescentLikelihood(1);

            double vector_prior = 0.0;

            double log_prior = p.getAllPriors();

            double log_posterior = log_prior + log_coalescent_likelihood + vector_prior;

            logf << "\t" << log_posterior;
            
            logf << "\t" << log_prior;

            logf << "\t" << vector_prior;

            logf << "\t" << log_coalescent_likelihood;

            double species_tree_height = p.getSpeciesTreeHeight();
            logf << "\t" << species_tree_height;
            
            _hpd_values.push_back(make_pair(log_posterior, species_tree_height));

            double species_tree_length = p.getSpeciesTreeLength();
            logf << "\t" << species_tree_length;

            double yule_model = p.getSpeciesTreePrior();
            logf << "\t" << yule_model;

            logf << "\t" << p.getPopMean() / 4.0; // beast uses Ne * u = theta / 4

            logf << "\t" << G::_lambda; // TODO: for now, not using estimate lambda option

            logf << endl;
            
            ofstream heightf;
            heightf.open("species_tree_heights.txt", std::ios::app);
            heightf << p.getSpeciesTreeHeight() << endl;
        }

        logf.close();
    }


    inline void Proj::writeDeepCoalescenceFile(vector<Particle> &v) {
        ofstream logf("deep_coalescences.txt");
        logf << "num deep coalescences = " << v[0].getNumDeepCoalescences() << endl;
        logf << "Maximum number of deep coalescences = " << v[0].getMaxDeepCoalescences() << endl;
        logf << "True species tree height = " << v[0].getSpeciesTreeHeight() << endl;
        
        ofstream spptreeheight("true_species_tree_height.txt");
        spptreeheight << v[0].getSpeciesTreeHeight() << endl;
        
        vector<double> gene_tree_heights = v[0].getGeneTreeHeights();
        for (unsigned l=0; l<G::_nloci; l++) {
            ofstream spptreeheight("true_gene_tree_height" + to_string(l+1) + ".txt");
            spptreeheight << gene_tree_heights[l] << endl;
        }
        
        // Expected height of species tree is (1/2 + 1/3 + ... + 1/nspecies)/lambda
        // Approximation to (1/2 + 1/3 + ... + 1/nspecies) is
        //   ln(nspecies) + 0.58 (Euler's constant)
        //   see http://www.dimostriamogoldbach.it/en/inverses-integers-sum/
        double expected_species_tree_height = 0.0;
        for (unsigned i = 2; i <= G::_sim_nspecies; i++) {
            expected_species_tree_height += 1.0/i;
        }
        expected_species_tree_height /= G::_lambda;

        logf << "Expected species tree height = " << expected_species_tree_height << endl;
        logf << "\n";
    }

    inline void Proj::writeThetaFile(vector<Particle> &v) {
        ofstream logf("simulated-thetas.txt");
        for (auto &p:v) {
            double theta_mean = p.getThetaMean();
            logf << "theta mean: " << theta_mean << "\n";
            vector<double> thetas = p.getThetaMap();
            unsigned count = 1;
            for (auto &t:thetas) {
                logf << "pop " << count << "\t" << t << "\n";
                count++;
            }
        }
        
        logf << "\n";
        logf << "\n";
        
        for (auto &p:v) {
            vector<double> thetas = p.getThetaMap();
            unsigned count = 1;
            for (auto &t:thetas) {
                logf << "pop / 4 " << count << "\t" << t/4 << "\n";
                count++;
            }
        }
    }

    inline void Proj::writeLoradFile(vector<Particle> &v) const {
        ofstream logf("params.log");
        logf << "iteration ";
        logf << "\t" << "likelihood ";
        for (int s=0; s<G::_nspecies-1; s++) {
            logf << "\t" << "species_increment";
        }
        logf << "\t" << "species_tree_prior";
        for (int g=1; g<G::_nloci+1; g++) {
            for (int i=1; i<G::_ntaxa; i++) {
                logf << "\t" << "gene_increment";
            }
        }
        logf << "\t" << "coalescent_likelihood";
        logf << endl;

        unsigned iter = 0;
        for (auto &p:v) {
            logf << iter;
            iter++;

            logf << "\t" << p.getLogLikelihood();

            double species_tree_prior = 0.0;

            for (unsigned g=0; g<G::_nloci+1; g++) {
                for (auto &b:p.getIncrementPriors(g)) {
                    logf << "\t" << b.first;
                    if (g == 0) { // species tree prior
                        species_tree_prior += b.second;
                    }
                    // no increment should be 0
                    assert (b.first > 0.0);
                }
                assert (species_tree_prior != 0.0);

                if (g == 0) {
                    logf << "\t" << species_tree_prior;
                }
            }

            double log_coalescent_likelihood = 0.0;
            for (unsigned g=1; g<G::_nloci+1; g++) {
                log_coalescent_likelihood += p.getCoalescentLikelihood(g);
            }
            logf << "\t" << log_coalescent_likelihood;

            logf << endl;
        }

        logf.close();
    }

    inline void Proj::writeLoradFileAfterSpeciesFiltering(vector<Particle> &v) const {
        ofstream logf("params.log");
        logf << "iteration ";
        for (int s=0; s<G::_nspecies-1; s++) {
            logf << "\t" << "species_increment";
        }
        logf << "\t" << "species_tree_prior";
        logf << "\t" << "coalescent_likelihood";
        logf << endl;

        unsigned iter = 0;
        for (auto &p:v) {
            logf << iter;
            iter++;
            
            for (auto &b:p.getIncrementPriors(0)) {
                logf << "\t" << b.first;
            }
            
            double species_tree_prior = p.getAllPriors(); // TODO: don't need this prior? included in coalescent likelihood?
            assert (species_tree_prior != 0.0);
            logf << "\t" << species_tree_prior;

            double log_coalescent_likelihood = p.getCoalescentLikelihood(0);
            logf << "\t" << log_coalescent_likelihood;

            logf << endl;
        }

        logf.close();
    }

    inline void Proj::saveSpeciesTreesAltHierarchical(vector<Particle> &v, unsigned group_number)  {
        string filename1 = "alt_species_trees.trees";

        assert (G::_start_mode_type != G::StartModeType::START_MODE_SIM);

        unsigned count = 0;
        // save all species trees
        std::ofstream treef;

        treef.open(filename1, std::ios_base::app);
        for (unsigned i=0; i<_second_level_indices_to_keep[group_number].size(); i++) {
            Particle p = v[_second_level_indices_to_keep[group_number][i]];
            treef << "  tree test = [&R] " << p.saveForestNewickAlt()  << ";\n";
            count++;
        }
        treef.close();
    }

    inline void Proj::saveSpeciesTreesHierarchical(vector<Particle> &v, string filename1, string filename2, unsigned group_number) {
        // save only unique species trees
        if (!G::_run_on_empty) {
            vector<vector<pair<double, double>>> unique_increments_and_priors;

            std::ofstream unique_treef;

            unique_treef.open(filename2, std::ios_base::app);
            
//            std::ofstream bhv_logf;
//            if (G::_hpd) {
//                bhv_logf.open("bhv_log.txt", std::ios_base::app);
//            }

            for (unsigned i=0; i<_second_level_indices_to_keep[group_number].size(); i++) {
                Particle p = v[_second_level_indices_to_keep[group_number][i]];
                vector<pair<double, double>> increments_and_priors = p.getSpeciesTreeIncrementPriors();
                
                if (G::_ruv) {
                    double species_tree_height = p.getSpeciesTreeHeight();
                    _species_tree_heights.push_back(species_tree_height);
                    pair<double, bool> rank = make_pair(species_tree_height, false);
                    _ranks.push_back(rank);
                }
                
                if (G::_bhv_reference != "" || G::_bhv_reference_path != ".") {
                    if (G::_bhv_reference_path != ".") {
                        G::_bhv_reference = readNewickFromFile(G::_bhv_reference_path);
                    }
                    _bhv_distances.push_back(p.calcBHVDistance());
                }
                
//                if (G::_hpd) {
//                    double freq = 1;
//                    bhv_logf << freq << "\t" << p.getCoalescentLikelihood(1) << "\t" << p.getAllPriors() << "\t" << p.saveForestNewick() << endl;
//                }
                
            if (G::_write_species_tree_file) {
                bool found = false;
                if(std::find(unique_increments_and_priors.begin(), unique_increments_and_priors.end(), increments_and_priors) != unique_increments_and_priors.end()) {
                    found = true;
                }
                if (!found) {
                    unique_increments_and_priors.push_back(increments_and_priors);
                    unique_treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
                }
            }
            unique_treef.close();
        }
//            if (G::_hpd) {
//                bhv_logf.close();
//            }
    }

        assert (G::_start_mode_type != G::StartModeType::START_MODE_SIM);

        if (G::_write_species_tree_file) {
            unsigned count = 0;
                // save all species trees
                std::ofstream treef;

                treef.open(filename1, std::ios_base::app);
            for (unsigned i=0; i<_second_level_indices_to_keep[group_number].size(); i++) {
                Particle p = v[_second_level_indices_to_keep[group_number][i]];
                    treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
                    count++;
                }
                treef.close();
        }
    }

    inline void Proj::saveSpeciesTreesAfterFirstRound(vector<Particle> &v) const {
        if (G::_sample_from_prior) {
            if (G::_write_species_tree_file) {
                ofstream treef("species_trees.trees");
                treef << "#nexus\n\n";
                treef << "begin trees;\n";
                for (auto &p:v) {
                    treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
                }
                treef << "end;\n";
                treef.close();
            }
        }
        else {
            if (G::_write_species_tree_file) {
                ofstream treef("species_trees_after_first_round.trees");
                treef << "#nexus\n\n";
                treef << "begin trees;\n";
                for (auto &p:v) {
                    treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
                }
                treef << "end;\n";
                treef.close();
            }
        }
    }

    inline void Proj::saveAllSpeciesTrees(vector<Particle> &v) const {
        
        // save all species trees
        ofstream treef("species_trees.trees");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:v) {
            treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
        }
        treef << "end;\n";
        treef.close();
    }

    inline void Proj::saveSpeciesTrees(vector<Particle> &v) const {
        // save only unique species trees
        if (!G::_run_on_empty) {
            vector<vector<pair<double, double>>> unique_increments_and_priors;

            ofstream unique_treef("unique_species_trees.trees");
            unique_treef << "#nexus\n\n";
            unique_treef << "begin trees;\n";
            for (auto &p:v) {
                vector<pair<double, double>> increments_and_priors = p.getSpeciesTreeIncrementPriors();
                bool found = false;
                if(std::find(unique_increments_and_priors.begin(), unique_increments_and_priors.end(), increments_and_priors) != unique_increments_and_priors.end()) {
                    found = true;
                }
                if (!found) {
                    unique_increments_and_priors.push_back(increments_and_priors);
                    unique_treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
                }
            }
            unique_treef << "end;\n";
            unique_treef.close();
        }

        if (G::_start_mode_type == G::StartModeType::START_MODE_SMC) {
            // save all species trees
            ofstream treef("species_trees.trees");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:v) {
                treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
            }
            treef << "end;\n";
            treef.close();
        }
        else {
            // save true species tree
            ofstream treef("true-species-tree.tre");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:v) {
                treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
            }
            treef << "end;\n";
            treef.close();
        }

        }

    inline void Proj::saveGeneTrees(vector<Particle> &v) const {
        if (G::_start_mode_type == G::StartModeType::START_MODE_SMC) {
            for (unsigned i=1; i<G::_nloci+1; i++) {
                string fname = "gene" + to_string(i) + ".trees";
                ofstream treef(fname);
                treef << "#nexus\n\n";
                treef << "begin trees;\n";
                for (auto &p:v) {
                    treef << "tree gene" << i << " = [&R] " << p.saveGeneNewick(i)  << ";\n";
                }
                treef << "end;\n";
            }
        }

        else {
            if (G::_save_gene_trees_separately) {
                for (auto &p:v) {
                    for (unsigned i=1; i<G::_nloci+1; i++) {
                        string fname = "gene" + to_string(i) + ".trees";
                        ofstream treef(fname);
                        treef << "#nexus\n\n";
                        treef << "begin trees;\n";
                        treef << "  tree gene" << i << " = [&R] " << p.saveGeneNewick(i)  << ";\n";
                        treef << endl;
                        treef << "end;\n";
                        treef.close();
                    }
                }
            }
            else {
                // save true species tree
                ofstream treef("true-gene-trees.tre");
                treef << "#nexus\n\n";
                treef << "begin trees;\n";
                for (auto &p:v) {
                        for (unsigned i=1; i<G::_nloci+1; i++) {
                            treef << "tree gene" << i << " = [&R] " << p.saveGeneNewick(i)  << ";\n";
                    }
                    treef << endl;
                }
                treef << "end;\n";
                treef.close();
            }
        }
    }


    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        std::vector<std::string> partition_subsets;
        boost::program_options::variables_map vm;
        boost::program_options::options_description desc("Allowed options");
        
#if defined(SPECIES_IN_CONF)
        vector<string> species_definitions;
#endif

        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile,d",  boost::program_options::value(&G::_data_file_name), "name of a data file in NEXUS format")
        ("subset",  boost::program_options::value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("gpu",           boost::program_options::value(&G::_use_gpu)->default_value(true), "use GPU if available")
        ("ambigmissing",  boost::program_options::value(&G::_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("nparticles",  boost::program_options::value(&G::_nparticles)->default_value(1000), "number of particles")
        ("seed,z", boost::program_options::value(&_random_seed)->default_value(1), "random seed")
        ("theta, t", boost::program_options::value(&G::_theta)->default_value(0.0), "theta")
        ("lambda", boost::program_options::value(&G::_lambda)->default_value(0.0), "speciation rate")
        ("model", boost::program_options::value(&G::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&G::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("nthreads",  boost::program_options::value(&G::_nthreads)->default_value(1), "number of threads for multi threading")
        ("run_on_empty", boost::program_options::value(&G::_run_on_empty)->default_value(false), "run program without data")
        ("run_on_empty_first_level_only", boost::program_options::value(&G::_run_on_empty_first_level_only)->default_value(false), "run program without data only for the first level; if setting this to true, also set run_on_empty = true")
        ("verbose", boost::program_options::value(&G::_verbose)->default_value(1), "set amount of output printed")
        ("save_memory", boost::program_options::value(&G::_save_memory)->default_value(false), "save memory at the expense of time")
        ("outgroup", boost::program_options::value(&G::_outgroup)->default_value("none"), "a string defining the outgroup")
        ("startmode", boost::program_options::value(&G::_start_mode)->default_value("smc"), "a string defining whether to simulate data or perform smc")
        ("nspecies", boost::program_options::value(&G::_sim_nspecies)->default_value(0), "number of species to simulate")
        ("ntaxaperspecies", boost::program_options::value(&G::_string_ntaxaperspecies)->default_value(""), "number of taxa per species to simulate")
        ("filename", boost::program_options::value(&G::_sim_file_name), "name of file to write simulated data to")
        ("particle_increase", boost::program_options::value(&G::_particle_increase)->default_value(1), "how much to increase particles for species filtering")
        ("thin", boost::program_options::value(&G::_thin)->default_value(1.0), "take this portion of particles for hierarchical species filtering")
        ("save_every", boost::program_options::value(&G::_save_every)->default_value(1.0), "take this portion of particles for output")
        ("save_gene_trees", boost::program_options::value(&G::_save_gene_trees)->default_value(true), "turn this off to not save gene trees and speed up program")
        ("gene_newicks", boost::program_options::value(&G::_gene_newicks_specified)->default_value(false), "set true if user is specifying gene tree files")
        ("ngenes", boost::program_options::value(&G::_ngenes_provided)->default_value(0), "number of gene newick files specified")
        ("theta_proposal_mean", boost::program_options::value(&G::_theta_proposal_mean)->default_value(0.0), "theta proposal mean")
        ("theta_prior_mean", boost::program_options::value(&G::_theta_prior_mean)->default_value(0.0), "theta prior mean")
        ("species_newick", boost::program_options::value(&G::_species_newick_name)->default_value("null"), "name of file containing species newick descriptions")
        ("fix_theta_for_simulations",  boost::program_options::value(&G::_fix_theta_for_simulations)->default_value(true), "set to true to fix one theta for all populations")
        ("fix_theta",  boost::program_options::value(&G::_fix_theta)->default_value(false), "set to true to fix one theta for all populations")
        ("relative_rates", boost::program_options::value(&G::_string_relative_rates)->default_value("null"))
        ("simoccupancy", boost::program_options::value(&Data::_occupancy)->default_value(1.0), "probability that any given taxon will have data for any given locus; 1-_occupancy is prob. all missing data for a taxon (used only if startmode is) 'sim'")
        ("simedgeratevar",   boost::program_options::value(&Forest::_edge_rate_variance)->default_value(0.0), "variance of lognormal relative rate distribution across edges in gene trees (used only if startmode is 'sim')")
        ("simasrvshape",   boost::program_options::value(&Forest::_asrv_shape)->default_value(G::_infinity), "Shape of gamma among-site rate heterogeneity within a locus (used only if startmode is 'sim')")
        ("simcomphet",   boost::program_options::value(&Forest::_comphet)->default_value(G::_infinity), "Dirichlet parameter governing compositional heterogeneity (default value results in compositional homogeneity (used only if startmode is 'sim')")
        ("save_gene_trees_separately", boost::program_options::value(&G::_save_gene_trees_separately)->default_value(false), "for simulations, save gene trees in separate files")
        ("newick_path", boost::program_options::value(&G::_newick_path)->default_value(""), "path to gene newicks are if starting from gene newicks and only performing SMC on second round")
        ("ngroups", boost::program_options::value(&G::_ngroups)->default_value(1), "number of populations")
        ("mcmc", boost::program_options::value(&G::_mcmc)->default_value(false), "use mcmc moves in analysis")
        ("sliding_window", boost::program_options::value(&G::_sliding_window)->default_value(0.05), "size of sliding window to use in mcmc analysis")
        ("n_mcmc_rounds", boost::program_options::value(&G::_n_mcmc_rounds)->default_value(1), "number of rounds to use for mcmc analysis")
        ("ruv", boost::program_options::value(&G::_ruv)->default_value(false), "calculate ranks for ruv")
        ("hpd", boost::program_options::value(&G::_hpd)->default_value(false), "calculate hpd intervals for coverage")
        ("ruv_first_level_species", boost::program_options::value(&G::_ruv_first_level_species)->default_value(false), "calculate species tree ranks for ruv after first round") // TODO: also turn on second_level = false
        ("hpd_first_level_species", boost::program_options::value(&G::_hpd_first_level_species)->default_value(false), "calculate hpd intervals for coverage after first level species filtering")
        ("bhv_reference", boost::program_options::value(&G::_bhv_reference)->default_value(""), "newick string to use for BHV distance reference")
        ("bhv_reference_path", boost::program_options::value(&G::_bhv_reference_path)->default_value("."), "path to use for BHV distance reference; can specify either bhv_reference or bhv_reference_path; bhv_reference string will take priority")
        ("write_species_tree_file", boost::program_options::value(&G::_write_species_tree_file)->default_value(true), "set to false to not write species tree newicks to a file - only use this option to turn on for RUV calculations when lots of trees will be saved")
        ("second_level", boost::program_options::value(&G::_second_level)->default_value(true), "set to false to not run second level")
        ("calc_bhv_distances_to_true_tree", boost::program_options::value(&G::_calc_bhv_distances_to_true_tree)->default_value(false), "set to true to calculate bhv distances between every sampled tree and the true tree")
        ("sample_from_prior", boost::program_options::value(&G::_sample_from_prior)->default_value(false), "sample species trees from prior")
        ("nloci_slow_rate", boost::program_options::value(&G::_nloci_slow_rate)->default_value(0), "for simulations - number of loci to simulate at slower rate")
        ("plus_G", boost::program_options::value(&G::_plus_G)->default_value(false), "+G rate het")
        ("plus_I", boost::program_options::value(&G::_plus_I)->default_value(false), "+I rate het")
        ("pinvar", boost::program_options::value(&G::_pinvar)->default_value(0.0), "pinvar")
        ("gamma_rate_var", boost::program_options::value(&G::_gamma_rate_var)->default_value(1000), "alpha value for +G rate het")
#if defined(SPECIES_IN_CONF)
        ("species", boost::program_options::value(&species_definitions), "a string defining a species, e.g. 'A:x,y,z' says that taxa x, y, and z are in species A")
#endif
        ;

        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        try {
            const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >("proj.conf", desc, false);
            boost::program_options::store(parsed, vm);
        }
        catch(boost::program_options::reading_file & x) {
            std::cout << "Note: configuration file (proj.conf) not found" << std::endl;
        }
        boost::program_options::notify(vm);

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            std::cout << desc << "\n";
            std::exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            std::cout << boost::str(boost::format("This is %s version %d.%d") % G::_program_name % G::_major_version % G::_minor_version) << std::endl;
            std::exit(1);
        }

        // If user specified --subset on command line, break specified partition subset
        // definition into name and character set string and add to _partition
        if (vm.count("subset") > 0) {
            _partition.reset(new Partition());
            for (auto  s : partition_subsets) {
                _partition->parseSubsetDefinition(s);
            }
        }
        
#if defined(SPECIES_IN_CONF)
        if (vm.count("species") > 0) {
            for (auto s : species_definitions) {
                parseSpeciesDefinition(s);
            }
        }
#endif

        // If user specified "base_frequencies" in conf file, convert them to a vector<double>
        if (vm.count("base_frequencies") > 0) {
            handleBaseFrequencies();
        }
        
        if (vm.count("relative_rates") > 0) {
            handleRelativeRates();
        }
            
        // if user specified "ntaxaperspecies" in conf file, convert them to a vector<unsigned>
        if (vm.count("ntaxaperspecies") > 0 && G::_start_mode == "sim") {
            handleNTaxaPerSpecies();
        }
        
        // if save_every > particle_increase, quit
        if (G::_save_every > G::_particle_increase) {
            throw XProj("particle_increase must be greater than or equal to save_every");
        }
        
        if (G::_model == "JC") {
            cout << "Setting kappa to 1.0 under JC model\n";
            cout << "Setting base frequencies equal under JC model\n";
            if (Forest::_kappa != 1.0) {
                cout << "\nIgnoring kappa under JC model\n";
            }
        }
        if (G::_start_mode == "sim") {
            if (G::_data_file_name != "") {
                cout << "\nIgnoring data file name for simulation\n";
            }
            if (G::_sim_nspecies == 0) {
                throw XProj("must specify number of species for which to simulate data");
            }
            
            if (G::_ntaxaperspecies.size() != 1 && G::_ntaxaperspecies.size() != G::_sim_nspecies) {
                throw XProj("ntaxaperspecies must be one number for all species or must match the total number of species specified; ex: ntaxaperspecies = 5 or ntaxaperspecies = 5, 2, 3 if nspecies = 3");
            }
            if (G::_run_on_empty) {
                cout << "\nIgnoring start_mode = run_on_empty and simulating data\n";
            }
            if (G::_sim_file_name == "") {
                throw XProj("must specify name of file to write simulated data to; ex. filename = sim.nex");
            }
            if (G::_theta == 0.0 && G::_theta_prior_mean == 0.0 && G::_theta_proposal_mean == 0.0) {
                throw XProj("must specify theta or theta proposal / prior mean for simulations");
            }
        }
        else {
            if (G::_data_file_name == "") {
                throw XProj("must specify name of data file if smc option is chosen; ex. data file = sim.nex");
            }
            if (G::_theta_prior_mean == 0.0 && G::_theta_proposal_mean > 0.0) {
                cout << boost::format("\nSetting theta prior mean equal to theta proposal mean of %d\n") % G::_theta_proposal_mean;
                G::_theta_prior_mean = G::_theta_proposal_mean;
            }
            // no proposal or prior mean if theta fixed
            else if (G::_theta_prior_mean > 0.0 && G::_theta_proposal_mean ==  0.0) {
                cout << boost::format("\nSetting theta proposal mean equal to theta prior mean of %d\n") % G::_theta_prior_mean;
                G::_theta_proposal_mean = G::_theta_prior_mean;
            }
            else if (G::_theta_prior_mean == 0.0 && G::_theta_proposal_mean == 0.0) {
                cout << boost::format("\nTheta mean of %d will be fixed for all particles; population sizes will all be drawn from the same theta\n") % G::_theta;
            }
        }
        
        if (G::_model == "JC") {
            G::_model_type = G::ModelType::MODEL_TYPE_JC;
        }
        else if (G::_model == "HKY") {
           G:: _model_type = G::ModelType::MODEL_TYPE_HKY;
        }
        else {
            throw XProj("must specify either JC or HKY for model type");
        }
        
        if (G::_start_mode == "smc") {
            G::_start_mode_type = G::StartModeType::START_MODE_SMC;
        }
        else if (G::_start_mode == "sim") {
            G::_start_mode_type = G::StartModeType::START_MODE_SIM;
        }
        else {
            throw XProj("must specify either JC or HKY for model type");
        }
        
        if (G::_species_newick_name != "null") {
            G::_species_newick_specified = true;
        }
        
        if (G::_proposal == "prior-prior") {
            G::_prior_prior = true;
        }
    }

    inline void Proj::checkOutgroupName() {
        bool found = false;
        for (auto &s:_taxon_map) {
            if (G::_outgroup == s.second) {
                found = true;
            }
        }
        if (!found) {
            throw XProj(format("outgroup name does not match any species name"));
        }
    }

    inline void Proj::handleNTaxaPerSpecies() {
        vector<string> temp;
        split(temp, G::_string_ntaxaperspecies, is_any_of(","));
        // iterate throgh temp
        if (temp[0] == "") {
            throw XProj("must specify number of taxa per species");
        }
        for (auto &i:temp) {
            double f = stof(i);
            G::_ntaxaperspecies.push_back(f);
        }
    }

    inline void Proj::handleBaseFrequencies() {
        vector <string> temp;
        split(temp, G::_string_base_frequencies, is_any_of(","));
        double sum = 0.0;
        // iterate through temp
        for (auto &i:temp) {
            double f = stof(i);
            G::_base_frequencies.push_back(f);
            sum +=f;
        }
        if (fabs(sum-1)>0.000001) {
            throw XProj(format("base frequencies (%s) don't add to 1")%G::_string_base_frequencies);
        }
        assert (fabs(sum-1) < 0.000001);
    }

    inline void Proj::handleRelativeRates() {
        vector <string> temp;
        assert (G::_double_relative_rates.size() == 0);
        split(temp, G::_string_relative_rates, is_any_of(","));
        if (G::_string_relative_rates == "null") {
            for (unsigned i=0; i<_partition->getNumSubsets(); i++) {
                G::_double_relative_rates.push_back(1.0); // if no relative rates provided, set them all to 1
            }
        }
        else {
        // iterate through temp
            for (auto &i:temp) {
                double f = stof(i);
                G::_double_relative_rates.push_back(f);
            }
        }
        
        if (G::_start_mode == "sim" && G::_nloci_slow_rate != 0) {
            // take first n loci and decrease relative rates
#if defined (INFO_TEST)
            // TODO: assuming 100 loci total
            for (unsigned n=0; n<100; n++) {
                if (n < 10) {
//                    G::_double_relative_rates[n] *= 0.0000842; // 0.001 unnormalized rate
                    G::_double_relative_rates[n] *= 0.0001;
                }
                else if (n < 20) {
//                    G::_double_relative_rates[n] *= 0.000842; // 0.01 rate
                    G::_double_relative_rates[n] *= 0.001;
                }
                else if (n < 30) {
//                    G::_double_relative_rates[n] *= 0.00842; // 0.1
                    G::_double_relative_rates[n] *= 0.01;
                }
                else if (n < 40) {
//                    G::_double_relative_rates[n] *= 0.0421; // 0.5
                    G::_double_relative_rates[n] *= 0.1;
                }
                else if (n < 60) {
//                    G::_double_relative_rates[n] *= 0.0842; // 1.0
                    G::_double_relative_rates[n] *= 1;
                }
                else if (n < 70) {
//                    G::_double_relative_rates[n] *= 0.09262; // 1.1
                    G::_double_relative_rates[n] *= 1.1;
                }
                else if (n < 80) {
//                    G::_double_relative_rates[n] *= 0.421; // 5.0
                    G::_double_relative_rates[n] *= 1.5;
                }
                else if (n < 90) {
//                    G::_double_relative_rates[n] *= 0.842; // 10.0
                    G::_double_relative_rates[n] *= 2.0;
                }
                else {
//                    G::_double_relative_rates[n] *= 8.42; // 100.0
                    G::_double_relative_rates[n] *= 3.3;
                }
            }
#else
            for (unsigned n=0; n < G::_nloci_slow_rate; n++) {
                G::_double_relative_rates[n] *= 0.01;
            }
#endif
        }
    }

    inline void Proj::summarizeData(Data::SharedPtr) {
        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();

        std::cout << "\nNumber of taxa: " << _data->getNumTaxa() << std::endl;
        std::cout << "Number of partition subsets: " << nsubsets << std::endl;
        std::cout << "Number of particles: " << G::_nparticles << std::endl;

        for (unsigned subset = 0; subset < nsubsets; subset++) {
            DataType dt = _partition->getDataTypeForSubset(subset);
            std::cout << "  Subset " << (subset+1) << " (" << _data->getSubsetName(subset) << ")" << std::endl;
            std::cout << "    data type: " << dt.getDataTypeAsString() << std::endl;
            std::cout << "    sites:     " << _data->calcSeqLenInSubset(subset) << std::endl;
            std::cout << "    patterns:  " << _data->getNumPatternsInSubset(subset) << std::endl;
        }
    }

    inline string Proj::readNewickFromFile(string file_name) {
        ifstream infile(file_name);
        string newick_string;
        string newick;
            
        int size_before = (int) newick.size();
        while (getline(infile, newick)) { // file newicks must start with the word "tree" two spaces from the left margin
            if (newick.find("tree") == 2) { // TODO: this only works if there are exactly 2 spaces before the newick - try starting at parenthesis
                size_t pos = newick.find("("); //find location of parenthesis
                newick.erase(0,pos); //delete everything prior to location found
                newick_string = newick;
                }
            }
            int size_after = (int) newick_string.size();
            if (size_before == size_after) {
                throw XProj("cannot find species newick file");
            }
        
        return newick_string;
    }

    inline string Proj::handleSpeciesNewick() {
        ifstream infile(G::_species_newick_name);
        string newick_string;
        string newick;
        int size_before = (int) newick.size();
        while (getline(infile, newick)) { // file newicks must start with the word "tree" two spaces from the left margin
            if (newick.find("tree") == 2) { // TODO: this only works if there are exactly 2 spaces before the newick - try starting at parenthesis
                size_t pos = newick.find("("); //find location of parenthesis
                newick.erase(0,pos); //delete everything prior to location found
                newick_string = newick;
                }
            }
            int size_after = (int) newick_string.size();
            if (size_before == size_after) {
                throw XProj("cannot find species newick file");
            }
        
        return newick_string;
    }

    inline void Proj::handleGeneNewicks() {
        G::_in_second_level = true;
        vector<vector<string>> newicks; // vector of vector of newicks, 1 vector per gene
        _first_line = true;
        if (G::_ngenes_provided == 0) {
            throw XProj("must specify number of genes in the conf file");
        }
        
        for (int i=1; i<G::_ngenes_provided+1; i++) {
            vector<string> current_gene_newicks;

            string file_name = "";
            if (G::_newick_path != ".") {
                file_name = G::_newick_path + "/" + "gene" + to_string(i) + ".trees"; // file must be named gene1.trees, gene2.trees, etc.
            }
            else {
                file_name = "gene" + to_string(i) + ".trees";
            }
            ifstream infile (file_name);
            string newick;
            unsigned size_before = (unsigned) current_gene_newicks.size();
            
            while (getline(infile, newick)) { // file newicks must start with the word "tree"
                if (current_gene_newicks.size() < G::_nparticles * G::_ngroups) { // stop adding newicks once the number of particles has been reached // TODO: add option to randomize this?
                if (newick.find("tree") == 2) { // TODO: this only works if there are exactly 2 spaces before the newick - try starting at parenthesis
                        size_t pos = newick.find("("); //find location of parenthesis
                        newick.erase(0,pos); //delete everything prior to location found
                        current_gene_newicks.push_back(newick);
                    }
                    else if (newick.find("tree") == 0) { // TODO: this only works if there are exactly 0 spaces before the newick - try starting at parenthesis
                            size_t pos = newick.find("("); //find location of parenthesis
                            newick.erase(0,pos); //delete everything prior to location found
                            current_gene_newicks.push_back(newick);
                        }
                }
            }
            int size_after = (int) current_gene_newicks.size();
            if (size_before == size_after) {
                string error_message = "cannot find gene " + to_string(i) + " file";
                throw XProj(error_message);
            }
            
            newicks.push_back(current_gene_newicks);
        }

        _data = Data::SharedPtr(new Data()); // don't need to set data b/c will not calculate a Felsenstein likelihood
        _data->setPartition(_partition);
        _data->getDataFromFile(G::_data_file_name);

        if (G::_verbose > 0) {
            summarizeData(_data);
        }
        createSpeciesMap(true); // use data
        G::_nspecies = (unsigned) G::_species_names.size();
        assert (G::_nspecies > 0);
        
        createSpeciesMapTyped(); // use data

        // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
        if (G::_outgroup != "none") {
            checkOutgroupName();
        }

        //set number of species to number in data file
        G::_ntaxa = _data->getNumTaxa();
        assert (G::_species_names.size() > 0);
        assert (G::_species_names_typed.size() > 0);
        G::_nloci = _data->getNumSubsets();
        rng.setSeed(_random_seed);
        rng_mcmc.setSeed(_random_seed + 2);
        
        buildSpeciesMap(/*taxa_from_data*/true);

        Particle p;
        initializeParticle(p); // initialize one particle and copy to all other particles
        
        vector<Particle> my_vec;
        my_vec.resize(G::_nparticles * G::_ngroups, p);
        
        // reset pointers so particles have their own pointers
        for (unsigned p=0; p<G::_nparticles * G::_ngroups; p++) {
            vector<Forest::SharedPtr> gfcpies;
            for (unsigned l=0; l<G::_nloci; l++) {
                gfcpies.push_back(Forest::SharedPtr(new Forest()));
            }
            my_vec[p].resetSubgroupPointers(gfcpies);
        }

        unsigned psuffix = 1;
        for (auto &p:my_vec) {
            p.setSeed(rng.randint(1,9999) + psuffix);
            psuffix += 2;
        }
        
        unsigned count = 0;
        for (auto &p:my_vec) {
            vector<string> particle_newicks;
            for (int i = 0; i<newicks.size(); i++) {
                // if number of particles > size of newicks, just use the same newick for each particle
                if (newicks[i].size() <= count) {
                    count = 0;
                }
                particle_newicks.push_back(newicks[i][count]);
            }
            assert (particle_newicks.size() == G::_nloci);
            p.processGeneNewicks(particle_newicks);
            count++;
        }

        cout << "\n";
        string filename1 = "species_trees.trees";
        string filename2 = "unique_species_trees.trees";
        string filename3 = "params-beast-comparison.log";
        if (filesystem::remove(filename1)) {
            ofstream speciestrf(filename1);
            speciestrf << "#nexus\n\n";
            speciestrf << "begin trees;\n";
            if (G::_verbose > 0) {
                cout << "existing file " << filename1 << " removed and replaced\n";
            }
        }
        else {
            ofstream speciestrf(filename1);
            speciestrf << "#nexus\n\n";
            speciestrf << "begin trees;\n";
            if (G::_verbose > 0) {
                cout << "created new file " << filename1 << "\n";
            }
        }
        if (filesystem::remove(filename2)) {
            ofstream uniquespeciestrf(filename2);
            uniquespeciestrf << "#nexus\n\n";
            uniquespeciestrf << "begin trees;\n";
            if (G::_verbose > 0) {
                cout << "existing file " << filename2 << " removed and replaced\n";
            }
        }
        else {
            ofstream uniquespeciestrf(filename2);
            uniquespeciestrf << "#nexus\n\n";
            uniquespeciestrf << "begin trees;\n";
            if (G::_verbose > 0) {
                cout << "created new file " << filename2 << "\n";
            }
        }
        if (filesystem::remove(filename3)) {
            ofstream paramsf(filename3);
            if (G::_verbose > 0) {
               cout << "existing file " << filename3 << " removed and replaced\n";
            }
        }
        else {
            ofstream paramsf(filename3);
            if (G::_verbose > 0) {
                cout << "created new file " << filename3 << "\n";
            }
        }
        
        cout << "\n";
        
        string altfname = "alt_species_trees.trees";
        if (filesystem::remove(altfname)) {
            if (G::_verbose > 0) {
               cout << "existing file " << altfname << " removed and replaced\n";
            }
            
            ofstream alttrf(altfname);

            alttrf << "#nexus\n\n";
            alttrf << "Begin trees;" << endl;
            string translate_block = my_vec[0].getTranslateBlock();
            alttrf << translate_block << endl;
        }
        else {
            ofstream alttrf(altfname);

            alttrf << "#nexus\n\n";
            alttrf << "Begin trees;" << endl;
            string translate_block = my_vec[0].getTranslateBlock();
            alttrf << translate_block << endl;
            
            if (G::_verbose > 0) {
                cout << "created new file " << altfname << "\n";
            }
        }

//        assert (G::_thin == 1.0);
//        unsigned ngroups = round(G::_nparticles * G::_ngroups * G::_thin);
//
//        assert(my_vec.size() == ngroups);
        
        secondLevel(my_vec);
        if (G::_save_gene_trees) {
            saveGeneTrees(my_vec);
        }
    }

    inline double Proj::getRunningSum(const vector<double> & log_weight_vec) const {
        double running_sum = 0.0;
        double log_particle_sum = 0.0;

        double log_max_weight = *max_element(log_weight_vec.begin(), log_weight_vec.end());
        for (auto & i:log_weight_vec) {
            running_sum += exp(i - log_max_weight);
        }
        log_particle_sum = log(running_sum) + log_max_weight;

        return log_particle_sum;
    }

    inline double Proj::computeEffectiveSampleSize(const vector<double> & probs) const {
        double ss = 0.0;
        for_each(probs.begin(), probs.end(), [&ss](double w){ss += w*w;});
        double ess = 1.0/ss;
        return ess;
    }

    inline unsigned Proj::multinomialDraw(const vector<double> & probs) {
        // Compute cumulative probababilities
        vector<double> cumprobs(probs.size());
        partial_sum(probs.begin(), probs.end(), cumprobs.begin());
        assert(fabs(*(cumprobs.rbegin()) - 1.0) < 0.0001);

        // Draw a Uniform(0,1) random deviate
        double u = rng.uniform();

        // Find first element in cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6 <-- u
        //
        // cumprobs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(cumprobs.begin(), cumprobs.end(), [u](double cumpr){return cumpr > u;});
        if (it == cumprobs.end()) {
            double last_cumprob = *(cumprobs.rbegin());
            throw XProj(format("G::multinomialDraw failed: u = %.9f, last cumprob = %.9f") % u % last_cumprob);
        }

        auto d = std::distance(cumprobs.begin(), it);
        assert(d >= 0);
        assert(d < probs.size());
        return (unsigned)d;
    }

    inline void Proj::filterParticlesRange(unsigned first, unsigned last, vector<Particle> &particles, unsigned g, vector<unsigned> particle_indices) {
        for (unsigned i=first; i<last; i++){
            // i is the group number
            unsigned start = i * G::_nparticles;
            unsigned end = start + (G::_nparticles) - 1;
            filterParticles(g, particles, particle_indices, start, end);
        }
    }

    inline void Proj::filterParticlesThreading(vector<Particle> &particles, unsigned g, vector<unsigned> particle_indices) {
        assert(G::_nthreads > 0);
        
        if (G::_nthreads == 1) {
            for (unsigned i=0; i<G::_ngroups; i++) {
                unsigned start = i * G::_nparticles;
                unsigned end = start + (G::_nparticles) - 1;
                            
                double ess = filterParticles(g, particles, particle_indices, start, end);
                
                if (G::_verbose > 1) {
                    cout << "   " << "ESS = " << ess << endl;
                }
            }
        }
        else {
          // divide up the particles as evenly as possible across threads
            unsigned first = 0;
            unsigned last = 0;
            unsigned stride = (G::_ngroups) / G::_nthreads; // divisor
            unsigned r = (G::_ngroups) % G::_nthreads; // remainder
            
            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;

            for (unsigned i=0; i<G::_nthreads; i++) {
                first = last;
                last = first + stride;
                
                if (r > 0) {
                    last += 1;
                    r -= 1;
                }
                
                if (last > (G::_ngroups)) {
                    last = (G::_ngroups);
                }
                
                threads.push_back(thread(&Proj::filterParticlesRange, this, first, last, std::ref(particles), g, particle_indices));
            }


          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }


    inline double Proj::filterParticles(unsigned step, vector<Particle> & particles, vector<unsigned> &particle_indices, unsigned start, unsigned end) {
        // Copy log weights for all bundles to prob vector
        
        vector<double> probs(G::_nparticles, 0.0);

        unsigned count_index = 0;
        for (unsigned p=start; p<end+1; p++) {
            unsigned particle_number = particle_indices[p];
            probs[count_index] = particles[particle_number].getLogWeight();
            count_index++;
        }

        assert (probs.size() == G::_nparticles);
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = getRunningSum(probs);

        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});

        // Compute component of the log marginal likelihood due to this step
        _log_marginal_likelihood += log_sum_weights - log(G::_nparticles);
        
//        if (step < G::_nloci) {
//            _log_marginal_likelihood += _starting_log_likelihood;
//        }

        double ess = 0.0;
        if (G::_verbose > 1) {
            // Compute effective sample size
            ess = computeEffectiveSampleSize(probs);
        }
        
//        if (ess < 200.0) {
        
#if defined(SYSTEMATIC_FILTERING)
        vector<unsigned> zeros;
        zeros.reserve(G::_nparticles);
        vector<unsigned> nonzeros;
        nonzeros.reserve(G::_nparticles);
        
        unsigned group_number = start / G::_nparticles;
        
        // Zero vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (G::_nparticles, 0);
        
        double cump = probs[0];
        double delta = _group_rng[group_number]->uniform() / G::_nparticles;
        unsigned c = (unsigned)(floor(1.0 + G::_nparticles*(cump - delta)));
        if (c > 0) {
            nonzeros.push_back(0);
        }
        else {
            zeros.push_back(0);
        }
        counts[0] = c;
        unsigned prev_cum_count = c;
        for (unsigned i = 1; i < G::_nparticles; ++i) {
            cump += probs[i];
            double cum_count = floor(1.0 + G::_nparticles*(cump - delta));
            if (cum_count > G::_nparticles) {
                cum_count = G::_nparticles;
            }
            unsigned c = (unsigned)cum_count - prev_cum_count;
            if (c > 0) {
                nonzeros.push_back(i);
            }
            else {
                zeros.push_back(i);
            }
            counts[i] = c;
            prev_cum_count = cum_count;
        }
        
        unsigned locus = particles[particle_indices[start]].getNextGene() - 1; // subtract 1 because vector of gene forests starts at 0
        assert (locus == particles[particle_indices[end]].getNextGene() - 1);
        // Create map (nonzero_map) in which the key for an element
        // is the memory address of a gene forest and
        // the value is a vector of indices of non-zero counts.
        // This map is used to determine which of the nonzeros
        // that need to be copied (last nonzero count for any
        // memory address does not need to be copied and can be
        // modified in place).
        map<const void *, list<unsigned> > nonzero_map;
        buildNonzeroMap(particles, locus, nonzero_map, nonzeros, particle_indices, start, end);
        
        // Example of following code that replaces dead
        // particles with copies of surviving particles:
        //             0  1  2  3  4  5  6  7  8  9
        // _counts  = {0, 2, 0, 0, 0, 8, 0, 0, 0, 0}  size = 10
        // zeros    = {0, 2, 3, 4, 6, 7, 8, 9}        size =  8
        // nonzeros = {1, 5}                          size =  2
        //
        //  next_zero   next_nonzero   k   copy action taken
        //  --------------------------------------------------------------
        //      0             0        0   _particles[1] --> _particles[0]
        //  --------------------------------------------------------------
        //      1             1        0   _particles[5] --> _particles[2]
        //      2             1        1   _particles[5] --> _particles[3]
        //      3             1        2   _particles[5] --> _particles[4]
        //      4             1        3   _particles[5] --> _particles[6]
        //      5             1        4   _particles[5] --> _particles[7]
        //      6             1        5   _particles[5] --> _particles[8]
        //      7             1        6   _particles[5] --> _particles[9]
        //  --------------------------------------------------------------
        unsigned next_zero = 0;
        unsigned next_nonzero = 0;
            while (next_nonzero < nonzeros.size()) {
                double index_survivor = nonzeros[next_nonzero];
                unsigned index_survivor_in_particles = particle_indices[index_survivor+start];
                
                if (!G::_mcmc) {
                    // for mcmc, wait to finalize joins
                    particles[index_survivor_in_particles].finalizeLatestJoin(locus, index_survivor_in_particles, nonzero_map);
                }
                unsigned ncopies = counts[index_survivor] - 1;
                for (unsigned k = 0; k < ncopies; k++) {
                    double index_nonsurvivor = zeros[next_zero++];
                    
                    // Replace non-survivor with copy of survivor
                    unsigned survivor_index_in_particles = particle_indices[index_survivor+start];
                    unsigned non_survivor_index_in_particles = particle_indices[index_nonsurvivor+start];
                    
                    particles[non_survivor_index_in_particles] = particles[survivor_index_in_particles];
                }
                
                ++next_nonzero;
            }
        
        return ess;
#else

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (G::_nparticles, 0);

        // Throw _nparticles darts
      for (unsigned i=start; i<end+1; i++) {
          unsigned group_number = start / G::_nparticles;
          double u = _group_rng[id_nuber]->uniform();
          auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
          assert(it != probs.end());
          unsigned which = (unsigned)std::distance(probs.begin(), it);
          counts[which]++;
        }
        

        // Copy particles

        bool copying_needed = true;

        // Locate first donor
        unsigned donor = 0;
        while (counts[donor] < 2) {
            donor++;
            if (donor >= counts.size()) {
                copying_needed = false; // all the particle counts are 1
                break;
            }
        }

        if (copying_needed) {
            // Locate first recipient
            unsigned recipient = 0;
            while (counts[recipient] != 0) {
                recipient++;
            }

            // Count number of cells with zero count that can serve as copy recipients
            unsigned nzeros = 0;
            for (unsigned i = 0; i < G::_nparticles; i++) {
                if (counts[i] == 0)
                    nzeros++;
            }

            while (nzeros > 0) {
                assert(donor < G::_nparticles);
                assert(recipient < G::_nparticles);

                // Copy donor to recipient
                unsigned recipient_index = particle_indices[recipient+start];
                unsigned donor_index = particle_indices[donor+start];

                particles[recipient_index] = particles[donor_index];

                counts[donor]--;
                counts[recipient]++;
                nzeros--;

                if (counts[donor] == 1) {
                    // Move donor to next slot with count > 1
                    donor++;
                    while (donor < G::_nparticles && counts[donor] < 2) {
                        donor++;
                    }
                }

                // Move recipient to next slot with count equal to 0
                recipient++;
                while (recipient < G::_nparticles && counts[recipient] > 0) {
                    recipient++;
                }
            }
      }
        return ess;
#endif

    }

    inline void Proj::buildNonzeroMap(vector<Particle> &particles, unsigned locus, map<const void *, list<unsigned> > & nonzero_map, const vector<unsigned> & nonzeros, vector<unsigned> particle_indices, unsigned start, unsigned end) {
        
        for (auto i : nonzeros) {
            unsigned particle_index = particle_indices[start + i];
            void * ptr = particles[particle_index].getGeneForestPtr(locus).get();

            if (nonzero_map.count(ptr) > 0) {
                nonzero_map[ptr].push_back(particle_index);
            }
            else {
                nonzero_map[ptr] = {particle_index};
            }
        }
    }

    inline double Proj::filterSpeciesParticles(unsigned step, vector<Particle> & particles, unsigned id_number) {
        unsigned nparticles = (unsigned) particles.size();
        assert (nparticles == G::_particle_increase);
        // Copy log weights for all bundles to prob vector
        vector<double> probs(G::_particle_increase, 0.0);
        
        for (unsigned p=0; p < nparticles; p++) {
            probs[p] = particles[p].getSpeciesLogWeight();
        }

        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = getRunningSum(probs);
        
        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});

    //        // Compute component of the log marginal likelihood due to this step
    //        _species_log_marginal_likelihood += log_sum_weights - log(nparticles);
    //        if (step == 0) {
    //            _species_log_marginal_likelihood += 0;
    //        }
        
        double ess = 0.0;
        if (G::_verbose > 1) {
        // Compute effective sample size
            ess = computeEffectiveSampleSize(probs);
        }
        
#if defined(SYSTEMATIC_FILTERING)
        vector<unsigned> zeros;
        zeros.reserve(G::_particle_increase);
        vector<unsigned> nonzeros;
        nonzeros.reserve(G::_particle_increase);

        // Zero vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (G::_particle_increase, 0);

        double cump = probs[0];

//        unsigned group_number = particles[0].getGroupNumber();
        double delta = _group_rng[id_number]->uniform() / G::_particle_increase;

        unsigned c = (unsigned)(floor(1.0 + G::_particle_increase*(cump - delta)));
        if (c > 0) {
            nonzeros.push_back(0);
        }
        else {
            zeros.push_back(0);
        }
        counts[0] = c;
        unsigned prev_cum_count = c;
        for (unsigned i = 1; i < G::_particle_increase; ++i) {
            cump += probs[i];
            double cum_count = floor(1.0 + G::_particle_increase*(cump - delta));
            if (cum_count > G::_particle_increase) {
                cum_count = G::_particle_increase;
            }
            unsigned c = (unsigned)cum_count - prev_cum_count;
            if (c > 0) {
                nonzeros.push_back(i);
            }
            else {
                zeros.push_back(i);
            }
            counts[i] = c;
            prev_cum_count = cum_count;
        }


        // Example of following code that replaces dead
        // particles with copies of surviving particles:
        //             0  1  2  3  4  5  6  7  8  9
        // _counts  = {0, 2, 0, 0, 0, 8, 0, 0, 0, 0}  size = 10
        // zeros    = {0, 2, 3, 4, 6, 7, 8, 9}        size =  8
        // nonzeros = {1, 5}                          size =  2
        //
        //  next_zero   next_nonzero   k   copy action taken
        //  --------------------------------------------------------------
        //      0             0        0   _particles[1] --> _particles[0]
        //  --------------------------------------------------------------
        //      1             1        0   _particles[5] --> _particles[2]
        //      2             1        1   _particles[5] --> _particles[3]
        //      3             1        2   _particles[5] --> _particles[4]
        //      4             1        3   _particles[5] --> _particles[6]
        //      5             1        4   _particles[5] --> _particles[7]
        //      6             1        5   _particles[5] --> _particles[8]
        //      7             1        6   _particles[5] --> _particles[9]
        //  --------------------------------------------------------------
        unsigned next_zero = 0;
        unsigned next_nonzero = 0;
        while (next_nonzero < nonzeros.size()) {
            double index_survivor = nonzeros[next_nonzero];
            unsigned ncopies = counts[index_survivor] - 1;
            for (unsigned k = 0; k < ncopies; k++) {
                double index_nonsurvivor = zeros[next_zero++];

                // Replace non-survivor with copy of survivor
                particles[index_nonsurvivor] = particles[index_survivor];
            }

            ++next_nonzero;
        }

        return ess;
#else

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (nparticles, 0);

        // Throw _nparticles darts
//        unsigned group_number = particles[0].getGroupNumber();
        
        for (unsigned i=0; i<nparticles; i++) {
            double u = _group_rng[id_number]->uniform();
            auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
            assert(it != probs.end());
            unsigned which = (unsigned)std::distance(probs.begin(), it);
            counts[which]++;
        }
        
        // Copy particles

      bool copying_needed = true;
      
        // Locate first donor
        unsigned donor = 0;
        while (counts[donor] < 2) {
            donor++;
            if (donor >= counts.size()) {
                copying_needed = false; // all the particle counts are 1
                break;
            }
        }
        
        if (copying_needed) {
            // Locate first recipient
            unsigned recipient = 0;
            while (counts[recipient] != 0) {
                recipient++;
            }

            // Count number of cells with zero count that can serve as copy recipients
            unsigned nzeros = 0;
            for (unsigned i = 0; i < nparticles; i++) {
                if (counts[i] == 0)
                    nzeros++;
            }

            while (nzeros > 0) {
                assert(donor < nparticles);
                assert(recipient < nparticles);

                // Copy donor to recipient
                particles[recipient] = particles[donor];

                counts[donor]--;
                counts[recipient]++;
                nzeros--;

                if (counts[donor] == 1) {
                    // Move donor to next slot with count > 1
                    donor++;
                    while (donor < nparticles && counts[donor] < 2) {
                        donor++;
                    }
                }

                // Move recipient to next slot with count equal to 0
                recipient++;
                while (recipient < nparticles && counts[recipient] > 0) {
                    recipient++;
                }
            }
        }
        
        return ess;
#endif
    }

    inline string Proj::inventName(unsigned k, bool lower_case) {
        // If   0 <= k < 26, returns A, B, ..., Z,
        // If  26 <= k < 702, returns AA, AB, ..., ZZ,
        // If 702 <= k < 18278, returns AAA, AAB, ..., ZZZ, and so on.
        //
        // For example, k = 19009 yields ABCD:
        // ABCD 19009 = 26 + 26*26 + 26*26*26 + 0*26*26*26 + 1*26*26 + 2*26 + 3
        //              <------- base ------>   ^first       ^second   ^third ^fourth
        // base = (26^4 - 1)/25 - 1 = 18278
        //   26^1 + 26^2 + 26^3 = 26^0 + 26^1 + 26^2 + 26^3 - 1 = (q^n - 1)/(q - 1) - 1, where q = 26, n = 4
        //   n = 1 + floor(log(19009)/log(26))
        // fourth = ((19009 - 18278                           )/26^0) % 26 = 3
        // third  = ((19009 - 18278 - 3*26^0                  )/26^1) % 26 = 2
        // second = ((19009 - 18278 - 3*26^0 - 2*26^1         )/26^2) % 26 = 1
        // first  = ((19009 - 18278 - 3*26^0 - 2*26^1 - 1*26^2)/26^3) % 26 = 0

        // Find how long a species name string must be
        double logibase26 = log(k)/log(26);
        unsigned n = 1 + (unsigned)floor(logibase26);
        vector<char> letters;
        unsigned base = (unsigned)((pow(26,n) - 1)/25.0 - 1);
        unsigned cum = 0;
        int ordA = (unsigned)(lower_case ? 'a' : 'A');
        for (unsigned i = 0; i < n; ++i) {
            unsigned ordi = (unsigned)((k - base - cum)/pow(26,i)) % 26;
            letters.push_back(char(ordA + ordi));
            cum += (unsigned)(ordi*pow(26,i));
        }
        string species_name(letters.rbegin(), letters.rend());
        return species_name;
    }

    inline void Proj::simSpeciesMap() {
        // nspecies is _sim_nspecies
        // ntaxa vector is _ntaxaperspecies
        unsigned count = 0;
        for (int s=0; s<G::_sim_nspecies; s++) {
            string species_name;
            species_name = inventName(s, false);
            for (int t=0; t<G::_ntaxaperspecies[s]; t++) {
                string taxon_name = inventName(count, true) + "^" + species_name;
                _taxon_map.insert({taxon_name, species_name});
                count++;
            }
            G::_species_names.push_back(species_name);
        } // TODO: if > 26 taxa, these wont' be grouped by species which will affect splits
    }

    inline void Proj::createSpeciesMap(bool use_data) {
        // this only works if names are in taxon^species format (no _)
        if (use_data) {
            G::_taxon_names = _data->getTaxonNames();
            for (auto &name:G::_taxon_names) {
                regex re(".+\\^(.+)");
                smatch match_obj;
                bool matched=regex_match(name, match_obj, re); //search name for regular expression, store result in match_obj
                cout << "taxon name: " << name << endl;
                if (matched) {
                    string species_name = match_obj[1];
                    string taxon_name = name;
                    if (find(G::_species_names.begin(), G::_species_names.end(), species_name) == G::_species_names.end()) {
                        G::_species_names.push_back(species_name);
                    }
                    _taxon_map[taxon_name]=species_name;
                }
            }
        }
        else {
            string prev_name = "";
            G::_taxon_names = _data->getTaxonNames();
            unsigned count = 0;
            for (auto &name:G::_taxon_names) {
                regex re(".+\\^(.+)");
                smatch match_obj;
                bool matched=regex_match(name, match_obj, re); //search name for regular expression, store result in match_obj
                cout << "taxon name: " << name << endl;
                if (matched) {
                    string species_name = match_obj[1];
                    string taxon_name = name;
                    _taxon_map[taxon_name]=G::_species_names[count];
                }
                if (count > 0) {
                    prev_name = G::_species_names[count-1];
                }
                if (prev_name != G::_species_names[count]) {
                    count++;
                }
            }
        }
    }

    inline void Proj::createSpeciesMapTyped() {
        assert (G::_nspecies > 0);
        G::species_t species_name = 1;
        G::species_t prev_species_name = 1;
        for (unsigned i=0; i<G::_nspecies; i++) {
            G::_species_names_typed.push_back(species_name);
            prev_species_name = species_name;
            species_name = prev_species_name + species_name;
        }
    }

    inline void Proj::showFinal(vector<Particle> my_vec) {
        for (auto &p:my_vec){
            p.showParticle();
        }

        double sum_h = 0.0;
        for (auto & p:my_vec) {
            double h = p.calcHeight();
            sum_h += h;
        }
        sum_h/=my_vec.size();
        cout << "mean height equals " << sum_h << endl;
        cout << "log marginal likelihood = " << _log_marginal_likelihood << endl;

#if defined (DRAW_NEW_THETA)
        cout << "different theta for each population in each particle " << endl;
#else
        cout << "theta = " << G::_theta << endl;
#endif

//        cout << "speciation rate = " << Forest::_lambda << endl;
    }

    inline void Proj::proposeSpeciesGroupRangeFaster(unsigned first, unsigned last, vector<Particle> &particles, unsigned ngroups, string filename1, string filename2, string filename3) {
        
        for (unsigned i=first; i<last; i++){
            vector<Particle> second_level_particles;
            
            unsigned group_number = _particle_indices_to_thin[i];
            unsigned id_number = i; // use for random seed, _second_level_particles_to_keep

            Particle p = particles[group_number];

            
            if (!G::_gene_newicks_specified) { // if starting from gene newicks, this is already built
                G::_generation = 0;
                p.resetSpecies();
                p.mapSpecies(_taxon_map);
            }

            second_level_particles.resize(G::_particle_increase, p);

            assert(second_level_particles.size() == G::_particle_increase);
        
            for (unsigned s=0; s<G::_nspecies-1; s++) {  // skip last round of filtering because weights are always 0
                
                // set particle random number seeds
                unsigned psuffix = 1;
                for (auto &p:second_level_particles) {
                    p.setSeed(_group_rng[id_number]->randint(1,9999) + psuffix);
                    psuffix += 2;
                }
                
                for (auto &p:second_level_particles) {
                    p.proposeSpeciationEvent();
                }
                
                filterSpeciesParticles(s, second_level_particles, id_number);
                G::_generation++;
            }
            
//            if (G::_save_every > 1.0) { // thin sample for output by taking a random sample
                unsigned sample_size = round (double (G::_particle_increase) / double(G::_save_every));
                if (sample_size == 0) {
                    cout << "\n";
                    cout << "current settings would save 0 species trees; saving every species tree\n";
                    cout << "\n";
                    sample_size = G::_particle_increase;
                }

                unsigned group_seed = _group_rng[id_number]->getSeed();
                

                unsigned count = 0;
                for (unsigned p = 0; p < G::_particle_increase; p++) {
                    _second_level_indices_to_keep[id_number].push_back(count);
                    count++;
                }
                
                std::shuffle(_second_level_indices_to_keep[id_number].begin(), _second_level_indices_to_keep[id_number].end(), std::default_random_engine(group_seed)); // shuffle particles using group seed
                
                // delete first (1-_thin) % of indices to keep
                _second_level_indices_to_keep[id_number].erase(next(_second_level_indices_to_keep[id_number].begin(), 0), next(_second_level_indices_to_keep[id_number].begin(), (G::_particle_increase-sample_size)));
                assert (_second_level_indices_to_keep[id_number].size() == sample_size);
//        }

            mtx.lock();
            saveSpeciesTreesHierarchical(second_level_particles, filename1, filename2, id_number);
            saveSpeciesTreesAltHierarchical(second_level_particles, id_number);
            _count++;
//            assert (i == group_number);
            if (G::_gene_newicks_specified) {
                writeParamsFileForBeastComparisonAfterSpeciesFilteringSpeciesOnly(second_level_particles, filename3, id_number);
            }
            else {
                writeParamsFileForBeastComparisonAfterSpeciesFiltering(second_level_particles, filename3, id_number);
            }
            mtx.unlock();
        }
    }

    inline void Proj::proposeSpeciesGroupsFaster(vector<Particle> &particles, unsigned ngroups, string filename1, string filename2, string filename3) {
        // ngroups = number of species SMCs to do (i.e. 100 particles for first round, thin = 1.0 means ngroups = 100 for this round)
        cout << "starting proposals second level" << endl;

        assert (G::_nthreads > 1);
        
        // divide up groups as evenly as possible across threads
        unsigned first = 0;
        unsigned incr = ngroups/G::_nthreads + (ngroups % G::_nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of groups
        unsigned last = incr;
        
        // need a vector of threads because we have to wait for each one to finish
        vector<thread> threads;

          while (true) {
          // create a thread to handle particles first through last - 1
            threads.push_back(thread(&Proj::proposeSpeciesGroupRangeFaster, this, first, last, std::ref(particles), ngroups, filename1, filename2, filename3));
          // update first and last
          first = last;
          last += incr;
          if (last > ngroups) {
            last = ngroups;
            }
          if (first >= ngroups) {
              break;
          }
        }

        // the join function causes this loop to pause until the ith thread finishes
        for (unsigned i = 0; i < threads.size(); i++) {
          threads[i].join();
        }
    }

    inline void Proj::initializeParticle(Particle &particle) {
        // set partials for first particle under save_memory setting for initial marginal likelihood calculation
        assert (G::_nthreads > 0);
        
        if (G::_plus_G) {
            double pinvar = 0.0;
            if (G::_plus_I) {
                pinvar = G::_pinvar;
            }
            // alpha and beta are shape and scale, respectively
            // mean = alpha * beta = 1.0
            double rate_variance = G::_gamma_rate_var;
            double alpha = 1 / G::_gamma_rate_var;
            double beta = rate_variance;
            double num_categ = 4.0; // TODO: fix this
            double mean_rate_variable_sites = 1.0;
            if (G::_plus_I) {
                mean_rate_variable_sites /= (1.0 - G::_pinvar);
            }
            double equal_prob = 1 / num_categ;
            
            boost::math::gamma_distribution<> my_gamma(alpha, beta);
            boost::math::gamma_distribution<> my_gamma_plus(alpha + 1.0, beta);

              double cum_upper        = 0.0;
              double cum_upper_plus   = 0.0;
              double upper            = 0.0;
              double cum_prob         = 0.0;
              for (unsigned i = 1; i <= num_categ; ++i) {
                  double cum_lower_plus       = cum_upper_plus;
                  double cum_lower            = cum_upper;
                  cum_prob                    += equal_prob;

                  if (i < num_categ) {
                      upper                   = boost::math::quantile(my_gamma, cum_prob);
                      cum_upper_plus          = boost::math::cdf(my_gamma_plus, upper);
                      cum_upper               = boost::math::cdf(my_gamma, upper);
                  }
                  else {
                      cum_upper_plus          = 1.0;
                      cum_upper               = 1.0;
                  }

                  double numer                = cum_upper_plus - cum_lower_plus;
                  double denom                = cum_upper - cum_lower;
                  double r_mean               = (denom > 0.0 ? (alpha*beta*numer/denom) : 0.0);
                  double mean = r_mean * mean_rate_variable_sites;
                  if ((mean - 0) < 0.001) {
                      mean = 0.001;
                  }
                  G::_gamma_rate_cat.push_back(mean);
//                  G::_gamma_rate_cat.push_back(1.0);
              }
        }
        else {
            double mean_rate_variable_sites = 1.0;
            if (G::_plus_I) {
                mean_rate_variable_sites /= (1.0 - G::_pinvar);
            }
            G::_gamma_rate_cat.push_back(mean_rate_variable_sites);
        }
        
        bool partials = true;
        if (G::_gene_newicks_specified) {
            partials = false;
            G::_save_memory = true;
        }
        
        if (G::_save_memory) {
            partials = false; // no partials needed for running on empty
        }
        
        particle.setData(_data, _taxon_map, partials);
        particle.mapSpecies(_taxon_map);
                
        if (!G::_gene_newicks_specified) { // if starting from gene newicks, don't calculate any likelihoods
            particle.calcGeneTreeLogLikelihoods();
        }
        
//        if (G::_species_newick_name != "null") {
//            string species_newick = handleSpeciesNewick();
//            particle.processSpeciesNewick(species_newick);
//        }
        
        if (G::_fix_theta) {
            particle.fixTheta();
        }
    }

    inline void Proj::proposeParticlesParallelizeByGroup(vector<vector<Particle>> &particles) {
        assert (G::_nthreads > 0);
                
        if (G::_nthreads == 1) {
            for (unsigned i=0; i<G::_ngroups; i++) {
                for (auto &p:particles[i]) {
                    p.proposal();
                }
            }
        }
        else {
            // run each group separately
            // divide up the groups as evenly as possible across threads
            unsigned first = 0;
            unsigned last = 0;
            unsigned stride = G::_ngroups / G::_nthreads; // divisor
            unsigned r = G::_ngroups % G::_nthreads; // remainder

            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;

            for (unsigned i=0; i<G::_nthreads; i++) {
                first = last;
                last = first + stride;

                if (r > 0) {
                    last += 1;
                    r -= 1;
                }

                if (last > G::_ngroups) {
                    last = G::_ngroups;
                }

                threads.push_back(thread(&Proj::proposeParticleGroupRange, this, first, last, std::ref(particles)));
            }


          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }

    inline void Proj::proposeParticlesSim(vector<Particle> &particles) {
        assert(G::_nthreads > 0); // don't thread for simulations
        for (auto & p : particles) {
            p.proposalSim();
        }
    }

    inline void Proj::mcmcMoves(vector<Particle> &particles, bool last_round) {
        assert(G::_nthreads > 0);
        if (G::_nthreads == 1) {
            unsigned count = 0;
            for (auto &p:particles) {
                unsigned prev_n_mcmc = G::_nmcmc_moves_accepted;
                p.proposeMCMCMove(last_round);
                
                bool accepted = false;
                if (G::_nmcmc_moves_accepted > prev_n_mcmc) {
                    accepted = true;
                }
                
                bool tune = false;
                
                double target_acceptance = 0.1;
                unsigned nattempts = count + 1;
                if (tune) {
                    // tune sliding window
                    double gamma_n = 10.0/(100.0 + nattempts);
                    if (accepted) {
                        G::_sliding_window *= 1.0 + gamma_n*(1.0 - target_acceptance)/(2.0*target_acceptance);
                    }
                    else {
                        G::_sliding_window *= 1.0 - gamma_n*0.5;
                    }

                    // Prevent run-away increases in boldness for low-information marginal densities
                    if (G::_sliding_window > 1000.0) {
                        G::_sliding_window = 1000.0;
                    }
                    // TODO: opposite problem - window gets too small?
                }
                count++;
            }
        }
        else { // TODO: no groups
            // divide up the particles as evenly as possible across threads
              unsigned first = 0;
              unsigned last = 0;
              unsigned stride = (G::_nparticles*G::_ngroups) / G::_nthreads; // divisor
              unsigned r = (G::_nparticles*G::_ngroups) % G::_nthreads; // remainder

            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;
            
            for (unsigned i=0; i<G::_nthreads; i++) {
                first = last;
                last = first + stride;
                
                if (r > 0) {
                    last += 1;
                    r -= 1;
                }
                
                if (last > (G::_nparticles*G::_ngroups)) {
                    last = (G::_nparticles*G::_ngroups);
                }
                
                threads.push_back(thread(&Proj::proposeMCMCMoveRange, this, first, last, std::ref(particles), last_round));
                // TODO: no tuning if parallelizing
            }


          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }

        }
    }

    inline void Proj::proposeParticles(vector<Particle> &particles) {
        assert(G::_nthreads > 0);
        if (G::_nthreads == 1) {
          for (auto & p : particles) {
              p.proposal();
          }
        }
        else {
          // divide up the particles as evenly as possible across threads
            unsigned first = 0;
            unsigned last = 0;
            unsigned stride = (G::_nparticles*G::_ngroups) / G::_nthreads; // divisor
            unsigned r = (G::_nparticles*G::_ngroups) % G::_nthreads; // remainder
            
            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;

            for (unsigned i=0; i<G::_nthreads; i++) {
                first = last;
                last = first + stride;
                
                if (r > 0) {
                    last += 1;
                    r -= 1;
                }
                
                if (last > (G::_nparticles*G::_ngroups)) {
                    last = (G::_nparticles*G::_ngroups);
                }
                
                threads.push_back(thread(&Proj::proposeParticleRange, this, first, last, std::ref(particles)));
            }


          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }

    inline void Proj::proposeParticleGroupRange(unsigned first, unsigned last, vector<vector<Particle>> &particles) {
        for (unsigned i=first; i<last; i++){
            for (auto &p:particles[i]) {
                p.proposal();
            }
        }
    }

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, vector<Particle> &particles) {
        for (unsigned i=first; i<last; i++){
            particles[i].proposal();
        }
    }

    inline void Proj::proposeMCMCMoveRange(unsigned first, unsigned last, vector<Particle> &particles, bool last_round) {
        for (unsigned i=first; i<last; i++){
            particles[i].proposeMCMCMove(last_round);
        }
    }

#if defined (DEBUG_MODE)
    inline void Proj::debugSpeciesTree(vector<Particle> &particles) {
        cout << "debugging species tree" << endl;
        for (auto &p:particles) {
            p.showSpeciesJoined();
            p.showSpeciesIncrement();
            p.showSpeciesTree();
            cout << " _______ " << endl;
        }
    }
#endif

    inline void Proj::writePaupFile(vector<Particle> particles, vector<string> taxpartition) {
        // Output a PAUP* command file for estimating the species tree using
        // svd quartets and qage
        cout << "  PAUP* commands saved in file \"svd-qage.nex\"\n";
        ofstream paupf("svd-qage.nex");
        paupf << "#NEXUS\n\n";
        paupf << "begin paup;\n";
        paupf << "  log start file=svdout.txt replace;\n";
        paupf << "  exe " + G::_sim_file_name + ";\n";
        paupf << "  taxpartition species (vector) = " << join(taxpartition," ") << ";\n";
        paupf << "svd taxpartition=species bootstrap nreps = 1000 treefile=test.tre replace;\n";
//        paupf << "  svd taxpartition=species;\n";
        paupf << "  roottrees;\n";
        paupf << "  qage taxpartition=species patprob=exactjc outUnits=substitutions treefile=svd.tre replace;\n";
        paupf << "  log stop;\n";
        paupf << "  quit;\n";
        paupf << "end;\n";
        paupf.close();
    }

    inline void Proj::simulateData() {
        cout << "\nSimulating data under multispecies coalescent model...\n" << endl;
        // set to one thread
        G::_nthreads = 1;
        rng.setSeed(_random_seed);

        unsigned ntaxa = 0;
        G::_nspecies = G::_sim_nspecies;
        
        if (G::_ntaxaperspecies.size() == 1) {
            ntaxa = G::_ntaxaperspecies[0];
            for (int i=0; i<G::_sim_nspecies-1; i++) {
                G::_ntaxaperspecies.push_back(ntaxa);
            }
            G::_ntaxa = ntaxa * G::_nspecies;
        }
        else {
            for (auto &t:G::_ntaxaperspecies) {
                ntaxa += t;
            }
            G::_ntaxa = ntaxa;
        }
        
        assert (G::_nspecies == (unsigned) G::_ntaxaperspecies.size());

        vector<Particle> sim_vec(1);
        sim_vec[0] = Particle();
        // set particle random number seed
        unsigned psuffix = 1;
        sim_vec[0].setSeed(rng.randint(1,9999) + psuffix);
        psuffix += 2;
        
        if (G::_plus_G) {
            double pinvar = 0.0;
            if (G::_plus_I) {
                pinvar = G::_pinvar;
            }
            // alpha and beta are shape and scale, respectively
            // mean = alpha * beta = 1.0
            double rate_variance = G::_gamma_rate_var;
            double alpha = 1 / G::_gamma_rate_var;
            double beta = rate_variance;
            double num_categ = 4.0; // TODO: fix this
            double mean_rate_variable_sites = 1.0;
            if (G::_plus_I) {
                mean_rate_variable_sites /= (1.0 - G::_pinvar);
            }
            double equal_prob = 1 / num_categ;
            
            boost::math::gamma_distribution<> my_gamma(alpha, beta);
            boost::math::gamma_distribution<> my_gamma_plus(alpha + 1.0, beta);

              double cum_upper        = 0.0;
              double cum_upper_plus   = 0.0;
              double upper            = 0.0;
              double cum_prob         = 0.0;
              for (unsigned i = 1; i <= num_categ; ++i) {
                  double cum_lower_plus       = cum_upper_plus;
                  double cum_lower            = cum_upper;
                  cum_prob                    += equal_prob;

                  if (i < num_categ) {
                      upper                   = boost::math::quantile(my_gamma, cum_prob);
                      cum_upper_plus          = boost::math::cdf(my_gamma_plus, upper);
                      cum_upper               = boost::math::cdf(my_gamma, upper);
                  }
                  else {
                      cum_upper_plus          = 1.0;
                      cum_upper               = 1.0;
                  }

                  double numer                = cum_upper_plus - cum_lower_plus;
                  double denom                = cum_upper - cum_lower;
                  double r_mean               = (denom > 0.0 ? (alpha*beta*numer/denom) : 0.0);
                  double mean = r_mean * mean_rate_variable_sites;
                  if ((mean - 0) < 0.001) {
                      mean = 0.001;
                  }
                  G::_gamma_rate_cat.push_back(mean);
//                  G::_gamma_rate_cat.push_back(1.0);
              }
        }
        else {
            double mean_rate_variable_sites = 1.0;
            if (G::_plus_I) {
                mean_rate_variable_sites /= (1.0 - G::_pinvar);
            }
            G::_gamma_rate_cat.push_back(mean_rate_variable_sites);
        }

        G::_run_on_empty = true;
        G::_proposal = "prior-prior";
        G::_prior_prior = true;

        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);

        // make up the species map
        simSpeciesMap();
        
#if defined (DRAW_NEW_THETA)
        updateSpeciesNames();
#endif
        
        createSpeciesMapTyped(); // use data

        vector<string> taxpartition;
        for (auto &t:_taxon_map) {
            taxpartition.push_back(t.second);
        }

        // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
        if (G::_outgroup != "none") {
            checkOutgroupName();
        }

        G::_nloci = _data->getNumSubsets();
        
        assert (ntaxa > 0);
        
        unsigned list_size = (G::_ntaxa-1)*G::_nloci;
        vector<unsigned> gene_order;
        
        unsigned count = 1;
        vector<pair<double, unsigned>> randomize;
        for (unsigned l=0; l<list_size; l++) {
            if (count == 1) {
                randomize.clear();
            }
            randomize.push_back(make_pair(rng.uniform(), count));
            count++;
            if (count > G::_nloci) {
                sort(randomize.begin(), randomize.end());
                for (auto &r:randomize) {
                    gene_order.push_back(r.second);
                }
                count = 1;
            }
        }
        
        assert (gene_order.size() == list_size);
        
        sim_vec[0].setGeneOrder(gene_order);

        sim_vec[0].setSimData(_data, _taxon_map, (unsigned) _taxon_map.size());

        sim_vec[0].mapSpecies(_taxon_map);

//        sim_vec[0].setNextSpeciesNumber(); // need to reset this now that number of species is known
        
        sim_vec[0].setNewTheta(G::_fix_theta_for_simulations);
                
        sim_vec[0].setNTaxaPerSpecies(G::_ntaxaperspecies);

        unsigned nsteps = (unsigned) (_taxon_map.size()-1)*G::_nloci;

        for (unsigned g=0; g<nsteps; g++){
            proposeParticlesSim(sim_vec);
            G::_generation++;
        }
        
        sim_vec[0].getNumDeepCoalescences();
        sim_vec[0].getMaxDeepCoalescences();

        cout << "\nBuilding species tree and associated gene trees....\n";
        G::_taxon_names.clear();
        for (auto &t:_taxon_map) {
            G::_taxon_names.push_back(t.first);
        }

        _data->setTaxonNames(G::_taxon_names);

        // Simulate sequence data
        cout << "\nSimulating sequence data....\n";
       vector<tuple<unsigned, unsigned, unsigned, unsigned>> sites_tuples = _partition->getSubsetRangeVect();
        vector<unsigned> sites_vector;
        for (auto &s:sites_tuples) {
            sites_vector.push_back(get<1>(s) - get<0>(s) + 1);
        }
        
        sim_vec[0].simulateData(sites_vector);
        
        _data->compressPatterns();
        _data->writeDataToFile(G::_sim_file_name);

        saveGeneTrees(sim_vec);
        saveSpeciesTrees(sim_vec);

        writePaupFile(sim_vec, taxpartition);
        writeDeepCoalescenceFile(sim_vec);
        writeThetaFile(sim_vec);
    }

    inline void Proj::secondLevel(vector<Particle> &particles) {
        if (G::_verbose > 1) {
            cout << "beginning second level SMC " << endl;
        }
        
        cout << "\n";
        string filename1 = "species_trees.trees";
        string filename2 = "unique_species_trees.trees";
        string filename3 = "params-beast-comparison.log";
        string altfname = "alt_species_trees.trees";

        cout << "\n";

        unsigned ngroups = round(G::_nparticles * G::_ngroups * G::_thin);
        if (ngroups == 0) {
            ngroups = 1;
            cout << "thin setting would result in 0 species groups; setting species groups to 1" << endl;
        }
        
        // set group rng
        unsigned psuffix = 1;
        _group_rng.resize(ngroups);
        psuffix = 1;
        for (auto &g:_group_rng) {
            g.reset(new Lot());
            g->setSeed(rng.randint(1,9999)+psuffix);
            psuffix += 2;
        }
        
        // save coal info for all existing gene pointers before starting any work to avoid timing issues
        for (auto &p:particles) {
            p.clearGeneForests();
        }
        
# if defined BEST_GENE_TREES
        vector<pair<unsigned, double>> particle_indices_and_posteriors;
        // TODO: try saving only particles with highest posteriors
        unsigned count = 0;
        for (auto &p:particles) {
            double log_coalescent_likelihood = 0.0;
            for (unsigned g=0; g<G::_nloci; g++) {
                log_coalescent_likelihood += p.getCoalescentLikelihood(g); // TODO: should this return 0 at this stage?
            }
            double log_likelihood = p.getLogLikelihood();;
            double log_prior = p.getAllPriorsFirstRound();
            
            double log_posterior = log_likelihood + log_prior + log_coalescent_likelihood;
            particle_indices_and_posteriors.push_back(make_pair(count, log_posterior));
            count++;
        }
        
        // sort particle_indices_and_posteriors by posterior
        std::sort(particle_indices_and_posteriors.begin(), particle_indices_and_posteriors.end(), [](const std::pair<unsigned, double>& a, const std::pair<unsigned, double>& b) {
            return b.second < a.second; // Compare based on the second element
        });
#endif
        
        _second_level_indices_to_keep.resize(ngroups);
        
        // TODO: trying multinomial resampling within each subgroup - this means thin will apply to each subgroup, not the entire thing
        unsigned ngroups_within_subgroup = round(G::_nparticles * G::_thin);
        if (ngroups_within_subgroup == 0) {
            ngroups_within_subgroup = 1;
        }
        
# if defined (BEST_GENE_TREES)
        // TODO: trying - just saving first (1 - thin)% of sorted particles
        unsigned index = 0;
        for (unsigned g = 0; g<G::_ngroups; g++) {
            for (unsigned i = 0; i <ngroups_within_subgroup; i++) {
                _particle_indices_to_thin.push_back(particle_indices_and_posteriors[index].first);
                index++;
            }
        }
#else
        for (unsigned g=0; g<G::_ngroups; g++) {
            for (unsigned i=0; i<ngroups_within_subgroup; i++) {
                unsigned start = g*G::_nparticles;
                unsigned end = start + (G::_nparticles) - 1;
                unsigned n = rng.randint(start, end);
//                unsigned n = rng.randint(0, ngroups - 1);
                _particle_indices_to_thin.push_back(n); // TODO: if particle indices are shuffled, is this different?
            }
        }
        
        // TODO: trying taking a random sample, not multinomial resampling
//            for (unsigned g=0; g<G::_ngroups; g++) {
//                vector<int> rand_numbers;
//                for (unsigned n=0; n<G::_nparticles-1; n++) {
//                    rand_numbers.push_back(n);
//                }
//                std::shuffle(rand_numbers.begin(), rand_numbers.end(), std::default_random_engine(_random_seed));
//                unsigned start = g*G::_nparticles;
//                for (unsigned i=0; i<ngroups_within_subgroup; i++) {
//                    _particle_indices_to_thin.push_back(rand_numbers[i] + start);
////                    unsigned start = g*G::_nparticles;
////                    unsigned end = start + (G::_nparticles) - 1;
////                    unsigned n = rng.randint(start, end);
////                    _particle_indices_to_thin.push_back(n); // TODO: if particle indices are shuffled, is this different?
//                }
//            }
#endif
        
        assert (_second_level_indices_to_keep.size() == ngroups);
        
        if (G::_nthreads == 1) {
            cout << "starting proposals second level" << endl;
            for (unsigned g=0; g<ngroups; g++) { // propose and filter for each particle saved from first round

//                unsigned group_number = g;
                unsigned group_number = _particle_indices_to_thin[g];
                unsigned id_number = g; // use for random seeds and in _second_level_indices_to_keep
                
                // grab a random index of a new particle
//                unsigned i = rng.randint(0, G::_nparticles - 1);
                Particle p = particles[group_number];
                
//                Particle p = particles[g];
                vector<Particle > second_level_particles;
                
                if (!G::_gene_newicks_specified) { // if starting from gene newicks, this is already built
                    G::_generation = 0;
                    p.resetSpecies();
                    p.mapSpecies(_taxon_map);
                }
            
                second_level_particles.resize(G::_particle_increase, p);
                                
                for (unsigned s=0; s<G::_nspecies-1; s++) {  // skip last round of filtering because weights are always 0
                    // set particle random number seeds
                    psuffix = 1;
                    for (auto &p:second_level_particles) {
                        p.setSeed(_group_rng[id_number]->randint(1,9999) + psuffix);
                        psuffix += 2;
                    }
                    
                    for (auto &p:second_level_particles) {
                        p.proposeSpeciationEvent();
                    }
                    
                    filterSpeciesParticles(s, second_level_particles, id_number);
                    
                    G::_generation++;
                
                }
                
//                if (G::_save_every > 1.0) { // thin sample for output by taking a random sample
                    unsigned sample_size = round (double (G::_particle_increase) / double(G::_save_every));
                    if (sample_size == 0) {
                        cout << "\n";
                        cout << "current settings would save 0 species trees; saving every species tree\n";
                        cout << "\n";
                        sample_size = G::_particle_increase;
                    }

                    unsigned group_seed = _group_rng[id_number]->getSeed();

                    unsigned count = 0;
                    for (unsigned p = 0; p < G::_particle_increase; p++) {
                        _second_level_indices_to_keep[id_number].push_back(count);
                        count++;
                    }
                    
                    std::shuffle(_second_level_indices_to_keep[id_number].begin(), _second_level_indices_to_keep[id_number].end(), std::default_random_engine(group_seed)); // shuffle particles using group seed
                    
                    // delete first (1-_thin) % of indices to keep
                    _second_level_indices_to_keep[id_number].erase(next(_second_level_indices_to_keep[id_number].begin(), 0), next(_second_level_indices_to_keep[id_number].begin(), (G::_particle_increase-sample_size)));
                    assert (_second_level_indices_to_keep[id_number].size() == sample_size);
//                }
                
                saveSpeciesTreesHierarchical(second_level_particles, filename1, filename2, id_number);
                saveSpeciesTreesAltHierarchical(second_level_particles, id_number);
                writeParamsFileForBeastComparisonAfterSpeciesFiltering(second_level_particles, filename3, id_number);
            }
//
            ofstream altfname_end;
            altfname_end.open("alt_species_trees.trees", std::ios::app);
            altfname_end << "end;" << endl;
            altfname_end.close();
            
            ofstream strees;
            strees.open("species_trees.trees", std::ios::app);
            strees << "end;" << endl;
            strees.close();
            
            ofstream u_strees;
            u_strees.open("unique_species_trees.trees", std::ios::app);
            u_strees << "end;" << endl;
            u_strees.close();
            
            string line;
            // For writing text file
            // Creating ofstream & ifstream class object
            ifstream in ("params-beast-comparison.log");
            ofstream f("params-beast-comparison-final.log");
            
            unsigned line_count = 0;

            while (!in.eof()) {
                string text;

                getline(in, text);

                if (line_count == 0) {
                    string add = "iter ";
                    text = add + text;
                }
                else {
                    if (text != "") {
                        string add = to_string(line_count);
                        text = add + text;
                    }
                }
                if (text != "") {
                    f << text << endl; // account for blank line at end of file
                }
                line_count++;
            }

            // remove existing params file and replace with copy
            char oldfname[] = "params-beast-comparison.log";
            char newfname[] = "params-beast-comparison-final.log";
            filesystem::remove(oldfname);
            std::rename(newfname, oldfname);
        }
        else {
            
            proposeSpeciesGroupsFaster(particles, ngroups, filename1, filename2, filename3);
            
            ofstream strees;
            strees.open("species_trees.trees", std::ios::app);
            strees << "end;" << endl;
            strees.close();
            
            ofstream u_strees;
            u_strees.open("unique_species_trees.trees", std::ios::app);
            u_strees << "end;" << endl;
            u_strees.close();

            string line;
            // For writing text file
            // Creating ofstream & ifstream class object
            ifstream in ("params-beast-comparison.log");
            ofstream f("params-beast-comparison-final.log");

            unsigned line_count = 0;

            while (!in.eof()) {
                string text;

                getline(in, text);

                if (line_count == 0) {
                    string add = "iter ";
                    text = add + text;
                }
                else {
                    if (text != "") {
                        string add = to_string(line_count);
                        text = add + text;
                    }
                }
                if (text != "") {
                    f << text << endl; // account for blank line at end of file
                }
                line_count++;
            }

            // remove existing params file and replace with copy
            char oldfname[] = "params-beast-comparison.log";
            char newfname[] = "params-beast-comparison-final.log";
            filesystem::remove(oldfname);
            std::rename(newfname, oldfname);
        }
        
        if (G::_ruv) {
            assert (_ranks.size() > 0);
            string sim_file_name;
            if (G::_newick_path != "") {
                sim_file_name = G::_newick_path + "/" + "true_species_tree_height.txt";
            }
            else {
                sim_file_name = "true_species_tree_height.txt";
            }
            string line;
            string height_as_string;
            ifstream infile(sim_file_name);
            while (getline(infile, line)) {
                height_as_string = line;
            }
            double true_species_tree_height = 0.0;
            if (height_as_string != "") {
                true_species_tree_height = std::stod(height_as_string);
            }
            _ranks.push_back(make_pair(true_species_tree_height, true));
            // sort ranks
            std::sort(_ranks.begin(), _ranks.end());
            
            // find rank of truth
            auto it = std::find_if(_ranks.begin(), _ranks.end(), [&](const pair<double, bool>& p) { return p.second == true;});
            unsigned index_value = (unsigned) std::distance(_ranks.begin(), it);
            
            // write rank value to file
            ofstream rankf("rank_species_tree_height.txt");
            rankf << "rank: " << index_value << endl;
        }
        
        if (G::_bhv_reference != "" || G::_bhv_reference_path != ".") {
            assert (_bhv_distances.size() > 0);
            string sim_file_name;
            if (G::_newick_path != "") {
                sim_file_name = G::_newick_path + "/" + "true-species-tree.tre";
            }
            else {
                sim_file_name = "true-species-tree.tre";
            }
            
            Particle p;
            string true_newick = readNewickFromFile(sim_file_name);
            
            double true_bhv = p.calcBHVDistanceTrueTree(true_newick);
            
            _bhv_distances.push_back(true_bhv);
            // sort distances
            std::sort(_bhv_distances.begin(), _bhv_distances.end());
            
            // find rank of truth
            auto it = std::find(_bhv_distances.begin(), _bhv_distances.end(), true_bhv);
            unsigned index_value = (unsigned) std::distance(_bhv_distances.begin(), it);
            
            // write rank value to file
            ofstream rankf("rank_bhv.txt");
            rankf << "rank: " << index_value << endl;
        }
        
        if (G::_hpd) {
           // species tree heights
            ofstream hpdf("hpd_heights.txt");
            hpdf << "min    " << "max " << endl;
            
            // sort hpd values largest to smallest
            std::sort(_hpd_values.begin(), _hpd_values.end());
            std::reverse(_hpd_values.begin(), _hpd_values.end());
            
            // take first 95% of values (round down to nearest integer)
            double total = size(_hpd_values);
            double ninety_five_index = floor(0.95*total);
            // TODO: double check this - should there be a log joining prob as part fo the species tree?
            if (ninety_five_index == 0) {
                ninety_five_index = 1;
            }
            
            vector<double> hpd_values_in_range;
            
            for (unsigned h=0; h<ninety_five_index; h++) {
                hpd_values_in_range.push_back(_hpd_values[h].second);
            }
            
            auto max = *std::max_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
            auto min = *std::min_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
            assert (min < max || min == max);
            
            // write min and max to file
            hpdf << min << "\t" << max << endl;
        }
        
        if (G::_hpd || G::_ruv) {
            // write mean species tree height to output file for validation
             double sum = accumulate(_species_tree_heights.begin(), _species_tree_heights.end(), 0.0);
             double mean = sum / _species_tree_heights.size();
             
             ofstream heightf("average_species_tree_height.txt");
             heightf << mean << endl;
        }
    }

    inline void Proj::buildSpeciesMap(bool taxa_from_data) {
        // Populates G::_species_names and G::_taxon_to_species
        // For example, given these taxon names
        //  t1^A, t2^A, t3^B, t4^B, t5^B, t6^C, t7^C, t8^D, t9^D, t10^D
        // G::_species_names = ["A", "B", "C", "D"]
        // G::_taxon_to_species["t1^A"]  = 0
        // G::_taxon_to_species["t2^A"]  = 0
        // G::_taxon_to_species["t3^B"]  = 1
        // G::_taxon_to_species["t4^B"]  = 1
        // G::_taxon_to_species["t5^B"]  = 1
        // G::_taxon_to_species["t6^C"]  = 2
        // G::_taxon_to_species["t7^C"]  = 2
        // G::_taxon_to_species["t8^D"]  = 3
        // G::_taxon_to_species["t9^D"]  = 3
        // G::_taxon_to_species["t10^D"] = 3
        // G::_taxon_to_species["A"]     = 0
        // G::_taxon_to_species["B"]     = 1
        // G::_taxon_to_species["C"]     = 2
        // G::_taxon_to_species["D"]     = 3
        unsigned nspecies = 0;
        map<string, unsigned> species_name_to_index;
        G::_taxon_to_species.clear();

//        output("\nMapping taxon names to species index:\n", 2);
        unsigned ntax = (unsigned)G::_taxon_names.size();
        assert(ntax > 0);
        if (taxa_from_data) {
            // Assume taxon names are already stored in _data object and no
            // species names have yet been stored
            G::_species_names.clear();
            assert(G::_taxon_names.size() > 0);
            for (auto & tname : G::_taxon_names) {
                string species_name = Node::taxonNameToSpeciesName(tname);
                unsigned species_index = ntax;
                if (species_name_to_index.find(species_name) == species_name_to_index.end()) {
                    // species_name not found
                    species_index = nspecies;
                    G::_species_names.push_back(species_name);
                    species_name_to_index[species_name] = nspecies++;
                } else {
                    // species_name found
                    species_index = species_name_to_index[species_name];
                }
                G::_taxon_to_species[tname] = species_index;
            }
        }
        else {
            // Assume G::_species_names and G::_taxon_names are already populated
            // and _data is null because no data is needed for this analysis ("spec" mode)
            
            // First build species_name_to_index from G::_species_names
            unsigned s = 0;
            for (auto & species_name : G::_species_names) {
                species_name_to_index[species_name] = s++;
            }
            
            // Now build G::_taxon_to_species from G::_taxon_names and species_name_to_index
            for (auto & tname : G::_taxon_names) {
                string species_name = Node::taxonNameToSpeciesName(tname);
                unsigned species_index = 0;
                if (species_name_to_index.count(species_name) == 0) {
                    // species_name not found
                    throw XProj(format("Expecting species name \"%s\" to be found in global species names vector but it was not") % species_name);
                }
                else {
                    // species_name found
                    species_index = species_name_to_index[species_name];
                }
                G::_taxon_to_species[tname] = species_index;
            }
        }
        
        for (auto & sname : G::_species_names) {
            unsigned species_index = 0;
            if (species_name_to_index.count(sname) == 0)
                throw XProj(format("Proj::buildSpeciesMap failed because key \"%s\" does not exist in species_name_to_index map") % sname);
            else {
                species_index = species_name_to_index.at(sname);
            }
            
            // Note: despite appearances, this next line does not
            // overwrite anything. We need to be able to map taxon
            // names in species trees to a species index as well as
            // taxon names in gene trees
            G::_taxon_to_species[sname] = species_index;
        }
    }

#if defined (DRAW_NEW_THETA)
    inline void Proj::updateSpeciesNames() {
        unsigned number = G::_nspecies;
        for (int i=0; i<G::_nspecies-1; i++) {
            string name = boost::str(boost::format("node-%d")%number);
            G::_species_names.push_back(name);
            number++;
        }
    }
#endif

#if defined(SPECIES_IN_CONF)
    inline void Proj::parseSpeciesDefinition(string s) {
        // Given these definitions in the conf file:
        //   species = A: a^A, b^A, c^A
        //   species = B: d^B, e^B, f^B
        //   species = C: g^C, h^C, i^C
        //   species = D: j^D, k^D, l^D
        //   species = E: m^E, n^E, o^E
        // This would be the result:
        //   G::_nspecies = 5
        //   G::_species_names = ["A", "B", "C", "D", "E"]
        //   G::_ntaxa = 15
        //   G::_taxon_names = [
        //      "a^A", "b^A", "c^A",
        //      "d^B", "e^B", "f^B",
        //      "g^C", "h^C", "i^C",
        //      "j^D", "k^D", "l^D",
        //      "m^E", "n^E", "o^E"]
        //   G::_taxon_to_species["a^A"] = 0
        //   G::_taxon_to_species["b^A"] = 0
        //   G::_taxon_to_species["c^A"] = 0
        //   G::_taxon_to_species["d^B"] = 1
        //   G::_taxon_to_species["e^B"] = 1
        //   G::_taxon_to_species["f^B"] = 1
        //   G::_taxon_to_species["g^C"] = 2
        //   G::_taxon_to_species["h^C"] = 2
        //   G::_taxon_to_species["i^C"] = 2
        //   G::_taxon_to_species["j^D"] = 3
        //   G::_taxon_to_species["k^D"] = 3
        //   G::_taxon_to_species["l^D"] = 3
        //   G::_taxon_to_species["m^E"] = 4
        //   G::_taxon_to_species["n^E"] = 4
        //   G::_taxon_to_species["o^E"] = 4
        //   G::_taxon_to_species["A"]   = 0
        //   G::_taxon_to_species["B"]   = 1
        //   G::_taxon_to_species["C"]   = 2
        //   G::_taxon_to_species["D"]   = 3
        //   G::_taxon_to_species["E"]   = 4
        
        vector<string> v;
        
        // First separate part before colon (stored in v[0])
        // from the part after colon (stored in v[1])
        split(v, s, boost::is_any_of(":"));
        if (v.size() != 2)
            throw XProj("Expecting exactly one colon in species definition");

        string species_name = v[0];
        string taxon_list = v[1];

        // Now separate the part after the colon at commas to
        // yield the taxa that are in that species
        boost::trim(taxon_list);
        split(v, taxon_list, boost::is_any_of(","));
        for_each(v.begin(), v.end(), [](string & s){boost::trim(s);});
        
        //output(format("Species \"%s\":\n") % species_name, LogCateg::ALWAYS);
        //for (string s : v) {
        //    output(format("   Taxon \"%s\"\n") % s, LogCateg::ALWAYS);
        //}
        
        G::_species_names.push_back(species_name);
        G::_nspecies = (unsigned)G::_species_names.size();
        unsigned i = (unsigned)(G::_nspecies - 1);
        G::_taxon_to_species[species_name] = i;
        for (auto t : v) {
            G::_taxon_names.push_back(t);
            G::_taxon_to_species[t] = i;
            G::_ntaxa++;
        }
    }
#endif

    inline void Proj::run() {
        G::_in_second_level = false;
        
        if (G::_calc_bhv_distances_to_true_tree) {
            assert (G::_newick_path != "");
            string true_spp_tree_file_name = "";
            if (G::_newick_path != ".") {
                true_spp_tree_file_name = G::_newick_path + "/" + "true-species-tree.tre";
            }
            else {
                true_spp_tree_file_name = "true-species-tree.tre";
            }
            
            string true_newick = readNewickFromFile(true_spp_tree_file_name);
        }
        output("Starting serial version...\n");
        if (G::_gene_newicks_specified) {
            if (G::_start_mode_type == G::StartModeType::START_MODE_SIM) {
                throw XProj("cannot specify gene newicks and simulations");
            }
            try {
//                if (G::_thin < 1.0) {
//                    cout << "thin setting will be set to 1.0 for gene newick start " << endl;
//                    G::_thin = 1.0;
//                }
                handleGeneNewicks();
            }
            catch (XProj & x) {
                std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
            }
        }
        
        else if (G::_start_mode_type == G::StartModeType::START_MODE_SIM) {
            if (G::_gene_newicks_specified) {
                throw XProj("cannot specify gene newicks and simulations");
            }
            
            _first_line = true;
            
            try {
                simulateData();
            }
            catch (XProj & x) {
                std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
            }
        }
        else {
            _first_line = true;
            if (G::_verbose > 0) {
                cout << "Starting..." << endl;
                cout << "Current working directory: " << filesystem::current_path() << endl;
                cout << "Random seed: " << _random_seed << endl;
#if defined (DRAW_NEW_THETA)
                cout << "drawing new theta for each particle " << endl;
#else
                cout << "Theta: " << G::_theta << endl;
#endif
                cout << "Number of threads: " << G::_nthreads << endl;
            }

            if (G::_run_on_empty) { // if running with no data, choose taxa to join at random
                G::_proposal = "prior-prior";
                G::_prior_prior = true;
            }

            try {
                if (G::_verbose > 0) {
                    cout << "\n*** Reading and storing the data in the file " << G::_data_file_name << endl;
                    cout << "data file name is " << G::_data_file_name << endl;
                }
                
                _data = Data::SharedPtr(new Data());
                _data->setPartition(_partition);
                _data->getDataFromFile(G::_data_file_name);

                if (G::_verbose > 0) {
                    summarizeData(_data);
                }
                
#if defined (SPECIES_IN_CONF)
//                createSpeciesMap(false);
#else
                createSpeciesMap(true);
#endif
                
                G::_nspecies = (unsigned) G::_species_names.size();
                createSpeciesMapTyped();

                // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
                if (G::_outgroup != "none") {
                    checkOutgroupName();
                }

#if defined(SPECIES_IN_CONF)
            if (G::_nspecies > 0) {
                // Species specified in the conf file
                // Check that taxon names are the same as those
                // in the data file
                _data->checkTaxonNames(G::_taxon_names);
                
                // set some global variables
                G::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(G::_taxon_names);
                assert (G::_species_names.size() > 0);
                
                assert (G::_ntaxa > 0);
                assert (G::_nspecies > 0);
            }
            else {
                // Copy taxon names to global variable _taxon_names
                G::_ntaxa = _data->getNumTaxa();
                _data->copyTaxonNames(G::_taxon_names);
                
                G::_nspecies = (unsigned) G::_species_names.size();
                assert (G::_nspecies > 0);
            }
#else
            // set some global variables
            G::_ntaxa = _data->getNumTaxa();
            assert (G::_species_names.size() > 0);
#endif
            G::_nloci = _data->getNumSubsets();
                
            assert (G::_nloci > 0);
            ps.setNLoci(G::_nloci);
            for (unsigned locus = 1; locus < G::_nloci+1; locus++) {
                // Set length of partials for gene g
                unsigned ncat = 1;
                if (G::_plus_G) {
                    ncat = 4; // TODO: don't hard code this
                }
                ps.setNElements(G::_nstates*_data->getNumPatternsInSubset(locus-1) * ncat, locus);
            }
                
            // set random number seed
            rng.setSeed(_random_seed);

            if (!G::_gene_newicks_specified) {
#if defined (SPECIES_IN_CONF)
//                    buildSpeciesMap(/*taxa_from_data*/false);

#else
                buildSpeciesMap(/*taxa_from_data*/true);
#endif
            }
                
            Particle p;
            initializeParticle(p); // initialize one particle and copy to all other particles
            
            // reset marginal likelihood
            _log_marginal_likelihood = 0.0;
            
            _starting_log_likelihoods = p.calcGeneTreeLogLikelihoods();
            
            for (unsigned g=0; g<G::_ngroups; g++) { // TODO: unsure if this is correct
                for (unsigned i=0; i<G::_nloci; i++) {
                    _log_marginal_likelihood += _starting_log_likelihoods[i];
                }
            }
            
            vector<Particle> my_vec;
            my_vec.resize(G::_nparticles * G::_ngroups, p);
                
            if (G::_sample_from_prior) {
                for (auto &p:my_vec) {
                    unsigned psuffix = 1;
                        p.setSeed(rng.randint(1,9999) + psuffix);
                        psuffix += 2;
                    }
                    for (auto &p:my_vec) {
                        p.buildEntireSpeciesTree();
                    }
                saveSpeciesTreesAfterFirstRound(my_vec);
                
            }

            else {
                unsigned psuffix = 1;
                for (auto &p:my_vec) {
                    p.setSeed(rng.randint(1,9999) + psuffix);
                    psuffix += 2;
                }
                
                if (G::_species_newick_specified) {
                    string species_newick = handleSpeciesNewick();
                    unsigned count = 0;
                    for (auto &p:my_vec) {
                        p.processSpeciesNewick(species_newick);
                        count++;
                    }
                }
                
                // set group rng
                _group_rng.resize(G::_ngroups);
                psuffix = 1;
                for (auto &g:_group_rng) {
                    g.reset(new Lot());
                    g->setSeed(rng.randint(1,9999)+psuffix);
                    psuffix += 2;
                }

    #if defined (DRAW_NEW_THETA)
                updateSpeciesNames();
                for (auto &p:my_vec) {
                    p.drawTheta();
                }
    #endif
                    
                // particle_indices holds the subgroup each particle is in
                unsigned total_n_particles = G::_nparticles * G::_ngroups;
                vector<unsigned> particle_indices(total_n_particles);

                // fill particle_indices with values starting from 0
                iota(particle_indices.begin(), particle_indices.end(), 0);
                    
                // if using subgroups, reset pointers so particles within a group have same gene forest pointers
                if (G::_ngroups > 1) {
                    for (unsigned i=0; i<G::_ngroups; i++) {
                        unsigned start = i * G::_nparticles;
                        unsigned end = start + (G::_nparticles) - 1;
                        vector<Forest::SharedPtr> gfcpies;
                        for (unsigned l=0; l<G::_nloci; l++) {
                            gfcpies.push_back(Forest::SharedPtr(new Forest()));
                        }
                        for (unsigned p=start; p<end+1; p++) {
                            my_vec[p].resetSubgroupPointers(gfcpies);
                        }
                    }
                }
                    //run through each generation of particles

                    unsigned nsteps = (G::_ntaxa-1)*G::_nloci;
                    
                if (G::_verbose > 1) {
                    cout << "step " << "\t" << "ESS before filtering " << "\t" << "n_unique particles before MCMC" <<  "\t" << "n_unique particles after MCMC" << endl;
                }
                    
                for (unsigned g=0; g<nsteps; g++) {
                    if (g == 0) {
                        // reset gene order
                        unsigned list_size = G::_nloci;
                        vector<vector<unsigned>> new_gene_order;
                        new_gene_order.resize(G::_ngroups);
                        for (auto n:new_gene_order) {
                            n.resize(G::_nloci);
                        }

                        for (unsigned n=0; n<G::_ngroups; n++) {
                            unsigned count = 1;
                            vector<pair<double, unsigned>> randomize;
                            for (unsigned l=0; l<list_size; l++) {
                                if (count == 1) {
                                    randomize.clear();
                                }
                                randomize.push_back(make_pair(rng.uniform(), count));
                                count++;
                                if (count > G::_nloci) {
                                    sort(randomize.begin(), randomize.end());
                                    for (auto &r:randomize) {
                                        new_gene_order[n].push_back(r.second);
                                    }
                                    count = 1;
                                }
                            }
                        }

                        // set gene order for first G::_nloci set of steps
                        unsigned ngroup = 0;
                        unsigned group_count = 0;
                        for (unsigned p=0; p<G::_nparticles*G::_ngroups; p++) {
                            my_vec[p].resetGeneOrder(g, new_gene_order[ngroup]);
                            if ((group_count+1)%G::_nparticles == 0) {
                                ngroup++;
                            }
                            group_count++;
                        }
                    }
                    if (G::_verbose == 1) {
                        cout << "starting step " << g << " of " << nsteps-1 << endl;
                    }

                    unsigned psuffix = 1;
                    if (g > 0) {
                        // set particle random number seeds
                        for (auto &p:my_vec) {
                            p.setSeed(rng.randint(1,9999) + psuffix);
                            psuffix += 2;
                        }
                    }
                    
                    //taxon joining and reweighting step
                    proposeParticles(my_vec);
                    bool filter = true;
                    
                    if (filter) {
                            
                        // TODO: can parallelize filtering by subgroup
                            if (G::_nthreads == 1) {
                            for (unsigned i=0; i<G::_ngroups; i++) {
                                unsigned start = i * G::_nparticles;
                                unsigned end = start + (G::_nparticles) - 1;

                                double ess = -1;
                                ess = filterParticles(g, my_vec, particle_indices, start, end);
                                
                                vector<double> weights_after_filtering(G::_nparticles);
                                
                                for (unsigned p=0; p<G::_nparticles; p++) {
                                    weights_after_filtering[p] = my_vec[p].getLogWeight();
                                }
                                
                                std::sort(weights_after_filtering.begin(), weights_after_filtering.end());
                                
                                double n_unique_particles_after_filtering = std::unique(weights_after_filtering.begin(), weights_after_filtering.end()) - weights_after_filtering.begin();
                                
                                
                                if (G::_verbose > 1) {
                                    cout << G::_generation  << "\t" << ess << "\t" << "\t" << "\t" << n_unique_particles_after_filtering << "\t";
                                    if (!G::_mcmc) {
                                        cout << endl;
                                    }
                                }
                            }
                        }
                        else {
                            filterParticlesThreading(my_vec, g, particle_indices);
                        }
                        
                        string filenamea = "params" + to_string(G::_generation) + "a";
                        
                        if (G::_generation == 0) {
                            string filename = "mcmc_moves_accepted.log";
                            if (filesystem::remove(filename)) {
                                ofstream mcmcfile(filename);
                                mcmcfile << "generation" << "\t" << "number of mcmc moves accepted" << "\t" << "proportion of mcmc moves accepted" << "\t" << "average log likelihood before mcmc" << "\t" << "average log likelihood after mcmc" << "\t" << "sliding window" << "\n";
                            }
                            else {
                                ofstream mcmcfile(filename);
                                mcmcfile << "generation" << "\t" << "number of mcmc moves accepted" << "\t" << "proportion of mcmc moves accepted" << "\t" << "average log likelihood before mcmc" << "\t" << "average log likelihood after mcmc" << "sliding window" <<"\n";
                            }
                        }
                                 
                        unsigned locus = my_vec[0].getNextGene() - 1;
                        
                        vector<double> log_likelihoods_before_mcmc;
                        for (auto &p:my_vec) {
                            log_likelihoods_before_mcmc.push_back(p.calcLogLikelihoodLocus(locus, false));
                        }
                        
                        double sum_log_likelihood_before_mcmc = std::accumulate(log_likelihoods_before_mcmc.begin(), log_likelihoods_before_mcmc.end(), 0);
                        double avg_log_likelihood_before_mcmc = sum_log_likelihood_before_mcmc / G::_nparticles;
                        
                        if (G::_mcmc) {
                            G::_nmcmc_moves_accepted = 0;
                            for (unsigned m=0; m<G::_n_mcmc_rounds; m++) {
                                bool last_round = false;
                                if (m == G::_n_mcmc_rounds - 1) {
                                    last_round = true;
                                }
                                mcmcMoves(my_vec, last_round);
                            }
                            
                            // TODO: no groups
                            unsigned locus = my_vec[0].getNextGene() - 1; // subtract 1 because vector of gene forests starts at 0
                            for (unsigned p=0; p<G::_nparticles; p++) {
                                // finalize join for every particle now
                                my_vec[p].finalizeLatestJoinMCMC(locus, p);
                            }

                            vector<double> log_likelihoods_after_mcmc;
                            for (auto &p:my_vec) {
                                log_likelihoods_after_mcmc.push_back(p.calcLogLikelihoodLocus(locus, true));
                            }
                            
                            double sum_log_likelihood_after_mcmc = std::accumulate(log_likelihoods_after_mcmc.begin(), log_likelihoods_after_mcmc.end(), 0);
                            double avg_log_likelihood_after_mcmc = sum_log_likelihood_after_mcmc / G::_nparticles;
                            
                            std::ofstream mcmcfile;

                            mcmcfile.open("mcmc_moves_accepted.log", std::ios_base::app); // append instead of overwrite
                            double proportion_accepted = (double) G::_nmcmc_moves_accepted / (G::_nparticles * G::_n_mcmc_rounds);
                            mcmcfile << G::_generation << "\t" << "\t" << G::_nmcmc_moves_accepted << "\t" << "\t" << proportion_accepted << "\t" << "\t" << "\t" << avg_log_likelihood_before_mcmc << "\t" << "\t" << "\t" << avg_log_likelihood_after_mcmc << "\t" << "\t" << "\t" << G::_sliding_window << "\n";
                            
                            vector<double> weights_after_mcmc(G::_nparticles);
                            
                            for (unsigned p=0; p<G::_nparticles; p++) {
                                weights_after_mcmc[p] = my_vec[p].getLogWeight();
                            }
                            
                            std::sort(weights_after_mcmc.begin(), weights_after_mcmc.end());
                            double n_unique_particles_after_mcmc = std::unique(weights_after_mcmc.begin(), weights_after_mcmc.end()) - weights_after_mcmc.begin();
                                    
                            if (G::_verbose > 1) {
                                cout << "\t" << "\t" << "\t" << n_unique_particles_after_mcmc << "\n";
                            }
                            
                            string filenameb = "params" + to_string(G::_generation) + "b";
                        }
                        else {
                            std::ofstream mcmcfile;

                            mcmcfile.open("mcmc_moves_accepted.log", std::ios_base::app); // append instead of overwrite

                            
                            mcmcfile << G::_generation << "\t" << "\t" << "N/A" << "\t" << "N/A" << "\t" << "\t" << "\t" << avg_log_likelihood_before_mcmc << "\t" << "\t" << "\t" << "N/A" << "\t" << "N/A" << "\n";
                        }
                    }
                                            
                        // only shuffle groups after all particles have coalesced each locus once, then reset gene order
                        if ((g+1)%G::_nloci == 0 && g != nsteps-1) {
                                // reset gene order
                            unsigned list_size = G::_nloci;
                            vector<vector<unsigned>> new_gene_order;
                            new_gene_order.resize(G::_ngroups);
                            for (auto n:new_gene_order) {
                                n.resize(G::_nloci);
                            }

                            for (unsigned n=0; n<G::_ngroups; n++) {
                                unsigned count = 1;
                                vector<pair<double, unsigned>> randomize;
                                for (unsigned l=0; l<list_size; l++) {
                                    if (count == 1) {
                                        randomize.clear();
                                    }
                                    randomize.push_back(make_pair(rng.uniform(), count));
                                    count++;
                                    if (count > G::_nloci) {
                                        sort(randomize.begin(), randomize.end());
                                        for (auto &r:randomize) {
                                            new_gene_order[n].push_back(r.second);
                                        }
                                        count = 1;
                                    }
                                }
                            }
                            
                            unsigned ngroup = 0;
                            unsigned group_count = 0;
                            for (unsigned p=0; p<G::_nparticles*G::_ngroups; p++) {
                                unsigned particle_number = particle_indices[p];
                                my_vec[particle_number].resetGeneOrder(g, new_gene_order[ngroup]);
                                if ((group_count+1)%G::_nparticles == 0) {
                                    ngroup++;
                                }
                                group_count++;
                            }
                        }
                    G::_generation++;
                    
    #if defined (VALGRIND)
                    VALGRIND_PRINTF("~~> post 1st-level step %d at time %d\n", g, (unsigned)clock()); // g = step
                    VALGRIND_MONITOR_COMMAND(str(format("detailed_snapshot stepsnaps-%d.txt") % g).c_str());
    #endif
                } // g loop
                    
                for (auto &p:my_vec) {
                    p.calcGeneTreeLengths();
                }
                    
                    if (G::_save_gene_trees) {
                        saveGeneTrees(my_vec);
                    }
                    if (G::_verbose > 0) {
                        cout << "\n";
                        cout << "marginal likelihood after combined filtering: " << _log_marginal_likelihood << endl;
                        cout << "\n";
                    }
                    
                     writePartialCountFile(my_vec);
                    
                    for (auto &p:my_vec) {
                        p.setSortedThetaVector();
                    }
                    
                    if (!G::_second_level) {
                        // remove 1 - thin % of particles
                        double n_elements_to_remove = (1 - G::_thin) * my_vec.size();
                        double n_elements_to_keep = my_vec.size() - n_elements_to_remove;
                        
                        // shuffle my_vec
                        shuffle(my_vec.begin(), my_vec.end(), std::default_random_engine(_random_seed));
                        
                        // keep only n_elements_to_keep
                        my_vec.resize(n_elements_to_keep);

                        
                        for (unsigned i=0; i<G::_nloci; i++) {
                            ofstream heightf;
                            string filename = "gene_tree_heights" + to_string(i+1) + ".txt";
                            if (filesystem::remove(filename)) {
                                ofstream heightf(filename);
                                if (G::_verbose > 0) {
                                   cout << "existing file " << filename << " removed and replaced\n";
                                }
                            }
                            else {
                                ofstream heightf(filename);
                                if (G::_verbose > 0) {
                                    cout << "created new file " << filename << "\n";
                                }
                            }
                        }
                        
                        writeParamsFileForBeastComparison(my_vec);
                    
                        if (G::_write_species_tree_file) {
                            saveAllSpeciesTrees(my_vec);
                        }
                        
                        if (G::_ruv) { // validation for gene trees only
                            assert (_gene_tree_ranks.size() == G::_nloci);
                            
                            for (unsigned l=0; l<G::_nloci; l++) {
                                string sim_file_name;
                                assert(G::_newick_path != "");
                                sim_file_name = G::_newick_path + "/" + "true_gene_tree_height" + to_string(l+1) + ".txt";
                                string line;
                                string height_as_string;
                                ifstream infile(sim_file_name);
                                while (getline(infile, line)) {
                                    height_as_string = line;
                                }
                                double true_gene_tree_height = 0.0;
                                assert (height_as_string != "");
                                true_gene_tree_height = std::stod(height_as_string);
                                _gene_tree_ranks[l].push_back(make_pair(true_gene_tree_height, true));
                                // sort ranks
                                std::sort(_gene_tree_ranks[l].begin(), _gene_tree_ranks[l].end());
                                
                                // find rank of truth
                                auto it = std::find_if(_gene_tree_ranks[l].begin(), _gene_tree_ranks[l].end(), [&](const pair<double, bool>& p) { return p.second == true;});
                                unsigned index_value = (unsigned) std::distance(_gene_tree_ranks[l].begin(), it);
                                
                                // write rank value to file
                                ofstream rankf("rank_gene_tree_height" + to_string(l+1) + ".txt");
                                rankf << "rank: " << index_value << endl;
                            }
                        }
                        if (G::_ruv_first_level_species) {
                            assert (_species_tree_ranks_after_first_round.size() > 0);
                            string sim_file_name;
                            assert(G::_newick_path != "");
                            sim_file_name = G::_newick_path + "/" + "true_species_tree_height.txt";
                            string line;
                            string height_as_string;
                            ifstream infile(sim_file_name);
                            while (getline(infile, line)) {
                                height_as_string = line;
                            }
                            double true_species_tree_height = 0.0;
                            assert (height_as_string != "");
                            true_species_tree_height = std::stod(height_as_string);
                            _species_tree_ranks_after_first_round.push_back(make_pair(true_species_tree_height, true));
                            // sort ranks
                            std::sort(_species_tree_ranks_after_first_round.begin(), _species_tree_ranks_after_first_round.end());
                            
                            // find rank of truth
                            auto it = std::find_if(_species_tree_ranks_after_first_round.begin(), _species_tree_ranks_after_first_round.end(), [&](const pair<double, bool>& p) { return p.second == true;});
                            unsigned index_value = (unsigned) std::distance(_species_tree_ranks_after_first_round.begin(), it);
                            
                            // write rank value to file
                            ofstream rankf("rank_species_tree_height_after_first_round.txt");
                            rankf << "rank: " << index_value << endl;
                        }
                        if (G::_bhv_reference != "" || G::_bhv_reference_path != ".") {
                            assert (_bhv_distances_genes.size() > 0);
                            for (unsigned l=0; l<G::_nloci; l++) {
                                string sim_file_name;
                                assert (G::_newick_path != "");
                                if (G::_newick_path == ".") {
                                    sim_file_name = "gene" + to_string(l+1) + ".trees";
                                }
                                else {
                                    sim_file_name = G::_newick_path + "/" + "gene" + to_string(l+1) + ".trees";
                                }
                                 
                                 Particle p;
                                 string true_newick = readNewickFromFile(sim_file_name);
                                 
                                 double true_bhv = p.calcBHVDistanceTrueTreeGene(true_newick);
                                 
                                _bhv_distances_genes[l].push_back(true_bhv);
    //                             // sort distances
                                 std::sort(_bhv_distances_genes[l].begin(), _bhv_distances_genes[l].end());
    
    //                             // find rank of truth
                                 auto it = std::find(_bhv_distances_genes[l].begin(), _bhv_distances_genes[l].end(), true_bhv);
                                 unsigned index_value = (unsigned) std::distance(_bhv_distances_genes[l].begin(), it);

                                 // write rank value to file
                                 ofstream rankf("rank_bhv" + to_string(l+1) + ".txt");
                                 rankf << "rank: " << index_value << endl;
                             }
                        }
                        
                        if (G::_hpd) {
                            assert (_hpd_values_genes.size() > 0);
                            for (unsigned l=0; l<G::_nloci; l++) {
                                 ofstream hpdf("hpd" + to_string(l+1) + ".txt");
                                 hpdf << "min    " << "max " << endl;
                                 
                                 // sort hpd values largest to smallest
                                 std::sort(_hpd_values_genes[l].begin(), _hpd_values_genes[l].end());
                                 std::reverse(_hpd_values_genes[l].begin(), _hpd_values_genes[l].end());
                                 
                                 // take first 95% of values (round down to nearest integer)
                                 double total = size(_hpd_values_genes[l]);
                                 double ninety_five_index = floor(0.95*total);
                                 
                                 if (ninety_five_index == 0) {
                                     ninety_five_index = 1;
                                 }
                                 
                                 vector<double> hpd_values_in_range;
                                 
                                 for (unsigned h=0; h<ninety_five_index; h++) {
                                     hpd_values_in_range.push_back(_hpd_values_genes[l][h].second);
                                 }
                                 
                                 auto max = *std::max_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
                                 auto min = *std::min_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
                                 assert (min < max || min == max);
                                 
                                 // write min and max to file
                                 hpdf << min << "\t" << max << endl;
                            }
                         }
                        
                        if (G::_hpd_first_level_species) {
                            ofstream hpdf("hpd_first_level_species.txt");
                            hpdf << "min    " << "max " << endl;
                            
                            // sort hpd values largest to smallest
                            std::sort(_hpd_first_level_values.begin(), _hpd_first_level_values.end());
                            std::reverse(_hpd_first_level_values.begin(), _hpd_first_level_values.end());
                            
                            // take first 95% of values (round down to nearest integer)
                            double total = size(_hpd_first_level_values);
                            double ninety_five_index = floor(0.95*total);
                            
                            if (ninety_five_index == 0) {
                                ninety_five_index = 1;
                            }
                            
                            vector<double> hpd_values_in_range;
                            
                            for (unsigned h=0; h<ninety_five_index; h++) {
                                hpd_values_in_range.push_back(_hpd_first_level_values[h].second);
                            }
                            
                            auto max = *std::max_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
                            auto min = *std::min_element(hpd_values_in_range.begin(), hpd_values_in_range.end());
                            assert (min < max || min == max);
                            
                            // write min and max to file
                            hpdf << min << "\t" << max << endl;
                        }
                        
                        if (G::_hpd || G::_ruv) {
                            // write mean gene tree heights to output file for validation
                            for (unsigned l=0; l<G::_nloci; l++) {
                                double sum = accumulate(_gene_tree_heights[l].begin(), _gene_tree_heights[l].end(), 0.0);
                                double mean = sum / _gene_tree_heights[l].size();
                                
                                ofstream heightf("average_gene_tree_height" + to_string(l+1) + ".txt");
                                heightf << mean << endl;
                            }
                        }
                        if (G::_hpd_first_level_species || G::_ruv_first_level_species) {
                            // write mean gene tree heights to output file for validation
                            double sum = accumulate(_species_tree_heights_after_first_round.begin(), _species_tree_heights_after_first_round.end(), 0.0);
                            double mean = sum / _species_tree_heights_after_first_round.size();

                            ofstream heightf("average_species_tree_height_after_first_round.txt");
                            heightf << mean << endl;
                        }
                    }

                    else {
    #if defined (HIERARCHICAL_FILTERING)
                        saveSpeciesTreesAfterFirstRound(my_vec);
                        
                        G::_in_second_level = true;
                        
                        cout << "\n";
                        string filename1 = "species_trees.trees";
                        string filename2 = "unique_species_trees.trees";
                        string filename3 = "params-beast-comparison.log";
                        if (filesystem::remove(filename1)) {
                            ofstream speciestrf(filename1);
                            speciestrf << "#nexus\n\n";
                            speciestrf << "begin trees;\n";
                            if (G::_verbose > 0) {
                                cout << "existing file " << filename1 << " removed and replaced\n";
                            }
                        }
                        else {
                            ofstream speciestrf(filename1);
                            speciestrf << "#nexus\n\n";
                            speciestrf << "begin trees;\n";
                            if (G::_verbose > 0) {
                                cout << "created new file " << filename1 << "\n";
                            }
                        }
                        if (filesystem::remove(filename2)) {
                            ofstream uniquespeciestrf(filename2);
                            uniquespeciestrf << "#nexus\n\n";
                            uniquespeciestrf << "begin trees;\n";
                            if (G::_verbose > 0) {
                                cout << "existing file " << filename2 << " removed and replaced\n";
                            }
                        }
                        else {
                            ofstream uniquespeciestrf(filename2);
                            uniquespeciestrf << "#nexus\n\n";
                            uniquespeciestrf << "begin trees;\n";
                            if (G::_verbose > 0) {
                                cout << "created new file " << filename2 << "\n";
                            }
                        }
                        
                        if (filesystem::remove(filename3)) {
                            ofstream paramsf(filename3);
                            if (G::_verbose > 0) {
                               cout << "existing file " << filename3 << " removed and replaced\n";
                            }
                        }
                        else {
                            ofstream paramsf(filename3);
                            if (G::_verbose > 0) {
                                cout << "created new file " << filename3 << "\n";
                            }
                        }
                        cout << "\n";

                        unsigned ngroups = round(G::_nparticles * G::_ngroups * G::_thin);
                        if (ngroups == 0) {
                            ngroups = 1;
                            cout << "thin setting would result in 0 species groups; setting species groups to 1" << endl;
                        }
                        
                        ofstream heightf;
                        if (filesystem::remove("species_heights.txt")) {
                            ofstream heightf("species_heights.txt");
                            if (G::_verbose > 0) {
                               cout << "existing file " << "species_heights.txt" << " removed and replaced\n";
                            }
                        }
                        else {
                            ofstream heightf("species_heights.txt");
                            if (G::_verbose > 0) {
                                cout << "created new file " << "species_heights.txt" << "\n";
                            }
                        }
                        
                        string altfname = "alt_species_trees.trees";
                        if (filesystem::remove(altfname)) {
                            
                            if (G::_verbose > 0) {
                               cout << "existing file " << altfname << " removed and replaced\n";
                            }
                            
                            ofstream alttrf(altfname);

                            alttrf << "#nexus\n\n";
                            alttrf << "Begin trees;" << endl;
                            string translate_block = my_vec[0].getTranslateBlock();
                            alttrf << translate_block << endl;
                        }
                        else {
                            ofstream alttrf(altfname);

                            alttrf << "#nexus\n\n";
                            alttrf << "Begin trees;" << endl;
                            string translate_block = my_vec[0].getTranslateBlock();
                            alttrf << translate_block << endl;
                            
                            if (G::_verbose > 0) {
                                cout << "created new file " << altfname << "\n";
                            }
                        }
                        cout << "\n";
                        
                        
                        heightf << "species_tree_height" << endl;
                        
                        secondLevel(my_vec);
    #endif
                    }
                }
            }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }
        }

        std::cout << "\nFinished!" << std::endl;
    }
}
