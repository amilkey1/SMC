#pragma once

#include <iostream>
#include "data.hpp"
#include "partition.hpp"
#include "stopwatch.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
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

using namespace std;
using namespace boost;
using namespace boost::algorithm;

#include "partial_store.hpp"
#include "g.hpp"
extern proj::PartialStore ps;
extern proj::Lot rng;
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
            void                saveSpeciesTreesHierarchical(vector<Particle> &v, string filename1, string filename2) const;
            void                saveGeneTrees(unsigned ngenes, vector<Particle> &v) const;
            void                writeLoradFile(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v) const;
            void                writeLoradFileAfterSpeciesFiltering(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v) const;
            void                writeDeepCoalescenceFile(vector<Particle> &v);
            void                writeThetaFile(vector<Particle> &v);
            void                writeParamsFileForBeastComparison (unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v) const;
            void                writeParamsFileForBeastComparisonAfterSpeciesFiltering (unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v, string filename, unsigned group_number);
            void                writeParamsFileForBeastComparisonAfterSpeciesFilteringSpeciesOnly(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v, string filename, unsigned group_number);
            void                writePartialCountFile(vector<Particle> &particles);
            void                createSpeciesMap(Data::SharedPtr);
            void                simSpeciesMap();
            string              inventName(unsigned k, bool lower_case);
            void                showFinal(vector<Particle> my_vec);
            void                proposeSpeciesParticleRange(unsigned first, unsigned last, vector<Particle> &particles);
            void                proposeSpeciesParticles(vector<Particle> &particles);
            void                proposeSpeciesGroups(vector<Particle> &particles, unsigned ngroups, string filename1, string filename2, string filename3, unsigned nsubsets, unsigned ntaxa);
            void                proposeSpeciesGroupRange(unsigned first, unsigned last, vector<Particle> &particles, unsigned ngroups, string filename1, string filename2, string filename3, unsigned nsubsets, unsigned ntaxa);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle> &particles);
            void                proposeParticleGroupRange(unsigned first, unsigned last, vector<vector<Particle>> &particles);
            void                proposeParticles(vector<Particle> &particles);
            void                proposeParticlesParallelizeByGroup(vector<vector<Particle>> &particles);
            void                simulateData();
            void                writePaupFile(vector<Particle> particles, vector<string> taxpartition);
            void                initializeParticles(vector<Particle> &particles);
            void                initializeParticleRange(unsigned first, unsigned last, vector<Particle> &particles);
            void                handleGeneNewicks();
            void                handleSpeciesNewick(vector<Particle> particles);
            double              filterParticles(unsigned step, vector<Particle> & particles, unsigned start, unsigned end);


            void                filterParticlesMixing(vector<Particle> & particles, unsigned step);
            unsigned            multinomialDraw(const vector<double> & probs);
            double              filterSpeciesParticles(unsigned step, vector<Particle> & particles);
            double              computeEffectiveSampleSize(const vector<double> & probs) const;
            void                saveSpeciesTreesAltHierarchical(vector<Particle> &v) const;


        private:

            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood = 0.0;
            double                      _log_species_tree_marginal_likelihood = 0.0;
            unsigned                    _random_seed;
            void                        summarizeData(Data::SharedPtr);
            unsigned                    setNumberTaxa(Data::SharedPtr);
            double                      getRunningSum(const vector<double> &) const;
            vector<string>              _species_names;
            map<string, string>         _taxon_map;
            void                        handleBaseFrequencies();
            void                        handleRelativeRates();
            void                        handleNTaxaPerSpecies();
            void                        checkOutgroupName();
            void                        debugSpeciesTree(vector<Particle> &particles);
            double                      _small_enough;
            bool                        _first_line;
            unsigned                    _count; // counter for params output file
            double                      _starting_log_likelihood;
            vector<Lot::SharedPtr>      _group_rng;

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
        partialf << "total times partials calculated: ";
        
        unsigned partial_count = 0;
            for (auto &p:particles) {
                partial_count += p.getPartialCount();
            }
        partialf << partial_count << "\n";
        partialf.close();
    }

    inline void Proj::writeParamsFileForBeastComparison(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v) const {
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

        for (int i=1; i<ngenes+1; i++) {
            logf << "\t" << "Tree.t:gene" + to_string(i) + "height";
            logf << "\t" << "Tree.t:gene" + to_string(i) + "treeLength";
        }

        logf << "\t" << "YuleModel.t:Species "; // this is the log probability of the species tree (multiply by log(3!) to get increment log prob)
        logf << "\t" << "popMean "; // this is psi in the InverseGamma(2,psi) distribution of popSize

        for (int i=0; i<(nspecies*2-1); i++) {
            logf << "\t" << "popSize." + to_string(i+1);
        }

        logf << "\t" << "speciationRate.t:Species ";

        for (int i=1; i<ngenes+1; i++) {
            logf << "\t" << "treeLikelihood:gene" + to_string(i);
        }
        for (int i=1; i<ngenes+1; i++) {
            logf << "\t" << "treePrior:gene" + to_string(i);
        }
        logf << endl;

        int iter = 0;
        
        for (auto &p:v) {
                logf << iter;
                iter++;

                double vector_prior = 0.0;
    #if defined DRAW_NEW_THETA
                vector<double> vector_priors = p.getVectorPrior();
                for (auto &v:vector_priors) {
                    vector_prior += v; // this is the InverseGamma(2, psi) prior on the 5 population sizes -- only for first round
                }
    #endif

                double log_coalescent_likelihood = 0.0;
                for (unsigned g=1; g<ngenes+1; g++) {
                    log_coalescent_likelihood += p.getCoalescentLikelihood(g);
                }

                double log_likelihood = p.getLogLikelihood();
                double log_prior = p.getAllPriorsFirstRound();

                double log_posterior = log_likelihood + log_prior + log_coalescent_likelihood + vector_prior;
                // no vector prior under Jones method

                logf << "\t" << log_posterior;

                logf << "\t" << log_likelihood;

                logf << "\t" << log_prior - log_coalescent_likelihood; // starbeast3 does not include coalescent likelihood in this prior


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

                for (int i=0; i<(nspecies*2-1); i++) {
    #if defined (DRAW_NEW_THETA)
                    vector<double> theta_vec = p.getThetaVector();
                    logf << "\t" << theta_vec[i] / 4.0;
    #else
                    logf << "\t" << Forest::_theta / 4.0; // all pop sizes are the same under this model, Ne*u = theta / 4?
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

    inline void Proj::writeParamsFileForBeastComparisonAfterSpeciesFiltering(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v, string filename, unsigned group_number) {
        // this function creates a params file that is comparable to output from starbeast3
        std::ofstream logf;

        logf.open(filename, std::ios_base::app);

        if (_first_line) {
            _first_line = false;
#if !defined PARALLELIZE_BY_GROUP
            logf << "iter ";
#endif
            logf << "\t" << "posterior ";
            logf << "\t" << "likelihood ";
            logf << "\t" << "prior ";
            logf << "\t" << "vectorPrior "; // log joint prior on population sizes (vectorPrior)
            logf << "\t" << "speciescoalescent ";
            logf << "\t" << "Tree.t:Species.height ";
            logf << "\t" << "Tree.t:Species.treeLength ";

            for (int i=1; i<ngenes+1; i++) {
                logf << "\t" << "Tree.t:gene" + to_string(i) + "height";
                logf << "\t" << "Tree.t:gene" + to_string(i) + "treeLength";
            }

            logf << "\t" << "YuleModel.t:Species "; // this is the log probability of the species tree (multiply by log(3!) to get increment log prob)
            logf << "\t" << "popMean "; // this is psi in the InverseGamma(2,psi) distribution of popSize

            for (int i=0; i<(nspecies*2-1); i++) {
                logf << "\t" << "popSize." + to_string(i+1);
            }

            logf << "\t" << "speciationRate.t:Species ";

            for (int i=1; i<ngenes+1; i++) {
                logf << "\t" << "treeLikelihood:gene" + to_string(i);
            }
            for (int i=1; i<ngenes+1; i++) {
                logf << "\t" << "treePrior:gene" + to_string(i);
            }
            logf << endl;
        }

        unsigned sample_size = round(double (G::_particle_increase) / double(G::_save_every) );
        if (sample_size == 0) {
            sample_size = G::_particle_increase;
        }
        
        for (auto &p:v) {
#if !defined PARALLELIZE_BY_GROUP
            unsigned iter = group_number * sample_size;
            logf << iter;
            iter++;
#endif

            double log_coalescent_likelihood = 0.0;
#if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
            log_coalescent_likelihood += p.getCoalescentLikelihood(1);
#else
            for (unsigned g=1; g<ngenes+1; g++) {
                log_coalescent_likelihood += p.getCoalescentLikelihood(g);
            }
#endif

            double vector_prior = 0.0;

            double log_likelihood = p.getLogLikelihood();
            double log_prior = p.getAllPriors();

            double log_posterior = log_likelihood + log_prior + log_coalescent_likelihood + vector_prior;

            logf << "\t" << log_posterior;

            logf << "\t" << log_likelihood;

            logf << "\t" << log_prior - log_coalescent_likelihood; // starbeast3 does not include coalescent likelihood in this prior

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

            for (int i=0; i<(nspecies*2-1); i++) {
#if defined (DRAW_NEW_THETA)
                if (!G::_gene_newicks_specified) {
                    vector<double> theta_vec = p.getThetaVector();
                    logf << "\t" << theta_vec[i] / 4.0;
                }
#else
                logf << "\t" << Forest::_theta / 4.0; // all pop sizes are the same under this model, Ne*u = theta / 4?
#endif
            }

            logf << "\t" << G::_lambda; // TODO: not estimating lambda for now

            vector<double> gene_tree_log_likelihoods = p.getGeneTreeLogLikelihoods();
            vector<double> gene_tree_priors = p.getGeneTreeCoalescentLikelihoods();

#if !defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
            double test = 0.0;
            for (auto &p:gene_tree_priors) {
                test += p;
            }
            assert (test == log_coalescent_likelihood);
            assert (gene_tree_log_likelihoods.size() == gene_tree_priors.size());
#endif

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

    inline void Proj::writeParamsFileForBeastComparisonAfterSpeciesFilteringSpeciesOnly(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v, string filename, unsigned group_number) {
        // this function creates a params file that is comparable to output from starbeast3
        std::ofstream logf;

        logf.open(filename, std::ios_base::app);

        // no gene tree parameters now
        if (_first_line) {
            _first_line = false;
    #if !defined PARALLELIZE_BY_GROUP
            logf << "iter ";
    #endif
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
        
        for (auto &p:v) {
    #if !defined PARALLELIZE_BY_GROUP
            unsigned iter = group_number * sample_size;
            logf << iter;
            iter++;
    #endif

            double log_coalescent_likelihood = 0.0;
    #if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
            log_coalescent_likelihood += p.getCoalescentLikelihood(1);
    #else
            for (unsigned g=1; g<ngenes+1; g++) {
                log_coalescent_likelihood += p.getCoalescentLikelihood(g);
            }
    #endif

            double vector_prior = 0.0;

            double log_prior = p.getAllPriors();

            double log_posterior = log_prior + log_coalescent_likelihood + vector_prior;

            logf << "\t" << log_posterior;
            
            logf << "\t" << log_prior;

            logf << "\t" << vector_prior;

            logf << "\t" << log_coalescent_likelihood;

            double species_tree_height = p.getSpeciesTreeHeight();
            logf << "\t" << species_tree_height;

            double species_tree_length = p.getSpeciesTreeLength();
            logf << "\t" << species_tree_length;

            double yule_model = p.getSpeciesTreePrior();
            logf << "\t" << yule_model;

            logf << "\t" << p.getPopMean() / 4.0; // beast uses Ne * u = theta / 4

            logf << "\t" << G::_lambda; // TODO: for now, not using estimate lambda option

            logf << endl;
        }

        logf.close();
    }


    inline void Proj::writeDeepCoalescenceFile(vector<Particle> &v) {
        ofstream logf("deep_coalescences.txt");
        logf << "num deep coalescences = " << v[0].getNumDeepCoalescences() << endl;
        logf << "Maximum number of deep coalescences = " << v[0].getMaxDeepCoalescences() << endl;
        logf << "True species tree height = " << v[0].getSpeciesTreeHeight() << endl;
        
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
    }

    inline void Proj::writeLoradFile(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v) const {
        ofstream logf("params.log");
        logf << "iteration ";
        logf << "\t" << "likelihood ";
        for (int s=0; s<nspecies-1; s++) {
            logf << "\t" << "species_increment";
        }
        logf << "\t" << "species_tree_prior";
        for (int g=1; g<ngenes+1; g++) {
            for (int i=1; i<ntaxa; i++) {
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

            for (unsigned g=0; g<ngenes+1; g++) {
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
            for (unsigned g=1; g<ngenes+1; g++) {
                log_coalescent_likelihood += p.getCoalescentLikelihood(g);
            }
            logf << "\t" << log_coalescent_likelihood;

            logf << endl;
        }

        logf.close();
    }

    inline void Proj::writeLoradFileAfterSpeciesFiltering(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle> &v) const {
        ofstream logf("params.log");
        logf << "iteration ";
        for (int s=0; s<nspecies-1; s++) {
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

#if defined (GRAHAM_JONES_COALESCENT_LIKELIHOOD)
            double log_coalescent_likelihood = p.getCoalescentLikelihood(0);
            logf << "\t" << log_coalescent_likelihood;
#else
            double log_coalescent_likelihood = 0.0;
            for (unsigned g=1; g<ngenes+1; g++) {
                log_coalescent_likelihood += p.getCoalescentLikelihood(g);
            }
            logf << "\t" << log_coalescent_likelihood;
#endif

            logf << endl;
        }

        logf.close();
    }

    inline void Proj::saveSpeciesTreesAltHierarchical(vector<Particle> &v) const {
            string filename1 = "alt_species_trees.trees";
            assert (G::_start_mode != "sim");

            unsigned count = 0;
            // save all species trees
            std::ofstream treef;

            treef.open(filename1, std::ios_base::app);
            for (auto &p:v) {
                treef << "  tree test = [&R] " << p.saveForestNewickAlt()  << ";\n";
                count++;
            }
            treef.close();
            }

    inline void Proj::saveSpeciesTreesHierarchical(vector<Particle> &v, string filename1, string filename2) const {
        // save only unique species trees
        if (!G::_run_on_empty) {
            vector<vector<pair<double, double>>> unique_increments_and_priors;

            std::ofstream unique_treef;

            unique_treef.open(filename2, std::ios_base::app);

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
            unique_treef.close();
        }

        assert (G::_start_mode != "sim");

        unsigned count = 0;
            // save all species trees
            std::ofstream treef;

            treef.open(filename1, std::ios_base::app);
            for (auto &p:v) {
                treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
                count++;
            }
            treef.close();
    }

    inline void Proj::saveSpeciesTreesAfterFirstRound(vector<Particle> &v) const {
        // save only unique species trees
        if (!G::_run_on_empty) {
            vector<vector<pair<double, double>>> unique_increments_and_priors;

            ofstream unique_treef("unique_species_trees_after_first_round.trees");
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

        if (G::_start_mode == "smc") {
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

    inline void Proj::saveGeneTrees(unsigned ngenes, vector<Particle> &v) const {
        if (G::_start_mode == "smc") {
            for (unsigned i=1; i<ngenes; i++) {
                string fname = "gene" + to_string(i) + ".trees";
                ofstream treef(fname);
                treef << "#nexus\n\n";
                treef << "begin trees;\n";
                for (auto &p:v) {
                    treef << "tree gene" << i << " = [&R] " << p.saveGeneNewick(i)  << ";\n";
                }
            }
        }

        else {
            if (G::_save_gene_trees_separately) {
                for (auto &p:v) {
                    for (unsigned i=1; i<ngenes+1; i++) {
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
                        for (unsigned i=1; i<ngenes+1; i++) {
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
        ("proposal",  boost::program_options::value(&G::_proposal)->default_value("prior-prior"), "a string defining a proposal (prior-prior or prior-post)")
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
        ("simasrvshape",   boost::program_options::value(&Forest::_asrv_shape)->default_value(Forest::_infinity), "Shape of gamma among-site rate heterogeneity within a locus (used only if startmode is 'sim')")
        ("simcomphet",   boost::program_options::value(&Forest::_comphet)->default_value(Forest::_infinity), "Dirichlet parameter governing compositional heterogeneity (default value results in compositional homogeneity (used only if startmode is 'sim')")
        ("save_gene_trees_separately", boost::program_options::value(&G::_save_gene_trees_separately)->default_value(false), "for simulations, save gene trees in separate files")
        ("newick_path", boost::program_options::value(&G::_newick_path)->default_value("."), "path to gene newicks are if starting from gene newicks and only performing SMC on second round")
        ("ngroups", boost::program_options::value(&G::_ngroups)->default_value(5), "number of populations")
        ("upgma", boost::program_options::value(&G::_upgma)->default_value(false), "set to true to use UPGMA completion")
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
            for (auto s : partition_subsets) {
                _partition->parseSubsetDefinition(s);
            }
        }

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

    inline unsigned Proj::setNumberTaxa(Data::SharedPtr) {
        unsigned ntaxa;
        ntaxa = _data->getNumTaxa();
        Forest::setNumTaxa(ntaxa);
        return ntaxa;
    }

    inline void Proj::handleSpeciesNewick(vector<Particle> particles) {
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
            while (getline(infile, newick)) {
                newick_string = newick;
                break;
            }
            int size_after = (int) newick.size();
            if (size_before == size_after) {
                throw XProj("cannot find gene newick file");
            }
            
            for (auto &p:particles) {
                p.processSpeciesNewick(newick_string);
                p.mapSpecies(_taxon_map, _species_names);
            }
    }

    inline void Proj::handleGeneNewicks() {
        vector<vector<string>> newicks; // vector of vector of newicks, 1 vector per gene
        _first_line = true;
        if (G::_ngenes_provided == 0) {
            throw XProj("must specify number of genes in the conf file");
        }
        
        for (int i=1; i<G::_ngenes_provided+1; i++) {
            vector<string> current_gene_newicks;

            string file_name = G::_newick_path + "/" + "gene" + to_string(i) + ".trees"; // file must be named gene1.trees, gene2.trees, etc.
            ifstream infile (file_name);
            string newick;
            unsigned size_before = (unsigned) current_gene_newicks.size();
            
            while (getline(infile, newick)) { // file newicks must start with the word "tree"
                if (current_gene_newicks.size() < G::_nparticles) { // stop adding newicks once the number of particles has been reached // TODO: add option to randomize this?
//                size_t found = newick.find("tree");
                if (newick.find("tree") == 2) { // TODO: this only works if there are exactly 2 spaces before the newick - try starting at parenthesis
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

        _data = Data::SharedPtr(new Data()); // TODO: don't need to set data b/c will not calculate a Felsenstein likelihood
        _data->setPartition(_partition);
        _data->getDataFromFile(G::_data_file_name);

        if (G::_verbose > 0) {
            summarizeData(_data);
        }
        createSpeciesMap(_data);

        // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
        if (G::_outgroup != "none") {
            checkOutgroupName();
        }

        //set number of species to number in data file
        unsigned ntaxa = setNumberTaxa(_data);
        unsigned nspecies = (unsigned) _species_names.size();
        Forest::setNumSpecies(nspecies);
        rng.setSeed(_random_seed);

//          create vector of particles

        unsigned nsubsets = _data->getNumSubsets();
        Particle::setNumSubsets(nsubsets);
        
        vector<Particle> my_vec;
        my_vec.resize(G::_nparticles);

        for (unsigned i=0; i<G::_nparticles; i++) {
            my_vec[i] = Particle();
        }

        initializeParticles(my_vec); // initialize in parallel with multithreading
        
        unsigned count = 0;
        for (auto &p:my_vec) {
            p.createSpeciesIndices();
            vector<string> particle_newicks;
            for (int i = 0; i<newicks.size(); i++) {
                // if number of particles > size of newicks, just use the same newick for each particle
                if (newicks[i].size() <= count) {
                    count = 0;
                }
                particle_newicks.push_back(newicks[i][count]);
            }
            assert (particle_newicks.size() == nsubsets);
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

        unsigned ngroups = round(G::_nparticles * G::_thin);
        if (ngroups == 0) {
            ngroups = 1;
            cout << "thin setting would result in 0 species groups; setting species groups to 1" << endl;
        }

        for (auto &p:my_vec) {
            // reset forest species partitions
            p.clearPartials(); // no more likelihood calculations
            p.resetSpecies();
            p.mapSpecies(_taxon_map, _species_names);
        }
        
        unsigned seed = rng.getSeed();
        std::shuffle(my_vec.begin(), my_vec.end(), std::default_random_engine(seed)); // shuffle particles using group seed
        
        // delete first (1-_thin) % of particles
        my_vec.erase(next(my_vec.begin(), 0), next(my_vec.begin(), (G::_nparticles-ngroups)));
        assert(my_vec.size() == ngroups);

        G::_nparticles = ngroups;

        vector<Particle> new_vec;

        G::_nparticles = G::_particle_increase; // TODO: G::_nparticles is static? be careful
        unsigned index = 0;
        
        bool parallelize_by_group = false;
        
        if (G::_nthreads > 1) {
#if defined PARALLELIZE_BY_GROUP
            parallelize_by_group = true;
#endif
        }
        
        // set group rng
        _group_rng.resize(ngroups);
        unsigned a=0;
        unsigned psuffix = 1;
        for (auto &g:_group_rng) {
            g.reset(new Lot());
            g->setSeed(rng.randint(1,9999)+psuffix);
            a++;
            psuffix += 2;
        }
        
        if (parallelize_by_group) {
        // don't bother with this if not multithreading
            
            unsigned group_number = 0;
            for (auto &p:my_vec) {
                p.setGroupNumber(group_number);
                group_number++;
            }
            
            proposeSpeciesGroups(my_vec, ngroups, filename1, filename2, filename3, nsubsets, ntaxa);
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
            
            unsigned group_number = 0;
            for (auto &p:my_vec) {
                p.setGroupNumber(group_number);
                group_number++;
            } // need to set group numbers because they are used in filtering
            
            for (unsigned a=0; a < ngroups; a++) {
//                        _log_species_tree_marginal_likelihood = 0.0; // TODO: for now, write lorad file for first set of species tree filtering and report the marginal likelihood for comparison
                                
                vector<Particle> use_vec;
                Particle p = my_vec[a];
                
                use_vec.resize(G::_particle_increase);
                
                fill(use_vec.begin(), use_vec.end(), p);
                
                assert(use_vec.size() == G::_particle_increase);

                index += G::_particle_increase;

                if (G::_verbose > 1) {
                    cout << "beginning species tree proposals for subset " << a+1 << endl;
                }
                for (unsigned s=0; s<nspecies-1; s++) {  // skip last round of filtering because weights are always 0
                    if (G::_verbose > 1) {
                        cout << "starting species step " << s+1 << " of " << nspecies-1 << endl;
                    }

                    // set particle random number seeds
                    unsigned group_number = use_vec[0].getGroupNumber();
                    unsigned psuffix = 1;
                    for (auto &p:use_vec) {
                        p.setSeed(_group_rng[group_number]->randint(1,9999) + psuffix);
                        psuffix += 2;
                    }

                    proposeSpeciesParticles(use_vec);
                    
                    
                    double ess = filterSpeciesParticles(s, use_vec);
                    
                    if (G::_verbose > 1) {
                        cout << "   " << "ESS = " << ess << endl;
                    }

                } // s loop

                if (G::_save_every > 1.0) { // thin sample for output by taking a random sample
                    unsigned sample_size = round (double (G::_particle_increase) / double(G::_save_every));
                    if (sample_size == 0) {
                        cout << "\n";
                        cout << "current settings would save 0 species trees; saving every species tree\n";
                        cout << "\n";
                        sample_size = G::_particle_increase;
                    }
                    
                    unsigned group_number = use_vec[0].getGroupNumber();
                    unsigned seed = _group_rng[group_number]->getSeed();
                    std::shuffle(use_vec.begin(), use_vec.end(), std::default_random_engine(seed)); // shuffle particles using group seed
                    
                    // delete first (1-_thin) % of particles
                    use_vec.erase(next(use_vec.begin(), 0), next(use_vec.begin(), (G::_particle_increase-sample_size)));
                    assert (use_vec.size() == sample_size);
                }

                saveSpeciesTreesHierarchical(use_vec, filename1, filename2);
                saveSpeciesTreesAltHierarchical(use_vec);
                writeParamsFileForBeastComparisonAfterSpeciesFilteringSpeciesOnly(nsubsets, nspecies, ntaxa, use_vec, filename3, a);
                if (a == 0) {
                    writeLoradFileAfterSpeciesFiltering(nsubsets, nspecies, ntaxa, use_vec); // testing the marginal likelihood by writing to file for lorad for first species group only
                    cout << "species tree log marginal likelihood is: " << _log_species_tree_marginal_likelihood << endl;
                }
            }

            std::ofstream treef;
            treef.open(filename1, std::ios_base::app);
            treef << "end;\n";
            treef.close();

            std::ofstream unique_treef;
            unique_treef.open(filename2, std::ios_base::app);
            unique_treef << "end;\n";
            unique_treef.close();
            
            // add iteration to params file
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


    inline void Proj::filterParticlesMixing(vector<Particle> & particles, unsigned step) {
//        string name = "mixing_count-" + to_string(step) + ".txt";
//        ofstream logf(name);
//        logf << "recipient_group" << "\t" << "donor_group" << endl;
        
        unsigned total_n_particles = G::_nparticles * G::_ngroups;
                
        double weight = (double) 1 / (G::_ngroups*G::_nparticles);
                
        // Copy log weights for all bundles to prob vector
        vector<double> probs(total_n_particles, weight);
        
        assert (particles.size() == total_n_particles);
                
        assert (probs.size() == total_n_particles);
        
          // Normalize log_weights to create discrete probability distribution
          double log_sum_weights = getRunningSum(probs);
                  
          transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});

          // Compute cumulative probabilities
          partial_sum(probs.begin(), probs.end(), probs.begin());

          // Initialize vector of counts storing number of darts hitting each particle
          vector<unsigned> counts (total_n_particles, 0);

          // Throw _nparticles darts
          for (unsigned i=0; i<total_n_particles; i++) {
              double u = rng.uniform();
              auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
              assert(it != probs.end());
              unsigned which = (unsigned)std::distance(probs.begin(), it);
              counts[which]++;
          } // TODO: with few groups, may wind up mostly resampling particles in the same pop
                  
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
              for (unsigned i = 0; i < total_n_particles; i++) {
                  if (counts[i] == 0)
                      nzeros++;
              }

              while (nzeros > 0) {
                  assert(donor < total_n_particles);
                  assert(recipient < total_n_particles);
                  
//                  unsigned recipient_group = recipient / G::_nparticles;
//                  unsigned recipient_index = recipient % G::_nparticles;
//
//                  unsigned donor_group = donor / G::_nparticles;
//                  unsigned donor_index = donor % G::_nparticles;
                  
//                  logf << recipient_group << "\t" << donor_group << endl;

                  // Copy donor to recipient
                  particles[recipient] = particles[donor];

                  counts[donor]--;
                  counts[recipient]++;
                  nzeros--;

                  if (counts[donor] == 1) {
                      // Move donor to next slot with count > 1
                      donor++;
                      while (donor < total_n_particles && counts[donor] < 2) {
                          donor++;
                      }
                  }

                  // Move recipient to next slot with count equal to 0
                  recipient++;
                  while (recipient < total_n_particles && counts[recipient] > 0) {
                      recipient++;
                  }
              }
        }
    }

    inline double Proj::filterParticles(unsigned step, vector<Particle> & particles, unsigned start, unsigned end) {
        // Copy log weights for all bundles to prob vector
      
        vector<double> probs(G::_nparticles, 0.0);
        
      unsigned count_index = 0;
      for (unsigned p=start; p<end+1; p++) {
          probs[count_index] = particles[p].getLogWeight();
          count_index++;
        }
      
      assert (probs.size() == G::_nparticles);
        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = getRunningSum(probs);
        
        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
        
        // Compute component of the log marginal likelihood due to this step
        _log_marginal_likelihood += log_sum_weights - log(G::_nparticles);
        if (step == 0) {
            _log_marginal_likelihood += _starting_log_likelihood;
        }
        
        double ess = 0.0;
        if (G::_verbose > 1) {
            // Compute effective sample size
            ess = computeEffectiveSampleSize(probs);
        }

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (G::_nparticles, 0);

        // Throw _nparticles darts
      for (unsigned i=start; i<end+1; i++) {
            double u = rng.uniform();
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
                particles[recipient+start] = particles[donor+start];

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

    }


//    inline double Proj::filterParticles(unsigned step, vector<Particle> & particles, unsigned start, unsigned end) {
//          // Copy log weights for all bundles to prob vector
//
//          vector<double> probs(G::_nparticles, 0.0);
//
//        unsigned count_index = 0;
//        for (unsigned p=start; p<end+1; p++) {
//            probs[count_index] = particles[p].getLogWeight();
//            count_index++;
//          }
//
//        assert (probs.size() == G::_nparticles);
//          // Normalize log_weights to create discrete probability distribution
//          double log_sum_weights = getRunningSum(probs);
//
//          transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});
//
//          // Compute component of the log marginal likelihood due to this step
//          _log_marginal_likelihood += log_sum_weights - log(G::_nparticles);
//          if (step == 0 && start == 0) {
//              _log_marginal_likelihood += _starting_log_likelihood;
//          }
//
//          double ess = 0.0;
//          if (G::_verbose > 1) {
//              // Compute effective sample size
//              ess = computeEffectiveSampleSize(probs);
//          }
//
//          // Compute cumulative probabilities
//          partial_sum(probs.begin(), probs.end(), probs.begin());
//
//          // Initialize vector of counts storing number of darts hitting each particle
//          vector<unsigned> counts (G::_nparticles, 0);
//
//          // Throw _nparticles darts
//        for (unsigned i=start; i<end+1; i++) {
//              double u = rng.uniform();
//              auto it = find_if(probs.begin(), probs.end(), [u](double cump){return cump > u;});
//              assert(it != probs.end());
//              unsigned which = (unsigned)std::distance(probs.begin(), it);
//              counts[which]++;
//          }
//
//          // Copy particles
//
//        bool copying_needed = true;
//
//          // Locate first donor
//          unsigned donor = 0;
//          while (counts[donor] < 2) {
//              donor++;
//              if (donor >= counts.size()) {
//                  copying_needed = false; // all the particle counts are 1
//                  break;
//              }
//          }
//
//        if (copying_needed) {
//              // Locate first recipient
//              unsigned recipient = 0;
//              while (counts[recipient] != 0) {
//                  recipient++;
//              }
//
//              // Count number of cells with zero count that can serve as copy recipients
//              unsigned nzeros = 0;
//              for (unsigned i = 0; i < G::_nparticles; i++) {
//                  if (counts[i] == 0)
//                      nzeros++;
//              }
//
//              while (nzeros > 0) {
//                  assert(donor < G::_nparticles);
//                  assert(recipient < G::_nparticles);
//
//                  // Copy donor to recipient
//                  particles[recipient+start] = particles[donor+start];
//
//                  counts[donor]--;
//                  counts[recipient]++;
//                  nzeros--;
//
//                  if (counts[donor] == 1) {
//                      // Move donor to next slot with count > 1
//                      donor++;
//                      while (donor < G::_nparticles && counts[donor] < 2) {
//                          donor++;
//                      }
//                  }
//
//                  // Move recipient to next slot with count equal to 0
//                  recipient++;
//                  while (recipient < G::_nparticles && counts[recipient] > 0) {
//                      recipient++;
//                  }
//              }
//        }
//          return ess;
//      }

    inline double Proj::filterSpeciesParticles(unsigned step, vector<Particle> & particles) {
        unsigned nparticles = (unsigned) particles.size();
        // Copy log weights for all bundles to prob vector
        vector<double> probs(nparticles, 0.0);
        
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
        
        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (nparticles, 0);

        // Throw _nparticles darts
        unsigned group_number = particles[0].getGroupNumber();
        
        for (unsigned i=0; i<nparticles; i++) {
//            double u = rng.uniform();
            double u = _group_rng[group_number]->uniform();
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
            _species_names.push_back(species_name);
        }
    }

    inline void Proj::createSpeciesMap(Data::SharedPtr d) {
        // this only works if names are in taxon^species format (no _)
        const vector<string> &names = d->getTaxonNames();
        for (auto &name:names) {
            regex re(".+\\^(.+)");
            smatch match_obj;
            bool matched=regex_match(name, match_obj, re); //search name for regular expression, store result in match_obj
            cout << "taxon name: " << name << endl;
            if (matched) {
                string species_name = match_obj[1];
                string taxon_name = name;
                if (find(_species_names.begin(), _species_names.end(), species_name) == _species_names.end()) {
                    _species_names.push_back(species_name);
                }
                _taxon_map[taxon_name]=species_name;
            }
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
        cout << "theta = " << Forest::_theta << endl;
#endif

//        cout << "speciation rate = " << Forest::_lambda << endl;
    }

    inline void Proj::proposeSpeciesGroupRange(unsigned first, unsigned last, vector<Particle> &particles, unsigned ngroups, string filename1, string filename2, string filename3, unsigned nsubsets, unsigned ntaxa) {
        unsigned nspecies = (unsigned) _species_names.size();
        
        for (unsigned i=first; i<last; i++){
            
            vector<Particle> use_vec;
            Particle p = particles[i];
            
            use_vec.resize(G::_particle_increase);
            
            fill(use_vec.begin(), use_vec.end(), p);

            assert(use_vec.size() == G::_particle_increase);

            assert(use_vec.size() == G::_particle_increase);

            if (G::_verbose > 1) {
                cout << "beginning species tree proposals for subset " << i+1 << endl;
            }
            for (unsigned s = 0; s < nspecies-1; s++) {  // skip last round of filtering because weights are always 0
                
                // set particle random number seeds
                unsigned psuffix = 1;
                unsigned group_number = use_vec[0].getGroupNumber();
                for (auto &p:use_vec) {
                    p.setSeed(_group_rng[group_number]->randint(1,9999) + psuffix);
                    psuffix += 2;
                }

                proposeSpeciesParticles(use_vec);

                double ess = filterSpeciesParticles(s, use_vec);
                if (G::_verbose > 1) {
                    cout << "   " << "ESS = " << ess << endl;
                }

            } // s loop
            
            if (G::_verbose > 1) {
                cout << "finished with species tree proposals for subset " << i+1 << endl;
            }
            
            if (G::_save_every > 1.0) { // thin sample for output by taking a random sample
                unsigned sample_size = round (double (G::_particle_increase) / double(G::_save_every));
                if (sample_size == 0) {
                    cout << "\n";
                    cout << "current settings would save 0 species trees; saving every species tree\n";
                    cout << "\n";
                    sample_size = G::_particle_increase;
                }

                unsigned group_number = use_vec[0].getGroupNumber();
               
                unsigned seed = _group_rng[group_number]->getSeed();
                std::shuffle(use_vec.begin(), use_vec.end(), std::default_random_engine(seed)); // shuffle particles using group seed

                // delete first (1-_thin) % of particles
                use_vec.erase(next(use_vec.begin(), 0), next(use_vec.begin(), (G::_particle_increase-sample_size)));
                assert (use_vec.size() == sample_size);
            }

            mtx.lock();
            saveSpeciesTreesHierarchical(use_vec, filename1, filename2);
            saveSpeciesTreesAltHierarchical(use_vec);
            _count++;
            if (G::_gene_newicks_specified) {
                writeParamsFileForBeastComparisonAfterSpeciesFilteringSpeciesOnly(nsubsets, nspecies, ntaxa, use_vec, filename3, i);
            }
            else {
                writeParamsFileForBeastComparisonAfterSpeciesFiltering(nsubsets, nspecies, ntaxa, use_vec, filename3, i);
            }
            mtx.unlock();
        }
    }

    inline void Proj::proposeSpeciesGroups(vector<Particle> &particles, unsigned ngroups, string filename1, string filename2, string filename3, unsigned nsubsets, unsigned ntaxa) {
        // ngroups = number of species SMCs to do (i.e. 100 particles for first round, thin = 1.0 means ngroups = 100 for this round)
        assert (G::_nthreads > 1);
        
        // divide up groups as evenly as possible across threads
        unsigned first = 0;
        unsigned incr = ngroups/G::_nthreads + (ngroups % G::_nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of groups
        unsigned last = incr;
        
        // need a vector of threads because we have to wait for each one to finish
        vector<thread> threads;

          while (true) {
          // create a thread to handle particles first through last - 1
            threads.push_back(thread(&Proj::proposeSpeciesGroupRange, this, first, last, std::ref(particles), ngroups, filename1, filename2, filename3, nsubsets, ntaxa));
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

    inline void Proj::proposeSpeciesParticles(vector<Particle> &particles) {
        assert(G::_nthreads > 0);
        if (G::_nthreads == 1) {
          for (auto & p : particles) {
              p.speciesOnlyProposal();
          }
        }
        else {
          // divide up the particles as evenly as possible across threads
          unsigned first = 0;
          unsigned incr = G::_nparticles/G::_nthreads + (G::_nparticles % G::_nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of particles
          unsigned last = incr;

          // need a vector of threads because we have to wait for each one to finish
          vector<thread> threads;

            while (true) {
            // create a thread to handle particles first through last - 1
              threads.push_back(thread(&Proj::proposeSpeciesParticleRange, this, first, last, std::ref(particles)));
            // update first and last
            first = last;
            last += incr;
            if (last > G::_nparticles) {
              last = G::_nparticles;
              }
            if (first>=G::_nparticles) {
                break;
            }
          }

          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }

    inline void Proj::initializeParticleRange(unsigned first, unsigned last, vector<Particle> &particles) {
        bool partials = false;

        for (unsigned i=first; i<last; i++){
            particles[i].setData(_data, _taxon_map, partials);
            partials = false;
            particles[i].mapSpecies(_taxon_map, _species_names);
            particles[i].setRelativeRatesByGene(G::_double_relative_rates);
        }
    }

    inline void Proj::initializeParticles(vector<Particle> &particles) {
        // set partials for first particle under save_memory setting for initial marginal likelihood calculation
        assert (G::_nthreads > 0);

        bool partials = true;
        if (G::_gene_newicks_specified) {
            partials = false;
            G::_save_memory = true;
        }

        if (G::_nthreads == 1) {
            for (auto & p:particles ) { // TODO: can initialize some of these things in parallel - probably not worth it
                p.setData(_data, _taxon_map, partials);
                partials = false;
                p.mapSpecies(_taxon_map, _species_names);
                
                p.setRelativeRatesByGene(G::_double_relative_rates);
                
            }
        }

        else {
            // always set partials for first particle under save memory setting for initial marginal likelihood calculation
            // for simplicity, do first particle separately under every setting
            bool partials = true;
            particles[0].setData(_data, _taxon_map, partials);
            particles[0].mapSpecies(_taxon_map, _species_names);
            particles[0].setRelativeRatesByGene(G::_double_relative_rates);

            if (particles.size() > 1) {
                // divide up the remaining particles as evenly as possible across threads
                unsigned first = 1;
                unsigned incr = (G::_nparticles*G::_ngroups-1) /G::_nthreads + ((G::_nparticles*G::_ngroups - 1) % G::_nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of particles

                unsigned last = incr;

                // need a vector of threads because we have to wait for each one to finish
                vector<thread> threads;

              while (true) {
                  // create a thread to handle particles first through last - 1
                    threads.push_back(thread(&Proj::initializeParticleRange, this, first, last, std::ref(particles)));
                  // update first and last
                  first = last;
                  last += incr;
                  
                  if (last > G::_nparticles*G::_ngroups) {
                    last = G::_nparticles*G::_ngroups;
                    }
                  if (first>=G::_nparticles*G::_ngroups) {
                      break;
                  }
            }

            // the join function causes this loop to pause until the ith thread finishes
            for (unsigned i = 0; i < threads.size(); i++) {
              threads[i].join();
            }
          }
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

    inline void Proj::proposeSpeciesParticleRange(unsigned first, unsigned last, vector<Particle> &particles) {
        for (unsigned i=first; i<last; i++){
            particles[i].speciesOnlyProposal();
        }
    }

    inline void Proj::debugSpeciesTree(vector<Particle> &particles) {
        cout << "debugging species tree" << endl;
        for (auto &p:particles) {
            p.showSpeciesJoined();
            p.showSpeciesIncrement();
            p.showSpeciesTree();
            cout << " _______ " << endl;
        }
    }

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
        
        if (G::_ntaxaperspecies.size() == 1) {
            ntaxa = G::_ntaxaperspecies[0];
            for (int i=0; i<G::_sim_nspecies-1; i++) {
                G::_ntaxaperspecies.push_back(ntaxa);
            }
        }
        else {
            for (auto &t:G::_ntaxaperspecies) {
                ntaxa += t;
            }
        }

        vector<Particle> sim_vec(1);
        sim_vec[0] = Particle();
        // set particle random number seed
        unsigned psuffix = 1;
        sim_vec[0].setSeed(rng.randint(1,9999) + psuffix);
        psuffix += 2;

        G::_run_on_empty = true;
        G::_proposal = "prior-prior";

        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);

        // make up the species map
        simSpeciesMap();

        vector<string> taxpartition;
        for (auto &t:_taxon_map) {
            taxpartition.push_back(t.second);
        }

        // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
        if (G::_outgroup != "none") {
            checkOutgroupName();
        }

        unsigned nsubsets = _data->getNumSubsets();
        Particle::setNumSubsets(nsubsets);
        
        assert (ntaxa > 0);
        
        unsigned list_size = (ntaxa-1)*nsubsets;
        vector<unsigned> gene_order;
        
        unsigned count = 1;
        vector<pair<double, unsigned>> randomize;
        for (unsigned l=0; l<list_size; l++) {
            if (count == 1) {
                randomize.clear();
            }
            randomize.push_back(make_pair(rng.uniform(), count));
            count++;
            if (count > nsubsets) {
                sort(randomize.begin(), randomize.end());
                for (auto &r:randomize) {
                    gene_order.push_back(r.second);
                }
                count = 1;
            }
        }
        
        assert (gene_order.size() == list_size);
        
        sim_vec[0].setGeneOrder(gene_order);

        sim_vec[0].setSimData(_data, _taxon_map, nsubsets, (unsigned) _taxon_map.size());

        sim_vec[0].mapSpecies(_taxon_map, _species_names);

        sim_vec[0].setNextSpeciesNumber(); // need to reset this now that number of species is known
        
        sim_vec[0].setNewTheta(G::_fix_theta_for_simulations);
                
        sim_vec[0].setNTaxaPerSpecies(G::_ntaxaperspecies);
        
        sim_vec[0].setRelativeRatesByGene(G::_double_relative_rates);

        unsigned nsteps = (unsigned) (_taxon_map.size()-1)*nsubsets;

        for (unsigned g=0; g<nsteps; g++){
            proposeParticles(sim_vec);
        }
        
        sim_vec[0].getNumDeepCoalescences();
        sim_vec[0].getMaxDeepCoalescences();

        cout << "\nBuilding species tree and associated gene trees....\n";
        vector<string> taxon_names;
        for (auto &t:_taxon_map) {
            taxon_names.push_back(t.first);
        }

        _data->setTaxonNames(taxon_names);

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

        saveGeneTrees(nsubsets, sim_vec);
        saveSpeciesTrees(sim_vec);

        writePaupFile(sim_vec, taxpartition);
        writeDeepCoalescenceFile(sim_vec);
        writeThetaFile(sim_vec);
    }

    inline void Proj::run() {
        if (G::_gene_newicks_specified) {
            if (G::_start_mode == "sim") {
                throw XProj("cannot specify gene newicks and simulations");
            }
            try {
                handleGeneNewicks();
            }
            catch (XProj & x) {
                std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
            }
        }
        
        else if (G::_start_mode == "sim") {
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
                cout << "Current working directory: " << boost::filesystem::current_path() << endl;
                cout << "Random seed: " << _random_seed << endl;
#if defined (DRAW_NEW_THETA)
                cout << "drawing new theta for each particle " << endl;
#else
                cout << "Theta: " << Forest::_theta << endl;
#endif
                cout << "Number of threads: " << G::_nthreads << endl;
            }

            if (G::_run_on_empty) { // if running with no data, choose taxa to join at random
                G::_proposal = "prior-prior";
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
                createSpeciesMap(_data);

                // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
                if (G::_outgroup != "none") {
                    checkOutgroupName();
                }

                //set number of species to number in data file
                unsigned ntaxa = setNumberTaxa(_data);
                unsigned nspecies = (unsigned) _species_names.size();
                Forest::setNumSpecies(nspecies);
                rng.setSeed(_random_seed);

    //          create vector of particles

                unsigned nsubsets = _data->getNumSubsets();
                Particle::setNumSubsets(nsubsets);
                
                vector<Particle> my_vec;
                my_vec.resize(G::_nparticles * G::_ngroups);
                
                for (unsigned p=0; p<my_vec.size(); p++) {
                    my_vec[p] = Particle();
                }

                initializeParticles(my_vec); // initialize in parallel with multithreading
      
                unsigned list_size = (ntaxa-1)*nsubsets;
                vector<unsigned> gene_order;
                
                unsigned count = 1;
                vector<pair<double, unsigned>> randomize;
                for (unsigned l=0; l<list_size; l++) {
                    if (count == 1) {
                        randomize.clear();
                    }
                    randomize.push_back(make_pair(rng.uniform(), count));
                    count++;
                    if (count > nsubsets) {
                        sort(randomize.begin(), randomize.end());
                        for (auto &r:randomize) {
                            gene_order.push_back(r.second);
                        }
                        count = 1;
                    }
                }
                
                assert (gene_order.size() == list_size);
                
                unsigned particle_num = 0;
                
                for (auto &p:my_vec) {
                    p.setGeneOrder(gene_order);
                    if (G::_upgma) {
                        if (particle_num == 0) {
                            p.calcStartingUPGMAMatrix();
                        }
                        else {
                            p.setStartingUPGMAMatrix(my_vec[0].getStartingUPGMAMatrix());
                        }
                        p.calcStartingRowCount();
                        particle_num++;
                    }
                }
                
//                if (_species_newick_name != "null") {
//                    handleSpeciesNewick(my_vec);
//                }
                
                // reset marginal likelihood
                _log_marginal_likelihood = 0.0;
                vector<double> starting_log_likelihoods = my_vec[0].calcGeneTreeLogLikelihoods(); // can't start at 0 because not every gene gets changed
                
                _starting_log_likelihood = 0.0;
                for (auto &l:starting_log_likelihoods) {
                    _starting_log_likelihood += l;
                }

                unsigned psuffix = 1;
                for (auto &p:my_vec) {
                        p.setSeed(rng.randint(1,9999) + psuffix);
                        psuffix += 2;
                    }

                    for (auto &p:my_vec) { // TODO: can parallelize this - is it worth it?
                        if (!G::_run_on_empty) {
                            p.setLogLikelihood(starting_log_likelihoods);
                        }
                        
                        if (G::_fix_theta) {
                            p.fixTheta();
                            p.setFixTheta(true);
                        }
                        
    #if defined (DRAW_NEW_THETA)
                        assert (!G::_fix_theta);
                        p.drawTheta();
    #endif
                    }
                
                if (G::_save_memory) {
                    my_vec[0].clearPartials(); // all other particles should have no partials
                }
                
                //run through each generation of particles

                    unsigned nsteps = (ntaxa-1)*nsubsets;
                
                    for (unsigned g=0; g<nsteps; g++){
                        if (G::_verbose > 0) {
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
                        
//                        StopWatch sw;
//                        sw.start();
                        proposeParticles(my_vec);
//                        double total_seconds = sw.stop();
//                        cout << "\nTotal time for proposal step " << g << " : " << total_seconds << endl;
//                        cout << total_seconds << endl;

                        bool filter = true;

                        if (G::_run_on_empty) {
                            filter = false;
                        }
                        
                        if (filter) {
                            
                            StopWatch sw;
                            
//                            sw.start();
                            for (unsigned i=0; i<G::_ngroups; i++) {
                                
                                unsigned start = i * G::_nparticles;
                                unsigned end = start + (G::_nparticles) - 1;
                                
                                double ess = filterParticles(g, my_vec, start, end); // filter only particles within a subgroup
                                
                                // TODO: can parallelize filtering now
                                
                                if (G::_verbose > 1) {
                                    cout << "\t" << "ESS is : " << ess << endl;
                                }
//                            double total_seconds = sw.stop();
//                            cout << "\nTotal time for filtering step " << g << " : " << total_seconds << endl;
//                            cout << total_seconds << endl;
                            }
                            
//                            saveAllSpeciesTrees(my_vec);
                            filterParticlesMixing(my_vec, g); // TODO: don't filter again at last step?
                            unsigned seed = rng.getSeed();
                            std::shuffle(my_vec.begin(), my_vec.end(), std::default_random_engine(seed)); // shuffle particles using group seed // TODO: trying shuffling for now to ensure particles are redistributed among populations
//                            saveAllSpeciesTrees(my_vec);
                        }
                } // g loop
                
                if (G::_save_gene_trees) {
                        for (int i=1; i<nsubsets+1; i++) {
                            saveGeneTrees(i, my_vec);
                        }
                }
                if (G::_verbose > 0) {
                    cout << "\n";
                    cout << "marginal likelihood after combined filtering: " << _log_marginal_likelihood << endl;
                    cout << "\n";
                }

#if !defined (HIERARCHICAL_FILTERING)
                writeParamsFileForBeastComparison(nsubsets, nspecies, ntaxa, my_vec);
                writePartialCountFile(my_vec);
                
                saveAllSpeciesTrees(my_vec);
#endif

#if defined (HIERARCHICAL_FILTERING)
                saveSpeciesTreesAfterFirstRound(my_vec);
                
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

                unsigned ngroups = round(_nparticles * _thin);
                if (ngroups == 0) {
                    ngroups = 1;
                    cout << "thin setting would result in 0 species groups; setting species groups to 1" << endl;
                }
                
                vector<unsigned> counts;
                for (unsigned index=0; index<my_vec.size(); index++) {
                  counts.push_back(index);
                }
                      
                unsigned seed = rng.getSeed();
                std::shuffle(counts.begin(), counts.end(), std::default_random_engine(seed));
                
                counts.erase(next(counts.begin(), 0), next(counts.begin(), (ngroups))); // choose what to delete - erase (thin) % of particle numbers
                sort (counts.begin(), counts.end(), greater<int>()); // sort highest to lowest for deletion of particles later

                // delete particles corresponding to those numbers
                for (unsigned c=0; c<counts.size(); c++) {
                    my_vec.erase(my_vec.begin() + counts[c]);
                }

                assert(my_vec.size() == ngroups);

                _nparticles = ngroups;

                for (auto &p:my_vec) {
                    // reset forest species partitions
                    p.clearPartials(); // no more likelihood calculations
                    p.resetSpecies();
                    p.mapSpecies(_taxon_map, _species_names);
                }

                vector<Particle> new_vec;

                _nparticles = _particle_increase;
                unsigned index = 0;
                
                bool parallelize_by_group = false;
            
                // set group rng
                _group_rng.resize(ngroups);
                unsigned a=0;
                for (auto &g:_group_rng) {
                    g.reset(new Lot());
                    g->setSeed(rng.randint(1,9999)+psuffix);
                    a++;
                    psuffix += 2;
                }
                
                if (_nthreads > 1) {
#if defined PARALLELIZE_BY_GROUP
                    parallelize_by_group = true;
#endif
                }
                
                if (parallelize_by_group) {
                // don't bother with this if not multithreading
                    unsigned group_number = 0;
                    for (auto &p:my_vec) {
                        p.setGroupNumber(group_number);
                        group_number++;
                    }
                    proposeSpeciesGroups(my_vec, ngroups, filename1, filename2, filename3, nsubsets, ntaxa);
                    
                    ofstream altfname;
                    altfname.open("alt_species_trees.trees", std::ios::app);
                    altfname << "end;" << endl;
                    altfname.close();
                    
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
                    unsigned group_number = 0;
                    for (auto &p:my_vec) {
                        p.setGroupNumber(group_number);
                        group_number++;
                    } // need to set group numbers because they are used in filtering
                    
                    for (unsigned a=0; a < ngroups; a++) {
//                        _log_species_tree_marginal_likelihood = 0.0; // for now, write lorad file for first set of species tree filtering and report the marginal likelihood for comparison
                        
                        vector<Particle> use_vec;
                        
                        use_vec.resize(_particle_increase);
                        
                        Particle chosen_particle = my_vec[a];
                                                
                        fill(use_vec.begin(), use_vec.end(), chosen_particle);
                        
                        assert(use_vec.size() == _particle_increase);

                        index += _particle_increase;

                        if (G::_verbose > 1) {
                            cout << "beginning species tree proposals for subset " << a+1 << endl;
                        }
                        for (unsigned s=0; s<nspecies-1; s++) {  // skip last round of filtering because weights are always 0
                            if (G::_verbose > 1) {
                                cout << "starting species step " << s+1 << " of " << nspecies-1 << endl;
                            }
                            
                            // set particle random number seeds
                            unsigned group_number = use_vec[0].getGroupNumber();
                            unsigned psuffix = 1;
                            for (auto &p:use_vec) {
                                p.setSeed(_group_rng[group_number]->randint(1,9999) + psuffix);
                                psuffix += 2;
                            }

                            proposeSpeciesParticles(use_vec);

                            double ess = filterSpeciesParticles(s, use_vec);
                            
                            if (G::_verbose > 1) {
                                cout << "   " << "ESS = " << ess << endl;
                            }

                        } // s loop

                        if (_save_every > 1.0) { // thin sample for output by taking a random sample
                            unsigned sample_size = round (double (_particle_increase) / double(_save_every));
                            if (sample_size == 0) {
                                cout << "\n";
                                cout << "current settings would save 0 species trees; saving every species tree\n";
                                cout << "\n";
                                sample_size = _particle_increase;
                            }

                            unsigned group_number = use_vec[0].getGroupNumber();
                            unsigned seed = _group_rng[group_number]->getSeed();
                            std::shuffle(use_vec.begin(), use_vec.end(), std::default_random_engine(seed)); // shuffle particles using group seed
                            
                            // delete first (1-_thin) % of particles
                            use_vec.erase(next(use_vec.begin(), 0), next(use_vec.begin(), (_particle_increase-sample_size)));
                             assert (use_vec.size() == sample_size);
                        }
                        
                        saveSpeciesTreesHierarchical(use_vec, filename1, filename2);
                        saveSpeciesTreesAltHierarchical(use_vec);
                        writeParamsFileForBeastComparisonAfterSpeciesFiltering(nsubsets, nspecies, ntaxa, use_vec, filename3, a);
                        if (a == 0) {
                            writeLoradFileAfterSpeciesFiltering(nsubsets, nspecies, ntaxa, use_vec); // testing the marginal likelihood by writing to file for lorad for first species group only
                            cout << "species tree log marginal likelihood is: " << _log_species_tree_marginal_likelihood << endl;
                        }
                    }

                    std::ofstream treef;
                    treef.open(filename1, std::ios_base::app);
                    treef << "end;\n";
                    treef.close();

                    std::ofstream unique_treef;
                    unique_treef.open(filename2, std::ios_base::app);
                    unique_treef << "end;\n";
                    unique_treef.close();
                    
                    std::ofstream alttrf;
                    alttrf.open("alt_species_trees.trees", std::ios_base::app);
                    alttrf << "end;\n";
                    alttrf.close();
                    
                    // add iterations to params file
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
#endif

            }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }
        }

        std::cout << "\nFinished!" << std::endl;
    }
}
