#pragma once

#include <iostream>
#include <fstream>
#include "data.hpp"
#include "partition.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "xproj.hpp"
#include "particle.hpp"
#include <vector>
#include <thread>
#include <boost/algorithm/string/split.hpp>
#include "conditionals.hpp"

using namespace std;
using namespace boost;
using namespace boost::algorithm;

#include "partial_store.hpp"
extern proj::PartialStore ps;

namespace proj {

    class Proj {
        public:
        
                                Proj();
                                ~Proj();

            void                clear();
            void                processCommandLineOptions(int argc, const char * argv[]);
            void                run();
            void                saveParticleWeights(vector<Particle::SharedPtr> &v) const;
            void                saveParticleLikelihoods(vector<Particle::SharedPtr> &v) const;

            void                normalizeWeights(vector<Particle::SharedPtr> & particles, string a, bool calc_marg_like);
            void                resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a);
            vector<int>         resampleSpeciesParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a);
            void                resetGeneParticles(vector<int> sel_indices, vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
            void                resetWeights(vector<Particle::SharedPtr> & particles, string a);
            void                createSpeciesMap(Data::SharedPtr);
            void                showFinal(vector<vector<Particle::SharedPtr>>);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, vector<pair<tuple<string, string, string>, double>> species_joined);
            void                proposeSpeciesParticleRange(unsigned first, unsigned last, vector<vector<Particle::SharedPtr>> &my_vec, unsigned s, unsigned nspecies, unsigned nsubsets);
            void                proposeSpeciesParticles(vector<vector<Particle::SharedPtr>> &my_vec, unsigned s, unsigned nspecies, unsigned nsubsets);
            void                proposeParticles(vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, Particle::SharedPtr species_tree_particle);
            void                growGeneTrees(vector<Particle::SharedPtr> &gene_particles, vector<Particle::SharedPtr> &my_vec_1_s, vector<Particle::SharedPtr> &my_vec_2_s, Particle::SharedPtr &species_tree_particle, unsigned ntaxa, unsigned ngenes, unsigned gene_number, unsigned iteration);
            void                growSpeciesTrees(vector<vector<Particle::SharedPtr>> &particles, vector<vector<Particle::SharedPtr>> &my_vec_1, vector<vector<Particle::SharedPtr>> &my_vec_2, unsigned ngenes, unsigned nspecies, unsigned nparticles);
            void                handleInitialNewicks(vector<vector<Particle::SharedPtr>> &particles, unsigned ngenes);
            void                saveAllHybridNodes(vector<Particle::SharedPtr> &v) const;
            void                writeGeneTreeFile();
            Particle::SharedPtr chooseTree(vector<Particle::SharedPtr> species_trees, string gene_or_species);
            void                writeLoradFile(vector<vector<Particle::SharedPtr>> my_vec, unsigned nparticles, unsigned nsubsets, unsigned nspecies, unsigned ntaxa);
            void                writeSpeciesTreeLoradFile(vector<Particle::SharedPtr> species_particles, unsigned nspecies);
            void                setStartingVariables();
            void                setUpInitialData();
            void                saveGeneAndSpeciesTrees(vector<Particle::SharedPtr> particles);
            void                updateTheta(Particle & species_particle, unsigned ntries, double delta, vector<Particle> & gene_particles);
            void                updateLambda(Particle & species_particle, unsigned ntries, double delta);
            void                estimateParameters(vector<vector<Particle::SharedPtr>> my_vec, Particle::SharedPtr species_tree_particle, unsigned ngenes);
            void                setUpForSpeciesFiltering(vector<vector<Particle::SharedPtr>> my_vec, unsigned ngenes, unsigned nparticles);
            double              calcLogSum(vector<double> vec);
            unsigned            multinomialDraw (const vector<double> & probs);
            void                proposeGeneParticlesFromPrior(vector<Particle::SharedPtr> &particles);
            void                proposeGeneParticlesFromPriorRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles);
            void                initializeTrees(vector<vector<Particle::SharedPtr>> particles, unsigned i, unsigned ngenes);
            void                sampleFromGeneTreePrior(vector<vector<Particle::SharedPtr>> &particles, unsigned ngenes, unsigned ntaxa, vector<vector<Particle::SharedPtr>> &my_vec_1, vector<vector<Particle::SharedPtr>> &my_vec_2);
            void                removeExtraParticles(vector<vector<Particle::SharedPtr>> &my_vec, vector<vector<Particle::SharedPtr>> &my_vec_1, vector<vector<Particle::SharedPtr>> &my_vec_2, unsigned nparticles, unsigned ngenes);
            void                resetGeneTreeWeights(vector<vector<Particle::SharedPtr>> &my_vec, unsigned ngenes);
            void                saveSelectedGeneTrees(vector<vector<Particle::SharedPtr>> &my_vec, unsigned ngenes);
        
        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            unsigned                    _nparticles;
            unsigned                    _random_seed;
            double                      _log_marginal_likelihood;
            double                      _species_tree_log_marginal_likelihood;
            bool                        _run_on_empty;
            bool                        _estimate_theta;
            bool                        _estimate_lambda;
            int                         _ntries_theta;
            int                         _ntries_lambda;

            static std::string          _program_name;
            static unsigned             _major_version;
            static unsigned             _minor_version;
            void                        summarizeData(Data::SharedPtr);
            unsigned                    setNumberTaxa(Data::SharedPtr);
            double                      getRunningSum(const vector<double> &) const;
            vector<string>              _species_names;
            map<string, string>         _taxon_map;
            vector<vector<Particle::SharedPtr>>            _accepted_particle_vec;
            unsigned                    _nthreads;
            double                      _hybridization_rate;
            void                        handleBaseFrequencies();
            void                        debugSpeciesTree(vector<Particle::SharedPtr> &particles);
            double                      _small_enough;
            int                         _niterations;
            string                      _species_newicks_name;
            string                      _gene_newicks_names;
            unsigned                    _species_particles_per_gene_particle;
            bool                        _sample_from_gene_tree_prior;
            bool                        _sample_from_species_tree_prior;
            bool                        _both;
            string                      _start;
            bool                        _use_first;
            bool                        _gene_first;
            bool                        _deconstruct;
    };

    inline Proj::Proj() {
//        std::cout << "Constructing a Proj" << std::endl;
        clear();
    }

    inline Proj::~Proj() {
//        std::cout << "Destroying a Proj" << std::endl;
    }

    inline void Proj::clear() {
        _data_file_name = "";
        _partition.reset(new Partition());
        _nparticles = 50000;
        _data = nullptr;
        _small_enough = 0.0000001;
    }

    inline void Proj::saveAllHybridNodes(vector<Particle::SharedPtr> &v) const {
        ofstream nodef("nodes.txt");
        for (auto &p:v) {
            nodef << "particle\n";
            nodef << p->saveHybridNodes()  << "\n";
        }
        nodef.close();
    }

    inline void Proj::saveGeneAndSpeciesTrees(vector<Particle::SharedPtr> particles) {
        // this function takes one species tree and associated gene trees and prints out the newicks + coalescent likelihood

        ofstream testf("coalescent-likelihood.txt");
        testf << "params : " << endl;
        testf << "\t" << "theta = " << Forest::_theta << endl;
        testf << "\t" << "lambda = " << Forest::_lambda << endl;
        
        testf << "coalescent likelihood = " << particles[0]->getCoalescentLikelihood() << endl;
        
        testf << "species tree: " << particles[0]->saveForestNewick() << endl;
        testf << endl;
        
        for (unsigned s=1; s < particles.size(); s++) {
            testf << "gene " << s << ": " << particles[s]->saveForestNewick() << endl;
            testf << endl;
        }
        
        testf.close();
        
    }

    inline void Proj::saveParticleWeights(vector<Particle::SharedPtr> &v) const {
        // this function saves particle weights + species tree newick
        
        ofstream weightf("species_trees.txt");
        weightf << "begin trees;" << endl;
        
        for (unsigned i=0; i < v.size(); i++) {
            weightf << "tree " << i << " = " << v[i]->getSpeciesNewick() << "; " << "\n";
        }
        
        weightf << "end trees;" << endl;
        weightf.close();
    }

    inline void Proj::saveParticleLikelihoods(vector<Particle::SharedPtr> &v) const {
        ofstream likelihoodf("likelihoods.txt");
        for (auto &p:v) {
            likelihoodf << "particle\n";
            likelihoodf << p->saveParticleLikelihoods() << "\n";
        }
        likelihoodf.close();
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        std::vector<std::string> partition_subsets;
        boost::program_options::variables_map vm;
        boost::program_options::options_description desc("Allowed options");
        
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile,d",  boost::program_options::value(&_data_file_name)->required(), "name of a data file in NEXUS format")
        ("subset",  boost::program_options::value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("nparticles",  boost::program_options::value(&_nparticles)->default_value(1000), "number of particles")
        ("seed,z", boost::program_options::value(&_random_seed)->default_value(1), "random seed")
        ("theta, t", boost::program_options::value(&Forest::_theta)->default_value(0.05), "theta")
        ("lambda", boost::program_options::value(&Forest::_lambda)->default_value(1), "lambda (speciation rate)")
        ("proposal",  boost::program_options::value(&Forest::_proposal)->default_value("prior-post"), "a string defining a proposal (prior-prior or prior-post or prior-post-ish)")
        ("model", boost::program_options::value(&Forest::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&Forest::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("nthreads",  boost::program_options::value(&_nthreads)->default_value(1.0), "number of threads for multi threading")
        ("migration_rate", boost::program_options::value(&Forest::_migration_rate)->default_value(0.0), "migration rate")
        ("hybridization_rate", boost::program_options::value(&Forest::_hybridization_rate)->default_value(0.0), "hybridization rate")
        ("run_on_empty", boost::program_options::value(&_run_on_empty)->default_value(false), "run with no data")
        ("niterations", boost::program_options::value(&_niterations)->default_value(1.0), "number of times to alternate between species and gene trees")
        ("gene_newicks", boost::program_options::value(&_gene_newicks_names)->default_value("null"), "name of file containing gene newick descriptions")
        ("species_newick", boost::program_options::value(&_species_newicks_name)->default_value("null"), "name of file containing species newick descriptions")
        ("species_particles_per_gene_particle", boost::program_options::value(&_species_particles_per_gene_particle)->default_value(1), "increase number of particles for species trees by this amount")
        ("outgroup", boost::program_options::value(&Forest::_outgroup)->default_value("null"), "specify outgroup in species tree")
        ("estimate_theta", boost::program_options::value(&_estimate_theta)->default_value(false), "estimate theta parameter")
        ("estimate_lambda", boost::program_options::value(&_estimate_lambda)->default_value(false), "estimate lambda parameter")
        ("ntries_theta", boost::program_options::value(&_ntries_theta)->default_value(100), "specify number of values of theta to try")
        ("ntries_lambda", boost::program_options::value(&_ntries_lambda)->default_value(100), "specify number of values of lambda to try")
        ("start_from_gene_tree_prior", boost::program_options::value(&_sample_from_gene_tree_prior)->default_value(false), "specify starting from gene tree prior")
        ("start_from_species_tree_prior", boost::program_options::value(&_sample_from_species_tree_prior)->default_value(false), "specify starting from species tree prior")
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
            std::cout << boost::str(boost::format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << std::endl;
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
    }

    inline void Proj::handleBaseFrequencies() {
        vector <string> temp;
        split(temp, Forest::_string_base_frequencies, is_any_of(","));
        double sum = 0.0;
        // iterate through temp
        for (auto &i:temp) {
            double f = stof(i);
            Forest::_base_frequencies.push_back(f);
            sum +=f;
        }
//        cout << sum-1 << endl;
        assert (fabs(sum-1) < 0.000001);
        if (fabs(sum-1)>0.000001) {
            throw XProj(format("base frequencies (%s) don't add to 1")%Forest::_string_base_frequencies);
        }
    }

    inline void Proj::summarizeData(Data::SharedPtr) {
        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        std::cout << "\nNumber of taxa: " << _data->getNumTaxa() << std::endl;
        std::cout << "Number of partition subsets: " << nsubsets << std::endl;
        std::cout << "Number of particles: " << _nparticles << std::endl;

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

    inline void Proj::normalizeWeights(vector<Particle::SharedPtr> & particles, string a, bool calc_marg_like) {
        unsigned nparticles = _nparticles;
        if (a == "s") {
            nparticles *= _species_particles_per_gene_particle;
        }
        unsigned i = 0;
        vector<double> log_weight_vec(nparticles);

        for (unsigned p=0; p<nparticles; p++ ) {
            log_weight_vec[i++] = particles[p]->getLogWeight(a);
        }

        double log_particle_sum = getRunningSum(log_weight_vec);

        for (unsigned p=0; p<nparticles; p++ ) {
            particles[p]->setLogWeight(particles[p]->getLogWeight(a) - log_particle_sum, a);
        }
        
        if (calc_marg_like) {
            _log_marginal_likelihood += log_particle_sum - log(_nparticles);
//            cout << setprecision(12) << "   " << _log_marginal_likelihood << endl;
        }
        else {
            _species_tree_log_marginal_likelihood += log_particle_sum - log(_nparticles);
//            cout << setprecision(12) << "   " << _species_tree_log_marginal_likelihood << endl;
        }
//        sort(particles.begin(), particles.end(), greater<Particle::SharedPtr>());
    }



    inline void Proj::writeGeneTreeFile() {
        // open log file
        ofstream genef("genetrees.txt");
        genef << "let gene_forests = [\n";
        for (unsigned p=0; p < _nparticles; p++) {
            genef << "particle " << endl;
            for (unsigned s=1; s < _accepted_particle_vec.size(); s++) {
                genef << "gene " << s << endl;
                vector<string> newicks = _accepted_particle_vec[s][p]->getGeneTreeNewicks();
                string newick = newicks[0];
                genef << newick << endl;
            }
        }
        genef << "];";
        genef.close();
    }

    inline void Proj::resetGeneParticles(vector<int> sel_indices, vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles) {
        
         unsigned nparticles = (unsigned)from_particles.size();
        
        for (unsigned i = 0; i < nparticles; i++) {
            int sel_index = sel_indices[i];
            assert(sel_index > -1);
            Particle::SharedPtr p0 = from_particles[sel_index];
            to_particles[i]=Particle::SharedPtr(new Particle(*p0));
            assert(nparticles == to_particles.size());
        }
    }

    inline double Proj::calcLogSum(vector<double> vec) {
        double running_sum = 0.0;
        double sum = 0.0;
        double log_max_weight = *max_element(vec.begin(), vec.end());
        for (auto & i:vec) {
            running_sum += exp(i - log_max_weight);
        }
        sum = log(running_sum) + log_max_weight;
        return sum;
    }

    inline void Proj::updateLambda(Particle & species_particle, unsigned ntries, double delta) {
        assert (_estimate_lambda);
        
        // Use multiple-try Metropolis to update lambda conditional on the
        // species forest defined in p. Uses the algorithm presented in
        // https://en.wikipedia.org/wiki/Multiple-try_Metropolis
        // assuming a symmetric proposal (so that w(x,y) = pi(x)).
        
        cout << "\nUpdating lambda...\n";

        // r is the rate of the lambda exponential prior
        double prior_rate = 1.0/Forest::_lambda_prior_mean;
        double log_prior_rate = log(prior_rate);
        double log_prior = log_prior_rate - prior_rate*Forest::_lambda;
        double log_likelihood = 0.0;

        // lambda0 is the current global lambda value
        double lambda0 = Forest::_lambda;
        
        // Sample ntries new values of lambda from symmetric proposal distribution
        // (window of width 2*delta centered on lambda0). Compute weights (coalescent
        // likelihood) for each proposed_lambdas value.
        vector<double> proposed_lambdas(ntries, 0.0);
        vector<double> logwstar(ntries, 0.0);
        for (unsigned i = 0; i < ntries; ++i) {
            double l = lambda0 - delta + 2.0*delta*rng.uniform();
            if (l < 0.0)
                l = -l;
            proposed_lambdas[i] = l;
            log_prior = log_prior_rate - prior_rate*l;
            log_likelihood =  species_particle.calcLogSpeciesTreeDensityGivenLambda(l);
            
            logwstar[i] = log_likelihood + log_prior;
        }
        
        // Compute log of the sum of the weights (this sum will form the
        // denominator of the acceptance ratio)
        double log_sum_denom_weights = calcLogSum(logwstar);

        // Normalize weights to create a discrete probability distribution
        vector<double> probs(ntries, 0.0);
        transform(logwstar.begin(), logwstar.end(), probs.begin(), [log_sum_denom_weights](double logw){
            return exp(logw - log_sum_denom_weights);
        });
        
        // Choose one lambda value from the probability distribution
        unsigned which = multinomialDraw(probs); // TODO: these must be on log scale
        double lambda_star = proposed_lambdas[which];
        
        // Sample ntries-1 new values of lambda from symmetric proposal distribution
        // (window of width 2*delta centered on lambda_star)
        log_prior = log_prior_rate - prior_rate*lambda0;
        log_likelihood =  species_particle.calcLogSpeciesTreeDensityGivenLambda(lambda0);
        logwstar[0] = log_likelihood + log_prior;
        for (unsigned i = 1; i < ntries; ++i) {
            double l = lambda_star - delta + 2.0*delta*rng.uniform();
            if (l < 0.0)
                l = -l;
            log_prior = log_prior_rate - prior_rate*l;
            log_likelihood =  species_particle.calcLogSpeciesTreeDensityGivenLambda(l);
            logwstar[i] = log_likelihood + log_prior;
        }
        
        // Compute log of the sum of the weights (this sum will form
        // the numerator of the acceptance ratio)
        double log_sum_numer_weights = calcLogSum(logwstar);

        // Compute acceptance ratio
        double logr = log_sum_numer_weights - log_sum_denom_weights;
        bool accept = true;
        if (logr < 0.0) {
            double logu = log(rng.uniform());
            accept = logu < logr;
        }
        if (accept) {
            Forest::_lambda = lambda_star;
            cout << str(format("  New lambda: %.5f\n") % Forest::_lambda);
        }
        else {
            cout << str(format("  Lambda unchanged: %.5f\n") % Forest::_lambda);
        }
    }

    inline void Proj::updateTheta(Particle & species_particle, unsigned ntries, double delta, vector<Particle> & gene_particles) {
        assert (_estimate_theta);
        // Use multiple-try Metropolis to update theta conditional on the gene forests
        // and species forest defined in p. Uses the algorithm presented in
        // https://en.wikipedia.org/wiki/Multiple-try_Metropolis
        // assuming a symmetric proposal (so that w(x,y) = pi(x)).
        
        cout << "\nUpdating theta...\n";

       // r is the rate of the theta exponential prior
       double prior_rate = 1.0/Forest::_theta_prior_mean;
       double log_prior_rate = log(prior_rate);
       double log_prior = log_prior_rate - prior_rate*Forest::_theta;

        // theta0 is the current global theta value
        double theta0 = Forest::_theta;
        assert (theta0 > 0.0);

        // Sample ntries new values of theta from symmetric proposal distribution
        // (window of width 2*delta centered on theta0). Compute weights (coalescent
        // likelihood) for each proposed_theta value.
        vector<double> proposed_thetas(ntries, 0.0);
        vector<double> logwstar(ntries, 0.0);
        
        vector<pair<tuple<string, string, string>, double>> species_info = species_particle.getSpeciesJoined();
        
        for (unsigned i = 0; i < ntries; ++i) {
            
            double q = theta0 - delta + 2.0*delta*rng.uniform();
            if (q < 0.0)
                q = -q;
            proposed_thetas[i] = q;
            double log_coalescent_likelihood = 0.0;
            for (unsigned s=0; s < gene_particles.size(); s++) {
                gene_particles[s].mapSpecies(_taxon_map, _species_names, s+1);
                gene_particles[s].refreshGeneTreePreorder();
                log_coalescent_likelihood += gene_particles[s].calcCoalLikeForNewTheta(q, species_info, false);
            }
            logwstar[i] = log_coalescent_likelihood + log_prior;
        }

        // Compute log of the sum of the weights (this sum will form the
        // numerator of the acceptance ratio)
        double log_sum_denom_weights = calcLogSum(logwstar);

        // Normalize weights to create a discrete probability distribution
        vector<double> probs(ntries, 0.0);
        transform(logwstar.begin(), logwstar.end(), probs.begin(), [log_sum_denom_weights](double logw){
            return exp(logw - log_sum_denom_weights);
        });

        // Choose one theta value from the probability distribution
        unsigned which = multinomialDraw(probs);
        double theta_star = proposed_thetas[which];

        // Sample ntries-1 new values of theta from symmetric proposal distribution
        // (window of width 2*delta centered on theta_star)
        log_prior = log_prior_rate - prior_rate * theta0;
        
        double log_coalescent_likelihood = 0.0;
        
        for (unsigned s=0; s < gene_particles.size(); s++) {
                gene_particles[s].mapSpecies(_taxon_map, _species_names, s+1);
                gene_particles[s].refreshGeneTreePreorder();
                log_coalescent_likelihood += gene_particles[s].calcCoalLikeForNewTheta(theta0, species_info, false);
            }
        
        logwstar[0] = log_coalescent_likelihood + log_prior;
        
        for (unsigned i=1; i<ntries; i++) {
            double q = theta_star - delta + 2.0 * delta * rng.uniform();
            if (q < 0.0) {
                q = -q;
            }
            
            log_prior = log_prior_rate - prior_rate * q;
            
            double log_coalescent_likelihood = 0.0;
            for (unsigned s=0; s < gene_particles.size(); s++) {
                gene_particles[s].mapSpecies(_taxon_map, _species_names, s+1);
                gene_particles[s].refreshGeneTreePreorder();
                log_coalescent_likelihood += gene_particles[s].calcCoalLikeForNewTheta(q, species_info, false);
                }

            logwstar[i] = log_coalescent_likelihood + log_prior;
        }

//         Compute log of the sum of the weights (this sum will form
        // the numerator of the acceptance ratio)
        double log_sum_numer_weights = calcLogSum(logwstar);

        // Compute acceptance ratio
        double logr = log_sum_numer_weights - log_sum_denom_weights;
        bool accept = true;
        if (logr < 0.0) {
            double logu = log(rng.uniform());
            accept = logu < logr;
        }
        if (accept) {
            Forest::_theta = theta_star;
            cout << str(format("\n*** new theta: %.5f\n") % theta_star);
        }
        else {
            cout << str(format("  Theta unchanged: %.5f\n") % Forest::_theta);
        }
    }

    inline unsigned Proj::multinomialDraw(const vector<double> & probs) {
        // Compute cumulative probababilities
        vector<double> cum_probs;
        
        cum_probs.resize(probs.size());
        partial_sum(probs.begin(), probs.end(), cum_probs.begin());
        
        // last element of cum_probs should hold 1
        if (abs(cum_probs.back() - 1.0) > _small_enough) {
            cout << "last element of cum_probs is " << cum_probs.back() << endl;
        }
        assert(abs((cum_probs.back() - 1.0)) < _small_enough);

        // Draw a Uniform(0,1) random deviate
        double u = rng.uniform();

        // Find first element in _cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6                <-- u
        //
        // cum_probs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(cum_probs.begin(), cum_probs.end(), [u](double cumpr){return cumpr > u;});
                
        assert(it != cum_probs.end());
        
        return (unsigned)std::distance(cum_probs.begin(), it);
    }

    inline Particle::SharedPtr Proj::chooseTree(vector<Particle::SharedPtr> particles, string gene_or_species) {
        normalizeWeights(particles, "gene_or_species", false);
        
        // get weights
        vector<double> log_weights;
        unsigned nparticles = _nparticles;
        if (gene_or_species == "s") {
            nparticles *= _species_particles_per_gene_particle;
        }
        for (unsigned p=0; p<nparticles; p++) {
            log_weights.push_back(particles[p]->getLogWeight(gene_or_species));
        }
        
//        // choose a random number [0,1]
        double u = rng.uniform();
        double cum_prob = 0.0;
        int index = 0.0;
        for (unsigned i=0; i < log_weights.size(); i++) {
            cum_prob += exp(log_weights[i]);
            if (u <= cum_prob) {
                index = i;
                break;
            }
        }
//        // return particle of choice
        return particles[index];
    }

    inline vector<int> Proj::resampleSpeciesParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a) {
        
        vector<int> sel_indices;
        unsigned nparticles = (unsigned)from_particles.size();
        vector<pair<double, double>> cum_probs;
            // Create vector of pairs p, with p.first = log weight and p.second = particle index
            cum_probs.resize(nparticles);
            unsigned i = 0;
            for(auto & p : from_particles) {
                cum_probs[i].first = p->getLogWeight(a);
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
            assert( fabs( 1.0 - cum_probs[_nparticles*_species_particles_per_gene_particle-1].first ) < Proj::_small_enough);

        
        // Draw new set of particles by sampling with replacement according to cum_probs
        to_particles.resize(_nparticles*_species_particles_per_gene_particle);
        for(unsigned i = 0; i < _nparticles*_species_particles_per_gene_particle; i++) {
        
            // Select a particle to copy to the ith slot in to_particles
            int sel_index = -1;
            double u = rng.uniform();
            for(unsigned j = 0; j < nparticles; j++) {
                if (u < cum_probs[j].first) {
                    sel_index = cum_probs[j].second;
                    break;
                }
            }
            assert(sel_index > -1);
            Particle::SharedPtr p0 = from_particles[sel_index];
            to_particles[i]=Particle::SharedPtr(new Particle(*p0));
            assert(nparticles == to_particles.size());
            
            sel_indices.push_back(sel_index);
        }
        return sel_indices;
    }

    inline void Proj::resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a) {
        
//        unsigned nparticles = (unsigned)from_particles.size();
        unsigned nparticles = _nparticles;
        
        vector<pair<double, double>> cum_probs;
            // Create vector of pairs p, with p.first = log weight and p.second = particle index
        cum_probs.resize(nparticles);
        unsigned i = 0;
        
        for (unsigned p=0; p < _nparticles; p++) {
                cum_probs[i].first = from_particles[p]->getLogWeight(a);
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
            assert( fabs( 1.0 - cum_probs[_nparticles-1].first ) < Proj::_small_enough);

        
        // Draw new set of particles by sampling with replacement according to cum_probs
//        to_particles.resize(_nparticles*_species_particles_per_gene_particle);
        to_particles.resize(_nparticles);
        for(unsigned i = 0; i < _nparticles; i++) {
        
            // Select a particle to copy to the ith slot in to_particles
            int sel_index = -1;
            double u = rng.uniform();
            for(unsigned j = 0; j < _nparticles; j++) {
                if (u < cum_probs[j].first) {
                    sel_index = cum_probs[j].second;
                    break;
                }
            }
            assert(sel_index > -1);
            Particle::SharedPtr p0 = from_particles[sel_index];
            to_particles[i]=Particle::SharedPtr(new Particle(*p0));
            
            assert(nparticles == to_particles.size());
        }
    }

    inline void Proj::resetWeights(vector<Particle::SharedPtr> & particles, string a) {
        unsigned nparticles = _nparticles;
        if (a == "s") {
            nparticles *= _species_particles_per_gene_particle;
        }
        double logw = -log(particles.size());

        for (unsigned p=0; p<nparticles; p++) {
            particles[p]->setLogWeight(logw, a);
        }
    }
    
    inline void Proj::createSpeciesMap(Data::SharedPtr d) {
        // TODO: this only works if names are in taxon^species format (no _)
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
    
    inline void Proj::showFinal(vector<vector<Particle::SharedPtr>> my_vec) {
        // this function displays the final species trees
//        for (int p=0; p<_nparticles; p++) {
//                my_vec[0][p]->showParticle();
//            }
        
        double sum_h = 0.0;
        for (auto & p:my_vec[0]) {
            double h = p->calcHeight();
            sum_h += h;
        }
        sum_h/=my_vec[0].size();
        cout << "species tree marg like: " << _species_tree_log_marginal_likelihood << endl;
        
        cout << "mean height equals " << sum_h << endl;
        cout << "log marginal likelihood = " << setprecision(12) << _log_marginal_likelihood << endl;
        cout << "theta = " << Forest::_theta << endl;
        cout << "speciation rate = " << Forest::_lambda << endl;
        cout << "hybridization rate = " << Forest::_hybridization_rate << endl;
    }

    inline void Proj::proposeParticles(vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, Particle::SharedPtr species_tree_particle) {
        assert(_nthreads > 0);
        vector<pair<tuple<string, string, string>, double>> species_joined = species_tree_particle->getSpeciesJoined();
        assert (species_joined.size() > 0);
        
        if (_nthreads == 1) {
            for (unsigned p=0; p<_nparticles; p++) {
                particles[p]->proposal(gene_trees_only, _deconstruct, species_joined);
          }
        }

        else {
          // divide up the particles as evenly as possible across threads
          unsigned first = 0;
          unsigned incr = _nparticles/_nthreads + (_nparticles % _nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of particles
          unsigned last = incr;

          // need a vector of threads because we have to wait for each one to finish
          vector<thread> threads;

            while (true) {
            // create a thread to handle particles first through last - 1
              threads.push_back(thread(&Proj::proposeParticleRange, this, first, last, std::ref(particles), gene_trees_only, a, species_joined));
            // update first and last
            first = last;
            last += incr;
            if (last > _nparticles) {
              last = _nparticles;
              }
            if (first>=_nparticles) {
                break;
            }
          }

          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, vector<pair<tuple<string, string, string>, double>> species_joined) {
        for (unsigned i=first; i<last; i++){
            particles[i]->proposal(gene_trees_only, _deconstruct, species_joined);
        }
    }

    inline void Proj::proposeSpeciesParticleRange(unsigned first, unsigned last, vector<vector<Particle::SharedPtr>> &my_vec, unsigned s, unsigned nspecies, unsigned nsubsets) {
            for (unsigned p=first; p<last; p++) {
                    tuple<string, string, string> species_joined = my_vec[0][p]->speciesTopologyProposal();

                    vector<double> max_depths; // this vector contains list of maximum depths for each gene tree

                    if (s < nspecies-1) {
                        for (unsigned j=1; j<nsubsets+1; j++) {
                            max_depths.push_back(my_vec[j][p]->calcConstrainedProposal(species_joined));
                        }

                        // now finish the species tree branch length proposal
                        my_vec[0][p]->speciesProposal(max_depths, species_joined);
                    }
                    
                    else {
                        for (unsigned j=1; j<nsubsets+1; j++) {
                            my_vec[j][p]->calcConstrainedProposal(species_joined); // update the gene tree species partitions
                            max_depths.push_back(0.0);
                        }
                        my_vec[0][p]->speciesProposal(max_depths, species_joined); // set last edge length of species tree to 0.0
                    }

                    double log_coalescent_likelihood = 0.0;

                    // calculate coalescent likelihood for each gene on each particle
                        for (unsigned j=1; j<nsubsets+1; j++) {
                            double last_edge_len = my_vec[0][p]->getLastEdgeLen();
                            double species_tree_height = my_vec[0][p]->getSpeciesTreeHeight();
                            log_coalescent_likelihood += my_vec[j][p]->calcGeneCoalescentLikelihood(last_edge_len, species_joined, species_tree_height);
                        }

                    my_vec[0][p]->calcSpeciesParticleWeight(log_coalescent_likelihood);
                } // p loop
    }

    inline void Proj::proposeSpeciesParticles( vector<vector<Particle::SharedPtr>> &my_vec, unsigned s, unsigned nspecies, unsigned nsubsets) {
        if (_nthreads == 1) {
            for (unsigned p=0; p<_nparticles*_species_particles_per_gene_particle; p++) {
                tuple<string, string, string> species_joined = my_vec[0][p]->speciesTopologyProposal();
                
                vector<double> max_depths; // this vector contains list of maximum depths for each gene tree
                
                if (s < nspecies-1) {
                    for (unsigned j=1; j<nsubsets+1; j++) {
                        
                        max_depths.push_back(my_vec[j][p]->calcConstrainedProposal(species_joined));
                    }
                    
                    // now finish the species tree branch length proposal
                    my_vec[0][p]->speciesProposal(max_depths, species_joined); // branch length proposal
                }
                else {
                    for (unsigned j=1; j<nsubsets+1; j++) {
                        my_vec[j][p]->calcConstrainedProposal(species_joined); // update the gene tree species partitions
                        max_depths.push_back(0.0);
                    }
                    my_vec[0][p]->speciesProposal(max_depths, species_joined); // set last edge length of species tree to 0.0
                }

                double log_coalescent_likelihood = 0.0;
                
                // calculate coalescent likelihood for each gene on each particle
                    for (unsigned j=1; j<nsubsets+1; j++) {
                        double last_edge_len = my_vec[0][p]->getLastEdgeLen();
                        double species_tree_height = my_vec[0][p]->getSpeciesTreeHeight();
                        log_coalescent_likelihood += my_vec[j][p]->calcGeneCoalescentLikelihood(last_edge_len, species_joined, species_tree_height);
                    }
                
                my_vec[0][p]->calcSpeciesParticleWeight(log_coalescent_likelihood);
            } // p loop
        }
            
        else { // multithreading
            // divide up the particles as evenly as possible across threads
            unsigned first = 0;
            unsigned nspecies_particles = _nparticles*_species_particles_per_gene_particle;
            unsigned incr = nspecies_particles/_nthreads + (nspecies_particles % _nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of particles
            unsigned last = incr;


            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;

            while (true) {
//                // create a thread to handle particles first through last - 1
                  threads.push_back(thread(&Proj::proposeSpeciesParticleRange, this, first, last, std::ref(my_vec), s, nspecies, nsubsets));
//                // update first and last
                first = last;
                last += incr;
                if (last > nspecies_particles) {
                  last = nspecies_particles;
                  }
                if (first>=nspecies_particles) {
                    break;
                }
          }

          // the join function causes this loop to pause until the ith thread finishes
              for (unsigned i = 0; i < threads.size(); i++) {
                threads[i].join();
              }
        }
    }

    inline void Proj::debugSpeciesTree(vector<Particle::SharedPtr> &particles) {
        cout << "debugging species tree" << endl;
        for (auto &p:particles) {
            p->showSpeciesIncrement();
            p->showSpeciesTree();
            cout << " _______ " << endl;
        }
    }

    inline void Proj::setUpInitialData() {
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_name);

        summarizeData(_data);
        createSpeciesMap(_data);

        //set number of species to number in data file
        setNumberTaxa(_data);
        
        
        Particle::_run_on_empty = false;
        if (_run_on_empty) {
            Particle::_run_on_empty = true;
        }
    }

    inline void Proj::proposeGeneParticlesFromPriorRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles) {
        
        for (unsigned i=first; i<last; i++){
            particles[i]->sampleGeneTreePrior();
        }
    }

    inline void Proj::proposeGeneParticlesFromPrior(vector<Particle::SharedPtr> &particles) {
        assert(_nthreads > 0);
        
        if (_nthreads == 1) {
            for (auto &p:particles) {
                p->sampleGeneTreePrior();
            }
        }

        else {
          // divide up the particles as evenly as possible across threads
          unsigned first = 0;
          unsigned incr = _nparticles/_nthreads + (_nparticles % _nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of particles
          unsigned last = incr;

          // need a vector of threads because we have to wait for each one to finish
          vector<thread> threads;

            while (true) {
            // create a thread to handle particles first through last - 1
                threads.push_back(thread(&Proj::proposeGeneParticlesFromPriorRange, this, first, last, std::ref(particles)));
                
            // update first and last
            first = last;
            last += incr;
            if (last > _nparticles) {
              last = _nparticles;
              }
            if (first>=_nparticles) {
                break;
            }
          }

          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
        
    }

    inline void Proj::handleInitialNewicks(vector<vector<Particle::SharedPtr>> &particles, unsigned ngenes) {
        vector<string> newicks;
        unsigned nparticles = _nparticles;
        
        if (_species_newicks_name != "null") {
            ifstream infile(_species_newicks_name);
            string newick;
            unsigned size_before = (int) newicks.size();
            while (getline(infile, newick)) {
                newicks.push_back(newick);
            }
            unsigned size_after = (int) newicks.size();
            if (size_before == size_after) {
                throw XProj("cannot find species newick file");
            }
        }
        
        if (_gene_newicks_names != "null") {
            if (_niterations == 1) {
                throw XProj(boost::str(boost::format("must specify more than 1 iteration if beginning from gene trees")));
            }
            ifstream infile(_gene_newicks_names);
            string newick;
            int size_before = (int) newicks.size();
            while (getline(infile, newick)) {
                newicks.push_back(newick);
            }
            int size_after = (int) newicks.size();
            if (size_before == size_after) {
                throw XProj("cannot find gene newick file");
            }
        }
        
        if (newicks.size() == 0 && !_sample_from_gene_tree_prior && !_sample_from_species_tree_prior) {
            throw XProj("Must specify gene newicks, species newicks, gene tree prior, or species tree prior");
        }
        
        _start = "species"; // start variable defines if program should start with gene or species trees
        
        for (unsigned s=0; s<ngenes+1; s++) {
            
            for (unsigned p=0; p<nparticles; p++) {
                particles[s][p]->setData(_data, _taxon_map, s);
                particles[s][p]->mapSpecies(_taxon_map, _species_names, s);
                particles[s][p]->setParticleGeneration(0);
                particles[s][p]->setLogLikelihood(0.0);
                particles[s][p]->setLogCoalescentLikelihood(0.0);
                particles[s][p]->setLogWeight(0.0, "g");
                particles[s][p]->setLogWeight(0.0, "s");
                
                // only sample 1 species tree, and use this tree for all the gene filtering
                if (s == 0 && _gene_newicks_names == "null" && p == 0) {
                    particles[0][p]->processSpeciesNewick(newicks, true); // if no newick specified, program will sample from species tree prior
                }
                
                if (s == 0 && _species_newicks_name != "null" && newicks.size() > 1) {
                    vector<string> species_newick;
                    species_newick.push_back(newicks[0]);
                    particles[0][p]->processSpeciesNewick(species_newick, false); // read in species newick
                    _both = true;
                }
            }
        }
        _gene_first = false;
        
        if (!_sample_from_gene_tree_prior) {
            
            if (_gene_newicks_names != "null") {
                if (_both) {
                    newicks.erase(newicks.begin()); // if specifying gene trees and species trees, erase the species newick because it's already been processed
                }
                assert (newicks.size() == ngenes);
                
                for (unsigned s=1; s<ngenes+1; s++) {
                    for (unsigned p=0; p<nparticles; p++) {
                        particles[s][p]->processGeneNewicks(newicks, s-1);
                        particles[s][p]->mapSpecies(_taxon_map, _species_names, s);
                        _start = "gene";
                        _gene_first = true;
                    }
                }
                
                _accepted_particle_vec = particles;
            }
        }
        
        // if both specified, calculate the coalescent likelihood and return it, ending the program
        if (_both) {
            cout << "...... calculating coalescent likelihood for specified trees ......" << endl;
            vector<pair<tuple<string, string, string>, double>> species_info = particles[0][0]->getSpeciesJoined();
            
            double log_coalescent_likelihood = 0.0;
            for (unsigned s=1; s<ngenes+1; s++) {
                    particles[s][0]->refreshGeneTreePreorder();
                    log_coalescent_likelihood += particles[s][0]->calcCoalLikeForNewTheta(Forest::_theta, species_info, _both);
                }
            cout << "log coalescent likelihood: " << log_coalescent_likelihood << endl;
            exit(0);
        }
    }

    inline void Proj::initializeTrees(vector<vector<Particle::SharedPtr>> particles, unsigned i, unsigned ngenes) {
        // reset trees for next iteration
        
        _species_tree_log_marginal_likelihood = 0.0;
        _log_marginal_likelihood = 0.0;
        
        unsigned nparticles = _nparticles;
        
        cout << "beginning iteration: " << i << endl;

        if (i > 0) {
            for (unsigned s=0; s<ngenes+1; s++) {
                for (unsigned p=0; p<nparticles; p++) {
                    particles[s][p]->setParticleGeneration(0);
                    if (s > 0) {
                        particles[s][p]->mapGeneTrees(_taxon_map, _species_names);
                        particles[s][p]->resetGeneIncrements();
                    }
                }
            }
        }
        
        // keep the species partition for the gene forests at this stage but clear the tree structure
         if (i == 1 && _gene_first == true) {
             for (unsigned s=1; s<ngenes+1; s++) {
                 // start at s=1 to only modify the gene trees
                 for (unsigned p=0; p<nparticles; p++) {
                     particles[s][p]->remakeGeneTrees(_taxon_map);
                     particles[s][p]->resetGeneTreePartials(_data, _taxon_map, s);
                 _deconstruct = false;
                 }
             }
         }
        
        if (i > 0) {
            _deconstruct = true;
        }
    }

    inline void Proj::sampleFromGeneTreePrior(vector<vector<Particle::SharedPtr>> &particles, unsigned ngenes, unsigned ntaxa, vector<vector<Particle::SharedPtr>> &my_vec_1, vector<vector<Particle::SharedPtr>> &my_vec_2) {
        
        _sample_from_gene_tree_prior = false;
        _start = "gene";
        for (unsigned s=1; s<ngenes+1; s++) {
            for (auto &p:particles[s]) {
                p->resetGeneTreePartials(_data, _taxon_map, s);
            }
        }
        for (unsigned g=0; g<ntaxa-1; g++) {
            cout << "prior generation " << g << endl;
            // filter particles within each gene
            
            for (unsigned s=1; s<ngenes+1; s++) { // skip species tree particles
                
                proposeGeneParticlesFromPrior(particles[s]);
                
                if (!_run_on_empty) {
                    bool calc_marg_like = true;
                    
                    normalizeWeights(particles[s], "g", calc_marg_like);
                    
                    double ess_inverse = 0.0;
                    
                    for (unsigned p=0; p<_nparticles; p++) {
                        ess_inverse += exp(2.0*particles[s][p]->getLogWeight("g"));
                    }

//                    double ess = 1.0/ess_inverse;
//                    cout << "   " << "ESS = " << ess << endl;
                 
                    resampleParticles(particles[s], _use_first ? my_vec_2[s]:my_vec_1[s], "g");
                    //if use_first is true, my_vec = my_vec_2
                    //if use_first is false, my_vec = my_vec_1
                    
                    particles[s] = _use_first ? my_vec_2[s]:my_vec_1[s];
                    // do not need to resample species trees; species tree will remain the same throughout all gene tree filtering
                    
                    assert(particles[s].size() == _nparticles);
                }
                //change use_first from true to false or false to true
                _use_first = !_use_first;
                if (g < ntaxa-2) {
                    resetWeights(particles[s], "g");
                }
                    assert (_accepted_particle_vec.size() == ngenes+1);
                    _accepted_particle_vec[s] = particles[s];
                } // s loop
        } // g loop
        writeGeneTreeFile();
    }

    inline void Proj::removeExtraParticles(vector<vector<Particle::SharedPtr>> &my_vec, vector<vector<Particle::SharedPtr>> &my_vec_1, vector<vector<Particle::SharedPtr>> &my_vec_2, unsigned nparticles, unsigned ngenes) {
        
        unsigned nparticles_to_remove = (nparticles*_species_particles_per_gene_particle) - nparticles;
        
        if (nparticles_to_remove > 0) {
            for (unsigned s=0; s<ngenes+1; s++) {
                my_vec[s].erase(my_vec[s].end() - nparticles_to_remove, my_vec[s].end());
                my_vec_2[s].erase(my_vec_2[s].end() - nparticles_to_remove, my_vec_2[s].end());
            }
        }
    }

    inline void Proj::resetGeneTreeWeights(vector<vector<Particle::SharedPtr>> &particles, unsigned ngenes) {
        // don't need to reset these variables for the first iteration
        for (unsigned s=1; s<ngenes+1; s++) {
            for (unsigned p=0; p<_nparticles; p++) {
                particles[s][p]->setLogLikelihood(0.0);
                particles[s][p]->setLogWeight(0.0, "g");
            }
        }
    }

    inline void Proj::saveSelectedGeneTrees(vector<vector<Particle::SharedPtr>> &particles, unsigned ngenes) {
        ofstream sel_gene_treesf("selected_gene_trees.txt");
        for (unsigned s=1; s<ngenes+1; s++) {
            sel_gene_treesf << "gene : " << s << " " << particles[s][0]->saveForestNewick();
        }
        sel_gene_treesf.close();
    }

    inline void Proj::growGeneTrees(vector<Particle::SharedPtr> &gene_particles, vector<Particle::SharedPtr> &my_vec_1_s, vector<Particle::SharedPtr> &my_vec_2_s, Particle::SharedPtr &species_tree_particle, unsigned ntaxa, unsigned ngenes, unsigned gene_number, unsigned iteration) {
        
        cout <<  "growing gene " << gene_number << endl;
        
        for (unsigned g=0; g<ntaxa-1; g++) {
            if (g == 0 && iteration > 0) {
                _deconstruct = true;
            }
//            cout << "generation " << g << endl;
            // filter particles within each gene
            
            bool gene_trees_only = true;

                proposeParticles(gene_particles, gene_trees_only, "g", species_tree_particle);

                if (!_run_on_empty) {
                    bool calc_marg_like = true;

                    normalizeWeights(gene_particles, "g", calc_marg_like);

                    double ess_inverse = 0.0;

                    for (unsigned p=0; p<_nparticles; p++) {
                        ess_inverse += exp(2.0*gene_particles[p]->getLogWeight("g"));
                    }

//                double ess = 1.0/ess_inverse;
//                cout << "   " << "ESS = " << ess << endl;

                    resampleParticles(gene_particles, _use_first ? my_vec_2_s:my_vec_1_s, "g");
                    //if use_first is true, my_vec = my_vec_2
                    //if use_first is false, my_vec = my_vec_1

                    gene_particles = _use_first ? my_vec_2_s:my_vec_1_s;
                    // do not need to resample species trees; species tree will remain the same throughout all gene tree filtering

                    assert(gene_particles.size() == _nparticles);
                }
                //change use_first from true to false or false to true
                _use_first = !_use_first;
                if (g < ntaxa-2) {
                    resetWeights(gene_particles, "g");
                }
                assert (_accepted_particle_vec.size() == ngenes+1);
                _accepted_particle_vec[gene_number] = gene_particles;
            _deconstruct = false;
        } // g loop
    }

    inline void Proj::estimateParameters(vector<vector<Particle::SharedPtr>> my_vec, Particle::SharedPtr species_tree_particle, unsigned ngenes) {
        if (_estimate_theta) {
            // try multiple theta values
            double delta_theta = 1.0;
            vector<Particle> gene_particles;
            for (unsigned s=1; s<ngenes+1; s++) {
                gene_particles.push_back(*my_vec[s][0]); // all gene trees are the same at this point, preserved from previous round of filtering
            }
            
            updateTheta(*species_tree_particle, _ntries_theta, delta_theta, gene_particles);
        }
        
        if (_estimate_lambda) {
            double delta_lambda = 50.0;
            updateLambda(*species_tree_particle, _ntries_lambda, delta_lambda);
        }
    }

    inline void Proj::setUpForSpeciesFiltering(vector<vector<Particle::SharedPtr>> particles, unsigned ngenes, unsigned nparticles) {
        _species_tree_log_marginal_likelihood = 0.0;
        
        for (unsigned p=0; p < particles[0].size(); p++) {
            particles[0][p]->mapSpecies(_taxon_map, _species_names, 0);
            particles[0][p]->resetSpecies();
        }
        
        for (unsigned s=1; s<ngenes+1; s++) {
            for (unsigned p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
                particles[s][p]->mapSpecies(_taxon_map, _species_names, s);
                particles[s][p]->refreshGeneTreePreorder();
                particles[s][p]->calcGeneTreeMinDepth(); // reset min depth vector for gene trees
                particles[s][p]->resetLogTopologyPrior();
            }
        }
        
        for (unsigned p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
            particles[0][p]->setLogLikelihood(0.0);
            particles[0][p]->setLogWeight(0.0, "g");
            particles[0][p]->setLogWeight(0.0, "s");
        }
    }

    inline void Proj::growSpeciesTrees(vector<vector<Particle::SharedPtr>> &particles, vector<vector<Particle::SharedPtr>> &particles_1, vector<vector<Particle::SharedPtr>> &particles_2, unsigned ngenes, unsigned nspecies, unsigned nparticles) {
        
        for (unsigned s=0; s<nspecies; s++) {
            cout << "beginning species tree proposals" << endl;
            
            if (s == 0) {
                for (unsigned j=1; j<ngenes+1; j++) {
                    // reset gene tree log coalescent likelihoods to 0
                    for (unsigned p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
                        particles[j][p]->setLogCoalescentLikelihood(0.0);
                    }
                }
            }
            
            proposeSpeciesParticles(particles, s, nspecies, ngenes);

            // filter - make sure all gene trees go along with correct species tree
            
            if (!_run_on_empty) {
                bool calc_marg_like = false;
                
                normalizeWeights(particles[0], "s", calc_marg_like);
                
                double ess_inverse = 0.0;
                
                for (auto & p:particles[0]) {
                    ess_inverse += exp(2.0*p->getLogWeight("s"));
                }

                double ess = 1.0/ess_inverse;
                cout << "   " << "ESS = " << ess << endl;
                vector<int> sel_indices = resampleSpeciesParticles(particles[0], _use_first ? particles_2[0]:particles_1[0], "s");
                //if use_first is true, my_vec = my_vec_2
                //if use_first is false, my_vec = my_vec_1
                
                particles[0] = _use_first ? particles_2[0]:particles_1[0];
                
                for (unsigned s=1; s<ngenes+1; s++) {
                    resetGeneParticles(sel_indices, particles[s], _use_first ? particles_2[s]:particles_1[s]);
                    particles[s] = _use_first ? particles_2[s]:particles_1[s];
                }
                
                //change use_first from true to false or false to true
                _use_first = !_use_first;
            }
            _accepted_particle_vec[0] = particles[0];
            _start = "species";

        } // s loop
    }

    inline void Proj::run() {
        cout << "Starting..." << endl;
        cout << "Current working directory: " << boost::filesystem::current_path() << endl;
        cout << "Random seed: " << _random_seed << endl;
        cout << "Theta: " << Forest::_theta << endl;
        cout << "Number of threads: " << _nthreads << endl;

        try {
            
            std::cout << "\n*** Reading and storing the data in the file " << _data_file_name << std::endl;
            std::cout << "data file name is " << _data_file_name << std::endl;
            
            setUpInitialData();
            
            unsigned nspecies = (unsigned) _species_names.size();
            Forest::setNumSpecies(nspecies);
            rng.setSeed(_random_seed);
            
            // create vector of particles
            unsigned nsubsets = _data->getNumSubsets();
            unsigned nparticles = _nparticles;
            
            Particle::setNumSubsets(nsubsets);
                            
            // set size of vectors (number of genes + 1 species tree)
            // my_vec contains a vector of vector of particles corresponding to species particles and gene particles
            vector<vector<Particle::SharedPtr>> my_vec_1(nsubsets+1, vector<Particle::SharedPtr> (nparticles)); // leave x% of vector empty for gene particles
            vector<vector<Particle::SharedPtr>> my_vec_2(nsubsets+1, vector<Particle::SharedPtr> (nparticles));
            vector<vector<Particle::SharedPtr>> &my_vec = my_vec_1;
            
            for (unsigned s=0; s<nsubsets+1; s++) {
                _accepted_particle_vec.push_back(vector<Particle::SharedPtr>(nparticles));
            }
            
            _log_marginal_likelihood = 0.0;
        
            for (unsigned s=0; s<nsubsets+1; s++) {
                unsigned nparticles = _nparticles;
                for (unsigned i=0; i<nparticles; i++) {
                    my_vec_1[s][i] = Particle::SharedPtr(new Particle); // TODO: wasteful that this makes species trees with extra lineages
                    my_vec_2[s][i] = Particle::SharedPtr(new Particle);
                }
            }

            _use_first = true;
            _both = false;
            
            handleInitialNewicks(my_vec, nsubsets);
            
            unsigned ntaxa = (unsigned) _taxon_map.size();
            _deconstruct = false;
            
            for (unsigned i=0; i<_niterations; i++) {
                initializeTrees(my_vec, i, nsubsets);
                
                if (_sample_from_gene_tree_prior) {
                    sampleFromGeneTreePrior(my_vec, nsubsets, ntaxa, my_vec_1, my_vec_2);
                }
                
                // my_vec[0] is the species tree particles
                // my_vec[1] is gene 1 particles
                // my_vec[2] is gene 2 particles
                // etc
                
                // filter gene trees
                if (_start == "species") {
                    // pick a species tree to use for all the gene trees for this step
                    
                    Particle::SharedPtr species_tree_particle;
                    
                    if (i > 0) {
                        species_tree_particle = chooseTree(my_vec[0], "s"); // pass in all the species trees
                        
                        if (_estimate_theta || _estimate_lambda) {
                            estimateParameters(my_vec, species_tree_particle, nsubsets);
                        }
            
                        // delete extra particles
                        removeExtraParticles(my_vec, my_vec_1, my_vec_2, nparticles, nsubsets);
                        
                    }
                    else {
                        species_tree_particle = my_vec[0][0];
                    }
                    
                    if (i > 0) {
                        // don't need to reset variables for first iteration
                        resetGeneTreeWeights(my_vec, nsubsets); // reset gene tree weights and likelihoods
                    }
                    
                    species_tree_particle->showParticle();
                    
                    for (unsigned s=1; s<nsubsets+1; s++) { // grow gene trees conditional on selected species tree
                        growGeneTrees(my_vec[s], my_vec_1[s], my_vec_2[s], species_tree_particle, ntaxa, nsubsets, s, i);
                    }
                }
                    
                // save gene tree variation for use in lorad file
                vector<vector<Particle>> variable_gene_trees(nsubsets, vector<Particle> (nparticles));

                for (unsigned s=1; s<nsubsets+1; s++) {
                    for (unsigned p=0; p<nparticles; p++) {
                        variable_gene_trees[s-1][p] = *my_vec[s][p];
                    }
                }
                    
                // choose one set of gene trees to use
                for (unsigned s=1; s<nsubsets+1; s++) {
                    Particle gene_x = *chooseTree(my_vec[s], "g");
                    for (unsigned p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
                        if (p<nparticles) {
                            *my_vec[s][p] = gene_x;
                        }
                        else {
                            // add in null particles here
                            my_vec[s].push_back(Particle::SharedPtr(new Particle));
                            *my_vec[s][p] = gene_x;
                            my_vec_2[s].push_back(Particle::SharedPtr(new Particle));
                        }
                    }
                }
                    
                // increase size of species vector
                for (unsigned p=nparticles; p<nparticles*_species_particles_per_gene_particle; p++) {
                    my_vec[0].push_back(Particle::SharedPtr(new Particle));
                    my_vec_2[0].push_back(Particle::SharedPtr(new Particle));
                }
                
                for (unsigned s=1; s<nsubsets+1; s++) {
                    for (unsigned p=0; p<nparticles; p++) {
                        *_accepted_particle_vec[s][p] = variable_gene_trees[s-1][p]; // preserve gene tree variation for output
                    }
                }
                    
                saveSelectedGeneTrees(my_vec, nsubsets);
                
                // build species trees
                
                if (i < _niterations-1) {
                    setUpForSpeciesFiltering(my_vec, nsubsets, nparticles);
                    growSpeciesTrees(my_vec, my_vec_1, my_vec_2, nsubsets, nspecies, nparticles); // grow and filter species trees conditional on selected gene trees
                    
                    if (i == _niterations - 2) {
                        writeSpeciesTreeLoradFile(my_vec[0], nspecies);
                    }
                    saveParticleWeights(my_vec[0]);
                }
            }
                                
            writeLoradFile(my_vec, nparticles, nsubsets, nspecies, ntaxa);
            
            writeGeneTreeFile();

//            saveAllHybridNodes(_accepted_particle_vec);
            showFinal(_accepted_particle_vec);
            cout << "marg like: " << setprecision(12) << _log_marginal_likelihood << endl;
        }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }

    inline void Proj::writeLoradFile(vector<vector<Particle::SharedPtr>> my_vec, unsigned nparticles, unsigned nsubsets, unsigned nspecies, unsigned ntaxa) {
        // open log file
        ofstream logf("params.log");
        
        double a = 0;
        unsigned col_count = 0;
        
        for (unsigned p=0; p<nparticles; p++) {
            
            vector<double> branch_length_vec;
            for (unsigned s=0; s<nsubsets+1; s++) {
                for (auto &b:my_vec[s][p]->getBranchLengths()) {
                    branch_length_vec.push_back(b);
                }
            }
            
            vector<double> prior_vec;
            for (unsigned s=0; s<nsubsets+1; s++) {
                for (auto &b:my_vec[s][p]->getBranchLengthPriors()) {
                    prior_vec.push_back(b);
                }
            }
            
            vector<double> gene_tree_log_like;
            for (unsigned s=1; s<nsubsets+1; s++) {
                for (auto &g:my_vec[s][p]->getGeneTreeLogLikelihoods()) {
                    gene_tree_log_like.push_back(g);
                }
            }
            
            vector<double> gene_tree_log_coalescent_like;
            for (unsigned s=1; s<nsubsets+1; s++) {
                for (auto &g:my_vec[s][p]->getGeneTreeLogCoalescentLikelihood()) {
                    gene_tree_log_coalescent_like.push_back(g);
                }
            }
            
            vector<double> log_topology_priors;
            for (unsigned s=0; s<nsubsets+1; s++) {
                for (auto &t:my_vec[s][p]->getTopologyPriors()) {
                    log_topology_priors.push_back(t);
                }
            }
            
            double species_tree_height = my_vec[0][p]->getSpeciesTreeHeight();
            
            assert(branch_length_vec.size() == prior_vec.size());
            
            double log_coalescent_likelihood = my_vec[0][p]->getCoalescentLikelihood();
            
            unsigned ngenes = nsubsets;
            
            if (col_count == 0) {
                logf << "iter" << "\t" << "theta";
                for (unsigned g=0; g<ngenes; g++) {
                    logf << "\t" << "gene_tree_log_like";
                }
                for (unsigned g=0; g<ngenes; g++) {
                    logf << "\t" << "gene_tree_log_coalescent_like";
                }
                for (unsigned i=0; i<nspecies-1; i++) {
                    logf << "\t" << "species_tree_increment" << "\t" << "increment_prior";
                }
                for (unsigned g=0; g<ngenes; g++) {
                    for (unsigned j=nspecies; j<nspecies+ntaxa-1; j++) {
                        logf << "\t" << "gene_tree_increment" << "\t" << "increment_prior";
                    }
                }
                logf << "\t" << "species_tree_topology_prior";
                for (unsigned g=0; g<ngenes; g++) {
                    logf << "\t" << "gene_tree_topology_prior";
                }
                logf << "\t" << "log_coal_like";
                
                logf << "\t" << "species_tree_height" << endl;
            }
            
            logf << a << "\t" << Forest::_theta;
            
            for (unsigned g=0; g<gene_tree_log_like.size(); g++) {
                logf << "\t" << setprecision(12) << gene_tree_log_like[g];
            }
            
            for (unsigned g=0; g<gene_tree_log_coalescent_like.size(); g++) {
                logf << "\t" << setprecision(12) << gene_tree_log_coalescent_like[g];
            }

            
            for (unsigned i=0; i<prior_vec.size(); i++) {
                logf << "\t" << setprecision(11) << branch_length_vec[i] << "\t" << prior_vec[i];
            }
            
            for (unsigned j=0; j < log_topology_priors.size(); j++) {
                logf << "\t" << setprecision(12) << log_topology_priors[j];
            }
            
            logf << "\t" << setprecision(12) << log_coalescent_likelihood;
            
            logf << "\t" << setprecision(12) << species_tree_height;
            
            logf << endl;
            a++;
            col_count++;
        }
        
        logf.close();
    }

    inline void Proj::writeSpeciesTreeLoradFile(vector<Particle::SharedPtr> species_particles, unsigned nspecies) {
        // open log file
        ofstream logf("species_params.log");
        
        double a = 0;
        unsigned col_count = 0;
        
        for (unsigned p=0; p<species_particles.size(); p++) {
            
            vector<double> branch_length_vec;
            for (auto &b:species_particles[p]->getBranchLengths()) {
                branch_length_vec.push_back(b);
            }
            
            vector<double> prior_vec;
            for (auto &b:species_particles[p]->getBranchLengthPriors()) {
                prior_vec.push_back(b);
            }
            
            vector<double> log_topology_priors;
            for (auto &t:species_particles[p]->getTopologyPriors()) {
                log_topology_priors.push_back(t);
            }
            
            double species_tree_height = species_particles[p]->getSpeciesTreeHeight();
            
            assert(branch_length_vec.size() == prior_vec.size());
            
            double log_coalescent_likelihood = species_particles[p]->getCoalescentLikelihood();
            
            if (col_count == 0) {
                logf << "iter" << "\t" << "theta";
                
                for (unsigned i=0; i<nspecies-1; i++) {
                    logf << "\t" << "species_tree_increment" << "\t" << "increment_prior";
                }
                
                logf << "\t" << "species_tree_topology_prior";
                
                logf << "\t" << "log_coal_like";
                
                logf << "\t" << "species_tree_height" << endl;
            }
            
            logf << a << "\t" << Forest::_theta;
            
            for (unsigned i=0; i<prior_vec.size(); i++) {
                logf << "\t" << setprecision(11) << branch_length_vec[i] << "\t" << prior_vec[i];
            }
            
            for (unsigned j=0; j < log_topology_priors.size(); j++) {
                logf << "\t" << setprecision(12) << log_topology_priors[j];
            }
            
            logf << "\t" << setprecision(12) << log_coalescent_likelihood;
            
            logf << "\t" << setprecision(12) << species_tree_height;
            
            logf << endl;
            a++;
            col_count++;
        }
        
        logf.close();
    }

}
