#pragma once

#include <iostream>
#include "data.hpp"
#include "partition.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "xproj.hpp"
#include "particle.hpp"
#include "bundle.hpp"
#include <vector>
#include <thread>
#include <boost/algorithm/string/split.hpp>
#include "conditionals.hpp"
#include <algorithm>
#include <random>
#include <cstdio>
#include <mutex>

using namespace std;
using namespace boost;
using namespace boost::algorithm;

#include "partial_store.hpp"
extern proj::PartialStore ps;
extern proj::Lot rng;

namespace proj {

    class Proj {
        public:

                                Proj();
                                ~Proj();

            void                clear();
            void                processCommandLineOptions(int argc, const char * argv[]);
            void                run();
            void                saveSpeciesTrees(vector<Bundle> &b) const;
            void                saveSpeciesTreesHierarchical(vector<Bundle> &b, string filename1) const;
            void                saveGeneTrees(unsigned ngenes, vector<Bundle> &v) const;
            void                saveGeneTree(unsigned gene_number, vector<Bundle> &b) const;
            void                resetWeights(vector<Bundle> & bundles);
            void                createSpeciesMap(Data::SharedPtr);
            void                showFinal(vector<Particle::SharedPtr> my_vec);
            void                initializeBundles(vector<Bundle> &particles);
            void                filterBundles(unsigned step, vector<Bundle> & bundles);
            void                runBundles(vector<Bundle> &b);
            void                runBundlesRange(unsigned first, unsigned last, vector<Bundle> &bundles);
            void                proposeSpeciesParticles(vector<Bundle> &b);
            void                proposeSpeciesParticleRange(unsigned first, unsigned last, vector<Bundle> &bundles);
            void                proposeSpeciesGroups(vector<Bundle> bundle_vec, unsigned ngroups, unsigned nsubsets, unsigned ntaxa, string filename1, string filename2);
            void                proposeSpeciesGroupRange(unsigned first, unsigned last, vector<Bundle> &bundles, unsigned ngroups, string filename1, unsigned nsubsets, unsigned ntaxa, string filename2);
            void                writeParamsFileForBeastComparison (unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Bundle> &b) const;
            void                writeParamsFileForBeastComparisonAfterSpeciesFiltering (unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Bundle> &b, string filename2, unsigned group_number);
        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood = 0.0;
            double                      _log_species_tree_marginal_likelihood = 0.0;
            bool                        _use_gpu;
            bool                        _ambig_missing;
            unsigned                    _nparticles;
            unsigned                    _random_seed;

            static string               _program_name;
            static unsigned             _major_version;
            static unsigned             _minor_version;
            static string               _start_mode;
            void                        summarizeData(Data::SharedPtr);
            unsigned                    setNumberTaxa(Data::SharedPtr);
            double                      getRunningSum(const vector<double> &) const;
            vector<string>              _species_names;
            map<string, string>         _taxon_map;
            unsigned                    _nthreads;
            void                        handleBaseFrequencies();
            void                        checkOutgroupName();
            void                        debugSpeciesTree(vector<Particle::SharedPtr> &particles);
            double                      _small_enough;
            unsigned                    _verbose;
            unsigned                    _particle_increase;
            double                      _thin;
            unsigned                    _save_every;
            bool                        _save_gene_trees;
            bool                        _first_line;
            unsigned                    _count; // counter for params output file
            bool                        _gene_newicks_specified;
            unsigned                    _ngenes_provided;
            string                      _species_newick_name;
            unsigned                    _nbundles;
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
        _use_gpu        = true;
        _ambig_missing  = true;
        _nparticles = 50000;
        _data = nullptr;
        _small_enough = 0.0000001;
        _nbundles = 1.0;
    }

    inline void Proj::saveSpeciesTreesHierarchical(vector<Bundle> &b, string filename1) const {
          assert (_start_mode != "sim");

          unsigned count = 0;
          // save all species trees
          std::ofstream treef;

          treef.open(filename1, std::ios_base::app);
          for (auto &p:b) {
              treef << "  tree test = [&R] " << p.saveSpeciesNewick()  << ";\n";
              count++;
          }
          treef.close();
    }

    inline void Proj::saveSpeciesTrees(vector<Bundle> &b) const {
        // save only unique species trees
        if (!Forest::_run_on_empty) {
            vector<vector<pair<double, double>>> unique_increments_and_priors;

            ofstream unique_treef("unique_species_trees_after_first_round.trees");
            unique_treef << "#nexus\n\n";
            unique_treef << "begin trees;\n";
            for (auto &p:b) {
                vector<pair<double, double>> increments_and_priors = p.getSpeciesTreeIncrementPriors();
                
                bool found = false;
                if(std::find(unique_increments_and_priors.begin(), unique_increments_and_priors.end(), increments_and_priors) != unique_increments_and_priors.end()) {
                    found = true;
                }
                if (!found) {
                    unique_increments_and_priors.push_back(increments_and_priors);
                    unique_treef << "  tree test = [&R] " << p.saveSpeciesNewick()  << ";\n";
                }
            }
            unique_treef << "end;\n";
            unique_treef.close();
        }

        if (_start_mode == "smc") {
            // save all species trees
            ofstream treef("species_trees_after_first_round.trees");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:b) {
                treef << "  tree test = [&R] " << p.saveSpeciesNewick()  << ";\n";
            }
            treef << "end;\n";
            treef.close();
        }
    }

    inline void Proj::writeParamsFileForBeastComparison(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Bundle> &bundles) const {
        
        // this function creates a params file that is comparable to output from starbeast3
        string filename1 = "params-beast-comparison-after-first-round.log";
        
        if (filesystem::remove(filename1)) {
            if (_verbose > 0) {
                cout << "existing file " << filename1 << " removed and replaced\n";
            }
        }
        
        ofstream logf(filename1);
        logf << "iter";
        for (unsigned i=1; i<ngenes+1; i++) {
            logf << "\t" << "Tree.t:gene" + to_string(i) + "likelihood";
        }

        logf << "\t" << "popMean";
        
        for (unsigned s=0; s<nspecies*2 - 1; s++) {
            logf << "\t" << "popMean." << s+1;
        }
        
        logf << endl;
        
        unsigned iter = 0;
        
        logf << iter << "\t";
        iter++;
        
        for (auto &b:bundles) {
            vector<double> log_likelihoods = b.getLogLikelihoods();
        
            unsigned gene_count = 0;
            
            for (unsigned p=0; p<log_likelihoods.size(); p++) {
                logf << log_likelihoods[p] << "\t";
                
                gene_count++;
                if (gene_count == ngenes) {
                    logf << b.getThetaMean() / 4.0 << "\t";
                    
#if defined (DRAW_NEW_THETA)
                    for (auto &t:b.getThetas()) {
                        logf << t / 4.0 << "\t"; // theta / 4 to be consistent with *beast3 output
                    }
#else
                    for (unsigned t=0; t<nspecies*2 - 1; t++) {
                        logf << Forest::_theta / 4 << "\t"; // theta / 4 to be consistent with *beast3 output
                    }
#endif
                    
                    logf << endl;
                    gene_count = 0;
                    
                    if (iter < _nparticles * _nbundles) { // don't add this for last one
                        logf << iter << "\t";
                        iter++;
                    }
                }
            }
        }
        
        logf.close();
    }

    inline void Proj::writeParamsFileForBeastComparisonAfterSpeciesFiltering(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Bundle> &bundles, string filename, unsigned group_number) {
        // save parameters of every gene and species particle (i.e. 1 bundle = n gene particles and 1 species particle - maybe not the best way to save because each species particle is represented particle_increase times
        
        std::ofstream logf;

        logf.open(filename, std::ios_base::app);
        
        if (_first_line) {
            _first_line = false;

#if !defined (PARALLELIZE_BY_GROUP)
            logf << "iter";
#endif
            if (_nthreads == 1) {
                logf << "iter";
            }
            
            for (unsigned i=1; i<ngenes+1; i++) {
                logf << "\t" << "Tree.t:gene" + to_string(i) + "likelihood";
            }

            logf << "\t" << "popMean";
            
            for (unsigned s=0; s<nspecies*2 - 1; s++) {
                logf << "\t" << "popMean." << s+1;
            }
            
            logf << endl;
            
        }
                
        unsigned iter = 0;
#if !defined (PARALLELIZE_BY_GROUP)
        unsigned sample_size = round(double (_particle_increase) / double(_save_every) );
        if (sample_size == 0) {
            sample_size = _particle_increase;
        }
        
        iter = group_number * sample_size;
#endif
        
        if (_nthreads == 1) {
            unsigned sample_size = round(double (_particle_increase) / double(_save_every) );
            if (sample_size == 0) {
                sample_size = _particle_increase;
            }
            
           iter = group_number * sample_size;
        }
        
        for (auto &b:bundles) {
#if !defined (PARALLELIZE_BY_GROUP)
            logf << iter << "\t";
            iter++;
#endif
            
            if (_nthreads == 1) {
                logf << iter << "\t";
                iter++;
            }
            
            vector<double> log_likelihoods = b.getLogLikelihoods();
        
            unsigned gene_count = 0;
            
            for (unsigned p=0; p<log_likelihoods.size(); p++) {
                logf << log_likelihoods[p] << "\t";
                
                gene_count++;
                if (gene_count == ngenes) {
                    logf << b.getThetaMean() / 4.0 << "\t";
                    
    #if defined (DRAW_NEW_THETA)
                    for (auto &t:b.getThetas()) {
                        logf << t / 4.0 << "\t"; // theta / 4 to be consistent with *beast3 output
                    }
    #else
                    for (unsigned t=0; t<nspecies*2 - 1; t++) {
                        logf << Forest::_theta / 4 << "\t"; // theta / 4 to be consistent with *beast3 output
                    }
    #endif
                    
                    logf << endl;
                    gene_count = 0;
                }
            }
        }
        
        logf.close();
    }


    inline void Proj::saveGeneTrees(unsigned ngenes, vector<Bundle> &b) const {
        if (_start_mode == "smc") {
            ofstream treef("gene_trees.trees");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:b) {
                    for (unsigned i=0; i<ngenes; i++) {
                        treef << "tree gene" << i << " = [&R] " << p.saveGeneNewick(i)  << ";\n";
                }
                treef << endl;
            }
            treef << "end;\n";
            treef.close();
        }

        else {
            // save true species tree
            ofstream treef("true-gene-trees.tre");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:b) {
                    for (unsigned i=0; i<ngenes; i++) {
                        treef << "tree gene" << i << " = [&R] " << p.saveGeneNewick(i)  << ";\n";
                }
                treef << endl;
            }
            treef << "end;\n";
            treef.close();
        }
    }

    inline void Proj::saveGeneTree(unsigned gene_number, vector<Bundle> &b) const {
        string name = "gene" + to_string(gene_number+1) + ".trees";
        ofstream treef(name);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:b) {
            treef << "  tree test = [&R] " << p.saveGeneNewick(gene_number)  << ";\n";
            treef << endl;
        }
        treef << "end;\n";

        treef.close();
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {
        std::vector<std::string> partition_subsets;
        boost::program_options::variables_map vm;
        boost::program_options::options_description desc("Allowed options");

        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile,d",  boost::program_options::value(&_data_file_name), "name of a data file in NEXUS format")
        ("subset",  boost::program_options::value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("gpu",           boost::program_options::value(&_use_gpu)->default_value(true), "use GPU if available")
        ("ambigmissing",  boost::program_options::value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("nparticles",  boost::program_options::value(&_nparticles)->default_value(1000), "number of particles")
        ("seed,z", boost::program_options::value(&_random_seed)->default_value(1), "random seed")
        ("theta, t", boost::program_options::value(&Forest::_theta)->default_value(0.0), "theta")
        ("lambda", boost::program_options::value(&Forest::_lambda)->default_value(1), "speciation rate")
        ("model", boost::program_options::value(&Forest::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&Forest::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("nthreads",  boost::program_options::value(&_nthreads)->default_value(1), "number of threads for multi threading")
        ("run_on_empty", boost::program_options::value(&Forest::_run_on_empty)->default_value(false), "run program without data")
        ("verbose", boost::program_options::value(&_verbose)->default_value(1), "set amount of output printed")
        ("save_memory", boost::program_options::value(&Forest::_save_memory)->default_value(false), "save memory at the expense of time")
        ("outgroup", boost::program_options::value(&Forest::_outgroup)->default_value("none"), "a string defining the outgroup")
        ("startmode", boost::program_options::value(&_start_mode)->default_value("smc"), "a string defining whether to simulate data or perform smc")
        ("particle_increase", boost::program_options::value(&_particle_increase)->default_value(1), "how much to increase particles for species filtering")
        ("thin", boost::program_options::value(&_thin)->default_value(1.0), "take this portion of particles for hierarchical species filtering")
        ("save_every", boost::program_options::value(&_save_every)->default_value(1.0), "take this portion of particles for output")
        ("save_gene_trees", boost::program_options::value(&_save_gene_trees)->default_value(true), "turn this off to not save gene trees and speed up program")
        ("gene_newicks", boost::program_options::value(&_gene_newicks_specified)->default_value(false), "set true if user is specifying gene tree files")
        ("ngenes", boost::program_options::value(&_ngenes_provided)->default_value(0), "number of gene newick files specified")
        ("theta_proposal_mean", boost::program_options::value(&Forest::_theta_proposal_mean)->default_value(0.0), "theta proposal mean")
        ("theta_prior_mean", boost::program_options::value(&Forest::_theta_prior_mean)->default_value(0.0), "theta prior mean")
        ("species_newick", boost::program_options::value(&_species_newick_name)->default_value("null"), "name of file containing species newick descriptions")
        ("nbundles", boost::program_options::value(&_nbundles)->default_value(1), "number of bundles")
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
        
        // if save_every > particle_increase, quit
        if (_save_every > _particle_increase) {
            throw XProj("particle_increase must be greater than or equal to save_every");
        }
        
        if (Forest::_model == "JC") {
            cout << "Setting kappa to 1.0 under JC model\n";
            cout << "Setting base frequencies equal under JC model\n";
            if (Forest::_kappa != 1.0) {
                cout << "\nIgnoring kappa under JC model\n";
            }
        }
        if (_start_mode == "sim") {
            if (_data_file_name != "") {
                cout << "\nIgnoring data file name for simulation\n";
            }

            if (Forest::_run_on_empty) {
                cout << "\nIgnoring start_mode = run_on_empty and simulating data\n";
            }

            if (Forest::_theta == 0.0 && Forest::_theta_prior_mean == 0.0 && Forest::_theta_proposal_mean == 0.0) {
                throw XProj("must specify theta or theta proposal / prior mean for simulations");
            }
        }
        else {
            if (_data_file_name == "") {
                throw XProj("must specify name of data file if smc option is chosen; ex. data file = sim.nex");
            }
            if (Forest::_theta_prior_mean == 0.0 && Forest::_theta_proposal_mean > 0.0) {
                cout << boost::format("\nSetting theta prior mean equal to theta proposal mean of %d\n") % Forest::_theta_proposal_mean;
                Forest::_theta_prior_mean = Forest::_theta_proposal_mean;
            }
            // no proposal or prior mean if theta fixed
            else if (Forest::_theta_prior_mean > 0.0 && Forest::_theta_proposal_mean ==  0.0) {
                cout << boost::format("\nSetting theta proposal mean equal to theta prior mean of %d\n") % Forest::_theta_prior_mean;
                Forest::_theta_proposal_mean = Forest::_theta_prior_mean;
            }
            else if (Forest::_theta_prior_mean == 0.0 && Forest::_theta_proposal_mean == 0.0) {
                cout << boost::format("\nTheta mean of %d will be fixed for all particles; population sizes will all be drawn from the same theta\n") % Forest::_theta;
            }
        }
    }

    inline void Proj::checkOutgroupName() {
        bool found = false;
        for (auto &s:_taxon_map) {
            if (Forest::_outgroup == s.second) {
                found = true;
            }
        }
        if (!found) {
            throw XProj(format("outgroup name does not match any species name"));
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
        if (fabs(sum-1)>0.000001) {
            throw XProj(format("base frequencies (%s) don't add to 1")%Forest::_string_base_frequencies);
        }
        assert (fabs(sum-1) < 0.000001);
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

    inline void Proj::proposeSpeciesGroups(vector<Bundle> bundles, unsigned ngroups, unsigned nsubsets, unsigned ntaxa, string filename1, string filename2) {
        // ngroups = number of species SMCs to do (i.e. 100 particles for first round, thin = 1.0 means ngroups = 100 for this round)
          assert (_nthreads > 1);
          
          // divide up groups as evenly as possible across threads
          unsigned first = 0;
          unsigned incr = ngroups/_nthreads + (ngroups % _nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of groups
          unsigned last = incr;
          
          // need a vector of threads because we have to wait for each one to finish
          vector<thread> threads;

            while (true) {
            // create a thread to handle particles first through last - 1
              threads.push_back(thread(&Proj::proposeSpeciesGroupRange, this, first, last, std::ref(bundles), ngroups, filename1, nsubsets, ntaxa, filename2));
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

    inline void Proj::proposeSpeciesGroupRange(unsigned first, unsigned last, vector<Bundle> &bundles, unsigned ngroups, string filename1, unsigned nsubsets, unsigned ntaxa, string filename2) {
        
        unsigned nspecies = (unsigned) _species_names.size();
           
           for (unsigned i=first; i<last; i++){
               vector<Bundle> use_vec;
               
               Bundle chosen_bundle = bundles[i]; // bundle to copy
               
               for (unsigned a=0; a<_particle_increase; a++) {
                   use_vec.push_back(*new Bundle(chosen_bundle));
               }

               assert(use_vec.size() == _particle_increase);

               if (_verbose > 0) {
                   cout << "beginning species tree proposals for subset " << i+1 << endl;
               }
               for (unsigned s = 0; s < nspecies-1; s++) {  // skip last round of filtering because weights are always 0

                   // set particle random number seeds
                   unsigned psuffix = 1;
                   for (auto &b:use_vec) {
                       b.setSeed(rng.randint(1,9999) + psuffix);
                       psuffix += 2;
                   }

                   proposeSpeciesParticles(use_vec);
                   filterBundles(s, use_vec);
                   resetWeights(use_vec);


               } // s loop
               
               if (_verbose > 0) {
                   cout << "finished with species tree proposals for subset " << i+1 << endl;
               }
               
               if (_save_every > 1.0) { // thin sample for output by taking a random sample
                   unsigned sample_size = round (double (_particle_increase) / double(_save_every));
                   if (sample_size == 0) {
                       cout << "\n";
                       cout << "current settings would save 0 species trees; saving every species tree\n";
                       cout << "\n";
                       sample_size = _particle_increase;
                   }

                   random_shuffle(use_vec.begin(), use_vec.end()); // shuffle particles, random_shuffle will always shuffle in same order
                   // delete first (1-_thin) % of particles
                   use_vec.erase(next(use_vec.begin(), 0), next(use_vec.begin(), (_particle_increase-sample_size)));
                   assert (use_vec.size() == sample_size);
               }

               mtx.lock(); // TODO: does this slow things down?
               saveSpeciesTreesHierarchical(use_vec, filename1);
               writeParamsFileForBeastComparisonAfterSpeciesFiltering(nsubsets, nspecies, ntaxa, use_vec, filename2, i);
               
               _count++;
               mtx.unlock();
           }
    }

    inline void Proj::proposeSpeciesParticles(vector<Bundle> &bundles) {
        assert(_nthreads > 0);
       if (_nthreads == 1) {
         for (auto & b : bundles) {
             b.proposeSpeciesParticles();
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
             threads.push_back(thread(&Proj::proposeSpeciesParticleRange, this, first, last, std::ref(bundles)));
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

    inline void Proj::proposeSpeciesParticleRange(unsigned first, unsigned last, vector<Bundle> &bundles) {
        for (unsigned i=first; i<last; i++){
              bundles[i].proposeSpeciesParticles();
          }
    }

    inline void Proj::runBundlesRange(unsigned first, unsigned last, vector<Bundle> &bundles) {
       for (unsigned i=first; i<last; i++){
           bundles[i].runBundle();
       }
   }

    inline void Proj::runBundles(vector<Bundle> &bundles) {
        assert(_nthreads > 0);
        if (_nthreads == 1) {
          for (auto & b : bundles) {
              b.runBundle();
          }
        }
        else {
          // divide up the particles as evenly as possible across threads
          unsigned first = 0;
          unsigned incr = _nbundles/_nthreads + (_nbundles % _nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of particles
          unsigned last = incr;

          // need a vector of threads because we have to wait for each one to finish
          vector<thread> threads;

            while (true) {
            // create a thread to handle particles first through last - 1
              threads.push_back(thread(&Proj::runBundlesRange, this, first, last, std::ref(bundles)));
            // update first and last
            first = last;
            last += incr;
            if (last > _nbundles) {
              last = _nbundles;
              }
            if (first>=_nbundles) {
                break;
            }
          }

          // the join function causes this loop to pause until the ith thread finishes
          for (unsigned i = 0; i < threads.size(); i++) {
            threads[i].join();
          }
        }
    }

    inline void Proj::filterBundles(unsigned step, vector<Bundle> & bundles) {
        // Copy log weights for all bundles to prob vector
        vector<double> probs(_nbundles, 0.0);
        
        for (unsigned b=0; b < _nbundles; b++) {
            probs[b] = bundles[b].getBundleLogWeight();
        }

        // Normalize log_weights to create discrete probability distribution
        double log_sum_weights = getRunningSum(probs);
        
        transform(probs.begin(), probs.end(), probs.begin(), [log_sum_weights](double logw){return exp(logw - log_sum_weights);});

        // Compute cumulative probabilities
        partial_sum(probs.begin(), probs.end(), probs.begin());

        // Initialize vector of counts storing number of darts hitting each particle
        vector<unsigned> counts (_nbundles, 0);

        // Throw _nparticles darts
        for (unsigned i=0; i<_nbundles; i++) {
            double u = rng.uniform();
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

        // Count number of cells with zero count that can serve as copy recipients
        unsigned nzeros = 0;
        for (unsigned i = 0; i < _nbundles; i++) {
            if (counts[i] == 0)
                nzeros++;
        }

        while (nzeros > 0) {
            assert(donor < _nbundles);
            assert(recipient < _nbundles);

            // Copy donor to recipient
            bundles[recipient] = bundles[donor];

            counts[donor]--;
            counts[recipient]++;
            nzeros--;

            if (counts[donor] == 1) {
                // Move donor to next slot with count > 1
                donor++;
                while (donor < _nbundles && counts[donor] < 2) {
                    donor++;
                }
            }

            // Move recipient to next slot with count equal to 0
            recipient++;
            while (recipient < _nbundles && counts[recipient] > 0) {
                recipient++;
            }
        }
    }

    inline void Proj::resetWeights(vector<Bundle> & bundles) {
        double logw = -log(bundles.size());
        for (auto & b : bundles) {
            b.setBundleLogWeight(logw);
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

    inline void Proj::showFinal(vector<Particle::SharedPtr> my_vec) {
        for (auto &p:my_vec){
            p->showParticle();
        }

        double sum_h = 0.0;
        for (auto & p:my_vec) {
            double h = p->calcHeight();
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

        cout << "speciation rate = " << Forest::_lambda << endl;
    }

    inline void Proj::initializeBundles(vector<Bundle> &bundles) {
        // set partials for first particle under save_memory setting for initial marginal likelihood calculation
        assert (_nthreads > 0);

        bool partials = false;
        if (_gene_newicks_specified) {
            partials = false;
            Forest::_save_memory = true;
        }

        for (auto & b:bundles ) { // TODO: can initialize some of these things in parallel?
            b.setNGeneParticles(_nparticles);
            b.setData(_data, _taxon_map, partials);
            partials = false;
            b.mapSpecies(_taxon_map, _species_names);
            // TODO: calculate likelihood for one particle, then delete those partials, set likelihood / weights for every other particle
            // TODO: can start weights at 0 because every gene will get changed?
        }
    }

    inline void Proj::debugSpeciesTree(vector<Particle::SharedPtr> &particles) {
        cout << "debugging species tree" << endl;
        for (auto &p:particles) {
            p->showSpeciesJoined();
            p->showSpeciesIncrement();
            p->showSpeciesTree();
            cout << " _______ " << endl;
        }
    }

    inline void Proj::run() {
            _first_line = true;
            if (_verbose > 0) {
                cout << "Starting..." << endl;
                cout << "Current working directory: " << boost::filesystem::current_path() << endl;
                cout << "Random seed: " << _random_seed << endl;
                cout << "Theta: " << Forest::_theta << endl;
                cout << "Number of threads: " << _nthreads << endl;
#if defined (DRAW_NEW_THETA)
                cout << "drawing new theta for each particle " << endl;
#else
                cout << "Theta: " << Forest::_theta << endl;
#endif
            }

            try {
                if (_verbose > 0) {
                    cout << "\n*** Reading and storing the data in the file " << _data_file_name << endl;
                    cout << "data file name is " << _data_file_name << endl;
                }
                
                _data = Data::SharedPtr(new Data());
                _data->setPartition(_partition);
                _data->getDataFromFile(_data_file_name);

                if (_verbose > 0) {
                    summarizeData(_data);
                }
                createSpeciesMap(_data);

                // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
                if (Forest::_outgroup != "none") {
                    checkOutgroupName();
                }

                //set number of species to number in data file
                unsigned ntaxa = setNumberTaxa(_data);
                unsigned nspecies = (unsigned) _species_names.size();
                Forest::setNumSpecies(nspecies);
                rng.setSeed(_random_seed);

    //          create vector of bundles
                unsigned nsubsets = _data->getNumSubsets();
                
                vector<Bundle> bundle_vec;
                bundle_vec.resize(_nbundles);

                Bundle::setNumGenes(nsubsets);

                
                initializeBundles(bundle_vec); // initialize in parallel with multithreading
                
#if defined (DRAW_NEW_THETA)
                unsigned psuffix = 1;
                for (auto &b:bundle_vec) {
                    b.setSeed(rng.randint(1,9999) + psuffix);
                    psuffix += 2;
                }
                
                for (auto &b:bundle_vec) {
                    b.drawTheta();
                }
#endif
                
                unsigned nsteps = ntaxa - 1;
                
                for (unsigned s=0; s<nsteps; s++) {
                    cout << "beginning step " << s << endl;
                    
                    unsigned psuffix = 1;
                    for (auto &b:bundle_vec) {
                        // set bundle random number seeds
                        b.setSeed(rng.randint(1,9999) + psuffix);
                        psuffix += 2;
                    }
                    
                    runBundles(bundle_vec);
                    
                    filterBundles(s, bundle_vec);
                    resetWeights(bundle_vec);
                
                }
                                
                saveSpeciesTrees(bundle_vec);
                if (_save_gene_trees) {
                    for (unsigned i=0; i<nsubsets; i++) {
                        saveGeneTree(i, bundle_vec); // TODO: for now, saving 1 gene per bundle
                    }
                }
                
                writeParamsFileForBeastComparison(nsubsets, nspecies, ntaxa, bundle_vec);

#if defined (HIERARCHICAL_FILTERING)
                string filename1 = "species_trees.trees";
                string filename2 = "params-beast-comparison.log";
                
                if (filesystem::remove(filename1)) {
                    ofstream speciestrf(filename1);
                    speciestrf << "#nexus\n\n";
                    speciestrf << "begin trees;\n";
                    if (_verbose > 0) {
                        cout << "existing file " << filename1 << " removed and replaced\n";
                    }
                }
                else {
                    ofstream speciestrf(filename1);
                    speciestrf << "#nexus\n\n";
                    speciestrf << "begin trees;\n";
                    if (_verbose > 0) {
                        cout << "created new file " << filename1 << "\n";
                    }
                }
                
                unsigned ngroups = round(_nbundles * _thin);
                if (ngroups == 0) {
                    ngroups = 1;
                    cout << "thin setting would result in 0 species groups; setting species groups to 1" << endl;
                }

                random_shuffle(bundle_vec.begin(), bundle_vec.end()); // shuffle particles, random_shuffle will always shuffle in same order
                // delete first (1-_thin) % of particles
                bundle_vec.erase(next(bundle_vec.begin(), 0), next(bundle_vec.begin(), (_nbundles-ngroups)));
                assert(bundle_vec.size() == ngroups);


                for (auto &b:bundle_vec) {
                    b.deleteExtraGeneParticles(); // delete all of the gene particles except the first one since we will only condition on one set per bundle
                    // reset forest species partitions
                    b.clearPartials(); // no more likelihood calculations
                    b.resetSpecies();
                    b.mapSpecies(_taxon_map, _species_names);
                }

                vector<Particle::SharedPtr> new_vec;

                _nparticles = _particle_increase;
                _nbundles = _particle_increase; // TODO: this might be necessary?
                unsigned index = 0;
                
                bool parallelize_by_group = false;
                
                if (_nthreads > 1) {
        #if defined PARALLELIZE_BY_GROUP
                    parallelize_by_group = true;
        #endif
                }
                
                if (parallelize_by_group) {
                    // don't bother with this if not multithreading
                        proposeSpeciesGroups(bundle_vec, ngroups, nsubsets, ntaxa, filename1, filename2);
                    
                    // formatting output files
                        ofstream strees;
                        strees.open("species_trees.trees", std::ios::app);
                        strees << "end;" << endl;
                        strees.close();

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
                                    string add = to_string(line_count) + "\t";
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
                    for (unsigned a=0; a < ngroups; a++) {
                        // TODO: for now, take 1 set of gene trees from each bundle
                        vector<Bundle> use_vec;
                        
                        Bundle chosen_bundle = bundle_vec[a]; // bundle to copy
                        for (unsigned i=0; i<_particle_increase; i++) {
                            use_vec.push_back(*new Bundle(chosen_bundle));
                        }

                        assert(use_vec.size() == _particle_increase);

                        index += _particle_increase;

                        if (_verbose > 0) {
                            cout << "beginning species tree proposals for subset " << a+1 << endl;
                        }
                        
                        for (unsigned s=0; s<nspecies-1; s++) {  // skip last round of filtering because weights are always 0
                            if (_verbose > 0) {
                                cout << "starting species step " << s+1 << " of " << nspecies-1 << endl;
                            }

                            // set particle random number seeds
                            unsigned psuffix = 1;
                            for (auto &b:use_vec) {
                                b.setSeed(rng.randint(1,9999) + psuffix);
                                psuffix += 2;
                            }
                            proposeSpeciesParticles(use_vec);
                            
                            filterBundles(s, use_vec);
                            resetWeights(use_vec);
                        } // s loop
                        
                        if (_save_every > 1.0) { // thin sample for output by taking a random sample
                            unsigned sample_size = round (double (_particle_increase) / double(_save_every));
                            if (sample_size == 0) {
                                cout << "\n";
                                cout << "current settings would save 0 species trees; saving every species tree\n";
                                cout << "\n";
                                sample_size = _particle_increase;
                            }

                            random_shuffle(use_vec.begin(), use_vec.end()); // shuffle particles, random_shuffle will always shuffle in same order
                            // delete first (1-_thin) % of particles
                            use_vec.erase(next(use_vec.begin(), 0), next(use_vec.begin(), (_particle_increase-sample_size)));
                            assert (use_vec.size() == sample_size);
                            
                        }
                        saveSpeciesTreesHierarchical(use_vec, filename1);
                        writeParamsFileForBeastComparisonAfterSpeciesFiltering(nsubsets, nspecies, ntaxa, use_vec, filename2, a);
                    }

                    // finish species tree file
                    std::ofstream treef;
                    treef.open(filename1, std::ios_base::app);
                    treef << "end;\n";
                    treef.close();
                }

#endif

            }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }
}
