#pragma once

#include <iostream>
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
extern proj::Lot rng;

namespace proj {

    class Proj {
        public:
        
                                Proj();
                                ~Proj();

            void                clear();
            void                processCommandLineOptions(int argc, const char * argv[]);
            void                run();
//            void                saveAllForests(const vector<Particle> &v) const ;
            void                saveAllForests(vector<Particle::SharedPtr> &v) const ;
            void                saveSpeciesTrees(vector<Particle::SharedPtr> &v) const;
            void                saveGeneTrees(unsigned ngenes, vector<Particle::SharedPtr> &v) const;
            void                saveGeneTree(unsigned gene_number, vector<Particle::SharedPtr> &v) const;
            void                writeLoradFile(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle::SharedPtr> &v) const;
            void                writeDeepCoalescenceFile(vector<Particle::SharedPtr> &v);
            void writeParamsFileForBeastComparison (unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle::SharedPtr> &v) const;
            void                normalizeWeights(vector<Particle::SharedPtr> & particles);
            void                normalizeSpeciesWeights(vector<Particle::SharedPtr> & particles);
            void                modifyWeights(vector<Particle::SharedPtr> & particles);
            void                correctWeights(vector<Particle::SharedPtr> & particles);
            void                resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
            void                resampleSpeciesParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
            void                resetWeights(vector<Particle::SharedPtr> & particles);
            void                resetSpeciesWeights(vector<Particle::SharedPtr> & particles);
            double              getWeightAverage(vector<double> log_weight_vec);
            void                createSpeciesMap(Data::SharedPtr);
            void                simSpeciesMap();
            string              inventName(unsigned k, bool lower_case);
            void                showFinal(vector<Particle::SharedPtr>);
            void                proposeSpeciesParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles);
            void                proposeSpeciesParticles(vector<Particle::SharedPtr> &particles);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles);
            void                proposeParticles(vector<Particle::SharedPtr> &particles);
            void                proposeParticleRangePriorPost(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles);
            void                proposeParticlesPriorPost(vector<Particle::SharedPtr> &particles);
            void                saveAllHybridNodes(vector<Particle::SharedPtr> &v) const;
            void                simulateData();
            void                writePaupFile(vector<Particle::SharedPtr> particles, vector<string> taxpartition);
            void                sanityChecks();

        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood = 0.0;
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
            void                        handleNTaxaPerSpecies();
            void                        checkOutgroupName();
            void                        debugSpeciesTree(vector<Particle::SharedPtr> &particles);
            double                      _small_enough;
            unsigned                    _verbose;
            double                      _phi;
            unsigned                    _sim_nspecies;
            vector<unsigned>            _ntaxaperspecies;
            string                      _string_ntaxaperspecies;
            string                      _sim_file_name;
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
        _phi = 1.0;
    }

//    inline void Proj::saveAllForests(const vector<Particle> &v) const {
inline void Proj::saveAllForests(vector<Particle::SharedPtr> &v) const {
        ofstream treef("forest.trees");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:v) {
            treef << "  tree test = [&R] " << p->saveForestNewick()  << ";\n";
        }
        treef << "end;\n";
        treef.close();
    }

    inline void Proj::writeParamsFileForBeastComparison(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle::SharedPtr> &v) const {
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
            
            double log_likelihood = p->getLogLikelihood();
            double log_prior = p->getAllPriors();
            
            double log_posterior = log_likelihood + log_prior;
            
            double log_coalescent_likelihood = 0.0;
            for (unsigned g=1; g<ngenes+1; g++) {
                log_coalescent_likelihood += p->getCoalescentLikelihood(g);
            }
            
            logf << "\t" << log_posterior;
            
            logf << "\t" << log_likelihood;
            
            logf << "\t" << log_prior - log_coalescent_likelihood; // starbeast3 does not include coalescent likelihood in this prior
            
            double vector_prior = 0.0;
            logf << "\t" << vector_prior; // log joint prior on population sizes (vectorPrior) - 0.0 for now since pop size is not a parameter
            
            logf << "\t" << log_coalescent_likelihood;
            
            double species_tree_height = p->getSpeciesTreeHeight();
            logf << "\t" << species_tree_height;
            
            double species_tree_length = p->getSpeciesTreeLength();
            logf << "\t" << species_tree_length;
            
            vector<double> gene_tree_heights = p->getGeneTreeHeights();
            vector<double> gene_tree_lengths = p->getGeneTreeLengths();
            assert (gene_tree_heights.size() == gene_tree_lengths.size());
            
            for (int i=0; i<gene_tree_heights.size(); i++) {
                logf << "\t" << gene_tree_heights[i];
                logf << "\t" << gene_tree_lengths[i];
            }
            
            double yule_model = p->getSpeciesTreePrior(); // TODO: unsure if this is correct
            logf << "\t" << yule_model;
            
            double pop_mean = 0.0; // setting this to 0.0 for now since pop size is not a parameter
            logf << "\t" << pop_mean;
            
            for (int i=0; i<(nspecies*2-1); i++) {
                logf << "\t" << Forest::_theta / 4.0; // all pop sizes are the same under this model, Ne*u = theta / 4?
//                logf << "\t" << p->getNewTheta() / 4.0;
            }
            
            logf << "\t" << Forest::_lambda;
            
            vector<double> gene_tree_log_likelihoods = p->getGeneTreeLogLikelihoods();
            vector<double> gene_tree_priors = p->getGeneTreePriors();
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

    inline void Proj::writeDeepCoalescenceFile(vector<Particle::SharedPtr> &v) {
        ofstream logf("deep_coalescences.txt");
        logf << "particle ";
        logf << "\t" << "num deep coalescences ";
        unsigned count = 0;
        for (auto &p:v) {
            logf << "\n" << count;
            logf << "\t" << p->getNumDeepCoalescences();
            count++;
        }
    }

    inline void Proj::writeLoradFile(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle::SharedPtr> &v) const {
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
            
#if defined (WEIGHT_CORRECTION)
            logf << "\t" << p->getLogLikelihood()*_phi;
#else
            logf << "\t" << p->getLogLikelihood();
#endif

            double species_tree_prior = 0.0;
            
            for (unsigned g=0; g<ngenes+1; g++) {
                for (auto &b:p->getIncrementPriors(g)) {
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
                log_coalescent_likelihood += p->getCoalescentLikelihood(g);
            }
            logf << "\t" << log_coalescent_likelihood;
            
            logf << endl;
        }
        
        logf.close();
    }

    inline void Proj::saveSpeciesTrees(vector<Particle::SharedPtr> &v) const {
        // save only unique species trees
        if (!Forest::_run_on_empty) {
            vector<vector<pair<double, double>>> unique_increments_and_priors;
            
            ofstream unique_treef("unique_species_trees.trees");
            unique_treef << "#nexus\n\n";
            unique_treef << "begin trees;\n";
            for (auto &p:v) {
                vector<pair<double, double>> increments_and_priors = p->getSpeciesTreeIncrementPriors();
                bool found = false;
                if(std::find(unique_increments_and_priors.begin(), unique_increments_and_priors.end(), increments_and_priors) != unique_increments_and_priors.end()) {
                    found = true;
                }
                if (!found) {
                    unique_increments_and_priors.push_back(increments_and_priors);
                    unique_treef << "  tree test = [&R] " << p->saveForestNewick()  << ";\n";
                }
            }
            unique_treef << "end;\n";
            unique_treef.close();
        }
        
        if (_start_mode == "smc") {
            // save all species trees
            ofstream treef("species_trees.trees");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:v) {
                treef << "  tree test = [&R] " << p->saveForestNewick()  << ";\n";
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
                treef << "  tree test = [&R] " << p->saveForestNewick()  << ";\n";
            }
            treef << "end;\n";
            treef.close();
        }
        
        }

    inline void Proj::saveGeneTrees(unsigned ngenes, vector<Particle::SharedPtr> &v) const {
        if (_start_mode == "smc") {
            ofstream treef("gene_trees.trees");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:v) {
                    for (int i=1; i<ngenes+1; i++) {
                        treef << "tree gene" << i << " = [&R] " << p->saveGeneNewick(i)  << ";\n";
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
            for (auto &p:v) {
                    for (int i=1; i<ngenes+1; i++) {
                        treef << "tree gene" << i << " = [&R] " << p->saveGeneNewick(i)  << ";\n";
                }
                treef << endl;
            }
            treef << "end;\n";
            treef.close();
        }
    }

    inline void Proj::saveGeneTree (unsigned gene_number, vector<Particle::SharedPtr> &v) const {
        string name = "gene" + to_string(gene_number) + ".trees";
        ofstream treef(name);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:v) {
            treef << "  tree test = [&R] " << p->saveGeneNewick(gene_number)  << ";\n";
            treef << endl;
        }
        treef << "end;\n";
        
        treef.close();
    }

    inline void Proj::saveAllHybridNodes(vector<Particle::SharedPtr> &v) const {
        ofstream nodef("nodes.txt");
        for (auto &p:v) {
            nodef << "particle\n";
            nodef << p->saveHybridNodes()  << "\n";
        }
        nodef.close();
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
        ("theta, t", boost::program_options::value(&Forest::_theta)->default_value(0.05), "theta")
        ("lambda", boost::program_options::value(&Forest::_lambda)->default_value(1), "speciation rate")
        ("proposal",  boost::program_options::value(&Forest::_proposal)->default_value("prior-prior"), "a string defining a proposal (prior-prior or prior-post)")
        ("model", boost::program_options::value(&Forest::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&Forest::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("nthreads",  boost::program_options::value(&_nthreads)->default_value(1), "number of threads for multi threading")
        ("migration_rate", boost::program_options::value(&Forest::_migration_rate)->default_value(0.0), "migration rate")
        ("hybridization_rate", boost::program_options::value(&Forest::_hybridization_rate)->default_value(0.0), "hybridization rate")
        ("run_on_empty", boost::program_options::value(&Particle::_run_on_empty)->default_value(false), "run program without data")
        ("verbose", boost::program_options::value(&_verbose)->default_value(1), "set amount of output printed")
        ("save_memory", boost::program_options::value(&Forest::_save_memory)->default_value(false), "save memory at the expense of time")
        ("phi", boost::program_options::value(&_phi)->default_value(1.0), "correct weights by this number")
        ("outgroup", boost::program_options::value(&Forest::_outgroup)->default_value("none"), "a string defining the outgroup")
        ("startmode", boost::program_options::value(&_start_mode)->default_value("smc"), "a string defining whether to simulate data or perform smc")
        ("nspecies", boost::program_options::value(&_sim_nspecies)->default_value(1), "number of species to simulate")
        ("ntaxaperspecies", boost::program_options::value(&_string_ntaxaperspecies)->default_value("1"), "number of taxa per species to simulate")
        ("filename", boost::program_options::value(&_sim_file_name), "name of file to write simulated data to")
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
        // if user specified "ntaxaperspecies" in conf file, convert them to a vector<unsigned>
        if (vm.count("base_frequencies") > 0) {
            handleNTaxaPerSpecies();
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
    
    inline void Proj::handleNTaxaPerSpecies() {
        vector<string> temp;
        split(temp, _string_ntaxaperspecies, is_any_of(","));
        // iterate throgh temp
        for (auto &i:temp) {
            double f = stof(i);
            _ntaxaperspecies.push_back(f);
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

    inline void Proj::modifyWeights(vector<Particle::SharedPtr> & particles) {
        for (auto &p:particles) {
            p->setLogWeight(p->getLogWeight()*_phi);
        }
    }

    inline void Proj::correctWeights(vector<Particle::SharedPtr> & particles) {
        for (auto &p:particles) {
            p->setLogWeight(p->getLogWeight()*(1 / _phi));
        }
    }

    inline void Proj::normalizeSpeciesWeights(vector<Particle::SharedPtr> & particles) {
        unsigned i = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            log_weight_vec[i++] = p->getSpeciesLogWeight();
        }

        double log_particle_sum = getRunningSum(log_weight_vec);
        
        double max = *max_element(std::begin(log_weight_vec), std::end(log_weight_vec));
        double min = *min_element(std::begin(log_weight_vec), std::end(log_weight_vec)); // C++11
        
        if (_verbose > 1) {
            cout << "\t" << "max weight = " << max << endl;;
            cout << "\t" << "min weight = " << min << endl;;
        }
//        for (auto &p:particles) {
//            p->showParticle();
//        }

        for (auto & p : particles) {
            p->setLogSpeciesWeight(p->getSpeciesLogWeight() - log_particle_sum);
        }
        
//        _log_marginal_likelihood += log_particle_sum - log(_nparticles);
    }

    inline void Proj::normalizeWeights(vector<Particle::SharedPtr> & particles) {
        unsigned i = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            log_weight_vec[i++] = p->getLogWeight();
        }

        double log_particle_sum = getRunningSum(log_weight_vec);
        
        double max = *max_element(std::begin(log_weight_vec), std::end(log_weight_vec));
        double min = *min_element(std::begin(log_weight_vec), std::end(log_weight_vec)); // C++11
        
        if (_verbose > 1) {
            cout << "\t" << "max weight = " << max << endl;;
            cout << "\t" << "min weight = " << min << endl;;
        }

        for (auto & p : particles) {
            p->setLogWeight(p->getLogWeight() - log_particle_sum);
        }
        
        _log_marginal_likelihood += log_particle_sum - log(_nparticles);
    }

    inline void Proj::resampleSpeciesParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles) {
        assert (from_particles.size() == _nparticles);
        assert (to_particles.size() == _nparticles);
        
        vector<pair<double, double>> cum_probs;
            // Create vector of pairs p, with p.first = log weight and p.second = particle index
        cum_probs.resize(_nparticles);
        unsigned i = 0;
        
        for (unsigned p=0; p < _nparticles; p++) {
            cum_probs[i].first = from_particles[p]->getSpeciesLogWeight();
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
            
            assert(_nparticles == to_particles.size());
            }
    }

    inline void Proj::resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles) {
         assert (from_particles.size() == _nparticles);
         assert (to_particles.size() == _nparticles);
         
         vector<pair<double, double>> cum_probs;
             // Create vector of pairs p, with p.first = log weight and p.second = particle index
         cum_probs.resize(_nparticles);
         unsigned i = 0;
         
         for (unsigned p=0; p < _nparticles; p++) {
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
             assert( fabs( 1.0 - cum_probs[_nparticles-1].first ) < Proj::_small_enough);

         
         // Draw new set of particles by sampling with replacement according to cum_probs
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
             
             assert(_nparticles == to_particles.size());
             }
        }

    inline void Proj::resetWeights(vector<Particle::SharedPtr> & particles) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p->setLogWeight(logw);
        }
    }

    inline void Proj::resetSpeciesWeights(vector<Particle::SharedPtr> & particles) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p->setLogSpeciesWeight(logw);
        }
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
        for (int s=0; s<_sim_nspecies; s++) {
            string species_name;
            species_name = inventName(s, false);
            for (int t=0; t<_ntaxaperspecies[s]; t++) {
                string taxon_name = inventName(count, true) + "^" + species_name;
                _taxon_map.insert({taxon_name, species_name});
                count++;
            }
            _species_names.push_back(species_name);
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
        cout << "theta = " << Forest::_theta << endl;
        cout << "speciation rate = " << Forest::_lambda << endl;
    }

    inline void Proj::proposeSpeciesParticles(vector<Particle::SharedPtr> &particles) {
        assert(_nthreads > 0);
        if (_nthreads == 1) {
          for (auto & p : particles) {
              p->speciesOnlyProposal();
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
              threads.push_back(thread(&Proj::proposeSpeciesParticleRange, this, first, last, std::ref(particles)));
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

    inline void Proj::proposeParticles(vector<Particle::SharedPtr> &particles) {
        assert(_nthreads > 0);
        if (_nthreads == 1) {
          for (auto & p : particles) {
              p->proposal();
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
              threads.push_back(thread(&Proj::proposeParticleRange, this, first, last, std::ref(particles)));
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

    inline void Proj::proposeParticlesPriorPost(vector<Particle::SharedPtr> &particles) {
        assert(_nthreads > 0);
        if (_nthreads == 1) {
          for (auto & p : particles) {
              p->proposalPriorPost();
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
              threads.push_back(thread(&Proj::proposeParticleRangePriorPost, this, first, last, std::ref(particles)));
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

    inline void Proj::proposeParticleRangePriorPost(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles) {
        for (unsigned i=first; i<last; i++){
            particles[i]->proposalPriorPost();
        }
    }

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles) {
        for (unsigned i=first; i<last; i++){
            particles[i]->proposal();
        }
    }

    inline void Proj::proposeSpeciesParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles) {
        for (unsigned i=first; i<last; i++){
            particles[i]->speciesOnlyProposal();
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

    inline void Proj::writePaupFile(vector<Particle::SharedPtr> particles, vector<string> taxpartition) {
        // Output a PAUP* command file for estimating the species tree using
        // svd quartets and qage
        cout << "  PAUP* commands saved in file \"svd-qage.nex\"\n";
        ofstream paupf("svd-qage.nex");
        paupf << "#NEXUS\n\n";
        paupf << "begin paup;\n";
        paupf << "  log start file=svdout.txt replace;\n";
        paupf << "  exe " + _sim_file_name + ";\n";
        paupf << "  taxpartition species (vector) = " << join(taxpartition," ") << ";\n";
        paupf << "  svd taxpartition=species;\n";
        paupf << "  roottrees;\n";
        paupf << "  qage taxpartition=species patprob=exactjc outUnits=substitutions treefile=svd.tre replace;\n";
        paupf << "  log stop;\n";
        paupf << "  quit;\n";
        paupf << "end;\n";
        paupf.close();
    }

    inline void Proj::simulateData() {
        cout << "\nSimulating data under multispecies coalescent model...\n" << endl;
        rng.setSeed(_random_seed);
        
        if (_ntaxaperspecies.size() == 1) {
            unsigned ntaxa = _ntaxaperspecies[0];
            for (int i=0; i<_sim_nspecies-1; i++) {
                _ntaxaperspecies.push_back(ntaxa);
            }
        }
        
        vector<Particle::SharedPtr> sim_vec(1);
        sim_vec[0] = Particle::SharedPtr(new Particle);
        
        // set particle randon number seed
        unsigned psuffix = 1;
        sim_vec[0]->setSeed(rng.randint(1,9999) + psuffix);
        psuffix += 2;

        Particle::_run_on_empty = true;
        Forest::_run_on_empty = true;
        Forest::_proposal = "prior-prior";
        
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        
        // make up the species map
        simSpeciesMap();
        
        vector<string> taxpartition;
        for (auto &t:_taxon_map) {
            taxpartition.push_back(t.second);
        }
        
        // if user specified an outgroup in conf file, check that the outgroup matches one of the species names
        if (Forest::_outgroup != "none") {
            checkOutgroupName();
        }
        
        unsigned nsubsets = _data->getNumSubsets();
        Particle::setNumSubsets(nsubsets);
        
        sim_vec[0]->setSimData(_data, _taxon_map, nsubsets, (unsigned) _taxon_map.size());
        
        sim_vec[0]->mapSpecies(_taxon_map, _species_names);

        unsigned nsteps = (unsigned) (_taxon_map.size()-1)*nsubsets;
        
        for (unsigned g=0; g<nsteps; g++){
            proposeParticles(sim_vec);
        }
        
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
        
        sim_vec[0]->simulateData(sites_vector);

        _data->compressPatterns();
        _data->writeDataToFile(_sim_file_name);
        
        saveGeneTrees(nsubsets, sim_vec);
        saveSpeciesTrees(sim_vec);
        
        writePaupFile(sim_vec, taxpartition);
        writeDeepCoalescenceFile(sim_vec);
    }

    inline void Proj::sanityChecks() {
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
            if (_ntaxaperspecies.size() != 1 && _ntaxaperspecies.size() != _sim_nspecies) {
                throw XProj("must specify number of taxa per species or one number if equal number of taxa per species");
            }
        }
        else {
            if (_data_file_name == "") {
                throw XProj("must specify name of data file if smc option is chosen");
            }
        }
    }

    inline void Proj::run() {
        sanityChecks();
        if (_start_mode == "sim") {
            try {
                simulateData();
            }
            catch (XProj & x) {
                std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
            }
        }
        else {
            if (_verbose > 0) {
                cout << "Starting..." << endl;
                cout << "Current working directory: " << boost::filesystem::current_path() << endl;
                cout << "Random seed: " << _random_seed << endl;
                cout << "Theta: " << Forest::_theta << endl;
                cout << "Number of threads: " << _nthreads << endl;
            }
            
            if (Particle::_run_on_empty) { // if running with no data, choose taxa to join at random
                Forest::_proposal = "prior-prior";
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

    //          create vector of particles
                unsigned nparticles = _nparticles;
                
                unsigned nsubsets = _data->getNumSubsets();
                Particle::setNumSubsets(nsubsets);
                
                vector<Particle::SharedPtr> my_vec_1(nparticles);
                vector<Particle::SharedPtr> my_vec_2(nparticles);
                vector<Particle::SharedPtr> &my_vec = my_vec_1;
                    
                for (unsigned i=0; i<nparticles; i++) {
                    my_vec_1[i] = Particle::SharedPtr(new Particle);
                    my_vec_2[i] = Particle::SharedPtr(new Particle);
                }

                bool use_first = true;
                bool partials = true;
                for (auto & p:my_vec ) {
                    p->setData(_data, _taxon_map, partials);
                    partials = false;
                    p->mapSpecies(_taxon_map, _species_names);
                }
                    
                // reset marginal likelihood
                _log_marginal_likelihood = 0.0;
                vector<double> starting_log_likelihoods = my_vec[0]->calcGeneTreeLogLikelihoods();
            
                for (auto &p:my_vec) {
                    p->setLogLikelihood(starting_log_likelihoods);
                    if (Forest::_save_memory) {
                        p->clearPartials();
                    }
                }
                    
        #if defined (WEIGHT_CORRECTION)
                            modifyWeights(my_vec);
        #endif
                        normalizeWeights(my_vec); // initialize marginal likelihood
                        
                        //run through each generation of particles
                    
                        unsigned nsteps = (ntaxa-1)*nsubsets;
                        
                        for (unsigned g=0; g<nsteps; g++){
                            if (_verbose > 0) {
                                cout << "starting step " << g << " of " << nsteps-1 << endl;
                            }
                            // set particle randon number seeds
                            unsigned psuffix = 1;
                            for (auto &p:my_vec) {
                                p->setSeed(rng.randint(1,9999) + psuffix);
                                psuffix += 2;
                            }
                            
                            //taxon joining and reweighting step
                            if (Forest::_proposal == "prior-prior") {
                                proposeParticles(my_vec);
                            }
                            else {
                                proposeParticlesPriorPost(my_vec);
                            }
                            
                            unsigned num_species_particles_proposed = 0;
                            
                            double ess_inverse = 0.0;
                            for (auto &p:my_vec) {
                                if (p->speciesJoinProposed()) {
                                    num_species_particles_proposed++;
                                }
                            }
                            
        #if defined (WEIGHT_CORRECTION)
                            modifyWeights(my_vec);
        #endif
                            normalizeWeights(my_vec);
                            
                            for (auto & p:my_vec) {
                                ess_inverse += exp(2.0*p->getLogWeight());
                            }
                            
                            double ess = 1.0/ess_inverse;
                            if (_verbose > 1) {
                                cout << "\t" << "ESS is : " << ess << endl;
                            }
                            
                            bool filter = true;
//                            if (ess < 2.0) {
//                                filter = false; // TODO: trying this
//                            }
                            if (Particle::_run_on_empty) {
                                filter = false;
                            }
                            if (filter) {
                                resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1);
                                //if use_first is true, my_vec = my_vec_2
                                //if use_first is false, my_vec = my_vec_1
                                
                                my_vec = use_first ? my_vec_2:my_vec_1;

                                //change use_first from true to false or false to true
                                use_first = !use_first;
                                
                                unsigned species_count = 0;
                                
                                for (auto &p:my_vec) {
                                    if (p->speciesJoinProposed()) {
                                        species_count++;
                                    }
                                }
                                
                                if (_verbose > 1) {
                                    cout << "\t" << "number of species join particles proposed = " << num_species_particles_proposed << endl;
                                    cout << "\t" << "number of species join particles accepted = " << species_count << endl;
                                }
                            resetWeights(my_vec);
                            }
                    } // g loop
                        
    //            saveAllHybridNodes(my_vec);
                
                saveSpeciesTrees(my_vec);
                for (int i=1; i < nsubsets+1; i++) {
                    saveGeneTree(i, my_vec);
                }
                
                if (Particle::_run_on_empty) { // make sure all gene tree log likelihoods are 0.0
                    vector<double> gene_tree_log_likelihoods;
                    gene_tree_log_likelihoods.resize(nsubsets);
                    for (int i=0; i<nsubsets; i++) {
                        gene_tree_log_likelihoods[i] = 0.0;
                    }
                    for (auto &p:my_vec) {
                        p->setLogLikelihood(gene_tree_log_likelihoods);
                    }
                }
                
                writeLoradFile(nsubsets, nspecies, ntaxa, my_vec);
                saveGeneTrees(nsubsets, my_vec);
                writeParamsFileForBeastComparison(nsubsets, nspecies, ntaxa, my_vec);
                
#if defined (EXTRA_SPECIES_SAMPLING)
                
                // TODO: do another round of species tree sampling
                for (auto &p:my_vec) {
                    // reset forest species partitions
                    p->resetSpecies();
                    p->mapSpecies(_taxon_map, _species_names);
//                    p->showParticle();
                }
                use_first = true;
                
                for (unsigned s=0; s<nspecies; s++){
                    cout << "beginning species tree proposals" << endl;
                    //taxon joining and reweighting step
                    
                    proposeSpeciesParticles(my_vec);
                    
                    normalizeSpeciesWeights(my_vec);
                    
                    double ess_inverse = 0.0;
                    
                    for (auto & p:my_vec) {
                        ess_inverse += exp(2.0*p->getSpeciesLogWeight());
                    }

                    double ess = 1.0/ess_inverse;
                    cout << "   " << "ESS = " << ess << endl;
                 
                    resampleSpeciesParticles(my_vec, use_first ? my_vec_2:my_vec_1);
                    //if use_first is true, my_vec = my_vec_2
                    //if use_first is false, my_vec = my_vec_1
                    
                    my_vec = use_first ? my_vec_2:my_vec_1;

                    //change use_first from true to false or false to true
                    use_first = !use_first;

                    resetSpeciesWeights(my_vec);
                } // s loop
                saveSpeciesTrees(my_vec);
//                for (auto &p:my_vec) {
//                    p->showParticle();
//                }
#endif
                
                if (_verbose > 0) {
                    cout << "marginal likelihood: " << _log_marginal_likelihood << endl;
                }
            }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }
        }

        std::cout << "\nFinished!" << std::endl;
    }
}

