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
            void                normalizeWeights(vector<Particle::SharedPtr> & particles);
            void                resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
            void                resetWeights(vector<Particle::SharedPtr> & particles);
            double              getWeightAverage(vector<double> log_weight_vec);
            void                createSpeciesMap(Data::SharedPtr);
//            void                showParticlesByWeight(vector<Particle::SharedPtr> my_vec);
            void                showFinal(vector<Particle::SharedPtr>);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles);
            void                proposeParticles(vector<Particle::SharedPtr> &particles);
            void                saveAllHybridNodes(vector<Particle::SharedPtr> &v) const;

        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood = 0.0;
            double                      _prev_log_marginal_likelihood = 0.0;
            bool                        _use_gpu;
            bool                        _ambig_missing;
            unsigned                    _nparticles;
            unsigned                    _random_seed;


            static std::string          _program_name;
            static unsigned             _major_version;
            static unsigned             _minor_version;
            void                        summarizeData(Data::SharedPtr);
            unsigned                    setNumberTaxa(Data::SharedPtr);
            double                      getRunningSum(const vector<double> &) const;
            vector<string>              _species_names;
            map<string, string>         _taxon_map;
            unsigned                    _nthreads;
            void                        handleBaseFrequencies();
            void                        debugSpeciesTree(vector<Particle::SharedPtr> &particles);
            void                        estimateTheta(vector<Particle::SharedPtr> &particles);
            void                        estimateSpeciationRate(vector<Particle::SharedPtr> &particles);
            double                      _small_enough;
            unsigned                    _verbose;
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

    inline void Proj::writeLoradFile(unsigned ngenes, unsigned nspecies, unsigned ntaxa, vector<Particle::SharedPtr> &v) const {
        ofstream logf("params.log");
        logf << "iteration ";
        logf << "\t" << "likelihood ";
        logf << "\t" << "species_tree_topology_prior ";
        for (int i=1; i<ngenes+1; i++) {
            logf << "\t" << "gene_tree_topology_prior ";
        }
        for (int s=0; s<nspecies-1; s++) {
            logf << "\t" << "species_increment";
            logf << "\t" << "species_increment_prior";
        }
        for (int g=1; g<ngenes+1; g++) {
            for (int i=1; i<ntaxa; i++) {
                logf << "\t" << "gene_increment";
                logf << "\t" << "log_gene_increment_prior";
            }
        }
        logf << "\t" << "coalescent_likelihood";
        logf << endl;
        
        unsigned iter = 0;
        for (auto &p:v) {
            logf << iter;
            iter++;
            
            logf << "\t" << p->getLogLikelihood();
            
            for (unsigned g=0; g<ngenes+1; g++) {
                logf << "\t" << p->getTopologyPrior(g);
            }

            for (unsigned g=0; g<ngenes+1; g++) {
                for (auto &b:p->getIncrementPriors(g)) {
                    logf << "\t" << b.first;
                    logf << "\t" << b.second;
                    // no increment or increment prior should be 0
                    assert (b.first > 0.0);
                    assert (b.second != 0.0);
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
            ofstream treef("species_trees.trees");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:v) {
                treef << "  tree test = [&R] " << p->saveForestNewick()  << ";\n";
            }
            treef << "end;\n";
            treef.close();
        }

    inline void Proj::saveGeneTrees(unsigned ngenes, vector<Particle::SharedPtr> &v) const {
        ofstream treef("gene_trees.txt");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:v) {
                for (int i=1; i<ngenes+1; i++) {
                    treef << "  gene " << i << " = [&R] " << p->saveGeneNewick(i)  << ";\n";
            }
            treef << endl;
        }
        treef << "end;\n";
        treef.close();
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
        ("datafile,d",  boost::program_options::value(&_data_file_name)->required(), "name of a data file in NEXUS format")
        ("subset",  boost::program_options::value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("gpu",           boost::program_options::value(&_use_gpu)->default_value(true), "use GPU if available")
        ("ambigmissing",  boost::program_options::value(&_ambig_missing)->default_value(true), "treat all ambiguities as missing data")
        ("nparticles",  boost::program_options::value(&_nparticles)->default_value(1000), "number of particles")
        ("seed,z", boost::program_options::value(&_random_seed)->default_value(1), "random seed")
        ("theta, t", boost::program_options::value(&Forest::_theta)->default_value(0.05), "theta")
        ("lambda", boost::program_options::value(&Forest::_lambda)->default_value(1), "speciation rate")
        ("proposal",  boost::program_options::value(&Forest::_proposal)->default_value("prior-post"), "a string defining a proposal (prior-prior or prior-post)")
        ("model", boost::program_options::value(&Forest::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&Forest::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("nthreads",  boost::program_options::value(&_nthreads)->default_value(1.0), "number of threads for multi threading")
        ("migration_rate", boost::program_options::value(&Forest::_migration_rate)->default_value(0.0), "migration rate")
        ("hybridization_rate", boost::program_options::value(&Forest::_hybridization_rate)->default_value(0.0), "hybridization rate")
        ("run_on_empty", boost::program_options::value(&Particle::_run_on_empty)->default_value(false), "run program without data")
        ("verbose", boost::program_options::value(&_verbose)->default_value(0), "set amount of output printed")
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

    inline void Proj::normalizeWeights(vector<Particle::SharedPtr> & particles) {
        unsigned i = 0;
        unsigned species_joins = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            if (p->speciesJoinProposed()) {
                species_joins++;
            }
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
//        sort(particles.begin(), particles.end(), greater<Particle>());
    }

    inline void Proj::resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles) {
         
    //        unsigned nparticles = (unsigned)from_particles.size();
         unsigned nparticles = _nparticles;
         assert (from_particles.size() == nparticles);
         assert (to_particles.size() == nparticles);
         
         vector<pair<double, double>> cum_probs;
             // Create vector of pairs p, with p.first = log weight and p.second = particle index
         cum_probs.resize(nparticles);
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
//             cout << "selected index is: " << sel_index;
             
             Particle::SharedPtr p0 = from_particles[sel_index];
             to_particles[i]=Particle::SharedPtr(new Particle(*p0));
             
             assert(nparticles == to_particles.size());
         }
     }
    inline void Proj::resetWeights(vector<Particle::SharedPtr> & particles) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p->setLogWeight(logw);
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

//        saveAllForests(my_vec);
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

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles) {
        for (unsigned i=first; i<last; i++){
            particles[i]->proposal();
        }
    }

//    inline void Proj::showParticlesByWeight(vector<Particle::SharedPtr> my_vec) {
//        vector <double> weights;
//
//        //create weight vector
//        for (auto & p:my_vec) {
//            weights.push_back(p->getLogWeight());
//        }
//
//        //sort particles by weight
//        sort(my_vec.begin(), my_vec.end(), greater<Particle>());
//
//        //print first particle
//        cout << "\n" << "Heaviest particle: ";
//        my_vec[0]->showParticle();
//    }

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
                cout << "data file name is " << _data_file_name << std::endl;
            }
            _data = Data::SharedPtr(new Data());
            _data->setPartition(_partition);
            _data->getDataFromFile(_data_file_name);

            if (_verbose > 0) {
                summarizeData(_data);
            }
            createSpeciesMap(_data);

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
                for (auto & p:my_vec ) {
                    p->setData(_data, _taxon_map);
                    p->mapSpecies(_taxon_map, _species_names);
                }
                
                _prev_log_marginal_likelihood = _log_marginal_likelihood;
                
                // reset marginal likelihood
                _log_marginal_likelihood = 0.0;
                for (auto &p:my_vec) {
                    p->calcLogLikelihood();
                }
                normalizeWeights(my_vec); // initialize marginal likelihood
                
                //run through each generation of particles
            
            unsigned nsteps = (ntaxa-1)*nsubsets-1+1;
                
//                    for (unsigned g=0; g<(ntaxa-1)*nsubsets+nspecies-1; g++){
                    for (unsigned g=0; g<nsteps; g++){

                    if (_verbose > 0) {
                        cout << "starting step " << g << " of " << nsteps-1 << endl;
                    }
                    //taxon joining and reweighting step
                    proposeParticles(my_vec);
                    
                    unsigned num_species_particles_proposed = 0;
                    
                    double ess_inverse = 0.0;
                    for (auto &p:my_vec) {
                        if (p->speciesJoinProposed()) {
                            num_species_particles_proposed++;
                        }
//                            p->showParticle();
                    }
                    normalizeWeights(my_vec);
                    
                    for (auto & p:my_vec) {
                        ess_inverse += exp(2.0*p->getLogWeight());
                    }
                    
                    double ess = 1.0/ess_inverse;
                    if (_verbose > 1) {
                        cout << "\t" << "ESS is : " << ess << endl;
                    }
                    
//                        if (ess < 100) {
//                        if (g % 2) { // resample every other generation
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
                    
                } // g loop
                    
            saveAllHybridNodes(my_vec);
            
            saveSpeciesTrees(my_vec);
            for (int i=1; i < nsubsets+1; i++) {
                saveGeneTree(i, my_vec);
            }
            writeLoradFile(nsubsets, nspecies, ntaxa, my_vec);
            
            if (_verbose > 0) {
                cout << "marginal likelihood: " << _log_marginal_likelihood << endl;
            }
        }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }
}

