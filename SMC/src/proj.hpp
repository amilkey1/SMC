#pragma once

#include <iostream>
#include "data.hpp"
#include "partition.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "xproj.hpp"
#include "particle.hpp"

using namespace std;
using namespace boost;

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
        void                saveAllForests(vector<Particle> &v) const ;

            void                normalizeWeights(vector<Particle> & particles);
            unsigned            chooseRandomParticle(vector<Particle> & particles, vector<double> & cum_prob);
            void                resampleParticles(vector<Particle> & from_particles, vector<Particle> & to_particles);
            void                resetWeights(vector<Particle> & particles);
            double              getWeightAverage(vector<double> log_weight_vec);
            void                createSpeciesMap(Data::SharedPtr);
            void                showParticlesByWeight(vector<Particle> my_vec);

        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood = 0.0;
            bool                        _use_gpu;
            bool                        _ambig_missing;
            unsigned                    _nparticles;
            unsigned                    _random_seed;


            static std::string           _program_name;
            static unsigned              _major_version;
            static unsigned              _minor_version;
            void                          summarizeData(Data::SharedPtr);
            unsigned                      setNumberTaxa(Data::SharedPtr);
            double                        getRunningSum(const vector<double> &) const;
            vector<string>                _species_names;
            map<string, string>           _taxon_map;

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
    }

//    inline void Proj::saveAllForests(const vector<Particle> &v) const {
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
        ("speciation_rate", boost::program_options::value(&Forest::_speciation_rate)->default_value(10), "speciation rate")
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
//        Forest::setNumSpecies(ntaxa);
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

    inline void Proj::normalizeWeights(vector<Particle> & particles) {
        unsigned i = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            log_weight_vec[i++] = p.getLogWeight();
        }

        double log_particle_sum = getRunningSum(log_weight_vec);

        for (auto & p : particles) {
            p.setLogWeight(p.getLogWeight() - log_particle_sum);
        }

        _log_marginal_likelihood += log_particle_sum - log(_nparticles);
        sort(particles.begin(), particles.end(), greater<Particle>());
    }

    inline unsigned Proj::chooseRandomParticle(vector<Particle> & particles, vector<double> & cum_probs) {
        int chosen_index = -1;
        unsigned nparticles = (unsigned)particles.size();
        double u = rng.uniform();
        double cum_prob = 0.0;
        
        for(unsigned j = 0; j < nparticles; j++) {
            if (cum_probs[j]<0.0){
                cum_prob += exp(particles[j].getLogWeight());
                cum_probs[j] = cum_prob;
            }
            else
                cum_prob = cum_probs[j];
            
            if (u < cum_prob) {
                chosen_index = j;
                break;
            }
        }
        assert(chosen_index > -1);
        return chosen_index;
    }

    inline void Proj::resampleParticles(vector<Particle> & from_particles, vector<Particle> & to_particles) {
        unsigned nparticles = (unsigned)from_particles.size();
        vector<double> cum_probs(nparticles, -1.0);
        
        // throw darts
        vector<unsigned> darts(nparticles, 0);
        unsigned max_j = 0;
        for(unsigned i = 0; i < nparticles; i++) {
            unsigned j = chooseRandomParticle(from_particles, cum_probs);
            if (j>max_j){
                max_j = j;
            }
            darts[j]++;
        }

        // create new particle vector
        unsigned m = 0;
        for (unsigned i = 0; i < nparticles; i++) {
            for (unsigned k = 0; k < darts[i]; k++) {
                to_particles[m++]=from_particles[i];
            }
        }
        assert(nparticles == to_particles.size());
    }

    inline void Proj::resetWeights(vector<Particle> & particles) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p.setLogWeight(logw);
        }
    }
    
    inline void Proj::createSpeciesMap(Data::SharedPtr d) {
        const vector<string> &names = d->getTaxonNames();
        for (auto &name:names) {
            regex re(".+\\^(.+)");
            smatch match_obj;
            bool matched=regex_match(name, match_obj, re); //search name for regular expression, store result in match_obj
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

    inline void Proj::showParticlesByWeight(vector<Particle> my_vec) {
        vector <double> weights;
        
        //create weight vector
        for (auto & p:my_vec) {
            weights.push_back(p.getLogWeight());
        }
        
        //sort particles by weight
        sort(my_vec.begin(), my_vec.end(), greater<Particle>());
        
        //print first particle
        cout << "\n" << "Heaviest particle: ";
        my_vec[0].showParticle();
    }

    inline void Proj::run() {
        std::cout << "Starting..." << std::endl;
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;
        std::cout << "Random seed: " << _random_seed << std::endl;
        std::cout << "Theta: " << Forest::_theta << std::endl;

        try {
            std::cout << "\n*** Reading and storing the data in the file " << _data_file_name << std::endl;
            std::cout << "data file name is " << _data_file_name << std::endl;
            _data = Data::SharedPtr(new Data());
            _data->setPartition(_partition);
            _data->getDataFromFile(_data_file_name);

            summarizeData(_data);
            createSpeciesMap(_data);

            //set number of species to number in data file
            unsigned ntaxa = setNumberTaxa(_data);
            unsigned nspecies = _species_names.size();
            Forest::setNumSpecies(nspecies);
            rng.setSeed(_random_seed);

//          create vector of particles
            unsigned nparticles = _nparticles;
            
            unsigned nsubsets = _data->getNumSubsets();
            Particle::setNumSubsets(nsubsets);
            
            vector<Particle> my_vec_1(nparticles);
            vector<Particle> my_vec_2(nparticles);
            vector<Particle> &my_vec = my_vec_1;
            bool use_first = true;
            for (auto & p:my_vec ) {
                p.setData(_data, _taxon_map);
                p.mapSpecies(_taxon_map, _species_names);
            }

            _log_marginal_likelihood = 0.0;
            //run through each generation of particles
            for (unsigned g=0; g<nspecies; g++){

                vector<double> log_weight_vec;
                double log_weight = 0.0;

                //taxon joining and reweighting step 
                for (auto & p:my_vec) {
                    log_weight = p.proposal();
                    log_weight_vec.push_back(log_weight);
//                    p.showParticle();
                }
                
                normalizeWeights(my_vec);
                resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1);

                //if use_first is true, my_vec = my_vec_2
                //if use_first is false, my_vec = my_vec_1
                
                my_vec = use_first ? my_vec_2:my_vec_1;

                //change use_first from true to false or false to true
                use_first = !use_first;

                // TODO: is this necessary?
                resetWeights(my_vec);
                
                for (auto &p:my_vec){
                    p.showParticle();
                }
            } // g loop

            double sum_h = 0.0;
            for (auto & p:my_vec) {
                double h = p.calcHeight();
                sum_h += h;
            }
            sum_h/=my_vec.size();
            cout << "mean height equals " << sum_h << endl;
            cout << "log marginal likelihood = " << _log_marginal_likelihood << endl;
            cout << "theta = " << Forest::_theta << endl;

            saveAllForests(my_vec);
//            showParticlesByWeight(my_vec);
            }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }
}
