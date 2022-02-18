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
            void                saveAllForests(const vector<Particle> &v) const ;

            void                normalizeWeights(vector<Particle> & particles);
            unsigned            chooseRandomParticle(vector<Particle> & particles);
            void                resampleParticles(vector<Particle> & from_particles, vector<Particle> & to_particles);
            void                resetWeights(vector<Particle> & particles);
            void                debugNormalizedWeights(const vector<Particle> & particles) const;
            void                calcMarginalLikelihood(double weight_average);
            double              getWeightAverage(vector<double> log_weight_vec);

        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _log_marginal_likelihood = 0.0;
            bool                        _use_gpu;
            bool                        _ambig_missing;
            unsigned                    _nparticles;
            unsigned                    _random_seed;
            double                      _theta;


            static std::string           _program_name;
            static unsigned              _major_version;
            static unsigned              _minor_version;
            void                          summarizeData(Data::SharedPtr);
            void                          printFirstParticle(vector<Particle>);
            unsigned                      setNumberSpecies(Data::SharedPtr);
            double                        getRunningSum(const vector<double> &) const;

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

    inline void Proj::saveAllForests(const vector<Particle> &v) const {
        ofstream treef("forest.trees");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:v) {
            treef << "  tree test = [&R] " << p.saveForestNewick()  << ";\n";
        }
        treef << "end;\n";
        treef.close();
    }

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {   ///begin_processCommandLineOptions
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
        ("theta, t", boost::program_options::value(&_theta)->default_value(0.05), "theta")
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

    inline void Proj::printFirstParticle(vector<Particle> my_vec) {
//            print out particles at the start
            cout << "\n Particles at the start: " << endl;
             for (auto & p:my_vec) {
                p.showParticle();
            }
    }

    inline unsigned Proj::setNumberSpecies(Data::SharedPtr) {
        unsigned nspecies;
        nspecies = _data->getNumTaxa();
        Forest::setNumSpecies(nspecies);
        return nspecies;
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

    inline void Proj::debugNormalizedWeights(const vector<Particle> & particles) const {
        cout << format("%12s %12s %12s\n") % "particle" % "weight" % "cumweight";
        double cumw = 0.0;
        unsigned i = 0;
        for (auto p : particles) {
            double w = exp(p.getLogWeight());
            cumw += w;
            cout << format("%12d %12.5f %12.5f\n") % (++i) % w % cumw;
        }
        cout << endl;
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
    }

    inline unsigned Proj::chooseRandomParticle(vector<Particle> & particles) {
        int chosen_index = -1;
        unsigned nparticles = (unsigned)particles.size();
        double u = rng.uniform();
        double cum_prob = 0.0;
        for(unsigned j = 0; j < nparticles; j++) {
            cum_prob += exp(particles[j].getLogWeight());
            if (u < cum_prob) {
                chosen_index = j;
                break;
            }
        }
        assert(chosen_index > -1);
        return chosen_index;
    }

    inline void Proj::resampleParticles(vector<Particle> & from_particles, vector<Particle> & to_particles) {
        // throw darts
        unsigned nparticles = (unsigned)from_particles.size();
        vector<double> darts(nparticles, 0.0);
        for(unsigned i = 0; i < nparticles; i++) {
            unsigned j = chooseRandomParticle(from_particles);
            darts[j]++;
        }

        // show darts
        double cumw = 0.0;
        unsigned cumd = 0;
//        cout << format("\n%12s %12s %12s\n") % "particle" % "weight" % "darts";
        for (unsigned i = 0; i < nparticles; i++) {
            double w = exp(from_particles[i].getLogWeight());
            cumw += w;
            cumd += darts[i];
//            cout << format("%12d %12.5f %12d\n") % i % w % darts[i];
        }
//        cout << format("%12s %12.5f %12d\n") % " " % cumw % cumd;

        // create new particle vector
        unsigned m = 0;
        for (unsigned i = 0; i < nparticles; i++) {
            for (unsigned k = 0; k < darts[i]; k++) {
                to_particles[m++]=from_particles[i];
            }
        }
        assert(nparticles == to_particles.size());

        // copy particles
//        copy(to_particles.begin(), to_particles.end(), from_particles.begin());
    }

    inline void Proj::resetWeights(vector<Particle> & particles) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p.setLogWeight(logw);
        }
    }

    inline void Proj::run() {
        std::cout << "Starting..." << std::endl;
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;
        std::cout << "Random seed: " << _random_seed << std::endl;
        std::cout << "Theta: " << _theta << std::endl;

        try {
            std::cout << "\n*** Reading and storing the data in the file " << _data_file_name << std::endl;
            std::cout << "data file name is " << _data_file_name << std::endl;
            _data = Data::SharedPtr(new Data());
            _data->setPartition(_partition);
            _data->getDataFromFile(_data_file_name);

            summarizeData(_data);
            ps.setnelements(4*_data->getNumPatterns());
            
            //set number of species to number in data file
            unsigned nspecies = setNumberSpecies(_data);
            rng.setSeed(_random_seed);
            
            Forest::setTheta(_theta);

//          create vector of particles
            unsigned nparticles = _nparticles;
            vector<Particle> my_vec_1(nparticles);
            vector<Particle> my_vec_2(nparticles);
            vector<Particle> &my_vec = my_vec_1;
            bool use_first = true;
            for (auto & p:my_vec ) {
                p.setData(_data);
            }

            _log_marginal_likelihood = 0.0;
            //run through each generation of particles
            for (unsigned g=0; g<nspecies-1; g++){
//                cout << "Generation " << g << endl;

                vector<double> log_weight_vec;
                double log_weight = 0.0;

                //taxon joining and reweighting step
                for (auto & p:my_vec) {
                    log_weight = p.proposal();
                    log_weight_vec.push_back(log_weight);
                }
                
                normalizeWeights(my_vec);
                
                resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1);
                
                //if use_first is true, my_vec = my_vec_2
                //if use_first if alse, my_vec = my_vec_1
                my_vec = use_first ? my_vec_2:my_vec_1;
                
                //change use_first from true to false or false to true
                use_first = !use_first;
                                
                resetWeights(my_vec);
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
            }


        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }

}
