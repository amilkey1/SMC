#pragma once    ///start

#include <iostream>
#include "data.hpp"
#include "likelihood.hpp"
#include "partition.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "xproj.hpp"
//#include "forest.hpp"
#include "particle.hpp"
using namespace std;

namespace proj {

    class Proj {
        public:
                                Proj();
                                ~Proj();

            void                clear();
            void                processCommandLineOptions(int argc, const char * argv[]);
            void                run();

        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            Likelihood::SharedPtr       _likelihood;
        
            bool                        _use_gpu;
            bool                        _ambig_missing;


            static std::string     _program_name;
            static unsigned        _major_version;
            static unsigned        _minor_version;

    };

    //end_class_declaration
    inline Proj::Proj() { ///begin_constructor
//        std::cout << "Constructing a Proj" << std::endl;
        clear();
    } //end_constructor

    inline Proj::~Proj() { ///begin_destructor
//        std::cout << "Destroying a Proj" << std::endl;
    } //end_destructor

    inline void Proj::clear() {    ///begin_clear
        _data_file_name = "";
        _partition.reset(new Partition());
        _use_gpu        = true;
        _ambig_missing  = true;
        _data = nullptr;
        _likelihood = nullptr;
    }   ///end_clear

    inline void Proj::processCommandLineOptions(int argc, const char * argv[]) {   ///begin_processCommandLineOptions
        std::vector<std::string> partition_subsets;
        boost::program_options::variables_map vm;
        boost::program_options::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("datafile,d",  boost::program_options::value(&_data_file_name)->required(), "name of a data file in NEXUS format")
        ("subset",  boost::program_options::value(&partition_subsets), "a string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
        ("gpu",           boost::program_options::value(&_use_gpu)->default_value(true),                "use GPU if available")
        ("ambigmissing",  boost::program_options::value(&_ambig_missing)->default_value(true),          "treat all ambiguities as missing data")
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
    }   ///end_processCommandLineOptions

    inline void Proj::run() {
        std::cout << "Starting..." << std::endl;
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;
        
        try {
            std::cout << "\n*** Reading and storing the data in the file " << _data_file_name << std::endl;
            std::cout << "data file name is " << _data_file_name << std::endl;
            _data = Data::SharedPtr(new Data());
            _data->setPartition(_partition);
            _data->getDataFromFile(_data_file_name);
            
            // Report information about data partition subsets
            unsigned nsubsets = _data->getNumSubsets();
            std::cout << "\nNumber of taxa: " << _data->getNumTaxa() << std::endl;
            
            std::cout << "Number of partition subsets: " << nsubsets << std::endl;
            
            for (unsigned subset = 0; subset < nsubsets; subset++) {
                DataType dt = _partition->getDataTypeForSubset(subset);
                std::cout << "  Subset " << (subset+1) << " (" << _data->getSubsetName(subset) << ")" << std::endl;
                std::cout << "    data type: " << dt.getDataTypeAsString() << std::endl;
                std::cout << "    sites:     " << _data->calcSeqLenInSubset(subset) << std::endl;
                std::cout << "    patterns:  " << _data->getNumPatternsInSubset(subset) << std::endl;
                }
            
            //set number of species to number in data file
            rng.setSeed(12345);
            unsigned nspecies;
            nspecies = _data->getNumTaxa();
            Forest::setNumSpecies(nspecies);
            
            //create vector of particles
            unsigned nparticles = 100;
            vector<Particle> my_vec(nparticles);
            for (auto & p:my_vec ) {
                p.setData(_data);
            }

            //print out particles at the start
//            cout << "\n Particles at the start: " << endl;
//            for (auto & p:my_vec ) {
//                p.showParticle();
//            }
            
            //run through each generation of particles
            for (unsigned g=0; g<nspecies-2; g++){
                vector<double> log_weight_vec;
                double running_sum = 0.0;
                double log_weight = 0.0;
                double log_particle_sum = 0.0;
                
                cout << "\n Particles after generation " << g << endl;
                
                //taxon joining and reweighting step
                for (auto & p:my_vec) {
                    log_weight = p.proposal();
                    log_weight_vec.push_back(log_weight);
                }
                
                double log_max_weight = *max_element(log_weight_vec.begin(), log_weight_vec.end());
                for (auto & i:log_weight_vec) {
                    running_sum += exp(i - log_max_weight);
                }
                log_particle_sum = log(running_sum) + log_max_weight;
                
                double ess_inverse = 0.0;
                
                for (auto & p:my_vec) {
                    //normalized weight is weight / particle sum = ln(weight)-ln(particle_sum)
                    p.setLogWeight(p.getLogWeight() - log_particle_sum);
//                    log_weight = p.getLogWeight() - log_particle_sum;
//                    ESS = 1/(weight^2)
                    ess_inverse += exp(2.0*p.getLogWeight());
                }
                double ess = 1.0/ess_inverse;
                cout << "ESS is " << ess << endl;
                
#if 0
                vector<unsigned> counts(my_vec.size(), 0);
                if (ess < 2.0) {
                    //sample particles
                    for(unsigned i=0; i<my_vec.size(); i++) {
                        double u = rng.uniform();
                        double cum_prob = 0.0;
                        unsigned j = 0;
                        for(auto & p:my_vec) {
                            cum_prob += exp(p.getLogWeight());
                            if (u < cum_prob) {
                                counts[j]++;
                                break;
                            }
                            j++;
                        }
                    }
                    //filter particles
                    double logNumParticles = log(my_vec.size());
                    vector<Particle> my_vec2;
                    for(unsigned i=0; i<my_vec.size(); i++) {
                        if (counts[i]>0) {
                            my_vec[i].setLogWeight(log(counts[i])-logNumParticles);
                            for(unsigned j=0; j<counts[i]; j++) {
                                my_vec2.push_back(my_vec[i]);
                            }
                        }
                    }
                    copy(my_vec2.begin(), my_vec2.end(), my_vec.begin());
                }
                
#endif
//                for (auto & p:my_vec) {
//                    p.showParticle(); //print all new particles
//                    double log_likelihood = p.calcLogLikelihood();
//                    p.savePaupFile("paup.nex", _data_file_name, "forest.tre", log_likelihood);
//                    p.saveForest("forest.tre");
//                }
            }
            double sum_h = 0.0;
            for (auto & p:my_vec) {
//               construct a vector of the new particles...?
                p.showParticle(); //print all new particles
                double log_likelihood = p.calcLogLikelihood();
                p.savePaupFile("paup.nex", _data_file_name, "forest.tre", log_likelihood);
                p.saveForest("forest.tre");
                double h = p.calcHeight();
                sum_h += h;
                            }
            sum_h/=my_vec.size();
            cout << "mean height equals " << sum_h << endl;
            }
        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }

} ///end
