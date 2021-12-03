#pragma once    ///start

#include <iostream>
#include "data.hpp"
//POL #include "likelihood.hpp"
#include "partition.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "xproj.hpp"
//#include "forest.hpp"
#include "particle.hpp"
using namespace std;
using namespace boost;  //POL

namespace proj {

    class Proj {
        public:
                                Proj();
                                ~Proj();

            void                clear();
            void                processCommandLineOptions(int argc, const char * argv[]);
            void                run();
            void                saveAllForests(const vector<Particle> &v) const ;
    
            void                normalizeWeights(vector<Particle> & particles);  //POL
            unsigned            chooseRandomParticle(vector<Particle> & particles);  //POL
            void                resampleParticles(vector<Particle> & particles);  //POL
            void                resetWeights(vector<Particle> & particles); //POL
            void                debugNormalizedWeights(const vector<Particle> & particles) const; //POL

        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            //Likelihood::SharedPtr       _likelihood;
        
            bool                        _use_gpu;
            bool                        _ambig_missing;


            static std::string     _program_name;
            static unsigned        _major_version;
            static unsigned        _minor_version;
            void                    summarizeData(Data::SharedPtr);
            void                    printFirstParticle(vector<Particle>);
            unsigned                setNumberSpecies(Data::SharedPtr);
            double                  getRunningSum(vector<double>);

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
        //POL _likelihood = nullptr;
    }   ///end_clear
///
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

    inline void Proj::summarizeData(Data::SharedPtr) {
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

    inline double Proj::getRunningSum(vector<double> log_weight_vec) {
        double running_sum = 0.0;
        double log_particle_sum = 0.0;
        
        double log_max_weight = *max_element(log_weight_vec.begin(), log_weight_vec.end());
        for (auto & i:log_weight_vec) {
            running_sum += exp(i - log_max_weight);
        }
        log_particle_sum = log(running_sum) + log_max_weight;
        
        return log_particle_sum;
    }

    //POL
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
    
    //POL
    inline void Proj::normalizeWeights(vector<Particle> & particles) { //POL
        unsigned i = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            log_weight_vec[i++] = p.getLogWeight();
        }

        double log_particle_sum = getRunningSum(log_weight_vec);

        for (auto & p : particles) {
            p.setLogWeight(p.getLogWeight() - log_particle_sum);
        }
    }
    
    //POL
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
    
    //POL
    inline void Proj::resampleParticles(vector<Particle> & particles) { //POL
        // throw darts
        unsigned nparticles = (unsigned)particles.size();
        vector<double> darts(nparticles, 0.0);
        for(unsigned i = 0; i < nparticles; i++) {
            unsigned j = chooseRandomParticle(particles);
            darts[j]++;
        }

#if 0
        // show darts
        double cumw = 0.0;
        unsigned cumd = 0;
        cout << format("\n%12s %12s %12s\n") % "particle" % "weight" % "darts";
        for (unsigned i = 0; i < nparticles; i++) {
            double w = exp(particles[i].getLogWeight());
            cumw += w;
            cumd += darts[i];
            cout << format("%12d %12.5f %12d\n") % i % w % darts[i];
        }
        cout << format("%12s %12.5f %12d\n") % " " % cumw % cumd;
#endif
        
        // create new particle vector
        vector<Particle> new_particles;
        for (unsigned i = 0; i < nparticles; i++) {
            for (unsigned k = 0; k < darts[i]; k++) {
                new_particles.push_back(particles[i]);
            }
        }
        assert(nparticles == new_particles.size());
        
        // copy particles
        copy(new_particles.begin(), new_particles.end(), particles.begin());
    }
    
    //POL
    inline void Proj::resetWeights(vector<Particle> & particles) { //POL
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p.setLogWeight(logw);
        }
    }
    
    inline void Proj::run() {
        std::cout << "Starting..." << std::endl;
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;
        
        try {
            std::cout << "\n*** Reading and storing the data in the file " << _data_file_name << std::endl;
            std::cout << "data file name is " << _data_file_name << std::endl;
            _data = Data::SharedPtr(new Data());
            _data->setPartition(_partition);
            _data->getDataFromFile(_data_file_name);
            
            summarizeData(_data);
            
            //set number of species to number in data file
            unsigned nspecies = setNumberSpecies(_data);
            rng.setSeed(5);
            
//          create vector of particles
            unsigned nparticles = 500;
            vector<Particle> my_vec(nparticles);
            for (auto & p:my_vec ) {
                p.setData(_data);
                p.setInitialWeight();
            }

//            printFirstParticle(my_vec);

            normalizeWeights(my_vec); //POL
            resampleParticles(my_vec); //POL
            resetWeights(my_vec); //POL

            //run through each generation of particles
            for (unsigned g=0; g<nspecies-2; g++){
                //POL cout << "Generation " << g << endl;
                
                vector<double> log_weight_vec;
                double log_weight = 0.0;
                
//                cout << "\n Particles after generation " << g << endl;
                
                //taxon joining and reweighting step
                for (auto & p:my_vec) {
                    log_weight = p.proposal();
                    log_weight_vec.push_back(log_weight);
                    //POL p.showParticle();
                }
                
#if 1
                normalizeWeights(my_vec); //POL
#else
                double log_particle_sum = getRunningSum(log_weight_vec);
                
                double ess_inverse = 0.0;
                
                for (auto & p:my_vec) {
                    //normalized weight is weight / particle sum = ln(weight)-ln(particle_sum)
                    p.setLogWeight(p.getLogWeight() - log_particle_sum);
                    double log_weight = p.getLogWeight();
//                    ESS = 1/(weight^2)
                    ess_inverse += exp(2.0*log_weight);
                    
//                    assert(ess_inverse > 0);
//                    p.showParticle();
                }
                
                double ess = 1.0/ess_inverse;
                cout << "ESS is " << ess << "  " << "g = " << g << endl;
                
                //debugNormalizedWeights(my_vec); //POL
#endif

#if 1
                resampleParticles(my_vec); //POL
                resetWeights(my_vec); //POL
#else
                if (ess < 100) {
                    vector<double> cum_probs(my_vec.size(), 0.0);
                    unsigned ndarts = (unsigned) my_vec.size();
                    //sample particles
                    for(unsigned i=0; i<ndarts; i++) {
                        double u = rng.uniform();
                        double cum_prob = 0.0;
                        unsigned j = 0;
                        for(auto & p:my_vec) {
                            cum_prob += exp(p.getLogWeight());
                            if (u < cum_prob) {
                                cum_probs[j]+=1.0;
                                break;
                            }
                            j++;
                        }
                    }
                    cum_probs[0]=cum_probs[0]/ndarts;
                    if (cum_probs[0]>0.0) {
                        my_vec[0].setLogWeight(log(cum_probs[0]));
                    }
                    for(unsigned i=1; i<my_vec.size(); i++) {
                        double w = cum_probs[i]/ndarts;
                        if (w>0.0){
                            my_vec[i].setLogWeight(log(w));
                        }
                        cum_probs[i]=cum_probs[i-1]+w;
//                        if (w>0.0) {
//                            cout << w << " | " << my_vec[i].firstPair() << endl;
//                            my_vec[i].showParticle();
//                        }
                    }
                    
                    //filter particles
//                    double logNumParticles = log(my_vec.size());
                    vector<Particle> my_vec2;
                    for(unsigned i=0; i<my_vec.size(); i++) {
                        double u = rng.uniform();
                        for(unsigned j=0; j<my_vec.size(); j++) {
                            if (u < cum_probs[j]) {
                                my_vec2.push_back(my_vec[j]);
                                break;
                            }
                        }
                    }
                    copy(my_vec2.begin(), my_vec2.end(), my_vec.begin());
                    
                    //reset weights
                    double log_w = -log(my_vec.size());
                    for(unsigned i=0; i<my_vec.size(); i++) {
                        //POL my_vec2[i].debugParticle(boost::str(boost::format("old %d") % i));
                        //POL my_vec[i].debugParticle(boost::str(boost::format("new %d") % i));
                        //POL cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
                        my_vec[i].setLogWeight(log_w);
                    }
                    cout << endl;
                } // if (ess < 100)...
#endif
            } // g loop
            double sum_h = 0.0;
            for (auto & p:my_vec) {
                double h = p.calcHeight();
                sum_h += h;
            }
            sum_h/=my_vec.size();
            cout << "mean height equals " << sum_h << endl;
            
            saveAllForests(my_vec);
            }
        
        
        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }

} ///end
