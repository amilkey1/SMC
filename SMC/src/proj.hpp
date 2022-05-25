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
            void                proposeTheta();
            void                proposeSpeciationRate();
            double              logThetaPrior(double theta);
            double              logSpeciationRatePrior(double speciation_rate);
            string              acceptTheta();
            string              acceptSpeciationRate();
            void                showFinal(vector<Particle>);
            void                tune(bool accepted);

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
            double                      _prev_theta = 0.0;
            double                      _prev_speciation_rate = 0.0;
            vector<Particle>            _accepted_particle_vec;
            vector<Particle>            _prev_particles;
            vector<double>              _theta_vector;
            vector<double>              _speciation_rate_vector;
            double                      _theta_accepted_number = 0.0;
            double                      _speciation_rate_accepted_number = 0.0;
            string                      _sample;
            int                         _nattempts = 0;
            bool                        _tuning = true;
            double                      _target_acceptance = 0.3;
            double                      _lambda;

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
        ("speciation_rate", boost::program_options::value(&Forest::_speciation_rate)->default_value(1), "speciation rate")
        ("proposal",  boost::program_options::value(&Forest::_proposal)->default_value("prior-post"), "a string defining a proposal (prior-prior or prior-post)")
        ("model", boost::program_options::value(&Forest::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
//        ("base_frequencies", boost::program_options::value(&Forest::_base_frequencies), "base frequencies A C G T")
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

    inline void Proj::tune(bool accepted) {
        _nattempts++;
        if (_tuning) {
            double gamma_n = 10.0/(100.0 + (double)_nattempts);
            if (accepted)
                _lambda *= 1.0 + gamma_n*(1.0 - _target_acceptance)/(2.0*_target_acceptance);
            else
                _lambda *= 1.0 - gamma_n*0.5;

            // Prevent run-away increases in boldness for low-information marginal densities
            if (_lambda > 1000.0)
                _lambda = 1000.0;
        }
    }

    inline void Proj::proposeTheta() {
        double u = rng.uniform();
        
        // TODO: window width?
//        double window_width = 0.4;
        double proposed_theta = _lambda*u+(Forest::_theta-_lambda/2.0);
        
        // make sure proposed theta is positive
        if (proposed_theta < 0.0) {
            proposed_theta*=-1;
        }
        
        _prev_theta = Forest::_theta;
        Forest::_theta = proposed_theta;
    }

    inline string Proj::acceptTheta() {
        double u = rng.uniform();
        double log_acceptance_ratio = (_log_marginal_likelihood+logThetaPrior(Forest::_theta))-(_prev_log_marginal_likelihood+logThetaPrior(_prev_theta));
        if (log(u) > log_acceptance_ratio){
            // reject proposed theta
            bool accepted = false;
            tune(accepted);
            return "reject";
        }
        else {
            bool accepted = true;
            tune(accepted);
            return "accept";
        }
    }

    inline void Proj::proposeSpeciationRate() {
        double u = rng.uniform();

        // TODO: window width?
//        _lambda = 25.0;
        double proposed_speciation_rate = _lambda*u+(Forest::_speciation_rate-_lambda/2.0);
        
        // make sure proposed speciation rate is positive
        if (proposed_speciation_rate < 0.0) {
            proposed_speciation_rate*=-1;
        }
        
        _prev_speciation_rate = Forest::_speciation_rate;
        Forest::_speciation_rate = proposed_speciation_rate;
    }

    inline string Proj::acceptSpeciationRate() {
        double u = rng.uniform();
        double log_acceptance_ratio = (_log_marginal_likelihood+logSpeciationRatePrior(Forest::_speciation_rate))-(_prev_log_marginal_likelihood+logSpeciationRatePrior(_prev_speciation_rate));
        if (log(u) > log_acceptance_ratio){
            bool accepted = false;
            tune(accepted);
            // reject proposed speciation rate
            return "reject";
        }
        else {
            bool accepted = true;
            tune(accepted);
            return "accept";
        }
    }

    inline double Proj::logThetaPrior(double theta) {
        double exponential_rate = -log(0.05);
        return (log(exponential_rate) - theta*exponential_rate);
    }

    inline double Proj::logSpeciationRatePrior(double speciation_rate) {
        double exponential_rate = -log(0.05)/400.0;
        return (log(exponential_rate) - speciation_rate*exponential_rate);
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
        // TODO: this only works if names are in taxon^species format (no _)
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
        cout << "theta = " << Forest::_theta << endl;
        cout << "speciation rate = " << Forest::_speciation_rate << endl;

        saveAllForests(my_vec);
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
            
            // sampling loop
                // z=0 -> theta loop
                // z=1 -> speciation rate loop
            for (int z=0; z<2; z++) {
                // set starting window size for each parameter
                if (z == 0) {
                    _lambda = 0.5;
                }
                else if (z == 1) {
                    _lambda = 100.0;
                }
            // theta or speciation rate loop
                for (int i=0; i<100; i++) {
                    vector<Particle> my_vec_1(nparticles);
                    vector<Particle> my_vec_2(nparticles);
                    vector<Particle> &my_vec = my_vec_1;
                    bool use_first = true;
                    for (auto & p:my_vec ) {
                        p.setData(_data, _taxon_map);
                        p.mapSpecies(_taxon_map, _species_names);
                    }
                    
                    _prev_log_marginal_likelihood = _log_marginal_likelihood;
                    
                    // reset marginal likelihood
                    _log_marginal_likelihood = 0.0;
                    //run through each generation of particles
                    for (unsigned g=0; g<nspecies; g++){

                        vector<double> log_weight_vec;
                        double log_weight = 0.0;

                        //taxon joining and reweighting step
                        for (auto & p:my_vec) {
                            log_weight = p.proposal();
                            log_weight_vec.push_back(log_weight);
                        }
                        
                        double ess_inverse = 0.0;
                        normalizeWeights(my_vec);
                        
                        for (auto & p:my_vec) {
                            ess_inverse += exp(2.0*p.getLogWeight());
                        }
                        
                        double ess = 1.0/ess_inverse;
//                        cout << "ESS is " << ess << " " << "g = " << g << endl;
                        
                        if (ess < 100) {
                            resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1);
                            //if use_first is true, my_vec = my_vec_2
                            //if use_first is false, my_vec = my_vec_1
                            
                            my_vec = use_first ? my_vec_2:my_vec_1;

                            //change use_first from true to false or false to true
                            use_first = !use_first;
                        }
                        resetWeights(my_vec);
                    } // g loop
                    
                    // theta proposal
                    if (z == 0) {
                        if (i == 0) {
                            _prev_particles=my_vec;
                        }
                        // propose new value of theta
                        if (_prev_log_marginal_likelihood != 0.0) {
                            cout << "\n" << "previous theta: " << _prev_theta << "\t" << "proposed theta: " << Forest::_theta << endl;
                            cout << "\t" << "proposed marg like: " << _log_marginal_likelihood;
                            cout << "\t" << "prev marg like: " << _prev_log_marginal_likelihood;
                            string outcome = acceptTheta();
                            if (outcome == "reject") {
                                my_vec = _prev_particles;
                                _log_marginal_likelihood = _prev_log_marginal_likelihood;
                                Forest::_theta = _prev_theta;
                                cout << "\n" << "REJECT" << endl;
                                cout << "lambda: " << _lambda << endl;
                                _theta_vector.push_back(Forest::_theta);
                            }
                            else {
                                _prev_particles = my_vec;
                                cout << "\n" << "ACCEPT" << endl;
                                cout << "lambda: " << _lambda << endl;
                                _theta_vector.push_back(Forest::_theta);
                                _theta_accepted_number++;
                            }
                        }
                        // propose new theta for all steps but last step
                        if (i<99) {
                            proposeTheta();
                        }
                        _accepted_particle_vec = my_vec;
                    }
                    if (z == 1) {
                        if (i == 0) {
                            _prev_particles = my_vec;
                            _prev_log_marginal_likelihood = 0.0;
                        }
                        // propose new value of speciation rate
                        if (_prev_log_marginal_likelihood != 0.0) {
                            cout << "\n" << "previous speciation rate: " << _prev_speciation_rate << "\t" << "proposed speciation rate: " << Forest::_speciation_rate << endl;
                            cout << "\t" << "proposed marg like: " << _log_marginal_likelihood;
                            cout << "\t" << "prev marg like: " << _prev_log_marginal_likelihood;
                            string outcome = acceptSpeciationRate();
                            if (outcome == "reject") {
                                my_vec = _prev_particles;
                                _log_marginal_likelihood = _prev_log_marginal_likelihood;
                                Forest::_speciation_rate = _prev_speciation_rate;
                                cout << "\n" << "REJECT" << endl;
                                cout << "lambda: " << _lambda << endl;
                                
                                _speciation_rate_vector.push_back(Forest::_speciation_rate);
                            }
                            else {
                                _prev_particles = my_vec;
                                _speciation_rate_vector.push_back(Forest::_speciation_rate);
                                _speciation_rate_accepted_number++;
                                cout << "\n" << "ACCEPT" << endl;
                                cout << "lambda: " << _lambda << endl;
                            }
                        }
                        if (i<99) {
                            proposeSpeciationRate();
                        }
                        _accepted_particle_vec = my_vec;
                    }
                } // i loop
            } // z loop
            showFinal(_accepted_particle_vec);
            
            cout << "________________________________________" << endl;
            
            for (auto &i:_theta_vector) {
                cout << i << "\n";
            }
            // sort sampled thetas by frequency
//            cout << "~theta, ~freq, " << endl;
//            sort(_theta_vector.begin(), _theta_vector.end(), greater<double>());
//            int b = 1;
//            for (int a = 0; a < _theta_vector.size(); a++) {
//                if (_theta_vector[a] == _theta_vector[a+1]) {
//                    b++;
//                }
//                else {
////                    cout << "theta: " << _theta_vector[a] << "\t" << "\t" << "frequency: " << b << endl;
//                    cout << _theta_vector[a] << " , " << b << " , " << endl;
//                    b = 1;
//                }
//                a++;
//            }
            
            cout << "\n" << "________________________________________" << endl;
            cout << "number of accepted theta proposals: " << _theta_accepted_number << endl;
            
            cout << "________________________________________" << endl;
            for (auto &i:_speciation_rate_vector) {
                cout << i << "\n";
            }
            // sort sampled speciation rates by frequency
//            cout << "~speciation_rate, ~freq, " << endl;
//            sort(_speciation_rate_vector.begin(), _speciation_rate_vector.end(), greater<double>());
//            int n = 1;
//            for (int m = 0; m < _speciation_rate_vector.size(); m++) {
//                if (_speciation_rate_vector[m] == _speciation_rate_vector[m+1]) {
//                    n++;
//                }
//                else {
////                    cout << "speciation rate: " << _speciation_rate_vector[m] << "\t" << "\t" << "frequency: " << n << endl;
//                    cout << _speciation_rate_vector[m] << " , " << n << " , " << endl;
//                    n = 1;
//                }
//                m++;
//            }
            
            cout << "\n" << "________________________________________" << endl;
            cout << "number of accepted speciation rate proposals: " << _speciation_rate_accepted_number << endl;
        }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }
}
