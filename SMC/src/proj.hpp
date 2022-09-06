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

            void                normalizeWeights(vector<Particle::SharedPtr> & particles);
            unsigned            chooseRandomParticle(vector<Particle::SharedPtr> & particles, vector<double> & cum_prob);
            void                resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
            void                resetWeights(vector<Particle::SharedPtr> & particles);
            double              getWeightAverage(vector<double> log_weight_vec);
            void                createSpeciesMap(Data::SharedPtr);
            void                showParticlesByWeight(vector<Particle::SharedPtr> my_vec);
            void                proposeTheta();
            void                proposeSpeciationRate();
            double              logThetaPrior(double theta);
            double              logSpeciationRatePrior(double speciation_rate);
            string              acceptTheta();
            string              acceptSpeciationRate();
            void                showFinal(vector<Particle::SharedPtr>);
            void                tune(bool accepted);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles);
            void                proposeParticles(vector<Particle::SharedPtr> &particles);
            void                printSpeciationRates();
            void                printThetas();
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
            double                      _prev_theta = 0.0;
            double                      _prev_speciation_rate = 0.0;
            vector<Particle::SharedPtr>            _accepted_particle_vec;
            vector<Particle::SharedPtr>            _prev_particles;
            vector<pair<double, double>>  _theta_vector;
            vector<pair<double, double>>  _speciation_rate_vector;
            double                      _theta_accepted_number = 0.0;
            double                      _speciation_rate_accepted_number = 0.0;
            int                         _nattempts = 0;
            bool                        _tuning = true;
            double                      _target_acceptance = 0.3;
            double                      _lambda;
            unsigned                    _nthreads;
            bool                        _estimate_theta;
            bool                        _estimate_speciation_rate;
            unsigned                    _nsamples; // number of total samples
            unsigned                    _sample = 0; // index of current sample
            void                        handleBaseFrequencies();
            void                        debugSpeciesTree(vector<Particle::SharedPtr> &particles);
            void                        estimateTheta(vector<Particle::SharedPtr> &particles);
            void                        estimateSpeciationRate(vector<Particle::SharedPtr> &particles);
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
        ("speciation_rate", boost::program_options::value(&Forest::_speciation_rate)->default_value(1), "speciation rate")
        ("proposal",  boost::program_options::value(&Forest::_proposal)->default_value("prior-post"), "a string defining a proposal (prior-prior or prior-post)")
        ("model", boost::program_options::value(&Forest::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&Forest::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("nthreads",  boost::program_options::value(&_nthreads)->default_value(1.0), "number of threads for multi threading")
        ("estimate_theta", boost::program_options::value(&_estimate_theta)->default_value(false), "bool: true if theta estimated, false if empirical")
        ("estimate_speciation_rate", boost::program_options::value(&_estimate_speciation_rate)->default_value(false), "bool: true if speciation rate estimated, false if empirical")
        ("nsamples", boost::program_options::value(&_nsamples)->default_value(1.0), "number of samples if parameters are being estimated")
        ("migration_rate", boost::program_options::value(&Forest::_migration_rate)->default_value(0.0), "migration rate")
        ("hybridization_rate", boost::program_options::value(&Forest::_hybridization_rate)->default_value(0.0), "hybridization rate")
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
        cout << sum-1 << endl;
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
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            log_weight_vec[i++] = p->getLogWeight();
        }

        double log_particle_sum = getRunningSum(log_weight_vec);

        for (auto & p : particles) {
            p->setLogWeight(p->getLogWeight() - log_particle_sum);
        }
        
        _log_marginal_likelihood += log_particle_sum - log(_nparticles);
        sort(particles.begin(), particles.end(), greater<Particle::SharedPtr>());
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

    inline void Proj::estimateTheta(vector<Particle::SharedPtr> &particles) {
        if (_sample == 0) {
            _theta_vector.push_back(make_pair(Forest::_theta, _log_marginal_likelihood));
            _prev_particles = particles;
            _prev_log_marginal_likelihood = 0.0;
        }
        // propose new value of theta
        if (_prev_log_marginal_likelihood != 0.0) {
//            cout << "\n" << "previous theta: " << _prev_theta << "\t" << "proposed theta: " << Forest::_theta << endl;
//            cout << "\t" << "proposed marg like: " << _log_marginal_likelihood;
//            cout << "\t" << "prev marg like: " << _prev_log_marginal_likelihood << endl;
            string outcome = acceptTheta();
            if (outcome == "reject") {
                particles = _prev_particles;
                _log_marginal_likelihood = _prev_log_marginal_likelihood;
                Forest::_theta = _prev_theta;
//                cout << "lambda: " << _lambda << endl;
                _theta_vector.push_back(make_pair(Forest::_theta, _log_marginal_likelihood));
            }
            else {
                _prev_particles = particles;
//                cout << "\n" << "ACCEPT" << endl;
//                cout << "lambda: " << _lambda << endl;
                _theta_vector.push_back(make_pair(Forest::_theta, _log_marginal_likelihood));
                _theta_accepted_number++;
            }
        }
        // propose new theta for all steps but last step
        if (_sample<_nsamples) {
            proposeTheta();
        }
        _accepted_particle_vec = particles;
    }

    inline void Proj::estimateSpeciationRate(vector<Particle::SharedPtr> &particles){
        if (_sample == 0) {
            _prev_particles = particles;
            _prev_log_marginal_likelihood = 0.0;
            _speciation_rate_vector.push_back(make_pair(Forest::_speciation_rate, _log_marginal_likelihood));
        }
        // propose new value of speciation rate
        if (_prev_log_marginal_likelihood != 0.0) {
//            cout << "\n" << "previous speciation rate: " << _prev_speciation_rate << "\t" << "proposed speciation rate: " << Forest::_speciation_rate << endl;
//            cout << "\t" << "proposed marg like: " << _log_marginal_likelihood;
//            cout << "\t" << "prev marg like: " << _prev_log_marginal_likelihood << endl;
            string outcome = acceptSpeciationRate();
            if (outcome == "reject") {
                particles = _prev_particles;
                _log_marginal_likelihood = _prev_log_marginal_likelihood;
                Forest::_speciation_rate = _prev_speciation_rate;
    //          cout << "\n" << "REJECT" << endl;
                _speciation_rate_vector.push_back(make_pair(Forest::_speciation_rate, _log_marginal_likelihood));
            }
            else {
                _prev_particles = particles;
                _speciation_rate_vector.push_back(make_pair(Forest::_speciation_rate, _log_marginal_likelihood));
                _speciation_rate_accepted_number++;
//                cout << "\n" << "ACCEPT" << endl;
//                cout << "lambda: " << _lambda << endl;
            }
        }
        if (_sample<_nsamples) {
            proposeSpeciationRate();
        }
        _accepted_particle_vec = particles;
    }

    inline void Proj::proposeTheta() {
        double u = rng.uniform();
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

    inline unsigned Proj::chooseRandomParticle(vector<Particle::SharedPtr> & particles, vector<double> & cum_probs) {
        int chosen_index = -1;
        unsigned nparticles = (unsigned)particles.size();
        double u = rng.uniform();
        double cum_prob = 0.0;
        
        for(unsigned j = 0; j < nparticles; j++) {
            if (cum_probs[j]<0.0){
                cum_prob += exp(particles[j]->getLogWeight());
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

    inline void Proj::resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles) {
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
                Particle::SharedPtr p0 = from_particles[i];
                // dereference p0, create a new particle, create new shared pointer to that particle
//                Particle::SharedPtr p = Particle::SharedPtr(new Particle(*p0));
                to_particles[m++]=Particle::SharedPtr(new Particle(*p0));
            }
        }
        assert(nparticles == to_particles.size());
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
        cout << "speciation rate = " << Forest::_speciation_rate << endl;

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

    inline void Proj::showParticlesByWeight(vector<Particle::SharedPtr> my_vec) {
        vector <double> weights;
        
        //create weight vector
        for (auto & p:my_vec) {
            weights.push_back(p->getLogWeight());
        }
        
        //sort particles by weight
        sort(my_vec.begin(), my_vec.end(), greater<Particle::SharedPtr>());
        
        //print first particle
        cout << "\n" << "Heaviest particle: ";
        my_vec[0]->showParticle();
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

    inline void Proj::printSpeciationRates() {
        // print sampled speciation rates and marginal likelihood in csv format
        cout << "theta, marginal_likelihood " << endl;
        for (auto &p:_speciation_rate_vector) {
            cout << p.first << ", " << p.second << endl;
        }
    }

    inline void Proj::printThetas() {
        // print sampled thetas and marginal likelihood in csv format
        cout << "theta, marginal_likelihood " << endl;
        for (auto &p:_theta_vector) {
            cout << p.first << ", " << p.second << endl;
        }
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
            _data = Data::SharedPtr(new Data());
            _data->setPartition(_partition);
            _data->getDataFromFile(_data_file_name);

            summarizeData(_data);
            createSpeciesMap(_data);

            //set number of species to number in data file
            setNumberTaxa(_data);
            unsigned nspecies = (unsigned) _species_names.size();
            Forest::setNumSpecies(nspecies);
            rng.setSeed(_random_seed);

//          create vector of particles
            unsigned nparticles = _nparticles;
            
            unsigned nsubsets = _data->getNumSubsets();
            Particle::setNumSubsets(nsubsets);
            
            // TODO: sampling will only work if the sampled params are only theta and / or speciation rate
            // estimate neither theta nor speciation rate OR estimate one but not the other
            double number_of_sampling_loops = 1.0;

            if (_estimate_theta == true && _estimate_speciation_rate == true) {
                // estimate both
                number_of_sampling_loops = 2.0;
            }
            
            // only go through z loop once if no sampling
            for (int z = 0; z<number_of_sampling_loops; z++) {
                
                // sampling both theta and speciation rate
                // set starting window size
                if (number_of_sampling_loops == 2.0) {
                    if (z == 0) {_lambda = 0.001;}
                    if (z == 1) {_lambda = 25.0;}
                }
                
                // sampling just theta
                // set starting window size
                if (_estimate_theta && !_estimate_speciation_rate){
                    _lambda = 0.001;
                }
                
                // sampling just speciation rate
                // set starting window size
                else if (_estimate_speciation_rate && !_estimate_theta) {
                    _lambda = 25.0;
                }
            // loop for number of samples (either theta or speciation rate)
                for (_sample=0; _sample<_nsamples; _sample++) {
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
                    //run through each generation of particles
                    
                    for (unsigned g=0; g<nspecies; g++){
                        //taxon joining and reweighting step
                        // save pointers to current particles ("parents" of particles about to be created)
//                        vector<Particle*> parent_particles;
//                        for (auto &p:my_vec) {
//                            parent_particles.push_back(&p);
//                        }
                        proposeParticles(my_vec);
                        
                        double ess_inverse = 0.0;
                        normalizeWeights(my_vec);
                        
                        for (auto & p:my_vec) {
                            ess_inverse += exp(2.0*p->getLogWeight());
                        }
                        
                        double ess = 1.0/ess_inverse;
                        
                        if (ess < 100) {
                            // save particle weights
//                            vector<double> particle_weights;
//                            for (auto &p:my_vec) {
//                                particle_weights.push_back(p.getLogWeight());
//                            }
                            
                            // save random seeds
//                            vector<map<int, vector<double>>> random_seeds;
//                            for (auto &p:my_vec) {
//                                random_seeds.push_back(p.getRandomSeeds());
//                            }
                            resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1);
                            //if use_first is true, my_vec = my_vec_2
                            //if use_first is false, my_vec = my_vec_1
                            
                            my_vec = use_first ? my_vec_2:my_vec_1;

                            //change use_first from true to false or false to true
                            use_first = !use_first;
                        }
                        resetWeights(my_vec);
                        _accepted_particle_vec = my_vec;
                        
//                        if (g == nspecies-1) {
//                            my_vec[0].showParticle();
//                        }
                    } // g loop
                    for (auto &p:my_vec) {
                        p->showHybridNodes();
                        saveAllHybridNodes(my_vec);
                    }
                    
                    if (number_of_sampling_loops == 2.0) {
                        if (z == 0) {estimateTheta(my_vec);}
                        if (z == 1) {estimateSpeciationRate(my_vec);}
                    }
                    
                    else if (number_of_sampling_loops == 1.0) {
                        if (_estimate_theta) {estimateTheta(my_vec);}
                        else if (_estimate_speciation_rate) {estimateSpeciationRate(my_vec);}
                    }
                } // _nsamples loop - number of samples
            } // z loop - theta or speciation rate
//            showFinal(_accepted_particle_vec);
            cout << "marginal likelihood: " << _log_marginal_likelihood << endl;
            
            if (_estimate_theta) {
                cout << "number of accepted theta proposals: " << _theta_accepted_number << endl;
                printThetas();
                cout << "Theta: " << _theta_vector[_nsamples-1].first << endl;
            }
            
            if (_estimate_speciation_rate) {
                cout << "number of accepted speciation rate proposals: " << _speciation_rate_accepted_number << endl;
                printSpeciationRates();
                cout << "Speciation rate: " << _speciation_rate_vector[_nsamples-1].first << endl;
            }
        }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }
}
