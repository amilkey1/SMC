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

//#include "node_manager.hpp"
//extern proj::NodeManager nm;

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
            void                saveSpeciesTrees();
            void                saveGeneTrees(unsigned ngenes);
            void                saveGeneTree(unsigned gene_number);
            void                writeLoradFile(unsigned ngenes, unsigned nspecies, unsigned ntaxa);
            void                normalizeWeights(vector<Particle::SharedPtr> & particles);
            unsigned            chooseRandomParticle(vector<Particle::SharedPtr> & particles, vector<double> & cum_prob);
            void                resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
            void                resetWeights(vector<Particle::SharedPtr> & particles);
            double              getWeightAverage(vector<double> log_weight_vec);
            void                createSpeciesMap(Data::SharedPtr);
//            void                showParticlesByWeight(vector<Particle::SharedPtr> my_vec);
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
            double                      _small_enough;
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

    inline void Proj::writeLoradFile(unsigned ngenes, unsigned nspecies, unsigned ntaxa) {
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
        for (auto &p:_accepted_particle_vec) {
            logf << iter;
            iter++;
            
            logf << "\t" << p->getLogLikelihood();
            
            for (unsigned g=0; g<ngenes+1; g++) {
                logf << "\t" << p->getTopologyPrior(g);
            }
//            vector<pair<double, double>> increments_and_priors;
            for (unsigned g=0; g<ngenes+1; g++) {
                for (auto &b:p->getIncrementPriors(g)) {
//                    increments_and_priors.push_back(b);
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

    inline void Proj::saveSpeciesTrees() {
            ofstream treef("species_trees.trees");
            treef << "#nexus\n\n";
            treef << "begin trees;\n";
            for (auto &p:_accepted_particle_vec) {
                treef << "  tree test = [&R] " << p->saveForestNewick()  << ";\n";
            }
            treef << "end;\n";
            treef.close();
        }

    inline void Proj::saveGeneTrees(unsigned ngenes) {
        ofstream treef("gene_trees.txt");
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:_accepted_particle_vec) {
                for (int i=1; i<ngenes+1; i++) {
                    treef << "  gene " << i << " = [&R] " << p->saveGeneNewick(i)  << ";\n";
            }
            treef << endl;
        }
        treef << "end;\n";
        treef.close();
    }

    inline void Proj::saveGeneTree (unsigned gene_number) {
        string name = "gene" + to_string(gene_number) + ".trees";
        ofstream treef(name);
        treef << "#nexus\n\n";
        treef << "begin trees;\n";
        for (auto &p:_accepted_particle_vec) {
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
        unsigned species_joins = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            if (p->speciesJoinProposed()) {
                species_joins++;
            }
            log_weight_vec[i++] = p->getLogWeight();
        }

        double log_particle_sum = getRunningSum(log_weight_vec);
        _log_marginal_likelihood += log_particle_sum - log(_nparticles);
        
#if defined(SPECIES_TREE_WEIGHT_AVERAGE)
        i = 0;
        species_joins = 0;
        double weight_average = log_particle_sum - _nparticles + species_joins;
        for (auto &p:particles) {
            if (p->speciesJoinProposed()) {
//                assert (p->getLogWeight() == 0);
                p->setLogWeight(weight_average);
                log_weight_vec[i] = p->getLogWeight();
            }
            i++;
        }
        
        log_particle_sum = getRunningSum(log_weight_vec);
#endif

        for (auto & p : particles) {
            p->setLogWeight(p->getLogWeight() - log_particle_sum);
        }
        
//        _log_marginal_likelihood += log_particle_sum - log(_nparticles);
//        sort(particles.begin(), particles.end(), greater<Particle>());
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
//
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
            unsigned ntaxa = setNumberTaxa(_data);
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
                    for (auto &p:my_vec) {
                        p->calcLogLikelihood();
                    }
                    normalizeWeights(my_vec); // initialize marginal likelihood
                    
                    //run through each generation of particles
                    
                    for (unsigned g=0; g<(ntaxa-1)*nsubsets+nspecies-1; g++){
                        cout << "starting step " << g << " of " << nspecies-1+(ntaxa-1)*nsubsets-1 << endl;
                        //taxon joining and reweighting step
                        proposeParticles(my_vec);
                        
                        unsigned num_species_particles_proposed = 0;
                        
                        double ess_inverse = 0.0;
                        for (auto &p:my_vec) {
                            if (p->speciesJoinProposed()) {
                                num_species_particles_proposed++;
                            }
//                            p->showSpeciesTree();
//                            p->showParticle();
                        }
                        normalizeWeights(my_vec);
                        
                        for (auto & p:my_vec) {
                            ess_inverse += exp(2.0*p->getLogWeight());
                        }
                        
                        double ess = 1.0/ess_inverse;
                        cout << "\t" << "ESS is : " << ess << endl;
                        
//                        if (ess < 100) {
                            resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1);
                            //if use_first is true, my_vec = my_vec_2
                            //if use_first is false, my_vec = my_vec_1
                            
                            my_vec = use_first ? my_vec_2:my_vec_1;

                            //change use_first from true to false or false to true
                            use_first = !use_first;
//                        }
                        unsigned species_count = 0;
                        
//                        for (auto &p:my_vec) {
//                            p->showParticle();
//                        }
                        
                        for (auto &p:my_vec) {
                            if (p->speciesJoinProposed()) {
                                species_count++;
//                                p->showParticle();
//                                p->showSpeciesTree();
//                                cout << "ESS for a species join is " << ess << endl;
//                                break;
                            }
                        }
                        cout << "\t" << "number of species join particles proposed = " << num_species_particles_proposed << endl;
                        cout << "\t" << "number of species join particles accepted = " << species_count << endl;
                        resetWeights(my_vec);
                        _accepted_particle_vec = my_vec;
                        
                    } // g loop
                    
                    saveAllHybridNodes(my_vec);
//                    for (auto &p:my_vec) {
//                        p.showHybridNodes();
//                    }
//                    for (auto &p:my_vec) {
//                        p->showParticle();
//                    }
                    
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
            saveSpeciesTrees();
//            saveGeneTrees(nsubsets);
            for (int i=1; i < nsubsets+1; i++) {
                saveGeneTree(i);
            }
            writeLoradFile(nsubsets, nspecies, ntaxa);
            
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

