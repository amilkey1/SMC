#pragma once

#include <iostream>
#include <fstream>
#include "data.hpp"
#include "partition.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "xproj.hpp"
#include "particle.hpp"
#include <vector>
#include <thread>
#include <boost/algorithm/string/split.hpp>
#include "tree_summary.hpp"

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
            void                saveAllForests(vector<Particle::SharedPtr> &v) const ;
            void                saveParticleWeights(vector<Particle::SharedPtr> &v) const;
            void                saveParticleLikelihoods(vector<Particle::SharedPtr> &v) const;

            void                normalizeWeights(vector<Particle::SharedPtr> & particles, string a, bool calc_marg_like);
            unsigned            chooseRandomParticle(vector<Particle::SharedPtr> & particles, vector<double> & cum_prob, string a);
            void                resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a);
            void                resetWeights(vector<Particle::SharedPtr> & particles, string a);
            double              getWeightAverage(vector<double> log_weight_vec);
            void                createSpeciesMap(Data::SharedPtr);
            void                showParticlesByWeight(vector<Particle::SharedPtr> my_vec, string a);
            void                proposeTheta();
            void                proposeSpeciationRate();
            void                proposeHybridizationRate();
            void                proposeParameters();
            double              logThetaPrior(double theta);
            double              logSpeciationRatePrior(double speciation_rate);
            double              logHybridizationRatePrior(double hybridization_rate);
            string              acceptParameters();
            void                showFinal(vector<Particle::SharedPtr>);
            double              tune(bool accepted, double lambda);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles, bool gene_trees_only, bool unconstrained, string a, bool deconstruct);
            void                proposeParticles(vector<Particle::SharedPtr> &particles, bool gene_trees_only, bool unconstrained, string a, bool deconstruct);
            void                printSpeciationRates();
            void                printThetas();
            void                saveAllHybridNodes(vector<Particle::SharedPtr> &v) const;
        vector< pair<double,unsigned> >                throwDarts(vector< pair<double,unsigned> > & cum_probs, vector<Particle::SharedPtr> &particles, string a);
            void                writeGeneTreeFile();
        
        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _prev_log_marginal_likelihood = 0.0;
            bool                        _use_gpu;
            bool                        _ambig_missing;
            unsigned                    _nparticles;
            unsigned                    _random_seed;
//            double                      _avg_marg_like;
            double                      _log_marginal_likelihood;
            double                      _gene_tree_log_marginal_likelihood;
            double                      _species_tree_log_marginal_likelihood;
            bool                        _run_on_empty;
            bool                        _build_species_tree_first;
            double                      _theta_prior;
            double                      _prev_theta_prior;
            double                      _speciation_rate_prior;
            double                      _prev_speciation_rate_prior;
            double                      _hybridization_rate_prior;
            double                      _prev_hybridization_rate_prior;

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
            double                      _prev_hybridization_rate = 0.0;
            vector<Particle::SharedPtr>            _accepted_particle_vec;
            vector<Particle::SharedPtr>            _prev_particles;
            vector<pair<double, double>>  _theta_vector;
            vector<pair<double, double>>  _speciation_rate_vector;
            vector<pair<double, double>> _hybridization_rate_vector;
            double                      _theta_accepted_number = 0.0;
            double                      _speciation_rate_accepted_number = 0.0;
            double                      _hybridization_rate_accepted_number = 0.0;
            int                         _nattempts = 0;
            bool                        _tuning = true;
            double                      _target_acceptance = 0.3;
            double                      _theta_lambda;
            double                      _speciation_rate_lambda;
            double                      _hybrid_rate_lambda;
            unsigned                    _nthreads;
            double                      _hybridization_rate;
            bool                        _estimate_theta;
            bool                        _estimate_speciation_rate;
            bool                        _estimate_hybridization_rate;
            unsigned                    _nsamples; // number of total samples
            unsigned                    _sample = 0; // index of current sample
            void                        handleBaseFrequencies();
            void                        debugSpeciesTree(vector<Particle::SharedPtr> &particles);
            void                        estimateParameters(vector<Particle::SharedPtr> &particles);
            double                      _small_enough;
            int                         _niterations;
            string                      _species_newicks_name;
            string                      _gene_newicks_names;
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

    inline void Proj::saveParticleWeights(vector<Particle::SharedPtr> &v) const {
        // this function saves particle weights + species tree newick in order of weight
        //sort particles by weight
        sort(v.begin(), v.end(), greater<Particle::SharedPtr>());
        
        ofstream weightf("weights.txt");
        weightf << "begin trees;" << endl;
//        for (auto &p:v) {
        for (int i=0; i<v.size(); i++) {
//            weightf << "particle\n";
//            weightf << p->getCoalescentLikelihood() << "\t" << p->saveParticleWeights() << "\t" << p->getSpeciesNewick() << "\n";
            weightf << "tree " << i << " = " << v[i]->getSpeciesNewick() << "; " << "\n";
        }
        weightf << "end trees;" << endl;
        weightf.close();
        
        ofstream alltreesf("alltrees.txt");
        for (int i=0; i<v.size(); i++) {
            alltreesf << "particle " << i << endl;
            alltreesf << "species tree = " << v[i]->getSpeciesNewick() << "; " << "\n";
            vector<string> gene_newicks = v[i]->getGeneTreeNewicks();
            for (int g=0; g<gene_newicks.size(); g++) {
                alltreesf << "gene tree " << g+1 << " = " << gene_newicks[g] << "\n";
            }
        }
    }

    inline void Proj::saveParticleLikelihoods(vector<Particle::SharedPtr> &v) const {
        ofstream likelihoodf("likelihoods.txt");
        for (auto &p:v) {
            likelihoodf << "particle\n";
            likelihoodf << p->saveParticleLikelihoods() << "\n";
        }
        likelihoodf.close();
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
        ("theta, t", boost::program_options::value(&Forest::_starting_theta)->default_value(0.05), "theta")
        ("speciation_rate", boost::program_options::value(&Forest::_speciation_rate)->default_value(1), "speciation rate")
        ("proposal",  boost::program_options::value(&Forest::_proposal)->default_value("prior-post"), "a string defining a proposal (prior-prior or prior-post or prior-post-ish)")
        ("model", boost::program_options::value(&Forest::_model)->default_value("JC"), "a string defining a substitution model")
        ("kappa",  boost::program_options::value(&Forest::_kappa)->default_value(1.0), "value of kappa")
        ("base_frequencies", boost::program_options::value(&Forest::_string_base_frequencies)->default_value("0.25, 0.25, 0.25, 0.25"), "string of base frequencies A C G T")
        ("nthreads",  boost::program_options::value(&_nthreads)->default_value(1.0), "number of threads for multi threading")
        ("estimate_theta", boost::program_options::value(&_estimate_theta)->default_value(false), "bool: true if theta estimated, false if empirical")
        ("estimate_speciation_rate", boost::program_options::value(&_estimate_speciation_rate)->default_value(false), "bool: true if speciation rate estimated, false if empirical")
        ("estimate_hybridization_rate", boost::program_options::value(&_estimate_hybridization_rate)->default_value(false), "bool: true if hybridization rate estimated, false if empirical")
        ("nsamples", boost::program_options::value(&_nsamples)->default_value(1.0), "number of samples if parameters are being estimated")
        ("migration_rate", boost::program_options::value(&Forest::_migration_rate)->default_value(0.0), "migration rate")
        ("hybridization_rate", boost::program_options::value(&Forest::_hybridization_rate)->default_value(0.0), "hybridization rate")
        ("run_on_empty", boost::program_options::value(&_run_on_empty)->default_value(false), "run with no data")
        ("build_species_tree_first", boost::program_options::value(&_build_species_tree_first)->default_value(true), "build the complete species tree before starting the gene trees")
        ("niterations", boost::program_options::value(&_niterations)->default_value(1.0), "number of times to alternate between species and gene trees")
        ("species_newicks", boost::program_options::value(&_species_newicks_name)->default_value("null"), "name of file containing species newick descriptions")
        ("gene_newicks", boost::program_options::value(&_gene_newicks_names)->default_value("null"), "name of file containing gene newick descriptions")
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
//        cout << sum-1 << endl;
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

    inline void Proj::normalizeWeights(vector<Particle::SharedPtr> & particles, string a, bool calc_marg_like) {
        unsigned i = 0;
        vector<double> log_weight_vec(particles.size());
        for (auto & p : particles) {
            log_weight_vec[i++] = p->getLogWeight(a);
        }

        double log_particle_sum = getRunningSum(log_weight_vec);

        for (auto & p : particles) {
            p->setLogWeight(p->getLogWeight(a) - log_particle_sum, a);
        }
        
        if (calc_marg_like) {
            _log_marginal_likelihood += log_particle_sum - log(_nparticles);
            cout << setprecision(12) << "   " << _log_marginal_likelihood << endl;
        }
        else {
            _species_tree_log_marginal_likelihood += log_particle_sum - log(_nparticles);
            cout << setprecision(12) << "   " << _species_tree_log_marginal_likelihood << endl;
        }
        sort(particles.begin(), particles.end(), greater<Particle::SharedPtr>());
    }

    inline double Proj::tune(bool accepted, double lambda) {
        _nattempts++;
        if (_tuning) {
            double gamma_n = 10.0/(100.0 + (double)_nattempts);
            if (accepted)
                lambda *= 1.0 + gamma_n*(1.0 - _target_acceptance)/(2.0*_target_acceptance);
            else
                lambda *= 1.0 - gamma_n*0.5;

            // Prevent run-away increases in boldness for low-information marginal densities
            if (lambda > 1000.0)
                lambda = 1000.0;
        }
        return lambda;
    }

    inline void Proj::proposeTheta() {
        double u = rng.uniform();
        double proposed_theta = _theta_lambda*u+(Forest::_starting_theta-_theta_lambda/2.0);

        // make sure proposed theta is positive
        if (proposed_theta < 0.0 ) {
            proposed_theta*=-1;
        }
        
        _prev_theta = Forest::_starting_theta;
        Forest::_starting_theta = proposed_theta;
    }

    inline void Proj::proposeHybridizationRate() {
        double u = rng.uniform();
        double proposed_hybrid_rate = _hybrid_rate_lambda*u+(Forest::_hybridization_rate-_hybrid_rate_lambda/2.0);
        
        // make sure proposed theta is positive
        if (proposed_hybrid_rate < 0.0) {
            proposed_hybrid_rate*=-1;
        }
        
        _prev_hybridization_rate = Forest::_hybridization_rate;
        Forest::_hybridization_rate = proposed_hybrid_rate;
    }

    inline void Proj::estimateParameters(vector<Particle::SharedPtr> &particles){
        if (_sample == 0) {
            _prev_particles = particles;
            _prev_log_marginal_likelihood = 0.0;
            _speciation_rate_vector.push_back(make_pair(Forest::_speciation_rate, _log_marginal_likelihood));
            _hybridization_rate_vector.push_back(make_pair(Forest::_hybridization_rate, _log_marginal_likelihood));
            _theta_vector.push_back(make_pair(Forest::_starting_theta, _log_marginal_likelihood));
        }
        // propose new value of all parameters being estimated
        if (_prev_log_marginal_likelihood != 0.0) {
                if (_estimate_speciation_rate) {
                    cout << "\n" << "previous speciation rate: " << _prev_speciation_rate << "\t" << "proposed speciation rate: " << Forest::_speciation_rate << endl;
                    cout << "speciation rate lambda: " << _speciation_rate_lambda << endl;
            }
            if (_estimate_theta) {
                cout << "\n" << "previous theta: " << _prev_theta << "\t" << "proposed theta: " << Forest::_starting_theta << endl;
                cout << "theta lambda: " << _theta_lambda << endl;
            }
            if (_estimate_hybridization_rate) {
                cout << "\n" << "previous hybridization rate: " << _prev_hybridization_rate << "\t" << "proposed hybridization rate: " << Forest::_hybridization_rate << endl;
                cout << "hybridization rate lambda: " << _hybrid_rate_lambda << endl;
            }
            string outcome = acceptParameters();
            if (outcome == "reject") {
                cout << "\n" << "REJECT" << endl;
                particles = _prev_particles;
                _log_marginal_likelihood = _prev_log_marginal_likelihood;
                Forest::_speciation_rate = _prev_speciation_rate;
                Forest::_hybridization_rate = _prev_hybridization_rate;
                Forest::_starting_theta = _prev_theta;
                _speciation_rate_vector.push_back(make_pair(Forest::_speciation_rate, _log_marginal_likelihood));
                _theta_vector.push_back(make_pair(Forest::_starting_theta, _log_marginal_likelihood));
                _hybridization_rate_vector.push_back(make_pair(Forest::_hybridization_rate, _log_marginal_likelihood));
            }
            else {
                cout << "\n" << "ACCEPT" << endl;
                _prev_particles = particles;
                _speciation_rate_vector.push_back(make_pair(Forest::_speciation_rate, _log_marginal_likelihood));
                _theta_vector.push_back(make_pair(Forest::_starting_theta, _log_marginal_likelihood));
                _hybridization_rate_vector.push_back(make_pair(Forest::_hybridization_rate, _log_marginal_likelihood));
                if (_estimate_theta) {_theta_accepted_number++;}
                if (_estimate_speciation_rate) {_speciation_rate_accepted_number++;}
                if (_estimate_hybridization_rate) {_hybridization_rate_accepted_number++;}
            }
        }
        if (_sample<_nsamples) {
            proposeParameters();
        }
        _accepted_particle_vec = particles;
    }



    inline string Proj::acceptParameters() {
        double u = rng.uniform();
        // TODO: check this is working when params not updated
        _prev_theta_prior = _theta_prior;
        _prev_speciation_rate_prior = _speciation_rate_prior;
        _prev_hybridization_rate_prior = _hybridization_rate_prior;
        
        if (_estimate_theta) {
            _theta_prior = logThetaPrior(Forest::_starting_theta);
        }
        if (_estimate_hybridization_rate) {
            _hybridization_rate_prior = logHybridizationRatePrior(Forest::_hybridization_rate);
        }
        if (_estimate_speciation_rate) {
            _speciation_rate_prior = logSpeciationRatePrior(Forest::_speciation_rate);
        }
        double log_acceptance_ratio = (_log_marginal_likelihood+_hybridization_rate_prior+_speciation_rate_prior+_theta_prior)-(_prev_log_marginal_likelihood+_prev_hybridization_rate_prior+_theta_prior+_speciation_rate_prior);
        if (log(u) > log_acceptance_ratio){
            // reject proposed theta
            bool accepted = false;
            if (_estimate_hybridization_rate) {
                _hybrid_rate_lambda = tune(accepted, _hybrid_rate_lambda);
            }
            if (_estimate_speciation_rate) {
                _speciation_rate_lambda = tune(accepted, _speciation_rate_lambda);
            }
            if (_estimate_theta) {
                _theta_lambda = tune(accepted, _theta_lambda);
            }
            return "reject";
        }
        else {
            bool accepted = true;
            if (_estimate_hybridization_rate) {
                _hybrid_rate_lambda = tune(accepted, _hybrid_rate_lambda);
            }
            if (_estimate_speciation_rate) {
                _speciation_rate_lambda = tune(accepted, _speciation_rate_lambda);
            }
            if (_estimate_theta) {
                _theta_lambda = tune(accepted, _theta_lambda);
            }
            return "accept";
        }
    }

    inline void Proj::proposeSpeciationRate() {
        double u = rng.uniform();
        double proposed_speciation_rate = _speciation_rate_lambda*u+(Forest::_speciation_rate-_speciation_rate_lambda/2.0);
        
        // make sure proposed speciation rate is positive
        if (proposed_speciation_rate < 0.0) {
            proposed_speciation_rate*=-1;
        }
        
        _prev_speciation_rate = Forest::_speciation_rate;
        Forest::_speciation_rate = proposed_speciation_rate;
    }

    inline void Proj::proposeParameters() {
        if (_estimate_speciation_rate) {
            proposeSpeciationRate();
        }
        if (_estimate_hybridization_rate) {
            proposeHybridizationRate();
        }
        if (_estimate_theta) {
            proposeTheta();
        }
    }

    inline double Proj::logThetaPrior(double theta) {
        double exponential_rate = -log(0.05);
        return (log(exponential_rate) - theta*exponential_rate);
    }

    inline double Proj::logSpeciationRatePrior(double speciation_rate) {
        double exponential_rate = -log(0.05)/100.0;
        return (log(exponential_rate) - speciation_rate*exponential_rate);
    }

    inline double Proj::logHybridizationRatePrior(double hybridization_rate) {
        double exponential_rate = -log(0.05);
        return (log(exponential_rate) - hybridization_rate*exponential_rate);
    }

    inline unsigned Proj::chooseRandomParticle(vector<Particle::SharedPtr> & particles, vector<double> & cum_probs, string a) {
        int chosen_index = -1;
        unsigned nparticles = (unsigned)particles.size();
        double u = rng.uniform();
        double cum_prob = 0.0;
        
        for(unsigned j = 0; j < nparticles; j++) {
            if (cum_probs[j]<0.0){
                cum_prob += exp(particles[j]->getLogWeight(a));
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

    inline vector< pair<double,unsigned> > Proj::throwDarts(vector< pair<double,unsigned> > & cum_probs, vector<Particle::SharedPtr> &particles, string a) {
        // Create vector of pairs p, with p.first = log weight and p.second = particle index
        cum_probs.resize(_nparticles);
        unsigned i = 0;
        for(auto & p : particles) {
            cum_probs[i].first = p->getLogWeight(a);
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
        
        return cum_probs;
    }

    inline void Proj::writeGeneTreeFile() {
        // open log file
        ofstream genef("genetrees.txt");
        genef << "let gene_forests = [\n";
        for (int j=0; j<_accepted_particle_vec.size(); j++) {
            vector<string> names = _accepted_particle_vec[j]->getGeneTreeNames(j);
            vector<string> newicks = _accepted_particle_vec[j]->getGeneTreeNewicks();
            // name and newick include quotes
            for (int i=0; i<names.size(); i++) {
                string name = names[i];
                string newick = newicks[i];
                genef << "{name:" << name << ", relrate:1.0, w:1.0, cumw:1.0, " << newick << "},\n";
            }
        }
        genef << "];";
        genef.close();
    }

    inline void Proj::resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a) {
        unsigned nparticles = (unsigned)from_particles.size();
        vector<pair<double, double>> cum_probs;
            // Create vector of pairs p, with p.first = log weight and p.second = particle index
            cum_probs.resize(nparticles);
            unsigned i = 0;
            for(auto & p : from_particles) {
                cum_probs[i].first = p->getLogWeight(a);
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
            assert(nparticles == to_particles.size());
        }
    }

    inline void Proj::resetWeights(vector<Particle::SharedPtr> & particles, string a) {
        double logw = -log(particles.size());
        for (auto & p : particles) {
            p->setLogWeight(logw, a);
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
//        cout << "mean height equals " << sum_h << endl;
////        cout << "log marginal likelihood = " << setprecision(12) << _log_marginal_likelihood << endl;
//        cout << "gene tree marg like: " << _gene_tree_log_marginal_likelihood << endl;
        cout << "species tree marg like: " << _species_tree_log_marginal_likelihood << endl;
//        cout << "starting theta = " << Forest::_starting_theta << endl;
//        cout << "speciation rate = " << Forest::_speciation_rate << endl;
//        cout << "hybridization rate = " << Forest::_hybridization_rate << endl;
        
        cout << "mean height equals " << sum_h << endl;
        cout << "log marginal likelihood = " << setprecision(12) << _log_marginal_likelihood << endl;
        cout << "starting theta = " << Forest::_starting_theta << endl;
        cout << "speciation rate = " << Forest::_speciation_rate << endl;
        cout << "hybridization rate = " << Forest::_hybridization_rate << endl;
    }

    inline void Proj::proposeParticles(vector<Particle::SharedPtr> &particles, bool gene_trees_only, bool unconstrained, string a, bool deconstruct) {
        assert(_nthreads > 0);
        if (_nthreads == 1) {
          for (auto & p : particles) {
              p->proposal(gene_trees_only, unconstrained, deconstruct);
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
              threads.push_back(thread(&Proj::proposeParticleRange, this, first, last, std::ref(particles), gene_trees_only, unconstrained, a, deconstruct));
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

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles, bool gene_trees_only, bool unconstrained, string a, bool deconstruct) {
        for (unsigned i=first; i<last; i++){
            particles[i]->proposal(gene_trees_only, unconstrained, deconstruct);
        }
    }

    inline void Proj::showParticlesByWeight(vector<Particle::SharedPtr> my_vec, string a) {
        vector <double> weights;
        
        //create weight vector
        for (auto & p:my_vec) {
            weights.push_back(p->getLogWeight(a));
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
        cout << "Starting Theta: " << Forest::_starting_theta << endl;
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
                
            // sampling both theta and speciation rate
            // set starting window size
            _prev_theta = Forest::_starting_theta;
            _prev_speciation_rate = Forest::_speciation_rate;
            _prev_hybridization_rate = Forest::_hybridization_rate;
            _theta_prior = logThetaPrior(_prev_theta);
            _speciation_rate_prior = logSpeciationRatePrior(_prev_speciation_rate);
            _hybridization_rate_prior = logHybridizationRatePrior(_prev_hybridization_rate);
            
            if (_estimate_theta) {_theta_lambda = 0.1;}
            if (_estimate_speciation_rate) {_speciation_rate_lambda = 5.0;}
            if (_estimate_hybridization_rate) {_hybrid_rate_lambda = 0.003;}
            
            // open log file
            ofstream logf("params.log");
            
        // loop for number of samples (either theta or speciation rate)
            for (_sample=0; _sample<_nsamples; _sample++) {
                cout << "sample: " << _sample << endl;
                vector<Particle::SharedPtr> my_vec_1(nparticles);
                vector<Particle::SharedPtr> my_vec_2(nparticles);
                vector<Particle::SharedPtr> &my_vec = my_vec_1;
                
                _prev_log_marginal_likelihood = _log_marginal_likelihood;
                _log_marginal_likelihood = 0.0;
                
                for (unsigned i=0; i<nparticles; i++) {
                    my_vec_1[i] = Particle::SharedPtr(new Particle);
                    my_vec_2[i] = Particle::SharedPtr(new Particle);
                }

                bool use_first = true;
                int h = 0;
                
                vector<string> newicks;
                if (_gene_newicks_names != "null" && _species_newicks_name != "null") {
                    throw XProj(boost::str(boost::format("cannot specify gene newicks and species newicks; choose one")));

                }
                
                if (_gene_newicks_names != "null") {
                    ifstream infile(_species_newicks_name);
                    string newick;
//                    vector<string> newicks;
                    while (getline(infile, newick)) {
                        newicks.push_back(newick);
                    }
                }
                else if (_species_newicks_name != "null") {
                    ifstream infile(_species_newicks_name);
                    string newick;
//                    vector<string> newicks;
                    while (getline(infile, newick)) {
                        newicks.push_back(newick);
                    }
                }
                
                    for (auto & p:my_vec ) {
                        p->setData(_data, _taxon_map);
                        p->mapSpecies(_taxon_map, _species_names);
                        p->setParticleGeneration(-1);
//                        if (newicks.size() > 0) {
                            if (_gene_newicks_names != "null") {
                                p->processGeneNewicks(newicks);
                            }
                            else {
                                p->processSpeciesNewick(newicks); // if no newick specified, program will sample from species tree prior
                            }
//                        }
                        p->mapSpecies(_taxon_map, _species_names);
                        if (!_run_on_empty) {
                            p->setParticleGeneration(0);
                            p->setLogLikelihood(0.0);
                            p->setLogWeight(0.0, "g");
                            p->setLogWeight(0.0, "s");
                        }
                        else {
                            p->setParticleGeneration(0);
                            p->setLogLikelihood(0.0);
                        }
                        h++;
                    }
                
                for (auto &p:my_vec) {
                    p->setRunOnEmpty(_run_on_empty);
                    p->setBuildEntireSpeciesTree(_build_species_tree_first);
                }
                
//                normalizeWeights(my_vec, -1, true);
                
                //run through each generation of particles
                int ntaxa = (int) _taxon_map.size();
                bool skip = false;
                bool deconstruct = true;
//                p->setLogCoalescentLikelihood(0.0);
                for (int i=0; i<_niterations; i++) {
                    _log_marginal_likelihood = 0.0;
                    cout << "beginning iteration: " << i << endl;
                    if (i > 0) {
                        for (auto &p:my_vec) {
                            p->setParticleGeneration(0);
                            p->setLogLikelihood(0.0);
//                            p->setLogCoalescentLikelihood(0.0);
                            p->setLogWeight(0.0, "g");
                            p->setLogWeight(0.0, "s");
//                            p->mapSpecies(_taxon_map, _species_names);
                            p->mapGeneTrees(_taxon_map, _species_names);
                            p->resetGeneIncrements();
                        }
                    }
                    // keep the species partition for the gene forests at this stage but clear the tree structure
                    if (!skip) {
                        bool no_newick = true;
                        if (!no_newick) {
                            if (i == 1) {
                                for (auto &p:my_vec) {
                                    p->remakeGeneTrees(_taxon_map);
                                    p->resetGeneTreePartials(_data, _taxon_map);
                                    deconstruct = false;
                                }
                            }
                        }
                    for (unsigned g=0; g<ntaxa-1; g++){
    //                    cout << "gen " << g << endl;
                        //taxon joining and reweighting step
                        
                        bool gene_trees_only = true;
                        bool unconstrained = true;
                        if (i > 0) {
                            unconstrained = false;
                        }
                        proposeParticles(my_vec, gene_trees_only, unconstrained, "g", deconstruct);
                        deconstruct = true;
                        
                        if (!_run_on_empty) {
                            bool calc_marg_like = true;
                            
                            normalizeWeights(my_vec, "g", calc_marg_like);
                            
                            double ess_inverse = 0.0;
                            
                            for (auto & p:my_vec) {
                                ess_inverse += exp(2.0*p->getLogWeight("g"));
                            }

                            double ess = 1.0/ess_inverse;
                            cout << "   " << "ESS = " << ess << endl;
                         
                            resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1, "g");
                            //if use_first is true, my_vec = my_vec_2
                            //if use_first is false, my_vec = my_vec_1
                            
                            my_vec = use_first ? my_vec_2:my_vec_1;

                            //change use_first from true to false or false to true
                            use_first = !use_first;
                            saveParticleWeights(my_vec);
                        }
                        resetWeights(my_vec, "g");
                        _accepted_particle_vec = my_vec;
                    } // g loop
                    }
                    skip = false;
                    // filter species trees now
                if (i < _niterations-1) {
                    _species_tree_log_marginal_likelihood = 0.0;
                    for (auto &p:my_vec) {
                        // reset forest species partitions
                        p->mapSpecies(_taxon_map, _species_names);
                        p->resetSpecies();
                    }
//                    _gene_tree_log_marginal_likelihood = _log_marginal_likelihood;
//                    _log_marginal_likelihood = 0.0;
                    
                        for (unsigned s=0; s<nspecies; s++){
                            cout << "beginning species tree proposals" << endl;
                            //taxon joining and reweighting step
                            
                            bool gene_trees_only = false;
                            bool unconstrained = false;
                            
                            proposeParticles(my_vec, gene_trees_only, unconstrained, "s", deconstruct);
                            
                            if (!_run_on_empty) {
                                bool calc_marg_like = false;
                                
                                normalizeWeights(my_vec, "s", calc_marg_like);
                                
                                double ess_inverse = 0.0;
                                
                                for (auto & p:my_vec) {
                                    ess_inverse += exp(2.0*p->getLogWeight("s"));
//                                    p->showParticle();
                                }

                                double ess = 1.0/ess_inverse;
                                cout << "   " << "ESS = " << ess << endl;
                             
                                resampleParticles(my_vec, use_first ? my_vec_2:my_vec_1, "s");
                                //if use_first is true, my_vec = my_vec_2
                                //if use_first is false, my_vec = my_vec_1
                                
                                my_vec = use_first ? my_vec_2:my_vec_1;

                                //change use_first from true to false or false to true
                                use_first = !use_first;
                                saveParticleWeights(my_vec);
                            }
                            resetWeights(my_vec, "s");
                            _accepted_particle_vec = my_vec;
                        } // s loop
                    }
                }
                    
//                }
                
//                cout << "gene tree marg like: " << _gene_tree_log_marginal_likelihood << endl;
//                cout << "species tree marg like: " << _log_marginal_likelihood << endl;
//                _log_marginal_likelihood = gene_tree_marg_like - _log_marginal_likelihood;
                for (auto &p:my_vec) {
                    p->getTopologyPriors();
                }
                
                cout << "\t" << "proposed marg like: " << _log_marginal_likelihood;
                cout << "\t" << "prev marg like: " << _prev_log_marginal_likelihood << endl;
                
                if (_estimate_theta || _estimate_speciation_rate || _estimate_hybridization_rate) {
                    estimateParameters(my_vec);
                }
                
                double a = 0;
                unsigned col_count = 0;
                for (auto &p:my_vec) {
                    
                    vector<double> branch_length_vec;
                    for (auto &b:p->getBranchLengths()) {
//                        branch_length_vec.push_back(log(b));
                        branch_length_vec.push_back(b);
                    }
                    
                    vector<double> prior_vec;
                    for (auto &b:p->getBranchLengthPriors()) {
                        prior_vec.push_back(b);
                    }
                    
                    vector<double> gene_tree_log_like;
                    for (auto &g:p->getGeneTreeLogLikelihoods()) {
                        gene_tree_log_like.push_back(g);
                    }
                    
                    vector<double> log_topology_priors;
                    for (auto &t:p->getTopologyPriors()) {
                        log_topology_priors.push_back(t);
                    }
                    
                    double log_coalescent_likelihood;
                    log_coalescent_likelihood = p->getCoalescentLikelihood();
                    
                    assert(branch_length_vec.size() == prior_vec.size());
                    int ngenes = p->getNGenes();
                    
                    if (col_count == 0) {
                        logf << "iter" << "\t" << "theta";
                        for (int g=0; g<ngenes; g++) {
                            logf << "\t" << "gene_tree_log_like";
                        }
                        for (int i=0; i<nspecies-1; i++) {
//                        for (int i = 0; i < branch_length_vec.size(); i++) {
                            logf << "\t" << "species_tree_increment" << "\t" << "increment_prior";
                        }
                        for (int g=0; g<ngenes; g++) {
                            for (int j=nspecies; j<nspecies+ntaxa-1; j++) {
                                logf << "\t" << "gene_tree_increment" << "\t" << "increment_prior";
                            }
                        }
                        logf << "\t" << "species_tree_topology_prior";
                        for (int g=0; g<ngenes; g++) {
                            logf << "\t" << "gene_tree_topology_prior";
                        }
                        logf << "\t" << "log_coal_like" << endl;
                    }
                    
                    logf << a << "\t" << Forest::_starting_theta;
                    double logJ = 0.0;
                    for (auto &m:branch_length_vec) {
                        logJ += branch_length_vec[m];
                    }
                    
//                    double species_logJ = 0.0;
//                    for (int s=0; s<nspecies-1; s++) {
//                        species_logJ += branch_length_vec[s];
//                    }
                    
//                    double gene_logJ = 0.0;
//                    for (int h=nspecies-1; h<nspecies+ntaxa-1; h++) {
//                        gene_logJ += branch_length_vec[h];
//                    }
                    
                    for (int g=0; g<gene_tree_log_like.size(); g++) {
                        logf << "\t" << setprecision(12) << gene_tree_log_like[g];
                    }

                    
                    for (int i=0; i<prior_vec.size(); i++) {
//                        logf << "\t" << setprecision(12) << branch_length_vec[i] << "\t" << prior_vec[i]+branch_length_vec[i];
                        logf << "\t" << setprecision(11) << branch_length_vec[i] << "\t" << prior_vec[i];
                    }
                    
                    for (int j=0; j<log_topology_priors.size(); j++) {
                        logf << "\t" << setprecision(12) << log_topology_priors[j];
                    }
                    
                    logf << "\t" << setprecision(12) << log_coalescent_likelihood;
                    
                    logf << endl;
                    a++;
                    col_count++;
                }
                
                if (_sample == _nsamples) {
                    logf.close();
                }
                } // _nsamples loop - number of samples
            writeGeneTreeFile();
//            saveParticleWeights(_accepted_particle_vec);
//            saveParticleLikelihoods(_accepted_particle_vec);
            
            if (_estimate_theta) {
                cout << "number of accepted theta proposals: " << _theta_accepted_number << endl;
                printThetas();
                cout << "\n" << "Theta: " << _theta_vector[_nsamples-1].first << endl;
            }
            
            if (_estimate_speciation_rate) {
                cout << "number of accepted speciation rate proposals: " << _speciation_rate_accepted_number << endl;
                printSpeciationRates();
                cout << "\n" << "Speciation rate: " << _speciation_rate_vector[_nsamples-1].first << endl;
            }
//            saveAllHybridNodes(_accepted_particle_vec);
            showFinal(_accepted_particle_vec);
            cout << "marg like: " << setprecision(12) << _log_marginal_likelihood << endl;
        }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }
}
