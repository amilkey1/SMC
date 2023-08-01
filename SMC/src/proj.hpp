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
            void                resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a);
            vector<int>         resampleSpeciesParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a);
            void                resetGeneParticles(vector<int> sel_indices, vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles);
            void                resetWeights(vector<Particle::SharedPtr> & particles, string a);
            void                createSpeciesMap(Data::SharedPtr);
            void                showFinal(vector<vector<Particle::SharedPtr>>);
            void                proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, bool deconstruct, vector<pair<tuple<string, string, string>, double>> species_joined);
            void                proposeSpeciesParticleRange(unsigned first, unsigned last, vector<vector<Particle::SharedPtr>> &my_vec, int s, int nspecies, int nsubsets);
            void                proposeSpeciesParticles( vector<vector<Particle::SharedPtr>> &my_vec, int s, int nspecies, int nsubsets);
            void                proposeParticles(vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, bool deconstruct, Particle::SharedPtr species_tree_particle);
            void                saveAllHybridNodes(vector<Particle::SharedPtr> &v) const;
            void                writeGeneTreeFile();
            Particle::SharedPtr chooseTree(vector<Particle::SharedPtr> species_trees, string gene_or_species);
            void                writeLoradFile(vector<vector<Particle::SharedPtr>> my_vec, int nparticles, int nsubsets, int nspecies, int ntaxa);
            void                setStartingVariables();
            void                setUpInitialData();
            void                saveGeneAndSpeciesTrees(Particle::SharedPtr species_particle, Particle::SharedPtr gene1, Particle::SharedPtr gene2, Particle::SharedPtr gene3);
        
        private:

            std::string                 _data_file_name;
            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            double                      _prev_log_marginal_likelihood = 0.0;
            bool                        _use_gpu;
            bool                        _ambig_missing;
            unsigned                    _nparticles;
            unsigned                    _random_seed;
            double                      _log_marginal_likelihood;
            double                      _gene_tree_log_marginal_likelihood;
            double                      _species_tree_log_marginal_likelihood;
            bool                        _run_on_empty;

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
            vector<vector<Particle::SharedPtr>>            _accepted_particle_vec;
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
            void                        handleBaseFrequencies();
            void                        debugSpeciesTree(vector<Particle::SharedPtr> &particles);
            double                      _small_enough;
            int                         _niterations;
            string                      _species_newicks_name;
            string                      _gene_newicks_names;
            int                         _species_particles_per_gene_particle;
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

    inline void Proj::saveGeneAndSpeciesTrees(Particle::SharedPtr species_particle, Particle::SharedPtr gene1, Particle::SharedPtr gene2, Particle::SharedPtr gene3) {
        // TODO: only works for one data set
        ofstream testf("coalescent-likelihood.txt");
        testf << "params : " << endl;
        testf << "\t" << "theta = 0.05" << endl;
        testf << "\t" << "lambda = 6.4" << endl;
        
        testf << "coalescent likelihood = " << species_particle->getCoalescentLikelihood() << endl;
        
        testf << "species tree: " << species_particle->saveForestNewick() << endl;
        testf << endl;
        testf << "gene1: " << gene1->saveForestNewick() << endl;
        testf << endl;
        testf << "gene2: " << gene2->saveForestNewick() << endl;
        testf << endl;
        testf << "gene3: " << gene3->saveForestNewick() << endl;
        
        testf.close();
        
    }

    inline void Proj::saveParticleWeights(vector<Particle::SharedPtr> &v) const {
        // this function saves particle weights + species tree newick
        
        ofstream weightf("weights.txt");
        weightf << "begin trees;" << endl;
//        for (auto &p:v) {
        for (int i=0; i<v.size(); i++) {
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
        ("migration_rate", boost::program_options::value(&Forest::_migration_rate)->default_value(0.0), "migration rate")
        ("hybridization_rate", boost::program_options::value(&Forest::_hybridization_rate)->default_value(0.0), "hybridization rate")
        ("run_on_empty", boost::program_options::value(&_run_on_empty)->default_value(false), "run with no data")
        ("niterations", boost::program_options::value(&_niterations)->default_value(1.0), "number of times to alternate between species and gene trees")
        ("gene_newicks", boost::program_options::value(&_gene_newicks_names)->default_value("null"), "name of file containing gene newick descriptions")
        ("species_newick", boost::program_options::value(&_species_newicks_name)->default_value("null"), "name of file containing species newick descriptions")
        ("species_particles_per_gene_particle", boost::program_options::value(&_species_particles_per_gene_particle)->default_value(1), "increase number of particles for species trees by this amount")
        ("outgroup", boost::program_options::value(&Forest::_outgroup)->default_value("null"), "specify outgroup in species tree")
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
        int nparticles = _nparticles;
        if (a == "s") {
            nparticles *= _species_particles_per_gene_particle;
        }
        unsigned i = 0;
        vector<double> log_weight_vec(nparticles);

        for (int p=0; p<nparticles; p++ ) {
            log_weight_vec[i++] = particles[p]->getLogWeight(a);
        }

        double log_particle_sum = getRunningSum(log_weight_vec);

        for (int p=0; p<nparticles; p++ ) {
            particles[p]->setLogWeight(particles[p]->getLogWeight(a) - log_particle_sum, a);
        }
        
        if (calc_marg_like) {
            _log_marginal_likelihood += log_particle_sum - log(_nparticles);
//            cout << setprecision(12) << "   " << _log_marginal_likelihood << endl;
        }
//        else {
//            _species_tree_log_marginal_likelihood += log_particle_sum - log(_nparticles);
//            cout << setprecision(12) << "   " << _species_tree_log_marginal_likelihood << endl;
//        }
//        sort(particles.begin(), particles.end(), greater<Particle::SharedPtr>());
    }



    inline void Proj::writeGeneTreeFile() {
        // open log file
        ofstream genef("genetrees.txt");
        genef << "let gene_forests = [\n";
        for (int p=0; p<_nparticles; p++) {
            genef << "particle " << endl;
            for (int s=1; s<_accepted_particle_vec.size(); s++) {
                genef << "gene " << s << endl;
                vector<string> newicks = _accepted_particle_vec[s][p]->getGeneTreeNewicks();
                string newick = newicks[0];
                genef << newick << endl;
            }
        }
        genef << "];";
        genef.close();
    }

    inline void Proj::resetGeneParticles(vector<int> sel_indices, vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles) {
        
         unsigned nparticles = (unsigned)from_particles.size();
        
        for (unsigned i = 0; i < nparticles; i++) {
            int sel_index = sel_indices[i];
            assert(sel_index > -1);
            Particle::SharedPtr p0 = from_particles[sel_index];
            to_particles[i]=Particle::SharedPtr(new Particle(*p0));
            assert(nparticles == to_particles.size());
        }
    }

    inline Particle::SharedPtr Proj::chooseTree(vector<Particle::SharedPtr> species_trees, string gene_or_species) {
        // get species tree weights
        vector<double> log_weights;
        int nparticles = _nparticles;
        if (gene_or_species == "s") {
            nparticles *= _species_particles_per_gene_particle;
        }
        for (int p=0; p<nparticles; p++) {
            log_weights.push_back(species_trees[p]->getLogWeight(gene_or_species));
        }
        
//        // choose a random number [0,1]
        double u = rng.uniform();
        double cum_prob = 0.0;
        int index = 0.0;
        for (int i=0; i < (int) log_weights.size(); i++) {
            cum_prob += exp(log_weights[i]);
            if (u <= cum_prob) {
                index = i;
                break;
            }
        }
//        // return particle of choice
        return species_trees[index];
    }

    inline vector<int> Proj::resampleSpeciesParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a) {
        
        vector<int> sel_indices;
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
            assert( fabs( 1.0 - cum_probs[_nparticles*_species_particles_per_gene_particle-1].first ) < Proj::_small_enough);

        
        // Draw new set of particles by sampling with replacement according to cum_probs
        to_particles.resize(_nparticles*_species_particles_per_gene_particle);
        for(unsigned i = 0; i < _nparticles*_species_particles_per_gene_particle; i++) {
        
            // Select a particle to copy to the ith slot in to_particles
            int sel_index = -1;
            double u = rng.uniform();
            for(unsigned j = 0; j < nparticles; j++) {
                if (u < cum_probs[j].first) {
                    sel_index = cum_probs[j].second;
                    break;
                }
            }
            assert(sel_index > -1);
            Particle::SharedPtr p0 = from_particles[sel_index];
            to_particles[i]=Particle::SharedPtr(new Particle(*p0));
            assert(nparticles == to_particles.size());
            
            sel_indices.push_back(sel_index);
        }
        return sel_indices;
    }

    inline void Proj::resampleParticles(vector<Particle::SharedPtr> & from_particles, vector<Particle::SharedPtr> & to_particles, string a) {
        
//        unsigned nparticles = (unsigned)from_particles.size();
        unsigned nparticles = _nparticles;
        
        vector<pair<double, double>> cum_probs;
            // Create vector of pairs p, with p.first = log weight and p.second = particle index
            cum_probs.resize(nparticles);
            unsigned i = 0;
//            for(auto & p : from_particles) {
        for (int p=0; p<_nparticles; p++) {
                cum_probs[i].first = from_particles[p]->getLogWeight(a);
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
            Particle::SharedPtr p0 = from_particles[sel_index];
            to_particles[i]=Particle::SharedPtr(new Particle(*p0));
            
            assert(nparticles == to_particles.size());
//            assert(nparticles*_species_particles_per_gene_particle == to_particles.size());
        }
    }

    inline void Proj::resetWeights(vector<Particle::SharedPtr> & particles, string a) {
        int nparticles = _nparticles;
        if (a == "s") {
            nparticles *= _species_particles_per_gene_particle;
        }
        double logw = -log(particles.size());
//        for (auto & p : particles) {
        for (int p=0; p<nparticles; p++) {
            particles[p]->setLogWeight(logw, a);
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
    
    inline void Proj::showFinal(vector<vector<Particle::SharedPtr>> my_vec) {
        // this function displays the final species trees
        for (int p=0; p<_nparticles; p++) {
                my_vec[0][p]->showParticle();
            }
        
        double sum_h = 0.0;
        for (auto & p:my_vec[0]) {
            double h = p->calcHeight();
            sum_h += h;
        }
        sum_h/=my_vec[0].size();
        cout << "species tree marg like: " << _species_tree_log_marginal_likelihood << endl;
        
        cout << "mean height equals " << sum_h << endl;
        cout << "log marginal likelihood = " << setprecision(12) << _log_marginal_likelihood << endl;
        cout << "starting theta = " << Forest::_starting_theta << endl;
        cout << "speciation rate = " << Forest::_speciation_rate << endl;
        cout << "hybridization rate = " << Forest::_hybridization_rate << endl;
    }

    inline void Proj::proposeParticles(vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, bool deconstruct, Particle::SharedPtr species_tree_particle) {
        assert(_nthreads > 0);
        vector<pair<tuple<string, string, string>, double>> species_joined = species_tree_particle->getSpeciesJoined();
        assert (species_joined.size() > 0);
        
        if (_nthreads == 1) {
            for (int p=0; p<_nparticles; p++) {
                particles[p]->proposal(gene_trees_only, deconstruct, species_joined);
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
              threads.push_back(thread(&Proj::proposeParticleRange, this, first, last, std::ref(particles), gene_trees_only, a, deconstruct, species_joined));
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

    inline void Proj::proposeParticleRange(unsigned first, unsigned last, vector<Particle::SharedPtr> &particles, bool gene_trees_only, string a, bool deconstruct, vector<pair<tuple<string, string, string>, double>> species_joined) {
        for (unsigned i=first; i<last; i++){
            particles[i]->proposal(gene_trees_only, deconstruct, species_joined);
        }
    }

    inline void Proj::proposeSpeciesParticleRange(unsigned first, unsigned last, vector<vector<Particle::SharedPtr>> &my_vec, int s, int nspecies, int nsubsets) {
            for (unsigned p=first; p<last; p++) {
                    tuple<string, string, string> species_joined = my_vec[0][p]->speciesTopologyProposal();

                    vector<double> max_depths; // this vector contains list of maximum depths for each gene tree

                    if (s < nspecies-1) {
                        for (int j=1; j<nsubsets+1; j++) {
                            max_depths.push_back(my_vec[j][p]->calcConstrainedProposal(species_joined));
                        }

                        // now finish the species tree branch length proposal
                        my_vec[0][p]->speciesProposal(max_depths, species_joined);
                    }
                    
                    else {
                        for (int j=1; j<nsubsets+1; j++) {
                            my_vec[j][p]->calcConstrainedProposal(species_joined); // update the gene tree species partitions
                            max_depths.push_back(0.0);
                        }
                        my_vec[0][p]->speciesProposal(max_depths, species_joined); // set last edge length of species tree to 0.0
                    }

                    double log_coalescent_likelihood = 0.0;

                    // calculate coalescent likelihood for each gene on each particle
                        for (int j=1; j<nsubsets+1; j++) {
                            double last_edge_len = my_vec[0][p]->getLastEdgeLen();
                            double species_tree_height = my_vec[0][p]->getSpeciesTreeHeight();
                            log_coalescent_likelihood += my_vec[j][p]->calcGeneCoalescentLikelihood(last_edge_len, species_joined, species_tree_height);
                        }

                    my_vec[0][p]->calcSpeciesParticleWeight(log_coalescent_likelihood);
                } // p loop
    }

    inline void Proj::proposeSpeciesParticles( vector<vector<Particle::SharedPtr>> &my_vec, int s, int nspecies, int nsubsets) {
        if (_nthreads == 1) {
            for (int p=0; p<_nparticles*_species_particles_per_gene_particle; p++) {
                tuple<string, string, string> species_joined = my_vec[0][p]->speciesTopologyProposal();
                
                vector<double> max_depths; // this vector contains list of maximum depths for each gene tree
                
                if (s < nspecies-1) {
                    for (int j=1; j<nsubsets+1; j++) {
                        
                        max_depths.push_back(my_vec[j][p]->calcConstrainedProposal(species_joined));
                    }
                    
                    // now finish the species tree branch length proposal
                    my_vec[0][p]->speciesProposal(max_depths, species_joined); // branch length proposal
                }
                else {
                    for (int j=1; j<nsubsets+1; j++) {
                        my_vec[j][p]->calcConstrainedProposal(species_joined); // update the gene tree species partitions
                        max_depths.push_back(0.0);
                    }
                    my_vec[0][p]->speciesProposal(max_depths, species_joined); // set last edge length of species tree to 0.0
                }

                double log_coalescent_likelihood = 0.0;
                
                // calculate coalescent likelihood for each gene on each particle
                    for (int j=1; j<nsubsets+1; j++) {
                        double last_edge_len = my_vec[0][p]->getLastEdgeLen();
                        double species_tree_height = my_vec[0][p]->getSpeciesTreeHeight();
                        log_coalescent_likelihood += my_vec[j][p]->calcGeneCoalescentLikelihood(last_edge_len, species_joined, species_tree_height);
                        
                        cout << log_coalescent_likelihood << endl;
                    }
                
                my_vec[0][p]->calcSpeciesParticleWeight(log_coalescent_likelihood);
            } // p loop
        }
            
        else { // multithreading
            // divide up the particles as evenly as possible across threads
            unsigned first = 0;
            int nspecies_particles = _nparticles*_species_particles_per_gene_particle;
            unsigned incr = nspecies_particles/_nthreads + (nspecies_particles % _nthreads != 0 ? 1:0); // adding 1 to ensure we don't have 1 dangling particle for odd number of particles
            unsigned last = incr;


            // need a vector of threads because we have to wait for each one to finish
            vector<thread> threads;

            while (true) {
//                // create a thread to handle particles first through last - 1
                  threads.push_back(thread(&Proj::proposeSpeciesParticleRange, this, first, last, std::ref(my_vec), s, nspecies, nsubsets));
//                // update first and last
                first = last;
                last += incr;
                if (last > nspecies_particles) {
                  last = nspecies_particles;
                  }
                if (first>=nspecies_particles) {
                    break;
                }
          }

          // the join function causes this loop to pause until the ith thread finishes
              for (unsigned i = 0; i < threads.size(); i++) {
                threads[i].join();
              }
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

    inline void Proj::setUpInitialData() {
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile(_data_file_name);

        summarizeData(_data);
        createSpeciesMap(_data);

        //set number of species to number in data file
        setNumberTaxa(_data);
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
            
            setUpInitialData();
            
            unsigned nspecies = (unsigned) _species_names.size();
            Forest::setNumSpecies(nspecies);
            rng.setSeed(_random_seed);
            
            // create vector of particles
            unsigned nsubsets = _data->getNumSubsets();
            unsigned nparticles = _nparticles;
            
            Particle::setNumSubsets(nsubsets);
                            
            // set size of vectors (number of genes + 1 species tree)
            // my_vec contains a vector of vector of particles corresponding to species particles and gene particles
            vector<vector<Particle::SharedPtr>> my_vec_1(nsubsets+1, vector<Particle::SharedPtr> (nparticles)); // leave x% of vector empty for gene particles
            vector<vector<Particle::SharedPtr>> my_vec_2(nsubsets+1, vector<Particle::SharedPtr> (nparticles));
            vector<vector<Particle::SharedPtr>> &my_vec = my_vec_1;
            
            for (int s=0; s<nsubsets+1; s++) {
                _accepted_particle_vec.push_back(vector<Particle::SharedPtr>(nparticles));
            }
            
            _prev_log_marginal_likelihood = _log_marginal_likelihood;
            _log_marginal_likelihood = 0.0;
        
            for (unsigned s=0; s<nsubsets+1; s++) {
                int nparticles = _nparticles;
                for (unsigned i=0; i<nparticles; i++) {
                    my_vec_1[s][i] = Particle::SharedPtr(new Particle); // TODO: wasteful that this makes species trees with extra lineages
                    my_vec_2[s][i] = Particle::SharedPtr(new Particle);
                }
            }

            bool use_first = true;
            
            vector<string> newicks;
            if (_gene_newicks_names != "null" && _species_newicks_name != "null") {
                throw XProj(boost::str(boost::format("cannot specify gene newicks and species newicks; choose one")));

            }
            
            if (_gene_newicks_names != "null") {
                if (_niterations == 1) {
                    throw XProj(boost::str(boost::format("must specify more than 1 iteration if beginning from gene trees")));
                }
                ifstream infile(_gene_newicks_names);
                string newick;
                while (getline(infile, newick)) {
                    newicks.push_back(newick);
                }
            }
            else if (_species_newicks_name != "null") {
                ifstream infile(_species_newicks_name);
                string newick;
                while (getline(infile, newick)) {
                    newicks.push_back(newick);
                }
            }
            
            string start = "species"; // start variable defines if program should start with gene or species trees
            
            for (int s=0; s<nsubsets+1; s++) {
                if (s == 0) {
                    nparticles = _nparticles;
                }
                else {
                    nparticles = _nparticles;
                }
                for (int p=0; p<nparticles; p++) {
                    my_vec[s][p]->setData(_data, _taxon_map, s);
                    my_vec[s][p]->mapSpecies(_taxon_map, _species_names, s);
                    my_vec[s][p]->setParticleGeneration(0);
                    my_vec[s][p]->setLogLikelihood(0.0);
                    my_vec[s][p]->setLogCoalescentLikelihood(0.0);
                    my_vec[s][p]->setLogWeight(0.0, "g");
                    my_vec[s][p]->setLogWeight(0.0, "s");
                    
                    // only sample 1 species tree, and use this tree for all the gene filtering
                    if (s == 0 && _gene_newicks_names == "null" && p == 0) {
                        my_vec[0][p]->processSpeciesNewick(newicks); // if no newick specified, program will sample from species tree prior
                    }
                }
            }
            bool gene_first = false;
            if (_gene_newicks_names != "null") {
                assert (newicks.size() == nsubsets);
                for (int s=1; s<nsubsets+1; s++) {
                    for (int p=0; p<nparticles; p++) {
                        my_vec[s][p]->processGeneNewicks(newicks, s-1);
                        my_vec[s][p]->mapSpecies(_taxon_map, _species_names, s);
                        start = "gene";
                        gene_first = true;
                    }
                }
                _accepted_particle_vec = my_vec;
            }
            
            if (_run_on_empty) {
                for (int s=0; s<nsubsets+1; s++) {
                    for (int p=0; p<nparticles; p++) {
                        my_vec[s][p]->setRunOnEmpty(_run_on_empty);
                    }
                }
            }
            
            int ntaxa = (int) _taxon_map.size();
            bool deconstruct = false;
            
            for (int i=0; i<_niterations; i++) {
                _log_marginal_likelihood = 0.0;
                cout << "beginning iteration: " << i << endl;
                
                if (i > 0) {
                    for (int s=0; s<nsubsets+1; s++) {
                        for (int p=0; p<nparticles; p++) {
                            my_vec[s][p]->setParticleGeneration(0);
                            if (s > 0) {
                                my_vec[s][p]->mapGeneTrees(_taxon_map, _species_names);
                                my_vec[s][p]->resetGeneIncrements();
                            }
                        }
                    }
                }
                
                // keep the species partition for the gene forests at this stage but clear the tree structure
                 if (i == 1 && gene_first == true) {
                     for (int s=1; s<nsubsets+1; s++) {
                         // start at s=1 to only modify the gene trees
                         for (int p=0; p<nparticles; p++) {
                         my_vec[s][p]->remakeGeneTrees(_taxon_map);
                         my_vec[s][p]->resetGeneTreePartials(_data, _taxon_map, s);
                         deconstruct = false;
                         }
                     }
                 }
            
                if (i > 0) {
                    deconstruct = true;
                }
                
                // my_vec[0] is the species tree particles
                // my_vec[1] is gene 1 particles
                // my_vec[2] is gene 2 particles
                // etc
                
                // filter gene trees
                if (start == "species") {
                    // pick a species tree to use for all the gene trees for this step
                    
                    Particle::SharedPtr species_tree_particle;
                    
                    if (i > 0) {
                        normalizeWeights(my_vec[0], "s", false);
                        species_tree_particle = chooseTree(my_vec[0], "s"); // pass in all the species trees
            
                        // delete extra particles
                        int nparticles_to_remove = (nparticles*_species_particles_per_gene_particle) - nparticles;
                        for (int s=0; s<nsubsets+1; s++) {
                            my_vec[s].erase(my_vec[s].end() - nparticles_to_remove, my_vec[s].end());
                            my_vec_2[s].erase(my_vec_2[s].end() - nparticles_to_remove, my_vec_2[s].end());
                        }
                        
                    }
                    else {
                        species_tree_particle = my_vec[0][0];
                    }
                    
                    if (i > 0) {
                        // don't need to reset these variables for the first ietration
                        for (int s=1; s<nsubsets+1; s++) {
                            for (int p=0; p<nparticles; p++) {
                                my_vec[s][p]->setLogLikelihood(0.0);
                                my_vec[s][p]->setLogWeight(0.0, "g");
                            }
                        }
                    }
                    
                    species_tree_particle->showParticle();
                                            
                    for (unsigned g=0; g<ntaxa-1; g++) {
                        cout << "generation " << g << endl;
                        // filter particles within each gene
                        
                        bool gene_trees_only = true;
                        
                        for (int s=1; s<nsubsets+1; s++) { // skip species tree particles
                            
                            proposeParticles(my_vec[s], gene_trees_only, "g", deconstruct, species_tree_particle);
                            
                            if (!_run_on_empty) {
                                bool calc_marg_like = true;
                                
                                normalizeWeights(my_vec[s], "g", calc_marg_like);
                                
                                double ess_inverse = 0.0;
                                
                                for (int p=0; p<_nparticles; p++) {
                                    ess_inverse += exp(2.0*my_vec[s][p]->getLogWeight("g"));
                                }

//                                    double ess = 1.0/ess_inverse;
//                                    cout << "   " << "ESS = " << ess << endl;
                             
                                resampleParticles(my_vec[s], use_first ? my_vec_2[s]:my_vec_1[s], "g");
                                //if use_first is true, my_vec = my_vec_2
                                //if use_first is false, my_vec = my_vec_1
                                
                                my_vec[s] = use_first ? my_vec_2[s]:my_vec_1[s];
                                // do not need to resample species trees; species tree will remain the same throughout all gene tree filtering
                                
                                assert(my_vec[s].size() == nparticles);
//                                    assert(my_vec[s].size() == nparticles*_species_particles_per_gene_particle);
                            }
                            //change use_first from true to false or false to true
                            use_first = !use_first;
                            if (g < ntaxa-2) {
                                resetWeights(my_vec[s], "g");
                            }
                                assert (_accepted_particle_vec.size() == nsubsets+1);
                                _accepted_particle_vec[s] = my_vec[s];
//                                    saveParticleWeights(my_vec[0]);
                            } // s loop
                        deconstruct = false;
                    } // g loop
                }
                
                    writeGeneTreeFile();
                    
                    // filter species trees now
                    
                    // save gene tree variation for use in lorad file
                    vector<vector<Particle>> variable_gene_trees(nsubsets, vector<Particle> (nparticles));

                    for (int s=1; s<nsubsets+1; s++) {
                        for (int p=0; p<nparticles; p++) {
                            variable_gene_trees[s-1][p] = *my_vec[s][p];
                        }
                    }
                    
                // choose one set of gene trees to use
                for (int s=1; s<nsubsets+1; s++) {
                    normalizeWeights(my_vec[s], "g", false);
                    Particle gene_x = *chooseTree(my_vec[s], "g");
                    for (int p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
                        if (p<nparticles) {
                            *my_vec[s][p] = gene_x;
                        }
                        else {
                            // add in null particles here
                            my_vec[s].push_back(Particle::SharedPtr(new Particle));
                            *my_vec[s][p] = gene_x;
                            my_vec_2[s].push_back(Particle::SharedPtr(new Particle));
                        }
                    }
                }
                    
                // increase size of species vector
                for (int p=nparticles; p<nparticles*_species_particles_per_gene_particle; p++) {
                    my_vec[0].push_back(Particle::SharedPtr(new Particle));
                    my_vec_2[0].push_back(Particle::SharedPtr(new Particle));
                }
                
                for (int s=1; s<nsubsets+1; s++) {
                    for (int p=0; p<nparticles; p++) {
                        *_accepted_particle_vec[s][p] = variable_gene_trees[s-1][p]; // preserve gene tree variation for output
                    }
                }
                    
                if (i < _niterations-1) {
                    _species_tree_log_marginal_likelihood = 0.0;
                    
                for (int p=0; p<my_vec[0].size(); p++) {
                        my_vec[0][p]->mapSpecies(_taxon_map, _species_names, 0);
                        my_vec[0][p]->resetSpecies();
                    }
                    
                    for (int s=1; s<nsubsets+1; s++) {
                        for (int p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
                            my_vec[s][p]->mapSpecies(_taxon_map, _species_names, s);
                            my_vec[s][p]->refreshGeneTreePreorder();
                            my_vec[s][p]->calcGeneTreeMinDepth(); // reset min depth vector for gene trees
                            my_vec[s][p]->resetLogTopologyPrior();
                        }
                    }
                    
                    // filter species trees
                    for (int p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
                        my_vec[0][p]->setLogLikelihood(0.0);
                        my_vec[0][p]->setLogWeight(0.0, "g");
                        my_vec[0][p]->setLogWeight(0.0, "s");
                    }
                    
                    for (unsigned s=0; s<nspecies; s++) {
                        
                        cout << "beginning species tree proposals" << endl;
                        //taxon joining and reweighting step
                        
                        if (s == 0) {
                            for (int j=1; j<nsubsets+1; j++) {
                                // reset gene tree log coalescent likelihoods to 0
                                for (int p=0; p<nparticles*_species_particles_per_gene_particle; p++) {
                                    my_vec[j][p]->setLogCoalescentLikelihood(0.0);
                                }
                            }
                        }
                        
                        proposeSpeciesParticles(my_vec, s, nspecies, nsubsets);

                        // filter - make sure all gene trees go along with correct species tree
                        
                        if (!_run_on_empty) {
                            bool calc_marg_like = false;
                            
                            normalizeWeights(my_vec[0], "s", calc_marg_like);
                            
                            double ess_inverse = 0.0;
                            
                            for (auto & p:my_vec[0]) {
                                ess_inverse += exp(2.0*p->getLogWeight("s"));
                            }

                            double ess = 1.0/ess_inverse;
                            cout << "   " << "ESS = " << ess << endl;
                            vector<int> sel_indices = resampleSpeciesParticles(my_vec[0], use_first ? my_vec_2[0]:my_vec_1[0], "s");
                            //if use_first is true, my_vec = my_vec_2
                            //if use_first is false, my_vec = my_vec_1
                            
                            my_vec[0] = use_first ? my_vec_2[0]:my_vec_1[0];
                            
                            for (int s=1; s<nsubsets+1; s++) {
                                resetGeneParticles(sel_indices, my_vec[s], use_first ? my_vec_2[s]:my_vec_1[s]);
                                my_vec[s] = use_first ? my_vec_2[s]:my_vec_1[s];
                            }
                            
                            //change use_first from true to false or false to true
                            use_first = !use_first;
                        }
                        _accepted_particle_vec[0] = my_vec[0];
                        start = "species";
                    } // s loop
                    saveParticleWeights(my_vec[0]);
                }
            }
                                
            writeLoradFile(my_vec, nparticles, nsubsets, nspecies, ntaxa);
            
            writeGeneTreeFile();
            // saveParticleWeights(_accepted_particle_vec);
            // saveParticleLikelihoods(_accepted_particle_vec);

//            saveAllHybridNodes(_accepted_particle_vec);
            showFinal(_accepted_particle_vec);
            cout << "marg like: " << setprecision(12) << _log_marginal_likelihood << endl;
        }

        catch (XProj & x) {
            std::cerr << "Proj encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }

    inline void Proj::writeLoradFile(vector<vector<Particle::SharedPtr>> my_vec, int nparticles, int nsubsets, int nspecies, int ntaxa) {
        // open log file
        ofstream logf("params.log");
        
        double a = 0;
        unsigned col_count = 0;
        
        for (int p=0; p<nparticles; p++) {
            
            vector<double> branch_length_vec;
            for (int s=0; s<nsubsets+1; s++) {
                for (auto &b:my_vec[s][p]->getBranchLengths()) {
                    branch_length_vec.push_back(b);
                }
            }
            
            vector<double> prior_vec;
            for (int s=0; s<nsubsets+1; s++) {
                for (auto &b:my_vec[s][p]->getBranchLengthPriors()) {
                    prior_vec.push_back(b);
                }
            }
            
            vector<double> gene_tree_log_like;
            for (int s=1; s<nsubsets+1; s++) {
                for (auto &g:my_vec[s][p]->getGeneTreeLogLikelihoods()) {
                    gene_tree_log_like.push_back(g);
                }
            }
            
            vector<double> gene_tree_log_coalescent_like;
            for (int s=1; s<nsubsets+1; s++) {
                for (auto &g:my_vec[s][p]->getGeneTreeLogCoalescentLikelihood()) {
                    gene_tree_log_coalescent_like.push_back(g);
                }
            }
            
            vector<double> log_topology_priors;
            for (int s=0; s<nsubsets+1; s++) {
                for (auto &t:my_vec[s][p]->getTopologyPriors()) {
                    log_topology_priors.push_back(t);
                }
            }
            
            double species_tree_height = my_vec[0][p]->getSpeciesTreeHeight();
            
            assert(branch_length_vec.size() == prior_vec.size());
            
            double log_coalescent_likelihood = my_vec[0][p]->getCoalescentLikelihood();
            
            int ngenes = nsubsets;
            
            if (col_count == 0) {
                logf << "iter" << "\t" << "theta";
                for (int g=0; g<ngenes; g++) {
                    logf << "\t" << "gene_tree_log_like";
                }
                for (int g=0; g<ngenes; g++) {
                    logf << "\t" << "gene_tree_log_coalescent_like";
                }
                for (int i=0; i<nspecies-1; i++) {
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
                logf << "\t" << "log_coal_like";
                
                logf << "\t" << "species_tree_height" << endl;
            }
            
            logf << a << "\t" << Forest::_starting_theta;
            
            for (int g=0; g<gene_tree_log_like.size(); g++) {
                logf << "\t" << setprecision(12) << gene_tree_log_like[g];
            }
            
            for (int g=0; g<gene_tree_log_coalescent_like.size(); g++) {
                logf << "\t" << setprecision(12) << gene_tree_log_coalescent_like[g];
            }

            
            for (int i=0; i<prior_vec.size(); i++) {
                logf << "\t" << setprecision(11) << branch_length_vec[i] << "\t" << prior_vec[i];
            }
            
            for (int j=0; j<log_topology_priors.size(); j++) {
                logf << "\t" << setprecision(12) << log_topology_priors[j];
            }
            
            logf << "\t" << setprecision(12) << log_coalescent_likelihood;
            
            logf << "\t" << setprecision(12) << species_tree_height;
            
            logf << endl;
            a++;
            col_count++;
        }
        
        logf.close();
    }

}
