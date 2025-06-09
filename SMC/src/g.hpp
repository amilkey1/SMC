#pragma once
#include "split.hpp"

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);
//extern proj::Lot::SharedPtr rng;

namespace proj {

    struct G {
        enum class ModelType {
            MODEL_TYPE_JC,
            MODEL_TYPE_HKY
        };
        
        enum class StartModeType {
            START_MODE_SIM,
            START_MODE_SMC
        };
        
        typedef unsigned long           species_t;
        static map<string, unsigned>    _taxon_to_species;
        
#if defined (LAZY_COPYING)
        typedef pair<unsigned, unsigned>                        uint_pair_t;
        typedef tuple<double,G::species_t,G::species_t>         merge_t; // <speciation height, spp 1, spp 2>
        typedef vector<merge_t>                                 merge_vect_t;
#endif
        
        // Program settings used in processCommandLineOptions
        static bool                 _use_gpu;
        static bool                 _ambig_missing;
        static string               _start_mode;
        static string               _program_name;
        static unsigned             _major_version;
        static unsigned             _minor_version;
        static string               _proposal;
        static bool                 _prior_prior;
        static string               _model;
        static string               _outgroup;
        static bool                 _run_on_empty;
        static bool                 _run_on_empty_first_level_only;
        static bool                 _save_memory;
        static unsigned             _nparticles;
        static vector<double>       _base_frequencies;
        static string               _string_base_frequencies;
        static string               _data_file_name;
        static unsigned             _nthreads;
        static unsigned             _verbose;
        static unsigned             _sim_nspecies;
        static string               _string_ntaxaperspecies;
        static string               _sim_file_name;
        static unsigned             _particle_increase;
        static double               _thin;
        static vector<unsigned>     _ntaxaperspecies;
        static unsigned             _save_every;
        static bool                 _save_gene_trees;
        static bool                 _gene_newicks_specified;
        static unsigned             _ngenes_provided;
        static string               _species_newick_name;
        static bool                 _species_newick_specified;
        static bool                 _fix_theta_for_simulations;
        static bool                 _fix_theta;
        static double               _theta;
        static double               _theta_proposal_mean;
        static double               _theta_prior_mean;
        static string               _string_relative_rates;
        static vector<double>       _double_relative_rates;
        static bool                 _save_gene_trees_separately;
        static string               _newick_path;
        static double               _lambda;
        static unsigned             _ngroups; // number of subpopulations for mixing SMC
        static bool                 _upgma;
        static unsigned             _nloci;
        static unsigned             _nspecies;
        static unsigned             _ntaxa;
        static double               _infinity;
        static unsigned             _nstates;
        static double               _ploidy;
        static double               _small_enough;
        static unsigned             _partial_count;
        static bool                 _mcmc;
        static double               _sliding_window;
        static unsigned             _nmcmc_moves_accepted;
        
        // other
        static vector<string>           _taxon_names;
        
        static vector<vector<double> > _dmatrix;
        static vector<Split>           _dmatrix_rows;
        
        static ModelType _model_type;
        static StartModeType _start_mode_type;
        
        static vector<string>       _species_names;
        
#if defined (LAZY_COPYING)
        static vector<G::species_t> _species_names_typed;
#endif
        
        static bool                 _in_second_level;
        static unsigned             _generation;
        
        // functions
        string inventName(unsigned k, bool lower_case);
        static unsigned multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs);
        
    };

    inline string G::inventName(unsigned k, bool lower_case) {
        // If   0 <= k < 26, returns A, B, ..., Z,
        // If  26 <= k < 702, returns AA, AB, ..., ZZ,
        // If 702 <= k < 18278, returns AAA, AAB, ..., ZZZ, and so on.
        //
        // For example, k = 19009 yields ABCD:
        // ABCD 19009 = 26 + 26*26 + 26*26*26 + 0*26*26*26 + 1*26*26 + 2*26 + 3
        //              <------- base ------>   ^first       ^second   ^third ^fourth
        // base = (26^4 - 1)/25 - 1 = 18278
        //   26^1 + 26^2 + 26^3 = 26^0 + 26^1 + 26^2 + 26^3 - 1 = (q^n - 1)/(q - 1) - 1, where q = 26, n = 4
        //   n = 1 + floor(log(19009)/log(26))
        // fourth = ((19009 - 18278                           )/26^0) % 26 = 3
        // third  = ((19009 - 18278 - 3*26^0                  )/26^1) % 26 = 2
        // second = ((19009 - 18278 - 3*26^0 - 2*26^1         )/26^2) % 26 = 1
        // first  = ((19009 - 18278 - 3*26^0 - 2*26^1 - 1*26^2)/26^3) % 26 = 0
                
        // Find how long a species name string must be
        double logibase26 = (k > 0 ? log(k)/log(26) : 0);
        unsigned n = 1 + (unsigned)floor(logibase26);
        vector<char> letters;
        unsigned base = (unsigned)((pow(26,n) - 1)/25.0 - 1);
        unsigned cum = 0;
        int ordA = (unsigned)(lower_case ? 'a' : 'A');
        for (unsigned i = 0; i < n; ++i) {
            unsigned ordi = (unsigned)((k - base - cum)/pow(26,i)) % 26;
            letters.push_back(char(ordA + ordi));
            cum += (unsigned)(ordi*pow(26,i));
        }
        string species_name(letters.rbegin(), letters.rend());
        return species_name;
    }

    inline unsigned G::multinomialDraw(Lot::SharedPtr lot, const vector<double> & probs) {
        // Compute cumulative probababilities
        vector<double> cumprobs(probs.size());
        partial_sum(probs.begin(), probs.end(), cumprobs.begin());
        assert(fabs(*(cumprobs.rbegin()) - 1.0) < 0.0001);

        // Draw a Uniform(0,1) random deviate
        double u = lot->uniform();

        // Find first element in cumprobs greater than u
        // e.g. probs = {0.2, 0.3, 0.4, 0.1}, u = 0.6, should return 2
        // because u falls in the third bin
        //
        //   |   0   |     1     |        2      | 3 | <-- bins
        //   |---+---+---+---+---+---+---+---+---+---|
        //   |       |           |   |           |   |
        //   0      0.2         0.5  |          0.9  1 <-- cumulative probabilities
        //                          0.6 <-- u
        //
        // cumprobs = {0.2, 0.5, 0.9, 1.0}, u = 0.6
        //               |         |
        //               begin()   it
        // returns 2 = 2 - 0
        auto it = find_if(cumprobs.begin(), cumprobs.end(), [u](double cumpr){return cumpr > u;});
        if (it == cumprobs.end()) {
            double last_cumprob = *(cumprobs.rbegin());
            throw XProj(format("G::multinomialDraw failed: u = %.9f, last cumprob = %.9f") % u % last_cumprob);
        }

        auto d = std::distance(cumprobs.begin(), it);
        assert(d >= 0);
        assert(d < probs.size());
        return (unsigned)d;
    }


}


