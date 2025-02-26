#pragma once

extern void output(string msg, unsigned level);
extern void output(format & fmt, unsigned level);
//extern proj::Lot::SharedPtr rng;

namespace proj {

    struct G {
        // Program settings used in processCommandLineOptions
        static bool                 _use_gpu;
        static bool                 _ambig_missing;
        static string               _start_mode;
        static string               _program_name;
        static unsigned             _major_version;
        static unsigned             _minor_version;
        static string               _proposal;
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
        static unsigned             _ngroups;
        static bool                 _upgma;
        
        // functions
        string inventName(unsigned k, bool lower_case);
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

}


