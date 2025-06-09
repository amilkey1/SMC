#include <iostream>
#include <vector>

#include "proj.hpp"
#include "conditionals.hpp"
#include "stopwatch.hpp"

#include "g.hpp"

// Initialize our random number generator here so it will be a global variable
#include "lot.hpp"
proj::Lot rng;

proj::StopWatch stopwatch;

#include "partial_store.hpp"
proj::PartialStore ps;

using namespace proj;
using namespace std;

// Include this header to enable macros that Valgrind recognizes
// but which are ignored if Valgrind is not being used
#include "valgrind.h"

int my_rank = 0;
#if defined(USING_MPI)
//int my_rank = 0;
int ntasks = 0;

void output(string msg) {
    if (my_rank == 0) {
        cout << msg;
    }
}
#else
void output(string msg) {
    cout << msg;
}
#endif

#include "data.hpp"
#include "forest.hpp"
#if defined(LAZY_COPYING)
#   include "forest-extension.hpp"
#endif
#include "species-forest.hpp"
#include "particle.hpp"

string      G::_program_name        = "smc";
unsigned    G::_major_version       = 1;
unsigned    G::_minor_version       = 0;
string      G::_start_mode          = "smc";
string      G::_proposal = "prior-prior";
bool        G::_prior_prior = true;
string      G::_model = "JC";
string      G::_outgroup = "none";
bool        G::_run_on_empty = false;
bool        G::_run_on_empty_first_level_only = false;
bool        G::_save_memory = false;
unsigned    G::_nparticles = 100;
unsigned    G::_nthreads = 1;
bool        G::_use_gpu = true;
bool        G::_ambig_missing = true;
unsigned    G::_verbose = 1;
unsigned    G::_sim_nspecies = 4;
string      G::_string_ntaxaperspecies = "";
string      G::_sim_file_name = "sim.nex";
unsigned    G::_particle_increase = 100;
double      G::_thin = 1.0;
vector<unsigned>    G::_ntaxaperspecies;
unsigned    G::_save_every = 1.0;
bool        G::_save_gene_trees = true;
bool        G::_gene_newicks_specified = false;
unsigned    G::_ngenes_provided = 2;
string      G::_species_newick_name = "";
bool        G::_species_newick_specified = false;
bool        G::_fix_theta_for_simulations = true;
bool        G::_fix_theta = true;
double      G::_theta = 0.1;
double      G::_theta_proposal_mean = 0.1;
double      G::_theta_prior_mean = 0.1;
string      G::_string_relative_rates = "";
vector<double>  G::_double_relative_rates;
string      G::_data_file_name = "sim.nex";
bool        G::_save_gene_trees_separately = false;
string      G::_newick_path = ".";
double      G::_lambda = 10.0;
unsigned    G::_ngroups = 1;
bool        G::_upgma = true;
double G::_infinity = numeric_limits<double>::infinity();
vector<double> G::_base_frequencies;
string      G::_string_base_frequencies;
vector<string>           G::_taxon_names;
unsigned    G::_nloci = 0;
unsigned    G::_nspecies = 0;
unsigned    G::_ntaxa = 0;
bool        G::_mcmc = false;
double      G::_sliding_window = 0.05;
unsigned    G::_nmcmc_moves_accepted = 0;

map<string, unsigned>    G::_taxon_to_species;
unsigned    G::_nstates = 4;
double      G::_ploidy = 2.0;
double      G::_small_enough = 0.0000001;

vector<vector<double> > G::_dmatrix;
vector<Split>           G::_dmatrix_rows;

G::ModelType G::_model_type;
G::StartModeType G::_start_mode_type;

const double Node::_smallest_edge_length=1.0e-12;
double Forest::_kappa = 1.0;
double Data::_occupancy = 1.0;
double Forest::_edge_rate_variance = 0.0;
double Forest::_asrv_shape = numeric_limits<double>::infinity();
double Forest::_comphet = numeric_limits<double>::infinity();

bool G::_in_second_level = false;
unsigned G::_generation = 0;

vector<string> G::_species_names;
unsigned G::_partial_count = 0;

#if defined (LAZY_COPYING)
vector<G::species_t> G::_species_names_typed;
#endif

GeneticCode::genetic_code_definitions_t GeneticCode::_definitions = {
                             // codon order is alphabetical: i.e. AAA, AAC, AAG, AAT, ACA, ..., TTT
    {"standard",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"vertmito",             "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"yeastmito",            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"moldmito",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"invertmito",           "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"ciliate",              "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"},
    {"echinomito",           "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"euplotid",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"},
    {"plantplastid",         "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"altyeast",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"ascidianmito",         "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"altflatwormmito",      "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"},
    {"blepharismamacro",     "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"},
    {"chlorophyceanmito",    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"},
    {"trematodemito",        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"scenedesmusmito",      "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"},
    {"thraustochytriummito", "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"}
};

int main(int argc, const char * argv[]) {
#if defined(USING_MPI)
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
//    std::cerr << my_rank << "    " << ntasks << endl;
#endif
    
    Proj proj;
    try {
        StopWatch sw;
        sw.start();
        proj.processCommandLineOptions(argc, argv);
        proj.run();
        double total_seconds = sw.stop();
        cout << "total time: " << total_seconds << endl;
    }
    catch(std::exception & x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    
#if defined(USING_MPI)
    MPI_Finalize();
#endif
    return 0;
}

