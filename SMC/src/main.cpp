#include <iostream>
#include <vector>
#include "proj.hpp"
#include "particle.hpp"
#include "conditionals.hpp"
#include "stopwatch.hpp"

// Initialize our random number generator here so it will be a global variable
#include "lot.hpp"
proj::Lot rng;

proj::StopWatch stopwatch;

#include "partial_store.hpp"
proj::PartialStore ps;

using namespace proj;
using namespace std;

#include "forest.hpp"
#include "data.hpp"
std::string  Proj::_program_name        = "proj";
unsigned     Proj::_major_version       = 1;
unsigned     Proj::_minor_version       = 0;
string       Proj::_start_mode          = "smc";

unsigned Forest::_nspecies = 4;
unsigned Forest::_ntaxa = 12;
double Forest::_theta = 0.05;
//double Forest::_lambda = 1;
unsigned Particle::_nsubsets = 1;
const double Node::_smallest_edge_length=1.0e-12;
string Forest::_proposal;
string Forest::_model;
string Forest::_outgroup;
double Forest::_kappa = 1.0;
vector<double> Forest::_base_frequencies;
string Forest::_string_base_frequencies;
bool Forest::_run_on_empty;
bool Forest::_run_on_empty_first_level_only;
bool Forest::_save_memory;
double Forest::_ploidy = 2.0;
double Forest::_theta_proposal_mean = 0.0;
double Forest::_theta_prior_mean = 0.0;
string Forest::_start_mode = "smc";
double Particle::_lambda_prior_mean = 0.0;
double Data::_occupancy = 1.0;
double Forest::_edge_rate_variance = 0.0;
double Forest::_asrv_shape = numeric_limits<double>::infinity();
double Forest::_infinity = numeric_limits<double>::infinity();
double Forest::_comphet = numeric_limits<double>::infinity();

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
    Proj proj;
    try {
        proj.processCommandLineOptions(argc, argv);
        proj.run();
    }
    catch(std::exception & x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    return 0;
}

