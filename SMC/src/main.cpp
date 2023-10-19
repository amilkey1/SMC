#include "conditionals.hpp"

#if defined(USING_SIGNPOSTS)
#include <os/log.h>
#include <os/signpost.h>
os_log_t log_handle;
os_signpost_id_t signpost_id;
#endif

#include <iostream>
#include <vector>
#include "proj.hpp"
#include "particle.hpp"

// Initialize our random number generator here so it will be a global variable
#include "lot.hpp"
proj::Lot rng;

#include "partial_store.hpp"
proj::PartialStore ps;

using namespace proj;

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

using namespace std;

#include "forest.hpp"
#include "data.hpp"
std::string  Proj::_program_name        = "proj";
unsigned     Proj::_major_version       = 1;
unsigned     Proj::_minor_version       = 0;

unsigned Forest::_nspecies = 4;
unsigned Forest::_ntaxa = 12;
double Forest::_theta = 0.05;
double Forest::_lambda = 1;
unsigned Particle::_nsubsets = 1;
const double Node::_smallest_edge_length=1.0e-12;
string Forest::_proposal;
string Forest::_model;
double Forest::_kappa = 1.0;
vector<double> Forest::_base_frequencies;
string Forest::_string_base_frequencies;
double Forest::_migration_rate;
double Forest::_hybridization_rate;
string Forest::_outgroup;
bool Particle::_run_on_empty;
double Forest::_theta_prior_mean = 1.0;
double Forest::_lambda_prior_mean = 1000.0;

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
    
#if defined(USING_SIGNPOSTS)
log_handle  = os_log_create("edu.uconn.eeb.phylogeny", OS_LOG_CATEGORY_POINTS_OF_INTEREST);
signpost_id = os_signpost_id_generate(log_handle);
assert(signpost_id != OS_SIGNPOST_ID_INVALID);
#endif
    
#if defined(USING_MPI)
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
//    std::cerr << my_rank << "    " << ntasks << endl;
#endif
    
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
    
//    ofstream memfile("allocs.txt");
//    ps.memoryReport(memfile);
//    memfile.close();
    
#if defined(USING_MPI)
    MPI_Finalize();
#endif
    return 0;
}
