//# define USING_MPI
//#if defined(USING_MPI)
//#   include <mpi.h>
//#endif

//#define SNAKE

//#define NDEBUG // disable asserts for efficiency / memory improvement

#define USE_TOTAL_RATE

//#define SIM_TEST // simulated example with 2 taxa, 2 species - 1 taxon per species
//#define SIM_TEST3 // simulated example with 3 species, 3 taxa

#define HIERARCHICAL_FILTERING

#define DRAW_NEW_THETA

//#define UNCONSTRAINED_PROPOSAL 

#define GRAHAM_JONES_COALESCENT_LIKELIHOOD // if defined, also define DRAW_NEW_THETA or turn off HIERARCHICAL_FILTERING

#define PARALLELIZE_BY_GROUP

//#define DEBUG_MODE

#define INV_GAMMA_PRIOR_TWO // turn this off if not defining DRAW_NEW_THETA

#define TESTING_UNEVEN_LIKELIHOOD_CORRECTION
