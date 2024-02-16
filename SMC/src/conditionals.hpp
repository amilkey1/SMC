//# define USING_MPI
//#if defined(USING_MPI)
//#   include <mpi.h>
//#endif

//#define SNAKE

//#define PRIOR_POST_ALL_PAIRS

//#define NDEBUG // disable asserts for efficiency / memory improvement

#define USE_TOTAL_RATE

//#define SAVE_UNIQUE_SPECIES_TREES

//#define SIM_TEST // simulated example with 2 taxa, 2 species - 1 taxon per species
//#define SIM_TEST3 // simulated example with 3 species, 3 taxa

#define WEIGHT_CORRECTION

#define EXTRA_SPECIES_SAMPLING // if defined, also define CONDITION_ON_ONE or HIERARCHICAL_FILTERING

// #define CONDITION_ON_ONE // do not uncomment both CONDITION_ON_ONE and HIERARCHICAL_FILTERING

#define HIERARCHICAL_FILTERING
