//# define USING_MPI
//#if defined(USING_MPI)
//#   include <mpi.h>
//#endif

//#define SNAKE

//#define NDEBUG // disable asserts for efficiency / memory improvement

#define USE_TOTAL_RATE

#define HIERARCHICAL_FILTERING

#define DRAW_NEW_THETA

//#define UNCONSTRAINED_PROPOSAL 

#define GRAHAM_JONES_COALESCENT_LIKELIHOOD // if defined, also define DRAW_NEW_THETA or turn off HIERARCHICAL_FILTERING
// do not define this if fix_theta defined in conf file

#define PARALLELIZE_BY_GROUP // do not turn this off if doing second level filtering - will be slower, and params file cannot be read by tracer

//#define DEBUG_MODE

//#define INV_GAMMA_PRIOR_TWO // turn this off if not defining DRAW_NEW_THETA

#define BUILD_UPGMA_TREE

//#define BUILD_UPGMA_TREE_CONSTRAINED // if this is turned on, also turn on BUILD_UPGMA_TREE

#define WEIGHT_MODIFIER // includes weight correction for differing theta prior mean and proposal mean

#define FASTER_UPGMA_TREE

//#define COAL_LIKE_TEST
//
// # define LIKELIHOOD_TEST // DO NOT USE
