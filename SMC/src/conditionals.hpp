//# define USING_MPI
//#if defined(USING_MPI)
//#   include <mpi.h>
//#endif

//#define NDEBUG // disable asserts for efficiency / memory improvement

#define USE_TOTAL_RATE

#define HIERARCHICAL_FILTERING

#define DRAW_NEW_THETA // TODO: FOR NOW, TURN THIS OFF IF MPI TURNED ON

//#define UNCONSTRAINED_PROPOSAL 

#define GRAHAM_JONES_COALESCENT_LIKELIHOOD // if defined, also define DRAW_NEW_THETA or turn off HIERARCHICAL_FILTERING
// do not define this if fix_theta defined in conf file

#define PARALLELIZE_BY_GROUP // do not turn this off if doing second level filtering - will be slower, and params file cannot be read by tracer

//#define DEBUG_MODE

//#define INV_GAMMA_PRIOR_TWO // turn this off if not defining DRAW_NEW_THETA

#define WEIGHT_MODIFIER // includes weight correction for differing theta prior mean and proposal mean

//#define COAL_LIKE_TEST // DO NOT USE
//
// # define LIKELIHOOD_TEST // DO NOT USE

#define FASTER_SECOND_LEVEL

//#define DEBUG_COALLIKE

//#define OLD_UPGMA

//#define UNUSED_FUNCTIONS

#define SYSTEMATIC_FILTERING

#define REUSE_PARTIALS
