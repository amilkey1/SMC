//# define USING_MPI
//#if defined(USING_MPI)
//#   include <mpi.h>
//#endif

//#define NDEBUG // disable asserts for efficiency / memory improvement

#define HIERARCHICAL_FILTERING

//#define UNCONSTRAINED_PROPOSAL 

//#define GRAHAM_JONES_COALESCENT_LIKELIHOOD // if defined, also define DRAW_NEW_THETA or turn off HIERARCHICAL_FILTERING

#define PARALLELIZE_BY_GROUP

//#define DEBUG_MODE

#define DRAW_NEW_THETA

#define INV_GAMMA_PRIOR_TWO // if DRAW_NEW_THETA turned off, INV_GAMMA_PRIOR_TWO needs to also be turned off

