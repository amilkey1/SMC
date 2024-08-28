//# define USING_MPI
//#if defined(USING_MPI)
//#   include <mpi.h>
//#endif

//#define NDEBUG // disable asserts for efficiency / memory improvement

#define HIERARCHICAL_FILTERING

#define DRAW_NEW_THETA

//#define UNCONSTRAINED_PROPOSAL 

#define GRAHAM_JONES_COALESCENT_LIKELIHOOD // if defined, also define DRAW_NEW_THETA or turn off HIERARCHICAL_FILTERING

#define PARALLELIZE_BY_GROUP

//#define DEBUG_MODE

#define INV_GAMMA_PRIOR_TWO // turn this off if not defining DRAW_NEW_THETA

//#define DELAY_FILTERING

#define CONSISTENCY_TEST
