//# define USING_MPI
//#if defined(USING_MPI)
//#   include <mpi.h>
//#endif

//#define NDEBUG // disable asserts for efficiency / memory improvement

#define USE_TOTAL_RATE

#define HIERARCHICAL_FILTERING

#define DRAW_NEW_THETA // TODO: FOR NOW, TURN THIS OFF IF MPI TURNED ON

//#define DEBUG_MODE

//#define INV_GAMMA_PRIOR_TWO // turn this off if not defining DRAW_NEW_THETA // TODO: this will not work with lazy copying

#define WEIGHT_MODIFIER // includes weight correction for differing theta prior mean and proposal mean

// # define LIKELIHOOD_TEST // DO NOT USE

//#define DEBUG_COALLIKE

#define SYSTEMATIC_FILTERING

#define REUSE_PARTIALS

#define UNROLL_LOOPS

#define LAZY_COPYING // TODO: only works with JC model

//#define UPGMA // turn on UPGMA in conf file but turn this off to avoid creating any members associated with UPGMA

//#define SPECIES_IN_CONF
