/*
 *         compute the duplex structure of two RNA strands,
 *              allowing only inter-strand base pairs.
 *       see cofold() for computing hybrid structures without
 *                           restriction.
 *
 *                           Ivo Hofacker
 *                        Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
// #include "utils/basic.h"
// #include "params/default.h"
// #include "fold_vars.h"
// #include "fold.h"
// #include "pair_mat.h"
// #include "params/basic.h"
// #include "alifold.h"
// #include "subopt.h"
// #include "loops/all.h"
// #include "loops/internal.h"
// #include "duplex.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define STACK_BULGE1  1     /* stacking energies for bulges of size 1 */
#define NEW_NINIO     1     /* new asymetry penalty */
#define MAXSECTORS    500   /* dimension for a backtrack array */
#define LOCALITY      0.    /* locality parameter for base-pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT
#define NONE -10000         /* score for forbidden pairs */
#define PUBLIC
#define PRIVATE static
#define MAXALPHA 20       /* maximal length of alphabet */

static short  alias[MAXALPHA + 1];
static int    pair[MAXALPHA + 1][MAXALPHA + 1];
/* rtype[pair[i][j]]:=pair[j][i] */
static int    rtype[8] = {
  0, 2, 1, 4, 3, 6, 5, 7
};
#define NBASES 8
static const char Law_and_Order[]         = "_ACGUTXKI";
static int        BP_pair[NBASES][NBASES] =
  /* _  A  C  G  U  X  K  I */
{ { 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 5, 0, 0, 5 },
  { 0, 0, 0, 1, 0, 0, 0, 0 },
  { 0, 0, 2, 0, 3, 0, 0, 0 },
  { 0, 6, 0, 4, 0, 0, 0, 6 },
  { 0, 0, 0, 0, 0, 0, 2, 0 },
  { 0, 0, 0, 0, 0, 1, 0, 0 },
  { 0, 6, 0, 0, 5, 0, 0, 0 } };

/** The gas constant */
#define GASCONST 1.98717  /* in [cal/K] */
/** 0 deg Celsius in Kelvin */
#define K0  273.15
/** Infinity as used in minimization routines */
#define INF 10000000 /* (INT_MAX/10) */

#define EMAX (INF/10)
/** forbidden */
#define FORBIDDEN 9999
/** bonus contribution */
#define BONUS 10000
/** The number of distinguishable base pairs */
#define NBPAIRS 7
/** The minimum loop length */
#define TURN 3
/** The maximum loop length */
#define MAXLOOP 30

#define   VRNA_GQUAD_MAX_STACK_SIZE     7
#define   VRNA_GQUAD_MIN_STACK_SIZE     2
#define   VRNA_GQUAD_MAX_LINKER_LENGTH  15
#define   VRNA_GQUAD_MIN_LINKER_LENGTH  1
#define   VRNA_GQUAD_MIN_BOX_SIZE       ((4 * VRNA_GQUAD_MIN_STACK_SIZE) + \
                                         (3 * VRNA_GQUAD_MIN_LINKER_LENGTH))
#define   VRNA_GQUAD_MAX_BOX_SIZE       ((4 * VRNA_GQUAD_MAX_STACK_SIZE) + \
                                         (3 * VRNA_GQUAD_MAX_LINKER_LENGTH))

#define VRNA_MODEL_DEFAULT_TEMPERATURE    37.0

#define VRNA_MODEL_DEFAULT_TEMPERATURE    37.0

/**
 *  @brief  Default scaling factor for partition function computations
 *
 *  @see  #vrna_exp_param_t.pf_scale, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_PF_SCALE       -1

/**
 *  @brief  Default scaling factor for absolute thermodynamic temperature in Boltzmann factors
 *
 *  @see    #vrna_exp_param_t.alpha, #vrna_md_t.betaScale, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BETA_SCALE     1.

/**
 *  @brief  Default dangling end model
 *
 *  @see  #vrna_md_t.dangles, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_DANGLES        2

/**
 *  @brief  Default model behavior for lookup of special tri-, tetra-, and hexa-loops
 *
 *  @see    #vrna_md_t.special_hp, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_SPECIAL_HP     1

/**
 *  @brief  Default model behavior for so-called 'lonely pairs'
 *
 *  @see    #vrna_md_t.noLP, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_LP          0

/**
 *  @brief  Default model behavior for G-U base pairs
 *
 *  @see    #vrna_md_t.noGU, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU          0

/**
 *  @brief  Default model behavior for G-U base pairs closing a loop
 *
 *  @see    #vrna_md_t.noGUclosure, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU_CLOSURE  0

/**
 *  @brief  Default model behavior to treat a molecule as a circular RNA (DNA)
 *  @see    #vrna_md_t.circ, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_CIRC           0

/**
 *  @brief  Default model behavior regarding the treatment of G-Quadruplexes
 *
 *  @see    #vrna_md_t.gquad, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_GQUAD          0

/**
 *  @brief  Default behavior of the model regarding unique multi-branch loop decomposition
 *
 *  @see    #vrna_md_t.uniq_ML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_UNIQ_ML        0

/**
 *  @brief  Default model behavior on which energy set to use
 *
 *  @see    #vrna_md_t.energy_set, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ENERGY_SET     0

/**
 *  @brief  Default model behavior with regards to backtracking of structures
 *  @see    #vrna_md_t.backtrack, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BACKTRACK      1

/**
 *  @brief  Default model behavior on what type of backtracking to perform
 *
 *  @see    #vrna_md_t.backtrack_type, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BACKTRACK_TYPE 'F'

/**
 *  @brief  Default model behavior with regards to computing base pair probabilities
 *
 *  @see    #vrna_md_t.compute_bpp, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_COMPUTE_BPP    1

/**
 *  @brief  Default model behavior for the allowed maximum base pair span
 *
 *  @see    #vrna_md_t.max_bp_span, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_MAX_BP_SPAN    -1

/**
 *  @brief  Default model behavior for the sliding window approach
 *
 *  @see    #vrna_md_t.window_size, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_WINDOW_SIZE    -1

/**
 *  @brief  Default model behavior on how to evaluate the energy contribution of multi-branch loops
 *
 *  @see    #vrna_md_t.logML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_LOG_ML         0

/**
 *  @brief  Default model behavior for consensus structure energy evaluation
 *  @see    #vrna_md_t.oldAliEn, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_OLD_EN     0

/**
 *  @brief  Default model behavior for consensus structure co-variance contribution assessment
 *
 *  @see    #vrna_md_t.ribo, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_RIBO       0

/**
 *  @brief  Default model behavior for weighting the co-variance score in consensus structure prediction
 *
 *  @see    #vrna_md_t.cv_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_CV_FACT    1.

/** @brief  Default model behavior for weighting the nucleotide conservation? in consensus structure prediction
 *
 *  @see    #vrna_md_t.nc_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_NC_FACT    1.


#define VRNA_MODEL_DEFAULT_PF_SMOOTH      1

/**
 *  @brief  Default model salt concentration (M)
 */
#define VRNA_MODEL_DEFAULT_SALT           1.021


/**
 *  @brief  Default model lower bound of multiloop size for salt correction fiting
 */
#define VRNA_MODEL_DEFAULT_SALT_MLLOWER    6


/**
 *  @brief  Default model upper bound of multiloop size for salt correction fiting
 */
#define VRNA_MODEL_DEFAULT_SALT_MLUPPER    24


/**
 *  @brief  Default model value to turn off user-provided salt correction for duplex initializtion
 */
#define VRNA_MODEL_DEFAULT_SALT_DPXINIT       99999

#define VRNA_MODEL_SALT_DPXINIT_FACT_RNA      -45.324
#define VRNA_MODEL_SALT_DPXINIT_FACT_DNA      -58.389


#define VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT  VRNA_MODEL_SALT_DPXINIT_FACT_RNA

/* Geometric parameters for RNA and DNA */

#define VRNA_MODEL_HELICAL_RISE_RNA   2.8
#define VRNA_MODEL_HELICAL_RISE_DNA   3.4
/**
 *  @brief  Default helical rise
 */
#define VRNA_MODEL_DEFAULT_HELICAL_RISE   VRNA_MODEL_HELICAL_RISE_RNA

#define VRNA_MODEL_BACKBONE_LENGTH_RNA   6.0
#define VRNA_MODEL_BACKBONE_LENGTH_DNA   6.76
/**
 *  @brief  Default backbone length
 */
#define VRNA_MODEL_DEFAULT_BACKBONE_LENGTH   VRNA_MODEL_BACKBONE_LENGTH_RNA


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#ifndef MAXALPHA
/**
 *  @brief Maximal length of alphabet
 */
#define MAXALPHA              20

#endif

#endif

#define   BP_REV_DEFAULT        { 0, 2, 1, 4, 3, 6, 5, 7 }

#define   BP_ALIAS_DEFAULT      { 0, 1, 2, 3, 4, 3, 2, 0 }

#define   BP_ENCODING_DEFAULT \
  /*  _  A  C  G  U  X  K  I */ \
  { { 0, 0, 0, 0, 0, 0, 0, 0 }, \
    { 0, 0, 0, 0, 5, 0, 0, 5 }, \
    { 0, 0, 0, 1, 0, 0, 0, 0 }, \
    { 0, 0, 2, 0, 3, 0, 0, 0 }, \
    { 0, 6, 0, 4, 0, 0, 0, 6 }, \
    { 0, 0, 0, 0, 0, 0, 2, 0 }, \
    { 0, 0, 0, 0, 0, 1, 0, 0 }, \
    { 0, 6, 0, 0, 5, 0, 0, 0 } }

#define   DM_DEFAULT \
  { { 0, 0, 0, 0, 0, 0, 0 }, /* hamming distance between pairs */ \
    { 0, 0, 2, 2, 1, 2, 2 } /* CG */, \
    { 0, 2, 0, 1, 2, 2, 2 } /* GC */, \
    { 0, 2, 1, 0, 2, 1, 2 } /* GU */, \
    { 0, 1, 2, 2, 0, 2, 1 } /* UG */, \
    { 0, 2, 2, 1, 2, 0, 2 } /* AU */, \
    { 0, 2, 2, 2, 1, 2, 0 } /* UA */ }

typedef double FLT_OR_DBL;


double          temperature     = VRNA_MODEL_DEFAULT_TEMPERATURE;
double          pf_scale        = VRNA_MODEL_DEFAULT_PF_SCALE;
int             dangles         = VRNA_MODEL_DEFAULT_DANGLES;
int             tetra_loop      = VRNA_MODEL_DEFAULT_SPECIAL_HP;
int             noLonelyPairs   = VRNA_MODEL_DEFAULT_NO_LP;
int             noGU            = 1;
int             no_closingGU    = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
int             circ            = VRNA_MODEL_DEFAULT_CIRC;
int             gquad           = VRNA_MODEL_DEFAULT_GQUAD;
int             uniq_ML         = VRNA_MODEL_DEFAULT_UNIQ_ML;
int             energy_set      = VRNA_MODEL_DEFAULT_ENERGY_SET;
int             do_backtrack    = VRNA_MODEL_DEFAULT_COMPUTE_BPP;
char            backtrack_type  = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE;
char            *nonstandards   = NULL;
int             max_bp_span     = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
int             oldAliEn        = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
int             ribo            = VRNA_MODEL_DEFAULT_ALI_RIBO;
double          cv_fact         = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
double          nc_fact         = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
int             logML           = VRNA_MODEL_DEFAULT_LOG_ML;

/* below are some more deprecated global symbols we need to get rid off */

int             james_rule        = 1;    /* interior loops of size 2 get energy 0.8Kcal and
                                           * no mismatches (no longer used) */
char            *RibosumFile      = NULL; /* TODO: compile ribosums into program
                                           * Warning: this variable will vanish */
int             csv               = 0;    /*generate comma seperated output*/
// vrna_bp_stack_t *base_pair        = NULL;
FLT_OR_DBL      *pr               = NULL; /* base pairing prob. matrix */
int             *iindx            = NULL; /* pr[i,j] -> pr[iindx[i]-j] */
int             fold_constrained  = 0;    /* fold with constraints */

double          salt             = VRNA_MODEL_DEFAULT_SALT;
int             saltDPXInit      = VRNA_MODEL_DEFAULT_SALT_DPXINIT;
float           saltDPXInitFact  = VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT;
float           helical_rise     = VRNA_MODEL_DEFAULT_HELICAL_RISE;
float           backbone_length  = VRNA_MODEL_DEFAULT_BACKBONE_LENGTH;

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

struct vrna_md_s {
  double  temperature;                      /**<  @brief  The temperature used to scale the thermodynamic parameters */
  double  betaScale;                        /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
  int     pf_smooth;                        /**<  @brief  A flat specifying whether energies in Boltzmann factors need to be smoothed */
  int     dangles;                          /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
                                             *
                                             *    If set to 0 no stabilizing energies are assigned to bases adjacent to
                                             *    helices in free ends and multiloops (so called dangling ends). Normally
                                             *    (dangles = 1) dangling end energies are assigned only to unpaired
                                             *    bases and a base cannot participate simultaneously in two dangling ends. In
                                             *    the partition function algorithm vrna_pf() these checks are neglected.
                                             *    To provide comparability between free energy minimization and partition function
                                             *    algorithms, the default setting is 2.
                                             *    This treatment of dangling ends gives more favorable energies to helices
                                             *    directly adjacent to one another, which can be beneficial since such
                                             *    helices often do engage in stabilizing interactions through co-axial
                                             *    stacking.\n
                                             *    If set to 3 co-axial stacking is explicitly included for
                                             *    adjacent helices in multiloops. The option affects only mfe folding
                                             *    and energy evaluation (vrna_mfe() and vrna_eval_structure()), as
                                             *    well as suboptimal folding (vrna_subopt()) via re-evaluation of energies.
                                             *    Co-axial stacking with one intervening mismatch is not considered so far.
                                             *    Note, that some function do not implement all dangle model but only a subset of
                                             *    (0,1,2,3). In particular, partition function algorithms can only handle
                                             *    0 and 2. Read the documentation of the particular recurrences or
                                             *    energy evaluation function for information about the provided dangle
                                             *    model.
                                             */
  int     special_hp;                       /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
  int     noLP;                             /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
  int     noGU;                             /**<  @brief  Do not allow GU pairs */
  int     noGUclosure;                      /**<  @brief  Do not allow loops to be closed by GU pair */
  int     logML;                            /**<  @brief  Use logarithmic scaling for multiloops */
  int     circ;                             /**<  @brief  Assume RNA to be circular instead of linear */
  int     gquad;                            /**<  @brief  Include G-quadruplexes in structure prediction */
  int     uniq_ML;                          /**<  @brief  Flag to ensure unique multi-branch loop decomposition during folding */
  int     energy_set;                       /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
  int     backtrack;                        /**<  @brief  Specifies whether or not secondary structures should be backtraced */
  char    backtrack_type;                   /**<  @brief  Specifies in which matrix to backtrack */
  int     compute_bpp;                      /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
  char    nonstandards[64];                 /**<  @brief  contains allowed non standard bases */
  int     max_bp_span;                      /**<  @brief  maximum allowed base pair span */

  int     min_loop_size;                    /**<  @brief  Minimum size of hairpin loops
                                             *
                                             *    The default value for this field is #TURN, however, it may
                                             *    be 0 in cofolding context.
                                             */
  int     window_size;                      /**<  @brief  Size of the sliding window for locally optimal structure prediction */
  int     oldAliEn;                         /**<  @brief  Use old alifold energy model */
  int     ribo;                             /**<  @brief  Use ribosum scoring table in alifold energy model */
  double  cv_fact;                          /**<  @brief  Co-variance scaling factor for consensus structure prediction */
  double  nc_fact;                          /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
  double  sfact;                            /**<  @brief  Scaling factor for partition function scaling */
  int     rtype[8];                         /**<  @brief  Reverse base pair type array */
  short   alias[MAXALPHA + 1];              /**<  @brief  alias of an integer nucleotide representation */
  int     pair[MAXALPHA + 1][MAXALPHA + 1]; /**<  @brief  Integer representation of a base pair */
  float   pair_dist[7][7];                  /**<  @brief  Base pair dissimilarity, a.k.a. distance matrix */
  double  salt;                             /**<  @brief  Salt (monovalent) concentration (M) in buffer */
  int     saltMLLower;                      /**<  @brief  Lower bound of multiloop size to use in loop salt correction linear fitting */
  int     saltMLUpper;                      /**<  @brief  Upper bound of multiloop size to use in loop salt correction linear fitting */
  int     saltDPXInit;                      /**<  @brief  User-provided salt correction for duplex initialization (in dcal/mol).
                                             *    If set to 99999 the default salt correction is used.
                                             *    If set to 0 there is no salt correction for duplex initialization.
                                             */
  float   saltDPXInitFact;                  /**<  @brief  */
  float   helical_rise;                     /**<  @brief  */
  float   backbone_length;                  /**<  @brief  */
};



typedef struct vrna_md_s vrna_md_t;


struct vrna_param_s {
  int       id;
  int       stack[NBPAIRS + 1][NBPAIRS + 1];
  int       hairpin[31];
  int       bulge[MAXLOOP + 1];
  int       internal_loop[MAXLOOP + 1];
  int       mismatchExt[NBPAIRS + 1][5][5];
  int       mismatchI[NBPAIRS + 1][5][5];
  int       mismatch1nI[NBPAIRS + 1][5][5];
  int       mismatch23I[NBPAIRS + 1][5][5];
  int       mismatchH[NBPAIRS + 1][5][5];
  int       mismatchM[NBPAIRS + 1][5][5];
  int       dangle5[NBPAIRS + 1][5];
  int       dangle3[NBPAIRS + 1][5];
  int       int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  int       int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  int       int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  int       ninio[5];
  double    lxc;
  int       MLbase;
  int       MLintern[NBPAIRS + 1];
  int       MLclosing;
  int       TerminalAU;
  int       DuplexInit;
  int       Tetraloop_E[200];
  char      Tetraloops[1401];
  int       Triloop_E[40];
  char      Triloops[241];
  int       Hexaloop_E[40];
  char      Hexaloops[1801];
  int       TripleC;
  int       MultipleCA;
  int       MultipleCB;
  int       gquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  int       gquadLayerMismatch;
  int       gquadLayerMismatchMax;

  double    temperature;      /**<  @brief  Temperature used for loop contribution scaling */

  vrna_md_t model_details;    /**<  @brief  Model details to be used in the recursions */
  char      param_file[256];  /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */
  int       SaltStack;
  int       SaltLoop[MAXLOOP + 2];
  double    SaltLoopDbl[MAXLOOP + 2];
  int       SaltMLbase;
  int       SaltMLintern;
  int       SaltMLclosing;
  int       SaltDPXInit;
};


typedef struct  vrna_param_s vrna_param_t;
PRIVATE vrna_md_t defaults = {
  VRNA_MODEL_DEFAULT_TEMPERATURE,
  1.,
  VRNA_MODEL_DEFAULT_PF_SMOOTH,
  VRNA_MODEL_DEFAULT_DANGLES,
  VRNA_MODEL_DEFAULT_SPECIAL_HP,
  VRNA_MODEL_DEFAULT_NO_LP,
  VRNA_MODEL_DEFAULT_NO_GU,
  VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,
  VRNA_MODEL_DEFAULT_LOG_ML,
  VRNA_MODEL_DEFAULT_CIRC,
  VRNA_MODEL_DEFAULT_GQUAD,
  VRNA_MODEL_DEFAULT_UNIQ_ML,
  VRNA_MODEL_DEFAULT_ENERGY_SET,
  VRNA_MODEL_DEFAULT_BACKTRACK,
  VRNA_MODEL_DEFAULT_BACKTRACK_TYPE,
  VRNA_MODEL_DEFAULT_COMPUTE_BPP,
  { 0 },
  VRNA_MODEL_DEFAULT_MAX_BP_SPAN,
  TURN,
  VRNA_MODEL_DEFAULT_WINDOW_SIZE,
  VRNA_MODEL_DEFAULT_ALI_OLD_EN,
  VRNA_MODEL_DEFAULT_ALI_RIBO,
  VRNA_MODEL_DEFAULT_ALI_CV_FACT,
  VRNA_MODEL_DEFAULT_ALI_NC_FACT,
  1.07,
  BP_REV_DEFAULT,
  BP_ALIAS_DEFAULT,
  BP_ENCODING_DEFAULT,
  DM_DEFAULT,
  VRNA_MODEL_DEFAULT_SALT,
  VRNA_MODEL_DEFAULT_SALT_MLLOWER,
  VRNA_MODEL_DEFAULT_SALT_MLUPPER,
  VRNA_MODEL_DEFAULT_SALT_DPXINIT,
  VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT,
  VRNA_MODEL_DEFAULT_HELICAL_RISE,
  VRNA_MODEL_DEFAULT_BACKBONE_LENGTH
};

PRIVATE vrna_param_t  *P = NULL;
PRIVATE int           **c = NULL;     /* energy array, given that i-j pair */
PRIVATE short         *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;
PRIVATE int           n1, n2;         /* sequence lengths */

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
 */
#pragma omp threadprivate(P, c, S1, SS1, S2, SS2, n1, n2)

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */





PRIVATE char *
backtrack(int i,
          int j);



/*
structure definitions


*/




typedef struct {
  int     i;
  int     j;
  int     end;
  char    *structure;
  double  energy;
  double  energy_backtrack;
  double  opening_backtrack_x;
  double  opening_backtrack_y;
  int     offset;
  double  dG1;
  double  dG2;
  double  ddG;
  int     tb;
  int     te;
  int     qb;
  int     qe;
} duplexT;

PRIVATE duplexT
duplexfold_cu(const char  *s1,
              const char  *s2,
              int         clean_up);








/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_nucleotide_encode(char       c,
                       vrna_md_t  *md)
{
  /* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */
  int code = -1;

  c = toupper(c);

  if (md) {
    if (md->energy_set > 0) {
      code = (int)(c - 'A') + 1;
    } else {
      const char *pos;
      pos = strchr(Law_and_Order, c);
      if (pos == NULL)
        code = 0;
      else
        code = (int)(pos - Law_and_Order);

      if (code > 5)
        code = 0;

      if (code > 4)
        code--;           /* make T and U equivalent */
    }
  }
  return code;
}


PUBLIC char
vrna_nucleotide_decode(int        enc,
                       vrna_md_t  *md)
{
  if (md) {
    if (md->energy_set > 0)
      return (char)enc + 'A' - 1;
    else
      return (char)Law_and_Order[enc];
  } else {
    return (char)0;
  }
}

PRIVATE const float
dm_default[7][7] = DM_DEFAULT;

PRIVATE void prepare_default_pairs(vrna_md_t *md)
{
  unsigned int i, j;

  for (i = 0; i < 5; i++)
    md->alias[i] = (short)i;

  md->alias[5]  = 3;          /* X <-> G */
  md->alias[6]  = 2;          /* K <-> C */
  md->alias[7]  = 0;          /* I <-> default base '@' */

  for (i = 0; i < NBASES; i++)
    for (j = 0; j < NBASES; j++)
      md->pair[i][j] = BP_pair[i][j];

  if (md->noGU)
    md->pair[3][4] = md->pair[4][3] = 0;


  if (md->nonstandards[0] != '\0') {
    /* allow nonstandard bp's (encoded by type=7) */
    for (i = 0; i < strlen(md->nonstandards); i += 2)
      md->pair[vrna_nucleotide_encode(md->nonstandards[i], md)]
      [vrna_nucleotide_encode(md->nonstandards[i + 1], md)] = 7;
  }
}


PRIVATE void
fill_pair_matrices(vrna_md_t *md)
{
  int i, j;

  /* nullify everything */
  for (i = 0; i <= MAXALPHA; i++)
    memset(md->pair[i], 0, (MAXALPHA + 1) * sizeof(int));

  memset(md->alias, 0, (MAXALPHA + 1) * sizeof(short));

  /* start setting actual base pair type encodings */
  switch (md->energy_set) {
    case  0:
      prepare_default_pairs(md);
      break;

    case 1:
      for (i = 1; i < MAXALPHA;) {
        md->alias[i++]  = 3;            /* A <-> G */
        md->alias[i++]  = 2;            /* B <-> C */
      }
      for (i = 1; i < MAXALPHA; i++) {
        md->pair[i][i + 1] = 2;             /* AB <-> GC */
        i++;
        md->pair[i][i - 1] = 1;             /* BA <-> CG */
      }

      break;

    case 2:
      for (i = 1; i < MAXALPHA;) {
        md->alias[i++]  = 1;            /* A <-> A*/
        md->alias[i++]  = 4;            /* B <-> U */
      }
      for (i = 1; i < MAXALPHA; i++) {
        md->pair[i][i + 1] = 5;             /* AB <-> AU */
        i++;
        md->pair[i][i - 1] = 6;             /* BA <-> UA */
      }

      break;

    case 3:
      for (i = 1; i < MAXALPHA - 2;) {
        md->alias[i++]  = 3;            /* A <-> G */
        md->alias[i++]  = 2;            /* B <-> C */
        md->alias[i++]  = 1;            /* C <-> A */
        md->alias[i++]  = 4;            /* D <-> U */
      }
      for (i = 1; i < MAXALPHA - 2; i++) {
        md->pair[i][i + 1] = 2;             /* AB <-> GC */
        i++;
        md->pair[i][i - 1] = 1;             /* BA <-> CG */
        i++;
        md->pair[i][i + 1] = 5;             /* CD <-> AU */
        i++;
        md->pair[i][i - 1] = 6;             /* DC <-> UA */
      }

      break;

    default:
      printf("vrna_md_update: "
                           "Unknown energy_set = %d. "
                           "Using defaults!",
                           md->energy_set);
      md->energy_set = 0;
      prepare_default_pairs(md);
      break;
  }

  /* set the reverse base pair types */
  for (i = 0; i <= MAXALPHA; i++)
    for (j = 0; j <= MAXALPHA; j++)
      md->rtype[md->pair[i][j]] = md->pair[j][i];

  /* handle special cases separately */
  md->rtype[0]  = 0;
  md->rtype[7]  = 7;

  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++)
      md->pair_dist[i][j] = dm_default[i][j];

  /* was used for energy_set == 0
   * for(i = 0; i < NBASES; i++)
   *  for(j = 0; j < NBASES; j++)
   *   md->rtype[md->pair[i][j]] = md->pair[j][i];
   */
}





PUBLIC void
set_model_details(vrna_md_t *md)
{
  if (md) {
    /* make sure there are no uninitialized data fields */
    memset(md, 0, sizeof(vrna_md_t));


    /* dangles: 用于处理悬挂末端的参数。
    special_hp: 特殊发夹的参数。
    noLP: 是否允许孤立碱基对。
    noGU: 是否允许 G-U 配对。
    noGUclosure: 是否允许 G-U 闭合。
    logML: 是否使用对数多分支环能量模型。
    gquad: 是否考虑 G 四链体结构。
    circ: 是否预测环状 RNA。
    uniq_ML: 是否使用唯一多分支环模型。
    compute_bpp: 是否计算碱基配对概率。
    backtrack: 回溯参数。
    backtrack_type: 回溯类型。
    energy_set: 能量集合。
    max_bp_span: 最大碱基对跨度。
    min_loop_size: 最小环大小。
    window_size: 窗口大小。
    temperature: 温度参数。
    betaScale: β 缩放因子。
    pf_smooth: 平滑因子。
    sfact: 标准因子。
    salt: 盐浓度参数。
    saltMLLower 和 saltMLUpper: 多分支环盐浓度下限和上限。
    saltDPXInit 和 saltDPXInitFact: 双链盐浓度初始化参数。
    helical_rise: 螺旋上升参数。
    backbone_length: 骨架长度参数。
    */
    md->dangles         = dangles;
    md->special_hp      = tetra_loop;
    md->noLP            = noLonelyPairs;
    md->noGU            = noGU;
    md->noGUclosure     = no_closingGU;
    md->logML           = logML;
    md->gquad           = gquad;
    md->circ            = circ;
    md->uniq_ML         = uniq_ML;
    md->compute_bpp     = do_backtrack;
    md->backtrack       = VRNA_MODEL_DEFAULT_BACKTRACK;
    md->backtrack_type  = backtrack_type;
    md->energy_set      = energy_set;
    md->max_bp_span     = max_bp_span;
    md->min_loop_size   = TURN;
    md->window_size     = VRNA_MODEL_DEFAULT_WINDOW_SIZE;
    md->oldAliEn        = oldAliEn;
    md->ribo            = ribo;
    md->cv_fact         = cv_fact;
    md->nc_fact         = nc_fact;
    md->temperature     = temperature;
    md->betaScale       = VRNA_MODEL_DEFAULT_BETA_SCALE;
    md->pf_smooth       = VRNA_MODEL_DEFAULT_PF_SMOOTH;
    md->sfact           = 1.07;
    md->salt            = defaults.salt;
    md->saltMLLower     = defaults.saltMLLower;
    md->saltMLUpper     = defaults.saltMLUpper;
    md->saltDPXInit     = defaults.saltDPXInit;
    md->saltDPXInitFact = defaults.saltDPXInitFact;
    md->helical_rise    = defaults.helical_rise;
    md->backbone_length = defaults.backbone_length;
    
    fill_pair_matrices(md);
  }
}

/* expn and kn definition*/
#define BIG 1.44115188075855872e+17
#define MAXLOG 7.09782712893383996843e+02
#define MACHEP 1.11022302462515654042e-16
#define MAXNUM 1.7976931348623158e+308
#define EUL 0.57721566490153286060
#define PI 3.14159265358979323846

double expn(int n, double x) {
    if (n < 0 || x < 0) {
        fprintf(stderr, "Invalid input: n and x must be non-negative\n");
        exit(EXIT_FAILURE);
    }

    if (x > MAXLOG) {
        return 0.0;
    }

    if (x == 0.0) {
        if (n < 2) {
            fprintf(stderr, "Singularity: E_n(x) is undefined for n < 2\n");
            exit(EXIT_FAILURE);
        } else {
            return 1.0 / (n - 1.0);
        }
    }

    if (n == 0) {
        return exp(-x) / x;
    }

    if (n > 5000) {
        double xk = x + n;
        double yk = 1.0 / (xk * xk);
        double t = n;
        double ans = yk * t * (6.0 * x * x - 8.0 * t * x + t * t);
        ans = yk * (ans + t * (t - 2.0 * x));
        ans = yk * (ans + t);
        ans = (ans + 1.0) * exp(-x) / xk;
        return ans;
    }

    if (x > 1.0) {
        double psi = -EUL - log(x);
        for (int i = 1; i < n; ++i) {
            psi += 1.0 / i;
        }

        double z = -x;
        double xk = 0.0;
        double yk = 1.0;
        double pk = 1.0 - n;
        double ans = n != 1 ? 1.0 / pk : 0.0;

        while (1) {
            xk += 1.0;
            yk *= z / xk;
            pk += 1.0;
            if (pk != 0.0) {
                ans += yk / pk;
            }
            double t = fabs(yk / ans);
            if (t <= MACHEP) {
                break;
            }
        }

        double t = (double)n;
        double r = n - 1;
        ans = (pow(z, r) * psi / tgamma(t)) - ans;
        return ans;
    }

    int k = 1;
    double pkm2 = 1.0;
    double qkm2 = x;
    double pkm1 = 1.0;
    double qkm1 = x + n;
    double ans = pkm1 / qkm1;

    while (1) {
        ++k;
        double yk = (k & 1) ? 1.0 : x;
        double xk = (k & 1) ? n + (k - 1) / 2 : k / 2;

        double pk = pkm1 * yk + pkm2 * xk;
        double qk = qkm1 * yk + qkm2 * xk;

        if (qk != 0) {
            double r = pk / qk;
            double t = fabs((ans - r) / r);
            ans = r;
            if (t <= MACHEP) {
                break;
            }
        } else {
            if (fabs(pk) > BIG) {
                pkm2 /= BIG;
                pkm1 /= BIG;
                qkm2 /= BIG;
                qkm1 /= BIG;
            }
        }
    }

    ans *= exp(-x);
    return ans;
}

double kn(int nn, double x) {
    int n = abs(nn);

    if (x <= 0.0) {
        if (x < 0.0) {
            fprintf(stderr, "kn: Domain error\n");
            exit(EXIT_FAILURE);
        } else {
            fprintf(stderr, "kn: Singularity\n");
            exit(EXIT_FAILURE);
        }
    }

    if (x > 9.55) {
        // Asymptotic expansion for large x
        double k = (double)n;
        double pn = 4.0 * k * k;
        double pk = 1.0;
        double z0 = 8.0 * x;
        double fn = 1.0;
        double t = 1.0;
        double s = t;
        double nkf = MAXNUM;
        int i = 0;

        while (1) {
            double z = pn - pk * pk;
            t *= z / (fn * z0);
            double nk1f = fabs(t);
            if (i >= n && nk1f > nkf) {
                break;
            }
            nkf = nk1f;
            s += t;
            fn += 1.0;
            pk += 2.0;
            ++i;
            if (fabs(t / s) <= MACHEP) {
                break;
            }
        }

        double ans = exp(-x) * sqrt(PI / (2.0 * x)) * s;
        return ans;
    }

    double ans = 0.0;
    double z0 = 0.25 * x * x;
    double fn = 1.0;
    double pn = 0.0;
    double zmn = 1.0;
    double tox = 2.0 / x;

    if (n > 0) {
        pn = -EUL;
        double k = 1.0;
        for (int i = 1; i < n; ++i) {
            pn += 1.0 / k;
            k += 1.0;
            fn *= k;
        }

        zmn = tox;

        if (n == 1) {
            ans = 1.0 / x;
        } else {
            double nk1f = fn / n;
            double kf = 1.0;
            double s = nk1f;
            double z = -z0;
            double zn = 1.0;
            for (int i = 1; i < n; ++i) {
                nk1f = nk1f / (n - i);
                kf *= i;
                zn *= z;
                double t = nk1f * zn / kf;
                s += t;
                if ((MAXNUM - fabs(t)) < fabs(s)) {
                    fprintf(stderr, "kn: Overflow\n");
                    exit(EXIT_FAILURE);
                }
                if (tox > 1.0 && (MAXNUM / tox) < zmn) {
                    fprintf(stderr, "kn: Overflow\n");
                    exit(EXIT_FAILURE);
                }
                zmn *= tox;
            }

            s *= 0.5;
            double t = fabs(s);
            if (zmn > 1.0 && (MAXNUM / zmn) < t) {
                fprintf(stderr, "kn: Overflow\n");
                exit(EXIT_FAILURE);
            }
            if (t > 1.0 && (MAXNUM / t) < zmn) {
                fprintf(stderr, "kn: Overflow\n");
                exit(EXIT_FAILURE);
            }
            ans = s * zmn;
        }
    }

    double tlg = 2.0 * log(0.5 * x);
    double pk = -EUL;
    double t;
    if (n == 0) {
        pn = pk;
        t = 1.0;
    } else {
        pn += 1.0 / n;
        t = 1.0 / fn;
    }

    double s = (pk + pn - tlg) * t;
    double k = 1.0;
    while (1) {
        t *= z0 / (k * (k + n));
        pk += 1.0 / k;
        pn += 1.0 / (k + n);
        s += (pk + pn - tlg) * t;
        k += 1.0;
        if (fabs(t / s) <= MACHEP) {
            break;
        }
    }

    s = 0.5 * s / zmn;
    if (n % 2 != 0) {
        s = -s;
    }
    ans += s;
    return ans;
}



#define EINVAL 22
#define ENOMEM 12
extern int *__errno_location (void) __THROW __attribute_const__;
# define errno (*__errno_location ())

#ifndef WITH_DMALLOC
/* include the following two functions only if not including <dmalloc.h> */
PUBLIC void *
vrna_alloc(unsigned size)
{
  void *pointer;

  if ((pointer = (void *)calloc(1, (size_t)size)) == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "vrna_alloc: requested size: %d\n", size);
      printf("Memory allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    printf("Memory allocation failure -> no memory");
  }

  return pointer;
}


PUBLIC void *
vrna_realloc(void     *p,
             unsigned size)
{
  if (p == NULL)
    return vrna_alloc(size);

  p = (void *)realloc(p, size);
  if (p == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "vrna_realloc: requested size: %d\n", size);
      printf("vrna_realloc allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    printf("vrna_realloc allocation failure -> no memory");
  }

  return p;
}


#endif

#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

PRIVATE INLINE double
epsilonr(double T)
{
	return 5321/T + 233.76 - 0.9297*T + 1.417*T*T/1000 - 0.8292*T*T*T/1000000;
}

PRIVATE INLINE double
bjerrum_length(double T)
{
	return 167100.052/(T*epsilonr(T));
};


PRIVATE INLINE double
ionic_strength(double rho)
{
  return rho;
};

PRIVATE INLINE double
kappa(double rho, double T)
{
  return sqrt(bjerrum_length(T)*ionic_strength(rho))/8.1284;
};

#define MIN2(A,B) ((A) < (B) ? (A) : (B))

PRIVATE double
approx_hyper(double y)
{
	double a, b, c;

	a = 1/(pow(y,6.)/pow(2*PI,6.)+1);
	b = pow(y, 4.)/(36*pow(PI, 4.)) - pow(y, 3.)/(24*PI*PI) + y*y/(2*PI*PI) - y/2;
	c = log(2*PI/y) - 1.96351;

	return a*b + (1-a)*c;
};



PRIVATE INLINE double
tau_ds(double T, double Helical_Rise)
{
  double bjerrum_length_inv, helical_rise_inv;

  bjerrum_length_inv = 1. / bjerrum_length(T);
  helical_rise_inv = 1. / Helical_Rise;
  return MIN2(helical_rise_inv, bjerrum_length_inv);
};


PRIVATE INLINE double
tau_ss(double T, double backbonelen)
{
  double bjerrum_length_inv;

  bjerrum_length_inv = 1/bjerrum_length(T);
  return MIN2(1/backbonelen, bjerrum_length_inv);
};

#define Eular_const 0.58
PRIVATE double
loop_salt_aux(double kmlss, int L, double T, double backbonelen)
{
	double a, b;

	a = (GASCONST / 1000.) * T * bjerrum_length(T) * L * backbonelen * tau_ss(T, backbonelen) * tau_ss(T, backbonelen);
	b = log(kmlss) - log(PI/2) + Eular_const + approx_hyper(kmlss) + 1/kmlss * (1 - exp(-kmlss) + kmlss*expn(1, kmlss));

	return a*b*100;
};

PUBLIC double
vrna_salt_loop(int L, double rho, double T, double backbonelen)
{
	if (L == 0)
		return 0;
	double correction, kmlss, kmlss_ref;
	
	kmlss_ref = kappa(VRNA_MODEL_DEFAULT_SALT, T) * L * backbonelen;
	kmlss = kappa(rho, T) * L * backbonelen;

	correction = loop_salt_aux(kmlss, L, T, backbonelen) - loop_salt_aux(kmlss_ref, L, T, backbonelen);

	return correction;
};
#define roundint(x) ((int) (x + 0.5 - (x<0)))
PUBLIC int
vrna_salt_loop_int(int L, double rho, double T, double backbonelen)
{
	double correction;

	correction = vrna_salt_loop(L, rho, T, backbonelen);
	return roundint(correction);
}



PRIVATE INLINE double
pairing_salt_const(double T, double Helical_Rise)
{
  return 2*(GASCONST / 1000.)*T*bjerrum_length(T)*Helical_Rise*tau_ds(T, Helical_Rise)*tau_ds(T, Helical_Rise);
};

#define Rods_dist 20.

PUBLIC int
vrna_salt_stack(double rho, double T, double hrise)
{
	double correction, kn_ref, kappa_ref;
	
	kappa_ref = kappa(VRNA_MODEL_DEFAULT_SALT, T);
	kn_ref = kn(0, Rods_dist*kappa_ref);
	correction = 100*pairing_salt_const(T, hrise)*(kn(0, Rods_dist*kappa(rho, T)) - kn_ref);
	return roundint(correction);
}

PUBLIC void
vrna_salt_ml(double saltLoop[], int lower, int upper, int* m, int* b)
{
	int sumx, sumxx;
	double y, sumy, sumyy, sumxy, denom, dm, db;

	sumx = sumxx = 0;
	sumy = sumyy = sumxy = 0.;

	for (int i=lower; i<=upper; i++)
	{ 
		sumx += i;
		sumxx += i * i;
		
		y = saltLoop[i];

		sumxy += i * y;
		sumy  += y;    
		sumyy += y*y;
	}

	denom = (double) ((upper-lower+1)*sumxx - sumx*sumx);
	dm = ((upper-lower+1) * sumxy  -  sumx * sumy) / denom;
	db = (sumy * sumxx  -  sumx * sumxy) / denom;

	*m = roundint(dm);
	*b = roundint(db);
	
}


PUBLIC vrna_md_t *
vrna_md_copy(vrna_md_t        *md_to,
             const vrna_md_t  *md_from)
{
  int       i;
  vrna_md_t *md;

  md = NULL;

  /* only process if md_from is non-NULL */
  if (md_from) {
    if (!md_to)
      /* create container to be filled */
      md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    else
      /* or directly write to target */
      md = md_to;

    /* check if not the same object */
    if (md_to != md_from) {
      /* copy simple members */
      memcpy(md, md_from, sizeof(vrna_md_t));
      /* copy arrays */
      memcpy(md->rtype, &(md_from->rtype[0]), 8 * sizeof(int));
      memcpy(md->alias, &(md_from->alias[0]), (MAXALPHA + 1) * sizeof(short));
      memcpy(md->nonstandards, &(md_from->nonstandards[0]), 64 * sizeof(char));
      /* copy matrices */
      for (i = 0; i <= MAXALPHA; i++)
        memcpy(md->pair[i], (md_from->pair[i]), (MAXALPHA + 1) * sizeof(int));
      /* copy pair dists */
      for (i = 0; i <= 6; i++)
        memcpy(md->pair_dist[i], (md_from->pair_dist[i]), 7 * sizeof(float));
    }
  }

  return md;
}

PUBLIC void
vrna_md_set_default(vrna_md_t *md)
{
  if (md) /* copy defaults */
    vrna_md_copy(md, &defaults);
}

/* Function to compute the salt correction for duplex initialization */
/* Fitted from 18 duplexes data (Chen & Znosko, 2013) */
PUBLIC int
vrna_salt_duplex_init(vrna_md_t *md_p)
{
  double    a, x, penalty;
  vrna_md_t md;

  if (md_p == NULL) {
    vrna_md_set_default(&md);
    md_p = &md;
  }

  if (md_p->saltDPXInit != 99999) {
    return md_p->saltDPXInit;
  } else {
    x = log(md_p->salt / VRNA_MODEL_DEFAULT_SALT);
    /* Converged duplex init correction */
    /* a = -1.25480589; */
    /* double b = -0.05306256; */
    /* int c = 160; */
    /* penalty = -exp(a*x*x+b*x+log(c)) + c; */
    /* a = -100.14040781; */
    return roundint(x * md_p->saltDPXInitFact);
  }
}

// #define lxc37  107.9
// #define VRNA_MODEL_DEFAULT_SALT  1.021
#define K0  273.15
#define Tmeasure  37 + K0
// #define MAX_NINIO  300
// #define GASCONST  1.98717
// #define GQuadLayerMismatchMax  1
// #define ninio37   60
// #define niniodH   320
// #define TripleC37   100
// #define TripleCdH   1860
// #define MultipleCA37   30
// #define MultipleCAdH   340
// #define MultipleCB37   160
// #define MultipleCBdH   760
// #define GQuadAlpha37   -1800
// #define GQuadAlphadH   -11934
// #define GQuadBeta37   1200
// #define GQuadBetadH   0
// #define GQuadLayerMismatch37     300
// #define GQuadLayerMismatchH      0
// #define GQuadLayerMismatchMax    1
// #define TerminalAU37   50
// #define TerminalAUdH   370
// #define DuplexInit37   410
// #define DuplexInitdH   360
PUBLIC double lxc37=107.856;
PUBLIC int ML_intern37=-90;
PUBLIC int ML_interndH=-220;
PUBLIC int ML_closing37=930;
PUBLIC int ML_closingdH=3000;
PUBLIC int ML_BASE37=0;
PUBLIC int ML_BASEdH=0;
PUBLIC int MAX_NINIO=300;
PUBLIC int ninio37=60;
PUBLIC int niniodH=320;
PUBLIC int TerminalAU37=50;
PUBLIC int TerminalAUdH=370;
PUBLIC int DuplexInit37=410;
PUBLIC int DuplexInitdH=360;
PUBLIC int TripleC37=100;
PUBLIC int TripleCdH=1860;
PUBLIC int MultipleCA37=30;
PUBLIC int MultipleCAdH=340;
PUBLIC int MultipleCB37=160;
PUBLIC int MultipleCBdH=760;

PUBLIC int GQuadAlpha37 = -1800;
PUBLIC int GQuadAlphadH = -11934;
PUBLIC int GQuadBeta37 = 1200;
PUBLIC int GQuadBetadH = 0;
PUBLIC int GQuadLayerMismatch37   = 300;
PUBLIC int GQuadLayerMismatchH    = 0;
PUBLIC int GQuadLayerMismatchMax  = 1;

PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] =
{{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
,{   INF,  -240,  -330,  -210,  -140,  -210,  -210,  -140}
,{   INF,  -330,  -340,  -250,  -150,  -220,  -240,  -150}
,{   INF,  -210,  -250,   130,   -50,  -140,  -130,   130}
,{   INF,  -140,  -150,   -50,    30,   -60,  -100,    30}
,{   INF,  -210,  -220,  -140,   -60,  -110,   -90,   -60}
,{   INF,  -210,  -240,  -130,  -100,   -90,  -130,   -90}
,{   INF,  -140,  -150,   130,    30,   -60,   -90,   130}};
PUBLIC int stackdH[NBPAIRS+1][NBPAIRS+1] =
{{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
,{   INF, -1060, -1340, -1210,  -560, -1050, -1040,  -560}
,{   INF, -1340, -1490, -1260,  -830, -1140, -1240,  -830}
,{   INF, -1210, -1260, -1460, -1350,  -880, -1280,  -880}
,{   INF,  -560,  -830, -1350,  -930,  -320,  -700,  -320}
,{   INF, -1050, -1140,  -880,  -320,  -940,  -680,  -320}
,{   INF, -1040, -1240, -1280,  -700,  -680,  -770,  -680}
,{   INF,  -560,  -830,  -880,  -320,  -320,  -680,  -320}};

PUBLIC int hairpin37[31] = {   INF,   INF,   INF,   540,   560,   570,   540,   600,   550,   640,   650,   660,   670,   680,   690,   690,   700,   710,   710,   720,   720,   730,   730,   740,   740,   750,   750,   750,   760,   760,   770};
PUBLIC int hairpindH[31] = {   INF,   INF,   INF,   130,   480,   360,  -290,   130,  -290,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500};
PUBLIC int bulge37[31] = {   INF,   380,   280,   320,   360,   400,   440,   460,   470,   480,   490,   500,   510,   520,   530,   540,   540,   550,   550,   560,   570,   570,   580,   580,   580,   590,   590,   600,   600,   600,   610};
PUBLIC int bulgedH[31] = {   INF,  1060,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710};
PUBLIC int internal_loop37[31] = {   INF,   INF,   100,   100,   110,   200,   200,   210,   230,   240,   250,   260,   270,   280,   290,   290,   300,   310,   310,   320,   330,   330,   340,   340,   350,   350,   350,   360,   360,   370,   370};
PUBLIC int internal_loopdH[31] = {   INF,   INF,   -720,   -720,  -720,  -680,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130};

PUBLIC int mismatchI37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -80,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -100,     0,  -100,     0}
 ,{     0,     0,     0,     0,   -60}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -80,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -100,     0,  -100,     0}
 ,{     0,     0,     0,     0,   -60}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }};
PUBLIC int mismatchIdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   280,     0,     0,   280,     0}
 ,{     0,     0,     0,  -340,     0}
 ,{     0,     0,     0,     0,     0}
 ,{   280,  -760,     0,   280,     0}
 ,{     0,     0,     0,     0,  -580}
 }
,{{   280,     0,     0,   280,     0}
 ,{     0,     0,     0,  -340,     0}
 ,{     0,     0,     0,     0,     0}
 ,{   280,  -760,     0,   280,     0}
 ,{     0,     0,     0,     0,  -580}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }};

PUBLIC int mismatchH37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   -80,  -100,  -110,  -100,   -80}
 ,{  -140,  -150,  -150,  -140,  -150}
 ,{   -80,  -100,  -110,  -100,   -80}
 ,{  -150,  -230,  -150,  -240,  -150}
 ,{  -100,  -100,  -140,  -100,  -210}
 }
,{{   -50,  -110,   -70,  -110,   -50}
 ,{  -110,  -110,  -150,  -130,  -150}
 ,{   -50,  -110,   -70,  -110,   -50}
 ,{  -150,  -250,  -150,  -220,  -150}
 ,{  -100,  -110,  -100,  -110,  -160}
 }
,{{    20,    20,   -20,   -10,   -20}
 ,{    20,    20,   -50,   -30,   -50}
 ,{   -10,   -10,   -20,   -10,   -20}
 ,{   -50,  -100,   -50,  -110,   -50}
 ,{   -10,   -10,   -30,   -10,  -100}
 }
,{{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -50,   -30,   -60,   -30}
 ,{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -90,   -30,  -110,   -30}
 ,{   -10,   -20,   -10,   -20,   -90}
 }
,{{   -10,   -10,   -20,   -10,   -20}
 ,{   -30,   -30,   -50,   -30,   -50}
 ,{   -10,   -10,   -20,   -10,   -20}
 ,{   -50,  -120,   -50,  -110,   -50}
 ,{   -10,   -10,   -30,   -10,  -120}
 }
,{{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -50,   -30,   -50,   -30}
 ,{     0,   -20,   -10,   -20,     0}
 ,{   -30,  -150,   -30,  -150,   -30}
 ,{   -10,   -20,   -10,   -20,   -90}
 }
,{{    20,    20,   -10,   -10,     0}
 ,{    20,    20,   -30,   -30,   -30}
 ,{     0,   -10,   -10,   -10,     0}
 ,{   -30,   -90,   -30,  -110,   -30}
 ,{   -10,   -10,   -10,   -10,   -90}
 }};
PUBLIC int mismatchHdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   560,  -570,   560,  -560,  -270}
 ,{  -560,  -910,  -560,  -560,  -560}
 ,{  -270,  -570,  -340,  -570,  -270}
 ,{   560, -1400,   560,  -920,  -560}
 ,{  -530,  -570,  -530,  -570, -1440}
 }
,{{    50,  -520,    50,  -560,  -400}
 ,{  -400,  -520,  -400,  -560,  -400}
 ,{    50,  -720,    50,  -720,  -420}
 ,{  -400, -1290,  -400,  -620,  -400}
 ,{   -30,  -720,   -30,  -720, -1080}
 }
,{{   970,   140,   970,   140,   570}
 ,{   570,    30,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,  -270,   570,    20,   570}
 ,{   830,   140,   830,   140,   -50}
 }
,{{   230,   100,   230,   220,   190}
 ,{  -110,  -110,  -260,  -520,  -260}
 ,{   190,   -60,  -140,   -60,   190}
 ,{   220,   100,  -260,   220,  -260}
 ,{   230,   -60,   230,   -60,   -70}
 }
,{{   970,   140,   970,   140,   570}
 ,{   570,   -20,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,  -520,   570,    20,   570}
 ,{   830,   140,   830,   140,  -380}
 }
,{{   230,   -30,   230,   -60,   190}
 ,{   -30,   -30,  -260,  -520,  -260}
 ,{   190,   -60,  -140,   -60,   190}
 ,{  -260,  -590,  -260,  -520,  -260}
 ,{   230,   -60,   230,   -60,   -70}
 }
,{{   970,   140,   970,   220,   570}
 ,{   570,    30,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,   100,   570,   220,   570}
 ,{   830,   140,   830,   140,   -50}
 }};

/* mismatch_multi */
PUBLIC int mismatchM37[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {   -50,  -110,   -50,  -140,   -70}
 ,{  -110,  -110,  -110,  -160,  -110}
 ,{   -70,  -150,   -70,  -150,  -100}
 ,{  -110,  -130,  -110,  -140,  -110}
 ,{   -50,  -150,   -50,  -150,   -70}
 },
 { /* GC.. */
  {   -80,  -140,   -80,  -140,  -100}
 ,{  -100,  -150,  -100,  -140,  -100}
 ,{  -110,  -150,  -110,  -150,  -140}
 ,{  -100,  -140,  -100,  -160,  -100}
 ,{   -80,  -150,   -80,  -150,  -120}
 },
 { /* GU.. */
  {   -50,   -80,   -50,   -50,   -50}
 ,{   -50,  -100,   -70,   -50,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,   -80,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UG.. */
  {   -30,   -30,   -60,   -60,   -60}
 ,{   -30,   -30,   -60,   -60,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,  -100,   -70,  -100,   -60}
 },
 { /* AU.. */
  {   -50,   -80,   -50,   -80,   -50}
 ,{   -70,  -100,   -70,  -110,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,  -120,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UA.. */
  {   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 },
 { /* NN.. */
  {   -30,   -30,   -50,   -50,   -50}
 ,{   -30,   -30,   -60,   -50,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -50,   -80,   -50,   -80,   -50}
 }};

/* mismatch_multi_enthalpies */
PUBLIC int mismatchMdH[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {    50,  -400,    50,  -400,   -30}
 ,{  -520,  -520,  -720,  -710,  -720}
 ,{    50,  -400,    50,  -400,   -30}
 ,{  -560,  -560,  -720,  -620,  -720}
 ,{  -400,  -400,  -420,  -400,  -500}
 },
 { /* GC.. */
  {  -270,  -560,  -270,  -560,  -530}
 ,{  -570,  -910,  -570,  -820,  -570}
 ,{  -340,  -560,  -340,  -560,  -530}
 ,{  -560,  -560,  -570,  -920,  -570}
 ,{  -270,  -560,  -270,  -560,  -860}
 },
 { /* GU.. */
  {   310,  -480,  -180,   310,   140}
 ,{   310,  -480,  -430,   310,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -150,  -890,  -430,  -150,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UG.. */
  {   600,   200,   600,   200,   460}
 ,{   -60,  -340,  -230,   -60,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,   160}
 },
 { /* AU.. */
  {   140,  -400,  -180,  -380,   140}
 ,{  -380,  -400,  -430,  -380,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -430,  -890,  -430,  -890,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UA.. */
  {   600,   200,   600,   200,   460}
 ,{  -230,  -390,  -230,  -310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,  -170}
 },
 { /* NN.. */
  {   600,   200,   600,   310,   460}
 ,{   310,  -340,  -230,   310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -150,  -350,  -230,  -150,  -230}
 ,{   200,   200,   -30,   200,   160}
 }};

PUBLIC int mismatch1nI37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }};
PUBLIC int mismatch1nIdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }};

PUBLIC int mismatch23I37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -50,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -110,     0,   -70,     0}
 ,{     0,     0,     0,     0,   -30}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -120,     0,   -70,     0}
 ,{     0,     0,     0,     0,   -30}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    20,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    20,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }};
PUBLIC int mismatch23IdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,  -570,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -860,     0,  -900,     0}
 ,{     0,     0,     0,     0,  -640}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0, -1090,     0,  -900,     0}
 ,{     0,     0,     0,     0,  -640}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -580,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   -60,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -580,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   -60,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }};

/* mismatch_exterior */
PUBLIC int mismatchExt37[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {   -50,  -110,   -50,  -140,   -70}
 ,{  -110,  -110,  -110,  -160,  -110}
 ,{   -70,  -150,   -70,  -150,  -100}
 ,{  -110,  -130,  -110,  -140,  -110}
 ,{   -50,  -150,   -50,  -150,   -70}
 },
 { /* GC.. */
  {   -80,  -140,   -80,  -140,  -100}
 ,{  -100,  -150,  -100,  -140,  -100}
 ,{  -110,  -150,  -110,  -150,  -140}
 ,{  -100,  -140,  -100,  -160,  -100}
 ,{   -80,  -150,   -80,  -150,  -120}
 },
 { /* GU.. */
  {   -50,   -80,   -50,   -50,   -50}
 ,{   -50,  -100,   -70,   -50,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,   -80,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UG.. */
  {   -30,   -30,   -60,   -60,   -60}
 ,{   -30,   -30,   -60,   -60,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,  -100,   -70,  -100,   -60}
 },
 { /* AU.. */
  {   -50,   -80,   -50,   -80,   -50}
 ,{   -70,  -100,   -70,  -110,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,  -120,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UA.. */
  {   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 },
 { /* NN.. */
  {   -30,   -30,   -50,   -50,   -50}
 ,{   -30,   -30,   -60,   -50,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -50,   -80,   -50,   -80,   -50}
 }};

/* mismatch_exterior_enthalpies */
PUBLIC int mismatchExtdH[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {    50,  -400,    50,  -400,   -30}
 ,{  -520,  -520,  -720,  -710,  -720}
 ,{    50,  -400,    50,  -400,   -30}
 ,{  -560,  -560,  -720,  -620,  -720}
 ,{  -400,  -400,  -420,  -400,  -500}
 },
 { /* GC.. */
  {  -270,  -560,  -270,  -560,  -530}
 ,{  -570,  -910,  -570,  -820,  -570}
 ,{  -340,  -560,  -340,  -560,  -530}
 ,{  -560,  -560,  -570,  -920,  -570}
 ,{  -270,  -560,  -270,  -560,  -860}
 },
 { /* GU.. */
  {   310,  -480,  -180,   310,   140}
 ,{   310,  -480,  -430,   310,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -150,  -890,  -430,  -150,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UG.. */
  {   600,   200,   600,   200,   460}
 ,{   -60,  -340,  -230,   -60,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,   160}
 },
 { /* AU.. */
  {   140,  -400,  -180,  -380,   140}
 ,{  -380,  -400,  -430,  -380,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -430,  -890,  -430,  -890,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UA.. */
  {   600,   200,   600,   200,   460}
 ,{  -230,  -390,  -230,  -310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,  -170}
 },
 { /* NN.. */
  {   600,   200,   600,   310,   460}
 ,{   310,  -340,  -230,   310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -150,  -350,  -230,  -150,  -230}
 ,{   200,   200,   -30,   200,   160}
 }};

/* dangle5 */
PUBLIC int dangle5_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -10,   -50,   -30,   -20,   -10},
/* GC */ {    -0,   -20,   -30,    -0,    -0},
/* GU */ {   -20,   -30,   -30,   -40,   -20},
/* UG */ {   -10,   -30,   -10,   -20,   -20},
/* AU */ {   -20,   -30,   -30,   -40,   -20},
/* UA */ {   -10,   -30,   -10,   -20,   -20},
/* NN */ {    -0,   -20,   -10,    -0,    -0}
};

/* dangle3 */
PUBLIC int dangle3_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -40,  -110,   -40,  -130,   -60},
/* GC */ {   -80,  -170,   -80,  -170,  -120},
/* GU */ {   -10,   -70,   -10,   -70,   -10},
/* UG */ {   -50,   -80,   -50,   -80,   -60},
/* AU */ {   -10,   -70,   -10,   -70,   -10},
/* UA */ {   -50,   -80,   -50,   -80,   -60},
/* NN */ {   -10,   -70,   -10,   -70,   -10}
};

/* dangle5_enthalpies */
PUBLIC int dangle5_dH[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   330,  -240,   330,    80,  -140},
/* GC */ {    70,  -160,    70,  -460,   -40},
/* GU */ {   310,   160,   220,    70,   310},
/* UG */ {   690,   -50,   690,    60,    60},
/* AU */ {   310,   160,   220,    70,   310},
/* UA */ {   690,   -50,   690,    60,    60},
/* NN */ {   690,   160,   690,    80,   310}
};

/* dangle3_enthalpies */
PUBLIC int dangle3_dH[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {  -280,  -740,  -280,  -640,  -360},
/* GC */ {  -410,  -900,  -410,  -860,  -750},
/* GU */ {   -70,  -570,   -70,  -580,  -220},
/* UG */ {   -90,  -490,   -90,  -550,  -230},
/* AU */ {   -70,  -570,   -70,  -580,  -220},
/* UA */ {   -90,  -490,   -90,  -550,  -230},
/* NN */ {   -70,  -490,   -70,  -550,  -220}
};

PUBLIC char Triloops[241] =
  "CAACG "
  "GUUAC "
;
PUBLIC int Triloop37[40] = {   680,   690};
PUBLIC int TriloopdH[40] = {  2370,  1080};

PUBLIC char Tetraloops[281] =
  "CAACGG "
  "CCAAGG "
  "CCACGG "
  "CCCAGG "
  "CCGAGG "
  "CCGCGG "
  "CCUAGG "
  "CCUCGG "
  "CUAAGG "
  "CUACGG "
  "CUCAGG "
  "CUCCGG "
  "CUGCGG "
  "CUUAGG "
  "CUUCGG "
  "CUUUGG "
;
PUBLIC int Tetraloop37[40] = {   550,   330,   370,   340,   350,   360,   370,   250,   360,   280,   370,   270,   280,   350,   370,   370};
PUBLIC int TetraloopdH[40] = {   690, -1030,  -330,  -890,  -660,  -750,  -350, -1390,  -760, -1070,  -660, -1290, -1070,  -620, -1530,  -680};

PUBLIC char Hexaloops[361] =
  "ACAGUACU "
  "ACAGUGAU "
  "ACAGUGCU "
  "ACAGUGUU "
;
PUBLIC int Hexaloop37[40] = {   280,   360,   290,   180};
PUBLIC int HexaloopdH[40] = { -1680, -1140, -1280, -1540};
#define saltT md->temperature+K0

/* long matrices*/
PUBLIC int int11_37[NBPAIRS+1][NBPAIRS+1][5][5] =
{{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{    90,    90,    50,    50,    50}
  ,{    90,    90,    50,    50,    50}
  ,{    50,    50,    50,    50,    50}
  ,{    50,    50,    50,  -140,    50}
  ,{    50,    50,    50,    50,    40}
  }
 ,{{    90,    90,    50,    50,    60}
  ,{    90,    90,   -40,    50,    50}
  ,{    60,    30,    50,    50,    60}
  ,{    50,   -10,    50,  -220,    50}
  ,{    50,    50,     0,    50,   -10}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,    60,    50,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   -20,   120,  -140,   120}
  ,{   120,   120,   100,   120,   110}
  }
 ,{{   220,   220,   170,   120,   120}
  ,{   220,   220,   130,   120,   120}
  ,{   170,   120,   170,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   110}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,    80}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   120}
  }
 ,{{   220,   220,   170,   120,   120}
  ,{   220,   220,   130,   120,   120}
  ,{   170,   120,   170,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   120}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{    90,    90,    60,    50,    50}
  ,{    90,    90,    30,   -10,    50}
  ,{    50,   -40,    50,    50,     0}
  ,{    50,    50,    50,  -220,    50}
  ,{    60,    50,    60,    50,   -10}
  }
 ,{{    80,    80,    50,    50,    50}
  ,{    80,    80,    50,    50,    50}
  ,{    50,    50,    50,    50,    50}
  ,{    50,    50,    50,  -230,    50}
  ,{    50,    50,    50,    50,   -60}
  }
 ,{{   190,   190,   120,   150,   150}
  ,{   190,   190,   120,   150,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   150,   120,   120,   120,   150}
  }
 ,{{   160,   160,   120,   120,   120}
  ,{   160,   160,   120,   100,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,    70}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,    80}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   120}
  }
 ,{{   190,   190,   120,   150,   150}
  ,{   190,   190,   120,   150,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   150,   120,   120,   120,   150}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,    60,   120,   -20,   120}
  ,{   120,    50,   120,   120,   100}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   110}
  }
 ,{{   190,   190,   120,   120,   150}
  ,{   190,   190,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   150,   150,   120,  -140,   120}
  ,{   150,   120,   120,   120,   150}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   120}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   120}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   220,   220,   170,   120,   120}
  ,{   220,   220,   120,   120,   120}
  ,{   170,   130,   170,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   110}
  }
 ,{{   160,   160,   120,   120,   120}
  ,{   160,   160,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   100,   120,  -140,   120}
  ,{   120,   120,   120,   120,    70}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 ,{{   220,   220,   190,   190,   190}
  ,{   220,   220,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,    80}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,    80}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   120}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   120}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   150}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   120}
  }
 ,{{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   120}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   150}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   170}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   220,   220,   170,   120,   120}
  ,{   220,   220,   120,   120,   120}
  ,{   170,   130,   170,   120,   120}
  ,{   120,   120,   120,  -140,   120}
  ,{   120,   120,   120,   120,   120}
  }
 ,{{   190,   190,   120,   120,   150}
  ,{   190,   190,   120,   120,   120}
  ,{   120,   120,   120,   120,   120}
  ,{   150,   150,   120,  -140,   120}
  ,{   150,   120,   120,   120,   150}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   220,   220,   190,   190,   190}
  ,{   220,   220,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   160}
  }
 ,{{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 ,{{   220,   220,   190,   190,   190}
  ,{   220,   220,   190,   190,   190}
  ,{   190,   190,   190,   190,   190}
  ,{   190,   190,   190,   -70,   190}
  ,{   190,   190,   190,   190,   190}
  }
 }};
PUBLIC int int11_dH[NBPAIRS+1][NBPAIRS+1][5][5] =
{{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1840, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  }
 ,{{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1840, -1050}
  ,{ -1050, -1050, -1050, -1050, -1390}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1840, -1050}
  ,{ -1050, -1050, -1050, -1050, -1390}
  }
 ,{{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1050, -1050}
  ,{ -1050, -1050, -1050, -1840, -1050}
  ,{ -1050, -1050, -1050, -1050, -1730}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550, -1230}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550, -1230}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550, -1230}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -730}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -730}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550, -1230}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -730}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -730}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 }
,{{{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  ,{   INF,   INF,   INF,   INF,   INF}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  }
 ,{{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550,  -550,  -550}
  ,{  -550,  -550,  -550, -1340,  -550}
  ,{  -550,  -550,  -550,  -550,  -890}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,  -390}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 ,{{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  ,{   -50,   -50,   -50,  -830,   -50}
  ,{   -50,   -50,   -50,   -50,   -50}
  }
 }};
PUBLIC int int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5] =
{{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   110,   230}
   ,{   230,   230,   230,   110,   230}
   ,{   230,   230,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   230,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   150,   150,   150,   150,   150}
   }
  }
 ,{{{   250,   250,   250,   230,   230}
   ,{   250,   250,   230,   230,   230}
   ,{   250,   230,   250,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   250,   250,   230,   230,   230}
   }
  ,{{   250,   250,   230,   110,   230}
   ,{   250,   250,   230,   110,   230}
   ,{   230,   230,   170,   110,   230}
   ,{   110,    80,   110,   110,   110}
   ,{   230,   230,   230,   110,   230}
   }
  ,{{   250,   250,   250,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   250,   230,   250,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   250,   250,   230,   230,   230}
   }
  ,{{   230,   170,   230,   110,   230}
   ,{   230,   170,   230,    80,   230}
   ,{   230,   110,   230,   110,   230}
   ,{   120,   120,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   220,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   170,   150,   170,   150,   140}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   250,   250,   230,   230,   230}
   ,{   250,   250,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   250,   250,   230,   230,   230}
   ,{   250,   250,   230,   210,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   120,   120,   110,   110,   110}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   190,   230,   230}
   }
  ,{{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   150,   150,   150,   150,   150}
   }
  }
 ,{{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   ,{   230,   230,   230,   230,   230}
   }
  ,{{   230,   110,   230,   110,   230}
   ,{   230,   110,   230,   110,   230}
   ,{   230,   110,   230,   110,   230}
   ,{   110,   110,   110,   110,   110}
   ,{   230,   110,   230,   110,   230}
   }
  ,{{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   230,   230,   230,   230,   150}
   ,{   150,   150,   150,   150,   150}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   250,   300,   210,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   120,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   190,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   250,   300,   210,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   120,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   190,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   250,   370,   210,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   120,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   190,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   300,   300,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   370,   370,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   ,{   300,   300,   300,   300,   300}
   }
  ,{{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   300,   190,   300,   190,   300}
   ,{   190,   190,   190,   190,   190}
   ,{   300,   190,   300,   190,   300}
   }
  ,{{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   300,   300,   300,   300,   220}
   ,{   220,   220,   220,   220,   220}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 ,{{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   ,{   370,   370,   370,   370,   370}
   }
  ,{{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   370,   260,   370,   260,   370}
   ,{   260,   260,   260,   260,   260}
   ,{   370,   260,   370,   260,   370}
   }
  ,{{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   370,   370,   370,   370,   300}
   ,{   300,   300,   300,   300,   300}
   }
  }
 }};
PUBLIC int int21_dH[NBPAIRS+1][NBPAIRS+1][5][5][5] =
{{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   }
  ,{{   350,   350,   350,  -230,   350}
   ,{   350,   350,   350,  -230,   350}
   ,{   350,   350,   350,  -230,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,   350,   350,  -230,   350}
   }
  ,{{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   }
  ,{{   350,  -230,   350,  -230,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,  -230,   350,  -230,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,  -230,   350,  -230,   350}
   }
  ,{{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{  -670,  -670,  -670,  -670,  -670}
   }
  }
 ,{{{   780,   640,   780,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   780,   350,   780,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   640,   640,   350,   350,   350}
   }
  ,{{   350,   350,   350,   250,   350}
   ,{   350,   260,   350,   250,   350}
   ,{   350,   350,  -250,  -230,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,   350,   350,  -230,   350}
   }
  ,{{   780,   640,   780,   350,   350}
   ,{   350,   160,   350,   350,   350}
   ,{   780,   350,   780,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   640,   640,   350,   350,   350}
   }
  ,{{   350,  -160,   350,  -230,   350}
   ,{   350,  -160,   350,  -410,   350}
   ,{   350,  -230,   350,  -230,   350}
   ,{  -230,  -310,  -230,  -230,  -230}
   ,{   350,  -230,   350,  -230,   350}
   }
  ,{{   580,   350,   580,   350,  -580}
   ,{   350,   350,   350,   350,  -670}
   ,{   580,   350,   580,   350,  -580}
   ,{   350,   350,   350,   350,  -670}
   ,{  -670,  -670,  -690,  -670,  -700}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   690,   690,   350,   350,   350}
   ,{   690,   690,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   }
  ,{{   690,   690,   350,   350,   350}
   ,{   690,   690,   350,   240,   350}
   ,{   350,   350,   350,   350,   350}
   ,{  -230,  -500,  -230,  -230,  -230}
   ,{   350,   350,   350,   350,   350}
   }
  ,{{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   130,   350,   350}
   }
  ,{{   350,  -230,   350,  -230,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,  -230,   350,  -230,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,  -230,   350,  -230,   350}
   }
  ,{{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{  -670,  -670,  -670,  -670,  -670}
   }
  }
 ,{{{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   }
  ,{{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,   350,   350,   350,   350}
   }
  ,{{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   ,{   350,   350,   350,   350,   350}
   }
  ,{{   350,  -230,   350,  -230,   350}
   ,{   350,  -230,   350,  -230,   350}
   ,{   350,  -230,   350,  -230,   350}
   ,{  -230,  -230,  -230,  -230,  -230}
   ,{   350,  -230,   350,  -230,   350}
   }
  ,{{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{   350,   350,   350,   350,  -670}
   ,{  -670,  -670,  -670,  -670,  -670}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   690,   850,   240,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,  -500,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   130,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   690,   850,   240,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,  -500,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   130,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,   690,  1350,   240,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,  -500,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,   130,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   850,   850,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{  1350,  1350,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 }
,{{{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  ,{{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   ,{   INF,   INF,   INF,   INF,   INF}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   ,{   850,   850,   850,   850,   850}
   }
  ,{{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   850,   280,   850,   280,   850}
   ,{   280,   280,   280,   280,   280}
   ,{   850,   280,   850,   280,   850}
   }
  ,{{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{   850,   850,   850,   850,  -160}
   ,{  -160,  -160,  -160,  -160,  -160}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 ,{{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   ,{  1350,  1350,  1350,  1350,  1350}
   }
  ,{{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{  1350,   780,  1350,   780,  1350}
   ,{   780,   780,   780,   780,   780}
   ,{  1350,   780,  1350,   780,  1350}
   }
  ,{{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{  1350,  1350,  1350,  1350,   340}
   ,{   340,   340,   340,   340,   340}
   }
  }
 }};
PUBLIC int int22_37[NBPAIRS+1][NBPAIRS+1][5][5][5][5] =
{{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   200,   160,   200,   150,   200}
    ,{   200,   160,   200,   150,   200}
    ,{   180,   140,   180,   140,   180}
    ,{   200,   160,   200,   150,   200}
    ,{   170,   130,   170,   120,   170}
    }
   ,{{   160,   120,   160,   110,   160}
    ,{   160,   120,   160,   110,   160}
    ,{   150,   110,   150,   110,   150}
    ,{   110,    20,   110,    20,    90}
    ,{   150,   110,   150,   110,   150}
    }
   ,{{   200,   160,   200,   150,   200}
    ,{   200,   160,   200,   150,   200}
    ,{   180,   140,   180,   140,   180}
    ,{   200,   160,   200,   150,   200}
    ,{   170,   130,   170,   120,   170}
    }
   ,{{   150,   110,   150,   110,   150}
    ,{   110,    20,   110,    20,    90}
    ,{   150,   110,   150,   110,   150}
    ,{    80,     0,    10,    80,    20}
    ,{   150,   110,   150,   110,   150}
    }
   ,{{   200,   160,   200,   150,   200}
    ,{   200,   160,   200,   150,   200}
    ,{   170,   130,   170,   120,   170}
    ,{   200,   160,   200,   150,   200}
    ,{   100,   100,    80,    30,    80}
    }
   }
  ,{{{   200,   160,   200,   110,   200}
    ,{   200,   160,   200,    60,   200}
    ,{   180,   140,   180,   110,   180}
    ,{   200,   160,   200,    60,   200}
    ,{   170,   130,   170,    90,   170}
    }
   ,{{   160,   120,   160,    20,   160}
    ,{   160,   120,   160,    20,   160}
    ,{   150,   110,   150,    20,   150}
    ,{    60,    20,    60,   -70,    60}
    ,{   150,   110,   150,    20,   150}
    }
   ,{{   200,   160,   200,   110,   200}
    ,{   200,   160,   200,    60,   200}
    ,{   180,   140,   180,   110,   180}
    ,{   200,   160,   200,    60,   200}
    ,{   170,   130,   170,    90,   170}
    }
   ,{{   150,   110,   150,    20,   150}
    ,{    60,    20,    60,   -70,    60}
    ,{   150,   110,   150,    20,   150}
    ,{    10,   -30,    10,     0,    10}
    ,{   150,   110,   150,    20,   150}
    }
   ,{{   200,   160,   200,    90,   200}
    ,{   200,   160,   200,    60,   200}
    ,{   170,   130,   170,    90,   170}
    ,{   200,   160,   200,    60,   200}
    ,{   100,   100,    80,   -50,    80}
    }
   }
  ,{{{   180,   150,   180,   150,   170}
    ,{   180,   150,   180,   150,   170}
    ,{   170,   140,   170,   140,   150}
    ,{   180,   150,   180,   150,   170}
    ,{   150,   120,   150,   120,   140}
    }
   ,{{   140,   110,   140,   110,   130}
    ,{   140,   110,   140,   110,   130}
    ,{   140,   110,   140,   110,   120}
    ,{   110,    20,   110,    20,    90}
    ,{   140,   110,   140,   110,   120}
    }
   ,{{   180,   150,   180,   150,   170}
    ,{   180,   150,   180,   150,   170}
    ,{   170,   140,   170,   140,   150}
    ,{   180,   150,   180,   150,   170}
    ,{   150,   120,   150,   120,   140}
    }
   ,{{   140,   110,   140,   110,   120}
    ,{   110,    20,   110,    20,    90}
    ,{   140,   110,   140,   110,   120}
    ,{   -10,   -40,   -10,   -40,   -20}
    ,{   140,   110,   140,   110,   120}
    }
   ,{{   180,   150,   180,   150,   170}
    ,{   180,   150,   180,   150,   170}
    ,{   150,   120,   150,   120,   140}
    ,{   180,   150,   180,   150,   170}
    ,{    60,    30,    60,    30,    50}
    }
   }
  ,{{{   200,   110,   200,    80,   200}
    ,{   200,    60,   200,    10,   200}
    ,{   180,   110,   180,   -10,   180}
    ,{   200,    60,   200,    80,   200}
    ,{   170,    90,   170,    20,   170}
    }
   ,{{   160,    20,   160,     0,   160}
    ,{   160,    20,   160,   -30,   160}
    ,{   150,    20,   150,   -40,   150}
    ,{    60,   -70,    60,     0,    60}
    ,{   150,    20,   150,   -40,   150}
    }
   ,{{   200,   110,   200,    10,   200}
    ,{   200,    60,   200,    10,   200}
    ,{   180,   110,   180,   -10,   180}
    ,{   200,    60,   200,    10,   200}
    ,{   170,    90,   170,   -20,   170}
    }
   ,{{   150,    20,   150,    80,   150}
    ,{    60,   -70,    60,     0,    60}
    ,{   150,    20,   150,   -40,   150}
    ,{    80,     0,    10,    80,    10}
    ,{   150,    20,   150,   -40,   150}
    }
   ,{{   200,    90,   200,    20,   200}
    ,{   200,    60,   200,    10,   200}
    ,{   170,    90,   170,   -20,   170}
    ,{   200,    60,   200,    10,   200}
    ,{    80,   -50,    80,    20,    80}
    }
   }
  ,{{{   170,   150,   170,   150,   100}
    ,{   170,   150,   170,   150,   100}
    ,{   150,   140,   150,   140,    60}
    ,{   170,   150,   170,   150,    80}
    ,{   140,   120,   140,   120,    50}
    }
   ,{{   130,   110,   130,   110,   100}
    ,{   130,   110,   130,   110,   100}
    ,{   120,   110,   120,   110,    30}
    ,{    90,    20,    90,    20,   -50}
    ,{   120,   110,   120,   110,    30}
    }
   ,{{   170,   150,   170,   150,    80}
    ,{   170,   150,   170,   150,    80}
    ,{   150,   140,   150,   140,    60}
    ,{   170,   150,   170,   150,    80}
    ,{   140,   120,   140,   120,    50}
    }
   ,{{   120,   110,   120,   110,    30}
    ,{    90,    20,    90,    20,   -50}
    ,{   120,   110,   120,   110,    30}
    ,{    20,   -40,   -20,   -40,    20}
    ,{   120,   110,   120,   110,    30}
    }
   ,{{   170,   150,   170,   150,    80}
    ,{   170,   150,   170,   150,    80}
    ,{   140,   120,   140,   120,    50}
    ,{   170,   150,   170,   150,    80}
    ,{    50,    30,    50,    30,   -40}
    }
   }
  }
 ,{{{{   220,   150,   220,   140,   170}
    ,{   220,   130,   220,   130,   170}
    ,{   150,   110,   150,   110,   150}
    ,{   140,   100,   140,   100,   140}
    ,{   170,   150,   150,   140,   170}
    }
   ,{{   220,   130,   220,   130,   170}
    ,{   220,   130,   220,   130,   170}
    ,{   150,   110,   150,   100,   150}
    ,{    70,   -30,    70,   -70,    50}
    ,{   150,   110,   150,   100,   150}
    }
   ,{{   190,   110,   190,   100,   170}
    ,{   190,   110,   190,   100,   140}
    ,{   150,   110,   150,   100,   150}
    ,{   140,   100,   140,   100,   140}
    ,{   170,   110,   150,   100,   170}
    }
   ,{{   150,   110,   150,   100,   150}
    ,{   140,    70,    70,   -10,   140}
    ,{   150,   110,   150,   100,   150}
    ,{    80,   -30,    10,    80,    70}
    ,{   150,   110,   150,   100,   150}
    }
   ,{{   150,   150,   150,   140,   150}
    ,{   140,   100,   140,   100,   140}
    ,{   150,   110,   150,   110,   150}
    ,{   140,   100,   140,   100,   140}
    ,{   150,   150,    70,   140,    70}
    }
   }
  ,{{{   170,   150,   150,    90,   170}
    ,{   170,   130,   140,    10,   170}
    ,{   150,   110,   150,    80,   150}
    ,{   140,   100,   140,    10,   140}
    ,{   150,   150,   150,    90,   150}
    }
   ,{{   170,   130,   150,    10,   170}
    ,{   170,   130,    60,     0,   170}
    ,{   150,   110,   150,   -70,   150}
    ,{    10,   -30,    10,  -160,   -30}
    ,{   150,   110,   150,    10,   150}
    }
   ,{{   150,   110,   150,    70,   150}
    ,{   140,   100,    50,  -100,   140}
    ,{   150,   110,   150,   -60,   150}
    ,{   140,   100,   140,    10,   140}
    ,{   150,   110,   150,    70,   150}
    }
   ,{{   150,   110,   150,    10,   150}
    ,{    40,    40,    30,   -70,    30}
    ,{   150,   110,   150,    10,   150}
    ,{    10,   -30,   -30,     0,    10}
    ,{   150,   110,   150,    10,   150}
    }
   ,{{   150,   150,   150,    90,   150}
    ,{   140,   100,   140,    10,   140}
    ,{   150,   110,   150,    80,   150}
    ,{   140,   100,   140,    10,   140}
    ,{   150,   150,     0,    90,    70}
    }
   }
  ,{{{   220,   130,   220,   130,   170}
    ,{   220,   130,   220,   130,   140}
    ,{   140,   110,   140,   110,   120}
    ,{   130,   100,   130,   100,   110}
    ,{   170,   100,   130,   100,   170}
    }
   ,{{   220,   130,   220,   130,   140}
    ,{   220,   130,   220,   130,   140}
    ,{   130,   100,   130,   100,   120}
    ,{    70,   -70,    70,   -70,     0}
    ,{   130,   100,   130,   100,   120}
    }
   ,{{   190,   110,   190,   100,   170}
    ,{   190,   110,   190,   100,   110}
    ,{   130,   100,   130,   100,   120}
    ,{   130,   100,   130,   100,   110}
    ,{   170,   100,   130,   100,   170}
    }
   ,{{   130,   100,   130,   100,   120}
    ,{    70,    70,    70,   -10,    60}
    ,{   130,   100,   130,   100,   120}
    ,{    20,   -40,   -10,   -40,    20}
    ,{   130,   100,   130,   100,   120}
    }
   ,{{   140,   110,   140,   110,   120}
    ,{   130,   100,   130,   100,   110}
    ,{   140,   110,   140,   110,   120}
    ,{   130,   100,   130,   100,   110}
    ,{    30,   -20,   -10,    30,    20}
    }
   }
  ,{{{   170,    90,   170,   140,   170}
    ,{   170,    70,   170,   -10,   170}
    ,{   150,    80,   150,   -40,   150}
    ,{   140,    10,   140,    80,   140}
    ,{   150,    90,   150,   140,   150}
    }
   ,{{   170,    10,   170,   -10,   170}
    ,{   170,   -20,   170,   -10,   170}
    ,{   150,   -40,   150,   -40,   150}
    ,{   -30,  -170,   -30,   -90,   -30}
    ,{   150,    10,   150,   -40,   150}
    }
   ,{{   150,    70,   150,    20,   150}
    ,{   140,    70,   140,   -50,   140}
    ,{   150,    70,   150,   -40,   150}
    ,{   140,    10,   140,   -50,   140}
    ,{   150,    70,   150,    20,   150}
    }
   ,{{   150,    10,   150,    80,   150}
    ,{    30,   -50,    30,   -30,    30}
    ,{   150,    10,   150,   -40,   150}
    ,{    80,   -30,    10,    80,    10}
    ,{   150,    10,   150,   -40,   150}
    }
   ,{{   150,    90,   150,   140,   150}
    ,{   140,    10,   140,   -50,   140}
    ,{   150,    80,   150,   -50,   150}
    ,{   140,    10,   140,   -50,   140}
    ,{   140,    90,    70,   140,    70}
    }
   }
  ,{{{   140,   130,   140,   130,   140}
    ,{   140,   130,   140,   130,   140}
    ,{   120,   110,   120,   110,    30}
    ,{   110,   100,   110,   100,    70}
    ,{   120,   100,   120,   100,    30}
    }
   ,{{   140,   130,   140,   130,   140}
    ,{   140,   130,   140,   130,   140}
    ,{   120,   100,   120,   100,    30}
    ,{    50,   -70,     0,   -70,    50}
    ,{   120,   100,   120,   100,    30}
    }
   ,{{   120,   100,   120,   100,    30}
    ,{   110,   100,   110,   100,    30}
    ,{   120,   100,   120,   100,    30}
    ,{   110,   100,   110,   100,    20}
    ,{   120,   100,   120,   100,    30}
    }
   ,{{   140,   100,   120,   100,   140}
    ,{   140,   -10,    50,   -10,   140}
    ,{   120,   100,   120,   100,    30}
    ,{    70,   -40,   -60,   -40,    70}
    ,{   120,   100,   120,   100,    30}
    }
   ,{{   120,   110,   120,   110,    30}
    ,{   110,   100,   110,   100,    20}
    ,{   120,   110,   120,   110,    30}
    ,{   110,   100,   110,   100,    20}
    ,{    40,    30,    40,    30,   -60}
    }
   }
  }
 ,{{{{   300,   290,   300,   260,   300}
    ,{   300,   270,   300,   260,   300}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   290,   290,   270,   220,   270}
    }
   ,{{   300,   270,   300,   260,   300}
    ,{   300,   270,   300,   260,   300}
    ,{   270,   230,   270,   220,   270}
    ,{   230,   150,   230,   140,   220}
    ,{   270,   230,   270,   220,   270}
    }
   ,{{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    }
   ,{{   270,   230,   270,   220,   270}
    ,{   270,   190,   270,   180,   260}
    ,{   270,   230,   270,   220,   270}
    ,{   210,   130,   140,   210,   150}
    ,{   270,   230,   270,   220,   270}
    }
   ,{{   290,   290,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   290,   290,   270,   220,   270}
    }
   }
  ,{{{   300,   290,   300,   190,   300}
    ,{   300,   270,   300,   170,   300}
    ,{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   290,   290,   270,   190,   270}
    }
   ,{{   300,   270,   300,   170,   300}
    ,{   300,   270,   300,   170,   300}
    ,{   270,   230,   270,   130,   270}
    ,{   190,   150,   190,    50,   190}
    ,{   270,   230,   270,   130,   270}
    }
   ,{{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   190,   270}
    }
   ,{{   270,   230,   270,   130,   270}
    ,{   230,   190,   230,    90,   230}
    ,{   270,   230,   270,   130,   270}
    ,{   140,   100,   140,   130,   140}
    ,{   270,   230,   270,   130,   270}
    }
   ,{{   290,   290,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   290,   290,   270,   130,   270}
    }
   }
  ,{{{   290,   260,   290,   260,   270}
    ,{   290,   260,   290,   260,   270}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   290,   260,   290,   260,   270}
    ,{   290,   260,   290,   260,   270}
    ,{   250,   220,   250,   220,   240}
    ,{   230,   140,   230,   140,   220}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   270,   220,   270,   220,   260}
    ,{   270,   180,   270,   180,   260}
    ,{   250,   220,   250,   220,   240}
    ,{   120,    90,   120,    90,   110}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    }
   }
  ,{{{   300,   190,   300,   210,   300}
    ,{   300,   170,   300,   170,   300}
    ,{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,   210,   270}
    ,{   270,   190,   270,   210,   270}
    }
   ,{{   300,   170,   300,   130,   300}
    ,{   300,   170,   300,   110,   300}
    ,{   270,   130,   270,    80,   270}
    ,{   190,    50,   190,   130,   190}
    ,{   270,   130,   270,    80,   270}
    }
   ,{{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   190,   270,    80,   270}
    }
   ,{{   270,   130,   270,   210,   270}
    ,{   230,    90,   230,   170,   230}
    ,{   270,   130,   270,    80,   270}
    ,{   210,   130,   140,   210,   140}
    ,{   270,   130,   270,    80,   270}
    }
   ,{{   270,   190,   270,   210,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   130,   270,   210,   270}
    }
   }
  ,{{{   270,   260,   270,   260,   240}
    ,{   270,   260,   270,   260,   240}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   270,   260,   270,   260,   240}
    ,{   270,   260,   270,   260,   240}
    ,{   240,   220,   240,   220,   150}
    ,{   220,   140,   220,   140,    70}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   260,   220,   260,   220,   150}
    ,{   260,   180,   260,   180,   110}
    ,{   240,   220,   240,   220,   150}
    ,{   150,    90,   110,    90,   150}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    }
   }
  }
 ,{{{{   310,   260,   310,   220,   300}
    ,{   310,   230,   310,   220,   300}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   260,   260,   240,   190,   240}
    }
   ,{{   240,   200,   240,   190,   240}
    ,{   200,   160,   200,   160,   200}
    ,{   240,   200,   240,   190,   240}
    ,{   150,    60,   150,    60,   130}
    ,{   240,   200,   240,   190,   240}
    }
   ,{{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    }
   ,{{   310,   230,   310,   220,   300}
    ,{   310,   230,   310,   220,   300}
    ,{   240,   200,   240,   190,   240}
    ,{   180,   100,   110,   180,   120}
    ,{   240,   200,   240,   190,   240}
    }
   ,{{   260,   260,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   260,   260,   240,   190,   240}
    }
   }
  ,{{{   270,   260,   270,   160,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   260,   260,   240,   160,   240}
    }
   ,{{   240,   200,   240,   100,   240}
    ,{   200,   160,   200,    70,   200}
    ,{   240,   200,   240,   100,   240}
    ,{   100,    60,   100,   -30,   100}
    ,{   240,   200,   240,   100,   240}
    }
   ,{{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   160,   240}
    }
   ,{{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   240,   200,   240,   100,   240}
    ,{   110,    70,   110,   100,   110}
    ,{   240,   200,   240,   100,   240}
    }
   ,{{   260,   260,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   260,   260,   240,   100,   240}
    }
   }
  ,{{{   310,   220,   310,   220,   300}
    ,{   310,   220,   310,   220,   300}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    }
   ,{{   220,   190,   220,   190,   210}
    ,{   190,   160,   190,   160,   170}
    ,{   220,   190,   220,   190,   210}
    ,{   150,    60,   150,    60,   130}
    ,{   220,   190,   220,   190,   210}
    }
   ,{{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    }
   ,{{   310,   220,   310,   220,   300}
    ,{   310,   220,   310,   220,   300}
    ,{   220,   190,   220,   190,   210}
    ,{    90,    60,    90,    60,    80}
    ,{   220,   190,   220,   190,   210}
    }
   ,{{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    }
   }
  ,{{{   270,   160,   270,   210,   270}
    ,{   270,   130,   270,   210,   270}
    ,{   240,   160,   240,    50,   240}
    ,{   240,   100,   240,   180,   240}
    ,{   240,   160,   240,   180,   240}
    }
   ,{{   240,   100,   240,    50,   240}
    ,{   200,    70,   200,    10,   200}
    ,{   240,   100,   240,    50,   240}
    ,{   100,   -30,   100,    40,   100}
    ,{   240,   100,   240,    50,   240}
    }
   ,{{   240,   160,   240,    50,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   240,   160,   240,    50,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   240,   160,   240,    50,   240}
    }
   ,{{   270,   130,   270,   210,   270}
    ,{   270,   130,   270,   210,   270}
    ,{   240,   100,   240,    50,   240}
    ,{   180,   100,   110,   180,   110}
    ,{   240,   100,   240,    50,   240}
    }
   ,{{   240,   160,   240,   180,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   240,   160,   240,    50,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   240,   100,   240,   180,   240}
    }
   }
  ,{{{   300,   220,   300,   220,   150}
    ,{   300,   220,   300,   220,   150}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    }
   ,{{   210,   190,   210,   190,   140}
    ,{   170,   160,   170,   160,   140}
    ,{   210,   190,   210,   190,   120}
    ,{   130,    60,   130,    60,   -10}
    ,{   210,   190,   210,   190,   120}
    }
   ,{{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    }
   ,{{   300,   220,   300,   220,   150}
    ,{   300,   220,   300,   220,   150}
    ,{   210,   190,   210,   190,   120}
    ,{   120,    60,    80,    60,   120}
    ,{   210,   190,   210,   190,   120}
    }
   ,{{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    }
   }
  }
 ,{{{{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   220,   180,   220,   170,   220}
    ,{   220,   180,   220,   180,   220}
    ,{   220,   180,   220,   170,   220}
    }
   ,{{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   210,   170,   210,   170,   210}
    ,{   160,    70,   160,    70,   140}
    ,{   210,   170,   210,   170,   210}
    }
   ,{{   220,   180,   220,   180,   220}
    ,{   220,   180,   220,   180,   220}
    ,{   220,   180,   220,   170,   220}
    ,{   220,   180,   220,   180,   220}
    ,{   220,   180,   220,   170,   220}
    }
   ,{{   230,   170,   230,   170,   210}
    ,{   230,   140,   230,   140,   210}
    ,{   210,   170,   210,   170,   210}
    ,{   130,    60,    60,   130,    70}
    ,{   210,   170,   210,   170,   210}
    }
   ,{{   220,   180,   220,   180,   220}
    ,{   220,   180,   220,   180,   220}
    ,{   220,   180,   220,   170,   220}
    ,{   220,   180,   220,   180,   220}
    ,{   150,   150,   130,    80,   130}
    }
   }
  ,{{{   240,   200,   240,   140,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   220,   180,   220,   140,   220}
    ,{   220,   180,   220,    90,   220}
    ,{   220,   180,   220,   140,   220}
    }
   ,{{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   210,   170,   210,    80,   210}
    ,{   110,    70,   110,   -20,   110}
    ,{   210,   170,   210,    80,   210}
    }
   ,{{   220,   180,   220,   140,   220}
    ,{   220,   180,   220,    90,   220}
    ,{   220,   180,   220,   140,   220}
    ,{   220,   180,   220,    90,   220}
    ,{   220,   180,   220,   140,   220}
    }
   ,{{   210,   170,   210,    80,   210}
    ,{   180,   140,   180,    50,   180}
    ,{   210,   170,   210,    80,   210}
    ,{    60,    20,    60,    60,    60}
    ,{   210,   170,   210,    80,   210}
    }
   ,{{   220,   180,   220,   140,   220}
    ,{   220,   180,   220,    90,   220}
    ,{   220,   180,   220,   140,   220}
    ,{   220,   180,   220,    90,   220}
    ,{   150,   150,   130,     0,   130}
    }
   }
  ,{{{   230,   190,   230,   190,   210}
    ,{   230,   190,   230,   190,   210}
    ,{   200,   170,   200,   170,   190}
    ,{   210,   180,   210,   180,   190}
    ,{   200,   170,   200,   170,   190}
    }
   ,{{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   200,   170,   200,   170,   180}
    ,{   160,    70,   160,    70,   140}
    ,{   200,   170,   200,   170,   180}
    }
   ,{{   210,   180,   210,   180,   190}
    ,{   210,   180,   210,   180,   190}
    ,{   200,   170,   200,   170,   190}
    ,{   210,   180,   210,   180,   190}
    ,{   200,   170,   200,   170,   190}
    }
   ,{{   230,   170,   230,   170,   210}
    ,{   230,   140,   230,   140,   210}
    ,{   200,   170,   200,   170,   180}
    ,{    50,    20,    50,    20,    30}
    ,{   200,   170,   200,   170,   180}
    }
   ,{{   210,   180,   210,   180,   190}
    ,{   210,   180,   210,   180,   190}
    ,{   200,   170,   200,   170,   190}
    ,{   210,   180,   210,   180,   190}
    ,{   110,    80,   110,    80,   100}
    }
   }
  ,{{{   240,   140,   240,   130,   240}
    ,{   240,   100,   240,   120,   240}
    ,{   220,   140,   220,    30,   220}
    ,{   220,    90,   220,   130,   220}
    ,{   220,   140,   220,    70,   220}
    }
   ,{{   240,   100,   240,    50,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   210,    80,   210,    20,   210}
    ,{   110,   -20,   110,    50,   110}
    ,{   210,    80,   210,    20,   210}
    }
   ,{{   220,   140,   220,    30,   220}
    ,{   220,    90,   220,    30,   220}
    ,{   220,   140,   220,    30,   220}
    ,{   220,    90,   220,    30,   220}
    ,{   220,   140,   220,    30,   220}
    }
   ,{{   210,    80,   210,   130,   210}
    ,{   180,    50,   180,   120,   180}
    ,{   210,    80,   210,    20,   210}
    ,{   130,    60,    60,   130,    60}
    ,{   210,    80,   210,    20,   210}
    }
   ,{{   220,   140,   220,    70,   220}
    ,{   220,    90,   220,    30,   220}
    ,{   220,   140,   220,    30,   220}
    ,{   220,    90,   220,    30,   220}
    ,{   130,     0,   130,    70,   130}
    }
   }
  ,{{{   210,   190,   210,   190,   180}
    ,{   210,   190,   210,   190,   180}
    ,{   190,   170,   190,   170,   100}
    ,{   190,   180,   190,   180,   100}
    ,{   190,   170,   190,   170,   100}
    }
   ,{{   210,   190,   210,   190,   180}
    ,{   210,   190,   210,   190,   180}
    ,{   180,   170,   180,   170,    90}
    ,{   140,    70,   140,    70,     0}
    ,{   180,   170,   180,   170,    90}
    }
   ,{{   190,   180,   190,   180,   100}
    ,{   190,   180,   190,   180,   100}
    ,{   190,   170,   190,   170,   100}
    ,{   190,   180,   190,   180,   100}
    ,{   190,   170,   190,   170,   100}
    }
   ,{{   210,   170,   210,   170,    90}
    ,{   210,   140,   210,   140,    60}
    ,{   180,   170,   180,   170,    90}
    ,{    70,    20,    30,    20,    70}
    ,{   180,   170,   180,   170,    90}
    }
   ,{{   190,   180,   190,   180,   100}
    ,{   190,   180,   190,   180,   100}
    ,{   190,   170,   190,   170,   100}
    ,{   190,   180,   190,   180,   100}
    ,{   100,    80,   100,    80,    10}
    }
   }
  }
 ,{{{{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    }
   ,{{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   190,   150,   190,   150,   190}
    ,{   180,    90,   180,    90,   160}
    ,{   190,   150,   190,   150,   190}
    }
   ,{{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    }
   ,{{   190,   150,   190,   150,   190}
    ,{   190,   100,   190,   100,   170}
    ,{   190,   150,   190,   150,   190}
    ,{   150,    80,    80,   150,    90}
    ,{   190,   150,   190,   150,   190}
    }
   ,{{   240,   200,   240,   190,   240}
    ,{   240,   200,   240,   190,   240}
    ,{   210,   170,   210,   160,   210}
    ,{   240,   200,   240,   190,   240}
    ,{   170,   170,   150,   110,   150}
    }
   }
  ,{{{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   160,   240}
    }
   ,{{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   190,   150,   190,    60,   190}
    ,{   130,    90,   130,     0,   130}
    ,{   190,   150,   190,    60,   190}
    }
   ,{{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   160,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   240,   200,   240,   160,   240}
    }
   ,{{   190,   150,   190,    80,   190}
    ,{   140,   100,   140,    10,   140}
    ,{   190,   150,   190,    60,   190}
    ,{    80,    40,    80,    80,    80}
    ,{   190,   150,   190,    60,   190}
    }
   ,{{   240,   200,   240,   130,   240}
    ,{   240,   200,   240,   100,   240}
    ,{   210,   170,   210,   130,   210}
    ,{   240,   200,   240,   100,   240}
    ,{   170,   170,   150,    20,   150}
    }
   }
  ,{{{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    }
   ,{{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   180,   150,   180,   150,   160}
    ,{   180,    90,   180,    90,   160}
    ,{   180,   150,   180,   150,   160}
    }
   ,{{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    }
   ,{{   190,   150,   190,   150,   170}
    ,{   190,   100,   190,   100,   170}
    ,{   180,   150,   180,   150,   160}
    ,{    70,    40,    70,    40,    50}
    ,{   180,   150,   180,   150,   160}
    }
   ,{{   220,   190,   220,   190,   210}
    ,{   220,   190,   220,   190,   210}
    ,{   190,   160,   190,   160,   180}
    ,{   220,   190,   220,   190,   210}
    ,{   140,   110,   140,   110,   120}
    }
   }
  ,{{{   240,   160,   240,   150,   240}
    ,{   240,   100,   240,    80,   240}
    ,{   240,   160,   240,    50,   240}
    ,{   240,   100,   240,   150,   240}
    ,{   240,   160,   240,    90,   240}
    }
   ,{{   240,   100,   240,    70,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   190,    60,   190,     0,   190}
    ,{   130,     0,   130,    70,   130}
    ,{   190,    60,   190,     0,   190}
    }
   ,{{   240,   160,   240,    50,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   240,   160,   240,    50,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   240,   160,   240,    50,   240}
    }
   ,{{   190,    80,   190,   150,   190}
    ,{   140,    10,   140,    80,   140}
    ,{   190,    60,   190,     0,   190}
    ,{   150,    80,    80,   150,    80}
    ,{   190,    60,   190,     0,   190}
    }
   ,{{   240,   130,   240,    90,   240}
    ,{   240,   100,   240,    50,   240}
    ,{   210,   130,   210,    20,   210}
    ,{   240,   100,   240,    50,   240}
    ,{   150,    20,   150,    90,   150}
    }
   }
  ,{{{   210,   190,   210,   190,   180}
    ,{   210,   190,   210,   190,   180}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    }
   ,{{   210,   190,   210,   190,   180}
    ,{   210,   190,   210,   190,   180}
    ,{   160,   150,   160,   150,    70}
    ,{   160,    90,   160,    90,    10}
    ,{   160,   150,   160,   150,    70}
    }
   ,{{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    }
   ,{{   170,   150,   170,   150,    90}
    ,{   170,   100,   170,   100,    20}
    ,{   160,   150,   160,   150,    70}
    ,{    90,    40,    50,    40,    90}
    ,{   160,   150,   160,   150,    70}
    }
   ,{{   210,   190,   210,   190,   120}
    ,{   210,   190,   210,   190,   120}
    ,{   180,   160,   180,   160,    90}
    ,{   210,   190,   210,   190,   120}
    ,{   120,   110,   120,   110,    30}
    }
   }
  }
 ,{{{{   310,   290,   310,   260,   300}
    ,{   310,   270,   310,   260,   300}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   290,   290,   270,   220,   270}
    }
   ,{{   300,   270,   300,   260,   300}
    ,{   300,   270,   300,   260,   300}
    ,{   270,   230,   270,   220,   270}
    ,{   230,   150,   230,   140,   220}
    ,{   270,   230,   270,   220,   270}
    }
   ,{{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    }
   ,{{   310,   230,   310,   220,   300}
    ,{   310,   230,   310,   220,   300}
    ,{   270,   230,   270,   220,   270}
    ,{   210,   130,   140,   210,   150}
    ,{   270,   230,   270,   220,   270}
    }
   ,{{   290,   290,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   270,   230,   270,   220,   270}
    ,{   290,   290,   270,   220,   270}
    }
   }
  ,{{{   300,   290,   300,   190,   300}
    ,{   300,   270,   300,   170,   300}
    ,{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   290,   290,   270,   190,   270}
    }
   ,{{   300,   270,   300,   170,   300}
    ,{   300,   270,   300,   170,   300}
    ,{   270,   230,   270,   130,   270}
    ,{   190,   150,   190,    50,   190}
    ,{   270,   230,   270,   130,   270}
    }
   ,{{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   190,   270}
    }
   ,{{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   140,   100,   140,   130,   140}
    ,{   270,   230,   270,   130,   270}
    }
   ,{{   290,   290,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   270,   230,   270,   190,   270}
    ,{   270,   230,   270,   130,   270}
    ,{   290,   290,   270,   130,   270}
    }
   }
  ,{{{   310,   260,   310,   260,   300}
    ,{   310,   260,   310,   260,   300}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   290,   260,   290,   260,   270}
    ,{   290,   260,   290,   260,   270}
    ,{   250,   220,   250,   220,   240}
    ,{   230,   140,   230,   140,   220}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   310,   220,   310,   220,   300}
    ,{   310,   220,   310,   220,   300}
    ,{   250,   220,   250,   220,   240}
    ,{   120,    90,   120,    90,   110}
    ,{   250,   220,   250,   220,   240}
    }
   ,{{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    ,{   250,   220,   250,   220,   240}
    }
   }
  ,{{{   300,   190,   300,   210,   300}
    ,{   300,   170,   300,   210,   300}
    ,{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,   210,   270}
    ,{   270,   190,   270,   210,   270}
    }
   ,{{   300,   170,   300,   130,   300}
    ,{   300,   170,   300,   110,   300}
    ,{   270,   130,   270,    80,   270}
    ,{   190,    50,   190,   130,   190}
    ,{   270,   130,   270,    80,   270}
    }
   ,{{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   190,   270,    80,   270}
    }
   ,{{   270,   130,   270,   210,   270}
    ,{   270,   130,   270,   210,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   210,   130,   140,   210,   140}
    ,{   270,   130,   270,    80,   270}
    }
   ,{{   270,   190,   270,   210,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   190,   270,    80,   270}
    ,{   270,   130,   270,    80,   270}
    ,{   270,   130,   270,   210,   270}
    }
   }
  ,{{{   300,   260,   300,   260,   240}
    ,{   300,   260,   300,   260,   240}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   270,   260,   270,   260,   240}
    ,{   270,   260,   270,   260,   240}
    ,{   240,   220,   240,   220,   150}
    ,{   220,   140,   220,   140,    70}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   300,   220,   300,   220,   150}
    ,{   300,   220,   300,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   150,    90,   110,    90,   150}
    ,{   240,   220,   240,   220,   150}
    }
   ,{{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    ,{   240,   220,   240,   220,   150}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   220,   220,   190,   150,   150}
    ,{   170,   170,   150,   150,   150}
    ,{   220,   220,   190,   130,   140}
    ,{   170,   170,   150,   150,   150}
    ,{   140,   140,   120,   140,   120}
    }
   ,{{   150,   130,   110,   110,   150}
    ,{   150,   130,   110,   110,   150}
    ,{   130,   130,   110,   100,   110}
    ,{    90,    10,    70,    10,    90}
    ,{   130,   130,   100,   100,   110}
    }
   ,{{   220,   220,   190,   150,   150}
    ,{   150,   150,   150,   150,   150}
    ,{   220,   220,   190,   130,   140}
    ,{   170,   170,   150,   150,   150}
    ,{   140,   140,   120,   120,   120}
    }
   ,{{   140,   130,   100,   100,   140}
    ,{    90,    10,    70,    10,    90}
    ,{   130,   130,   100,   100,   110}
    ,{   140,   -10,    20,    80,   140}
    ,{   130,   130,   100,   100,   110}
    }
   ,{{   170,   170,   170,   150,   150}
    ,{   170,   170,   150,   150,   150}
    ,{   170,   140,   170,   120,   120}
    ,{   170,   170,   150,   150,   150}
    ,{   140,   140,    30,   140,    30}
    }
   }
  ,{{{   220,   220,   190,   140,   140}
    ,{   170,   170,   140,    40,   140}
    ,{   220,   220,   190,    70,   130}
    ,{   170,   170,   140,    30,   140}
    ,{   140,   140,   110,   140,   110}
    }
   ,{{   130,   130,   110,    70,   100}
    ,{   130,   130,   100,    40,   100}
    ,{   130,   130,   110,    70,   100}
    ,{    70,   -20,    70,   -50,    10}
    ,{   130,   130,   100,   -10,   100}
    }
   ,{{   220,   220,   190,    70,   140}
    ,{   140,    60,    50,    30,   140}
    ,{   220,   220,   190,    70,   130}
    ,{   170,   170,   140,    30,   140}
    ,{   140,   140,   110,    50,   110}
    }
   ,{{   130,   130,   100,   -10,   100}
    ,{    10,     0,  -100,   -70,    10}
    ,{   130,   130,   100,   -10,   100}
    ,{   -10,   -10,   -50,   -30,   -50}
    ,{   130,   130,   100,   -10,   100}
    }
   ,{{   170,   170,   140,   140,   140}
    ,{   170,   170,   140,    30,   140}
    ,{   140,   140,   110,    60,   110}
    ,{   170,   170,   140,    30,   140}
    ,{   140,   140,    30,   140,    20}
    }
   }
  ,{{{   150,   150,   150,   150,   150}
    ,{   150,   150,   150,   150,   150}
    ,{   140,   130,   130,   130,   140}
    ,{   150,   150,   150,   150,   150}
    ,{   120,   120,   120,   120,   120}
    }
   ,{{   110,   110,   110,   110,   110}
    ,{   110,   110,   110,   110,   110}
    ,{   110,   100,   100,   100,   110}
    ,{    80,   -40,    70,    10,    80}
    ,{   110,   100,   100,   100,   110}
    }
   ,{{   150,   150,   150,   150,   150}
    ,{   150,   150,   150,   150,   150}
    ,{   140,   130,   130,   130,   140}
    ,{   150,   150,   150,   150,   150}
    ,{   120,   120,   120,   120,   120}
    }
   ,{{   110,   100,   100,   100,   110}
    ,{    80,   -70,   -60,    10,    80}
    ,{   110,   100,   100,   100,   110}
    ,{   -40,   -40,   -40,   -40,   -50}
    ,{   110,   100,   100,   100,   110}
    }
   ,{{   150,   150,   150,   150,   150}
    ,{   150,   150,   150,   150,   150}
    ,{   120,   120,   120,   120,   120}
    ,{   150,   150,   150,   150,   150}
    ,{    30,    30,    30,    30,    30}
    }
   }
  ,{{{   140,    70,   140,    80,   140}
    ,{   140,    10,   140,    10,   140}
    ,{   130,    70,   130,    20,   130}
    ,{   140,   -30,   140,    80,   140}
    ,{   110,    50,   110,    70,   110}
    }
   ,{{   100,   -30,   100,   -30,   100}
    ,{   100,   -30,   100,   -30,   100}
    ,{   100,   -70,   100,   -40,   100}
    ,{    10,  -170,    10,   -30,    10}
    ,{   100,   -70,   100,   -40,   100}
    }
   ,{{   140,    70,   140,    10,   140}
    ,{   140,    10,   140,   -30,   140}
    ,{   130,    70,   130,   -10,   130}
    ,{   140,   -30,   140,    10,   140}
    ,{   110,     0,   110,   -60,   110}
    }
   ,{{   100,   -70,   100,    80,   100}
    ,{    10,  -160,    10,     0,    10}
    ,{   100,   -70,   100,   -40,   100}
    ,{    80,   -90,   -50,    80,   -50}
    ,{   100,   -70,   100,   -40,   100}
    }
   ,{{   140,    50,   140,    70,   140}
    ,{   140,   -30,   140,    10,   140}
    ,{   110,     0,   110,    20,   110}
    ,{   140,   -30,   140,    10,   140}
    ,{    70,    50,    20,    70,    20}
    }
   }
  ,{{{   170,   150,   170,   150,   150}
    ,{   150,   150,   150,   150,   150}
    ,{   170,   130,   170,   130,    30}
    ,{   150,   150,   150,   150,   140}
    ,{   120,   120,   120,   120,    40}
    }
   ,{{   150,   110,   110,   110,   150}
    ,{   150,   110,   110,   110,   150}
    ,{   100,   100,   100,   100,   -20}
    ,{    90,    10,    70,    10,    90}
    ,{   100,   100,   100,   100,    30}
    }
   ,{{   150,   150,   150,   150,    70}
    ,{   150,   150,   150,   150,     0}
    ,{   130,   130,   130,   130,   -10}
    ,{   150,   150,   150,   150,    70}
    ,{   120,   120,   120,   120,    40}
    }
   ,{{   140,   100,   100,   100,   140}
    ,{    90,    10,    70,    10,    90}
    ,{   100,   100,   100,   100,    30}
    ,{   140,   -40,    20,   -40,   140}
    ,{   100,   100,   100,   100,    30}
    }
   ,{{   170,   150,   170,   150,    70}
    ,{   150,   150,   150,   150,    70}
    ,{   170,   120,   170,   120,    20}
    ,{   150,   150,   150,   150,    70}
    ,{    30,    30,    30,    30,   -60}
    }
   }
  }
 ,{{{{   150,   150,   120,   120,   130}
    ,{   150,   150,   120,   120,   130}
    ,{   130,   130,   100,   100,   110}
    ,{   120,   120,    90,    90,   100}
    ,{   120,   120,   100,   100,   100}
    }
   ,{{   150,   150,   120,   120,   130}
    ,{   150,   150,   120,   120,   130}
    ,{   120,   120,   100,   100,   100}
    ,{   -10,   -50,   -20,   -80,   -10}
    ,{   120,   120,   100,   100,   100}
    }
   ,{{   120,   120,   100,   100,   100}
    ,{   120,   120,    90,    90,   100}
    ,{   120,   120,   100,   100,   100}
    ,{   120,   120,    90,    90,   100}
    ,{   120,   120,   100,   100,   100}
    }
   ,{{   120,   120,   100,   100,   100}
    ,{    50,    10,    50,   -10,    50}
    ,{   120,   120,   100,   100,   100}
    ,{    80,   -20,   -40,    80,    10}
    ,{   120,   120,   100,   100,   100}
    }
   ,{{   130,   130,   100,   100,   110}
    ,{   120,   120,    90,    90,   100}
    ,{   130,   130,   100,   100,   110}
    ,{   120,   120,    90,    90,   100}
    ,{   110,   110,    20,    20,    30}
    }
   }
  ,{{{   150,   150,   120,    50,   120}
    ,{   150,   150,   120,    10,   120}
    ,{   130,   130,   100,    50,   100}
    ,{   120,   120,    90,   -20,    90}
    ,{   120,   120,    90,    50,    90}
    }
   ,{{   150,   150,   120,    10,   120}
    ,{   150,   150,   120,    10,   120}
    ,{   120,   120,    90,   -10,    90}
    ,{   -50,   -50,   -80,  -190,   -80}
    ,{   120,   120,    90,   -10,    90}
    }
   ,{{   120,   120,    90,    50,    90}
    ,{   120,   120,    90,   -20,    90}
    ,{   120,   120,    90,    50,    90}
    ,{   120,   120,    90,   -20,    90}
    ,{   120,   120,    90,    50,    90}
    }
   ,{{   120,   120,    90,   -10,    90}
    ,{    10,    10,   -20,  -130,   -20}
    ,{   120,   120,    90,   -10,    90}
    ,{   -20,   -20,   -50,   -20,   -50}
    ,{   120,   120,    90,   -10,    90}
    }
   ,{{   130,   130,   100,    50,   100}
    ,{   120,   120,    90,   -20,    90}
    ,{   130,   130,   100,    50,   100}
    ,{   120,   120,    90,   -20,    90}
    ,{   110,   110,    20,   -90,    20}
    }
   }
  ,{{{   130,   120,   120,   120,   130}
    ,{   130,   120,   120,   120,   130}
    ,{   110,   100,   100,   100,   110}
    ,{   100,    90,    90,    90,   100}
    ,{   100,   100,   100,   100,   100}
    }
   ,{{   130,   120,   120,   120,   130}
    ,{   130,   120,   120,   120,   130}
    ,{   100,   100,   100,   100,   100}
    ,{   -10,   -80,   -20,   -80,   -10}
    ,{   100,   100,   100,   100,   100}
    }
   ,{{   100,   100,   100,   100,   100}
    ,{   100,    90,    90,    90,   100}
    ,{   100,   100,   100,   100,   100}
    ,{   100,    90,    90,    90,   100}
    ,{   100,   100,   100,   100,   100}
    }
   ,{{   100,   100,   100,   100,   100}
    ,{    50,   -10,    50,   -10,    50}
    ,{   100,   100,   100,   100,   100}
    ,{   -40,   -40,   -40,   -40,   -40}
    ,{   100,   100,   100,   100,   100}
    }
   ,{{   110,   100,   100,   100,   110}
    ,{   100,    90,    90,    90,   100}
    ,{   110,   100,   100,   100,   110}
    ,{   100,    90,    90,    90,   100}
    ,{    30,    20,    20,    20,    30}
    }
   }
  ,{{{   120,   -10,   120,    80,   120}
    ,{   120,   -50,   120,   -20,   120}
    ,{   100,   -10,   100,   -40,   100}
    ,{    90,   -80,    90,    80,    90}
    ,{    90,   -20,    90,    10,    90}
    }
   ,{{   120,   -50,   120,   -20,   120}
    ,{   120,   -50,   120,   -20,   120}
    ,{    90,   -80,    90,   -40,    90}
    ,{   -80,  -260,   -80,   -90,   -80}
    ,{    90,   -80,    90,   -40,    90}
    }
   ,{{    90,   -20,    90,   -40,    90}
    ,{    90,   -80,    90,   -50,    90}
    ,{    90,   -20,    90,   -40,    90}
    ,{    90,   -80,    90,   -50,    90}
    ,{    90,   -20,    90,   -40,    90}
    }
   ,{{    90,   -80,    90,    80,    90}
    ,{   -20,  -190,   -20,   -20,   -20}
    ,{    90,   -80,    90,   -40,    90}
    ,{    80,   -90,   -50,    80,   -50}
    ,{    90,   -80,    90,   -40,    90}
    }
   ,{{   100,   -10,   100,    10,   100}
    ,{    90,   -80,    90,   -50,    90}
    ,{   100,   -10,   100,   -40,   100}
    ,{    90,   -80,    90,   -50,    90}
    ,{    20,  -150,    20,    10,    20}
    }
   }
  ,{{{   120,   120,   120,   120,   110}
    ,{   120,   120,   120,   120,   110}
    ,{   100,   100,   100,   100,    30}
    ,{    90,    90,    90,    90,    20}
    ,{   100,   100,   100,   100,    20}
    }
   ,{{   120,   120,   120,   120,   110}
    ,{   120,   120,   120,   120,   110}
    ,{   100,   100,   100,   100,    20}
    ,{   -20,   -80,   -20,   -80,  -150}
    ,{   100,   100,   100,   100,    20}
    }
   ,{{   100,   100,   100,   100,    20}
    ,{    90,    90,    90,    90,    20}
    ,{   100,   100,   100,   100,    20}
    ,{    90,    90,    90,    90,    20}
    ,{   100,   100,   100,   100,    20}
    }
   ,{{   100,   100,   100,   100,    20}
    ,{    50,   -10,    50,   -10,   -90}
    ,{   100,   100,   100,   100,    20}
    ,{    10,   -40,   -40,   -40,    10}
    ,{   100,   100,   100,   100,    20}
    }
   ,{{   100,   100,   100,   100,    30}
    ,{    90,    90,    90,    90,    20}
    ,{   100,   100,   100,   100,    30}
    ,{    90,    90,    90,    90,    20}
    ,{    20,    20,    20,    20,   -50}
    }
   }
  }
 ,{{{{   300,   300,   250,   250,   260}
    ,{   280,   280,   250,   250,   260}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   300,   300,   220,   220,   220}
    }
   ,{{   280,   280,   250,   250,   260}
    ,{   280,   280,   250,   250,   260}
    ,{   240,   240,   220,   220,   220}
    ,{   200,   160,   200,   140,   200}
    ,{   240,   240,   220,   220,   220}
    }
   ,{{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    }
   ,{{   240,   240,   240,   220,   240}
    ,{   240,   200,   240,   180,   240}
    ,{   240,   240,   220,   220,   220}
    ,{   210,   110,    90,   210,   140}
    ,{   240,   240,   220,   220,   220}
    }
   ,{{   300,   300,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   300,   300,   220,   220,   220}
    }
   }
  ,{{{   300,   300,   250,   160,   250}
    ,{   280,   280,   250,   140,   250}
    ,{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   300,   300,   210,   160,   210}
    }
   ,{{   280,   280,   250,   140,   250}
    ,{   280,   280,   250,   140,   250}
    ,{   240,   240,   210,   100,   210}
    ,{   160,   160,   130,    20,   130}
    ,{   240,   240,   210,   100,   210}
    }
   ,{{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   160,   210}
    }
   ,{{   240,   240,   210,   100,   210}
    ,{   200,   200,   170,    60,   170}
    ,{   240,   240,   210,   100,   210}
    ,{   110,   110,    80,   100,    80}
    ,{   240,   240,   210,   100,   210}
    }
   ,{{   300,   300,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   300,   300,   210,   100,   210}
    }
   }
  ,{{{   260,   250,   250,   250,   260}
    ,{   260,   250,   250,   250,   260}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   260,   250,   250,   250,   260}
    ,{   260,   250,   250,   250,   260}
    ,{   220,   220,   220,   220,   220}
    ,{   200,   140,   200,   140,   200}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   240,   220,   240,   220,   240}
    ,{   240,   180,   240,   180,   240}
    ,{   220,   220,   220,   220,   220}
    ,{    90,    90,    90,    90,    90}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    }
   }
  ,{{{   250,   100,   250,   210,   250}
    ,{   250,    70,   250,   170,   250}
    ,{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,   210,   210}
    ,{   210,   100,   210,   210,   210}
    }
   ,{{   250,    70,   250,   130,   250}
    ,{   250,    70,   250,   110,   250}
    ,{   210,    40,   210,    80,   210}
    ,{   130,   -40,   130,   130,   130}
    ,{   210,    40,   210,    80,   210}
    }
   ,{{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,   100,   210,    80,   210}
    }
   ,{{   210,    40,   210,   210,   210}
    ,{   170,     0,   170,   170,   170}
    ,{   210,    40,   210,    80,   210}
    ,{   210,    40,    80,   210,    80}
    ,{   210,    40,   210,    80,   210}
    }
   ,{{   210,   100,   210,   210,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,    40,   210,   210,   210}
    }
   }
  ,{{{   250,   250,   250,   250,   240}
    ,{   250,   250,   250,   250,   240}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   250,   250,   250,   250,   240}
    ,{   250,   250,   250,   250,   240}
    ,{   220,   220,   220,   220,   140}
    ,{   200,   140,   200,   140,    60}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   240,   220,   240,   220,   140}
    ,{   240,   180,   240,   180,   100}
    ,{   220,   220,   220,   220,   140}
    ,{   140,    90,    90,    90,   140}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    }
   }
  }
 ,{{{{   280,   270,   280,   220,   280}
    ,{   280,   240,   280,   220,   280}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   270,   270,   190,   190,   190}
    }
   ,{{   210,   210,   190,   190,   190}
    ,{   190,   190,   150,   150,   160}
    ,{   210,   210,   190,   190,   190}
    ,{   120,    80,   110,    50,   120}
    ,{   210,   210,   190,   190,   190}
    }
   ,{{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    }
   ,{{   280,   240,   280,   220,   280}
    ,{   280,   240,   280,   220,   280}
    ,{   210,   210,   190,   190,   190}
    ,{   180,    80,    60,   180,   110}
    ,{   210,   210,   190,   190,   190}
    }
   ,{{   270,   270,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   270,   270,   190,   190,   190}
    }
   }
  ,{{{   270,   270,   210,   130,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   270,   270,   180,   130,   180}
    }
   ,{{   210,   210,   180,    70,   180}
    ,{   190,   190,   150,    40,   150}
    ,{   210,   210,   180,    70,   180}
    ,{    80,    80,    50,   -60,    50}
    ,{   210,   210,   180,    70,   180}
    }
   ,{{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,   130,   180}
    }
   ,{{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   210,   210,   180,    70,   180}
    ,{    80,    80,    50,    70,    50}
    ,{   210,   210,   180,    70,   180}
    }
   ,{{   270,   270,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   270,   270,   180,    70,   180}
    }
   }
  ,{{{   280,   220,   280,   220,   280}
    ,{   280,   220,   280,   220,   280}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   160,   150,   150,   150,   160}
    ,{   190,   190,   190,   190,   190}
    ,{   120,    50,   110,    50,   120}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   280,   220,   280,   220,   280}
    ,{   280,   220,   280,   220,   280}
    ,{   190,   190,   190,   190,   190}
    ,{    60,    60,    60,    60,    60}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    }
   }
  ,{{{   210,    70,   210,   210,   210}
    ,{   210,    40,   210,   210,   210}
    ,{   180,    70,   180,    50,   180}
    ,{   180,    10,   180,   180,   180}
    ,{   180,    70,   180,   180,   180}
    }
   ,{{   180,    10,   180,    50,   180}
    ,{   150,   -20,   150,    10,   150}
    ,{   180,    10,   180,    50,   180}
    ,{    50,  -120,    50,    40,    50}
    ,{   180,    10,   180,    50,   180}
    }
   ,{{   180,    70,   180,    50,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   180,    70,   180,    50,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   180,    70,   180,    50,   180}
    }
   ,{{   210,    40,   210,   210,   210}
    ,{   210,    40,   210,   210,   210}
    ,{   180,    10,   180,    50,   180}
    ,{   180,    10,    50,   180,    50}
    ,{   180,    10,   180,    50,   180}
    }
   ,{{   180,    70,   180,   180,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   180,    70,   180,    50,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   180,    10,   180,   180,   180}
    }
   }
  ,{{{   280,   220,   280,   220,   140}
    ,{   280,   220,   280,   220,   140}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    }
   ,{{   190,   190,   190,   190,   140}
    ,{   150,   150,   150,   150,   140}
    ,{   190,   190,   190,   190,   110}
    ,{   110,    50,   110,    50,   -20}
    ,{   190,   190,   190,   190,   110}
    }
   ,{{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    }
   ,{{   280,   220,   280,   220,   140}
    ,{   280,   220,   280,   220,   140}
    ,{   190,   190,   190,   190,   110}
    ,{   110,    60,    60,    60,   110}
    ,{   190,   190,   190,   190,   110}
    }
   ,{{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    }
   }
  }
 ,{{{{   210,   210,   190,   190,   200}
    ,{   210,   210,   190,   190,   200}
    ,{   190,   190,   170,   170,   170}
    ,{   200,   200,   170,   170,   180}
    ,{   190,   190,   170,   170,   170}
    }
   ,{{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   190,   190,   160,   160,   170}
    ,{   130,    90,   120,    60,   130}
    ,{   190,   190,   160,   160,   170}
    }
   ,{{   200,   200,   170,   170,   180}
    ,{   200,   200,   170,   170,   180}
    ,{   190,   190,   170,   170,   170}
    ,{   200,   200,   170,   170,   180}
    ,{   190,   190,   170,   170,   170}
    }
   ,{{   200,   190,   190,   160,   200}
    ,{   200,   160,   190,   130,   200}
    ,{   190,   190,   160,   160,   170}
    ,{   130,    40,    10,   130,    70}
    ,{   190,   190,   160,   160,   170}
    }
   ,{{   200,   200,   170,   170,   180}
    ,{   200,   200,   170,   170,   180}
    ,{   190,   190,   170,   170,   170}
    ,{   200,   200,   170,   170,   180}
    ,{   160,   160,    80,    80,    80}
    }
   }
  ,{{{   210,   210,   180,   110,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   190,   190,   160,   110,   160}
    ,{   200,   200,   170,    60,   170}
    ,{   190,   190,   160,   110,   160}
    }
   ,{{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   190,   190,   160,    50,   160}
    ,{    90,    90,    60,   -50,    60}
    ,{   190,   190,   160,    50,   160}
    }
   ,{{   200,   200,   170,   110,   170}
    ,{   200,   200,   170,    60,   170}
    ,{   190,   190,   160,   110,   160}
    ,{   200,   200,   170,    60,   170}
    ,{   190,   190,   160,   110,   160}
    }
   ,{{   190,   190,   160,    50,   160}
    ,{   160,   160,   130,    20,   130}
    ,{   190,   190,   160,    50,   160}
    ,{    40,    40,    10,    30,    10}
    ,{   190,   190,   160,    50,   160}
    }
   ,{{   200,   200,   170,   110,   170}
    ,{   200,   200,   170,    60,   170}
    ,{   190,   190,   160,   110,   160}
    ,{   200,   200,   170,    60,   170}
    ,{   160,   160,    70,   -30,    70}
    }
   }
  ,{{{   200,   190,   190,   190,   200}
    ,{   200,   190,   190,   190,   200}
    ,{   170,   170,   170,   170,   170}
    ,{   180,   170,   170,   170,   180}
    ,{   170,   170,   170,   170,   170}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   170,   160,   160,   160,   170}
    ,{   130,    60,   120,    60,   130}
    ,{   170,   160,   160,   160,   170}
    }
   ,{{   180,   170,   170,   170,   180}
    ,{   180,   170,   170,   170,   180}
    ,{   170,   170,   170,   170,   170}
    ,{   180,   170,   170,   170,   180}
    ,{   170,   170,   170,   170,   170}
    }
   ,{{   200,   160,   190,   160,   200}
    ,{   200,   130,   190,   130,   200}
    ,{   170,   160,   160,   160,   170}
    ,{    20,    10,    10,    10,    20}
    ,{   170,   160,   160,   160,   170}
    }
   ,{{   180,   170,   170,   170,   180}
    ,{   180,   170,   170,   170,   180}
    ,{   170,   170,   170,   170,   170}
    ,{   180,   170,   170,   170,   180}
    ,{    80,    80,    80,    80,    80}
    }
   }
  ,{{{   180,    50,   180,   130,   180}
    ,{   180,    10,   180,   120,   180}
    ,{   160,    50,   160,    30,   160}
    ,{   170,     0,   170,   130,   170}
    ,{   160,    50,   160,    70,   160}
    }
   ,{{   180,    10,   180,    50,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   160,   -10,   160,    20,   160}
    ,{    60,  -110,    60,    50,    60}
    ,{   160,   -10,   160,    20,   160}
    }
   ,{{   170,    50,   170,    30,   170}
    ,{   170,     0,   170,    30,   170}
    ,{   160,    50,   160,    30,   160}
    ,{   170,     0,   170,    30,   170}
    ,{   160,    50,   160,    30,   160}
    }
   ,{{   160,   -10,   160,   130,   160}
    ,{   130,   -40,   130,   120,   130}
    ,{   160,   -10,   160,    20,   160}
    ,{   130,   -30,    10,   130,    10}
    ,{   160,   -10,   160,    20,   160}
    }
   ,{{   170,    50,   170,    70,   170}
    ,{   170,     0,   170,    30,   170}
    ,{   160,    50,   160,    30,   160}
    ,{   170,     0,   170,    30,   170}
    ,{    70,  -100,    70,    70,    70}
    }
   }
  ,{{{   190,   190,   190,   190,   170}
    ,{   190,   190,   190,   190,   170}
    ,{   170,   170,   170,   170,    90}
    ,{   170,   170,   170,   170,   100}
    ,{   170,   170,   170,   170,    90}
    }
   ,{{   190,   190,   190,   190,   170}
    ,{   190,   190,   190,   190,   170}
    ,{   160,   160,   160,   160,    90}
    ,{   120,    60,   120,    60,   -10}
    ,{   160,   160,   160,   160,    90}
    }
   ,{{   170,   170,   170,   170,   100}
    ,{   170,   170,   170,   170,   100}
    ,{   170,   170,   170,   170,    90}
    ,{   170,   170,   170,   170,   100}
    ,{   170,   170,   170,   170,    90}
    }
   ,{{   190,   160,   190,   160,    90}
    ,{   190,   130,   190,   130,    60}
    ,{   160,   160,   160,   160,    90}
    ,{    70,    10,    10,    10,    70}
    ,{   160,   160,   160,   160,    90}
    }
   ,{{   170,   170,   170,   170,   100}
    ,{   170,   170,   170,   170,   100}
    ,{   170,   170,   170,   170,    90}
    ,{   170,   170,   170,   170,   100}
    ,{    80,    80,    80,    80,     0}
    }
   }
  }
 ,{{{{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    }
   ,{{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   170,   170,   140,   140,   150}
    ,{   150,   110,   140,    80,   150}
    ,{   170,   170,   140,   140,   150}
    }
   ,{{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    }
   ,{{   170,   170,   150,   150,   160}
    ,{   160,   120,   150,    90,   160}
    ,{   170,   170,   140,   140,   150}
    ,{   150,    60,    30,   150,    90}
    ,{   170,   170,   140,   140,   150}
    }
   ,{{   210,   210,   190,   190,   190}
    ,{   210,   210,   190,   190,   190}
    ,{   180,   180,   160,   160,   160}
    ,{   210,   210,   190,   190,   190}
    ,{   190,   190,   100,   100,   110}
    }
   }
  ,{{{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,   130,   180}
    }
   ,{{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   170,   170,   140,    30,   140}
    ,{   110,   110,    80,   -30,    80}
    ,{   170,   170,   140,    30,   140}
    }
   ,{{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,   130,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   210,   210,   180,   130,   180}
    }
   ,{{   170,   170,   140,    50,   140}
    ,{   120,   120,    90,   -20,    90}
    ,{   170,   170,   140,    30,   140}
    ,{    60,    60,    30,    50,    30}
    ,{   170,   170,   140,    30,   140}
    }
   ,{{   210,   210,   180,   100,   180}
    ,{   210,   210,   180,    70,   180}
    ,{   180,   180,   150,   100,   150}
    ,{   210,   210,   180,    70,   180}
    ,{   190,   190,   100,   -10,   100}
    }
   }
  ,{{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   150,   140,   140,   140,   150}
    ,{   150,    80,   140,    80,   150}
    ,{   150,   140,   140,   140,   150}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   160,   140,   150,   140,   160}
    ,{   160,    90,   150,    90,   160}
    ,{   150,   140,   140,   140,   150}
    ,{    40,    30,    30,    30,    40}
    ,{   150,   140,   140,   140,   150}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   160,   160,   160,   160,   160}
    ,{   190,   190,   190,   190,   190}
    ,{   110,   100,   100,   100,   110}
    }
   }
  ,{{{   180,    70,   180,   150,   180}
    ,{   180,    10,   180,    80,   180}
    ,{   180,    70,   180,    50,   180}
    ,{   180,    10,   180,   150,   180}
    ,{   180,    70,   180,    90,   180}
    }
   ,{{   180,    10,   180,    70,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   140,   -30,   140,     0,   140}
    ,{    80,   -90,    80,    70,    80}
    ,{   140,   -30,   140,     0,   140}
    }
   ,{{   180,    70,   180,    50,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   180,    70,   180,    50,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   180,    70,   180,    50,   180}
    }
   ,{{   150,   -10,   140,   150,   140}
    ,{    90,   -80,    90,    80,    90}
    ,{   140,   -30,   140,     0,   140}
    ,{   150,   -10,    30,   150,    30}
    ,{   140,   -30,   140,     0,   140}
    }
   ,{{   180,    40,   180,    90,   180}
    ,{   180,    10,   180,    50,   180}
    ,{   150,    40,   150,    20,   150}
    ,{   180,    10,   180,    50,   180}
    ,{   100,   -70,   100,    90,   100}
    }
   }
  ,{{{   190,   190,   190,   190,   170}
    ,{   190,   190,   190,   190,   170}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    }
   ,{{   190,   190,   190,   190,   170}
    ,{   190,   190,   190,   190,   170}
    ,{   140,   140,   140,   140,    70}
    ,{   140,    80,   140,    80,    10}
    ,{   140,   140,   140,   140,    70}
    }
   ,{{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    }
   ,{{   150,   140,   150,   140,    90}
    ,{   150,    90,   150,    90,    20}
    ,{   140,   140,   140,   140,    70}
    ,{    90,    30,    30,    30,    90}
    ,{   140,   140,   140,   140,    70}
    }
   ,{{   190,   190,   190,   190,   110}
    ,{   190,   190,   190,   190,   110}
    ,{   160,   160,   160,   160,    80}
    ,{   190,   190,   190,   190,   110}
    ,{   100,   100,   100,   100,    30}
    }
   }
  }
 ,{{{{   300,   300,   280,   250,   280}
    ,{   280,   280,   280,   250,   280}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   300,   300,   220,   220,   220}
    }
   ,{{   280,   280,   250,   250,   260}
    ,{   280,   280,   250,   250,   260}
    ,{   240,   240,   220,   220,   220}
    ,{   200,   160,   200,   140,   200}
    ,{   240,   240,   220,   220,   220}
    }
   ,{{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    }
   ,{{   280,   240,   280,   220,   280}
    ,{   280,   240,   280,   220,   280}
    ,{   240,   240,   220,   220,   220}
    ,{   210,   110,    90,   210,   140}
    ,{   240,   240,   220,   220,   220}
    }
   ,{{   300,   300,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   240,   240,   220,   220,   220}
    ,{   300,   300,   220,   220,   220}
    }
   }
  ,{{{   300,   300,   250,   160,   250}
    ,{   280,   280,   250,   140,   250}
    ,{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   300,   300,   210,   160,   210}
    }
   ,{{   280,   280,   250,   140,   250}
    ,{   280,   280,   250,   140,   250}
    ,{   240,   240,   210,   100,   210}
    ,{   160,   160,   130,    20,   130}
    ,{   240,   240,   210,   100,   210}
    }
   ,{{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   160,   210}
    }
   ,{{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   110,   110,    80,   100,    80}
    ,{   240,   240,   210,   100,   210}
    }
   ,{{   300,   300,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   240,   240,   210,   160,   210}
    ,{   240,   240,   210,   100,   210}
    ,{   300,   300,   210,   140,   210}
    }
   }
  ,{{{   280,   250,   280,   250,   280}
    ,{   280,   250,   280,   250,   280}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   260,   250,   250,   250,   260}
    ,{   260,   250,   250,   250,   260}
    ,{   220,   220,   220,   220,   220}
    ,{   200,   140,   200,   140,   200}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   280,   220,   280,   220,   280}
    ,{   280,   220,   280,   220,   280}
    ,{   220,   220,   220,   220,   220}
    ,{    90,    90,    90,    90,    90}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    ,{   220,   220,   220,   220,   220}
    }
   }
  ,{{{   250,   100,   250,   210,   250}
    ,{   250,    70,   250,   210,   250}
    ,{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,   210,   210}
    ,{   210,   100,   210,   210,   210}
    }
   ,{{   250,    70,   250,   130,   250}
    ,{   250,    70,   250,   110,   250}
    ,{   210,    40,   210,    80,   210}
    ,{   130,   -40,   130,   130,   130}
    ,{   210,    40,   210,    80,   210}
    }
   ,{{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,   100,   210,    80,   210}
    }
   ,{{   210,    40,   210,   210,   210}
    ,{   210,    40,   210,   210,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,    40,    80,   210,    80}
    ,{   210,    40,   210,    80,   210}
    }
   ,{{   210,   100,   210,   210,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,   100,   210,    80,   210}
    ,{   210,    40,   210,    80,   210}
    ,{   210,    50,   210,   210,   210}
    }
   }
  ,{{{   280,   250,   280,   250,   240}
    ,{   280,   250,   280,   250,   240}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   250,   250,   250,   250,   240}
    ,{   250,   250,   250,   250,   240}
    ,{   220,   220,   220,   220,   140}
    ,{   200,   140,   200,   140,    90}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   280,   220,   280,   220,   140}
    ,{   280,   220,   280,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   140,    90,    90,    90,   140}
    ,{   220,   220,   220,   220,   140}
    }
   ,{{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    ,{   220,   220,   220,   220,   140}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   300,   300,   270,   270,   290}
    ,{   300,   300,   270,   270,   290}
    ,{   290,   290,   250,   270,   250}
    ,{   300,   300,   270,   270,   270}
    ,{   270,   270,   240,   260,   240}
    }
   ,{{   290,   270,   230,   230,   290}
    ,{   290,   270,   230,   230,   290}
    ,{   260,   260,   220,   220,   220}
    ,{   190,   170,   190,   130,   190}
    ,{   260,   260,   220,   220,   220}
    }
   ,{{   300,   300,   270,   270,   270}
    ,{   300,   300,   270,   270,   270}
    ,{   290,   290,   250,   270,   250}
    ,{   300,   300,   270,   270,   270}
    ,{   270,   270,   240,   260,   240}
    }
   ,{{   260,   260,   220,   220,   220}
    ,{   190,   170,   190,   130,   190}
    ,{   260,   260,   220,   220,   220}
    ,{   210,   130,    80,   210,   210}
    ,{   260,   260,   220,   220,   220}
    }
   ,{{   300,   300,   270,   270,   270}
    ,{   300,   300,   270,   270,   270}
    ,{   270,   270,   240,   260,   240}
    ,{   300,   300,   270,   270,   270}
    ,{   240,   240,   150,   150,   150}
    }
   }
  ,{{{   300,   300,   270,   270,   270}
    ,{   300,   300,   270,   230,   270}
    ,{   290,   290,   250,   270,   250}
    ,{   300,   300,   270,   230,   270}
    ,{   270,   270,   240,   260,   240}
    }
   ,{{   270,   270,   230,   190,   230}
    ,{   270,   270,   230,   190,   230}
    ,{   260,   260,   220,   180,   220}
    ,{   170,   170,   130,    90,   130}
    ,{   260,   260,   220,   180,   220}
    }
   ,{{   300,   300,   270,   270,   270}
    ,{   300,   300,   270,   230,   270}
    ,{   290,   290,   250,   270,   250}
    ,{   300,   300,   270,   230,   270}
    ,{   270,   270,   240,   260,   240}
    }
   ,{{   260,   260,   220,   180,   220}
    ,{   170,   170,   130,    90,   130}
    ,{   260,   260,   220,   180,   220}
    ,{   170,   110,    80,   170,    80}
    ,{   260,   260,   220,   180,   220}
    }
   ,{{   300,   300,   270,   260,   270}
    ,{   300,   300,   270,   230,   270}
    ,{   270,   270,   240,   260,   240}
    ,{   300,   300,   270,   230,   270}
    ,{   240,   240,   150,   110,   150}
    }
   }
  ,{{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   230,   230,   230,   230,   230}
    ,{   230,   230,   230,   230,   230}
    ,{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   190}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   190}
    ,{   220,   220,   220,   220,   220}
    ,{    80,    80,    80,    80,    80}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    ,{   270,   270,   270,   270,   270}
    ,{   150,   150,   150,   150,   150}
    }
   }
  ,{{{   270,   230,   270,   210,   270}
    ,{   270,   190,   270,   140,   270}
    ,{   250,   230,   250,   120,   250}
    ,{   270,   190,   270,   210,   270}
    ,{   240,   220,   240,   150,   240}
    }
   ,{{   230,   150,   230,   130,   230}
    ,{   230,   150,   230,   100,   230}
    ,{   220,   140,   220,    90,   220}
    ,{   130,    50,   130,   130,   130}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   270,   230,   270,   140,   270}
    ,{   270,   190,   270,   140,   270}
    ,{   250,   230,   250,   120,   250}
    ,{   270,   190,   270,   140,   270}
    ,{   240,   220,   240,   110,   240}
    }
   ,{{   220,   140,   220,   210,   220}
    ,{   130,    50,   130,   130,   130}
    ,{   220,   140,   220,    90,   220}
    ,{   210,   130,    80,   210,    80}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   270,   220,   270,   150,   270}
    ,{   270,   190,   270,   140,   270}
    ,{   240,   220,   240,   110,   240}
    ,{   270,   190,   270,   140,   270}
    ,{   150,    70,   150,   150,   150}
    }
   }
  ,{{{   290,   270,   270,   270,   290}
    ,{   290,   270,   270,   270,   290}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   290,   230,   230,   230,   290}
    ,{   290,   230,   230,   230,   290}
    ,{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   130}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   130}
    ,{   220,   220,   220,   220,   220}
    ,{   210,    80,    80,    80,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    ,{   270,   270,   270,   270,   270}
    ,{   150,   150,   150,   150,   150}
    }
   }
  }
 ,{{{{   300,   280,   240,   240,   300}
    ,{   300,   280,   240,   240,   300}
    ,{   260,   260,   220,   240,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   250,   250,   220,   240,   220}
    }
   ,{{   300,   280,   240,   240,   300}
    ,{   300,   280,   240,   240,   300}
    ,{   250,   250,   220,   220,   220}
    ,{   100,    70,   100,    40,   100}
    ,{   250,   250,   220,   220,   220}
    }
   ,{{   250,   250,   220,   240,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   250,   250,   220,   240,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   250,   250,   220,   240,   220}
    }
   ,{{   250,   250,   220,   220,   220}
    ,{   160,   140,   160,   100,   160}
    ,{   250,   250,   220,   220,   220}
    ,{   210,   130,    80,   210,   210}
    ,{   250,   250,   220,   220,   220}
    }
   ,{{   260,   260,   220,   240,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   260,   260,   220,   240,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   240,   240,   140,   140,   140}
    }
   }
  ,{{{   280,   280,   240,   240,   240}
    ,{   280,   280,   240,   200,   240}
    ,{   260,   260,   220,   240,   220}
    ,{   250,   250,   210,   170,   210}
    ,{   250,   250,   220,   240,   220}
    }
   ,{{   280,   280,   240,   200,   240}
    ,{   280,   280,   240,   200,   240}
    ,{   250,   250,   220,   180,   220}
    ,{    70,    70,    40,     0,    40}
    ,{   250,   250,   220,   180,   220}
    }
   ,{{   250,   250,   220,   240,   220}
    ,{   250,   250,   210,   170,   210}
    ,{   250,   250,   220,   240,   220}
    ,{   250,   250,   210,   170,   210}
    ,{   250,   250,   220,   240,   220}
    }
   ,{{   250,   250,   220,   180,   220}
    ,{   140,   140,   100,    60,   100}
    ,{   250,   250,   220,   180,   220}
    ,{   170,   110,    80,   170,    80}
    ,{   250,   250,   220,   180,   220}
    }
   ,{{   260,   260,   220,   240,   220}
    ,{   250,   250,   210,   170,   210}
    ,{   260,   260,   220,   240,   220}
    ,{   250,   250,   210,   170,   210}
    ,{   240,   240,   140,   100,   140}
    }
   }
  ,{{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   220,   220,   220,   220,   220}
    ,{   100,    40,   100,    40,   100}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   160,   100,   160,   100,   160}
    ,{   220,   220,   220,   220,   220}
    ,{    80,    80,    80,    80,    80}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   140,   140,   140,   140,   140}
    }
   }
  ,{{{   240,   200,   240,   210,   240}
    ,{   240,   160,   240,   110,   240}
    ,{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,   210,   210}
    ,{   220,   200,   220,   140,   220}
    }
   ,{{   240,   160,   240,   110,   240}
    ,{   240,   160,   240,   110,   240}
    ,{   220,   140,   220,    90,   220}
    ,{    40,   -40,    40,    40,    40}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   220,   200,   220,    90,   220}
    }
   ,{{   220,   140,   220,   210,   220}
    ,{   100,    20,   100,   100,   100}
    ,{   220,   140,   220,    90,   220}
    ,{   210,   130,    80,   210,    80}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   220,   200,   220,   140,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   140,    60,   140,   140,   140}
    }
   }
  ,{{{   300,   240,   240,   240,   300}
    ,{   300,   240,   240,   240,   300}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   300,   240,   240,   240,   300}
    ,{   300,   240,   240,   240,   300}
    ,{   220,   220,   220,   220,   220}
    ,{   100,    40,   100,    40,    40}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   160,   100,   160,   100,   100}
    ,{   220,   220,   220,   220,   220}
    ,{   210,    80,    80,    80,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   140,   140,   140,   140,   140}
    }
   }
  }
 ,{{{{   430,   430,   370,   370,   430}
    ,{   430,   410,   370,   370,   430}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   360,   340}
    }
   ,{{   430,   410,   370,   370,   430}
    ,{   430,   410,   370,   370,   430}
    ,{   370,   370,   340,   340,   340}
    ,{   320,   290,   320,   260,   320}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   360,   340}
    }
   ,{{   370,   370,   360,   340,   360}
    ,{   360,   330,   360,   300,   360}
    ,{   370,   370,   340,   340,   340}
    ,{   340,   260,   210,   340,   340}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   430,   430,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   340,   340}
    }
   }
  ,{{{   430,   430,   370,   360,   370}
    ,{   410,   410,   370,   330,   370}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   430,   430,   340,   360,   340}
    }
   ,{{   410,   410,   370,   330,   370}
    ,{   410,   410,   370,   330,   370}
    ,{   370,   370,   340,   300,   340}
    ,{   290,   290,   260,   220,   260}
    ,{   370,   370,   340,   300,   340}
    }
   ,{{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   360,   340}
    }
   ,{{   370,   370,   340,   300,   340}
    ,{   330,   330,   300,   260,   300}
    ,{   370,   370,   340,   300,   340}
    ,{   300,   240,   210,   300,   210}
    ,{   370,   370,   340,   300,   340}
    }
   ,{{   430,   430,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   430,   430,   340,   300,   340}
    }
   }
  ,{{{   370,   370,   370,   370,   370}
    ,{   370,   370,   370,   370,   370}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   370,   370,   370,   370,   370}
    ,{   370,   370,   370,   370,   370}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   320}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   360,   340,   360,   340,   360}
    ,{   360,   300,   360,   300,   360}
    ,{   340,   340,   340,   340,   340}
    ,{   210,   210,   210,   210,   210}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  ,{{{   370,   320,   370,   340,   370}
    ,{   370,   290,   370,   300,   370}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   340,   320,   340,   340,   340}
    }
   ,{{   370,   290,   370,   260,   370}
    ,{   370,   290,   370,   240,   370}
    ,{   340,   260,   340,   210,   340}
    ,{   260,   180,   260,   260,   260}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    }
   ,{{   340,   260,   340,   340,   340}
    ,{   300,   220,   300,   300,   300}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   210,   340,   210}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   340,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    }
   }
  ,{{{   430,   370,   370,   370,   430}
    ,{   430,   370,   370,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   430,   370,   370,   370,   430}
    ,{   430,   370,   370,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   260}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   360,   340,   360,   340,   340}
    ,{   360,   300,   360,   300,   300}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   210,   210,   210,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  }
 ,{{{{   400,   400,   400,   360,   400}
    ,{   400,   370,   400,   360,   400}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   400,   400,   310,   330,   310}
    }
   ,{{   360,   360,   310,   360,   330}
    ,{   360,   360,   270,   360,   330}
    ,{   340,   340,   310,   310,   310}
    ,{   230,   220,   230,   170,   230}
    ,{   340,   340,   310,   310,   310}
    }
   ,{{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   330,   310}
    }
   ,{{   400,   370,   400,   340,   400}
    ,{   400,   370,   400,   340,   400}
    ,{   340,   340,   310,   310,   310}
    ,{   310,   230,   180,   310,   310}
    ,{   340,   340,   310,   310,   310}
    }
   ,{{   400,   400,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   400,   400,   310,   310,   310}
    }
   }
  ,{{{   400,   400,   340,   360,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   400,   400,   310,   330,   310}
    }
   ,{{   360,   360,   310,   360,   310}
    ,{   360,   360,   270,   360,   270}
    ,{   340,   340,   310,   270,   310}
    ,{   220,   220,   170,   130,   170}
    ,{   340,   340,   310,   270,   310}
    }
   ,{{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   330,   310}
    }
   ,{{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   340,   340,   310,   270,   310}
    ,{   270,   210,   180,   270,   180}
    ,{   340,   340,   310,   270,   310}
    }
   ,{{   400,   400,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   400,   400,   310,   270,   310}
    }
   }
  ,{{{   400,   340,   400,   340,   400}
    ,{   400,   340,   400,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   270,   270,   270,   270,   270}
    ,{   310,   310,   310,   310,   310}
    ,{   230,   170,   230,   170,   230}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   400,   340,   400,   340,   400}
    ,{   400,   340,   400,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   180,   180,   180,   180,   180}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  ,{{{   340,   290,   340,   340,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   310,   310}
    ,{   310,   290,   310,   310,   310}
    }
   ,{{   310,   230,   310,   180,   310}
    ,{   270,   190,   270,   140,   270}
    ,{   310,   230,   310,   180,   310}
    ,{   170,    20,   170,   170,   170}
    ,{   310,   230,   310,   180,   310}
    }
   ,{{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    }
   ,{{   340,   260,   340,   340,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   230,   180,   310,   180}
    ,{   310,   230,   310,   180,   310}
    }
   ,{{   310,   290,   310,   310,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   230,   310,   310,   310}
    }
   }
  ,{{{   400,   340,   400,   340,   340}
    ,{   400,   340,   400,   340,   340}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   330,   310,   310,   310,   330}
    ,{   330,   270,   270,   270,   330}
    ,{   310,   310,   310,   310,   310}
    ,{   230,   170,   230,   170,   170}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   400,   340,   400,   340,   340}
    ,{   400,   340,   400,   340,   340}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   180,   180,   180,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  }
 ,{{{{   370,   340,   310,   310,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   320,   320,   290,   310,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   320,   320,   290,   310,   290}
    }
   ,{{   370,   340,   310,   310,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   320,   320,   280,   280,   280}
    ,{   240,   220,   240,   180,   240}
    ,{   320,   320,   280,   280,   280}
    }
   ,{{   330,   330,   290,   310,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   320,   320,   290,   310,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   320,   320,   290,   310,   290}
    }
   ,{{   320,   320,   310,   280,   310}
    ,{   310,   290,   310,   250,   310}
    ,{   320,   320,   280,   280,   280}
    ,{   260,   180,   130,   260,   260}
    ,{   320,   320,   280,   280,   280}
    }
   ,{{   330,   330,   290,   310,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   320,   320,   290,   310,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   290,   290,   200,   200,   200}
    }
   }
  ,{{{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   320,   320,   290,   310,   290}
    ,{   330,   330,   290,   250,   290}
    ,{   320,   320,   290,   310,   290}
    }
   ,{{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   320,   320,   280,   240,   280}
    ,{   220,   220,   180,   140,   180}
    ,{   320,   320,   280,   240,   280}
    }
   ,{{   330,   330,   290,   310,   290}
    ,{   330,   330,   290,   250,   290}
    ,{   320,   320,   290,   310,   290}
    ,{   330,   330,   290,   250,   290}
    ,{   320,   320,   290,   310,   290}
    }
   ,{{   320,   320,   280,   240,   280}
    ,{   290,   290,   250,   210,   250}
    ,{   320,   320,   280,   240,   280}
    ,{   220,   170,   130,   220,   130}
    ,{   320,   320,   280,   240,   280}
    }
   ,{{   330,   330,   290,   310,   290}
    ,{   330,   330,   290,   250,   290}
    ,{   320,   320,   290,   310,   290}
    ,{   330,   330,   290,   250,   290}
    ,{   290,   290,   200,   160,   200}
    }
   }
  ,{{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   240,   180,   240,   180,   240}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   310,   280,   310,   280,   310}
    ,{   310,   250,   310,   250,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   130,   130,   130,   130,   130}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   200,   200,   200,   200,   200}
    }
   }
  ,{{{   310,   270,   310,   260,   310}
    ,{   310,   230,   310,   250,   310}
    ,{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   260,   290}
    ,{   290,   270,   290,   200,   290}
    }
   ,{{   310,   230,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   280,   200,   280,   150,   280}
    ,{   180,   100,   180,   180,   180}
    ,{   280,   200,   280,   150,   280}
    }
   ,{{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   290,   270,   290,   160,   290}
    }
   ,{{   280,   200,   280,   260,   280}
    ,{   250,   170,   250,   250,   250}
    ,{   280,   200,   280,   150,   280}
    ,{   260,   180,   130,   260,   130}
    ,{   280,   200,   280,   150,   280}
    }
   ,{{   290,   270,   290,   200,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   200,   120,   200,   200,   200}
    }
   }
  ,{{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   280,   280,   280,   280,   280}
    ,{   240,   180,   240,   180,   180}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   310,   280,   310,   280,   280}
    ,{   310,   250,   310,   250,   250}
    ,{   280,   280,   280,   280,   280}
    ,{   260,   130,   130,   130,   260}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   200,   200,   200,   200,   200}
    }
   }
  }
 ,{{{{   370,   340,   310,   330,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   330,   310}
    }
   ,{{   370,   340,   310,   310,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   300,   300,   260,   260,   260}
    ,{   260,   240,   260,   200,   260}
    ,{   300,   300,   260,   260,   260}
    }
   ,{{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   330,   310}
    }
   ,{{   300,   300,   270,   280,   280}
    ,{   270,   250,   270,   210,   270}
    ,{   300,   300,   260,   260,   260}
    ,{   280,   200,   150,   280,   280}
    ,{   300,   300,   260,   260,   260}
    }
   ,{{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   310,   310,   280,   300,   280}
    ,{   340,   340,   310,   310,   310}
    ,{   320,   320,   220,   220,   220}
    }
   }
  ,{{{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   330,   310}
    }
   ,{{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   300,   300,   260,   220,   260}
    ,{   240,   240,   200,   160,   200}
    ,{   300,   300,   260,   220,   260}
    }
   ,{{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   330,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   340,   340,   310,   330,   310}
    }
   ,{{   300,   300,   260,   240,   260}
    ,{   250,   250,   210,   170,   210}
    ,{   300,   300,   260,   220,   260}
    ,{   240,   190,   150,   240,   150}
    ,{   300,   300,   260,   220,   260}
    }
   ,{{   340,   340,   310,   300,   310}
    ,{   340,   340,   310,   270,   310}
    ,{   310,   310,   280,   300,   280}
    ,{   340,   340,   310,   270,   310}
    ,{   320,   320,   220,   180,   220}
    }
   }
  ,{{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   200,   260,   200,   260}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   270,   260,   270,   260,   270}
    ,{   270,   210,   270,   210,   270}
    ,{   260,   260,   260,   260,   260}
    ,{   150,   150,   150,   150,   150}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   310,   310,   310,   310,   310}
    ,{   220,   220,   220,   220,   220}
    }
   }
  ,{{{   310,   290,   310,   280,   310}
    ,{   310,   230,   310,   210,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   280,   310}
    ,{   310,   290,   310,   220,   310}
    }
   ,{{   310,   230,   310,   200,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   260,   180,   260,   130,   260}
    ,{   200,   120,   200,   200,   200}
    ,{   260,   180,   260,   130,   260}
    }
   ,{{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    }
   ,{{   280,   200,   260,   280,   260}
    ,{   210,   130,   210,   210,   210}
    ,{   260,   180,   260,   130,   260}
    ,{   280,   200,   150,   280,   150}
    ,{   260,   180,   260,   130,   260}
    }
   ,{{   310,   260,   310,   220,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   280,   260,   280,   150,   280}
    ,{   310,   230,   310,   180,   310}
    ,{   220,   140,   220,   220,   220}
    }
   }
  ,{{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   200,   260,   200,   200}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   280,   260,   270,   260,   280}
    ,{   270,   210,   270,   210,   210}
    ,{   260,   260,   260,   260,   260}
    ,{   280,   150,   150,   150,   280}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   310,   310,   310,   310,   310}
    ,{   220,   220,   220,   220,   220}
    }
   }
  }
 ,{{{{   430,   430,   400,   370,   430}
    ,{   430,   410,   400,   370,   430}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   360,   340}
    }
   ,{{   430,   410,   370,   370,   430}
    ,{   430,   410,   370,   370,   430}
    ,{   370,   370,   340,   340,   340}
    ,{   320,   290,   320,   260,   320}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   360,   340}
    }
   ,{{   400,   370,   400,   340,   400}
    ,{   400,   370,   400,   340,   400}
    ,{   370,   370,   340,   340,   340}
    ,{   340,   260,   210,   340,   340}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   430,   430,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   340,   340}
    }
   }
  ,{{{   430,   430,   370,   360,   370}
    ,{   410,   410,   370,   360,   370}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   430,   430,   340,   360,   340}
    }
   ,{{   410,   410,   370,   360,   370}
    ,{   410,   410,   370,   360,   370}
    ,{   370,   370,   340,   300,   340}
    ,{   290,   290,   260,   220,   260}
    ,{   370,   370,   340,   300,   340}
    }
   ,{{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   360,   340}
    }
   ,{{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   300,   240,   210,   300,   210}
    ,{   370,   370,   340,   300,   340}
    }
   ,{{   430,   430,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   370,   340,   300,   340}
    ,{   430,   430,   340,   300,   340}
    }
   }
  ,{{{   400,   370,   400,   370,   400}
    ,{   400,   370,   400,   370,   400}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   370,   370,   370,   370,   370}
    ,{   370,   370,   370,   370,   370}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   320}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   400,   340,   400,   340,   400}
    ,{   400,   340,   400,   340,   400}
    ,{   340,   340,   340,   340,   340}
    ,{   210,   210,   210,   210,   210}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  ,{{{   370,   320,   370,   340,   370}
    ,{   370,   290,   370,   340,   370}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   340,   320,   340,   340,   340}
    }
   ,{{   370,   290,   370,   260,   370}
    ,{   370,   290,   370,   240,   370}
    ,{   340,   260,   340,   210,   340}
    ,{   260,   180,   260,   260,   260}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    }
   ,{{   340,   260,   340,   340,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   210,   340,   210}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   340,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    }
   }
  ,{{{   430,   370,   400,   370,   430}
    ,{   430,   370,   400,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   430,   370,   370,   370,   430}
    ,{   430,   370,   370,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   260}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   400,   340,   400,   340,   340}
    ,{   400,   340,   400,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   210,   210,   210,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   310,   240,   240,   310,   260}
    ,{   270,   240,   240,   270,   260}
    ,{   310,   220,   220,   310,   220}
    ,{   270,   240,   240,   270,   240}
    ,{   300,   210,   210,   300,   210}
    }
   ,{{   260,   200,   200,   230,   260}
    ,{   260,   200,   200,   230,   260}
    ,{   220,   190,   190,   220,   190}
    ,{   160,   100,   160,   130,   160}
    ,{   220,   190,   190,   220,   190}
    }
   ,{{   310,   240,   240,   310,   240}
    ,{   270,   240,   240,   270,   240}
    ,{   310,   220,   220,   310,   220}
    ,{   270,   240,   240,   270,   240}
    ,{   300,   210,   210,   300,   210}
    }
   ,{{   220,   190,   190,   220,   190}
    ,{   160,   100,   160,   130,   160}
    ,{   220,   190,   190,   220,   190}
    ,{   210,    50,    50,   210,   180}
    ,{   220,   190,   190,   220,   190}
    }
   ,{{   300,   240,   240,   300,   240}
    ,{   270,   240,   240,   270,   240}
    ,{   300,   210,   210,   300,   210}
    ,{   270,   240,   240,   270,   240}
    ,{   150,   140,   120,   150,   120}
    }
   }
  ,{{{   310,   200,   240,   310,   240}
    ,{   270,   200,   240,   270,   240}
    ,{   310,   190,   220,   310,   220}
    ,{   270,   200,   240,   270,   240}
    ,{   300,   170,   210,   300,   210}
    }
   ,{{   230,   160,   200,   230,   200}
    ,{   230,   160,   200,   230,   200}
    ,{   220,   160,   190,   220,   190}
    ,{   130,    70,   100,   130,   100}
    ,{   220,   160,   190,   220,   190}
    }
   ,{{   310,   200,   240,   310,   240}
    ,{   270,   200,   240,   270,   240}
    ,{   310,   190,   220,   310,   220}
    ,{   270,   200,   240,   270,   240}
    ,{   300,   170,   210,   300,   210}
    }
   ,{{   220,   160,   190,   220,   190}
    ,{   130,    70,   100,   130,   100}
    ,{   220,   160,   190,   220,   190}
    ,{   210,    10,    50,   210,    50}
    ,{   220,   160,   190,   220,   190}
    }
   ,{{   300,   200,   240,   300,   240}
    ,{   270,   200,   240,   270,   240}
    ,{   300,   170,   210,   300,   210}
    ,{   270,   200,   240,   270,   240}
    ,{   150,   140,   120,   150,   120}
    }
   }
  ,{{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   220,   220,   220,   220,   220}
    ,{   240,   240,   240,   240,   240}
    ,{   210,   210,   210,   210,   210}
    }
   ,{{   200,   200,   200,   200,   200}
    ,{   200,   200,   200,   200,   200}
    ,{   190,   190,   190,   190,   190}
    ,{   160,   100,   160,   100,   160}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   220,   220,   220,   220,   220}
    ,{   240,   240,   240,   240,   240}
    ,{   210,   210,   210,   210,   210}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   160,   100,   160,   100,   160}
    ,{   190,   190,   190,   190,   190}
    ,{    50,    50,    50,    50,    50}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   210,   210,   210,   210,   210}
    ,{   240,   240,   240,   240,   240}
    ,{   120,   120,   120,   120,   120}
    }
   }
  ,{{{   240,   150,   240,   180,   240}
    ,{   240,   100,   240,   110,   240}
    ,{   220,   150,   220,    90,   220}
    ,{   240,   100,   240,   180,   240}
    ,{   210,   130,   210,   120,   210}
    }
   ,{{   200,    60,   200,   100,   200}
    ,{   200,    60,   200,    70,   200}
    ,{   190,    60,   190,    60,   190}
    ,{   100,   -30,   100,   100,   100}
    ,{   190,    60,   190,    60,   190}
    }
   ,{{   240,   150,   240,   110,   240}
    ,{   240,   100,   240,   110,   240}
    ,{   220,   150,   220,    90,   220}
    ,{   240,   100,   240,   110,   240}
    ,{   210,   130,   210,    80,   210}
    }
   ,{{   190,    60,   190,   180,   190}
    ,{   100,   -30,   100,   100,   100}
    ,{   190,    60,   190,    60,   190}
    ,{   180,    40,    50,   180,    50}
    ,{   190,    60,   190,    60,   190}
    }
   ,{{   240,   130,   240,   120,   240}
    ,{   240,   100,   240,   110,   240}
    ,{   210,   130,   210,    80,   210}
    ,{   240,   100,   240,   110,   240}
    ,{   120,   -10,   120,   120,   120}
    }
   }
  ,{{{   260,   240,   240,   240,   260}
    ,{   260,   240,   240,   240,   260}
    ,{   220,   220,   220,   220,   220}
    ,{   240,   240,   240,   240,   240}
    ,{   210,   210,   210,   210,   210}
    }
   ,{{   260,   200,   200,   200,   260}
    ,{   260,   200,   200,   200,   260}
    ,{   190,   190,   190,   190,   190}
    ,{   160,   100,   160,   100,   100}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   220,   220,   220,   220,   220}
    ,{   240,   240,   240,   240,   240}
    ,{   210,   210,   210,   210,   210}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   160,   100,   160,   100,   100}
    ,{   190,   190,   190,   190,   190}
    ,{   180,    50,    50,    50,   180}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   210,   210,   210,   210,   210}
    ,{   240,   240,   240,   240,   240}
    ,{   120,   120,   120,   120,   120}
    }
   }
  }
 ,{{{{   280,   210,   210,   280,   270}
    ,{   270,   210,   210,   240,   270}
    ,{   280,   190,   190,   280,   190}
    ,{   210,   180,   180,   210,   180}
    ,{   280,   190,   190,   280,   190}
    }
   ,{{   270,   210,   210,   240,   270}
    ,{   270,   210,   210,   240,   270}
    ,{   220,   190,   190,   220,   190}
    ,{    70,    10,    70,    40,    70}
    ,{   220,   190,   190,   220,   190}
    }
   ,{{   280,   190,   190,   280,   190}
    ,{   210,   180,   180,   210,   180}
    ,{   280,   190,   190,   280,   190}
    ,{   210,   180,   180,   210,   180}
    ,{   280,   190,   190,   280,   190}
    }
   ,{{   220,   190,   190,   220,   190}
    ,{   130,    70,   130,   100,   130}
    ,{   220,   190,   190,   220,   190}
    ,{   210,    50,    50,   210,   180}
    ,{   220,   190,   190,   220,   190}
    }
   ,{{   280,   190,   190,   280,   190}
    ,{   210,   180,   180,   210,   180}
    ,{   280,   190,   190,   280,   190}
    ,{   210,   180,   180,   210,   180}
    ,{   140,   140,   110,   140,   110}
    }
   }
  ,{{{   280,   190,   210,   280,   210}
    ,{   240,   190,   210,   240,   210}
    ,{   280,   160,   190,   280,   190}
    ,{   210,   150,   180,   210,   180}
    ,{   280,   150,   190,   280,   190}
    }
   ,{{   240,   190,   210,   240,   210}
    ,{   240,   190,   210,   240,   210}
    ,{   220,   150,   190,   220,   190}
    ,{    40,   -20,    10,    40,    10}
    ,{   220,   150,   190,   220,   190}
    }
   ,{{   280,   150,   190,   280,   190}
    ,{   210,   150,   180,   210,   180}
    ,{   280,   150,   190,   280,   190}
    ,{   210,   150,   180,   210,   180}
    ,{   280,   150,   190,   280,   190}
    }
   ,{{   220,   150,   190,   220,   190}
    ,{   100,    40,    70,   100,    70}
    ,{   220,   150,   190,   220,   190}
    ,{   210,    10,    50,   210,    50}
    ,{   220,   150,   190,   220,   190}
    }
   ,{{   280,   160,   190,   280,   190}
    ,{   210,   150,   180,   210,   180}
    ,{   280,   160,   190,   280,   190}
    ,{   210,   150,   180,   210,   180}
    ,{   140,   140,   110,   140,   110}
    }
   }
  ,{{{   210,   210,   210,   210,   210}
    ,{   210,   210,   210,   210,   210}
    ,{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   210,   210,   210,   210,   210}
    ,{   210,   210,   210,   210,   210}
    ,{   190,   190,   190,   190,   190}
    ,{    70,    10,    70,    10,    70}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   130,    70,   130,    70,   130}
    ,{   190,   190,   190,   190,   190}
    ,{    50,    50,    50,    50,    50}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   110,   110,   110,   110,   110}
    }
   }
  ,{{{   210,   120,   210,   180,   210}
    ,{   210,    80,   210,    80,   210}
    ,{   190,   120,   190,    60,   190}
    ,{   180,    50,   180,   180,   180}
    ,{   190,   110,   190,   110,   190}
    }
   ,{{   210,    80,   210,    80,   210}
    ,{   210,    80,   210,    80,   210}
    ,{   190,    50,   190,    60,   190}
    ,{    10,  -120,    10,    10,    10}
    ,{   190,    50,   190,    60,   190}
    }
   ,{{   190,   110,   190,    60,   190}
    ,{   180,    50,   180,    50,   180}
    ,{   190,   110,   190,    60,   190}
    ,{   180,    50,   180,    50,   180}
    ,{   190,   110,   190,    60,   190}
    }
   ,{{   190,    50,   190,   180,   190}
    ,{    70,   -60,    70,    70,    70}
    ,{   190,    50,   190,    60,   190}
    ,{   180,    40,    50,   180,    50}
    ,{   190,    50,   190,    60,   190}
    }
   ,{{   190,   120,   190,   110,   190}
    ,{   180,    50,   180,    50,   180}
    ,{   190,   120,   190,    60,   190}
    ,{   180,    50,   180,    50,   180}
    ,{   110,   -20,   110,   110,   110}
    }
   }
  ,{{{   270,   210,   210,   210,   270}
    ,{   270,   210,   210,   210,   270}
    ,{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   270,   210,   210,   210,   270}
    ,{   270,   210,   210,   210,   270}
    ,{   190,   190,   190,   190,   190}
    ,{    70,    10,    70,    10,    10}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   130,    70,   130,    70,    70}
    ,{   190,   190,   190,   190,   190}
    ,{   180,    50,    50,    50,   180}
    ,{   190,   190,   190,   190,   190}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   190,   190,   190,   190,   190}
    ,{   180,   180,   180,   180,   180}
    ,{   110,   110,   110,   110,   110}
    }
   }
  }
 ,{{{{   400,   360,   340,   400,   400}
    ,{   400,   360,   340,   370,   400}
    ,{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   330,   310,   400,   310}
    }
   ,{{   400,   360,   340,   370,   400}
    ,{   400,   360,   340,   370,   400}
    ,{   340,   310,   310,   340,   310}
    ,{   290,   230,   290,   260,   290}
    ,{   340,   310,   310,   340,   310}
    }
   ,{{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   310,   310,   400,   310}
    }
   ,{{   360,   360,   330,   340,   330}
    ,{   360,   360,   330,   300,   330}
    ,{   340,   310,   310,   340,   310}
    ,{   340,   180,   180,   340,   310}
    ,{   340,   310,   310,   340,   310}
    }
   ,{{   400,   330,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   340,   330,   310,   340,   310}
    }
   }
  ,{{{   400,   360,   340,   400,   340}
    ,{   370,   360,   340,   370,   340}
    ,{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   330,   310,   400,   310}
    }
   ,{{   370,   360,   340,   370,   340}
    ,{   370,   360,   340,   370,   340}
    ,{   340,   270,   310,   340,   310}
    ,{   260,   190,   230,   260,   230}
    ,{   340,   270,   310,   340,   310}
    }
   ,{{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   270,   310,   400,   310}
    }
   ,{{   360,   360,   310,   340,   310}
    ,{   360,   360,   270,   300,   270}
    ,{   340,   270,   310,   340,   310}
    ,{   340,   140,   180,   340,   180}
    ,{   340,   270,   310,   340,   310}
    }
   ,{{   400,   330,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   340,   330,   310,   340,   310}
    }
   }
  ,{{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   310,   310,   310,   310,   310}
    ,{   290,   230,   290,   230,   290}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   330,   310,   330,   310,   330}
    ,{   330,   270,   330,   270,   330}
    ,{   310,   310,   310,   310,   310}
    ,{   180,   180,   180,   180,   180}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  ,{{{   340,   230,   340,   310,   340}
    ,{   340,   220,   340,   270,   340}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   310,   310}
    ,{   310,   230,   310,   310,   310}
    }
   ,{{   340,   220,   340,   230,   340}
    ,{   340,   220,   340,   210,   340}
    ,{   310,   170,   310,   180,   310}
    ,{   230,    20,   230,   230,   230}
    ,{   310,   170,   310,   180,   310}
    }
   ,{{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    }
   ,{{   310,   170,   310,   310,   310}
    ,{   270,   130,   270,   270,   270}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   170,   180,   310,   180}
    ,{   310,   170,   310,   180,   310}
    }
   ,{{   310,   230,   310,   310,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   170,   310,   310,   310}
    }
   }
  ,{{{   400,   340,   340,   340,   400}
    ,{   400,   340,   340,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   400,   340,   340,   340,   400}
    ,{   400,   340,   340,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   290,   230,   290,   230,   230}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   330,   310,   330,   310,   310}
    ,{   330,   270,   330,   270,   270}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   180,   180,   180,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  }
 ,{{{{   370,   310,   370,   370,   370}
    ,{   370,   310,   370,   340,   370}
    ,{   370,   280,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   370,   300,   280,   370,   280}
    }
   ,{{   310,   280,   280,   310,   300}
    ,{   300,   240,   240,   270,   300}
    ,{   310,   280,   280,   310,   280}
    ,{   200,   140,   200,   170,   200}
    ,{   310,   280,   280,   310,   280}
    }
   ,{{   370,   280,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   370,   280,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   370,   280,   280,   370,   280}
    }
   ,{{   370,   310,   370,   340,   370}
    ,{   370,   310,   370,   340,   370}
    ,{   310,   280,   280,   310,   280}
    ,{   310,   150,   150,   310,   280}
    ,{   310,   280,   280,   310,   280}
    }
   ,{{   370,   300,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   370,   280,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   310,   300,   280,   310,   280}
    }
   }
  ,{{{   370,   300,   310,   370,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   300,   280,   370,   280}
    }
   ,{{   310,   240,   280,   310,   280}
    ,{   270,   210,   240,   270,   240}
    ,{   310,   240,   280,   310,   280}
    ,{   170,   110,   140,   170,   140}
    ,{   310,   240,   280,   310,   280}
    }
   ,{{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   240,   280,   370,   280}
    }
   ,{{   340,   270,   310,   340,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   310,   240,   280,   310,   280}
    ,{   310,   110,   150,   310,   150}
    ,{   310,   240,   280,   310,   280}
    }
   ,{{   370,   300,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   310,   300,   280,   310,   280}
    }
   }
  ,{{{   370,   310,   370,   310,   370}
    ,{   370,   310,   370,   310,   370}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   240,   240,   240,   240,   240}
    ,{   280,   280,   280,   280,   280}
    ,{   200,   140,   200,   140,   200}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   370,   310,   370,   310,   370}
    ,{   370,   310,   370,   310,   370}
    ,{   280,   280,   280,   280,   280}
    ,{   150,   150,   150,   150,   150}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   }
  ,{{{   310,   200,   310,   310,   310}
    ,{   310,   170,   310,   310,   310}
    ,{   280,   200,   280,   150,   280}
    ,{   280,   140,   280,   280,   280}
    ,{   280,   200,   280,   280,   280}
    }
   ,{{   280,   140,   280,   150,   280}
    ,{   240,   110,   240,   110,   240}
    ,{   280,   140,   280,   150,   280}
    ,{   140,    10,   140,   140,   140}
    ,{   280,   140,   280,   150,   280}
    }
   ,{{   280,   200,   280,   150,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   280,   200,   280,   150,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   280,   200,   280,   150,   280}
    }
   ,{{   310,   170,   310,   310,   310}
    ,{   310,   170,   310,   310,   310}
    ,{   280,   140,   280,   150,   280}
    ,{   280,   140,   150,   280,   150}
    ,{   280,   140,   280,   150,   280}
    }
   ,{{   280,   200,   280,   280,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   280,   200,   280,   150,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   280,   140,   280,   280,   280}
    }
   }
  ,{{{   370,   310,   370,   310,   310}
    ,{   370,   310,   370,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   300,   280,   280,   280,   300}
    ,{   300,   240,   240,   240,   300}
    ,{   280,   280,   280,   280,   280}
    ,{   200,   140,   200,   140,   140}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   370,   310,   370,   310,   310}
    ,{   370,   310,   370,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   150,   150,   150,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   }
  }
 ,{{{{   350,   280,   280,   350,   340}
    ,{   340,   280,   280,   310,   340}
    ,{   350,   260,   260,   350,   260}
    ,{   290,   260,   260,   290,   260}
    ,{   350,   260,   260,   350,   260}
    }
   ,{{   340,   280,   280,   310,   340}
    ,{   340,   280,   280,   310,   340}
    ,{   280,   250,   250,   280,   250}
    ,{   210,   150,   210,   180,   210}
    ,{   280,   250,   250,   280,   250}
    }
   ,{{   350,   260,   260,   350,   260}
    ,{   290,   260,   260,   290,   260}
    ,{   350,   260,   260,   350,   260}
    ,{   290,   260,   260,   290,   260}
    ,{   350,   260,   260,   350,   260}
    }
   ,{{   280,   250,   280,   280,   280}
    ,{   280,   220,   280,   250,   280}
    ,{   280,   250,   250,   280,   250}
    ,{   260,   100,   100,   260,   230}
    ,{   280,   250,   250,   280,   250}
    }
   ,{{   350,   260,   260,   350,   260}
    ,{   290,   260,   260,   290,   260}
    ,{   350,   260,   260,   350,   260}
    ,{   290,   260,   260,   290,   260}
    ,{   200,   190,   170,   200,   170}
    }
   }
  ,{{{   350,   240,   280,   350,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   350,   220,   260,   350,   260}
    ,{   290,   230,   260,   290,   260}
    ,{   350,   220,   260,   350,   260}
    }
   ,{{   310,   240,   280,   310,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   280,   220,   250,   280,   250}
    ,{   180,   120,   150,   180,   150}
    ,{   280,   220,   250,   280,   250}
    }
   ,{{   350,   230,   260,   350,   260}
    ,{   290,   230,   260,   290,   260}
    ,{   350,   220,   260,   350,   260}
    ,{   290,   230,   260,   290,   260}
    ,{   350,   220,   260,   350,   260}
    }
   ,{{   280,   220,   250,   280,   250}
    ,{   250,   190,   220,   250,   220}
    ,{   280,   220,   250,   280,   250}
    ,{   260,    70,   100,   260,   100}
    ,{   280,   220,   250,   280,   250}
    }
   ,{{   350,   230,   260,   350,   260}
    ,{   290,   230,   260,   290,   260}
    ,{   350,   220,   260,   350,   260}
    ,{   290,   230,   260,   290,   260}
    ,{   200,   190,   170,   200,   170}
    }
   }
  ,{{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   250,   250,   250,   250,   250}
    ,{   210,   150,   210,   150,   210}
    ,{   250,   250,   250,   250,   250}
    }
   ,{{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   280,   250,   280,   250,   280}
    ,{   280,   220,   280,   220,   280}
    ,{   250,   250,   250,   250,   250}
    ,{   100,   100,   100,   100,   100}
    ,{   250,   250,   250,   250,   250}
    }
   ,{{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   170,   170,   170,   170,   170}
    }
   }
  ,{{{   280,   180,   280,   230,   280}
    ,{   280,   140,   280,   220,   280}
    ,{   260,   180,   260,   130,   260}
    ,{   260,   130,   260,   230,   260}
    ,{   260,   180,   260,   170,   260}
    }
   ,{{   280,   140,   280,   150,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   250,   120,   250,   120,   250}
    ,{   150,    20,   150,   150,   150}
    ,{   250,   120,   250,   120,   250}
    }
   ,{{   260,   180,   260,   130,   260}
    ,{   260,   130,   260,   130,   260}
    ,{   260,   180,   260,   130,   260}
    ,{   260,   130,   260,   130,   260}
    ,{   260,   180,   260,   130,   260}
    }
   ,{{   250,   120,   250,   230,   250}
    ,{   220,    90,   220,   220,   220}
    ,{   250,   120,   250,   120,   250}
    ,{   230,   100,   100,   230,   100}
    ,{   250,   120,   250,   120,   250}
    }
   ,{{   260,   180,   260,   170,   260}
    ,{   260,   130,   260,   130,   260}
    ,{   260,   180,   260,   130,   260}
    ,{   260,   130,   260,   130,   260}
    ,{   170,    30,   170,   170,   170}
    }
   }
  ,{{{   340,   280,   280,   280,   340}
    ,{   340,   280,   280,   280,   340}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   340,   280,   280,   280,   340}
    ,{   340,   280,   280,   280,   340}
    ,{   250,   250,   250,   250,   250}
    ,{   210,   150,   210,   150,   150}
    ,{   250,   250,   250,   250,   250}
    }
   ,{{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   280,   250,   280,   250,   250}
    ,{   280,   220,   280,   220,   220}
    ,{   250,   250,   250,   250,   250}
    ,{   230,   100,   100,   100,   230}
    ,{   250,   250,   250,   250,   250}
    }
   ,{{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   260,   260,   260,   260}
    ,{   170,   170,   170,   170,   170}
    }
   }
  }
 ,{{{{   370,   280,   280,   370,   340}
    ,{   340,   280,   280,   310,   340}
    ,{   370,   280,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   370,   280,   280,   370,   280}
    }
   ,{{   340,   280,   280,   310,   340}
    ,{   340,   280,   280,   310,   340}
    ,{   260,   230,   230,   260,   230}
    ,{   230,   170,   230,   200,   230}
    ,{   260,   230,   230,   260,   230}
    }
   ,{{   370,   280,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   370,   280,   280,   370,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   370,   280,   280,   370,   280}
    }
   ,{{   280,   230,   240,   280,   250}
    ,{   240,   180,   240,   210,   240}
    ,{   260,   230,   230,   260,   230}
    ,{   280,   120,   120,   280,   250}
    ,{   260,   230,   230,   260,   230}
    }
   ,{{   340,   280,   280,   340,   280}
    ,{   310,   280,   280,   310,   280}
    ,{   340,   250,   250,   340,   250}
    ,{   310,   280,   280,   310,   280}
    ,{   220,   220,   190,   220,   190}
    }
   }
  ,{{{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   240,   280,   370,   280}
    }
   ,{{   310,   240,   280,   310,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   260,   200,   230,   260,   230}
    ,{   200,   140,   170,   200,   170}
    ,{   260,   200,   230,   260,   230}
    }
   ,{{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   240,   280,   370,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   370,   240,   280,   370,   280}
    }
   ,{{   280,   200,   230,   280,   230}
    ,{   210,   150,   180,   210,   180}
    ,{   260,   200,   230,   260,   230}
    ,{   280,    90,   120,   280,   120}
    ,{   260,   200,   230,   260,   230}
    }
   ,{{   340,   240,   280,   340,   280}
    ,{   310,   240,   280,   310,   280}
    ,{   340,   210,   250,   340,   250}
    ,{   310,   240,   280,   310,   280}
    ,{   220,   220,   190,   220,   190}
    }
   }
  ,{{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   230,   230,   230,   230,   230}
    ,{   230,   170,   230,   170,   230}
    ,{   230,   230,   230,   230,   230}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   240,   230,   240,   230,   240}
    ,{   240,   180,   240,   180,   240}
    ,{   230,   230,   230,   230,   230}
    ,{   120,   120,   120,   120,   120}
    ,{   230,   230,   230,   230,   230}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   250,   250,   250,   250,   250}
    ,{   280,   280,   280,   280,   280}
    ,{   190,   190,   190,   190,   190}
    }
   }
  ,{{{   280,   200,   280,   250,   280}
    ,{   280,   140,   280,   180,   280}
    ,{   280,   200,   280,   150,   280}
    ,{   280,   140,   280,   250,   280}
    ,{   280,   200,   280,   190,   280}
    }
   ,{{   280,   140,   280,   170,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   230,   100,   230,   100,   230}
    ,{   170,    40,   170,   170,   170}
    ,{   230,   100,   230,   100,   230}
    }
   ,{{   280,   200,   280,   150,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   280,   200,   280,   150,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   280,   200,   280,   150,   280}
    }
   ,{{   250,   120,   230,   250,   230}
    ,{   180,    50,   180,   180,   180}
    ,{   230,   100,   230,   100,   230}
    ,{   250,   120,   120,   250,   120}
    ,{   230,   100,   230,   100,   230}
    }
   ,{{   280,   170,   280,   190,   280}
    ,{   280,   140,   280,   150,   280}
    ,{   250,   170,   250,   120,   250}
    ,{   280,   140,   280,   150,   280}
    ,{   190,    60,   190,   190,   190}
    }
   }
  ,{{{   340,   280,   280,   280,   340}
    ,{   340,   280,   280,   280,   340}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   340,   280,   280,   280,   340}
    ,{   340,   280,   280,   280,   340}
    ,{   230,   230,   230,   230,   230}
    ,{   230,   170,   230,   170,   170}
    ,{   230,   230,   230,   230,   230}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   250,   230,   240,   230,   250}
    ,{   240,   180,   240,   180,   180}
    ,{   230,   230,   230,   230,   230}
    ,{   250,   120,   120,   120,   250}
    ,{   230,   230,   230,   230,   230}
    }
   ,{{   280,   280,   280,   280,   280}
    ,{   280,   280,   280,   280,   280}
    ,{   250,   250,   250,   250,   250}
    ,{   280,   280,   280,   280,   280}
    ,{   190,   190,   190,   190,   190}
    }
   }
  }
 ,{{{{   400,   360,   370,   400,   400}
    ,{   400,   360,   370,   370,   400}
    ,{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   330,   310,   400,   310}
    }
   ,{{   400,   360,   340,   370,   400}
    ,{   400,   360,   340,   370,   400}
    ,{   340,   310,   310,   340,   310}
    ,{   290,   230,   290,   260,   290}
    ,{   340,   310,   310,   340,   310}
    }
   ,{{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   310,   310,   400,   310}
    }
   ,{{   370,   360,   370,   340,   370}
    ,{   370,   360,   370,   340,   370}
    ,{   340,   310,   310,   340,   310}
    ,{   340,   180,   180,   340,   310}
    ,{   340,   310,   310,   340,   310}
    }
   ,{{   400,   330,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   400,   310,   310,   400,   310}
    ,{   340,   310,   310,   340,   310}
    ,{   340,   330,   310,   340,   310}
    }
   }
  ,{{{   400,   360,   340,   400,   340}
    ,{   370,   360,   340,   370,   340}
    ,{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   330,   310,   400,   310}
    }
   ,{{   370,   360,   340,   370,   340}
    ,{   370,   360,   340,   370,   340}
    ,{   340,   270,   310,   340,   310}
    ,{   260,   190,   230,   260,   230}
    ,{   340,   270,   310,   340,   310}
    }
   ,{{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   270,   310,   400,   310}
    }
   ,{{   360,   360,   310,   340,   310}
    ,{   360,   360,   310,   340,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   340,   140,   180,   340,   180}
    ,{   340,   270,   310,   340,   310}
    }
   ,{{   400,   330,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   400,   270,   310,   400,   310}
    ,{   340,   270,   310,   340,   310}
    ,{   340,   330,   310,   340,   310}
    }
   }
  ,{{{   370,   340,   370,   340,   370}
    ,{   370,   340,   370,   340,   370}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   310,   310,   310,   310,   310}
    ,{   290,   230,   290,   230,   290}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   370,   310,   370,   310,   370}
    ,{   370,   310,   370,   310,   370}
    ,{   310,   310,   310,   310,   310}
    ,{   180,   180,   180,   180,   180}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  ,{{{   340,   230,   340,   310,   340}
    ,{   340,   220,   340,   310,   340}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   310,   310}
    ,{   310,   230,   310,   310,   310}
    }
   ,{{   340,   220,   340,   230,   340}
    ,{   340,   220,   340,   210,   340}
    ,{   310,   170,   310,   180,   310}
    ,{   230,    40,   230,   230,   230}
    ,{   310,   170,   310,   180,   310}
    }
   ,{{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    }
   ,{{   310,   170,   310,   310,   310}
    ,{   310,   170,   310,   310,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   170,   180,   310,   180}
    ,{   310,   170,   310,   180,   310}
    }
   ,{{   310,   230,   310,   310,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   170,   310,   180,   310}
    ,{   310,   170,   310,   310,   310}
    }
   }
  ,{{{   400,   340,   370,   340,   400}
    ,{   400,   340,   370,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   400,   340,   340,   340,   400}
    ,{   400,   340,   340,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   290,   230,   290,   230,   230}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   370,   310,   370,   310,   310}
    ,{   370,   310,   370,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   180,   180,   180,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   240,   240,   220,   230,   220}
    ,{   240,   240,   220,   210,   220}
    ,{   230,   220,   210,   230,   210}
    ,{   240,   240,   220,   210,   220}
    ,{   210,   210,   190,   210,   190}
    }
   ,{{   200,   200,   180,   170,   180}
    ,{   200,   200,   180,   170,   180}
    ,{   190,   190,   180,   170,   180}
    ,{   140,   100,   140,    80,   140}
    ,{   190,   190,   180,   170,   180}
    }
   ,{{   240,   240,   220,   230,   220}
    ,{   240,   240,   220,   210,   220}
    ,{   230,   220,   210,   230,   210}
    ,{   240,   240,   220,   210,   220}
    ,{   210,   210,   190,   210,   190}
    }
   ,{{   190,   190,   180,   170,   180}
    ,{   140,   100,   140,    80,   140}
    ,{   190,   190,   180,   170,   180}
    ,{   130,    50,    30,   130,    70}
    ,{   190,   190,   180,   170,   180}
    }
   ,{{   240,   240,   220,   210,   220}
    ,{   240,   240,   220,   210,   220}
    ,{   210,   210,   190,   210,   190}
    ,{   240,   240,   220,   210,   220}
    ,{   180,   180,   100,    90,   100}
    }
   }
  ,{{{   240,   240,   220,   230,   220}
    ,{   240,   240,   220,   180,   220}
    ,{   230,   220,   210,   230,   210}
    ,{   240,   240,   220,   180,   220}
    ,{   210,   210,   190,   210,   190}
    }
   ,{{   200,   200,   180,   140,   180}
    ,{   200,   200,   180,   140,   180}
    ,{   190,   190,   180,   140,   180}
    ,{   100,   100,    90,    50,    90}
    ,{   190,   190,   180,   140,   180}
    }
   ,{{   240,   240,   220,   230,   220}
    ,{   240,   240,   220,   180,   220}
    ,{   230,   220,   210,   230,   210}
    ,{   240,   240,   220,   180,   220}
    ,{   210,   210,   190,   210,   190}
    }
   ,{{   190,   190,   180,   140,   180}
    ,{   100,   100,    90,    50,    90}
    ,{   190,   190,   180,   140,   180}
    ,{   120,    50,    30,   120,    30}
    ,{   190,   190,   180,   140,   180}
    }
   ,{{   240,   240,   220,   210,   220}
    ,{   240,   240,   220,   180,   220}
    ,{   210,   210,   190,   210,   190}
    ,{   240,   240,   220,   180,   220}
    ,{   180,   180,   100,    60,   100}
    }
   }
  ,{{{   220,   210,   220,   210,   220}
    ,{   220,   210,   220,   210,   220}
    ,{   200,   200,   200,   200,   200}
    ,{   220,   210,   220,   210,   220}
    ,{   190,   180,   190,   180,   190}
    }
   ,{{   180,   170,   180,   170,   180}
    ,{   180,   170,   180,   170,   180}
    ,{   170,   170,   170,   170,   170}
    ,{   140,    80,   140,    80,   140}
    ,{   170,   170,   170,   170,   170}
    }
   ,{{   220,   210,   220,   210,   220}
    ,{   220,   210,   220,   210,   220}
    ,{   200,   200,   200,   200,   200}
    ,{   220,   210,   220,   210,   220}
    ,{   190,   180,   190,   180,   190}
    }
   ,{{   170,   170,   170,   170,   170}
    ,{   140,    80,   140,    80,   140}
    ,{   170,   170,   170,   170,   170}
    ,{    30,    20,    30,    20,    30}
    ,{   170,   170,   170,   170,   170}
    }
   ,{{   220,   210,   220,   210,   220}
    ,{   220,   210,   220,   210,   220}
    ,{   190,   180,   190,   180,   190}
    ,{   220,   210,   220,   210,   220}
    ,{   100,    90,   100,    90,   100}
    }
   }
  ,{{{   220,   160,   220,   130,   220}
    ,{   220,   110,   220,    60,   220}
    ,{   210,   160,   210,    50,   210}
    ,{   220,   110,   220,   130,   220}
    ,{   190,   140,   190,    70,   190}
    }
   ,{{   180,    70,   180,    60,   180}
    ,{   180,    70,   180,    20,   180}
    ,{   180,    70,   180,    20,   180}
    ,{    90,   -20,    90,    60,    90}
    ,{   180,    70,   180,    20,   180}
    }
   ,{{   220,   160,   220,    60,   220}
    ,{   220,   110,   220,    60,   220}
    ,{   210,   160,   210,    50,   210}
    ,{   220,   110,   220,    60,   220}
    ,{   190,   140,   190,    30,   190}
    }
   ,{{   180,    70,   180,   130,   180}
    ,{    90,   -20,    90,    60,    90}
    ,{   180,    70,   180,    20,   180}
    ,{   130,    50,    30,   130,    30}
    ,{   180,    70,   180,    20,   180}
    }
   ,{{   220,   140,   220,    70,   220}
    ,{   220,   110,   220,    60,   220}
    ,{   190,   140,   190,    30,   190}
    ,{   220,   110,   220,    60,   220}
    ,{   100,     0,   100,    70,   100}
    }
   }
  ,{{{   220,   210,   220,   210,   150}
    ,{   220,   210,   220,   210,   150}
    ,{   200,   200,   200,   200,   110}
    ,{   220,   210,   220,   210,   130}
    ,{   190,   180,   190,   180,   100}
    }
   ,{{   180,   170,   180,   170,   150}
    ,{   180,   170,   180,   170,   150}
    ,{   170,   170,   170,   170,    80}
    ,{   140,    80,   140,    80,     0}
    ,{   170,   170,   170,   170,    80}
    }
   ,{{   220,   210,   220,   210,   130}
    ,{   220,   210,   220,   210,   130}
    ,{   200,   200,   200,   200,   110}
    ,{   220,   210,   220,   210,   130}
    ,{   190,   180,   190,   180,   100}
    }
   ,{{   170,   170,   170,   170,    80}
    ,{   140,    80,   140,    80,     0}
    ,{   170,   170,   170,   170,    80}
    ,{    70,    20,    30,    20,    70}
    ,{   170,   170,   170,   170,    80}
    }
   ,{{   220,   210,   220,   210,   130}
    ,{   220,   210,   220,   210,   130}
    ,{   190,   180,   190,   180,   100}
    ,{   220,   210,   220,   210,   130}
    ,{   100,    90,   100,    90,    10}
    }
   }
  }
 ,{{{{   210,   210,   200,   200,   200}
    ,{   210,   210,   200,   190,   200}
    ,{   200,   190,   180,   200,   180}
    ,{   180,   180,   170,   160,   170}
    ,{   190,   190,   170,   190,   170}
    }
   ,{{   210,   210,   200,   190,   200}
    ,{   210,   210,   200,   190,   200}
    ,{   190,   190,   170,   160,   170}
    ,{    50,    10,    50,   -10,    50}
    ,{   190,   190,   170,   160,   170}
    }
   ,{{   190,   190,   170,   190,   170}
    ,{   180,   180,   170,   160,   170}
    ,{   190,   190,   170,   190,   170}
    ,{   180,   180,   170,   160,   170}
    ,{   190,   190,   170,   190,   170}
    }
   ,{{   190,   190,   170,   160,   170}
    ,{   110,    70,   110,    50,   110}
    ,{   190,   190,   170,   160,   170}
    ,{   130,    50,    30,   130,    70}
    ,{   190,   190,   170,   160,   170}
    }
   ,{{   200,   190,   180,   200,   180}
    ,{   180,   180,   170,   160,   170}
    ,{   200,   190,   180,   200,   180}
    ,{   180,   180,   170,   160,   170}
    ,{   170,   170,   100,    90,   100}
    }
   }
  ,{{{   210,   210,   200,   200,   200}
    ,{   210,   210,   200,   160,   200}
    ,{   200,   190,   180,   200,   180}
    ,{   180,   180,   170,   130,   170}
    ,{   190,   190,   170,   190,   170}
    }
   ,{{   210,   210,   200,   160,   200}
    ,{   210,   210,   200,   160,   200}
    ,{   190,   190,   170,   130,   170}
    ,{    10,    10,     0,   -40,     0}
    ,{   190,   190,   170,   130,   170}
    }
   ,{{   190,   190,   170,   190,   170}
    ,{   180,   180,   170,   130,   170}
    ,{   190,   190,   170,   190,   170}
    ,{   180,   180,   170,   130,   170}
    ,{   190,   190,   170,   190,   170}
    }
   ,{{   190,   190,   170,   130,   170}
    ,{    70,    70,    60,    20,    60}
    ,{   190,   190,   170,   130,   170}
    ,{   120,    50,    30,   120,    30}
    ,{   190,   190,   170,   130,   170}
    }
   ,{{   200,   190,   180,   200,   180}
    ,{   180,   180,   170,   130,   170}
    ,{   200,   190,   180,   200,   180}
    ,{   180,   180,   170,   130,   170}
    ,{   170,   170,   100,    60,   100}
    }
   }
  ,{{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   170,   170,   170,   170,   170}
    ,{   160,   160,   160,   160,   160}
    ,{   170,   160,   170,   160,   170}
    }
   ,{{   190,   190,   190,   190,   190}
    ,{   190,   190,   190,   190,   190}
    ,{   170,   160,   170,   160,   170}
    ,{    50,   -10,    50,   -10,    50}
    ,{   170,   160,   170,   160,   170}
    }
   ,{{   170,   160,   170,   160,   170}
    ,{   160,   160,   160,   160,   160}
    ,{   170,   160,   170,   160,   170}
    ,{   160,   160,   160,   160,   160}
    ,{   170,   160,   170,   160,   170}
    }
   ,{{   170,   160,   170,   160,   170}
    ,{   110,    50,   110,    50,   110}
    ,{   170,   160,   170,   160,   170}
    ,{    30,    20,    30,    20,    30}
    ,{   170,   160,   170,   160,   170}
    }
   ,{{   170,   170,   170,   170,   170}
    ,{   160,   160,   160,   160,   160}
    ,{   170,   170,   170,   170,   170}
    ,{   160,   160,   160,   160,   160}
    ,{    90,    90,    90,    90,    90}
    }
   }
  ,{{{   200,   130,   200,   130,   200}
    ,{   200,    90,   200,    40,   200}
    ,{   180,   130,   180,    20,   180}
    ,{   170,    60,   170,   130,   170}
    ,{   170,   120,   170,    70,   170}
    }
   ,{{   200,    90,   200,    40,   200}
    ,{   200,    90,   200,    40,   200}
    ,{   170,    60,   170,    10,   170}
    ,{     0,  -110,     0,   -30,     0}
    ,{   170,    60,   170,    10,   170}
    }
   ,{{   170,   120,   170,    10,   170}
    ,{   170,    60,   170,    10,   170}
    ,{   170,   120,   170,    10,   170}
    ,{   170,    60,   170,    10,   170}
    ,{   170,   120,   170,    10,   170}
    }
   ,{{   170,    60,   170,   130,   170}
    ,{    60,   -50,    60,    30,    60}
    ,{   170,    60,   170,    10,   170}
    ,{   130,    50,    30,   130,    30}
    ,{   170,    60,   170,    10,   170}
    }
   ,{{   180,   130,   180,    70,   180}
    ,{   170,    60,   170,    10,   170}
    ,{   180,   130,   180,    20,   180}
    ,{   170,    60,   170,    10,   170}
    ,{   100,   -10,   100,    70,   100}
    }
   }
  ,{{{   190,   190,   190,   190,   160}
    ,{   190,   190,   190,   190,   160}
    ,{   170,   170,   170,   170,    80}
    ,{   160,   160,   160,   160,    70}
    ,{   170,   160,   170,   160,    80}
    }
   ,{{   190,   190,   190,   190,   160}
    ,{   190,   190,   190,   190,   160}
    ,{   170,   160,   170,   160,    80}
    ,{    50,   -10,    50,   -10,  -100}
    ,{   170,   160,   170,   160,    80}
    }
   ,{{   170,   160,   170,   160,    80}
    ,{   160,   160,   160,   160,    70}
    ,{   170,   160,   170,   160,    80}
    ,{   160,   160,   160,   160,    70}
    ,{   170,   160,   170,   160,    80}
    }
   ,{{   170,   160,   170,   160,    80}
    ,{   110,    50,   110,    50,   -30}
    ,{   170,   160,   170,   160,    80}
    ,{    70,    20,    30,    20,    70}
    ,{   170,   160,   170,   160,    80}
    }
   ,{{   170,   170,   170,   170,    80}
    ,{   160,   160,   160,   160,    70}
    ,{   170,   170,   170,   170,    80}
    ,{   160,   160,   160,   160,    70}
    ,{    90,    90,    90,    90,     0}
    }
   }
  }
 ,{{{{   370,   370,   330,   320,   330}
    ,{   340,   340,   330,   320,   330}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   370,   370,   290,   310,   290}
    }
   ,{{   340,   340,   330,   320,   330}
    ,{   340,   340,   330,   320,   330}
    ,{   310,   310,   290,   280,   290}
    ,{   270,   230,   270,   200,   270}
    ,{   310,   310,   290,   280,   290}
    }
   ,{{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   310,   310,   290,   310,   290}
    }
   ,{{   310,   310,   310,   280,   310}
    ,{   310,   270,   310,   240,   310}
    ,{   310,   310,   290,   280,   290}
    ,{   260,   180,   160,   260,   200}
    ,{   310,   310,   290,   280,   290}
    }
   ,{{   370,   370,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   370,   370,   290,   280,   290}
    }
   }
  ,{{{   370,   370,   330,   310,   330}
    ,{   340,   340,   330,   290,   330}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   370,   370,   290,   310,   290}
    }
   ,{{   340,   340,   330,   290,   330}
    ,{   340,   340,   330,   290,   330}
    ,{   310,   310,   290,   250,   290}
    ,{   230,   230,   210,   170,   210}
    ,{   310,   310,   290,   250,   290}
    }
   ,{{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   310,   290}
    }
   ,{{   310,   310,   290,   250,   290}
    ,{   270,   270,   250,   210,   250}
    ,{   310,   310,   290,   250,   290}
    ,{   250,   180,   160,   250,   160}
    ,{   310,   310,   290,   250,   290}
    }
   ,{{   370,   370,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   370,   370,   290,   250,   290}
    }
   }
  ,{{{   320,   320,   320,   320,   320}
    ,{   320,   320,   320,   320,   320}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   320,   320,   320,   320,   320}
    ,{   320,   320,   320,   320,   320}
    ,{   290,   280,   290,   280,   290}
    ,{   270,   200,   270,   200,   270}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   310,   280,   310,   280,   310}
    ,{   310,   240,   310,   240,   310}
    ,{   290,   280,   290,   280,   290}
    ,{   160,   150,   160,   150,   160}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    }
   }
  ,{{{   330,   240,   330,   260,   330}
    ,{   330,   220,   330,   220,   330}
    ,{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   260,   290}
    ,{   290,   240,   290,   260,   290}
    }
   ,{{   330,   220,   330,   180,   330}
    ,{   330,   220,   330,   170,   330}
    ,{   290,   180,   290,   130,   290}
    ,{   210,   100,   210,   180,   210}
    ,{   290,   180,   290,   130,   290}
    }
   ,{{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   240,   290,   130,   290}
    }
   ,{{   290,   180,   290,   260,   290}
    ,{   250,   140,   250,   220,   250}
    ,{   290,   180,   290,   130,   290}
    ,{   260,   180,   160,   260,   160}
    ,{   290,   180,   290,   130,   290}
    }
   ,{{   290,   240,   290,   260,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   180,   290,   260,   290}
    }
   }
  ,{{{   320,   320,   320,   320,   290}
    ,{   320,   320,   320,   320,   290}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   320,   320,   320,   320,   290}
    ,{   320,   320,   320,   320,   290}
    ,{   290,   280,   290,   280,   200}
    ,{   270,   200,   270,   200,   120}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   310,   280,   310,   280,   200}
    ,{   310,   240,   310,   240,   160}
    ,{   290,   280,   290,   280,   200}
    ,{   200,   150,   160,   150,   200}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    }
   }
  }
 ,{{{{   350,   340,   350,   280,   350}
    ,{   350,   310,   350,   280,   350}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   340,   340,   260,   280,   260}
    }
   ,{{   280,   280,   260,   250,   260}
    ,{   240,   240,   230,   220,   230}
    ,{   280,   280,   260,   250,   260}
    ,{   180,   140,   180,   120,   180}
    ,{   280,   280,   260,   250,   260}
    }
   ,{{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   280,   260}
    }
   ,{{   350,   310,   350,   280,   350}
    ,{   350,   310,   350,   280,   350}
    ,{   280,   280,   260,   250,   260}
    ,{   230,   150,   130,   230,   170}
    ,{   280,   280,   260,   250,   260}
    }
   ,{{   340,   340,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   340,   340,   260,   250,   260}
    }
   }
  ,{{{   340,   340,   290,   280,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   340,   340,   260,   280,   260}
    }
   ,{{   280,   280,   260,   220,   260}
    ,{   240,   240,   230,   190,   230}
    ,{   280,   280,   260,   220,   260}
    ,{   140,   140,   130,    90,   130}
    ,{   280,   280,   260,   220,   260}
    }
   ,{{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   280,   260}
    }
   ,{{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   280,   280,   260,   220,   260}
    ,{   220,   150,   130,   220,   130}
    ,{   280,   280,   260,   220,   260}
    }
   ,{{   340,   340,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   340,   340,   260,   220,   260}
    }
   }
  ,{{{   350,   280,   350,   280,   350}
    ,{   350,   280,   350,   280,   350}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    }
   ,{{   260,   250,   260,   250,   260}
    ,{   220,   220,   220,   220,   220}
    ,{   260,   250,   260,   250,   260}
    ,{   180,   120,   180,   120,   180}
    ,{   260,   250,   260,   250,   260}
    }
   ,{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    }
   ,{{   350,   280,   350,   280,   350}
    ,{   350,   280,   350,   280,   350}
    ,{   260,   250,   260,   250,   260}
    ,{   130,   120,   130,   120,   130}
    ,{   260,   250,   260,   250,   260}
    }
   ,{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    }
   }
  ,{{{   290,   210,   290,   260,   290}
    ,{   290,   180,   290,   260,   290}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   150,   260,   230,   260}
    ,{   260,   210,   260,   230,   260}
    }
   ,{{   260,   150,   260,   100,   260}
    ,{   230,   120,   230,    70,   230}
    ,{   260,   150,   260,   100,   260}
    ,{   130,    20,   130,   100,   130}
    ,{   260,   150,   260,   100,   260}
    }
   ,{{   260,   210,   260,   100,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   260,   210,   260,   100,   260}
    }
   ,{{   290,   180,   290,   260,   290}
    ,{   290,   180,   290,   260,   290}
    ,{   260,   150,   260,   100,   260}
    ,{   230,   150,   130,   230,   130}
    ,{   260,   150,   260,   100,   260}
    }
   ,{{   260,   210,   260,   230,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   260,   150,   260,   230,   260}
    }
   }
  ,{{{   350,   280,   350,   280,   200}
    ,{   350,   280,   350,   280,   200}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    }
   ,{{   260,   250,   260,   250,   190}
    ,{   220,   220,   220,   220,   190}
    ,{   260,   250,   260,   250,   170}
    ,{   180,   120,   180,   120,    30}
    ,{   260,   250,   260,   250,   170}
    }
   ,{{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    }
   ,{{   350,   280,   350,   280,   200}
    ,{   350,   280,   350,   280,   200}
    ,{   260,   250,   260,   250,   170}
    ,{   170,   120,   130,   120,   170}
    ,{   260,   250,   260,   250,   170}
    }
   ,{{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    }
   }
  }
 ,{{{{   280,   280,   260,   260,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   260,   260,   240,   260,   240}
    ,{   260,   260,   250,   240,   250}
    ,{   260,   260,   240,   260,   240}
    }
   ,{{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   250,   250,   240,   230,   240}
    ,{   190,   150,   190,   130,   190}
    ,{   250,   250,   240,   230,   240}
    }
   ,{{   260,   260,   250,   260,   250}
    ,{   260,   260,   250,   240,   250}
    ,{   260,   260,   240,   260,   240}
    ,{   260,   260,   250,   240,   250}
    ,{   260,   260,   240,   260,   240}
    }
   ,{{   260,   250,   260,   230,   260}
    ,{   260,   220,   260,   200,   260}
    ,{   250,   250,   240,   230,   240}
    ,{   190,   110,    90,   190,   120}
    ,{   250,   250,   240,   230,   240}
    }
   ,{{   260,   260,   250,   260,   250}
    ,{   260,   260,   250,   240,   250}
    ,{   260,   260,   240,   260,   240}
    ,{   260,   260,   250,   240,   250}
    ,{   230,   230,   150,   140,   150}
    }
   }
  ,{{{   280,   280,   260,   260,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   260,   260,   240,   260,   240}
    ,{   260,   260,   250,   210,   250}
    ,{   260,   260,   240,   260,   240}
    }
   ,{{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   250,   250,   240,   200,   240}
    ,{   150,   150,   140,   100,   140}
    ,{   250,   250,   240,   200,   240}
    }
   ,{{   260,   260,   250,   260,   250}
    ,{   260,   260,   250,   210,   250}
    ,{   260,   260,   240,   260,   240}
    ,{   260,   260,   250,   210,   250}
    ,{   260,   260,   240,   260,   240}
    }
   ,{{   250,   250,   240,   200,   240}
    ,{   220,   220,   210,   170,   210}
    ,{   250,   250,   240,   200,   240}
    ,{   180,   100,    90,   180,    90}
    ,{   250,   250,   240,   200,   240}
    }
   ,{{   260,   260,   250,   260,   250}
    ,{   260,   260,   250,   210,   250}
    ,{   260,   260,   240,   260,   240}
    ,{   260,   260,   250,   210,   250}
    ,{   230,   230,   150,   110,   150}
    }
   }
  ,{{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   240,   230,   240,   230,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   240,   230,   240,   230,   240}
    }
   ,{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   230,   230,   230,   230,   230}
    ,{   190,   130,   190,   130,   190}
    ,{   230,   230,   230,   230,   230}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   240,   230,   240,   230,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   240,   230,   240,   230,   240}
    }
   ,{{   260,   230,   260,   230,   260}
    ,{   260,   200,   260,   200,   260}
    ,{   230,   230,   230,   230,   230}
    ,{    80,    80,    80,    80,    80}
    ,{   230,   230,   230,   230,   230}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   240,   230,   240,   230,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   150,   140,   150,   140,   150}
    }
   }
  ,{{{   260,   190,   260,   190,   260}
    ,{   260,   150,   260,   180,   260}
    ,{   240,   190,   240,    80,   240}
    ,{   250,   140,   250,   190,   250}
    ,{   240,   190,   240,   120,   240}
    }
   ,{{   260,   150,   260,   110,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   240,   130,   240,    80,   240}
    ,{   140,    30,   140,   110,   140}
    ,{   240,   130,   240,    80,   240}
    }
   ,{{   250,   190,   250,    90,   250}
    ,{   250,   140,   250,    90,   250}
    ,{   240,   190,   240,    80,   240}
    ,{   250,   140,   250,    90,   250}
    ,{   240,   190,   240,    80,   240}
    }
   ,{{   240,   130,   240,   190,   240}
    ,{   210,   100,   210,   180,   210}
    ,{   240,   130,   240,    80,   240}
    ,{   190,   110,    90,   190,    90}
    ,{   240,   130,   240,    80,   240}
    }
   ,{{   250,   190,   250,   120,   250}
    ,{   250,   140,   250,    90,   250}
    ,{   240,   190,   240,    80,   240}
    ,{   250,   140,   250,    90,   250}
    ,{   150,    40,   150,   120,   150}
    }
   }
  ,{{{   260,   250,   260,   250,   230}
    ,{   260,   250,   260,   250,   230}
    ,{   240,   230,   240,   230,   150}
    ,{   240,   240,   240,   240,   150}
    ,{   240,   230,   240,   230,   150}
    }
   ,{{   260,   250,   260,   250,   230}
    ,{   260,   250,   260,   250,   230}
    ,{   230,   230,   230,   230,   140}
    ,{   190,   130,   190,   130,    40}
    ,{   230,   230,   230,   230,   140}
    }
   ,{{   240,   240,   240,   240,   150}
    ,{   240,   240,   240,   240,   150}
    ,{   240,   230,   240,   230,   150}
    ,{   240,   240,   240,   240,   150}
    ,{   240,   230,   240,   230,   150}
    }
   ,{{   260,   230,   260,   230,   140}
    ,{   260,   200,   260,   200,   110}
    ,{   230,   230,   230,   230,   140}
    ,{   120,    80,    80,    80,   120}
    ,{   230,   230,   230,   230,   140}
    }
   ,{{   240,   240,   240,   240,   150}
    ,{   240,   240,   240,   240,   150}
    ,{   240,   230,   240,   230,   150}
    ,{   240,   240,   240,   240,   150}
    ,{   150,   140,   150,   140,    60}
    }
   }
  }
 ,{{{{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   280,   260}
    }
   ,{{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   230,   230,   220,   210,   220}
    ,{   210,   170,   210,   150,   210}
    ,{   230,   230,   220,   210,   220}
    }
   ,{{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   280,   260}
    }
   ,{{   230,   230,   220,   210,   220}
    ,{   220,   180,   220,   160,   220}
    ,{   230,   230,   220,   210,   220}
    ,{   210,   130,   110,   210,   140}
    ,{   230,   230,   220,   210,   220}
    }
   ,{{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   250,   260}
    ,{   250,   250,   230,   250,   230}
    ,{   280,   280,   260,   250,   260}
    ,{   250,   250,   180,   170,   180}
    }
   }
  ,{{{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   280,   260}
    }
   ,{{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   230,   230,   220,   180,   220}
    ,{   170,   170,   160,   120,   160}
    ,{   230,   230,   220,   180,   220}
    }
   ,{{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   280,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   280,   280,   260,   280,   260}
    }
   ,{{   230,   230,   220,   200,   220}
    ,{   180,   180,   170,   130,   170}
    ,{   230,   230,   220,   180,   220}
    ,{   200,   120,   110,   200,   110}
    ,{   230,   230,   220,   180,   220}
    }
   ,{{   280,   280,   260,   250,   260}
    ,{   280,   280,   260,   220,   260}
    ,{   250,   250,   230,   250,   230}
    ,{   280,   280,   260,   220,   260}
    ,{   250,   250,   180,   140,   180}
    }
   }
  ,{{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    }
   ,{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   210,   210,   210,   210,   210}
    ,{   210,   150,   210,   150,   210}
    ,{   210,   210,   210,   210,   210}
    }
   ,{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    }
   ,{{   220,   210,   220,   210,   220}
    ,{   220,   160,   220,   160,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   100,   100,   100,   100,   100}
    ,{   210,   210,   210,   210,   210}
    }
   ,{{   260,   250,   260,   250,   260}
    ,{   260,   250,   260,   250,   260}
    ,{   230,   220,   230,   220,   230}
    ,{   260,   250,   260,   250,   260}
    ,{   170,   170,   170,   170,   170}
    }
   }
  ,{{{   260,   210,   260,   210,   260}
    ,{   260,   150,   260,   140,   260}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   150,   260,   210,   260}
    ,{   260,   210,   260,   150,   260}
    }
   ,{{   260,   150,   260,   130,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   220,   110,   220,    60,   220}
    ,{   160,    50,   160,   130,   160}
    ,{   220,   110,   220,    60,   220}
    }
   ,{{   260,   210,   260,   100,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   260,   210,   260,   100,   260}
    }
   ,{{   220,   130,   220,   210,   220}
    ,{   170,    60,   170,   140,   170}
    ,{   220,   110,   220,    60,   220}
    ,{   210,   130,   110,   210,   110}
    ,{   220,   110,   220,    60,   220}
    }
   ,{{   260,   180,   260,   150,   260}
    ,{   260,   150,   260,   100,   260}
    ,{   230,   180,   230,    70,   230}
    ,{   260,   150,   260,   100,   260}
    ,{   180,    70,   180,   150,   180}
    }
   }
  ,{{{   260,   250,   260,   250,   230}
    ,{   260,   250,   260,   250,   230}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    }
   ,{{   260,   250,   260,   250,   230}
    ,{   260,   250,   260,   250,   230}
    ,{   210,   210,   210,   210,   120}
    ,{   210,   150,   210,   150,    60}
    ,{   210,   210,   210,   210,   120}
    }
   ,{{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    }
   ,{{   220,   210,   220,   210,   140}
    ,{   220,   160,   220,   160,    70}
    ,{   210,   210,   210,   210,   120}
    ,{   140,   100,   100,   100,   140}
    ,{   210,   210,   210,   210,   120}
    }
   ,{{   260,   250,   260,   250,   170}
    ,{   260,   250,   260,   250,   170}
    ,{   230,   220,   230,   220,   140}
    ,{   260,   250,   260,   250,   170}
    ,{   170,   170,   170,   170,    80}
    }
   }
  }
 ,{{{{   370,   370,   350,   320,   350}
    ,{   350,   340,   350,   320,   350}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   370,   370,   290,   310,   290}
    }
   ,{{   340,   340,   330,   320,   330}
    ,{   340,   340,   330,   320,   330}
    ,{   310,   310,   290,   280,   290}
    ,{   270,   230,   270,   200,   270}
    ,{   310,   310,   290,   280,   290}
    }
   ,{{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   310,   310,   290,   310,   290}
    }
   ,{{   350,   310,   350,   280,   350}
    ,{   350,   310,   350,   280,   350}
    ,{   310,   310,   290,   280,   290}
    ,{   260,   180,   160,   260,   200}
    ,{   310,   310,   290,   280,   290}
    }
   ,{{   370,   370,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   280,   290}
    ,{   370,   370,   290,   280,   290}
    }
   }
  ,{{{   370,   370,   330,   310,   330}
    ,{   340,   340,   330,   290,   330}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   370,   370,   290,   310,   290}
    }
   ,{{   340,   340,   330,   290,   330}
    ,{   340,   340,   330,   290,   330}
    ,{   310,   310,   290,   250,   290}
    ,{   230,   230,   210,   170,   210}
    ,{   310,   310,   290,   250,   290}
    }
   ,{{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   310,   290}
    }
   ,{{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   250,   180,   160,   250,   160}
    ,{   310,   310,   290,   250,   290}
    }
   ,{{   370,   370,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   310,   310,   290,   310,   290}
    ,{   310,   310,   290,   250,   290}
    ,{   370,   370,   290,   250,   290}
    }
   }
  ,{{{   350,   320,   350,   320,   350}
    ,{   350,   320,   350,   320,   350}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   320,   320,   320,   320,   320}
    ,{   320,   320,   320,   320,   320}
    ,{   290,   280,   290,   280,   290}
    ,{   270,   200,   270,   200,   270}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   350,   280,   350,   280,   350}
    ,{   350,   280,   350,   280,   350}
    ,{   290,   280,   290,   280,   290}
    ,{   160,   150,   160,   150,   160}
    ,{   290,   280,   290,   280,   290}
    }
   ,{{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    ,{   290,   280,   290,   280,   290}
    }
   }
  ,{{{   330,   240,   330,   260,   330}
    ,{   330,   220,   330,   260,   330}
    ,{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   260,   290}
    ,{   290,   240,   290,   260,   290}
    }
   ,{{   330,   220,   330,   180,   330}
    ,{   330,   220,   330,   170,   330}
    ,{   290,   180,   290,   130,   290}
    ,{   210,   100,   210,   180,   210}
    ,{   290,   180,   290,   130,   290}
    }
   ,{{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   240,   290,   130,   290}
    }
   ,{{   290,   180,   290,   260,   290}
    ,{   290,   180,   290,   260,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   260,   180,   160,   260,   160}
    ,{   290,   180,   290,   130,   290}
    }
   ,{{   290,   240,   290,   260,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   240,   290,   130,   290}
    ,{   290,   180,   290,   130,   290}
    ,{   290,   180,   290,   260,   290}
    }
   }
  ,{{{   350,   320,   350,   320,   290}
    ,{   350,   320,   350,   320,   290}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   320,   320,   320,   320,   290}
    ,{   320,   320,   320,   320,   290}
    ,{   290,   280,   290,   280,   200}
    ,{   270,   200,   270,   200,   120}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   350,   280,   350,   280,   200}
    ,{   350,   280,   350,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   200,   150,   160,   150,   200}
    ,{   290,   280,   290,   280,   200}
    }
   ,{{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    ,{   290,   280,   290,   280,   200}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   240,   240,   240,   190,   240}
    ,{   240,   240,   240,   190,   240}
    ,{   220,   220,   220,   190,   220}
    ,{   240,   240,   240,   190,   240}
    ,{   210,   210,   210,   170,   210}
    }
   ,{{   200,   200,   200,   150,   200}
    ,{   200,   200,   200,   150,   200}
    ,{   190,   190,   190,   150,   190}
    ,{   160,   100,   160,    80,   130}
    ,{   190,   190,   190,   150,   190}
    }
   ,{{   240,   240,   240,   190,   240}
    ,{   240,   240,   240,   190,   240}
    ,{   220,   220,   220,   190,   220}
    ,{   240,   240,   240,   190,   240}
    ,{   210,   210,   210,   170,   210}
    }
   ,{{   190,   190,   190,   150,   190}
    ,{   160,   100,   160,    80,   130}
    ,{   190,   190,   190,   150,   190}
    ,{   150,    70,    50,   150,    90}
    ,{   190,   190,   190,   150,   190}
    }
   ,{{   240,   240,   240,   190,   240}
    ,{   240,   240,   240,   190,   240}
    ,{   210,   210,   210,   170,   210}
    ,{   240,   240,   240,   190,   240}
    ,{   180,   180,   120,    90,   120}
    }
   }
  ,{{{   240,   240,   240,   190,   240}
    ,{   240,   240,   240,   140,   240}
    ,{   220,   220,   220,   190,   220}
    ,{   240,   240,   240,   140,   240}
    ,{   210,   210,   210,   170,   210}
    }
   ,{{   200,   200,   200,   100,   200}
    ,{   200,   200,   200,   100,   200}
    ,{   190,   190,   190,   100,   190}
    ,{   100,   100,   100,    10,   100}
    ,{   190,   190,   190,   100,   190}
    }
   ,{{   240,   240,   240,   190,   240}
    ,{   240,   240,   240,   140,   240}
    ,{   220,   220,   220,   190,   220}
    ,{   240,   240,   240,   140,   240}
    ,{   210,   210,   210,   170,   210}
    }
   ,{{   190,   190,   190,   100,   190}
    ,{   100,   100,   100,    10,   100}
    ,{   190,   190,   190,   100,   190}
    ,{    80,    50,    50,    80,    50}
    ,{   190,   190,   190,   100,   190}
    }
   ,{{   240,   240,   240,   170,   240}
    ,{   240,   240,   240,   140,   240}
    ,{   210,   210,   210,   170,   210}
    ,{   240,   240,   240,   140,   240}
    ,{   180,   180,   120,    20,   120}
    }
   }
  ,{{{   240,   190,   240,   190,   210}
    ,{   240,   190,   240,   190,   210}
    ,{   220,   180,   220,   180,   190}
    ,{   240,   190,   240,   190,   210}
    ,{   210,   160,   210,   160,   180}
    }
   ,{{   200,   150,   200,   150,   170}
    ,{   200,   150,   200,   150,   170}
    ,{   190,   150,   190,   150,   160}
    ,{   160,    60,   160,    60,   130}
    ,{   190,   150,   190,   150,   160}
    }
   ,{{   240,   190,   240,   190,   210}
    ,{   240,   190,   240,   190,   210}
    ,{   220,   180,   220,   180,   190}
    ,{   240,   190,   240,   190,   210}
    ,{   210,   160,   210,   160,   180}
    }
   ,{{   190,   150,   190,   150,   160}
    ,{   160,    60,   160,    60,   130}
    ,{   190,   150,   190,   150,   160}
    ,{    50,     0,    50,     0,    20}
    ,{   190,   150,   190,   150,   160}
    }
   ,{{   240,   190,   240,   190,   210}
    ,{   240,   190,   240,   190,   210}
    ,{   210,   160,   210,   160,   180}
    ,{   240,   190,   240,   190,   210}
    ,{   120,    70,   120,    70,    90}
    }
   }
  ,{{{   240,   180,   240,   150,   240}
    ,{   240,   130,   240,    80,   240}
    ,{   220,   180,   220,    70,   220}
    ,{   240,   130,   240,   150,   240}
    ,{   210,   160,   210,    90,   210}
    }
   ,{{   200,    90,   200,    80,   200}
    ,{   200,    90,   200,    40,   200}
    ,{   190,    90,   190,    40,   190}
    ,{   100,     0,   100,    80,   100}
    ,{   190,    90,   190,    40,   190}
    }
   ,{{   240,   180,   240,    80,   240}
    ,{   240,   130,   240,    80,   240}
    ,{   220,   180,   220,    70,   220}
    ,{   240,   130,   240,    80,   240}
    ,{   210,   160,   210,    50,   210}
    }
   ,{{   190,    90,   190,   150,   190}
    ,{   100,     0,   100,    80,   100}
    ,{   190,    90,   190,    40,   190}
    ,{   150,    70,    50,   150,    50}
    ,{   190,    90,   190,    40,   190}
    }
   ,{{   240,   160,   240,    90,   240}
    ,{   240,   130,   240,    80,   240}
    ,{   210,   160,   210,    50,   210}
    ,{   240,   130,   240,    80,   240}
    ,{   120,    10,   120,    90,   120}
    }
   }
  ,{{{   240,   190,   240,   190,   170}
    ,{   240,   190,   240,   190,   170}
    ,{   220,   180,   220,   180,   140}
    ,{   240,   190,   240,   190,   150}
    ,{   210,   160,   210,   160,   120}
    }
   ,{{   200,   150,   200,   150,   170}
    ,{   200,   150,   200,   150,   170}
    ,{   190,   150,   190,   150,   110}
    ,{   160,    60,   160,    60,    20}
    ,{   190,   150,   190,   150,   110}
    }
   ,{{   240,   190,   240,   190,   150}
    ,{   240,   190,   240,   190,   150}
    ,{   220,   180,   220,   180,   140}
    ,{   240,   190,   240,   190,   150}
    ,{   210,   160,   210,   160,   120}
    }
   ,{{   190,   150,   190,   150,   110}
    ,{   160,    60,   160,    60,    20}
    ,{   190,   150,   190,   150,   110}
    ,{    90,     0,    50,     0,    90}
    ,{   190,   150,   190,   150,   110}
    }
   ,{{   240,   190,   240,   190,   150}
    ,{   240,   190,   240,   190,   150}
    ,{   210,   160,   210,   160,   120}
    ,{   240,   190,   240,   190,   150}
    ,{   120,    70,   120,    70,    30}
    }
   }
  }
 ,{{{{   210,   210,   210,   170,   210}
    ,{   210,   210,   210,   170,   210}
    ,{   190,   190,   190,   160,   190}
    ,{   180,   180,   180,   150,   180}
    ,{   190,   190,   190,   150,   190}
    }
   ,{{   210,   210,   210,   170,   210}
    ,{   210,   210,   210,   170,   210}
    ,{   190,   190,   190,   140,   190}
    ,{    70,    10,    70,   -10,    40}
    ,{   190,   190,   190,   140,   190}
    }
   ,{{   190,   190,   190,   150,   190}
    ,{   180,   180,   180,   140,   180}
    ,{   190,   190,   190,   150,   190}
    ,{   180,   180,   180,   140,   180}
    ,{   190,   190,   190,   150,   190}
    }
   ,{{   190,   190,   190,   150,   190}
    ,{   130,    70,   130,    50,   100}
    ,{   190,   190,   190,   140,   190}
    ,{   150,    70,    50,   150,    90}
    ,{   190,   190,   190,   140,   190}
    }
   ,{{   190,   190,   190,   160,   190}
    ,{   180,   180,   180,   140,   180}
    ,{   190,   190,   190,   160,   190}
    ,{   180,   180,   180,   140,   180}
    ,{   170,   170,   110,    90,   110}
    }
   }
  ,{{{   210,   210,   210,   160,   210}
    ,{   210,   210,   210,   120,   210}
    ,{   190,   190,   190,   160,   190}
    ,{   180,   180,   180,    90,   180}
    ,{   190,   190,   190,   150,   190}
    }
   ,{{   210,   210,   210,   120,   210}
    ,{   210,   210,   210,   120,   210}
    ,{   190,   190,   190,    90,   190}
    ,{    10,    10,    10,   -80,    10}
    ,{   190,   190,   190,    90,   190}
    }
   ,{{   190,   190,   190,   150,   190}
    ,{   180,   180,   180,    90,   180}
    ,{   190,   190,   190,   150,   190}
    ,{   180,   180,   180,    90,   180}
    ,{   190,   190,   190,   150,   190}
    }
   ,{{   190,   190,   190,    90,   190}
    ,{    70,    70,    70,   -20,    70}
    ,{   190,   190,   190,    90,   190}
    ,{    80,    50,    50,    80,    50}
    ,{   190,   190,   190,    90,   190}
    }
   ,{{   190,   190,   190,   160,   190}
    ,{   180,   180,   180,    90,   180}
    ,{   190,   190,   190,   160,   190}
    ,{   180,   180,   180,    90,   180}
    ,{   170,   170,   110,    20,   110}
    }
   }
  ,{{{   210,   170,   210,   170,   180}
    ,{   210,   170,   210,   170,   180}
    ,{   190,   150,   190,   150,   160}
    ,{   180,   140,   180,   140,   150}
    ,{   190,   140,   190,   140,   160}
    }
   ,{{   210,   170,   210,   170,   180}
    ,{   210,   170,   210,   170,   180}
    ,{   190,   140,   190,   140,   160}
    ,{    70,   -30,    70,   -30,    40}
    ,{   190,   140,   190,   140,   160}
    }
   ,{{   190,   140,   190,   140,   160}
    ,{   180,   140,   180,   140,   150}
    ,{   190,   140,   190,   140,   160}
    ,{   180,   140,   180,   140,   150}
    ,{   190,   140,   190,   140,   160}
    }
   ,{{   190,   140,   190,   140,   160}
    ,{   130,    30,   130,    30,   100}
    ,{   190,   140,   190,   140,   160}
    ,{    50,     0,    50,     0,    20}
    ,{   190,   140,   190,   140,   160}
    }
   ,{{   190,   150,   190,   150,   160}
    ,{   180,   140,   180,   140,   150}
    ,{   190,   150,   190,   150,   160}
    ,{   180,   140,   180,   140,   150}
    ,{   110,    70,   110,    70,    80}
    }
   }
  ,{{{   210,   150,   210,   150,   210}
    ,{   210,   110,   210,    60,   210}
    ,{   190,   150,   190,    40,   190}
    ,{   180,    80,   180,   150,   180}
    ,{   190,   140,   190,    90,   190}
    }
   ,{{   210,   110,   210,    60,   210}
    ,{   210,   110,   210,    60,   210}
    ,{   190,    80,   190,    30,   190}
    ,{    10,   -90,    10,   -10,    10}
    ,{   190,    80,   190,    30,   190}
    }
   ,{{   190,   140,   190,    30,   190}
    ,{   180,    80,   180,    30,   180}
    ,{   190,   140,   190,    30,   190}
    ,{   180,    80,   180,    30,   180}
    ,{   190,   140,   190,    30,   190}
    }
   ,{{   190,    80,   190,   150,   190}
    ,{    70,   -30,    70,    50,    70}
    ,{   190,    80,   190,    30,   190}
    ,{   150,    70,    50,   150,    50}
    ,{   190,    80,   190,    30,   190}
    }
   ,{{   190,   150,   190,    90,   190}
    ,{   180,    80,   180,    30,   180}
    ,{   190,   150,   190,    40,   190}
    ,{   180,    80,   180,    30,   180}
    ,{   110,    10,   110,    90,   110}
    }
   }
  ,{{{   210,   170,   210,   170,   190}
    ,{   210,   170,   210,   170,   190}
    ,{   190,   150,   190,   150,   110}
    ,{   180,   140,   180,   140,   100}
    ,{   190,   140,   190,   140,   100}
    }
   ,{{   210,   170,   210,   170,   190}
    ,{   210,   170,   210,   170,   190}
    ,{   190,   140,   190,   140,   100}
    ,{    70,   -30,    70,   -30,   -70}
    ,{   190,   140,   190,   140,   100}
    }
   ,{{   190,   140,   190,   140,   100}
    ,{   180,   140,   180,   140,   100}
    ,{   190,   140,   190,   140,   100}
    ,{   180,   140,   180,   140,   100}
    ,{   190,   140,   190,   140,   100}
    }
   ,{{   190,   140,   190,   140,   100}
    ,{   130,    30,   130,    30,   -10}
    ,{   190,   140,   190,   140,   100}
    ,{    90,     0,    50,     0,    90}
    ,{   190,   140,   190,   140,   100}
    }
   ,{{   190,   150,   190,   150,   110}
    ,{   180,   140,   180,   140,   100}
    ,{   190,   150,   190,   150,   110}
    ,{   180,   140,   180,   140,   100}
    ,{   110,    70,   110,    70,    30}
    }
   }
  }
 ,{{{{   370,   370,   340,   300,   340}
    ,{   340,   340,   340,   300,   340}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   280,   310}
    ,{   370,   370,   310,   280,   310}
    }
   ,{{   340,   340,   340,   300,   340}
    ,{   340,   340,   340,   300,   340}
    ,{   310,   310,   310,   260,   310}
    ,{   290,   230,   290,   200,   260}
    ,{   310,   310,   310,   260,   310}
    }
   ,{{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   310,   310,   310,   270,   310}
    }
   ,{{   330,   310,   330,   280,   310}
    ,{   330,   270,   330,   240,   300}
    ,{   310,   310,   310,   260,   310}
    ,{   280,   200,   180,   280,   220}
    ,{   310,   310,   310,   260,   310}
    }
   ,{{   370,   370,   310,   280,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   370,   370,   310,   280,   310}
    }
   }
  ,{{{   370,   370,   340,   270,   340}
    ,{   340,   340,   340,   250,   340}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   370,   370,   310,   270,   310}
    }
   ,{{   340,   340,   340,   250,   340}
    ,{   340,   340,   340,   250,   340}
    ,{   310,   310,   310,   210,   310}
    ,{   230,   230,   230,   130,   230}
    ,{   310,   310,   310,   210,   310}
    }
   ,{{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   270,   310}
    }
   ,{{   310,   310,   310,   210,   310}
    ,{   270,   270,   270,   170,   270}
    ,{   310,   310,   310,   210,   310}
    ,{   210,   180,   180,   210,   180}
    ,{   310,   310,   310,   210,   310}
    }
   ,{{   370,   370,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   370,   370,   310,   210,   310}
    }
   }
  ,{{{   340,   300,   340,   300,   310}
    ,{   340,   300,   340,   300,   310}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   340,   300,   340,   300,   310}
    ,{   340,   300,   340,   300,   310}
    ,{   310,   260,   310,   260,   280}
    ,{   290,   180,   290,   180,   260}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   330,   260,   330,   260,   300}
    ,{   330,   220,   330,   220,   300}
    ,{   310,   260,   310,   260,   280}
    ,{   180,   130,   180,   130,   150}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    }
   }
  ,{{{   340,   260,   340,   280,   340}
    ,{   340,   240,   340,   240,   340}
    ,{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   280,   310}
    ,{   310,   260,   310,   280,   310}
    }
   ,{{   340,   240,   340,   200,   340}
    ,{   340,   240,   340,   190,   340}
    ,{   310,   200,   310,   150,   310}
    ,{   230,   120,   230,   200,   230}
    ,{   310,   200,   310,   150,   310}
    }
   ,{{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   260,   310,   150,   310}
    }
   ,{{   310,   200,   310,   280,   310}
    ,{   270,   160,   270,   240,   270}
    ,{   310,   200,   310,   150,   310}
    ,{   280,   200,   180,   280,   180}
    ,{   310,   200,   310,   150,   310}
    }
   ,{{   310,   260,   310,   280,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   200,   310,   280,   310}
    }
   }
  ,{{{   340,   300,   340,   300,   320}
    ,{   340,   300,   340,   300,   320}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   340,   300,   340,   300,   320}
    ,{   340,   300,   340,   300,   320}
    ,{   310,   260,   310,   260,   220}
    ,{   290,   180,   290,   180,   140}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   330,   260,   330,   260,   220}
    ,{   330,   220,   330,   220,   180}
    ,{   310,   260,   310,   260,   220}
    ,{   220,   130,   180,   130,   220}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    }
   }
  }
 ,{{{{   370,   340,   370,   280,   340}
    ,{   370,   310,   370,   280,   340}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   250,   280}
    ,{   340,   340,   280,   250,   280}
    }
   ,{{   280,   280,   280,   230,   280}
    ,{   240,   240,   240,   200,   240}
    ,{   280,   280,   280,   230,   280}
    ,{   200,   140,   200,   120,   170}
    ,{   280,   280,   280,   230,   280}
    }
   ,{{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   240,   280}
    }
   ,{{   370,   310,   370,   280,   340}
    ,{   370,   310,   370,   280,   340}
    ,{   280,   280,   280,   230,   280}
    ,{   250,   170,   150,   250,   190}
    ,{   280,   280,   280,   230,   280}
    }
   ,{{   340,   340,   280,   250,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   340,   340,   280,   250,   280}
    }
   }
  ,{{{   340,   340,   310,   240,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   340,   340,   280,   240,   280}
    }
   ,{{   280,   280,   280,   180,   280}
    ,{   240,   240,   240,   150,   240}
    ,{   280,   280,   280,   180,   280}
    ,{   140,   140,   140,    50,   140}
    ,{   280,   280,   280,   180,   280}
    }
   ,{{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   240,   280}
    }
   ,{{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   280,   280,   280,   180,   280}
    ,{   180,   150,   150,   180,   150}
    ,{   280,   280,   280,   180,   280}
    }
   ,{{   340,   340,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   340,   340,   280,   180,   280}
    }
   }
  ,{{{   370,   260,   370,   260,   340}
    ,{   370,   260,   370,   260,   340}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   240,   200,   240,   200,   210}
    ,{   280,   230,   280,   230,   250}
    ,{   200,   100,   200,   100,   170}
    ,{   280,   230,   280,   230,   250}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    }
   ,{{   370,   260,   370,   260,   340}
    ,{   370,   260,   370,   260,   340}
    ,{   280,   230,   280,   230,   250}
    ,{   150,   100,   150,   100,   120}
    ,{   280,   230,   280,   230,   250}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    }
   }
  ,{{{   310,   230,   310,   280,   310}
    ,{   310,   200,   310,   280,   310}
    ,{   280,   230,   280,   120,   280}
    ,{   280,   170,   280,   250,   280}
    ,{   280,   230,   280,   250,   280}
    }
   ,{{   280,   170,   280,   120,   280}
    ,{   240,   140,   240,    90,   240}
    ,{   280,   170,   280,   120,   280}
    ,{   140,    40,   140,   120,   140}
    ,{   280,   170,   280,   120,   280}
    }
   ,{{   280,   230,   280,   120,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   280,   230,   280,   120,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   280,   230,   280,   120,   280}
    }
   ,{{   310,   200,   310,   280,   310}
    ,{   310,   200,   310,   280,   310}
    ,{   280,   170,   280,   120,   280}
    ,{   250,   170,   150,   250,   150}
    ,{   280,   170,   280,   120,   280}
    }
   ,{{   280,   230,   280,   250,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   280,   230,   280,   120,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   280,   170,   280,   250,   280}
    }
   }
  ,{{{   370,   260,   370,   260,   220}
    ,{   370,   260,   370,   260,   220}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    }
   ,{{   280,   230,   280,   230,   220}
    ,{   240,   200,   240,   200,   220}
    ,{   280,   230,   280,   230,   190}
    ,{   200,   100,   200,   100,    60}
    ,{   280,   230,   280,   230,   190}
    }
   ,{{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    }
   ,{{   370,   260,   370,   260,   220}
    ,{   370,   260,   370,   260,   220}
    ,{   280,   230,   280,   230,   190}
    ,{   190,   100,   150,   100,   190}
    ,{   280,   230,   280,   230,   190}
    }
   ,{{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    }
   }
  }
 ,{{{{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    }
   ,{{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   250,   250,   250,   210,   250}
    ,{   210,   150,   210,   130,   180}
    ,{   250,   250,   250,   210,   250}
    }
   ,{{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    }
   ,{{   280,   250,   280,   210,   250}
    ,{   280,   220,   280,   200,   250}
    ,{   250,   250,   250,   210,   250}
    ,{   210,   130,   100,   210,   150}
    ,{   250,   250,   250,   210,   250}
    }
   ,{{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   230,   230,   170,   140,   170}
    }
   }
  ,{{{   280,   280,   280,   220,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   170,   260}
    ,{   260,   260,   260,   220,   260}
    }
   ,{{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   250,   250,   250,   160,   250}
    ,{   150,   150,   150,    60,   150}
    ,{   250,   250,   250,   160,   250}
    }
   ,{{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   170,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   170,   260}
    ,{   260,   260,   260,   220,   260}
    }
   ,{{   250,   250,   250,   160,   250}
    ,{   220,   220,   220,   130,   220}
    ,{   250,   250,   250,   160,   250}
    ,{   140,   100,   100,   140,   100}
    ,{   250,   250,   250,   160,   250}
    }
   ,{{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   170,   260}
    ,{   260,   260,   260,   220,   260}
    ,{   260,   260,   260,   170,   260}
    ,{   230,   230,   170,    70,   170}
    }
   }
  ,{{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   260,   210,   260,   210,   230}
    ,{   260,   220,   260,   220,   230}
    ,{   260,   210,   260,   210,   230}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   250,   210,   250,   210,   220}
    ,{   210,   110,   210,   110,   180}
    ,{   250,   210,   250,   210,   220}
    }
   ,{{   260,   220,   260,   220,   230}
    ,{   260,   220,   260,   220,   230}
    ,{   260,   210,   260,   210,   230}
    ,{   260,   220,   260,   220,   230}
    ,{   260,   210,   260,   210,   230}
    }
   ,{{   280,   210,   280,   210,   250}
    ,{   280,   180,   280,   180,   250}
    ,{   250,   210,   250,   210,   220}
    ,{   100,    60,   100,    60,    70}
    ,{   250,   210,   250,   210,   220}
    }
   ,{{   260,   220,   260,   220,   230}
    ,{   260,   220,   260,   220,   230}
    ,{   260,   210,   260,   210,   230}
    ,{   260,   220,   260,   220,   230}
    ,{   170,   120,   170,   120,   140}
    }
   }
  ,{{{   280,   210,   280,   210,   280}
    ,{   280,   170,   280,   200,   280}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   160,   260,   210,   260}
    ,{   260,   210,   260,   140,   260}
    }
   ,{{   280,   170,   280,   130,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   250,   150,   250,   100,   250}
    ,{   150,    50,   150,   130,   150}
    ,{   250,   150,   250,   100,   250}
    }
   ,{{   260,   210,   260,   110,   260}
    ,{   260,   160,   260,   110,   260}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   160,   260,   110,   260}
    ,{   260,   210,   260,   100,   260}
    }
   ,{{   250,   150,   250,   210,   250}
    ,{   220,   120,   220,   200,   220}
    ,{   250,   150,   250,   100,   250}
    ,{   210,   130,   100,   210,   100}
    ,{   250,   150,   250,   100,   250}
    }
   ,{{   260,   210,   260,   140,   260}
    ,{   260,   160,   260,   110,   260}
    ,{   260,   210,   260,   100,   260}
    ,{   260,   160,   260,   110,   260}
    ,{   170,    60,   170,   140,   170}
    }
   }
  ,{{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   260,   210,   260,   210,   170}
    ,{   260,   220,   260,   220,   180}
    ,{   260,   210,   260,   210,   170}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   250,   210,   250,   210,   170}
    ,{   210,   110,   210,   110,    70}
    ,{   250,   210,   250,   210,   170}
    }
   ,{{   260,   220,   260,   220,   180}
    ,{   260,   220,   260,   220,   180}
    ,{   260,   210,   260,   210,   170}
    ,{   260,   220,   260,   220,   180}
    ,{   260,   210,   260,   210,   170}
    }
   ,{{   280,   210,   280,   210,   170}
    ,{   280,   180,   280,   180,   140}
    ,{   250,   210,   250,   210,   170}
    ,{   150,    60,   100,    60,   150}
    ,{   250,   210,   250,   210,   170}
    }
   ,{{   260,   220,   260,   220,   180}
    ,{   260,   220,   260,   220,   180}
    ,{   260,   210,   260,   210,   170}
    ,{   260,   220,   260,   220,   180}
    ,{   170,   120,   170,   120,    80}
    }
   }
  }
 ,{{{{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   240,   280}
    }
   ,{{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   230,   230,   230,   190,   230}
    ,{   230,   170,   230,   150,   200}
    ,{   230,   230,   230,   190,   230}
    }
   ,{{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   240,   280}
    }
   ,{{   240,   230,   240,   230,   230}
    ,{   240,   180,   240,   160,   210}
    ,{   230,   230,   230,   190,   230}
    ,{   230,   150,   120,   230,   170}
    ,{   230,   230,   230,   190,   230}
    }
   ,{{   280,   280,   280,   230,   280}
    ,{   280,   280,   280,   230,   280}
    ,{   250,   250,   250,   210,   250}
    ,{   280,   280,   280,   230,   280}
    ,{   250,   250,   190,   170,   190}
    }
   }
  ,{{{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   240,   280}
    }
   ,{{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   230,   230,   230,   140,   230}
    ,{   170,   170,   170,    80,   170}
    ,{   230,   230,   230,   140,   230}
    }
   ,{{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   240,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   280,   280,   280,   240,   280}
    }
   ,{{   230,   230,   230,   160,   230}
    ,{   180,   180,   180,    90,   180}
    ,{   230,   230,   230,   140,   230}
    ,{   160,   120,   120,   160,   120}
    ,{   230,   230,   230,   140,   230}
    }
   ,{{   280,   280,   280,   210,   280}
    ,{   280,   280,   280,   180,   280}
    ,{   250,   250,   250,   210,   250}
    ,{   280,   280,   280,   180,   280}
    ,{   250,   250,   190,   100,   190}
    }
   }
  ,{{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   230,   190,   230,   190,   200}
    ,{   230,   130,   230,   130,   200}
    ,{   230,   190,   230,   190,   200}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    }
   ,{{   240,   190,   240,   190,   210}
    ,{   240,   140,   240,   140,   210}
    ,{   230,   190,   230,   190,   200}
    ,{   120,    80,   120,    80,    90}
    ,{   230,   190,   230,   190,   200}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   250,   200,   250,   200,   220}
    ,{   280,   230,   280,   230,   250}
    ,{   190,   150,   190,   150,   160}
    }
   }
  ,{{{   280,   230,   280,   230,   280}
    ,{   280,   170,   280,   160,   280}
    ,{   280,   230,   280,   120,   280}
    ,{   280,   170,   280,   230,   280}
    ,{   280,   230,   280,   170,   280}
    }
   ,{{   280,   170,   280,   150,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   230,   130,   230,    80,   230}
    ,{   170,    70,   170,   150,   170}
    ,{   230,   130,   230,    80,   230}
    }
   ,{{   280,   230,   280,   120,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   280,   230,   280,   120,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   280,   230,   280,   120,   280}
    }
   ,{{   230,   150,   230,   230,   230}
    ,{   180,    80,   180,   160,   180}
    ,{   230,   130,   230,    80,   230}
    ,{   230,   150,   120,   230,   120}
    ,{   230,   130,   230,    80,   230}
    }
   ,{{   280,   200,   280,   170,   280}
    ,{   280,   170,   280,   120,   280}
    ,{   250,   200,   250,    90,   250}
    ,{   280,   170,   280,   120,   280}
    ,{   190,    90,   190,   170,   190}
    }
   }
  ,{{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    }
   ,{{   280,   230,   280,   230,   250}
    ,{   280,   230,   280,   230,   250}
    ,{   230,   190,   230,   190,   150}
    ,{   230,   130,   230,   130,    90}
    ,{   230,   190,   230,   190,   150}
    }
   ,{{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    }
   ,{{   240,   190,   240,   190,   170}
    ,{   240,   140,   240,   140,   100}
    ,{   230,   190,   230,   190,   150}
    ,{   170,    80,   120,    80,   170}
    ,{   230,   190,   230,   190,   150}
    }
   ,{{   280,   230,   280,   230,   190}
    ,{   280,   230,   280,   230,   190}
    ,{   250,   200,   250,   200,   160}
    ,{   280,   230,   280,   230,   190}
    ,{   190,   150,   190,   150,   110}
    }
   }
  }
 ,{{{{   370,   370,   370,   300,   340}
    ,{   370,   340,   370,   300,   340}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   280,   310}
    ,{   370,   370,   310,   280,   310}
    }
   ,{{   340,   340,   340,   300,   340}
    ,{   340,   340,   340,   300,   340}
    ,{   310,   310,   310,   260,   310}
    ,{   290,   230,   290,   200,   260}
    ,{   310,   310,   310,   260,   310}
    }
   ,{{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   310,   310,   310,   270,   310}
    }
   ,{{   370,   310,   370,   280,   340}
    ,{   370,   310,   370,   280,   340}
    ,{   310,   310,   310,   260,   310}
    ,{   280,   200,   180,   280,   220}
    ,{   310,   310,   310,   260,   310}
    }
   ,{{   370,   370,   310,   280,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   260,   310}
    ,{   370,   370,   310,   280,   310}
    }
   }
  ,{{{   370,   370,   340,   270,   340}
    ,{   340,   340,   340,   250,   340}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   370,   370,   310,   270,   310}
    }
   ,{{   340,   340,   340,   250,   340}
    ,{   340,   340,   340,   250,   340}
    ,{   310,   310,   310,   210,   310}
    ,{   230,   230,   230,   130,   230}
    ,{   310,   310,   310,   210,   310}
    }
   ,{{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   270,   310}
    }
   ,{{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   210,   180,   180,   210,   180}
    ,{   310,   310,   310,   210,   310}
    }
   ,{{   370,   370,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   310,   310,   310,   270,   310}
    ,{   310,   310,   310,   210,   310}
    ,{   370,   370,   310,   210,   310}
    }
   }
  ,{{{   370,   300,   370,   300,   340}
    ,{   370,   300,   370,   300,   340}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   340,   300,   340,   300,   310}
    ,{   340,   300,   340,   300,   310}
    ,{   310,   260,   310,   260,   280}
    ,{   290,   180,   290,   180,   260}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   370,   260,   370,   260,   340}
    ,{   370,   260,   370,   260,   340}
    ,{   310,   260,   310,   260,   280}
    ,{   180,   130,   180,   130,   150}
    ,{   310,   260,   310,   260,   280}
    }
   ,{{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    ,{   310,   260,   310,   260,   280}
    }
   }
  ,{{{   340,   260,   340,   280,   340}
    ,{   340,   240,   340,   280,   340}
    ,{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   280,   310}
    ,{   310,   260,   310,   280,   310}
    }
   ,{{   340,   240,   340,   200,   340}
    ,{   340,   240,   340,   190,   340}
    ,{   310,   200,   310,   150,   310}
    ,{   230,   120,   230,   200,   230}
    ,{   310,   200,   310,   150,   310}
    }
   ,{{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   260,   310,   150,   310}
    }
   ,{{   310,   200,   310,   280,   310}
    ,{   310,   200,   310,   280,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   280,   200,   180,   280,   180}
    ,{   310,   200,   310,   150,   310}
    }
   ,{{   310,   260,   310,   280,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   260,   310,   150,   310}
    ,{   310,   200,   310,   150,   310}
    ,{   310,   200,   310,   280,   310}
    }
   }
  ,{{{   370,   300,   370,   300,   320}
    ,{   370,   300,   370,   300,   320}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   340,   300,   340,   300,   320}
    ,{   340,   300,   340,   300,   320}
    ,{   310,   260,   310,   260,   220}
    ,{   290,   180,   290,   180,   140}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   370,   260,   370,   260,   220}
    ,{   370,   260,   370,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   220,   130,   180,   130,   220}
    ,{   310,   260,   310,   260,   220}
    }
   ,{{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    ,{   310,   260,   310,   260,   220}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   310,   300,   270,   310,   290}
    ,{   300,   300,   270,   270,   290}
    ,{   310,   290,   250,   310,   250}
    ,{   300,   300,   270,   270,   270}
    ,{   300,   270,   240,   300,   240}
    }
   ,{{   290,   270,   230,   230,   290}
    ,{   290,   270,   230,   230,   290}
    ,{   260,   260,   220,   220,   220}
    ,{   190,   170,   190,   130,   190}
    ,{   260,   260,   220,   220,   220}
    }
   ,{{   310,   300,   270,   310,   270}
    ,{   300,   300,   270,   270,   270}
    ,{   310,   290,   250,   310,   250}
    ,{   300,   300,   270,   270,   270}
    ,{   300,   270,   240,   300,   240}
    }
   ,{{   260,   260,   220,   220,   220}
    ,{   190,   170,   190,   130,   190}
    ,{   260,   260,   220,   220,   220}
    ,{   210,   130,    80,   210,   210}
    ,{   260,   260,   220,   220,   220}
    }
   ,{{   300,   300,   270,   300,   270}
    ,{   300,   300,   270,   270,   270}
    ,{   300,   270,   240,   300,   240}
    ,{   300,   300,   270,   270,   270}
    ,{   240,   240,   150,   150,   150}
    }
   }
  ,{{{   310,   300,   270,   310,   270}
    ,{   300,   300,   270,   270,   270}
    ,{   310,   290,   250,   310,   250}
    ,{   300,   300,   270,   270,   270}
    ,{   300,   270,   240,   300,   240}
    }
   ,{{   270,   270,   230,   230,   230}
    ,{   270,   270,   230,   230,   230}
    ,{   260,   260,   220,   220,   220}
    ,{   170,   170,   130,   130,   130}
    ,{   260,   260,   220,   220,   220}
    }
   ,{{   310,   300,   270,   310,   270}
    ,{   300,   300,   270,   270,   270}
    ,{   310,   290,   250,   310,   250}
    ,{   300,   300,   270,   270,   270}
    ,{   300,   270,   240,   300,   240}
    }
   ,{{   260,   260,   220,   220,   220}
    ,{   170,   170,   130,   130,   130}
    ,{   260,   260,   220,   220,   220}
    ,{   210,   110,    80,   210,    80}
    ,{   260,   260,   220,   220,   220}
    }
   ,{{   300,   300,   270,   300,   270}
    ,{   300,   300,   270,   270,   270}
    ,{   300,   270,   240,   300,   240}
    ,{   300,   300,   270,   270,   270}
    ,{   240,   240,   150,   150,   150}
    }
   }
  ,{{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   230,   230,   230,   230,   230}
    ,{   230,   230,   230,   230,   230}
    ,{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   190}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   190}
    ,{   220,   220,   220,   220,   220}
    ,{    80,    80,    80,    80,    80}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    ,{   270,   270,   270,   270,   270}
    ,{   150,   150,   150,   150,   150}
    }
   }
  ,{{{   270,   230,   270,   210,   270}
    ,{   270,   190,   270,   140,   270}
    ,{   250,   230,   250,   120,   250}
    ,{   270,   190,   270,   210,   270}
    ,{   240,   220,   240,   150,   240}
    }
   ,{{   230,   150,   230,   130,   230}
    ,{   230,   150,   230,   100,   230}
    ,{   220,   140,   220,    90,   220}
    ,{   130,    50,   130,   130,   130}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   270,   230,   270,   140,   270}
    ,{   270,   190,   270,   140,   270}
    ,{   250,   230,   250,   120,   250}
    ,{   270,   190,   270,   140,   270}
    ,{   240,   220,   240,   110,   240}
    }
   ,{{   220,   140,   220,   210,   220}
    ,{   130,    50,   130,   130,   130}
    ,{   220,   140,   220,    90,   220}
    ,{   210,   130,    80,   210,    80}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   270,   220,   270,   150,   270}
    ,{   270,   190,   270,   140,   270}
    ,{   240,   220,   240,   110,   240}
    ,{   270,   190,   270,   140,   270}
    ,{   150,    70,   150,   150,   150}
    }
   }
  ,{{{   290,   270,   270,   270,   290}
    ,{   290,   270,   270,   270,   290}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   290,   230,   230,   230,   290}
    ,{   290,   230,   230,   230,   290}
    ,{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   130}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   250,   250,   250,   250,   250}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   190,   130,   190,   130,   130}
    ,{   220,   220,   220,   220,   220}
    ,{   210,    80,    80,    80,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   270,   270,   270,   270,   270}
    ,{   270,   270,   270,   270,   270}
    ,{   240,   240,   240,   240,   240}
    ,{   270,   270,   270,   270,   270}
    ,{   150,   150,   150,   150,   150}
    }
   }
  }
 ,{{{{   300,   280,   240,   280,   300}
    ,{   300,   280,   240,   240,   300}
    ,{   280,   260,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   250,   220,   280,   220}
    }
   ,{{   300,   280,   240,   240,   300}
    ,{   300,   280,   240,   240,   300}
    ,{   250,   250,   220,   220,   220}
    ,{   100,    70,   100,    40,   100}
    ,{   250,   250,   220,   220,   220}
    }
   ,{{   280,   250,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   250,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   250,   220,   280,   220}
    }
   ,{{   250,   250,   220,   220,   220}
    ,{   160,   140,   160,   100,   160}
    ,{   250,   250,   220,   220,   220}
    ,{   210,   130,    80,   210,   210}
    ,{   250,   250,   220,   220,   220}
    }
   ,{{   280,   260,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   260,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   240,   240,   140,   140,   140}
    }
   }
  ,{{{   280,   280,   240,   280,   240}
    ,{   280,   280,   240,   240,   240}
    ,{   280,   260,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   250,   220,   280,   220}
    }
   ,{{   280,   280,   240,   240,   240}
    ,{   280,   280,   240,   240,   240}
    ,{   250,   250,   220,   220,   220}
    ,{    70,    70,    40,    40,    40}
    ,{   250,   250,   220,   220,   220}
    }
   ,{{   280,   250,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   250,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   250,   220,   280,   220}
    }
   ,{{   250,   250,   220,   220,   220}
    ,{   140,   140,   100,   100,   100}
    ,{   250,   250,   220,   220,   220}
    ,{   210,   110,    80,   210,    80}
    ,{   250,   250,   220,   220,   220}
    }
   ,{{   280,   260,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   280,   260,   220,   280,   220}
    ,{   250,   250,   210,   210,   210}
    ,{   240,   240,   140,   140,   140}
    }
   }
  ,{{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   220,   220,   220,   220,   220}
    ,{   100,    40,   100,    40,   100}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   160,   100,   160,   100,   160}
    ,{   220,   220,   220,   220,   220}
    ,{    80,    80,    80,    80,    80}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   140,   140,   140,   140,   140}
    }
   }
  ,{{{   240,   200,   240,   210,   240}
    ,{   240,   160,   240,   110,   240}
    ,{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,   210,   210}
    ,{   220,   200,   220,   140,   220}
    }
   ,{{   240,   160,   240,   110,   240}
    ,{   240,   160,   240,   110,   240}
    ,{   220,   140,   220,    90,   220}
    ,{    40,   -40,    40,    40,    40}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   220,   200,   220,    90,   220}
    }
   ,{{   220,   140,   220,   210,   220}
    ,{   100,    20,   100,   100,   100}
    ,{   220,   140,   220,    90,   220}
    ,{   210,   130,    80,   210,    80}
    ,{   220,   140,   220,    90,   220}
    }
   ,{{   220,   200,   220,   140,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   220,   200,   220,    90,   220}
    ,{   210,   130,   210,    80,   210}
    ,{   140,    90,   140,   140,   140}
    }
   }
  ,{{{   300,   240,   240,   240,   300}
    ,{   300,   240,   240,   240,   300}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   300,   240,   240,   240,   300}
    ,{   300,   240,   240,   240,   300}
    ,{   220,   220,   220,   220,   220}
    ,{   100,    40,   100,    40,    50}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   160,   100,   160,   100,   140}
    ,{   220,   220,   220,   220,   220}
    ,{   210,    80,    80,    80,   210}
    ,{   220,   220,   220,   220,   220}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   220,   220,   220,   220,   220}
    ,{   210,   210,   210,   210,   210}
    ,{   140,   140,   140,   140,   140}
    }
   }
  }
 ,{{{{   430,   430,   370,   400,   430}
    ,{   430,   410,   370,   370,   430}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   400,   340}
    }
   ,{{   430,   410,   370,   370,   430}
    ,{   430,   410,   370,   370,   430}
    ,{   370,   370,   340,   340,   340}
    ,{   320,   290,   320,   260,   320}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    }
   ,{{   370,   370,   360,   340,   360}
    ,{   360,   360,   360,   300,   360}
    ,{   370,   370,   340,   340,   340}
    ,{   340,   260,   210,   340,   340}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   430,   430,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   340,   340}
    }
   }
  ,{{{   430,   430,   370,   400,   370}
    ,{   410,   410,   370,   370,   370}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   400,   340}
    }
   ,{{   410,   410,   370,   370,   370}
    ,{   410,   410,   370,   370,   370}
    ,{   370,   370,   340,   340,   340}
    ,{   290,   290,   260,   260,   260}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    }
   ,{{   370,   370,   340,   340,   340}
    ,{   360,   360,   300,   300,   300}
    ,{   370,   370,   340,   340,   340}
    ,{   340,   240,   210,   340,   210}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   430,   430,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   340,   340}
    }
   }
  ,{{{   370,   370,   370,   370,   370}
    ,{   370,   370,   370,   370,   370}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   370,   370,   370,   370,   370}
    ,{   370,   370,   370,   370,   370}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   320}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   360,   340,   360,   340,   360}
    ,{   360,   300,   360,   300,   360}
    ,{   340,   340,   340,   340,   340}
    ,{   210,   210,   210,   210,   210}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  ,{{{   370,   320,   370,   340,   370}
    ,{   370,   290,   370,   300,   370}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   340,   320,   340,   340,   340}
    }
   ,{{   370,   290,   370,   260,   370}
    ,{   370,   290,   370,   240,   370}
    ,{   340,   260,   340,   210,   340}
    ,{   260,   180,   260,   260,   260}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    }
   ,{{   340,   260,   340,   340,   340}
    ,{   300,   220,   300,   300,   300}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   210,   340,   210}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   340,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    }
   }
  ,{{{   430,   370,   370,   370,   430}
    ,{   430,   370,   370,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   430,   370,   370,   370,   430}
    ,{   430,   370,   370,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   260}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   360,   340,   360,   340,   340}
    ,{   360,   300,   360,   300,   300}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   210,   210,   210,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  }
 ,{{{{   400,   400,   400,   370,   400}
    ,{   400,   370,   400,   360,   400}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   400,   400,   310,   370,   310}
    }
   ,{{   360,   360,   310,   360,   330}
    ,{   360,   360,   270,   360,   330}
    ,{   340,   340,   310,   310,   310}
    ,{   230,   220,   230,   170,   230}
    ,{   340,   340,   310,   310,   310}
    }
   ,{{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    }
   ,{{   400,   370,   400,   340,   400}
    ,{   400,   370,   400,   340,   400}
    ,{   340,   340,   310,   310,   310}
    ,{   310,   230,   180,   310,   310}
    ,{   340,   340,   310,   310,   310}
    }
   ,{{   400,   400,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   400,   400,   310,   310,   310}
    }
   }
  ,{{{   400,   400,   340,   370,   340}
    ,{   370,   370,   340,   360,   340}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   400,   400,   310,   370,   310}
    }
   ,{{   360,   360,   310,   360,   310}
    ,{   360,   360,   270,   360,   270}
    ,{   340,   340,   310,   310,   310}
    ,{   220,   220,   170,   170,   170}
    ,{   340,   340,   310,   310,   310}
    }
   ,{{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    }
   ,{{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   340,   340,   310,   310,   310}
    ,{   310,   210,   180,   310,   180}
    ,{   340,   340,   310,   310,   310}
    }
   ,{{   400,   400,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   400,   400,   310,   310,   310}
    }
   }
  ,{{{   400,   340,   400,   340,   400}
    ,{   400,   340,   400,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   270,   270,   270,   270,   270}
    ,{   310,   310,   310,   310,   310}
    ,{   230,   170,   230,   170,   230}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   400,   340,   400,   340,   400}
    ,{   400,   340,   400,   340,   400}
    ,{   310,   310,   310,   310,   310}
    ,{   180,   180,   180,   180,   180}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  ,{{{   340,   290,   340,   340,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   310,   310}
    ,{   310,   290,   310,   310,   310}
    }
   ,{{   310,   230,   310,   180,   310}
    ,{   270,   190,   270,   140,   270}
    ,{   310,   230,   310,   180,   310}
    ,{   170,    40,   170,   170,   170}
    ,{   310,   230,   310,   180,   310}
    }
   ,{{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    }
   ,{{   340,   260,   340,   340,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   230,   180,   310,   180}
    ,{   310,   230,   310,   180,   310}
    }
   ,{{   310,   290,   310,   310,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   230,   310,   310,   310}
    }
   }
  ,{{{   400,   340,   400,   340,   340}
    ,{   400,   340,   400,   340,   340}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   330,   310,   310,   310,   330}
    ,{   330,   270,   270,   270,   330}
    ,{   310,   310,   310,   310,   310}
    ,{   230,   170,   230,   170,   170}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   400,   340,   400,   340,   340}
    ,{   400,   340,   400,   340,   340}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   180,   180,   180,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   }
  }
 ,{{{{   370,   340,   310,   350,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   350,   320,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    }
   ,{{   370,   340,   310,   310,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   320,   320,   280,   280,   280}
    ,{   240,   220,   240,   180,   240}
    ,{   320,   320,   280,   280,   280}
    }
   ,{{   350,   330,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    }
   ,{{   320,   320,   310,   280,   310}
    ,{   310,   290,   310,   250,   310}
    ,{   320,   320,   280,   280,   280}
    ,{   260,   180,   130,   260,   260}
    ,{   320,   320,   280,   280,   280}
    }
   ,{{   350,   330,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   290,   290,   200,   200,   200}
    }
   }
  ,{{{   350,   340,   310,   350,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   350,   320,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    }
   ,{{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   320,   320,   280,   280,   280}
    ,{   220,   220,   180,   180,   180}
    ,{   320,   320,   280,   280,   280}
    }
   ,{{   350,   330,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    }
   ,{{   320,   320,   280,   280,   280}
    ,{   290,   290,   250,   250,   250}
    ,{   320,   320,   280,   280,   280}
    ,{   260,   170,   130,   260,   130}
    ,{   320,   320,   280,   280,   280}
    }
   ,{{   350,   330,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   350,   320,   290,   350,   290}
    ,{   330,   330,   290,   290,   290}
    ,{   290,   290,   200,   200,   200}
    }
   }
  ,{{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   240,   180,   240,   180,   240}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   310,   280,   310,   280,   310}
    ,{   310,   250,   310,   250,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   130,   130,   130,   130,   130}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   200,   200,   200,   200,   200}
    }
   }
  ,{{{   310,   270,   310,   260,   310}
    ,{   310,   230,   310,   250,   310}
    ,{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   260,   290}
    ,{   290,   270,   290,   200,   290}
    }
   ,{{   310,   230,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   280,   200,   280,   150,   280}
    ,{   180,   100,   180,   180,   180}
    ,{   280,   200,   280,   150,   280}
    }
   ,{{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   290,   270,   290,   160,   290}
    }
   ,{{   280,   200,   280,   260,   280}
    ,{   250,   170,   250,   250,   250}
    ,{   280,   200,   280,   150,   280}
    ,{   260,   180,   130,   260,   130}
    ,{   280,   200,   280,   150,   280}
    }
   ,{{   290,   270,   290,   200,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   290,   270,   290,   160,   290}
    ,{   290,   210,   290,   160,   290}
    ,{   200,   120,   200,   200,   200}
    }
   }
  ,{{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   280,   280,   280,   280,   280}
    ,{   240,   180,   240,   180,   180}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    }
   ,{{   310,   280,   310,   280,   280}
    ,{   310,   250,   310,   250,   250}
    ,{   280,   280,   280,   280,   280}
    ,{   260,   130,   130,   130,   260}
    ,{   280,   280,   280,   280,   280}
    }
   ,{{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   290,   290,   290,   290,   290}
    ,{   200,   200,   200,   200,   200}
    }
   }
  }
 ,{{{{   370,   340,   310,   370,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    }
   ,{{   370,   340,   310,   310,   370}
    ,{   370,   340,   310,   310,   370}
    ,{   300,   300,   260,   260,   260}
    ,{   260,   240,   260,   200,   260}
    ,{   300,   300,   260,   260,   260}
    }
   ,{{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    }
   ,{{   300,   300,   270,   280,   280}
    ,{   270,   250,   270,   210,   270}
    ,{   300,   300,   260,   260,   260}
    ,{   280,   200,   150,   280,   280}
    ,{   300,   300,   260,   260,   260}
    }
   ,{{   340,   340,   310,   340,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   310,   280,   340,   280}
    ,{   340,   340,   310,   310,   310}
    ,{   320,   320,   220,   220,   220}
    }
   }
  ,{{{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    }
   ,{{   340,   340,   310,   310,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   300,   300,   260,   260,   260}
    ,{   240,   240,   200,   200,   200}
    ,{   300,   300,   260,   260,   260}
    }
   ,{{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   370,   340,   310,   370,   310}
    }
   ,{{   300,   300,   260,   280,   260}
    ,{   250,   250,   210,   210,   210}
    ,{   300,   300,   260,   260,   260}
    ,{   280,   190,   150,   280,   150}
    ,{   300,   300,   260,   260,   260}
    }
   ,{{   340,   340,   310,   340,   310}
    ,{   340,   340,   310,   310,   310}
    ,{   340,   310,   280,   340,   280}
    ,{   340,   340,   310,   310,   310}
    ,{   320,   320,   220,   220,   220}
    }
   }
  ,{{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   200,   260,   200,   260}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   270,   260,   270,   260,   270}
    ,{   270,   210,   270,   210,   270}
    ,{   260,   260,   260,   260,   260}
    ,{   150,   150,   150,   150,   150}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   310,   310,   310,   310,   310}
    ,{   220,   220,   220,   220,   220}
    }
   }
  ,{{{   310,   290,   310,   280,   310}
    ,{   310,   230,   310,   210,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   280,   310}
    ,{   310,   290,   310,   220,   310}
    }
   ,{{   310,   230,   310,   200,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   260,   180,   260,   130,   260}
    ,{   200,   120,   200,   200,   200}
    ,{   260,   180,   260,   130,   260}
    }
   ,{{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   310,   290,   310,   180,   310}
    }
   ,{{   280,   200,   260,   280,   260}
    ,{   210,   130,   210,   210,   210}
    ,{   260,   180,   260,   130,   260}
    ,{   280,   200,   150,   280,   150}
    ,{   260,   180,   260,   130,   260}
    }
   ,{{   310,   260,   310,   220,   310}
    ,{   310,   230,   310,   180,   310}
    ,{   280,   260,   280,   150,   280}
    ,{   310,   230,   310,   180,   310}
    ,{   220,   140,   220,   220,   220}
    }
   }
  ,{{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   370,   310,   310,   310,   370}
    ,{   370,   310,   310,   310,   370}
    ,{   260,   260,   260,   260,   260}
    ,{   260,   200,   260,   200,   200}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    }
   ,{{   280,   260,   270,   260,   280}
    ,{   270,   210,   270,   210,   210}
    ,{   260,   260,   260,   260,   260}
    ,{   280,   150,   150,   150,   280}
    ,{   260,   260,   260,   260,   260}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{   280,   280,   280,   280,   280}
    ,{   310,   310,   310,   310,   310}
    ,{   220,   220,   220,   220,   220}
    }
   }
  }
 ,{{{{   430,   430,   400,   400,   430}
    ,{   430,   410,   400,   370,   430}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   400,   340}
    }
   ,{{   430,   410,   370,   370,   430}
    ,{   430,   410,   370,   370,   430}
    ,{   370,   370,   340,   340,   340}
    ,{   320,   290,   320,   260,   320}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    }
   ,{{   400,   370,   400,   340,   400}
    ,{   400,   370,   400,   340,   400}
    ,{   370,   370,   340,   340,   340}
    ,{   340,   260,   210,   340,   340}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   430,   430,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   340,   340}
    }
   }
  ,{{{   430,   430,   370,   400,   370}
    ,{   410,   410,   370,   370,   370}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   400,   340}
    }
   ,{{   410,   410,   370,   370,   370}
    ,{   410,   410,   370,   370,   370}
    ,{   370,   370,   340,   340,   340}
    ,{   290,   290,   260,   260,   260}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    }
   ,{{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   340,   240,   210,   340,   210}
    ,{   370,   370,   340,   340,   340}
    }
   ,{{   430,   430,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   400,   370,   340,   400,   340}
    ,{   370,   370,   340,   340,   340}
    ,{   430,   430,   340,   340,   340}
    }
   }
  ,{{{   400,   370,   400,   370,   400}
    ,{   400,   370,   400,   370,   400}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   370,   370,   370,   370,   370}
    ,{   370,   370,   370,   370,   370}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   320}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   400,   340,   400,   340,   400}
    ,{   400,   340,   400,   340,   400}
    ,{   340,   340,   340,   340,   340}
    ,{   210,   210,   210,   210,   210}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  ,{{{   370,   320,   370,   340,   370}
    ,{   370,   290,   370,   340,   370}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   340,   320,   340,   340,   340}
    }
   ,{{   370,   290,   370,   260,   370}
    ,{   370,   290,   370,   240,   370}
    ,{   340,   260,   340,   210,   340}
    ,{   260,   180,   260,   260,   260}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    }
   ,{{   340,   260,   340,   340,   340}
    ,{   340,   260,   340,   340,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   210,   340,   210}
    ,{   340,   260,   340,   210,   340}
    }
   ,{{   340,   320,   340,   340,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   320,   340,   210,   340}
    ,{   340,   260,   340,   210,   340}
    ,{   340,   260,   340,   340,   340}
    }
   }
  ,{{{   430,   370,   400,   370,   430}
    ,{   430,   370,   400,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   430,   370,   370,   370,   430}
    ,{   430,   370,   370,   370,   430}
    ,{   340,   340,   340,   340,   340}
    ,{   320,   260,   320,   260,   260}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   400,   340,   400,   340,   340}
    ,{   400,   340,   400,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   210,   210,   210,   340}
    ,{   340,   340,   340,   340,   340}
    }
   ,{{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    ,{   340,   340,   340,   340,   340}
    }
   }
  }
 }};
PUBLIC int int22_dH[NBPAIRS+1][NBPAIRS+1][5][5][5][5] =
{{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{    80,  -120,    30,    80,    80}
    ,{    30,  -310,  -170,    30,  -110}
    ,{    80,  -230,  -110,    80,   -60}
    ,{    80,  -120,    30,    30,    80}
    ,{   -30,  -340,  -220,   -30,  -170}
    }
   ,{{  -120,  -460,  -290,  -120,  -230}
    ,{  -120,  -460,  -310,  -120,  -260}
    ,{  -430,  -770,  -620,  -430,  -570}
    ,{  -230,  -670,  -290,  -980,  -230}
    ,{  -430,  -770,  -620,  -430,  -570}
    }
   ,{{    30,  -290,  -170,    30,  -110}
    ,{    30,  -310,  -170,    30,  -110}
    ,{    20,  -290,  -170,    20,  -120}
    ,{    30,  -310,  -170,    30,  -110}
    ,{   -30,  -340,  -220,   -30,  -170}
    }
   ,{{    80,  -120,    30,  -430,    80}
    ,{  -520,  -960,  -580, -1270,  -520}
    ,{  -430,  -770,  -620,  -430,  -570}
    ,{    80,  -120,    30,  -430,    80}
    ,{  -430,  -770,  -620,  -430,  -570}
    }
   ,{{    80,  -230,  -110,    80,   -60}
    ,{    30,  -310,  -170,    30,  -110}
    ,{    80,  -230,  -110,    80,   -60}
    ,{    30,  -310,  -170,    30,  -110}
    ,{  -860,  -860,  -960, -1410,  -900}
    }
   }
  ,{{{    30,  -120,    30,  -520,    30}
    ,{  -170,  -310,  -170,  -810,  -170}
    ,{  -110,  -260,  -110,  -520,  -110}
    ,{    30,  -120,    30,  -810,    30}
    ,{  -220,  -370,  -220,  -630,  -220}
    }
   ,{{  -310,  -460,  -310,  -960,  -310}
    ,{  -310,  -460,  -310,  -960,  -310}
    ,{  -620,  -770,  -620, -1270,  -620}
    ,{  -530,  -670,  -530, -1170,  -530}
    ,{  -620,  -770,  -620, -1270,  -620}
    }
   ,{{  -170,  -310,  -170,  -580,  -170}
    ,{  -170,  -310,  -170,  -810,  -170}
    ,{  -170,  -320,  -170,  -580,  -170}
    ,{  -170,  -310,  -170,  -810,  -170}
    ,{  -220,  -370,  -220,  -630,  -220}
    }
   ,{{    30,  -120,    30, -1270,    30}
    ,{  -810,  -960,  -810, -1460,  -810}
    ,{  -620,  -770,  -620, -1270,  -620}
    ,{    30,  -120,    30, -1870,    30}
    ,{  -620,  -770,  -620, -1270,  -620}
    }
   ,{{  -110,  -260,  -110,  -520,  -110}
    ,{  -170,  -310,  -170,  -810,  -170}
    ,{  -110,  -260,  -110,  -520,  -110}
    ,{  -170,  -310,  -170,  -810,  -170}
    ,{  -860,  -860,  -960, -1600,  -960}
    }
   }
  ,{{{    80,  -430,    20,  -430,    80}
    ,{  -110,  -620,  -170,  -620,  -110}
    ,{   -60,  -570,  -120,  -570,   -60}
    ,{    80,  -430,    20,  -430,    80}
    ,{  -170,  -680,  -230,  -680,  -170}
    }
   ,{{  -230,  -770,  -290,  -770,  -230}
    ,{  -260,  -770,  -320,  -770,  -260}
    ,{  -570, -1080,  -630, -1080,  -570}
    ,{  -230,  -980,  -290,  -980,  -230}
    ,{  -570, -1080,  -630, -1080,  -570}
    }
   ,{{  -110,  -620,  -170,  -620,  -110}
    ,{  -110,  -620,  -170,  -620,  -110}
    ,{  -120,  -630,  -180,  -630,  -120}
    ,{  -110,  -620,  -170,  -620,  -110}
    ,{  -170,  -680,  -230,  -680,  -170}
    }
   ,{{    80,  -430,    20,  -430,    80}
    ,{  -520, -1270,  -580, -1270,  -520}
    ,{  -570, -1080,  -630, -1080,  -570}
    ,{    80,  -430,    20,  -430,    80}
    ,{  -570, -1080,  -630, -1080,  -570}
    }
   ,{{   -60,  -570,  -120,  -570,   -60}
    ,{  -110,  -620,  -170,  -620,  -110}
    ,{   -60,  -570,  -120,  -570,   -60}
    ,{  -110,  -620,  -170,  -620,  -110}
    ,{  -900, -1410,  -960, -1410,  -900}
    }
   }
  ,{{{    80,  -230,    30,    80,    30}
    ,{    30,  -530,  -170,    30,  -170}
    ,{    80,  -230,  -110,    80,  -110}
    ,{    30,  -530,    30,    30,    30}
    ,{   -30,  -340,  -220,   -30,  -220}
    }
   ,{{  -120,  -670,  -310,  -120,  -310}
    ,{  -120,  -670,  -310,  -120,  -310}
    ,{  -430,  -980,  -620,  -430,  -620}
    ,{  -530,  -890,  -530, -1580,  -530}
    ,{  -430,  -980,  -620,  -430,  -620}
    }
   ,{{    30,  -290,  -170,    30,  -170}
    ,{    30,  -530,  -170,    30,  -170}
    ,{    20,  -290,  -170,    20,  -170}
    ,{    30,  -530,  -170,    30,  -170}
    ,{   -30,  -340,  -220,   -30,  -220}
    }
   ,{{    30,  -980,    30,  -430,    30}
    ,{  -810, -1170,  -810, -1870,  -810}
    ,{  -430,  -980,  -620,  -430,  -620}
    ,{    30, -1580,    30, -2280,    30}
    ,{  -430,  -980,  -620,  -430,  -620}
    }
   ,{{    80,  -230,  -110,    80,  -110}
    ,{    30,  -530,  -170,    30,  -170}
    ,{    80,  -230,  -110,    80,  -110}
    ,{    30,  -530,  -170,    30,  -170}
    ,{  -960, -1320,  -960, -2010,  -960}
    }
   }
  ,{{{   -30,  -430,   -30,  -430,  -860}
    ,{  -220,  -620,  -220,  -620,  -860}
    ,{  -170,  -570,  -170,  -570,  -900}
    ,{   -30,  -430,   -30,  -430,  -960}
    ,{  -280,  -680,  -280,  -680, -1010}
    }
   ,{{  -340,  -770,  -340,  -770,  -860}
    ,{  -370,  -770,  -370,  -770,  -860}
    ,{  -680, -1080,  -680, -1080, -1410}
    ,{  -340,  -980,  -340,  -980, -1320}
    ,{  -680, -1080,  -680, -1080, -1410}
    }
   ,{{  -220,  -620,  -220,  -620,  -960}
    ,{  -220,  -620,  -220,  -620,  -960}
    ,{  -230,  -630,  -230,  -630,  -960}
    ,{  -220,  -620,  -220,  -620,  -960}
    ,{  -280,  -680,  -280,  -680, -1010}
    }
   ,{{   -30,  -430,   -30,  -430, -1410}
    ,{  -630, -1270,  -630, -1270, -1600}
    ,{  -680, -1080,  -680, -1080, -1410}
    ,{   -30,  -430,   -30,  -430, -2010}
    ,{  -680, -1080,  -680, -1080, -1410}
    }
   ,{{  -170,  -570,  -170,  -570,  -900}
    ,{  -220,  -620,  -220,  -620,  -960}
    ,{  -170,  -570,  -170,  -570,  -900}
    ,{  -220,  -620,  -220,  -620,  -960}
    ,{ -1010, -1410, -1010, -1410, -1750}
    }
   }
  }
 ,{{{{   540,   180,    30,   540,   180}
    ,{    10,  -580,  -150,    10,   -90}
    ,{   540,  -350,  -600,   540,  -540}
    ,{   180,   180,    30,  -320,   180}
    ,{   -90,  -740,   -90,  -260,  -540}
    }
   ,{{   -90,  -350,  -150,  -100,   -90}
    ,{   -90,  -580,  -150,  -200,   -90}
    ,{  -100,  -350,  -600,  -100,  -540}
    ,{  -630, -1790,  -630, -1790, -1040}
    ,{  -400,  -740,  -600,  -400,  -540}
    }
   ,{{   540,  -660,  -510,   540,  -400}
    ,{    10,  -660,  -510,    10,  -400}
    ,{   540,  -940,  -820,   540,  -760}
    ,{  -320,  -660,  -510,  -320,  -460}
    ,{  -260,  -940,  -820,  -260,  -550}
    }
   ,{{   180,   180,    30,  -400,   180}
    ,{  -500, -1070,  -500, -1080,  -570}
    ,{  -400,  -740,  -600,  -400,  -540}
    ,{   180,   180,    30,  -430,   180}
    ,{  -400,  -740,  -600,  -400,  -540}
    }
   ,{{   -90,  -660,   -90,  -210,  -460}
    ,{  -320,  -660,  -510,  -320,  -460}
    ,{  -210, -1250, -1130,  -210, -1070}
    ,{  -320,  -660,  -510,  -320,  -460}
    ,{   -90,  -830,   -90,  -810,  -800}
    }
   }
  ,{{{   540,   180,   -90,   540,    30}
    ,{    10,  -580,  -220,    10,  -150}
    ,{   540,  -740,  -600,   540,  -600}
    ,{   180,   180,  -390, -1160,    30}
    ,{   -90,  -740,   -90,  -810,  -600}
    }
   ,{{  -100,  -580,  -220,  -100,  -150}
    ,{  -150,  -580,  -220,  -970,  -150}
    ,{  -100,  -740,  -600,  -100,  -600}
    ,{ -1340, -2010, -1650, -1980, -1340}
    ,{  -600,  -740,  -600, -1240,  -600}
    }
   ,{{   540,  -660,  -510,   540,  -510}
    ,{    10,  -660, -1150,    10,  -510}
    ,{   540,  -960,  -820,   540,  -820}
    ,{  -510,  -660,  -510, -1160,  -510}
    ,{  -820,  -960,  -820, -1220,  -820}
    }
   ,{{   180,   180,  -390, -1240,    30}
    ,{  -860, -1340,  -860, -2450,  -860}
    ,{  -600,  -740,  -600, -1240,  -600}
    ,{   180,   180,  -390, -1870,    30}
    ,{  -600,  -740,  -600, -1240,  -600}
    }
   ,{{   -90,  -660,   -90,  -810,  -510}
    ,{  -510,  -660,  -510, -1160,  -510}
    ,{ -1130, -1270, -1130, -1530, -1130}
    ,{  -510,  -660,  -510, -1160,  -510}
    ,{   -90, -1240,   -90,  -810,  -800}
    }
   }
  ,{{{   180,  -430,    20,  -430,   180}
    ,{   -90,  -600,  -500,  -600,   -90}
    ,{  -540, -1050,  -600, -1050,  -540}
    ,{   180,  -430,    20,  -430,   180}
    ,{  -540,  -830,  -600, -1050,  -540}
    }
   ,{{   -90,  -600,  -600,  -600,   -90}
    ,{   -90,  -600, -1070,  -600,   -90}
    ,{  -540, -1050,  -600, -1050,  -540}
    ,{  -630, -1790,  -630, -1790, -1040}
    ,{  -540, -1050,  -600, -1050,  -540}
    }
   ,{{  -460,  -970,  -520,  -970,  -460}
    ,{  -460,  -970,  -750,  -970,  -460}
    ,{  -760, -1270,  -820, -1270,  -760}
    ,{  -460,  -970,  -520,  -970,  -460}
    ,{  -550, -1270,  -820, -1270,  -550}
    }
   ,{{   180,  -430,    20,  -430,   180}
    ,{  -500, -1070,  -500, -1320,  -570}
    ,{  -540, -1050,  -600, -1050,  -540}
    ,{   180,  -430,    20,  -430,   180}
    ,{  -540, -1050,  -600, -1050,  -540}
    }
   ,{{  -460,  -830,  -520,  -970,  -460}
    ,{  -460,  -970,  -520,  -970,  -460}
    ,{ -1070, -1580, -1130, -1580, -1070}
    ,{  -460,  -970,  -520,  -970,  -460}
    ,{  -830,  -830, -1710, -1260, -1460}
    }
   }
  ,{{{    30,  -350,    30,  -200,    30}
    ,{  -150,  -870,  -150,  -200,  -150}
    ,{  -210,  -350,  -600,  -210,  -600}
    ,{    30,  -870,    30,  -320,    30}
    ,{  -260,  -940,  -600,  -260,  -600}
    }
   ,{{  -150,  -350,  -150,  -200,  -150}
    ,{  -150, -1600,  -150,  -200,  -150}
    ,{  -350,  -350,  -600,  -440,  -600}
    ,{ -1340, -3070, -1340, -2390, -1340}
    ,{  -400,  -960,  -600,  -400,  -600}
    }
   ,{{  -260,  -870,  -510,  -260,  -510}
    ,{  -320, -1110,  -510,  -320,  -510}
    ,{  -620,  -940,  -820,  -620,  -820}
    ,{  -320,  -870,  -510,  -320,  -510}
    ,{  -260,  -940,  -820,  -260,  -820}
    }
   ,{{    30,  -960,    30,  -400,    30}
    ,{  -860, -1880,  -860, -1080,  -860}
    ,{  -400,  -960,  -600,  -400,  -600}
    ,{    30, -1370,    30, -2280,    30}
    ,{  -400,  -960,  -600,  -400,  -600}
    }
   ,{{  -210,  -870,  -510,  -210,  -510}
    ,{  -320,  -870,  -510,  -320,  -510}
    ,{  -210, -1250, -1130,  -210, -1130}
    ,{  -320,  -870,  -510,  -320,  -510}
    ,{  -800, -1360,  -800, -1550,  -800}
    }
   }
  ,{{{  -200,  -430,  -200,  -430,  -230}
    ,{  -200,  -600,  -200,  -600,  -400}
    ,{  -650, -1050,  -650, -1050, -1390}
    ,{  -230,  -430,  -570,  -430,  -230}
    ,{  -650, -1050,  -650, -1050, -1390}
    }
   ,{{  -200,  -600,  -200,  -600, -1390}
    ,{  -200,  -600,  -200,  -600, -1490}
    ,{  -650, -1050,  -650, -1050, -1390}
    ,{ -1150, -1790, -1150, -1790, -1520}
    ,{  -650, -1050,  -650, -1050, -1390}
    }
   ,{{  -400,  -970,  -570,  -970,  -400}
    ,{  -400,  -970,  -570,  -970,  -400}
    ,{  -870, -1270,  -870, -1270, -1610}
    ,{  -570,  -970,  -570,  -970, -1300}
    ,{  -870, -1270,  -870, -1270, -1610}
    }
   ,{{  -230,  -430,  -650,  -430,  -230}
    ,{ -1300, -1320, -1750, -1320, -1300}
    ,{  -650, -1050,  -650, -1050, -1390}
    ,{  -230,  -430,  -880,  -430,  -230}
    ,{  -650, -1050,  -650, -1050, -1390}
    }
   ,{{  -570,  -970,  -570,  -970, -1300}
    ,{  -570,  -970,  -570,  -970, -1300}
    ,{ -1180, -1580, -1180, -1580, -1920}
    ,{  -570,  -970,  -570,  -970, -1300}
    ,{  -860, -1260,  -860, -1260, -2350}
    }
   }
  }
 ,{{{{   240,    40,   190,  -270,   240}
    ,{  -590, -1030,  -650,  -870,  -590}
    ,{  -870, -1180, -1060,  -870, -1010}
    ,{   240,    40,   190,  -270,   240}
    ,{  -870,  -970, -1060,  -870, -1010}
    }
   ,{{  -780, -1210,  -840,  -870,  -780}
    ,{ -1050, -1370, -1240, -1050, -1190}
    ,{  -870, -1210, -1060,  -870, -1010}
    ,{  -780, -1220,  -840, -1530,  -780}
    ,{  -870, -1210, -1060,  -870, -1010}
    }
   ,{{  -870, -1180, -1060,  -870, -1010}
    ,{  -870, -1210, -1060,  -870, -1010}
    ,{  -870, -1180, -1060,  -870, -1010}
    ,{  -870, -1210, -1060,  -870, -1010}
    ,{  -870, -1180, -1060,  -870, -1010}
    }
   ,{{   240,    40,   190,  -270,   240}
    ,{  -590, -1030,  -650, -1340,  -590}
    ,{  -870, -1210, -1060,  -870, -1010}
    ,{   240,    40,   190,  -270,   240}
    ,{  -870, -1210, -1060,  -870, -1010}
    }
   ,{{  -870,  -970, -1060,  -870, -1010}
    ,{  -870, -1210, -1060,  -870, -1010}
    ,{  -870, -1180, -1060,  -870, -1010}
    ,{  -870, -1210, -1060,  -870, -1010}
    ,{  -970,  -970, -1060, -1520, -1010}
    }
   }
  ,{{{   190,    40,   190, -1470,   190}
    ,{  -890, -1030,  -890, -1530,  -890}
    ,{ -1060, -1210, -1060, -1470, -1060}
    ,{   190,    40,   190, -1710,   190}
    ,{  -970,  -970, -1060, -1470, -1060}
    }
   ,{{ -1060, -1210, -1060, -1710, -1060}
    ,{ -1240, -1370, -1240, -1890, -1240}
    ,{ -1060, -1210, -1060, -1710, -1060}
    ,{ -1080, -1220, -1080, -1720, -1080}
    ,{ -1060, -1210, -1060, -1710, -1060}
    }
   ,{{ -1060, -1210, -1060, -1470, -1060}
    ,{ -1060, -1210, -1060, -1710, -1060}
    ,{ -1060, -1210, -1060, -1470, -1060}
    ,{ -1060, -1210, -1060, -1710, -1060}
    ,{ -1060, -1210, -1060, -1470, -1060}
    }
   ,{{   190,    40,   190, -1530,   190}
    ,{  -890, -1030,  -890, -1530,  -890}
    ,{ -1060, -1210, -1060, -1710, -1060}
    ,{   190,    40,   190, -1710,   190}
    ,{ -1060, -1210, -1060, -1710, -1060}
    }
   ,{{  -970,  -970, -1060, -1470, -1060}
    ,{ -1060, -1210, -1060, -1710, -1060}
    ,{ -1060, -1210, -1060, -1470, -1060}
    ,{ -1060, -1210, -1060, -1710, -1060}
    ,{  -970,  -970, -1060, -1710, -1060}
    }
   }
  ,{{{   240,  -270,   180,  -270,   240}
    ,{  -590, -1340,  -650, -1340,  -590}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{   240,  -270,   180,  -270,   240}
    ,{ -1010, -1520, -1070, -1520, -1010}
    }
   ,{{  -780, -1520,  -840, -1520,  -780}
    ,{ -1190, -1700, -1250, -1700, -1190}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{  -780, -1530,  -840, -1530,  -780}
    ,{ -1010, -1520, -1070, -1520, -1010}
    }
   ,{{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    }
   ,{{   240,  -270,   180,  -270,   240}
    ,{  -590, -1340,  -650, -1340,  -590}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{   240,  -270,   180,  -270,   240}
    ,{ -1010, -1520, -1070, -1520, -1010}
    }
   ,{{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    ,{ -1010, -1520, -1070, -1520, -1010}
    }
   }
  ,{{{   190, -1180,   190,  -870,   190}
    ,{  -870, -1250,  -890,  -870,  -890}
    ,{  -870, -1180, -1060,  -870, -1060}
    ,{   190, -1420,   190,  -870,   190}
    ,{  -870, -1180, -1060,  -870, -1060}
    }
   ,{{  -870, -1420, -1060,  -870, -1060}
    ,{ -1050, -1600, -1240, -1050, -1240}
    ,{  -870, -1420, -1060,  -870, -1060}
    ,{ -1080, -1440, -1080, -2130, -1080}
    ,{  -870, -1420, -1060,  -870, -1060}
    }
   ,{{  -870, -1180, -1060,  -870, -1060}
    ,{  -870, -1420, -1060,  -870, -1060}
    ,{  -870, -1180, -1060,  -870, -1060}
    ,{  -870, -1420, -1060,  -870, -1060}
    ,{  -870, -1180, -1060,  -870, -1060}
    }
   ,{{   190, -1250,   190,  -870,   190}
    ,{  -890, -1250,  -890, -1940,  -890}
    ,{  -870, -1420, -1060,  -870, -1060}
    ,{   190, -1420,   190, -2120,   190}
    ,{  -870, -1420, -1060,  -870, -1060}
    }
   ,{{  -870, -1180, -1060,  -870, -1060}
    ,{  -870, -1420, -1060,  -870, -1060}
    ,{  -870, -1180, -1060,  -870, -1060}
    ,{  -870, -1420, -1060,  -870, -1060}
    ,{ -1060, -1420, -1060, -2120, -1060}
    }
   }
  ,{{{   130,  -270,   130,  -270, -1680}
    ,{  -700, -1340,  -700, -1340, -1680}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{   130,  -270,   130,  -270, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    }
   ,{{  -890, -1520,  -890, -1520, -1790}
    ,{ -1300, -1700, -1300, -1700, -1790}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{  -890, -1530,  -890, -1530, -1870}
    ,{ -1120, -1520, -1120, -1520, -1850}
    }
   ,{{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    }
   ,{{   130,  -270,   130,  -270, -1680}
    ,{  -700, -1340,  -700, -1340, -1680}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{   130,  -270,   130,  -270, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    }
   ,{{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    ,{ -1120, -1520, -1120, -1520, -1850}
    }
   }
  }
 ,{{{{   800,   600,   740,   290,   800}
    ,{   200,  -140,     0,   200,    50}
    ,{  -310,  -630,  -510,  -310,  -450}
    ,{   800,   600,   740,   290,   800}
    ,{  -310,  -410,  -510,  -310,  -450}
    }
   ,{{   200,  -140,     0,   200,    50}
    ,{   200,  -140,     0,   200,    50}
    ,{  -310,  -650,  -510,  -310,  -450}
    ,{  -550,  -990,  -610, -1300,  -550}
    ,{  -310,  -650,  -510,  -310,  -450}
    }
   ,{{  -310,  -630,  -510,  -310,  -450}
    ,{  -310,  -650,  -510,  -310,  -450}
    ,{  -310,  -630,  -510,  -310,  -450}
    ,{  -310,  -650,  -510,  -310,  -450}
    ,{  -310,  -630,  -510,  -310,  -450}
    }
   ,{{   800,   600,   740,   290,   800}
    ,{  -720, -1160,  -780, -1470,  -720}
    ,{  -310,  -650,  -510,  -310,  -450}
    ,{   800,   600,   740,   290,   800}
    ,{  -310,  -650,  -510,  -310,  -450}
    }
   ,{{  -310,  -410,  -510,  -310,  -450}
    ,{  -310,  -650,  -510,  -310,  -450}
    ,{  -310,  -630,  -510,  -310,  -450}
    ,{  -310,  -650,  -510,  -310,  -450}
    ,{  -410,  -410,  -510,  -960,  -450}
    }
   }
  ,{{{   740,   600,   740,  -640,   740}
    ,{     0,  -140,     0,  -640,     0}
    ,{  -510,  -650,  -510,  -910,  -510}
    ,{   740,   600,   740, -1150,   740}
    ,{  -410,  -410,  -510,  -910,  -510}
    }
   ,{{     0,  -140,     0,  -640,     0}
    ,{     0,  -140,     0,  -640,     0}
    ,{  -510,  -650,  -510, -1150,  -510}
    ,{  -850,  -990,  -850, -1490,  -850}
    ,{  -510,  -650,  -510, -1150,  -510}
    }
   ,{{  -510,  -650,  -510,  -910,  -510}
    ,{  -510,  -650,  -510, -1150,  -510}
    ,{  -510,  -650,  -510,  -910,  -510}
    ,{  -510,  -650,  -510, -1150,  -510}
    ,{  -510,  -650,  -510,  -910,  -510}
    }
   ,{{   740,   600,   740, -1150,   740}
    ,{ -1020, -1160, -1020, -1660, -1020}
    ,{  -510,  -650,  -510, -1150,  -510}
    ,{   740,   600,   740, -1150,   740}
    ,{  -510,  -650,  -510, -1150,  -510}
    }
   ,{{  -410,  -410,  -510,  -910,  -510}
    ,{  -510,  -650,  -510, -1150,  -510}
    ,{  -510,  -650,  -510,  -910,  -510}
    ,{  -510,  -650,  -510, -1150,  -510}
    ,{  -410,  -410,  -510, -1150,  -510}
    }
   }
  ,{{{   800,   290,   740,   290,   800}
    ,{    50,  -450,     0,  -450,    50}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{   800,   290,   740,   290,   800}
    ,{  -450,  -960,  -510,  -960,  -450}
    }
   ,{{    50,  -450,     0,  -450,    50}
    ,{    50,  -450,     0,  -450,    50}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{  -550, -1300,  -610, -1300,  -550}
    ,{  -450,  -960,  -510,  -960,  -450}
    }
   ,{{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    }
   ,{{   800,   290,   740,   290,   800}
    ,{  -720, -1470,  -780, -1470,  -720}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{   800,   290,   740,   290,   800}
    ,{  -450,  -960,  -510,  -960,  -450}
    }
   ,{{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    ,{  -450,  -960,  -510,  -960,  -450}
    }
   }
  ,{{{   740,  -360,   740,   200,   740}
    ,{   200,  -360,     0,   200,     0}
    ,{  -310,  -630,  -510,  -310,  -510}
    ,{   740,  -870,   740,  -310,   740}
    ,{  -310,  -630,  -510,  -310,  -510}
    }
   ,{{   200,  -360,     0,   200,     0}
    ,{   200,  -360,     0,   200,     0}
    ,{  -310,  -870,  -510,  -310,  -510}
    ,{  -850, -1210,  -850, -1900,  -850}
    ,{  -310,  -870,  -510,  -310,  -510}
    }
   ,{{  -310,  -630,  -510,  -310,  -510}
    ,{  -310,  -870,  -510,  -310,  -510}
    ,{  -310,  -630,  -510,  -310,  -510}
    ,{  -310,  -870,  -510,  -310,  -510}
    ,{  -310,  -630,  -510,  -310,  -510}
    }
   ,{{   740,  -870,   740,  -310,   740}
    ,{ -1020, -1380, -1020, -2070, -1020}
    ,{  -310,  -870,  -510,  -310,  -510}
    ,{   740,  -870,   740, -1560,   740}
    ,{  -310,  -870,  -510,  -310,  -510}
    }
   ,{{  -310,  -630,  -510,  -310,  -510}
    ,{  -310,  -870,  -510,  -310,  -510}
    ,{  -310,  -630,  -510,  -310,  -510}
    ,{  -310,  -870,  -510,  -310,  -510}
    ,{  -510,  -870,  -510, -1560,  -510}
    }
   }
  ,{{{   690,   290,   690,   290,  -550}
    ,{   -50,  -450,   -50,  -450,  -550}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{   690,   290,   690,   290, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    }
   ,{{   -50,  -450,   -50,  -450,  -550}
    ,{   -50,  -450,   -50,  -450,  -550}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{  -660, -1300,  -660, -1300, -1640}
    ,{  -560,  -960,  -560,  -960, -1300}
    }
   ,{{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    }
   ,{{   690,   290,   690,   290, -1300}
    ,{  -830, -1470,  -830, -1470, -1810}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{   690,   290,   690,   290, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    }
   ,{{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    ,{  -560,  -960,  -560,  -960, -1300}
    }
   }
  }
 ,{{{{  1170,   970,  1120,   780,  1170}
    ,{   780,   440,   580,   780,   640}
    ,{   480,   170,   280,   480,   340}
    ,{  1170,   970,  1120,   660,  1170}
    ,{   480,   170,   280,   480,   340}
    }
   ,{{   780,   440,   580,   780,   640}
    ,{   780,   440,   580,   780,   640}
    ,{   470,   130,   270,   470,   330}
    ,{  -510,  -950,  -570, -1260,  -510}
    ,{   470,   130,   270,   470,   330}
    }
   ,{{   490,   170,   290,   490,   340}
    ,{   490,   140,   290,   490,   340}
    ,{   480,   170,   280,   480,   340}
    ,{   490,   140,   290,   490,   340}
    ,{   480,   170,   280,   480,   340}
    }
   ,{{  1170,   970,  1120,   660,  1170}
    ,{  -330,  -770,  -390, -1080,  -330}
    ,{   470,   130,   270,   470,   330}
    ,{  1170,   970,  1120,   660,  1170}
    ,{   470,   130,   270,   470,   330}
    }
   ,{{   490,   170,   290,   490,   340}
    ,{   490,   140,   290,   490,   340}
    ,{   480,   170,   280,   480,   340}
    ,{   490,   140,   290,   490,   340}
    ,{  -600,  -600,  -690, -1150,  -640}
    }
   }
  ,{{{  1120,   970,  1120,   -60,  1120}
    ,{   580,   440,   580,   -60,   580}
    ,{   280,   140,   280,  -120,   280}
    ,{  1120,   970,  1120,  -350,  1120}
    ,{   280,   140,   280,  -120,   280}
    }
   ,{{   580,   440,   580,   -60,   580}
    ,{   580,   440,   580,   -60,   580}
    ,{   270,   130,   270,  -370,   270}
    ,{  -800,  -950,  -800, -1450,  -800}
    ,{   270,   130,   270,  -370,   270}
    }
   ,{{   290,   140,   290,  -120,   290}
    ,{   290,   140,   290,  -350,   290}
    ,{   280,   140,   280,  -120,   280}
    ,{   290,   140,   290,  -350,   290}
    ,{   280,   140,   280,  -120,   280}
    }
   ,{{  1120,   970,  1120,  -370,  1120}
    ,{  -620,  -770,  -620, -1270,  -620}
    ,{   270,   130,   270,  -370,   270}
    ,{  1120,   970,  1120,  -780,  1120}
    ,{   270,   130,   270,  -370,   270}
    }
   ,{{   290,   140,   290,  -120,   290}
    ,{   290,   140,   290,  -350,   290}
    ,{   280,   140,   280,  -120,   280}
    ,{   290,   140,   290,  -350,   290}
    ,{  -600,  -600,  -690, -1340,  -690}
    }
   }
  ,{{{  1170,   660,  1110,   660,  1170}
    ,{   640,   130,   580,   130,   640}
    ,{   340,  -170,   280,  -170,   340}
    ,{  1170,   660,  1110,   660,  1170}
    ,{   340,  -170,   280,  -170,   340}
    }
   ,{{   640,   130,   580,   130,   640}
    ,{   640,   130,   580,   130,   640}
    ,{   330,  -180,   270,  -180,   330}
    ,{  -510, -1260,  -570, -1260,  -510}
    ,{   330,  -180,   270,  -180,   330}
    }
   ,{{   340,  -160,   280,  -160,   340}
    ,{   340,  -160,   280,  -160,   340}
    ,{   340,  -170,   280,  -170,   340}
    ,{   340,  -160,   280,  -160,   340}
    ,{   340,  -170,   280,  -170,   340}
    }
   ,{{  1170,   660,  1110,   660,  1170}
    ,{  -330, -1080,  -390, -1080,  -330}
    ,{   330,  -180,   270,  -180,   330}
    ,{  1170,   660,  1110,   660,  1170}
    ,{   330,  -180,   270,  -180,   330}
    }
   ,{{   340,  -160,   280,  -160,   340}
    ,{   340,  -160,   280,  -160,   340}
    ,{   340,  -170,   280,  -170,   340}
    ,{   340,  -160,   280,  -160,   340}
    ,{  -640, -1150,  -700, -1150,  -640}
    }
   }
  ,{{{  1120,   220,  1120,   780,  1120}
    ,{   780,   220,   580,   780,   580}
    ,{   480,   170,   280,   480,   280}
    ,{  1120,   -70,  1120,   490,  1120}
    ,{   480,   170,   280,   480,   280}
    }
   ,{{   780,   220,   580,   780,   580}
    ,{   780,   220,   580,   780,   580}
    ,{   470,   -80,   270,   470,   270}
    ,{  -800, -1160,  -800, -1860,  -800}
    ,{   470,   -80,   270,   470,   270}
    }
   ,{{   490,   170,   290,   490,   290}
    ,{   490,   -70,   290,   490,   290}
    ,{   480,   170,   280,   480,   280}
    ,{   490,   -70,   290,   490,   290}
    ,{   480,   170,   280,   480,   280}
    }
   ,{{  1120,   -80,  1120,   470,  1120}
    ,{  -620,  -980,  -620, -1680,  -620}
    ,{   470,   -80,   270,   470,   270}
    ,{  1120,  -490,  1120, -1190,  1120}
    ,{   470,   -80,   270,   470,   270}
    }
   ,{{   490,   170,   290,   490,   290}
    ,{   490,   -70,   290,   490,   290}
    ,{   480,   170,   280,   480,   280}
    ,{   490,   -70,   290,   490,   290}
    ,{  -690, -1050,  -690, -1750,  -690}
    }
   }
  ,{{{  1060,   660,  1060,   660,    40}
    ,{   530,   130,   530,   130,    40}
    ,{   230,  -170,   230,  -170,  -500}
    ,{  1060,   660,  1060,   660,  -500}
    ,{   230,  -170,   230,  -170,  -500}
    }
   ,{{   530,   130,   530,   130,    40}
    ,{   530,   130,   530,   130,    40}
    ,{   220,  -180,   220,  -180,  -510}
    ,{  -620, -1260,  -620, -1260, -1590}
    ,{   220,  -180,   220,  -180,  -510}
    }
   ,{{   230,  -160,   230,  -160,  -500}
    ,{   230,  -160,   230,  -160,  -500}
    ,{   230,  -170,   230,  -170,  -500}
    ,{   230,  -160,   230,  -160,  -500}
    ,{   230,  -170,   230,  -170,  -500}
    }
   ,{{  1060,   660,  1060,   660,  -510}
    ,{  -440, -1080,  -440, -1080, -1410}
    ,{   220,  -180,   220,  -180,  -510}
    ,{  1060,   660,  1060,   660,  -920}
    ,{   220,  -180,   220,  -180,  -510}
    }
   ,{{   230,  -160,   230,  -160,  -500}
    ,{   230,  -160,   230,  -160,  -500}
    ,{   230,  -170,   230,  -170,  -500}
    ,{   230,  -160,   230,  -160,  -500}
    ,{  -750, -1150,  -750, -1150, -1480}
    }
   }
  }
 ,{{{{  1350,  1160,  1300,   850,  1350}
    ,{   850,   500,   650,   850,   700}
    ,{   720,   400,   520,   720,   570}
    ,{  1350,  1160,  1300,   850,  1350}
    ,{   590,   270,   390,   590,   440}
    }
   ,{{   850,   500,   650,   850,   700}
    ,{   850,   500,   650,   850,   700}
    ,{   570,   220,   370,   570,   420}
    ,{  -460,  -900,  -520, -1210,  -460}
    ,{   570,   220,   370,   570,   420}
    }
   ,{{   720,   400,   520,   720,   570}
    ,{   720,   370,   520,   720,   570}
    ,{   720,   400,   520,   720,   570}
    ,{   720,   370,   520,   720,   570}
    ,{   590,   270,   390,   590,   440}
    }
   ,{{  1350,  1160,  1300,   850,  1350}
    ,{  -760, -1200,  -820, -1510,  -760}
    ,{   570,   220,   370,   570,   420}
    ,{  1350,  1160,  1300,   850,  1350}
    ,{   570,   220,   370,   570,   420}
    }
   ,{{   720,   370,   520,   720,   570}
    ,{   720,   370,   520,   720,   570}
    ,{   280,   -40,    80,   280,   130}
    ,{   720,   370,   520,   720,   570}
    ,{  -320,  -320,  -420,  -870,  -360}
    }
   }
  ,{{{  1300,  1160,  1300,   120,  1300}
    ,{   650,   500,   650,     0,   650}
    ,{   520,   370,   520,   120,   520}
    ,{  1300,  1160,  1300,  -120,  1300}
    ,{   390,   240,   390,   -10,   390}
    }
   ,{{   650,   500,   650,     0,   650}
    ,{   650,   500,   650,     0,   650}
    ,{   370,   220,   370,  -270,   370}
    ,{  -750,  -900,  -750, -1400,  -750}
    ,{   370,   220,   370,  -270,   370}
    }
   ,{{   520,   370,   520,   120,   520}
    ,{   520,   370,   520,  -120,   520}
    ,{   520,   370,   520,   120,   520}
    ,{   520,   370,   520,  -120,   520}
    ,{   390,   240,   390,   -10,   390}
    }
   ,{{  1300,  1160,  1300,  -270,  1300}
    ,{ -1050, -1200, -1050, -1700, -1050}
    ,{   370,   220,   370,  -270,   370}
    ,{  1300,  1160,  1300,  -590,  1300}
    ,{   370,   220,   370,  -270,   370}
    }
   ,{{   520,   370,   520,  -120,   520}
    ,{   520,   370,   520,  -120,   520}
    ,{    80,   -60,    80,  -320,    80}
    ,{   520,   370,   520,  -120,   520}
    ,{  -320,  -320,  -420, -1060,  -420}
    }
   }
  ,{{{  1350,   850,  1290,   850,  1350}
    ,{   700,   190,   640,   190,   700}
    ,{   570,    60,   510,    60,   570}
    ,{  1350,   850,  1290,   850,  1350}
    ,{   440,   -60,   380,   -60,   440}
    }
   ,{{   700,   190,   640,   190,   700}
    ,{   700,   190,   640,   190,   700}
    ,{   420,   -80,   360,   -80,   420}
    ,{  -460, -1210,  -520, -1210,  -460}
    ,{   420,   -80,   360,   -80,   420}
    }
   ,{{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   440,   -60,   380,   -60,   440}
    }
   ,{{  1350,   850,  1290,   850,  1350}
    ,{  -760, -1510,  -820, -1510,  -760}
    ,{   420,   -80,   360,   -80,   420}
    ,{  1350,   850,  1290,   850,  1350}
    ,{   420,   -80,   360,   -80,   420}
    }
   ,{{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   130,  -370,    70,  -370,   130}
    ,{   570,    60,   510,    60,   570}
    ,{  -360,  -870,  -420,  -870,  -360}
    }
   }
  ,{{{  1300,   400,  1300,   850,  1300}
    ,{   850,   290,   650,   850,   650}
    ,{   720,   400,   520,   720,   520}
    ,{  1300,   160,  1300,   720,  1300}
    ,{   590,   270,   390,   590,   390}
    }
   ,{{   850,   290,   650,   850,   650}
    ,{   850,   290,   650,   850,   650}
    ,{   570,    10,   370,   570,   370}
    ,{  -750, -1110,  -750, -1810,  -750}
    ,{   570,    10,   370,   570,   370}
    }
   ,{{   720,   400,   520,   720,   520}
    ,{   720,   160,   520,   720,   520}
    ,{   720,   400,   520,   720,   520}
    ,{   720,   160,   520,   720,   520}
    ,{   590,   270,   390,   590,   390}
    }
   ,{{  1300,    10,  1300,   570,  1300}
    ,{ -1050, -1410, -1050, -2110, -1050}
    ,{   570,    10,   370,   570,   370}
    ,{  1300,  -310,  1300, -1000,  1300}
    ,{   570,    10,   370,   570,   370}
    }
   ,{{   720,   160,   520,   720,   520}
    ,{   720,   160,   520,   720,   520}
    ,{   280,   -40,    80,   280,    80}
    ,{   720,   160,   520,   720,   520}
    ,{  -420,  -780,  -420, -1470,  -420}
    }
   }
  ,{{{  1250,   850,  1250,   850,   100}
    ,{   590,   190,   590,   190,   100}
    ,{   460,    60,   460,    60,  -270}
    ,{  1250,   850,  1250,   850,  -270}
    ,{   330,   -60,   330,   -60,  -400}
    }
   ,{{   590,   190,   590,   190,   100}
    ,{   590,   190,   590,   190,   100}
    ,{   310,   -80,   310,   -80,  -420}
    ,{  -570, -1210,  -570, -1210, -1540}
    ,{   310,   -80,   310,   -80,  -420}
    }
   ,{{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{   330,   -60,   330,   -60,  -400}
    }
   ,{{  1250,   850,  1250,   850,  -420}
    ,{  -870, -1510,  -870, -1510, -1840}
    ,{   310,   -80,   310,   -80,  -420}
    ,{  1250,   850,  1250,   850,  -740}
    ,{   310,   -80,   310,   -80,  -420}
    }
   ,{{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{    20,  -370,    20,  -370,  -710}
    ,{   460,    60,   460,    60,  -270}
    ,{  -470,  -870,  -470,  -870, -1210}
    }
   }
  }
 ,{{{{  1350,  1160,  1300,   850,  1350}
    ,{   850,   500,   650,   850,   700}
    ,{   720,   400,   520,   720,   570}
    ,{  1350,  1160,  1300,   850,  1350}
    ,{   590,   270,   390,   590,   440}
    }
   ,{{   850,   500,   650,   850,   700}
    ,{   850,   500,   650,   850,   700}
    ,{   570,   220,   370,   570,   420}
    ,{  -230,  -670,  -290,  -980,  -230}
    ,{   570,   220,   370,   570,   420}
    }
   ,{{   720,   400,   520,   720,   570}
    ,{   720,   370,   520,   720,   570}
    ,{   720,   400,   520,   720,   570}
    ,{   720,   370,   520,   720,   570}
    ,{   590,   270,   390,   590,   440}
    }
   ,{{  1350,  1160,  1300,   850,  1350}
    ,{  -330,  -770,  -390, -1080,  -330}
    ,{   570,   220,   370,   570,   420}
    ,{  1350,  1160,  1300,   850,  1350}
    ,{   570,   220,   370,   570,   420}
    }
   ,{{   720,   370,   520,   720,   570}
    ,{   720,   370,   520,   720,   570}
    ,{   480,   170,   280,   480,   340}
    ,{   720,   370,   520,   720,   570}
    ,{   -90,  -320,   -90,  -810,  -360}
    }
   }
  ,{{{  1300,  1160,  1300,   540,  1300}
    ,{   650,   500,   650,    10,   650}
    ,{   540,   370,   520,   540,   520}
    ,{  1300,  1160,  1300,  -120,  1300}
    ,{   390,   240,   390,   -10,   390}
    }
   ,{{   650,   500,   650,     0,   650}
    ,{   650,   500,   650,     0,   650}
    ,{   370,   220,   370,  -100,   370}
    ,{  -530,  -670,  -530, -1170,  -530}
    ,{   370,   220,   370,  -270,   370}
    }
   ,{{   540,   370,   520,   540,   520}
    ,{   520,   370,   520,    10,   520}
    ,{   540,   370,   520,   540,   520}
    ,{   520,   370,   520,  -120,   520}
    ,{   390,   240,   390,   -10,   390}
    }
   ,{{  1300,  1160,  1300,  -270,  1300}
    ,{  -620,  -770,  -620, -1270,  -620}
    ,{   370,   220,   370,  -270,   370}
    ,{  1300,  1160,  1300,  -590,  1300}
    ,{   370,   220,   370,  -270,   370}
    }
   ,{{   520,   370,   520,  -120,   520}
    ,{   520,   370,   520,  -120,   520}
    ,{   280,   140,   280,  -120,   280}
    ,{   520,   370,   520,  -120,   520}
    ,{   -90,  -320,   -90,  -810,  -420}
    }
   }
  ,{{{  1350,   850,  1290,   850,  1350}
    ,{   700,   190,   640,   190,   700}
    ,{   570,    60,   510,    60,   570}
    ,{  1350,   850,  1290,   850,  1350}
    ,{   440,   -60,   380,   -60,   440}
    }
   ,{{   700,   190,   640,   190,   700}
    ,{   700,   190,   640,   190,   700}
    ,{   420,   -80,   360,   -80,   420}
    ,{  -230,  -980,  -290,  -980,  -230}
    ,{   420,   -80,   360,   -80,   420}
    }
   ,{{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   440,   -60,   380,   -60,   440}
    }
   ,{{  1350,   850,  1290,   850,  1350}
    ,{  -330, -1070,  -390, -1080,  -330}
    ,{   420,   -80,   360,   -80,   420}
    ,{  1350,   850,  1290,   850,  1350}
    ,{   420,   -80,   360,   -80,   420}
    }
   ,{{   570,    60,   510,    60,   570}
    ,{   570,    60,   510,    60,   570}
    ,{   340,  -170,   280,  -170,   340}
    ,{   570,    60,   510,    60,   570}
    ,{  -360,  -830,  -420,  -870,  -360}
    }
   }
  ,{{{  1300,   400,  1300,   850,  1300}
    ,{   850,   290,   650,   850,   650}
    ,{   720,   400,   520,   720,   520}
    ,{  1300,   160,  1300,   720,  1300}
    ,{   590,   270,   390,   590,   390}
    }
   ,{{   850,   290,   650,   850,   650}
    ,{   850,   290,   650,   850,   650}
    ,{   570,    10,   370,   570,   370}
    ,{  -530,  -890,  -530, -1580,  -530}
    ,{   570,    10,   370,   570,   370}
    }
   ,{{   720,   400,   520,   720,   520}
    ,{   720,   160,   520,   720,   520}
    ,{   720,   400,   520,   720,   520}
    ,{   720,   160,   520,   720,   520}
    ,{   590,   270,   390,   590,   390}
    }
   ,{{  1300,    10,  1300,   570,  1300}
    ,{  -620,  -980,  -620, -1080,  -620}
    ,{   570,    10,   370,   570,   370}
    ,{  1300,  -310,  1300, -1000,  1300}
    ,{   570,    10,   370,   570,   370}
    }
   ,{{   720,   170,   520,   720,   520}
    ,{   720,   160,   520,   720,   520}
    ,{   480,   170,   280,   480,   280}
    ,{   720,   160,   520,   720,   520}
    ,{  -420,  -780,  -420, -1470,  -420}
    }
   }
  ,{{{  1250,   850,  1250,   850,   100}
    ,{   590,   190,   590,   190,   100}
    ,{   460,    60,   460,    60,  -270}
    ,{  1250,   850,  1250,   850,  -230}
    ,{   330,   -60,   330,   -60,  -400}
    }
   ,{{   590,   190,   590,   190,   100}
    ,{   590,   190,   590,   190,   100}
    ,{   310,   -80,   310,   -80,  -420}
    ,{  -340,  -980,  -340,  -980, -1320}
    ,{   310,   -80,   310,   -80,  -420}
    }
   ,{{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{   330,   -60,   330,   -60,  -400}
    }
   ,{{  1250,   850,  1250,   850,  -230}
    ,{  -440, -1080,  -440, -1080, -1300}
    ,{   310,   -80,   310,   -80,  -420}
    ,{  1250,   850,  1250,   850,  -230}
    ,{   310,   -80,   310,   -80,  -420}
    }
   ,{{   460,    60,   460,    60,  -270}
    ,{   460,    60,   460,    60,  -270}
    ,{   230,  -170,   230,  -170,  -500}
    ,{   460,    60,   460,    60,  -270}
    ,{  -470,  -870,  -470,  -870, -1210}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   540,   -90,   540,   180,   -90}
    ,{   540,  -100,   540,   180,   -90}
    ,{   180,   -90,  -460,   180,  -460}
    ,{    30,  -150,  -260,    30,  -210}
    ,{  -200,  -200,  -400,  -230,  -570}
    }
   ,{{   180,  -350,  -660,   180,  -660}
    ,{   180,  -580,  -660,   180,  -660}
    ,{  -430,  -600,  -970,  -430,  -830}
    ,{  -350,  -350,  -870,  -960,  -870}
    ,{  -430,  -600,  -970,  -430,  -970}
    }
   ,{{    30,  -150,  -510,    30,   -90}
    ,{   -90,  -220,  -510,  -390,   -90}
    ,{    20,  -600,  -520,    20,  -520}
    ,{    30,  -150,  -510,    30,  -510}
    ,{  -200,  -200,  -570,  -650,  -570}
    }
   ,{{   540,  -100,   540,  -400,  -210}
    ,{   540,  -100,   540, -1240,  -810}
    ,{  -430,  -600,  -970,  -430,  -970}
    ,{  -200,  -200,  -260,  -400,  -210}
    ,{  -430,  -600,  -970,  -430,  -970}
    }
   ,{{   180,   -90,  -400,   180,  -460}
    ,{    30,  -150,  -510,    30,  -510}
    ,{   180,   -90,  -460,   180,  -460}
    ,{    30,  -150,  -510,    30,  -510}
    ,{  -230, -1390,  -400,  -230, -1300}
    }
   }
  ,{{{    10,   -90,    10,  -500,  -320}
    ,{    10,  -150,    10,  -860,  -510}
    ,{   -90,   -90,  -460,  -500,  -460}
    ,{  -150,  -150,  -320,  -860,  -320}
    ,{  -200,  -200,  -400, -1300,  -570}
    }
   ,{{  -580,  -580,  -660, -1070,  -660}
    ,{  -580,  -580,  -660, -1340,  -660}
    ,{  -600,  -600,  -970, -1070,  -970}
    ,{  -870, -1600, -1110, -1880,  -870}
    ,{  -600,  -600,  -970, -1320,  -970}
    }
   ,{{  -150,  -150,  -510,  -500,  -510}
    ,{  -220,  -220, -1150,  -860,  -510}
    ,{  -500, -1070,  -750,  -500,  -520}
    ,{  -150,  -150,  -510,  -860,  -510}
    ,{  -200,  -200,  -570, -1750,  -570}
    }
   ,{{    10,  -200,    10, -1080,  -320}
    ,{    10,  -970,    10, -2450, -1160}
    ,{  -600,  -600,  -970, -1320,  -970}
    ,{  -200,  -200,  -320, -1080,  -320}
    ,{  -600,  -600,  -970, -1320,  -970}
    }
   ,{{   -90,   -90,  -400,  -570,  -460}
    ,{  -150,  -150,  -510,  -860,  -510}
    ,{   -90,   -90,  -460,  -570,  -460}
    ,{  -150,  -150,  -510,  -860,  -510}
    ,{  -400, -1490,  -400, -1300, -1300}
    }
   }
  ,{{{   540,  -100,   540,  -400,  -210}
    ,{   540,  -100,   540,  -600, -1130}
    ,{  -540,  -540,  -760,  -540, -1070}
    ,{  -210,  -350,  -620,  -400,  -210}
    ,{  -650,  -650,  -870,  -650, -1180}
    }
   ,{{  -350,  -350,  -940,  -740, -1250}
    ,{  -740,  -740,  -960,  -740, -1270}
    ,{ -1050, -1050, -1270, -1050, -1580}
    ,{  -350,  -350,  -940,  -960, -1250}
    ,{ -1050, -1050, -1270, -1050, -1580}
    }
   ,{{  -600,  -600,  -820,  -600, -1130}
    ,{  -600,  -600,  -820,  -600, -1130}
    ,{  -600,  -600,  -820,  -600, -1130}
    ,{  -600,  -600,  -820,  -600, -1130}
    ,{  -650,  -650,  -870,  -650, -1180}
    }
   ,{{   540,  -100,   540,  -400,  -210}
    ,{   540,  -100,   540, -1240, -1530}
    ,{ -1050, -1050, -1270, -1050, -1580}
    ,{  -210,  -440,  -620,  -400,  -210}
    ,{ -1050, -1050, -1270, -1050, -1580}
    }
   ,{{  -540,  -540,  -760,  -540, -1070}
    ,{  -600,  -600,  -820,  -600, -1130}
    ,{  -540,  -540,  -760,  -540, -1070}
    ,{  -600,  -600,  -820,  -600, -1130}
    ,{ -1390, -1390, -1610, -1390, -1920}
    }
   }
  ,{{{   180,  -630,  -320,   180,  -320}
    ,{   180, -1340,  -510,   180,  -510}
    ,{   180,  -630,  -460,   180,  -460}
    ,{    30, -1340,  -320,    30,  -320}
    ,{  -230, -1150,  -570,  -230,  -570}
    }
   ,{{   180, -1790,  -660,   180,  -660}
    ,{   180, -2010,  -660,   180,  -660}
    ,{  -430, -1790,  -970,  -430,  -970}
    ,{  -870, -3070,  -870, -1370,  -870}
    ,{  -430, -1790,  -970,  -430,  -970}
    }
   ,{{    30,  -630,  -510,    30,  -510}
    ,{  -390, -1650,  -510,  -390,  -510}
    ,{    20,  -630,  -520,    20,  -520}
    ,{    30, -1340,  -510,    30,  -510}
    ,{  -570, -1150,  -570,  -880,  -570}
    }
   ,{{  -320, -1790,  -320,  -430,  -320}
    ,{ -1160, -1980, -1160, -1870, -1160}
    ,{  -430, -1790,  -970,  -430,  -970}
    ,{  -320, -2390,  -320, -2280,  -320}
    ,{  -430, -1790,  -970,  -430,  -970}
    }
   ,{{   180, -1040,  -460,   180,  -460}
    ,{    30, -1340,  -510,    30,  -510}
    ,{   180, -1040,  -460,   180,  -460}
    ,{    30, -1340,  -510,    30,  -510}
    ,{  -230, -1520, -1300,  -230, -1300}
    }
   }
  ,{{{   -90,  -400,  -260,  -400,   -90}
    ,{   -90,  -600,  -820,  -600,   -90}
    ,{  -540,  -540,  -550,  -540,  -830}
    ,{  -260,  -400,  -260,  -400,  -800}
    ,{  -650,  -650,  -870,  -650,  -860}
    }
   ,{{  -740,  -740,  -940,  -740,  -830}
    ,{  -740,  -740,  -960,  -740, -1240}
    ,{  -830, -1050, -1270, -1050,  -830}
    ,{  -940,  -960,  -940,  -960, -1360}
    ,{ -1050, -1050, -1270, -1050, -1260}
    }
   ,{{   -90,  -600,  -820,  -600,   -90}
    ,{   -90,  -600,  -820,  -600,   -90}
    ,{  -600,  -600,  -820,  -600, -1710}
    ,{  -600,  -600,  -820,  -600,  -800}
    ,{  -650,  -650,  -870,  -650,  -860}
    }
   ,{{  -260,  -400,  -260,  -400,  -810}
    ,{  -810, -1240, -1220, -1240,  -810}
    ,{ -1050, -1050, -1270, -1050, -1260}
    ,{  -260,  -400,  -260,  -400, -1550}
    ,{ -1050, -1050, -1270, -1050, -1260}
    }
   ,{{  -540,  -540,  -550,  -540,  -800}
    ,{  -600,  -600,  -820,  -600,  -800}
    ,{  -540,  -540,  -550,  -540, -1460}
    ,{  -600,  -600,  -820,  -600,  -800}
    ,{ -1390, -1390, -1610, -1390, -2350}
    }
   }
  }
 ,{{{{    50,    50,  -320,    50,  -320}
    ,{    50,  -130,  -490,    50,  -490}
    ,{  -400,  -580,  -940,  -400,  -940}
    ,{    50,    50,  -320,  -320,  -320}
    ,{  -400,  -540,  -940,  -400,  -940}
    }
   ,{{    50,  -130,  -490,    50,  -490}
    ,{    50,  -130,  -490,    50,  -490}
    ,{  -400,  -580,  -940,  -400,  -940}
    ,{ -1320, -1320, -1680, -1770, -1680}
    ,{  -400,  -580,  -940,  -400,  -940}
    }
   ,{{  -320,  -490,  -860,  -320,  -860}
    ,{  -320,  -490,  -860,  -320,  -860}
    ,{  -620,  -800, -1160,  -620, -1160}
    ,{  -320,  -490,  -860,  -320,  -860}
    ,{  -620,  -800, -1160,  -620, -1160}
    }
   ,{{    50,    50,  -320,  -400,  -320}
    ,{  -840,  -840, -1210, -1290, -1210}
    ,{  -400,  -580,  -940,  -400,  -940}
    ,{    50,    50,  -320,  -400,  -320}
    ,{  -400,  -580,  -940,  -400,  -940}
    }
   ,{{  -320,  -490,  -860,  -320,  -860}
    ,{  -320,  -490,  -860,  -320,  -860}
    ,{  -930, -1110, -1470,  -930, -1470}
    ,{  -320,  -490,  -860,  -320,  -860}
    ,{  -540,  -540, -1150, -1230, -1150}
    }
   }
  ,{{{    50,    50,  -320,  -840,  -320}
    ,{  -130,  -130,  -490,  -840,  -490}
    ,{  -580,  -580,  -940, -1270,  -940}
    ,{    50,    50,  -320, -1210,  -320}
    ,{  -540,  -540,  -940, -1270,  -940}
    }
   ,{{  -130,  -130,  -490,  -840,  -490}
    ,{  -130,  -130,  -490,  -840,  -490}
    ,{  -580,  -580,  -940, -1290,  -940}
    ,{ -1320, -1320, -1680, -2030, -1680}
    ,{  -580,  -580,  -940, -1290,  -940}
    }
   ,{{  -490,  -490,  -860, -1210,  -860}
    ,{  -490,  -490,  -860, -1210,  -860}
    ,{  -800,  -800, -1160, -1270, -1160}
    ,{  -490,  -490,  -860, -1210,  -860}
    ,{  -800,  -800, -1160, -1270, -1160}
    }
   ,{{    50,    50,  -320, -1290,  -320}
    ,{  -840,  -840, -1210, -1560, -1210}
    ,{  -580,  -580,  -940, -1290,  -940}
    ,{    50,    50,  -320, -1920,  -320}
    ,{  -580,  -580,  -940, -1290,  -940}
    }
   ,{{  -490,  -490,  -860, -1210,  -860}
    ,{  -490,  -490,  -860, -1210,  -860}
    ,{ -1110, -1110, -1470, -1580, -1470}
    ,{  -490,  -490,  -860, -1210,  -860}
    ,{  -540,  -540, -1150, -1500, -1150}
    }
   }
  ,{{{  -400,  -400,  -620,  -400,  -930}
    ,{  -580,  -580,  -800,  -580, -1110}
    ,{ -1030, -1030, -1250, -1030, -1560}
    ,{  -400,  -400,  -620,  -400,  -930}
    ,{ -1030, -1030, -1250, -1030, -1560}
    }
   ,{{  -580,  -580,  -800,  -580, -1110}
    ,{  -580,  -580,  -800,  -580, -1110}
    ,{ -1030, -1030, -1250, -1030, -1560}
    ,{ -1750, -1770, -1750, -1770, -2060}
    ,{ -1030, -1030, -1250, -1030, -1560}
    }
   ,{{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{ -1250, -1250, -1470, -1250, -1780}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{ -1250, -1250, -1470, -1250, -1780}
    }
   ,{{  -400,  -400,  -620,  -400,  -930}
    ,{ -1270, -1290, -1270, -1290, -1580}
    ,{ -1030, -1030, -1250, -1030, -1560}
    ,{  -400,  -400,  -620,  -400,  -930}
    ,{ -1030, -1030, -1250, -1030, -1560}
    }
   ,{{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{ -1560, -1560, -1780, -1560, -2090}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{ -1230, -1230, -1450, -1230, -1760}
    }
   }
  ,{{{    50, -1320,  -320,    50,  -320}
    ,{    50, -1320,  -490,    50,  -490}
    ,{  -400, -1750,  -940,  -400,  -940}
    ,{  -320, -1680,  -320,  -320,  -320}
    ,{  -400, -1750,  -940,  -400,  -940}
    }
   ,{{    50, -1320,  -490,    50,  -490}
    ,{    50, -1320,  -490,    50,  -490}
    ,{  -400, -1770,  -940,  -400,  -940}
    ,{ -1680, -2510, -1680, -2390, -1680}
    ,{  -400, -1770,  -940,  -400,  -940}
    }
   ,{{  -320, -1680,  -860,  -320,  -860}
    ,{  -320, -1680,  -860,  -320,  -860}
    ,{  -620, -1750, -1160,  -620, -1160}
    ,{  -320, -1680,  -860,  -320,  -860}
    ,{  -620, -1750, -1160,  -620, -1160}
    }
   ,{{  -320, -1770,  -320,  -400,  -320}
    ,{ -1210, -2030, -1210, -1920, -1210}
    ,{  -400, -1770,  -940,  -400,  -940}
    ,{  -320, -2390,  -320, -2280,  -320}
    ,{  -400, -1770,  -940,  -400,  -940}
    }
   ,{{  -320, -1680,  -860,  -320,  -860}
    ,{  -320, -1680,  -860,  -320,  -860}
    ,{  -930, -2060, -1470,  -930, -1470}
    ,{  -320, -1680,  -860,  -320,  -860}
    ,{ -1150, -1970, -1150, -1860, -1150}
    }
   }
  ,{{{  -400,  -400,  -620,  -400,  -540}
    ,{  -540,  -580,  -800,  -580,  -540}
    ,{ -1030, -1030, -1250, -1030, -1230}
    ,{  -400,  -400,  -620,  -400, -1150}
    ,{ -1030, -1030, -1250, -1030, -1230}
    }
   ,{{  -540,  -580,  -800,  -580,  -540}
    ,{  -540,  -580,  -800,  -580,  -540}
    ,{ -1030, -1030, -1250, -1030, -1230}
    ,{ -1750, -1770, -1750, -1770, -1970}
    ,{ -1030, -1030, -1250, -1030, -1230}
    }
   ,{{  -940,  -940, -1160,  -940, -1150}
    ,{  -940,  -940, -1160,  -940, -1150}
    ,{ -1250, -1250, -1470, -1250, -1450}
    ,{  -940,  -940, -1160,  -940, -1150}
    ,{ -1250, -1250, -1470, -1250, -1450}
    }
   ,{{  -400,  -400,  -620,  -400, -1230}
    ,{ -1270, -1290, -1270, -1290, -1500}
    ,{ -1030, -1030, -1250, -1030, -1230}
    ,{  -400,  -400,  -620,  -400, -1860}
    ,{ -1030, -1030, -1250, -1030, -1230}
    }
   ,{{  -940,  -940, -1160,  -940, -1150}
    ,{  -940,  -940, -1160,  -940, -1150}
    ,{ -1560, -1560, -1780, -1560, -1760}
    ,{  -940,  -940, -1160,  -940, -1150}
    ,{ -1230, -1230, -1450, -1230, -1440}
    }
   }
  }
 ,{{{{   210,   210,  -160,  -240,  -160}
    ,{  -870,  -870, -1230,  -870, -1230}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{   210,   210,  -160,  -240,  -160}
    ,{  -800,  -800, -1410,  -870, -1410}
    }
   ,{{  -870, -1040, -1410,  -870, -1410}
    ,{ -1050, -1220, -1590, -1050, -1590}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{ -1060, -1060, -1420, -1510, -1420}
    ,{  -870, -1040, -1410,  -870, -1410}
    }
   ,{{  -870, -1040, -1410,  -870, -1410}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{  -870, -1040, -1410,  -870, -1410}
    }
   ,{{   210,   210,  -160,  -240,  -160}
    ,{  -870,  -870, -1230, -1320, -1230}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{   210,   210,  -160,  -240,  -160}
    ,{  -870, -1040, -1410,  -870, -1410}
    }
   ,{{  -800,  -800, -1410,  -870, -1410}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{  -870, -1040, -1410,  -870, -1410}
    ,{  -800,  -800, -1410, -1490, -1410}
    }
   }
  ,{{{   210,   210,  -160, -1520,  -160}
    ,{  -870,  -870, -1230, -1580, -1230}
    ,{ -1040, -1040, -1410, -1520, -1410}
    ,{   210,   210,  -160, -1760,  -160}
    ,{  -800,  -800, -1410, -1520, -1410}
    }
   ,{{ -1040, -1040, -1410, -1760, -1410}
    ,{ -1220, -1220, -1590, -1940, -1590}
    ,{ -1040, -1040, -1410, -1760, -1410}
    ,{ -1060, -1060, -1420, -1770, -1420}
    ,{ -1040, -1040, -1410, -1760, -1410}
    }
   ,{{ -1040, -1040, -1410, -1520, -1410}
    ,{ -1040, -1040, -1410, -1760, -1410}
    ,{ -1040, -1040, -1410, -1520, -1410}
    ,{ -1040, -1040, -1410, -1760, -1410}
    ,{ -1040, -1040, -1410, -1520, -1410}
    }
   ,{{   210,   210,  -160, -1580,  -160}
    ,{  -870,  -870, -1230, -1580, -1230}
    ,{ -1040, -1040, -1410, -1760, -1410}
    ,{   210,   210,  -160, -1760,  -160}
    ,{ -1040, -1040, -1410, -1760, -1410}
    }
   ,{{  -800,  -800, -1410, -1520, -1410}
    ,{ -1040, -1040, -1410, -1760, -1410}
    ,{ -1040, -1040, -1410, -1520, -1410}
    ,{ -1040, -1040, -1410, -1760, -1410}
    ,{  -800,  -800, -1410, -1760, -1410}
    }
   }
  ,{{{  -240,  -240,  -460,  -240,  -770}
    ,{ -1300, -1320, -1300, -1320, -1610}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{  -240,  -240,  -460,  -240,  -770}
    ,{ -1490, -1490, -1710, -1490, -2020}
    }
   ,{{ -1490, -1490, -1490, -1490, -1800}
    ,{ -1670, -1670, -1890, -1670, -2200}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1510, -1490, -1510, -1800}
    ,{ -1490, -1490, -1710, -1490, -2020}
    }
   ,{{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    }
   ,{{  -240,  -240,  -460,  -240,  -770}
    ,{ -1300, -1320, -1300, -1320, -1610}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{  -240,  -240,  -460,  -240,  -770}
    ,{ -1490, -1490, -1710, -1490, -2020}
    }
   ,{{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    ,{ -1490, -1490, -1710, -1490, -2020}
    }
   }
  ,{{{  -160, -1990,  -160,  -870,  -160}
    ,{  -870, -2060, -1230,  -870, -1230}
    ,{  -870, -1990, -1410,  -870, -1410}
    ,{  -160, -2230,  -160,  -870,  -160}
    ,{  -870, -1990, -1410,  -870, -1410}
    }
   ,{{  -870, -2230, -1410,  -870, -1410}
    ,{ -1050, -2410, -1590, -1050, -1590}
    ,{  -870, -2230, -1410,  -870, -1410}
    ,{ -1420, -2250, -1420, -2130, -1420}
    ,{  -870, -2230, -1410,  -870, -1410}
    }
   ,{{  -870, -1990, -1410,  -870, -1410}
    ,{  -870, -2230, -1410,  -870, -1410}
    ,{  -870, -1990, -1410,  -870, -1410}
    ,{  -870, -2230, -1410,  -870, -1410}
    ,{  -870, -1990, -1410,  -870, -1410}
    }
   ,{{  -160, -2060,  -160,  -870,  -160}
    ,{ -1230, -2060, -1230, -1940, -1230}
    ,{  -870, -2230, -1410,  -870, -1410}
    ,{  -160, -2230,  -160, -2120,  -160}
    ,{  -870, -2230, -1410,  -870, -1410}
    }
   ,{{  -870, -1990, -1410,  -870, -1410}
    ,{  -870, -2230, -1410,  -870, -1410}
    ,{  -870, -1990, -1410,  -870, -1410}
    ,{  -870, -2230, -1410,  -870, -1410}
    ,{ -1410, -2230, -1410, -2120, -1410}
    }
   }
  ,{{{  -240,  -240,  -460,  -240, -1520}
    ,{ -1300, -1320, -1300, -1320, -1520}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{  -240,  -240,  -460,  -240, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    }
   ,{{ -1490, -1490, -1490, -1490, -1640}
    ,{ -1640, -1670, -1890, -1670, -1640}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1510, -1490, -1510, -1710}
    ,{ -1490, -1490, -1710, -1490, -1700}
    }
   ,{{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    }
   ,{{  -240,  -240,  -460,  -240, -1520}
    ,{ -1300, -1320, -1300, -1320, -1520}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{  -240,  -240,  -460,  -240, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    }
   ,{{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    ,{ -1490, -1490, -1710, -1490, -1700}
    }
   }
  }
 ,{{{{   760,   760,   400,   310,   400}
    ,{   200,  -430,  -340,   200,  -340}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{   760,   760,   400,   310,   400}
    ,{  -250,  -250,  -850,  -310,  -850}
    }
   ,{{   200,  -430,  -340,   200,  -340}
    ,{   200,  -430,  -340,   200,  -340}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{  -830,  -830, -1190, -1280, -1190}
    ,{  -310,  -490,  -850,  -310,  -850}
    }
   ,{{  -310,  -490,  -850,  -310,  -850}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{  -310,  -490,  -850,  -310,  -850}
    }
   ,{{   760,   760,   400,   310,   400}
    ,{ -1000, -1000, -1360, -1450, -1360}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{   760,   760,   400,   310,   400}
    ,{  -310,  -490,  -850,  -310,  -850}
    }
   ,{{  -250,  -250,  -850,  -310,  -850}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{  -310,  -490,  -850,  -310,  -850}
    ,{  -250,  -250,  -850,  -940,  -850}
    }
   }
  ,{{{   760,   760,   400,  -690,   400}
    ,{  -340,  -490,  -340,  -690,  -340}
    ,{  -490,  -490,  -850,  -960,  -850}
    ,{   760,   760,   400, -1200,   400}
    ,{  -250,  -250,  -850,  -960,  -850}
    }
   ,{{  -340,  -490,  -340,  -690,  -340}
    ,{  -340, -2040,  -340,  -690,  -340}
    ,{  -490,  -490,  -850, -1200,  -850}
    ,{  -830,  -830, -1190, -1540, -1190}
    ,{  -490,  -490,  -850, -1200,  -850}
    }
   ,{{  -490,  -490,  -850,  -960,  -850}
    ,{  -490,  -490,  -850, -1200,  -850}
    ,{  -490,  -490,  -850,  -960,  -850}
    ,{  -490,  -490,  -850, -1200,  -850}
    ,{  -490,  -490,  -850,  -960,  -850}
    }
   ,{{   760,   760,   400, -1200,   400}
    ,{ -1000, -1000, -1360, -1710, -1360}
    ,{  -490,  -490,  -850, -1200,  -850}
    ,{   760,   760,   400, -1200,   400}
    ,{  -490,  -490,  -850, -1200,  -850}
    }
   ,{{  -250,  -250,  -850,  -960,  -850}
    ,{  -490,  -490,  -850, -1200,  -850}
    ,{  -490,  -490,  -850,  -960,  -850}
    ,{  -490,  -490,  -850, -1200,  -850}
    ,{  -250,  -250,  -850, -1200,  -850}
    }
   }
  ,{{{   310,   310,    90,   310,  -220}
    ,{  -430,  -430,  -650,  -430,  -960}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{   310,   310,    90,   310,  -220}
    ,{  -940,  -940, -1160,  -940, -1470}
    }
   ,{{  -430,  -430,  -650,  -430,  -960}
    ,{  -430,  -430,  -650,  -430,  -960}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{ -1260, -1280, -1260, -1280, -1570}
    ,{  -940,  -940, -1160,  -940, -1470}
    }
   ,{{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    }
   ,{{   310,   310,    90,   310,  -220}
    ,{ -1430, -1450, -1430, -1450, -1740}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{   310,   310,    90,   310,  -220}
    ,{  -940,  -940, -1160,  -940, -1470}
    }
   ,{{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    ,{  -940,  -940, -1160,  -940, -1470}
    }
   }
  ,{{{   400, -1170,   400,   200,   400}
    ,{   200, -1170,  -340,   200,  -340}
    ,{  -310, -1440,  -850,  -310,  -850}
    ,{   400, -1680,   400,  -310,   400}
    ,{  -310, -1440,  -850,  -310,  -850}
    }
   ,{{   200, -1170,  -340,   200,  -340}
    ,{   200, -1170,  -340,   200,  -340}
    ,{  -310, -1680,  -850,  -310,  -850}
    ,{ -1190, -2020, -1190, -1900, -1190}
    ,{  -310, -1680,  -850,  -310,  -850}
    }
   ,{{  -310, -1440,  -850,  -310,  -850}
    ,{  -310, -1680,  -850,  -310,  -850}
    ,{  -310, -1440,  -850,  -310,  -850}
    ,{  -310, -1680,  -850,  -310,  -850}
    ,{  -310, -1440,  -850,  -310,  -850}
    }
   ,{{   400, -1680,   400,  -310,   400}
    ,{ -1360, -2190, -1360, -2070, -1360}
    ,{  -310, -1680,  -850,  -310,  -850}
    ,{   400, -1680,   400, -1560,   400}
    ,{  -310, -1680,  -850,  -310,  -850}
    }
   ,{{  -310, -1440,  -850,  -310,  -850}
    ,{  -310, -1680,  -850,  -310,  -850}
    ,{  -310, -1440,  -850,  -310,  -850}
    ,{  -310, -1680,  -850,  -310,  -850}
    ,{  -850, -1680,  -850, -1560,  -850}
    }
   }
  ,{{{   310,   310,    90,   310,  -390}
    ,{  -390,  -430,  -650,  -430,  -390}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{   310,   310,    90,   310, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    }
   ,{{  -390,  -430,  -650,  -430,  -390}
    ,{  -390,  -430,  -650,  -430,  -390}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{ -1260, -1280, -1260, -1280, -1480}
    ,{  -940,  -940, -1160,  -940, -1140}
    }
   ,{{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    }
   ,{{   310,   310,    90,   310, -1140}
    ,{ -1430, -1450, -1430, -1450, -1650}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{   310,   310,    90,   310, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    }
   ,{{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    ,{  -940,  -940, -1160,  -940, -1140}
    }
   }
  }
 ,{{{{  1140,  1140,   770,   780,   770}
    ,{   780,   600,   240,   780,   240}
    ,{   480,   300,   -60,   480,   -60}
    ,{  1140,  1140,   770,   690,   770}
    ,{   480,   300,   -60,   480,   -60}
    }
   ,{{   780,   600,   240,   780,   240}
    ,{   780,   600,   240,   780,   240}
    ,{   470,   290,   -70,   470,   -70}
    ,{  -780,  -780, -1150, -1230, -1150}
    ,{   470,   290,   -70,   470,   -70}
    }
   ,{{   490,   310,   -50,   490,   -50}
    ,{   490,   310,   -50,   490,   -50}
    ,{   480,   300,   -60,   480,   -60}
    ,{   490,   310,   -50,   490,   -50}
    ,{   480,   300,   -60,   480,   -60}
    }
   ,{{  1140,  1140,   770,   690,   770}
    ,{  -600,  -600,  -970, -1050,  -970}
    ,{   470,   290,   -70,   470,   -70}
    ,{  1140,  1140,   770,   690,   770}
    ,{   470,   290,   -70,   470,   -70}
    }
   ,{{   490,   310,   -50,   490,   -50}
    ,{   490,   310,   -50,   490,   -50}
    ,{   480,   300,   -60,   480,   -60}
    ,{   490,   310,   -50,   490,   -50}
    ,{  -430,  -430, -1040, -1120, -1040}
    }
   }
  ,{{{  1140,  1140,   770,  -110,   770}
    ,{   600,   600,   240,  -110,   240}
    ,{   300,   300,   -60,  -170,   -60}
    ,{  1140,  1140,   770,  -400,   770}
    ,{   300,   300,   -60,  -170,   -60}
    }
   ,{{   600,   600,   240,  -110,   240}
    ,{   600,   600,   240,  -110,   240}
    ,{   290,   290,   -70,  -420,   -70}
    ,{  -780,  -780, -1150, -1500, -1150}
    ,{   290,   290,   -70,  -420,   -70}
    }
   ,{{   310,   310,   -50,  -170,   -50}
    ,{   310,   310,   -50,  -400,   -50}
    ,{   300,   300,   -60,  -170,   -60}
    ,{   310,   310,   -50,  -400,   -50}
    ,{   300,   300,   -60,  -170,   -60}
    }
   ,{{  1140,  1140,   770,  -420,   770}
    ,{  -600,  -600,  -970, -1320,  -970}
    ,{   290,   290,   -70,  -420,   -70}
    ,{  1140,  1140,   770,  -830,   770}
    ,{   290,   290,   -70,  -420,   -70}
    }
   ,{{   310,   310,   -50,  -170,   -50}
    ,{   310,   310,   -50,  -400,   -50}
    ,{   300,   300,   -60,  -170,   -60}
    ,{   310,   310,   -50,  -400,   -50}
    ,{  -430,  -430, -1040, -1390, -1040}
    }
   }
  ,{{{   690,   690,   470,   690,   160}
    ,{   150,   150,   -60,   150,  -370}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{   690,   690,   470,   690,   160}
    ,{  -140,  -140,  -360,  -140,  -670}
    }
   ,{{   150,   150,   -60,   150,  -370}
    ,{   150,   150,   -60,   150,  -370}
    ,{  -150,  -150,  -370,  -150,  -680}
    ,{ -1210, -1230, -1210, -1230, -1520}
    ,{  -150,  -150,  -370,  -150,  -680}
    }
   ,{{  -140,  -140,  -360,  -140,  -670}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{  -140,  -140,  -360,  -140,  -670}
    }
   ,{{   690,   690,   470,   690,   160}
    ,{ -1030, -1050, -1030, -1050, -1340}
    ,{  -150,  -150,  -370,  -150,  -680}
    ,{   690,   690,   470,   690,   160}
    ,{  -150,  -150,  -370,  -150,  -680}
    }
   ,{{  -140,  -140,  -360,  -140,  -670}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{ -1120, -1120, -1340, -1120, -1650}
    }
   }
  ,{{{   780,  -580,   770,   780,   770}
    ,{   780,  -580,   240,   780,   240}
    ,{   480,  -640,   -60,   480,   -60}
    ,{   770,  -880,   770,   490,   770}
    ,{   480,  -640,   -60,   480,   -60}
    }
   ,{{   780,  -580,   240,   780,   240}
    ,{   780,  -580,   240,   780,   240}
    ,{   470,  -890,   -70,   470,   -70}
    ,{ -1150, -1970, -1150, -1860, -1150}
    ,{   470,  -890,   -70,   470,   -70}
    }
   ,{{   490,  -640,   -50,   490,   -50}
    ,{   490,  -880,   -50,   490,   -50}
    ,{   480,  -640,   -60,   480,   -60}
    ,{   490,  -880,   -50,   490,   -50}
    ,{   480,  -640,   -60,   480,   -60}
    }
   ,{{   770,  -890,   770,   470,   770}
    ,{  -970, -1790,  -970, -1680,  -970}
    ,{   470,  -890,   -70,   470,   -70}
    ,{   770, -1300,   770, -1190,   770}
    ,{   470,  -890,   -70,   470,   -70}
    }
   ,{{   490,  -640,   -50,   490,   -50}
    ,{   490,  -880,   -50,   490,   -50}
    ,{   480,  -640,   -60,   480,   -60}
    ,{   490,  -880,   -50,   490,   -50}
    ,{ -1040, -1860, -1040, -1750, -1040}
    }
   }
  ,{{{   690,   690,   470,   690,   190}
    ,{   190,   150,   -60,   150,   190}
    ,{  -140,  -140,  -360,  -140,  -350}
    ,{   690,   690,   470,   690,  -340}
    ,{  -140,  -140,  -360,  -140,  -350}
    }
   ,{{   190,   150,   -60,   150,   190}
    ,{   190,   150,   -60,   150,   190}
    ,{  -150,  -150,  -370,  -150,  -360}
    ,{ -1210, -1230, -1210, -1230, -1440}
    ,{  -150,  -150,  -370,  -150,  -360}
    }
   ,{{  -140,  -140,  -360,  -140,  -340}
    ,{  -140,  -140,  -360,  -140,  -340}
    ,{  -140,  -140,  -360,  -140,  -350}
    ,{  -140,  -140,  -360,  -140,  -340}
    ,{  -140,  -140,  -360,  -140,  -350}
    }
   ,{{   690,   690,   470,   690,  -360}
    ,{ -1030, -1050, -1030, -1050, -1260}
    ,{  -150,  -150,  -370,  -150,  -360}
    ,{   690,   690,   470,   690,  -770}
    ,{  -150,  -150,  -370,  -150,  -360}
    }
   ,{{  -140,  -140,  -360,  -140,  -340}
    ,{  -140,  -140,  -360,  -140,  -340}
    ,{  -140,  -140,  -360,  -140,  -350}
    ,{  -140,  -140,  -360,  -140,  -340}
    ,{ -1120, -1120, -1340, -1120, -1330}
    }
   }
  }
 ,{{{{  1320,  1320,   960,   870,   960}
    ,{   850,   670,   300,   850,   300}
    ,{   720,   540,   170,   720,   170}
    ,{  1320,  1320,   960,   870,   960}
    ,{   590,   410,    40,   590,    40}
    }
   ,{{   850,   670,   300,   850,   300}
    ,{   850,   670,   300,   850,   300}
    ,{   570,   390,    20,   570,    20}
    ,{  -730,  -730, -1100, -1180, -1100}
    ,{   570,   390,    20,   570,    20}
    }
   ,{{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   590,   410,    40,   590,    40}
    }
   ,{{  1320,  1320,   960,   870,   960}
    ,{ -1030, -1030, -1400, -1480, -1400}
    ,{   570,   390,    20,   570,    20}
    ,{  1320,  1320,   960,   870,   960}
    ,{   570,   390,    20,   570,    20}
    }
   ,{{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   280,   100,  -260,   280,  -260}
    ,{   720,   540,   170,   720,   170}
    ,{  -160,  -160,  -760,  -850,  -760}
    }
   }
  ,{{{  1320,  1320,   960,    70,   960}
    ,{   670,   670,   300,   -40,   300}
    ,{   540,   540,   170,    70,   170}
    ,{  1320,  1320,   960,  -170,   960}
    ,{   410,   410,    40,   -60,    40}
    }
   ,{{   670,   670,   300,   -40,   300}
    ,{   670,   670,   300,   -40,   300}
    ,{   390,   390,    20,  -320,    20}
    ,{  -730,  -730, -1100, -1450, -1100}
    ,{   390,   390,    20,  -320,    20}
    }
   ,{{   540,   540,   170,    70,   170}
    ,{   540,   540,   170,  -170,   170}
    ,{   540,   540,   170,    70,   170}
    ,{   540,   540,   170,  -170,   170}
    ,{   410,   410,    40,   -60,    40}
    }
   ,{{  1320,  1320,   960,  -320,   960}
    ,{ -1030, -1030, -1400, -1750, -1400}
    ,{   390,   390,    20,  -320,    20}
    ,{  1320,  1320,   960,  -640,   960}
    ,{   390,   390,    20,  -320,    20}
    }
   ,{{   540,   540,   170,  -170,   170}
    ,{   540,   540,   170,  -170,   170}
    ,{   100,   100,  -260,  -370,  -260}
    ,{   540,   540,   170,  -170,   170}
    ,{  -160,  -160,  -760, -1110,  -760}
    }
   }
  ,{{{   870,   870,   650,   870,   340}
    ,{   220,   220,     0,   220,  -310}
    ,{    90,    90,  -130,    90,  -440}
    ,{   870,   870,   650,   870,   340}
    ,{   -40,   -40,  -260,   -40,  -570}
    }
   ,{{   220,   220,     0,   220,  -310}
    ,{   220,   220,     0,   220,  -310}
    ,{   -60,   -60,  -280,   -60,  -590}
    ,{ -1160, -1180, -1160, -1180, -1470}
    ,{   -60,   -60,  -280,   -60,  -590}
    }
   ,{{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{   -40,   -40,  -260,   -40,  -570}
    }
   ,{{   870,   870,   650,   870,   340}
    ,{ -1460, -1480, -1460, -1480, -1770}
    ,{   -60,   -60,  -280,   -60,  -590}
    ,{   870,   870,   650,   870,   340}
    ,{   -60,   -60,  -280,   -60,  -590}
    }
   ,{{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{  -350,  -350,  -570,  -350,  -880}
    ,{    90,    90,  -130,    90,  -440}
    ,{  -850,  -850, -1070,  -850, -1380}
    }
   }
  ,{{{   960,  -410,   960,   850,   960}
    ,{   850,  -520,   300,   850,   300}
    ,{   720,  -410,   170,   720,   170}
    ,{   960,  -650,   960,   720,   960}
    ,{   590,  -540,    40,   590,    40}
    }
   ,{{   850,  -520,   300,   850,   300}
    ,{   850,  -520,   300,   850,   300}
    ,{   570,  -800,    20,   570,    20}
    ,{ -1100, -1920, -1100, -1810, -1100}
    ,{   570,  -800,    20,   570,    20}
    }
   ,{{   720,  -410,   170,   720,   170}
    ,{   720,  -650,   170,   720,   170}
    ,{   720,  -410,   170,   720,   170}
    ,{   720,  -650,   170,   720,   170}
    ,{   590,  -540,    40,   590,    40}
    }
   ,{{   960,  -800,   960,   570,   960}
    ,{ -1400, -2220, -1400, -2110, -1400}
    ,{   570,  -800,    20,   570,    20}
    ,{   960, -1120,   960, -1000,   960}
    ,{   570,  -800,    20,   570,    20}
    }
   ,{{   720,  -650,   170,   720,   170}
    ,{   720,  -650,   170,   720,   170}
    ,{   280,  -850,  -260,   280,  -260}
    ,{   720,  -650,   170,   720,   170}
    ,{  -760, -1590,  -760, -1470,  -760}
    }
   }
  ,{{{   870,   870,   650,   870,   250}
    ,{   250,   220,     0,   220,   250}
    ,{    90,    90,  -130,    90,  -110}
    ,{   870,   870,   650,   870,  -110}
    ,{   -40,   -40,  -260,   -40,  -240}
    }
   ,{{   250,   220,     0,   220,   250}
    ,{   250,   220,     0,   220,   250}
    ,{   -60,   -60,  -280,   -60,  -260}
    ,{ -1160, -1180, -1160, -1180, -1390}
    ,{   -60,   -60,  -280,   -60,  -260}
    }
   ,{{    90,    90,  -130,    90,  -110}
    ,{    90,    90,  -130,    90,  -110}
    ,{    90,    90,  -130,    90,  -110}
    ,{    90,    90,  -130,    90,  -110}
    ,{   -40,   -40,  -260,   -40,  -240}
    }
   ,{{   870,   870,   650,   870,  -260}
    ,{ -1460, -1480, -1460, -1480, -1690}
    ,{   -60,   -60,  -280,   -60,  -260}
    ,{   870,   870,   650,   870,  -580}
    ,{   -60,   -60,  -280,   -60,  -260}
    }
   ,{{    90,    90,  -130,    90,  -110}
    ,{    90,    90,  -130,    90,  -110}
    ,{  -350,  -350,  -570,  -350,  -550}
    ,{    90,    90,  -130,    90,  -110}
    ,{  -850,  -850, -1070,  -850, -1050}
    }
   }
  }
 ,{{{{  1320,  1320,   960,   870,   960}
    ,{   850,   670,   540,   850,   300}
    ,{   720,   540,   170,   720,   170}
    ,{  1320,  1320,   960,   870,   960}
    ,{   590,   410,    40,   590,    40}
    }
   ,{{   850,   670,   300,   850,   300}
    ,{   850,   670,   300,   850,   300}
    ,{   570,   390,    20,   570,    20}
    ,{  -350,  -350,  -870,  -960,  -870}
    ,{   570,   390,    20,   570,    20}
    }
   ,{{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   590,   410,    40,   590,    40}
    }
   ,{{  1320,  1320,   960,   870,   960}
    ,{   540,  -100,   540, -1050,  -810}
    ,{   570,   390,    20,   570,    20}
    ,{  1320,  1320,   960,   870,   960}
    ,{   570,   390,    20,   570,    20}
    }
   ,{{   720,   540,   170,   720,   170}
    ,{   720,   540,   170,   720,   170}
    ,{   480,   300,   -60,   480,   -60}
    ,{   720,   540,   170,   720,   170}
    ,{  -160,  -160,  -400,  -230,  -760}
    }
   }
  ,{{{  1320,  1320,   960,    70,   960}
    ,{   670,   670,   300,   -40,   300}
    ,{   540,   540,   170,    70,   170}
    ,{  1320,  1320,   960,  -170,   960}
    ,{   410,   410,    40,   -60,    40}
    }
   ,{{   670,   670,   300,   -40,   300}
    ,{   670,   670,   300,   -40,   300}
    ,{   390,   390,    20,  -320,    20}
    ,{  -730,  -730, -1100, -1450,  -870}
    ,{   390,   390,    20,  -320,    20}
    }
   ,{{   540,   540,   170,    70,   170}
    ,{   540,   540,   170,  -170,   170}
    ,{   540,   540,   170,    70,   170}
    ,{   540,   540,   170,  -170,   170}
    ,{   410,   410,    40,   -60,    40}
    }
   ,{{  1320,  1320,   960,  -320,   960}
    ,{    10,  -600,    10, -1320,  -970}
    ,{   390,   390,    20,  -320,    20}
    ,{  1320,  1320,   960,  -640,   960}
    ,{   390,   390,    20,  -320,    20}
    }
   ,{{   540,   540,   170,  -170,   170}
    ,{   540,   540,   170,  -170,   170}
    ,{   300,   300,   -60,  -170,   -60}
    ,{   540,   540,   170,  -170,   170}
    ,{  -160,  -160,  -400, -1110,  -760}
    }
   }
  ,{{{   870,   870,   650,   870,   340}
    ,{   540,   220,   540,   220,  -310}
    ,{    90,    90,  -130,    90,  -440}
    ,{   870,   870,   650,   870,   340}
    ,{   -40,   -40,  -260,   -40,  -570}
    }
   ,{{   220,   220,     0,   220,  -310}
    ,{   220,   220,     0,   220,  -310}
    ,{   -60,   -60,  -280,   -60,  -590}
    ,{  -350,  -350,  -940,  -960, -1250}
    ,{   -60,   -60,  -280,   -60,  -590}
    }
   ,{{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{   -40,   -40,  -260,   -40,  -570}
    }
   ,{{   870,   870,   650,   870,   340}
    ,{   540,  -100,   540, -1050, -1340}
    ,{   -60,   -60,  -280,   -60,  -590}
    ,{   870,   870,   650,   870,   340}
    ,{   -60,   -60,  -280,   -60,  -590}
    }
   ,{{    90,    90,  -130,    90,  -440}
    ,{    90,    90,  -130,    90,  -440}
    ,{  -140,  -140,  -360,  -140,  -670}
    ,{    90,    90,  -130,    90,  -440}
    ,{  -850,  -850, -1070,  -850, -1380}
    }
   }
  ,{{{   960,  -410,   960,   850,   960}
    ,{   850,  -520,   300,   850,   300}
    ,{   720,  -410,   170,   720,   170}
    ,{   960,  -650,   960,   720,   960}
    ,{   590,  -540,    40,   590,    40}
    }
   ,{{   850,  -520,   300,   850,   300}
    ,{   850,  -520,   300,   850,   300}
    ,{   570,  -800,    20,   570,    20}
    ,{  -870, -1920,  -870, -1370,  -870}
    ,{   570,  -800,    20,   570,    20}
    }
   ,{{   720,  -410,   170,   720,   170}
    ,{   720,  -650,   170,   720,   170}
    ,{   720,  -410,   170,   720,   170}
    ,{   720,  -650,   170,   720,   170}
    ,{   590,  -540,    40,   590,    40}
    }
   ,{{   960,  -800,   960,   570,   960}
    ,{  -970, -1790,  -970, -1680,  -970}
    ,{   570,  -800,    20,   570,    20}
    ,{   960, -1120,   960, -1000,   960}
    ,{   570,  -800,    20,   570,    20}
    }
   ,{{   720,  -640,   170,   720,   170}
    ,{   720,  -650,   170,   720,   170}
    ,{   480,  -640,   -60,   480,   -60}
    ,{   720,  -650,   170,   720,   170}
    ,{  -230, -1520,  -760,  -230,  -760}
    }
   }
  ,{{{   870,   870,   650,   870,   250}
    ,{   250,   220,     0,   220,   250}
    ,{    90,    90,  -130,    90,  -110}
    ,{   870,   870,   650,   870,  -110}
    ,{   -40,   -40,  -260,   -40,  -240}
    }
   ,{{   250,   220,     0,   220,   250}
    ,{   250,   220,     0,   220,   250}
    ,{   -60,   -60,  -280,   -60,  -260}
    ,{  -940,  -960,  -940,  -960, -1360}
    ,{   -60,   -60,  -280,   -60,  -260}
    }
   ,{{    90,    90,  -130,    90,   -90}
    ,{    90,    90,  -130,    90,   -90}
    ,{    90,    90,  -130,    90,  -110}
    ,{    90,    90,  -130,    90,  -110}
    ,{   -40,   -40,  -260,   -40,  -240}
    }
   ,{{   870,   870,   650,   870,  -260}
    ,{  -810, -1050, -1030, -1050,  -810}
    ,{   -60,   -60,  -280,   -60,  -260}
    ,{   870,   870,   650,   870,  -580}
    ,{   -60,   -60,  -280,   -60,  -260}
    }
   ,{{    90,    90,  -130,    90,  -110}
    ,{    90,    90,  -130,    90,  -110}
    ,{  -140,  -140,  -360,  -140,  -350}
    ,{    90,    90,  -130,    90,  -110}
    ,{  -850,  -850, -1070,  -850, -1050}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   240,  -780,  -870,   240,  -870}
    ,{   190, -1060, -1060,   190,  -970}
    ,{   240,  -780, -1010,   240, -1010}
    ,{   190,  -870,  -870,   190,  -870}
    ,{   130,  -890, -1120,   130, -1120}
    }
   ,{{    40, -1210, -1180,    40,  -970}
    ,{    40, -1210, -1210,    40,  -970}
    ,{  -270, -1520, -1520,  -270, -1520}
    ,{ -1180, -1420, -1180, -1250, -1180}
    ,{  -270, -1520, -1520,  -270, -1520}
    }
   ,{{   190,  -840, -1060,   190, -1060}
    ,{   190, -1060, -1060,   190, -1060}
    ,{   180,  -840, -1070,   180, -1070}
    ,{   190, -1060, -1060,   190, -1060}
    ,{   130,  -890, -1120,   130, -1120}
    }
   ,{{  -270,  -870,  -870,  -270,  -870}
    ,{ -1470, -1710, -1470, -1530, -1470}
    ,{  -270, -1520, -1520,  -270, -1520}
    ,{  -870,  -870,  -870,  -870,  -870}
    ,{  -270, -1520, -1520,  -270, -1520}
    }
   ,{{   240,  -780, -1010,   240, -1010}
    ,{   190, -1060, -1060,   190, -1060}
    ,{   240,  -780, -1010,   240, -1010}
    ,{   190, -1060, -1060,   190, -1060}
    ,{ -1680, -1790, -1850, -1680, -1850}
    }
   }
  ,{{{  -590, -1050,  -870,  -590,  -870}
    ,{  -890, -1240, -1060,  -890, -1060}
    ,{  -590, -1190, -1010,  -590, -1010}
    ,{  -870, -1050,  -870,  -890,  -870}
    ,{  -700, -1300, -1120,  -700, -1120}
    }
   ,{{ -1030, -1370, -1210, -1030, -1210}
    ,{ -1030, -1370, -1210, -1030, -1210}
    ,{ -1340, -1700, -1520, -1340, -1520}
    ,{ -1250, -1600, -1420, -1250, -1420}
    ,{ -1340, -1700, -1520, -1340, -1520}
    }
   ,{{  -650, -1240, -1060,  -650, -1060}
    ,{  -890, -1240, -1060,  -890, -1060}
    ,{  -650, -1250, -1070,  -650, -1070}
    ,{  -890, -1240, -1060,  -890, -1060}
    ,{  -700, -1300, -1120,  -700, -1120}
    }
   ,{{  -870, -1050,  -870, -1340,  -870}
    ,{ -1530, -1890, -1710, -1530, -1710}
    ,{ -1340, -1700, -1520, -1340, -1520}
    ,{  -870, -1050,  -870, -1940,  -870}
    ,{ -1340, -1700, -1520, -1340, -1520}
    }
   ,{{  -590, -1190, -1010,  -590, -1010}
    ,{  -890, -1240, -1060,  -890, -1060}
    ,{  -590, -1190, -1010,  -590, -1010}
    ,{  -890, -1240, -1060,  -890, -1060}
    ,{ -1680, -1790, -1850, -1680, -1850}
    }
   }
  ,{{{  -870,  -870,  -870,  -870,  -870}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1010, -1010, -1010, -1010, -1010}
    ,{  -870,  -870,  -870,  -870,  -870}
    ,{ -1120, -1120, -1120, -1120, -1120}
    }
   ,{{ -1180, -1210, -1180, -1210, -1180}
    ,{ -1210, -1210, -1210, -1210, -1210}
    ,{ -1520, -1520, -1520, -1520, -1520}
    ,{ -1180, -1420, -1180, -1420, -1180}
    ,{ -1520, -1520, -1520, -1520, -1520}
    }
   ,{{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1070, -1070, -1070, -1070, -1070}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1120, -1120, -1120, -1120, -1120}
    }
   ,{{  -870,  -870,  -870,  -870,  -870}
    ,{ -1470, -1710, -1470, -1710, -1470}
    ,{ -1520, -1520, -1520, -1520, -1520}
    ,{  -870,  -870,  -870,  -870,  -870}
    ,{ -1520, -1520, -1520, -1520, -1520}
    }
   ,{{ -1010, -1010, -1010, -1010, -1010}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1010, -1010, -1010, -1010, -1010}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1850, -1850, -1850, -1850, -1850}
    }
   }
  ,{{{   240,  -780,  -870,   240,  -870}
    ,{   190, -1080, -1060,   190, -1060}
    ,{   240,  -780, -1010,   240, -1010}
    ,{   190, -1080,  -870,   190,  -870}
    ,{   130,  -890, -1120,   130, -1120}
    }
   ,{{    40, -1220, -1210,    40, -1210}
    ,{    40, -1220, -1210,    40, -1210}
    ,{  -270, -1530, -1520,  -270, -1520}
    ,{ -1420, -1440, -1420, -1420, -1420}
    ,{  -270, -1530, -1520,  -270, -1520}
    }
   ,{{   190,  -840, -1060,   190, -1060}
    ,{   190, -1080, -1060,   190, -1060}
    ,{   180,  -840, -1070,   180, -1070}
    ,{   190, -1080, -1060,   190, -1060}
    ,{   130,  -890, -1120,   130, -1120}
    }
   ,{{  -270, -1530,  -870,  -270,  -870}
    ,{ -1710, -1720, -1710, -1710, -1710}
    ,{  -270, -1530, -1520,  -270, -1520}
    ,{  -870, -2130,  -870, -2120,  -870}
    ,{  -270, -1530, -1520,  -270, -1520}
    }
   ,{{   240,  -780, -1010,   240, -1010}
    ,{   190, -1080, -1060,   190, -1060}
    ,{   240,  -780, -1010,   240, -1010}
    ,{   190, -1080, -1060,   190, -1060}
    ,{ -1850, -1870, -1850, -1850, -1850}
    }
   }
  ,{{{  -870,  -870,  -870,  -870,  -970}
    ,{  -970, -1060, -1060, -1060,  -970}
    ,{ -1010, -1010, -1010, -1010, -1010}
    ,{  -870,  -870,  -870,  -870, -1060}
    ,{ -1120, -1120, -1120, -1120, -1120}
    }
   ,{{  -970, -1210, -1180, -1210,  -970}
    ,{  -970, -1210, -1210, -1210,  -970}
    ,{ -1520, -1520, -1520, -1520, -1520}
    ,{ -1180, -1420, -1180, -1420, -1420}
    ,{ -1520, -1520, -1520, -1520, -1520}
    }
   ,{{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1070, -1070, -1070, -1070, -1070}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1120, -1120, -1120, -1120, -1120}
    }
   ,{{  -870,  -870,  -870,  -870, -1520}
    ,{ -1470, -1710, -1470, -1710, -1710}
    ,{ -1520, -1520, -1520, -1520, -1520}
    ,{  -870,  -870,  -870,  -870, -2120}
    ,{ -1520, -1520, -1520, -1520, -1520}
    }
   ,{{ -1010, -1010, -1010, -1010, -1010}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1010, -1010, -1010, -1010, -1010}
    ,{ -1060, -1060, -1060, -1060, -1060}
    ,{ -1850, -1850, -1850, -1850, -1850}
    }
   }
  }
 ,{{{{   210,  -870,  -870,   210,  -800}
    ,{   210, -1040, -1040,   210,  -800}
    ,{  -240, -1490, -1490,  -240, -1490}
    ,{  -160,  -870,  -870,  -160,  -870}
    ,{  -240, -1490, -1490,  -240, -1490}
    }
   ,{{   210, -1040, -1040,   210,  -800}
    ,{   210, -1040, -1040,   210,  -800}
    ,{  -240, -1490, -1490,  -240, -1490}
    ,{ -1990, -2230, -1990, -2060, -1990}
    ,{  -240, -1490, -1490,  -240, -1490}
    }
   ,{{  -160, -1410, -1410,  -160, -1410}
    ,{  -160, -1410, -1410,  -160, -1410}
    ,{  -460, -1490, -1710,  -460, -1710}
    ,{  -160, -1410, -1410,  -160, -1410}
    ,{  -460, -1490, -1710,  -460, -1710}
    }
   ,{{  -240,  -870,  -870,  -240,  -870}
    ,{ -1520, -1760, -1520, -1580, -1520}
    ,{  -240, -1490, -1490,  -240, -1490}
    ,{  -870,  -870,  -870,  -870,  -870}
    ,{  -240, -1490, -1490,  -240, -1490}
    }
   ,{{  -160, -1410, -1410,  -160, -1410}
    ,{  -160, -1410, -1410,  -160, -1410}
    ,{  -770, -1800, -2020,  -770, -2020}
    ,{  -160, -1410, -1410,  -160, -1410}
    ,{ -1520, -1640, -1700, -1520, -1700}
    }
   }
  ,{{{  -870, -1050,  -870,  -870,  -870}
    ,{  -870, -1220, -1040,  -870, -1040}
    ,{ -1300, -1670, -1490, -1300, -1490}
    ,{  -870, -1050,  -870, -1230,  -870}
    ,{ -1300, -1640, -1490, -1300, -1490}
    }
   ,{{  -870, -1220, -1040,  -870, -1040}
    ,{  -870, -1220, -1040,  -870, -1040}
    ,{ -1320, -1670, -1490, -1320, -1490}
    ,{ -2060, -2410, -2230, -2060, -2230}
    ,{ -1320, -1670, -1490, -1320, -1490}
    }
   ,{{ -1230, -1590, -1410, -1230, -1410}
    ,{ -1230, -1590, -1410, -1230, -1410}
    ,{ -1300, -1890, -1710, -1300, -1710}
    ,{ -1230, -1590, -1410, -1230, -1410}
    ,{ -1300, -1890, -1710, -1300, -1710}
    }
   ,{{  -870, -1050,  -870, -1320,  -870}
    ,{ -1580, -1940, -1760, -1580, -1760}
    ,{ -1320, -1670, -1490, -1320, -1490}
    ,{  -870, -1050,  -870, -1940,  -870}
    ,{ -1320, -1670, -1490, -1320, -1490}
    }
   ,{{ -1230, -1590, -1410, -1230, -1410}
    ,{ -1230, -1590, -1410, -1230, -1410}
    ,{ -1610, -2200, -2020, -1610, -2020}
    ,{ -1230, -1590, -1410, -1230, -1410}
    ,{ -1520, -1640, -1700, -1520, -1700}
    }
   }
  ,{{{  -870,  -870,  -870,  -870,  -870}
    ,{ -1040, -1040, -1040, -1040, -1040}
    ,{ -1490, -1490, -1490, -1490, -1490}
    ,{  -870,  -870,  -870,  -870,  -870}
    ,{ -1490, -1490, -1490, -1490, -1490}
    }
   ,{{ -1040, -1040, -1040, -1040, -1040}
    ,{ -1040, -1040, -1040, -1040, -1040}
    ,{ -1490, -1490, -1490, -1490, -1490}
    ,{ -1990, -2230, -1990, -2230, -1990}
    ,{ -1490, -1490, -1490, -1490, -1490}
    }
   ,{{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1710, -1710, -1710, -1710, -1710}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1710, -1710, -1710, -1710, -1710}
    }
   ,{{  -870,  -870,  -870,  -870,  -870}
    ,{ -1520, -1760, -1520, -1760, -1520}
    ,{ -1490, -1490, -1490, -1490, -1490}
    ,{  -870,  -870,  -870,  -870,  -870}
    ,{ -1490, -1490, -1490, -1490, -1490}
    }
   ,{{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -2020, -2020, -2020, -2020, -2020}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1700, -1700, -1700, -1700, -1700}
    }
   }
  ,{{{   210, -1060,  -870,   210,  -870}
    ,{   210, -1060, -1040,   210, -1040}
    ,{  -240, -1490, -1490,  -240, -1490}
    ,{  -160, -1420,  -870,  -160,  -870}
    ,{  -240, -1490, -1490,  -240, -1490}
    }
   ,{{   210, -1060, -1040,   210, -1040}
    ,{   210, -1060, -1040,   210, -1040}
    ,{  -240, -1510, -1490,  -240, -1490}
    ,{ -2230, -2250, -2230, -2230, -2230}
    ,{  -240, -1510, -1490,  -240, -1490}
    }
   ,{{  -160, -1420, -1410,  -160, -1410}
    ,{  -160, -1420, -1410,  -160, -1410}
    ,{  -460, -1490, -1710,  -460, -1710}
    ,{  -160, -1420, -1410,  -160, -1410}
    ,{  -460, -1490, -1710,  -460, -1710}
    }
   ,{{  -240, -1510,  -870,  -240,  -870}
    ,{ -1760, -1770, -1760, -1760, -1760}
    ,{  -240, -1510, -1490,  -240, -1490}
    ,{  -870, -2130,  -870, -2120,  -870}
    ,{  -240, -1510, -1490,  -240, -1490}
    }
   ,{{  -160, -1420, -1410,  -160, -1410}
    ,{  -160, -1420, -1410,  -160, -1410}
    ,{  -770, -1800, -2020,  -770, -2020}
    ,{  -160, -1420, -1410,  -160, -1410}
    ,{ -1700, -1710, -1700, -1700, -1700}
    }
   }
  ,{{{  -800,  -870,  -870,  -870,  -800}
    ,{  -800, -1040, -1040, -1040,  -800}
    ,{ -1490, -1490, -1490, -1490, -1490}
    ,{  -870,  -870,  -870,  -870, -1410}
    ,{ -1490, -1490, -1490, -1490, -1490}
    }
   ,{{  -800, -1040, -1040, -1040,  -800}
    ,{  -800, -1040, -1040, -1040,  -800}
    ,{ -1490, -1490, -1490, -1490, -1490}
    ,{ -1990, -2230, -1990, -2230, -2230}
    ,{ -1490, -1490, -1490, -1490, -1490}
    }
   ,{{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1710, -1710, -1710, -1710, -1710}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1710, -1710, -1710, -1710, -1710}
    }
   ,{{  -870,  -870,  -870,  -870, -1490}
    ,{ -1520, -1760, -1520, -1760, -1760}
    ,{ -1490, -1490, -1490, -1490, -1490}
    ,{  -870,  -870,  -870,  -870, -2120}
    ,{ -1490, -1490, -1490, -1490, -1490}
    }
   ,{{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -2020, -2020, -2020, -2020, -2020}
    ,{ -1410, -1410, -1410, -1410, -1410}
    ,{ -1700, -1700, -1700, -1700, -1700}
    }
   }
  }
 ,{{{{  -710,  -710,  -710,  -710,  -710}
    ,{  -710, -1780, -1540,  -710, -1540}
    ,{  -710, -1730, -1960,  -710, -1960}
    ,{  -710,  -710,  -710,  -710,  -710}
    ,{  -710, -1730, -1960,  -710, -1960}
    }
   ,{{  -710, -1960, -1730,  -710, -1730}
    ,{  -890, -2140, -2140,  -890, -1900}
    ,{  -710, -1960, -1960,  -710, -1960}
    ,{ -1730, -1970, -1730, -1800, -1730}
    ,{  -710, -1960, -1960,  -710, -1960}
    }
   ,{{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1960, -1960,  -710, -1960}
    ,{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1960, -1960,  -710, -1960}
    ,{  -710, -1730, -1960,  -710, -1960}
    }
   ,{{  -710,  -710,  -710,  -710,  -710}
    ,{ -1540, -1780, -1540, -1610, -1540}
    ,{  -710, -1960, -1960,  -710, -1960}
    ,{  -710,  -710,  -710,  -710,  -710}
    ,{  -710, -1960, -1960,  -710, -1960}
    }
   ,{{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1960, -1960,  -710, -1960}
    ,{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1960, -1960,  -710, -1960}
    ,{ -1780, -1900, -1960, -1780, -1960}
    }
   }
  ,{{{  -710,  -890,  -710, -1540,  -710}
    ,{ -1610, -1960, -1780, -1610, -1780}
    ,{ -1540, -2140, -1960, -1540, -1960}
    ,{  -710,  -890,  -710, -1780,  -710}
    ,{ -1540, -1900, -1960, -1540, -1960}
    }
   ,{{ -1780, -2140, -1960, -1780, -1960}
    ,{ -1960, -2320, -2140, -1960, -2140}
    ,{ -1780, -2140, -1960, -1780, -1960}
    ,{ -1800, -2150, -1970, -1800, -1970}
    ,{ -1780, -2140, -1960, -1780, -1960}
    }
   ,{{ -1540, -2140, -1960, -1540, -1960}
    ,{ -1780, -2140, -1960, -1780, -1960}
    ,{ -1540, -2140, -1960, -1540, -1960}
    ,{ -1780, -2140, -1960, -1780, -1960}
    ,{ -1540, -2140, -1960, -1540, -1960}
    }
   ,{{  -710,  -890,  -710, -1610,  -710}
    ,{ -1610, -1960, -1780, -1610, -1780}
    ,{ -1780, -2140, -1960, -1780, -1960}
    ,{  -710,  -890,  -710, -1780,  -710}
    ,{ -1780, -2140, -1960, -1780, -1960}
    }
   ,{{ -1540, -1900, -1960, -1540, -1960}
    ,{ -1780, -2140, -1960, -1780, -1960}
    ,{ -1540, -2140, -1960, -1540, -1960}
    ,{ -1780, -2140, -1960, -1780, -1960}
    ,{ -1780, -1900, -1960, -1780, -1960}
    }
   }
  ,{{{  -710,  -710,  -710,  -710,  -710}
    ,{ -1540, -1780, -1540, -1780, -1540}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{  -710,  -710,  -710,  -710,  -710}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{ -1730, -1960, -1730, -1960, -1730}
    ,{ -2140, -2140, -2140, -2140, -2140}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1730, -1970, -1730, -1970, -1730}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{  -710,  -710,  -710,  -710,  -710}
    ,{ -1540, -1780, -1540, -1780, -1540}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{  -710,  -710,  -710,  -710,  -710}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   }
  ,{{{  -710, -1730,  -710,  -710,  -710}
    ,{  -710, -1800, -1780,  -710, -1780}
    ,{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1970,  -710,  -710,  -710}
    ,{  -710, -1730, -1960,  -710, -1960}
    }
   ,{{  -710, -1970, -1960,  -710, -1960}
    ,{  -890, -2150, -2140,  -890, -2140}
    ,{  -710, -1970, -1960,  -710, -1960}
    ,{ -1970, -1990, -1970, -1970, -1970}
    ,{  -710, -1970, -1960,  -710, -1960}
    }
   ,{{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1970, -1960,  -710, -1960}
    ,{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1970, -1960,  -710, -1960}
    ,{  -710, -1730, -1960,  -710, -1960}
    }
   ,{{  -710, -1800,  -710,  -710,  -710}
    ,{ -1780, -1800, -1780, -1780, -1780}
    ,{  -710, -1970, -1960,  -710, -1960}
    ,{  -710, -1970,  -710, -1960,  -710}
    ,{  -710, -1970, -1960,  -710, -1960}
    }
   ,{{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1970, -1960,  -710, -1960}
    ,{  -710, -1730, -1960,  -710, -1960}
    ,{  -710, -1970, -1960,  -710, -1960}
    ,{ -1960, -1970, -1960, -1960, -1960}
    }
   }
  ,{{{  -710,  -710,  -710,  -710, -1780}
    ,{ -1540, -1780, -1540, -1780, -1780}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{  -710,  -710,  -710,  -710, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{ -1730, -1960, -1730, -1960, -1900}
    ,{ -1900, -2140, -2140, -2140, -1900}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1730, -1970, -1730, -1970, -1970}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{  -710,  -710,  -710,  -710, -1780}
    ,{ -1540, -1780, -1540, -1780, -1780}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{  -710,  -710,  -710,  -710, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   ,{{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    ,{ -1960, -1960, -1960, -1960, -1960}
    }
   }
  }
 ,{{{{   360,   -70,  -150,   360,  -150}
    ,{   360,   -70,  -890,   360,  -650}
    ,{  -150, -1180, -1400,  -150, -1400}
    ,{  -150,  -150,  -150,  -150,  -150}
    ,{  -150, -1180, -1400,  -150, -1400}
    }
   ,{{   360,   -70,  -890,   360,  -650}
    ,{   360,   -70,  -890,   360,  -650}
    ,{  -150, -1400, -1400,  -150, -1400}
    ,{ -1500, -1600, -1500, -1570, -1500}
    ,{  -150, -1400, -1400,  -150, -1400}
    }
   ,{{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1400, -1400,  -150, -1400}
    ,{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1400, -1400,  -150, -1400}
    ,{  -150, -1180, -1400,  -150, -1400}
    }
   ,{{  -150,  -150,  -150,  -150,  -150}
    ,{ -1670, -1910, -1670, -1740, -1670}
    ,{  -150, -1400, -1400,  -150, -1400}
    ,{  -150,  -150,  -150,  -150,  -150}
    ,{  -150, -1400, -1400,  -150, -1400}
    }
   ,{{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1400, -1400,  -150, -1400}
    ,{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1400, -1400,  -150, -1400}
    ,{ -1230, -1340, -1400, -1230, -1400}
    }
   }
  ,{{{   -30,   -70,  -150,   -30,  -150}
    ,{   -30,   -70,  -890,   -30,  -890}
    ,{  -990, -1580, -1400,  -990, -1400}
    ,{  -150,  -330,  -150, -1230,  -150}
    ,{  -990, -1340, -1400,  -990, -1400}
    }
   ,{{   -30,   -70,  -890,   -30,  -890}
    ,{   -30,   -70,  -890,   -30,  -890}
    ,{ -1230, -1580, -1400, -1230, -1400}
    ,{ -1570, -1600, -1740, -1570, -1740}
    ,{ -1230, -1580, -1400, -1230, -1400}
    }
   ,{{  -990, -1580, -1400,  -990, -1400}
    ,{ -1230, -1580, -1400, -1230, -1400}
    ,{  -990, -1580, -1400,  -990, -1400}
    ,{ -1230, -1580, -1400, -1230, -1400}
    ,{  -990, -1580, -1400,  -990, -1400}
    }
   ,{{  -150,  -330,  -150, -1230,  -150}
    ,{ -1740, -2090, -1910, -1740, -1910}
    ,{ -1230, -1580, -1400, -1230, -1400}
    ,{  -150,  -330,  -150, -1230,  -150}
    ,{ -1230, -1580, -1400, -1230, -1400}
    }
   ,{{  -990, -1340, -1400,  -990, -1400}
    ,{ -1230, -1580, -1400, -1230, -1400}
    ,{  -990, -1580, -1400,  -990, -1400}
    ,{ -1230, -1580, -1400, -1230, -1400}
    ,{ -1230, -1340, -1400, -1230, -1400}
    }
   }
  ,{{{  -150,  -150,  -150,  -150,  -150}
    ,{  -890,  -890,  -890,  -890,  -890}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150,  -150}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{  -890,  -890,  -890,  -890,  -890}
    ,{  -890,  -890,  -890,  -890,  -890}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1500, -1740, -1500, -1740, -1500}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{  -150,  -150,  -150,  -150,  -150}
    ,{ -1670, -1910, -1670, -1910, -1670}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150,  -150}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   }
  ,{{{   360,  -910,  -150,   360,  -150}
    ,{   360,  -910,  -890,   360,  -890}
    ,{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1420,  -150,  -150,  -150}
    ,{  -150, -1180, -1400,  -150, -1400}
    }
   ,{{   360,  -910,  -890,   360,  -890}
    ,{   360,  -910,  -890,   360,  -890}
    ,{  -150, -1420, -1400,  -150, -1400}
    ,{ -1740, -3040, -1740, -1740, -1740}
    ,{  -150, -1420, -1400,  -150, -1400}
    }
   ,{{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1420, -1400,  -150, -1400}
    ,{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1420, -1400,  -150, -1400}
    ,{  -150, -1180, -1400,  -150, -1400}
    }
   ,{{  -150, -1420,  -150,  -150,  -150}
    ,{ -1910, -1930, -1910, -1910, -1910}
    ,{  -150, -1420, -1400,  -150, -1400}
    ,{  -150, -1420,  -150, -1400,  -150}
    ,{  -150, -1420, -1400,  -150, -1400}
    }
   ,{{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1420, -1400,  -150, -1400}
    ,{  -150, -1180, -1400,  -150, -1400}
    ,{  -150, -1420, -1400,  -150, -1400}
    ,{ -1400, -1420, -1400, -1400, -1400}
    }
   }
  ,{{{  -150,  -150,  -150,  -150,  -650}
    ,{  -650,  -890,  -890,  -890,  -650}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{  -650,  -890,  -890,  -890,  -650}
    ,{  -650,  -890,  -890,  -890,  -650}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1500, -1740, -1500, -1740, -1740}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{  -150,  -150,  -150,  -150, -1400}
    ,{ -1670, -1910, -1670, -1910, -1910}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   }
  }
 ,{{{{   940,   220,   220,   940,   220}
    ,{   940,  -310,  -310,   940,   -70}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   650,   220,   220,   650,   220}
    ,{   640,  -380,  -610,   640,  -610}
    }
   ,{{   940,  -310,  -310,   940,   -70}
    ,{   940,  -310,  -310,   940,   -70}
    ,{   630,  -620,  -620,   630,  -620}
    ,{ -1460, -1700, -1460, -1520, -1460}
    ,{   630,  -620,  -620,   630,  -620}
    }
   ,{{   650,  -380,  -600,   650,  -600}
    ,{   650,  -600,  -600,   650,  -600}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   650,  -600,  -600,   650,  -600}
    ,{   640,  -380,  -610,   640,  -610}
    }
   ,{{   630,   220,   220,   630,   220}
    ,{ -1280, -1520, -1280, -1340, -1280}
    ,{   630,  -620,  -620,   630,  -620}
    ,{   220,   220,   220,   220,   220}
    ,{   630,  -620,  -620,   630,  -620}
    }
   ,{{   650,  -380,  -600,   650,  -600}
    ,{   650,  -600,  -600,   650,  -600}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   650,  -600,  -600,   650,  -600}
    ,{ -1410, -1530, -1590, -1410, -1590}
    }
   }
  ,{{{   220,    40,   220,  -130,   220}
    ,{  -130,  -490,  -310,  -130,  -310}
    ,{  -190,  -790,  -610,  -190,  -610}
    ,{   220,    40,   220,  -430,   220}
    ,{  -190,  -790,  -610,  -190,  -610}
    }
   ,{{  -130,  -490,  -310,  -130,  -310}
    ,{  -130,  -490,  -310,  -130,  -310}
    ,{  -440,  -800,  -620,  -440,  -620}
    ,{ -1520, -1880, -1700, -1520, -1700}
    ,{  -440,  -800,  -620,  -440,  -620}
    }
   ,{{  -190,  -780,  -600,  -190,  -600}
    ,{  -430,  -780,  -600,  -430,  -600}
    ,{  -190,  -790,  -610,  -190,  -610}
    ,{  -430,  -780,  -600,  -430,  -600}
    ,{  -190,  -790,  -610,  -190,  -610}
    }
   ,{{   220,    40,   220,  -440,   220}
    ,{ -1340, -1700, -1520, -1340, -1520}
    ,{  -440,  -800,  -620,  -440,  -620}
    ,{   220,    40,   220,  -850,   220}
    ,{  -440,  -800,  -620,  -440,  -620}
    }
   ,{{  -190,  -780,  -600,  -190,  -600}
    ,{  -430,  -780,  -600,  -430,  -600}
    ,{  -190,  -790,  -610,  -190,  -610}
    ,{  -430,  -780,  -600,  -430,  -600}
    ,{ -1410, -1530, -1590, -1410, -1590}
    }
   }
  ,{{{   220,   220,   220,   220,   220}
    ,{  -310,  -310,  -310,  -310,  -310}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{   220,   220,   220,   220,   220}
    ,{  -610,  -610,  -610,  -610,  -610}
    }
   ,{{  -310,  -310,  -310,  -310,  -310}
    ,{  -310,  -310,  -310,  -310,  -310}
    ,{  -620,  -620,  -620,  -620,  -620}
    ,{ -1460, -1700, -1460, -1700, -1460}
    ,{  -620,  -620,  -620,  -620,  -620}
    }
   ,{{  -600,  -600,  -600,  -600,  -600}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{  -610,  -610,  -610,  -610,  -610}
    }
   ,{{   220,   220,   220,   220,   220}
    ,{ -1280, -1520, -1280, -1520, -1280}
    ,{  -620,  -620,  -620,  -620,  -620}
    ,{   220,   220,   220,   220,   220}
    ,{  -620,  -620,  -620,  -620,  -620}
    }
   ,{{  -600,  -600,  -600,  -600,  -600}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{ -1590, -1590, -1590, -1590, -1590}
    }
   }
  ,{{{   940,  -320,   220,   940,   220}
    ,{   940,  -320,  -310,   940,  -310}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   650,  -620,   220,   650,   220}
    ,{   640,  -380,  -610,   640,  -610}
    }
   ,{{   940,  -320,  -310,   940,  -310}
    ,{   940,  -320,  -310,   940,  -310}
    ,{   630,  -630,  -620,   630,  -620}
    ,{ -1700, -1710, -1700, -1700, -1700}
    ,{   630,  -630,  -620,   630,  -620}
    }
   ,{{   650,  -380,  -600,   650,  -600}
    ,{   650,  -620,  -600,   650,  -600}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   650,  -620,  -600,   650,  -600}
    ,{   640,  -380,  -610,   640,  -610}
    }
   ,{{   630,  -630,   220,   630,   220}
    ,{ -1520, -1530, -1520, -1520, -1520}
    ,{   630,  -630,  -620,   630,  -620}
    ,{   220, -1040,   220, -1030,   220}
    ,{   630,  -630,  -620,   630,  -620}
    }
   ,{{   650,  -380,  -600,   650,  -600}
    ,{   650,  -620,  -600,   650,  -600}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   650,  -620,  -600,   650,  -600}
    ,{ -1590, -1600, -1590, -1590, -1590}
    }
   }
  ,{{{   220,   220,   220,   220,   -70}
    ,{   -70,  -310,  -310,  -310,   -70}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{   220,   220,   220,   220,  -600}
    ,{  -610,  -610,  -610,  -610,  -610}
    }
   ,{{   -70,  -310,  -310,  -310,   -70}
    ,{   -70,  -310,  -310,  -310,   -70}
    ,{  -620,  -620,  -620,  -620,  -620}
    ,{ -1460, -1700, -1460, -1700, -1700}
    ,{  -620,  -620,  -620,  -620,  -620}
    }
   ,{{  -600,  -600,  -600,  -600,  -600}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{  -610,  -610,  -610,  -610,  -610}
    }
   ,{{   220,   220,   220,   220,  -620}
    ,{ -1280, -1520, -1280, -1520, -1520}
    ,{  -620,  -620,  -620,  -620,  -620}
    ,{   220,   220,   220,   220, -1030}
    ,{  -620,  -620,  -620,  -620,  -620}
    }
   ,{{  -600,  -600,  -600,  -600,  -600}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{  -600,  -600,  -600,  -600,  -600}
    ,{ -1590, -1590, -1590, -1590, -1590}
    }
   }
  }
 ,{{{{  1010,   410,   410,  1010,   410}
    ,{  1010,  -240,  -240,  1010,     0}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,   410,   410,   880,   410}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{  1010,  -240,  -240,  1010,     0}
    ,{  1010,  -240,  -240,  1010,     0}
    ,{   730,  -520,  -520,   730,  -520}
    ,{ -1410, -1650, -1410, -1470, -1410}
    ,{   730,  -520,  -520,   730,  -520}
    }
   ,{{   880,  -150,  -370,   880,  -370}
    ,{   880,  -370,  -370,   880,  -370}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,  -370,  -370,   880,  -370}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{   730,   410,   410,   730,   410}
    ,{ -1710, -1950, -1710, -1770, -1710}
    ,{   730,  -520,  -520,   730,  -520}
    ,{   410,   410,   410,   410,   410}
    ,{   730,  -520,  -520,   730,  -520}
    }
   ,{{   880,  -370,  -370,   880,  -370}
    ,{   880,  -370,  -370,   880,  -370}
    ,{   440,  -590,  -810,   440,  -810}
    ,{   880,  -370,  -370,   880,  -370}
    ,{ -1140, -1250, -1310, -1140, -1310}
    }
   }
  ,{{{   410,   230,   410,    40,   410}
    ,{   -70,  -420,  -240,   -70,  -240}
    ,{    40,  -550,  -370,    40,  -370}
    ,{   410,   230,   410,  -200,   410}
    ,{   -90,  -680,  -500,   -90,  -500}
    }
   ,{{   -70,  -420,  -240,   -70,  -240}
    ,{   -70,  -420,  -240,   -70,  -240}
    ,{  -350,  -700,  -520,  -350,  -520}
    ,{ -1470, -1830, -1650, -1470, -1650}
    ,{  -350,  -700,  -520,  -350,  -520}
    }
   ,{{    40,  -550,  -370,    40,  -370}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{    40,  -550,  -370,    40,  -370}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{   -90,  -680,  -500,   -90,  -500}
    }
   ,{{   410,   230,   410,  -350,   410}
    ,{ -1770, -2130, -1950, -1770, -1950}
    ,{  -350,  -700,  -520,  -350,  -520}
    ,{   410,   230,   410,  -670,   410}
    ,{  -350,  -700,  -520,  -350,  -520}
    }
   ,{{  -200,  -550,  -370,  -200,  -370}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{  -400,  -990,  -810,  -400,  -810}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{ -1140, -1250, -1310, -1140, -1310}
    }
   }
  ,{{{   410,   410,   410,   410,   410}
    ,{  -240,  -240,  -240,  -240,  -240}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{   410,   410,   410,   410,   410}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{  -240,  -240,  -240,  -240,  -240}
    ,{  -240,  -240,  -240,  -240,  -240}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{ -1410, -1650, -1410, -1650, -1410}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{   410,   410,   410,   410,   410}
    ,{ -1710, -1950, -1710, -1950, -1710}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{   410,   410,   410,   410,   410}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -810,  -810,  -810,  -810,  -810}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{ -1310, -1310, -1310, -1310, -1310}
    }
   }
  ,{{{  1010,  -150,   410,  1010,   410}
    ,{  1010,  -260,  -240,  1010,  -240}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,  -390,   410,   880,   410}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{  1010,  -260,  -240,  1010,  -240}
    ,{  1010,  -260,  -240,  1010,  -240}
    ,{   730,  -540,  -520,   730,  -520}
    ,{ -1650, -1660, -1650, -1650, -1650}
    ,{   730,  -540,  -520,   730,  -520}
    }
   ,{{   880,  -150,  -370,   880,  -370}
    ,{   880,  -390,  -370,   880,  -370}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,  -390,  -370,   880,  -370}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{   730,  -540,   410,   730,   410}
    ,{ -1950, -1960, -1950, -1950, -1950}
    ,{   730,  -540,  -520,   730,  -520}
    ,{   410,  -860,   410,  -840,   410}
    ,{   730,  -540,  -520,   730,  -520}
    }
   ,{{   880,  -390,  -370,   880,  -370}
    ,{   880,  -390,  -370,   880,  -370}
    ,{   440,  -590,  -810,   440,  -810}
    ,{   880,  -390,  -370,   880,  -370}
    ,{ -1310, -1330, -1310, -1310, -1310}
    }
   }
  ,{{{   410,   410,   410,   410,     0}
    ,{     0,  -240,  -240,  -240,     0}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{   410,   410,   410,   410,  -370}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{     0,  -240,  -240,  -240,     0}
    ,{     0,  -240,  -240,  -240,     0}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{ -1410, -1650, -1410, -1650, -1650}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{   410,   410,   410,   410,  -520}
    ,{ -1710, -1950, -1710, -1950, -1950}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{   410,   410,   410,   410,  -840}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -810,  -810,  -810,  -810,  -810}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{ -1310, -1310, -1310, -1310, -1310}
    }
   }
  }
 ,{{{{  1010,   410,   410,  1010,   410}
    ,{  1010,   -70,  -240,  1010,     0}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,   410,   410,   880,   410}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{  1010,   -70,  -240,  1010,     0}
    ,{  1010,   -70,  -240,  1010,     0}
    ,{   730,  -520,  -520,   730,  -520}
    ,{ -1180, -1420, -1180, -1250, -1180}
    ,{   730,  -520,  -520,   730,  -520}
    }
   ,{{   880,  -150,  -370,   880,  -370}
    ,{   880,  -370,  -370,   880,  -370}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,  -370,  -370,   880,  -370}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{   730,   410,   410,   730,   410}
    ,{ -1280, -1520, -1280, -1340, -1280}
    ,{   730,  -520,  -520,   730,  -520}
    ,{   410,   410,   410,   410,   410}
    ,{   730,  -520,  -520,   730,  -520}
    }
   ,{{   880,  -370,  -370,   880,  -370}
    ,{   880,  -370,  -370,   880,  -370}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   880,  -370,  -370,   880,  -370}
    ,{ -1140, -1250, -1310, -1140, -1310}
    }
   }
  ,{{{   410,   230,   410,    40,   410}
    ,{   -30,   -70,  -240,   -30,  -240}
    ,{    40,  -550,  -370,    40,  -370}
    ,{   410,   230,   410,  -200,   410}
    ,{   -90,  -680,  -500,   -90,  -500}
    }
   ,{{   -30,   -70,  -240,   -30,  -240}
    ,{   -30,   -70,  -240,   -30,  -240}
    ,{  -350,  -700,  -520,  -350,  -520}
    ,{ -1250, -1600, -1420, -1250, -1420}
    ,{  -350,  -700,  -520,  -350,  -520}
    }
   ,{{    40,  -550,  -370,    40,  -370}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{    40,  -550,  -370,    40,  -370}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{   -90,  -680,  -500,   -90,  -500}
    }
   ,{{   410,   230,   410,  -350,   410}
    ,{ -1340, -1700, -1520, -1340, -1520}
    ,{  -350,  -700,  -520,  -350,  -520}
    ,{   410,   230,   410,  -670,   410}
    ,{  -350,  -700,  -520,  -350,  -520}
    }
   ,{{  -190,  -550,  -370,  -190,  -370}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{  -190,  -790,  -610,  -190,  -610}
    ,{  -200,  -550,  -370,  -200,  -370}
    ,{ -1140, -1250, -1310, -1140, -1310}
    }
   }
  ,{{{   410,   410,   410,   410,   410}
    ,{  -240,  -240,  -240,  -240,  -240}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{   410,   410,   410,   410,   410}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{  -240,  -240,  -240,  -240,  -240}
    ,{  -240,  -240,  -240,  -240,  -240}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{ -1180, -1420, -1180, -1420, -1180}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{   410,   410,   410,   410,   410}
    ,{ -1280, -1520, -1280, -1520, -1280}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{   410,   410,   410,   410,   410}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{ -1310, -1310, -1310, -1310, -1310}
    }
   }
  ,{{{  1010,  -150,   410,  1010,   410}
    ,{  1010,  -260,  -240,  1010,  -240}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,  -390,   410,   880,   410}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{  1010,  -260,  -240,  1010,  -240}
    ,{  1010,  -260,  -240,  1010,  -240}
    ,{   730,  -540,  -520,   730,  -520}
    ,{ -1420, -1440, -1420, -1420, -1420}
    ,{   730,  -540,  -520,   730,  -520}
    }
   ,{{   880,  -150,  -370,   880,  -370}
    ,{   880,  -390,  -370,   880,  -370}
    ,{   880,  -150,  -370,   880,  -370}
    ,{   880,  -390,  -370,   880,  -370}
    ,{   750,  -280,  -500,   750,  -500}
    }
   ,{{   730,  -540,   410,   730,   410}
    ,{ -1520, -1530, -1520, -1520, -1520}
    ,{   730,  -540,  -520,   730,  -520}
    ,{   410,  -860,   410,  -840,   410}
    ,{   730,  -540,  -520,   730,  -520}
    }
   ,{{   880,  -380,  -370,   880,  -370}
    ,{   880,  -390,  -370,   880,  -370}
    ,{   640,  -380,  -610,   640,  -610}
    ,{   880,  -390,  -370,   880,  -370}
    ,{ -1310, -1330, -1310, -1310, -1310}
    }
   }
  ,{{{   410,   410,   410,   410,     0}
    ,{     0,  -240,  -240,  -240,     0}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{   410,   410,   410,   410,  -370}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{     0,  -240,  -240,  -240,     0}
    ,{     0,  -240,  -240,  -240,     0}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{ -1180, -1420, -1180, -1420, -1420}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -500,  -500,  -500,  -500,  -500}
    }
   ,{{   410,   410,   410,   410,  -520}
    ,{ -1280, -1520, -1280, -1520, -1520}
    ,{  -520,  -520,  -520,  -520,  -520}
    ,{   410,   410,   410,   410,  -840}
    ,{  -520,  -520,  -520,  -520,  -520}
    }
   ,{{  -370,  -370,  -370,  -370,  -370}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{  -610,  -610,  -610,  -610,  -610}
    ,{  -370,  -370,  -370,  -370,  -370}
    ,{ -1310, -1310, -1310, -1310, -1310}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{   800,   200,  -310,   800,  -310}
    ,{   740,     0,  -510,   740,  -410}
    ,{   800,    50,  -450,   800,  -450}
    ,{   740,   200,  -310,   740,  -310}
    ,{   690,   -50,  -560,   690,  -560}
    }
   ,{{   600,  -140,  -630,   600,  -410}
    ,{   600,  -140,  -650,   600,  -410}
    ,{   290,  -450,  -960,   290,  -960}
    ,{  -360,  -360,  -630,  -870,  -630}
    ,{   290,  -450,  -960,   290,  -960}
    }
   ,{{   740,     0,  -510,   740,  -510}
    ,{   740,     0,  -510,   740,  -510}
    ,{   740,     0,  -510,   740,  -510}
    ,{   740,     0,  -510,   740,  -510}
    ,{   690,   -50,  -560,   690,  -560}
    }
   ,{{   290,   200,  -310,   290,  -310}
    ,{  -640,  -640,  -910, -1150,  -910}
    ,{   290,  -450,  -960,   290,  -960}
    ,{   200,   200,  -310,  -310,  -310}
    ,{   290,  -450,  -960,   290,  -960}
    }
   ,{{   800,    50,  -450,   800,  -450}
    ,{   740,     0,  -510,   740,  -510}
    ,{   800,    50,  -450,   800,  -450}
    ,{   740,     0,  -510,   740,  -510}
    ,{  -550,  -550, -1300, -1300, -1300}
    }
   }
  ,{{{   200,   200,  -310,  -720,  -310}
    ,{     0,     0,  -510, -1020,  -510}
    ,{    50,    50,  -450,  -720,  -450}
    ,{   200,   200,  -310, -1020,  -310}
    ,{   -50,   -50,  -560,  -830,  -560}
    }
   ,{{  -140,  -140,  -650, -1160,  -650}
    ,{  -140,  -140,  -650, -1160,  -650}
    ,{  -450,  -450,  -960, -1470,  -960}
    ,{  -360,  -360,  -870, -1380,  -870}
    ,{  -450,  -450,  -960, -1470,  -960}
    }
   ,{{     0,     0,  -510,  -780,  -510}
    ,{     0,     0,  -510, -1020,  -510}
    ,{     0,     0,  -510,  -780,  -510}
    ,{     0,     0,  -510, -1020,  -510}
    ,{   -50,   -50,  -560,  -830,  -560}
    }
   ,{{   200,   200,  -310, -1470,  -310}
    ,{  -640,  -640, -1150, -1660, -1150}
    ,{  -450,  -450,  -960, -1470,  -960}
    ,{   200,   200,  -310, -2070,  -310}
    ,{  -450,  -450,  -960, -1470,  -960}
    }
   ,{{    50,    50,  -450,  -720,  -450}
    ,{     0,     0,  -510, -1020,  -510}
    ,{    50,    50,  -450,  -720,  -450}
    ,{     0,     0,  -510, -1020,  -510}
    ,{  -550,  -550, -1300, -1810, -1300}
    }
   }
  ,{{{  -310,  -310,  -310,  -310,  -310}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -450,  -450,  -450,  -450,  -450}
    ,{  -310,  -310,  -310,  -310,  -310}
    ,{  -560,  -560,  -560,  -560,  -560}
    }
   ,{{  -630,  -650,  -630,  -650,  -630}
    ,{  -650,  -650,  -650,  -650,  -650}
    ,{  -960,  -960,  -960,  -960,  -960}
    ,{  -630,  -870,  -630,  -870,  -630}
    ,{  -960,  -960,  -960,  -960,  -960}
    }
   ,{{  -510,  -510,  -510,  -510,  -510}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -560,  -560,  -560,  -560,  -560}
    }
   ,{{  -310,  -310,  -310,  -310,  -310}
    ,{  -910, -1150,  -910, -1150,  -910}
    ,{  -960,  -960,  -960,  -960,  -960}
    ,{  -310,  -310,  -310,  -310,  -310}
    ,{  -960,  -960,  -960,  -960,  -960}
    }
   ,{{  -450,  -450,  -450,  -450,  -450}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -450,  -450,  -450,  -450,  -450}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{ -1300, -1300, -1300, -1300, -1300}
    }
   }
  ,{{{   800,  -550,  -310,   800,  -310}
    ,{   740,  -850,  -510,   740,  -510}
    ,{   800,  -550,  -450,   800,  -450}
    ,{   740,  -850,  -310,   740,  -310}
    ,{   690,  -660,  -560,   690,  -560}
    }
   ,{{   600,  -990,  -650,   600,  -650}
    ,{   600,  -990,  -650,   600,  -650}
    ,{   290, -1300,  -960,   290,  -960}
    ,{  -870, -1210,  -870,  -870,  -870}
    ,{   290, -1300,  -960,   290,  -960}
    }
   ,{{   740,  -610,  -510,   740,  -510}
    ,{   740,  -850,  -510,   740,  -510}
    ,{   740,  -610,  -510,   740,  -510}
    ,{   740,  -850,  -510,   740,  -510}
    ,{   690,  -660,  -560,   690,  -560}
    }
   ,{{   290, -1300,  -310,   290,  -310}
    ,{ -1150, -1490, -1150, -1150, -1150}
    ,{   290, -1300,  -960,   290,  -960}
    ,{  -310, -1900,  -310, -1560,  -310}
    ,{   290, -1300,  -960,   290,  -960}
    }
   ,{{   800,  -550,  -450,   800,  -450}
    ,{   740,  -850,  -510,   740,  -510}
    ,{   800,  -550,  -450,   800,  -450}
    ,{   740,  -850,  -510,   740,  -510}
    ,{ -1300, -1640, -1300, -1300, -1300}
    }
   }
  ,{{{  -310,  -310,  -310,  -310,  -410}
    ,{  -410,  -510,  -510,  -510,  -410}
    ,{  -450,  -450,  -450,  -450,  -450}
    ,{  -310,  -310,  -310,  -310,  -510}
    ,{  -560,  -560,  -560,  -560,  -560}
    }
   ,{{  -410,  -650,  -630,  -650,  -410}
    ,{  -410,  -650,  -650,  -650,  -410}
    ,{  -960,  -960,  -960,  -960,  -960}
    ,{  -630,  -870,  -630,  -870,  -870}
    ,{  -960,  -960,  -960,  -960,  -960}
    }
   ,{{  -510,  -510,  -510,  -510,  -510}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -560,  -560,  -560,  -560,  -560}
    }
   ,{{  -310,  -310,  -310,  -310,  -960}
    ,{  -910, -1150,  -910, -1150, -1150}
    ,{  -960,  -960,  -960,  -960,  -960}
    ,{  -310,  -310,  -310,  -310, -1560}
    ,{  -960,  -960,  -960,  -960,  -960}
    }
   ,{{  -450,  -450,  -450,  -450,  -450}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{  -450,  -450,  -450,  -450,  -450}
    ,{  -510,  -510,  -510,  -510,  -510}
    ,{ -1300, -1300, -1300, -1300, -1300}
    }
   }
  }
 ,{{{{   760,   200,  -310,   760,  -250}
    ,{   760,  -340,  -490,   760,  -250}
    ,{   310,  -430,  -940,   310,  -940}
    ,{   400,   200,  -310,   400,  -310}
    ,{   310,  -390,  -940,   310,  -940}
    }
   ,{{   760,  -430,  -490,   760,  -250}
    ,{   760,  -490,  -490,   760,  -250}
    ,{   310,  -430,  -940,   310,  -940}
    ,{ -1170, -1170, -1440, -1680, -1440}
    ,{   310,  -430,  -940,   310,  -940}
    }
   ,{{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{    90,  -650, -1160,    90, -1160}
    ,{   400,  -340,  -850,   400,  -850}
    ,{    90,  -650, -1160,    90, -1160}
    }
   ,{{   310,   200,  -310,   310,  -310}
    ,{  -690,  -690,  -960, -1200,  -960}
    ,{   310,  -430,  -940,   310,  -940}
    ,{   200,   200,  -310,  -310,  -310}
    ,{   310,  -430,  -940,   310,  -940}
    }
   ,{{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{  -220,  -960, -1470,  -220, -1470}
    ,{   400,  -340,  -850,   400,  -850}
    ,{  -390,  -390, -1140, -1140, -1140}
    }
   }
  ,{{{   200,   200,  -310, -1000,  -310}
    ,{  -340,  -340,  -490, -1000,  -490}
    ,{  -430,  -430,  -940, -1430,  -940}
    ,{   200,   200,  -310, -1360,  -310}
    ,{  -390,  -390,  -940, -1430,  -940}
    }
   ,{{  -430,  -430,  -490, -1000,  -490}
    ,{  -490, -2040,  -490, -1000,  -490}
    ,{  -430,  -430,  -940, -1450,  -940}
    ,{ -1170, -1170, -1680, -2190, -1680}
    ,{  -430,  -430,  -940, -1450,  -940}
    }
   ,{{  -340,  -340,  -850, -1360,  -850}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -650,  -650, -1160, -1430, -1160}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -650,  -650, -1160, -1430, -1160}
    }
   ,{{   200,   200,  -310, -1450,  -310}
    ,{  -690,  -690, -1200, -1710, -1200}
    ,{  -430,  -430,  -940, -1450,  -940}
    ,{   200,   200,  -310, -2070,  -310}
    ,{  -430,  -430,  -940, -1450,  -940}
    }
   ,{{  -340,  -340,  -850, -1360,  -850}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -960,  -960, -1470, -1740, -1470}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -390,  -390, -1140, -1650, -1140}
    }
   }
  ,{{{  -310,  -310,  -310,  -310,  -310}
    ,{  -490,  -490,  -490,  -490,  -490}
    ,{  -940,  -940,  -940,  -940,  -940}
    ,{  -310,  -310,  -310,  -310,  -310}
    ,{  -940,  -940,  -940,  -940,  -940}
    }
   ,{{  -490,  -490,  -490,  -490,  -490}
    ,{  -490,  -490,  -490,  -490,  -490}
    ,{  -940,  -940,  -940,  -940,  -940}
    ,{ -1440, -1680, -1440, -1680, -1440}
    ,{  -940,  -940,  -940,  -940,  -940}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1160, -1160, -1160, -1160, -1160}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1160, -1160, -1160, -1160, -1160}
    }
   ,{{  -310,  -310,  -310,  -310,  -310}
    ,{  -960, -1200,  -960, -1200,  -960}
    ,{  -940,  -940,  -940,  -940,  -940}
    ,{  -310,  -310,  -310,  -310,  -310}
    ,{  -940,  -940,  -940,  -940,  -940}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1470, -1470, -1470, -1470, -1470}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1140, -1140, -1140, -1140, -1140}
    }
   }
  ,{{{   760,  -830,  -310,   760,  -310}
    ,{   760,  -830,  -490,   760,  -490}
    ,{   310, -1260,  -940,   310,  -940}
    ,{   400, -1190,  -310,   400,  -310}
    ,{   310, -1260,  -940,   310,  -940}
    }
   ,{{   760,  -830,  -490,   760,  -490}
    ,{   760,  -830,  -490,   760,  -490}
    ,{   310, -1280,  -940,   310,  -940}
    ,{ -1680, -2020, -1680, -1680, -1680}
    ,{   310, -1280,  -940,   310,  -940}
    }
   ,{{   400, -1190,  -850,   400,  -850}
    ,{   400, -1190,  -850,   400,  -850}
    ,{    90, -1260, -1160,    90, -1160}
    ,{   400, -1190,  -850,   400,  -850}
    ,{    90, -1260, -1160,    90, -1160}
    }
   ,{{   310, -1280,  -310,   310,  -310}
    ,{ -1200, -1540, -1200, -1200, -1200}
    ,{   310, -1280,  -940,   310,  -940}
    ,{  -310, -1900,  -310, -1560,  -310}
    ,{   310, -1280,  -940,   310,  -940}
    }
   ,{{   400, -1190,  -850,   400,  -850}
    ,{   400, -1190,  -850,   400,  -850}
    ,{  -220, -1570, -1470,  -220, -1470}
    ,{   400, -1190,  -850,   400,  -850}
    ,{ -1140, -1480, -1140, -1140, -1140}
    }
   }
  ,{{{  -250,  -310,  -310,  -310,  -250}
    ,{  -250,  -490,  -490,  -490,  -250}
    ,{  -940,  -940,  -940,  -940,  -940}
    ,{  -310,  -310,  -310,  -310,  -850}
    ,{  -940,  -940,  -940,  -940,  -940}
    }
   ,{{  -250,  -490,  -490,  -490,  -250}
    ,{  -250,  -490,  -490,  -490,  -250}
    ,{  -940,  -940,  -940,  -940,  -940}
    ,{ -1440, -1680, -1440, -1680, -1680}
    ,{  -940,  -940,  -940,  -940,  -940}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1160, -1160, -1160, -1160, -1160}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1160, -1160, -1160, -1160, -1160}
    }
   ,{{  -310,  -310,  -310,  -310,  -940}
    ,{  -960, -1200,  -960, -1200, -1200}
    ,{  -940,  -940,  -940,  -940,  -940}
    ,{  -310,  -310,  -310,  -310, -1560}
    ,{  -940,  -940,  -940,  -940,  -940}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1470, -1470, -1470, -1470, -1470}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{ -1140, -1140, -1140, -1140, -1140}
    }
   }
  }
 ,{{{{   360,   360,  -150,  -150,  -150}
    ,{   -30,   -30,  -990,  -150,  -990}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{   360,   360,  -150,  -150,  -150}
    ,{  -150,  -650, -1400,  -150, -1400}
    }
   ,{{   -70,   -70, -1180,  -150, -1180}
    ,{   -70,   -70, -1580,  -330, -1340}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{  -910,  -910, -1180, -1420, -1180}
    ,{  -150,  -890, -1400,  -150, -1400}
    }
   ,{{  -150,  -890, -1400,  -150, -1400}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{  -150,  -890, -1400,  -150, -1400}
    }
   ,{{   360,   360,  -150,  -150,  -150}
    ,{   -30,   -30,  -990, -1230,  -990}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{   360,   360,  -150,  -150,  -150}
    ,{  -150,  -890, -1400,  -150, -1400}
    }
   ,{{  -150,  -650, -1400,  -150, -1400}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{  -150,  -890, -1400,  -150, -1400}
    ,{  -650,  -650, -1400, -1400, -1400}
    }
   }
  ,{{{   360,   360,  -150, -1670,  -150}
    ,{   -30,   -30, -1230, -1740, -1230}
    ,{  -890,  -890, -1400, -1670, -1400}
    ,{   360,   360,  -150, -1910,  -150}
    ,{  -650,  -650, -1400, -1670, -1400}
    }
   ,{{   -70,   -70, -1400, -1910, -1400}
    ,{   -70,   -70, -1580, -2090, -1580}
    ,{  -890,  -890, -1400, -1910, -1400}
    ,{  -910,  -910, -1420, -1930, -1420}
    ,{  -890,  -890, -1400, -1910, -1400}
    }
   ,{{  -890,  -890, -1400, -1670, -1400}
    ,{  -890,  -890, -1400, -1910, -1400}
    ,{  -890,  -890, -1400, -1670, -1400}
    ,{  -890,  -890, -1400, -1910, -1400}
    ,{  -890,  -890, -1400, -1670, -1400}
    }
   ,{{   360,   360,  -150, -1740,  -150}
    ,{   -30,   -30, -1230, -1740, -1230}
    ,{  -890,  -890, -1400, -1910, -1400}
    ,{   360,   360,  -150, -1910,  -150}
    ,{  -890,  -890, -1400, -1910, -1400}
    }
   ,{{  -650,  -650, -1400, -1670, -1400}
    ,{  -890,  -890, -1400, -1910, -1400}
    ,{  -890,  -890, -1400, -1670, -1400}
    ,{  -890,  -890, -1400, -1910, -1400}
    ,{  -650,  -650, -1400, -1910, -1400}
    }
   }
  ,{{{  -150,  -150,  -150,  -150,  -150}
    ,{  -990, -1230,  -990, -1230,  -990}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150,  -150}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1180, -1400, -1180, -1400, -1180}
    ,{ -1580, -1580, -1580, -1580, -1580}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1180, -1420, -1180, -1420, -1180}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{  -150,  -150,  -150,  -150,  -150}
    ,{  -990, -1230,  -990, -1230,  -990}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150,  -150}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   }
  ,{{{  -150, -1500,  -150,  -150,  -150}
    ,{  -150, -1570, -1230,  -150, -1230}
    ,{  -150, -1500, -1400,  -150, -1400}
    ,{  -150, -1740,  -150,  -150,  -150}
    ,{  -150, -1500, -1400,  -150, -1400}
    }
   ,{{  -150, -1600, -1400,  -150, -1400}
    ,{  -330, -1600, -1580,  -330, -1580}
    ,{  -150, -1740, -1400,  -150, -1400}
    ,{ -1420, -3040, -1420, -1420, -1420}
    ,{  -150, -1740, -1400,  -150, -1400}
    }
   ,{{  -150, -1500, -1400,  -150, -1400}
    ,{  -150, -1740, -1400,  -150, -1400}
    ,{  -150, -1500, -1400,  -150, -1400}
    ,{  -150, -1740, -1400,  -150, -1400}
    ,{  -150, -1500, -1400,  -150, -1400}
    }
   ,{{  -150, -1570,  -150,  -150,  -150}
    ,{ -1230, -1570, -1230, -1230, -1230}
    ,{  -150, -1740, -1400,  -150, -1400}
    ,{  -150, -1740,  -150, -1400,  -150}
    ,{  -150, -1740, -1400,  -150, -1400}
    }
   ,{{  -150, -1500, -1400,  -150, -1400}
    ,{  -150, -1740, -1400,  -150, -1400}
    ,{  -150, -1500, -1400,  -150, -1400}
    ,{  -150, -1740, -1400,  -150, -1400}
    ,{ -1400, -1740, -1400, -1400, -1400}
    }
   }
  ,{{{  -150,  -150,  -150,  -150, -1230}
    ,{  -990, -1230,  -990, -1230, -1230}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1180, -1400, -1180, -1400, -1340}
    ,{ -1340, -1580, -1580, -1580, -1340}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1180, -1420, -1180, -1420, -1420}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{  -150,  -150,  -150,  -150, -1230}
    ,{  -990, -1230,  -990, -1230, -1230}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{  -150,  -150,  -150,  -150, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   ,{{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    ,{ -1400, -1400, -1400, -1400, -1400}
    }
   }
  }
 ,{{{{   910,   910,   400,   910,   400}
    ,{   910,   170,  -340,   910,  -100}
    ,{   400,  -340,  -850,   400,  -850}
    ,{   910,   910,   400,   400,   400}
    ,{   400,  -100,  -850,   400,  -850}
    }
   ,{{   910,   170,  -340,   910,  -100}
    ,{   910,   170,  -340,   910,  -100}
    ,{   400,  -340,  -850,   400,  -850}
    ,{  -680,  -680,  -950, -1190,  -950}
    ,{   400,  -340,  -850,   400,  -850}
    }
   ,{{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    }
   ,{{   910,   910,   400,   400,   400}
    ,{  -850,  -850, -1120, -1360, -1120}
    ,{   400,  -340,  -850,   400,  -850}
    ,{   910,   910,   400,   400,   400}
    ,{   400,  -340,  -850,   400,  -850}
    }
   ,{{   400,  -100,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{   400,  -340,  -850,   400,  -850}
    ,{  -100,  -100,  -850,  -850,  -850}
    }
   }
  ,{{{   910,   910,   400,  -850,   400}
    ,{   170,   170,  -340,  -850,  -340}
    ,{  -340,  -340,  -850, -1120,  -850}
    ,{   910,   910,   400, -1360,   400}
    ,{  -100,  -100,  -850, -1120,  -850}
    }
   ,{{   170,   170,  -340,  -850,  -340}
    ,{   170,   170,  -340,  -850,  -340}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -680,  -680, -1190, -1700, -1190}
    ,{  -340,  -340,  -850, -1360,  -850}
    }
   ,{{  -340,  -340,  -850, -1120,  -850}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -340,  -340,  -850, -1120,  -850}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -340,  -340,  -850, -1120,  -850}
    }
   ,{{   910,   910,   400, -1360,   400}
    ,{  -850,  -850, -1360, -1870, -1360}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{   910,   910,   400, -1360,   400}
    ,{  -340,  -340,  -850, -1360,  -850}
    }
   ,{{  -100,  -100,  -850, -1120,  -850}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -340,  -340,  -850, -1120,  -850}
    ,{  -340,  -340,  -850, -1360,  -850}
    ,{  -100,  -100,  -850, -1360,  -850}
    }
   }
  ,{{{   400,   400,   400,   400,   400}
    ,{  -340,  -340,  -340,  -340,  -340}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{   400,   400,   400,   400,   400}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{  -340,  -340,  -340,  -340,  -340}
    ,{  -340,  -340,  -340,  -340,  -340}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -950, -1190,  -950, -1190,  -950}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{   400,   400,   400,   400,   400}
    ,{ -1120, -1360, -1120, -1360, -1120}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{   400,   400,   400,   400,   400}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   }
  ,{{{   910,  -680,   400,   910,   400}
    ,{   910,  -680,  -340,   910,  -340}
    ,{   400,  -950,  -850,   400,  -850}
    ,{   400, -1190,   400,   400,   400}
    ,{   400,  -950,  -850,   400,  -850}
    }
   ,{{   910,  -680,  -340,   910,  -340}
    ,{   910,  -680,  -340,   910,  -340}
    ,{   400, -1190,  -850,   400,  -850}
    ,{ -1190, -1530, -1190, -1190, -1190}
    ,{   400, -1190,  -850,   400,  -850}
    }
   ,{{   400,  -950,  -850,   400,  -850}
    ,{   400, -1190,  -850,   400,  -850}
    ,{   400,  -950,  -850,   400,  -850}
    ,{   400, -1190,  -850,   400,  -850}
    ,{   400,  -950,  -850,   400,  -850}
    }
   ,{{   400, -1190,   400,   400,   400}
    ,{ -1360, -1700, -1360, -1360, -1360}
    ,{   400, -1190,  -850,   400,  -850}
    ,{   400, -1190,   400,  -850,   400}
    ,{   400, -1190,  -850,   400,  -850}
    }
   ,{{   400,  -950,  -850,   400,  -850}
    ,{   400, -1190,  -850,   400,  -850}
    ,{   400,  -950,  -850,   400,  -850}
    ,{   400, -1190,  -850,   400,  -850}
    ,{  -850, -1190,  -850,  -850,  -850}
    }
   }
  ,{{{   400,   400,   400,   400,  -100}
    ,{  -100,  -340,  -340,  -340,  -100}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{   400,   400,   400,   400,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{  -100,  -340,  -340,  -340,  -100}
    ,{  -100,  -340,  -340,  -340,  -100}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -950, -1190,  -950, -1190, -1190}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{   400,   400,   400,   400,  -850}
    ,{ -1120, -1360, -1120, -1360, -1360}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{   400,   400,   400,   400,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   ,{{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    ,{  -850,  -850,  -850,  -850,  -850}
    }
   }
  }
 ,{{{{  1490,  1280,   780,  1490,   780}
    ,{  1490,   750,   240,  1490,   480}
    ,{  1200,   450,   -50,  1200,   -50}
    ,{  1280,  1280,   780,  1200,   780}
    ,{  1200,   450,   -50,  1200,   -50}
    }
   ,{{  1490,   750,   240,  1490,   480}
    ,{  1490,   750,   240,  1490,   480}
    ,{  1190,   440,   -60,  1190,   -60}
    ,{  -630,  -630,  -900, -1140,  -900}
    ,{  1190,   440,   -60,  1190,   -60}
    }
   ,{{  1200,   460,   -50,  1200,   -50}
    ,{  1200,   460,   -50,  1200,   -50}
    ,{  1200,   450,   -50,  1200,   -50}
    ,{  1200,   460,   -50,  1200,   -50}
    ,{  1200,   450,   -50,  1200,   -50}
    }
   ,{{  1280,  1280,   780,  1190,   780}
    ,{  -450,  -450,  -720,  -960,  -720}
    ,{  1190,   440,   -60,  1190,   -60}
    ,{  1280,  1280,   780,   780,   780}
    ,{  1190,   440,   -60,  1190,   -60}
    }
   ,{{  1200,   460,   -50,  1200,   -50}
    ,{  1200,   460,   -50,  1200,   -50}
    ,{  1200,   450,   -50,  1200,   -50}
    ,{  1200,   460,   -50,  1200,   -50}
    ,{  -280,  -280, -1030, -1030, -1030}
    }
   }
  ,{{{  1280,  1280,   780,  -260,   780}
    ,{   750,   750,   240,  -260,   240}
    ,{   450,   450,   -50,  -320,   -50}
    ,{  1280,  1280,   780,  -560,   780}
    ,{   450,   450,   -50,  -320,   -50}
    }
   ,{{   750,   750,   240,  -260,   240}
    ,{   750,   750,   240,  -260,   240}
    ,{   440,   440,   -60,  -570,   -60}
    ,{  -630,  -630, -1140, -1650, -1140}
    ,{   440,   440,   -60,  -570,   -60}
    }
   ,{{   460,   460,   -50,  -320,   -50}
    ,{   460,   460,   -50,  -560,   -50}
    ,{   450,   450,   -50,  -320,   -50}
    ,{   460,   460,   -50,  -560,   -50}
    ,{   450,   450,   -50,  -320,   -50}
    }
   ,{{  1280,  1280,   780,  -570,   780}
    ,{  -450,  -450,  -960, -1470,  -960}
    ,{   440,   440,   -60,  -570,   -60}
    ,{  1280,  1280,   780,  -980,   780}
    ,{   440,   440,   -60,  -570,   -60}
    }
   ,{{   460,   460,   -50,  -320,   -50}
    ,{   460,   460,   -50,  -560,   -50}
    ,{   450,   450,   -50,  -320,   -50}
    ,{   460,   460,   -50,  -560,   -50}
    ,{  -280,  -280, -1030, -1540, -1030}
    }
   }
  ,{{{   780,   780,   780,   780,   780}
    ,{   240,   240,   240,   240,   240}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   780,   780,   780,   780,   780}
    ,{   -50,   -50,   -50,   -50,   -50}
    }
   ,{{   240,   240,   240,   240,   240}
    ,{   240,   240,   240,   240,   240}
    ,{   -60,   -60,   -60,   -60,   -60}
    ,{  -900, -1140,  -900, -1140,  -900}
    ,{   -60,   -60,   -60,   -60,   -60}
    }
   ,{{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    }
   ,{{   780,   780,   780,   780,   780}
    ,{  -720,  -960,  -720,  -960,  -720}
    ,{   -60,   -60,   -60,   -60,   -60}
    ,{   780,   780,   780,   780,   780}
    ,{   -60,   -60,   -60,   -60,   -60}
    }
   ,{{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{ -1030, -1030, -1030, -1030, -1030}
    }
   }
  ,{{{  1490,   -90,   780,  1490,   780}
    ,{  1490,   -90,   240,  1490,   240}
    ,{  1200,  -150,   -50,  1200,   -50}
    ,{  1200,  -390,   780,  1200,   780}
    ,{  1200,  -150,   -50,  1200,   -50}
    }
   ,{{  1490,   -90,   240,  1490,   240}
    ,{  1490,   -90,   240,  1490,   240}
    ,{  1190,  -400,   -60,  1190,   -60}
    ,{ -1140, -1480, -1140, -1140, -1140}
    ,{  1190,  -400,   -60,  1190,   -60}
    }
   ,{{  1200,  -150,   -50,  1200,   -50}
    ,{  1200,  -390,   -50,  1200,   -50}
    ,{  1200,  -150,   -50,  1200,   -50}
    ,{  1200,  -390,   -50,  1200,   -50}
    ,{  1200,  -150,   -50,  1200,   -50}
    }
   ,{{  1190,  -400,   780,  1190,   780}
    ,{  -960, -1300,  -960,  -960,  -960}
    ,{  1190,  -400,   -60,  1190,   -60}
    ,{   780,  -810,   780,  -470,   780}
    ,{  1190,  -400,   -60,  1190,   -60}
    }
   ,{{  1200,  -150,   -50,  1200,   -50}
    ,{  1200,  -390,   -50,  1200,   -50}
    ,{  1200,  -150,   -50,  1200,   -50}
    ,{  1200,  -390,   -50,  1200,   -50}
    ,{ -1030, -1370, -1030, -1030, -1030}
    }
   }
  ,{{{   780,   780,   780,   780,   480}
    ,{   480,   240,   240,   240,   480}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   780,   780,   780,   780,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    }
   ,{{   480,   240,   240,   240,   480}
    ,{   480,   240,   240,   240,   480}
    ,{   -60,   -60,   -60,   -60,   -60}
    ,{  -900, -1140,  -900, -1140, -1140}
    ,{   -60,   -60,   -60,   -60,   -60}
    }
   ,{{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    }
   ,{{   780,   780,   780,   780,   -60}
    ,{  -720,  -960,  -720,  -960,  -960}
    ,{   -60,   -60,   -60,   -60,   -60}
    ,{   780,   780,   780,   780,  -470}
    ,{   -60,   -60,   -60,   -60,   -60}
    }
   ,{{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{ -1030, -1030, -1030, -1030, -1030}
    }
   }
  }
 ,{{{{  1560,  1470,   960,  1560,   960}
    ,{  1560,   820,   310,  1560,   550}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1470,  1470,   960,  1430,   960}
    ,{  1300,   560,    50,  1300,    50}
    }
   ,{{  1560,   820,   310,  1560,   550}
    ,{  1560,   820,   310,  1560,   550}
    ,{  1280,   540,    30,  1280,    30}
    ,{  -580,  -580,  -850, -1090,  -850}
    ,{  1280,   540,    30,  1280,    30}
    }
   ,{{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1300,   560,    50,  1300,    50}
    }
   ,{{  1470,  1470,   960,  1280,   960}
    ,{  -880,  -880, -1150, -1390, -1150}
    ,{  1280,   540,    30,  1280,    30}
    ,{  1470,  1470,   960,   960,   960}
    ,{  1280,   540,    30,  1280,    30}
    }
   ,{{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{   990,   250,  -260,   990,  -260}
    ,{  1430,   690,   180,  1430,   180}
    ,{   -10,   -10,  -760,  -760,  -760}
    }
   }
  ,{{{  1470,  1470,   960,   -90,   960}
    ,{   820,   820,   310,  -200,   310}
    ,{   690,   690,   180,   -90,   180}
    ,{  1470,  1470,   960,  -330,   960}
    ,{   560,   560,    50,  -220,    50}
    }
   ,{{   820,   820,   310,  -200,   310}
    ,{   820,   820,   310,  -200,   310}
    ,{   540,   540,    30,  -480,    30}
    ,{  -580,  -580, -1090, -1600, -1090}
    ,{   540,   540,    30,  -480,    30}
    }
   ,{{   690,   690,   180,   -90,   180}
    ,{   690,   690,   180,  -330,   180}
    ,{   690,   690,   180,   -90,   180}
    ,{   690,   690,   180,  -330,   180}
    ,{   560,   560,    50,  -220,    50}
    }
   ,{{  1470,  1470,   960,  -480,   960}
    ,{  -880,  -880, -1390, -1900, -1390}
    ,{   540,   540,    30,  -480,    30}
    ,{  1470,  1470,   960,  -800,   960}
    ,{   540,   540,    30,  -480,    30}
    }
   ,{{   690,   690,   180,  -330,   180}
    ,{   690,   690,   180,  -330,   180}
    ,{   250,   250,  -260,  -530,  -260}
    ,{   690,   690,   180,  -330,   180}
    ,{   -10,   -10,  -760, -1270,  -760}
    }
   }
  ,{{{   960,   960,   960,   960,   960}
    ,{   310,   310,   310,   310,   310}
    ,{   180,   180,   180,   180,   180}
    ,{   960,   960,   960,   960,   960}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{    30,    30,    30,    30,    30}
    ,{  -850, -1090,  -850, -1090,  -850}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   960,   960,   960,   960,   960}
    ,{ -1150, -1390, -1150, -1390, -1150}
    ,{    30,    30,    30,    30,    30}
    ,{   960,   960,   960,   960,   960}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{  -260,  -260,  -260,  -260,  -260}
    ,{   180,   180,   180,   180,   180}
    ,{  -760,  -760,  -760,  -760,  -760}
    }
   }
  ,{{{  1560,    80,   960,  1560,   960}
    ,{  1560,   -30,   310,  1560,   310}
    ,{  1430,    80,   180,  1430,   180}
    ,{  1430,  -160,   960,  1430,   960}
    ,{  1300,   -50,    50,  1300,    50}
    }
   ,{{  1560,   -30,   310,  1560,   310}
    ,{  1560,   -30,   310,  1560,   310}
    ,{  1280,  -310,    30,  1280,    30}
    ,{ -1090, -1430, -1090, -1090, -1090}
    ,{  1280,  -310,    30,  1280,    30}
    }
   ,{{  1430,    80,   180,  1430,   180}
    ,{  1430,  -160,   180,  1430,   180}
    ,{  1430,    80,   180,  1430,   180}
    ,{  1430,  -160,   180,  1430,   180}
    ,{  1300,   -50,    50,  1300,    50}
    }
   ,{{  1280,  -310,   960,  1280,   960}
    ,{ -1390, -1730, -1390, -1390, -1390}
    ,{  1280,  -310,    30,  1280,    30}
    ,{   960,  -630,   960,  -290,   960}
    ,{  1280,  -310,    30,  1280,    30}
    }
   ,{{  1430,  -160,   180,  1430,   180}
    ,{  1430,  -160,   180,  1430,   180}
    ,{   990,  -360,  -260,   990,  -260}
    ,{  1430,  -160,   180,  1430,   180}
    ,{  -760, -1100,  -760,  -760,  -760}
    }
   }
  ,{{{   960,   960,   960,   960,   550}
    ,{   550,   310,   310,   310,   550}
    ,{   180,   180,   180,   180,   180}
    ,{   960,   960,   960,   960,   180}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   550,   310,   310,   310,   550}
    ,{   550,   310,   310,   310,   550}
    ,{    30,    30,    30,    30,    30}
    ,{  -850, -1090,  -850, -1090, -1090}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   960,   960,   960,   960,    30}
    ,{ -1150, -1390, -1150, -1390, -1390}
    ,{    30,    30,    30,    30,    30}
    ,{   960,   960,   960,   960,  -290}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{  -260,  -260,  -260,  -260,  -260}
    ,{   180,   180,   180,   180,   180}
    ,{  -760,  -760,  -760,  -760,  -760}
    }
   }
  }
 ,{{{{  1560,  1470,   960,  1560,   960}
    ,{  1560,   820,   310,  1560,   550}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1470,  1470,   960,  1430,   960}
    ,{  1300,   560,    50,  1300,    50}
    }
   ,{{  1560,   820,   310,  1560,   550}
    ,{  1560,   820,   310,  1560,   550}
    ,{  1280,   540,    30,  1280,    30}
    ,{  -360,  -360,  -630,  -870,  -630}
    ,{  1280,   540,    30,  1280,    30}
    }
   ,{{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1300,   560,    50,  1300,    50}
    }
   ,{{  1470,  1470,   960,  1280,   960}
    ,{   -30,   -30,  -720,  -960,  -720}
    ,{  1280,   540,    30,  1280,    30}
    ,{  1470,  1470,   960,   960,   960}
    ,{  1280,   540,    30,  1280,    30}
    }
   ,{{  1430,   690,   180,  1430,   180}
    ,{  1430,   690,   180,  1430,   180}
    ,{  1200,   450,   -50,  1200,   -50}
    ,{  1430,   690,   180,  1430,   180}
    ,{   -10,   -10,  -760,  -760,  -760}
    }
   }
  ,{{{  1470,  1470,   960,   -90,   960}
    ,{   820,   820,   310,  -200,   310}
    ,{   690,   690,   180,   -90,   180}
    ,{  1470,  1470,   960,  -330,   960}
    ,{   560,   560,    50,  -220,    50}
    }
   ,{{   820,   820,   310,  -200,   310}
    ,{   820,   820,   310,  -200,   310}
    ,{   540,   540,    30,  -480,    30}
    ,{  -360,  -360,  -870, -1380,  -870}
    ,{   540,   540,    30,  -480,    30}
    }
   ,{{   690,   690,   180,   -90,   180}
    ,{   690,   690,   180,  -330,   180}
    ,{   690,   690,   180,   -90,   180}
    ,{   690,   690,   180,  -330,   180}
    ,{   560,   560,    50,  -220,    50}
    }
   ,{{  1470,  1470,   960,  -480,   960}
    ,{   -30,   -30,  -960, -1470,  -960}
    ,{   540,   540,    30,  -480,    30}
    ,{  1470,  1470,   960,  -800,   960}
    ,{   540,   540,    30,  -480,    30}
    }
   ,{{   690,   690,   180,  -320,   180}
    ,{   690,   690,   180,  -330,   180}
    ,{   450,   450,   -50,  -320,   -50}
    ,{   690,   690,   180,  -330,   180}
    ,{   -10,   -10,  -760, -1270,  -760}
    }
   }
  ,{{{   960,   960,   960,   960,   960}
    ,{   310,   310,   310,   310,   310}
    ,{   180,   180,   180,   180,   180}
    ,{   960,   960,   960,   960,   960}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   310,   310,   310,   310,   310}
    ,{   310,   310,   310,   310,   310}
    ,{    30,    30,    30,    30,    30}
    ,{  -630,  -870,  -630,  -870,  -630}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   960,   960,   960,   960,   960}
    ,{  -720,  -960,  -720,  -960,  -720}
    ,{    30,    30,    30,    30,    30}
    ,{   960,   960,   960,   960,   960}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   180,   180,   180,   180,   180}
    ,{  -760,  -760,  -760,  -760,  -760}
    }
   }
  ,{{{  1560,    80,   960,  1560,   960}
    ,{  1560,   -30,   310,  1560,   310}
    ,{  1430,    80,   180,  1430,   180}
    ,{  1430,  -160,   960,  1430,   960}
    ,{  1300,   -50,    50,  1300,    50}
    }
   ,{{  1560,   -30,   310,  1560,   310}
    ,{  1560,   -30,   310,  1560,   310}
    ,{  1280,  -310,    30,  1280,    30}
    ,{  -870, -1210,  -870,  -870,  -870}
    ,{  1280,  -310,    30,  1280,    30}
    }
   ,{{  1430,    80,   180,  1430,   180}
    ,{  1430,  -160,   180,  1430,   180}
    ,{  1430,    80,   180,  1430,   180}
    ,{  1430,  -160,   180,  1430,   180}
    ,{  1300,   -50,    50,  1300,    50}
    }
   ,{{  1280,  -310,   960,  1280,   960}
    ,{  -960, -1300,  -960,  -960,  -960}
    ,{  1280,  -310,    30,  1280,    30}
    ,{   960,  -630,   960,  -290,   960}
    ,{  1280,  -310,    30,  1280,    30}
    }
   ,{{  1430,  -150,   180,  1430,   180}
    ,{  1430,  -160,   180,  1430,   180}
    ,{  1200,  -150,   -50,  1200,   -50}
    ,{  1430,  -160,   180,  1430,   180}
    ,{  -760, -1100,  -760,  -760,  -760}
    }
   }
  ,{{{   960,   960,   960,   960,   550}
    ,{   550,   310,   310,   310,   550}
    ,{   180,   180,   180,   180,   180}
    ,{   960,   960,   960,   960,   180}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   550,   310,   310,   310,   550}
    ,{   550,   310,   310,   310,   550}
    ,{    30,    30,    30,    30,    30}
    ,{  -630,  -870,  -630,  -870,  -870}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{    50,    50,    50,    50,    50}
    }
   ,{{   960,   960,   960,   960,    30}
    ,{  -720,  -960,  -720,  -960,  -960}
    ,{    30,    30,    30,    30,    30}
    ,{   960,   960,   960,   960,  -290}
    ,{    30,    30,    30,    30,    30}
    }
   ,{{   180,   180,   180,   180,   180}
    ,{   180,   180,   180,   180,   180}
    ,{   -50,   -50,   -50,   -50,   -50}
    ,{   180,   180,   180,   180,   180}
    ,{  -760,  -760,  -760,  -760,  -760}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{  1170,   780,   490,  1170,   490}
    ,{  1120,   580,   290,  1120,   290}
    ,{  1170,   640,   340,  1170,   340}
    ,{  1120,   780,   490,  1120,   490}
    ,{  1060,   530,   230,  1060,   230}
    }
   ,{{   970,   440,   170,   970,   170}
    ,{   970,   440,   140,   970,   140}
    ,{   660,   130,  -160,   660,  -160}
    ,{   220,   220,   170,   -80,   170}
    ,{   660,   130,  -160,   660,  -160}
    }
   ,{{  1120,   580,   290,  1120,   290}
    ,{  1120,   580,   290,  1120,   290}
    ,{  1110,   580,   280,  1110,   280}
    ,{  1120,   580,   290,  1120,   290}
    ,{  1060,   530,   230,  1060,   230}
    }
   ,{{   780,   780,   490,   660,   490}
    ,{   -60,   -60,  -120,  -370,  -120}
    ,{   660,   130,  -160,   660,  -160}
    ,{   780,   780,   490,   470,   490}
    ,{   660,   130,  -160,   660,  -160}
    }
   ,{{  1170,   640,   340,  1170,   340}
    ,{  1120,   580,   290,  1120,   290}
    ,{  1170,   640,   340,  1170,   340}
    ,{  1120,   580,   290,  1120,   290}
    ,{    40,    40,  -500,  -510,  -500}
    }
   }
  ,{{{   780,   780,   490,  -330,   490}
    ,{   580,   580,   290,  -620,   290}
    ,{   640,   640,   340,  -330,   340}
    ,{   780,   780,   490,  -620,   490}
    ,{   530,   530,   230,  -440,   230}
    }
   ,{{   440,   440,   140,  -770,   140}
    ,{   440,   440,   140,  -770,   140}
    ,{   130,   130,  -160, -1080,  -160}
    ,{   220,   220,   -70,  -980,   -70}
    ,{   130,   130,  -160, -1080,  -160}
    }
   ,{{   580,   580,   290,  -390,   290}
    ,{   580,   580,   290,  -620,   290}
    ,{   580,   580,   280,  -390,   280}
    ,{   580,   580,   290,  -620,   290}
    ,{   530,   530,   230,  -440,   230}
    }
   ,{{   780,   780,   490, -1080,   490}
    ,{   -60,   -60,  -350, -1270,  -350}
    ,{   130,   130,  -160, -1080,  -160}
    ,{   780,   780,   490, -1680,   490}
    ,{   130,   130,  -160, -1080,  -160}
    }
   ,{{   640,   640,   340,  -330,   340}
    ,{   580,   580,   290,  -620,   290}
    ,{   640,   640,   340,  -330,   340}
    ,{   580,   580,   290,  -620,   290}
    ,{    40,    40,  -500, -1410,  -500}
    }
   }
  ,{{{   480,   470,   480,   470,   480}
    ,{   280,   270,   280,   270,   280}
    ,{   340,   330,   340,   330,   340}
    ,{   480,   470,   480,   470,   480}
    ,{   230,   220,   230,   220,   230}
    }
   ,{{   170,   130,   170,   130,   170}
    ,{   140,   130,   140,   130,   140}
    ,{  -170,  -180,  -170,  -180,  -170}
    ,{   170,   -80,   170,   -80,   170}
    ,{  -170,  -180,  -170,  -180,  -170}
    }
   ,{{   280,   270,   280,   270,   280}
    ,{   280,   270,   280,   270,   280}
    ,{   280,   270,   280,   270,   280}
    ,{   280,   270,   280,   270,   280}
    ,{   230,   220,   230,   220,   230}
    }
   ,{{   480,   470,   480,   470,   480}
    ,{  -120,  -370,  -120,  -370,  -120}
    ,{  -170,  -180,  -170,  -180,  -170}
    ,{   480,   470,   480,   470,   480}
    ,{  -170,  -180,  -170,  -180,  -170}
    }
   ,{{   340,   330,   340,   330,   340}
    ,{   280,   270,   280,   270,   280}
    ,{   340,   330,   340,   330,   340}
    ,{   280,   270,   280,   270,   280}
    ,{  -500,  -510,  -500,  -510,  -500}
    }
   }
  ,{{{  1170,  -510,   490,  1170,   490}
    ,{  1120,  -800,   290,  1120,   290}
    ,{  1170,  -510,   340,  1170,   340}
    ,{  1120,  -800,   490,  1120,   490}
    ,{  1060,  -620,   230,  1060,   230}
    }
   ,{{   970,  -950,   140,   970,   140}
    ,{   970,  -950,   140,   970,   140}
    ,{   660, -1260,  -160,   660,  -160}
    ,{   -70, -1160,   -70,  -490,   -70}
    ,{   660, -1260,  -160,   660,  -160}
    }
   ,{{  1120,  -570,   290,  1120,   290}
    ,{  1120,  -800,   290,  1120,   290}
    ,{  1110,  -570,   280,  1110,   280}
    ,{  1120,  -800,   290,  1120,   290}
    ,{  1060,  -620,   230,  1060,   230}
    }
   ,{{   660, -1260,   490,   660,   490}
    ,{  -350, -1450,  -350,  -780,  -350}
    ,{   660, -1260,  -160,   660,  -160}
    ,{   490, -1860,   490, -1190,   490}
    ,{   660, -1260,  -160,   660,  -160}
    }
   ,{{  1170,  -510,   340,  1170,   340}
    ,{  1120,  -800,   290,  1120,   290}
    ,{  1170,  -510,   340,  1170,   340}
    ,{  1120,  -800,   290,  1120,   290}
    ,{  -500, -1590,  -500,  -920,  -500}
    }
   }
  ,{{{   480,   470,   480,   470,  -600}
    ,{   280,   270,   280,   270,  -600}
    ,{   340,   330,   340,   330,  -640}
    ,{   480,   470,   480,   470,  -690}
    ,{   230,   220,   230,   220,  -750}
    }
   ,{{   170,   130,   170,   130,  -600}
    ,{   140,   130,   140,   130,  -600}
    ,{  -170,  -180,  -170,  -180, -1150}
    ,{   170,   -80,   170,   -80, -1050}
    ,{  -170,  -180,  -170,  -180, -1150}
    }
   ,{{   280,   270,   280,   270,  -690}
    ,{   280,   270,   280,   270,  -690}
    ,{   280,   270,   280,   270,  -700}
    ,{   280,   270,   280,   270,  -690}
    ,{   230,   220,   230,   220,  -750}
    }
   ,{{   480,   470,   480,   470, -1150}
    ,{  -120,  -370,  -120,  -370, -1340}
    ,{  -170,  -180,  -170,  -180, -1150}
    ,{   480,   470,   480,   470, -1750}
    ,{  -170,  -180,  -170,  -180, -1150}
    }
   ,{{   340,   330,   340,   330,  -640}
    ,{   280,   270,   280,   270,  -690}
    ,{   340,   330,   340,   330,  -640}
    ,{   280,   270,   280,   270,  -690}
    ,{  -500,  -510,  -500,  -510, -1480}
    }
   }
  }
 ,{{{{  1140,   780,   490,  1140,   490}
    ,{  1140,   600,   310,  1140,   310}
    ,{   690,   150,  -140,   690,  -140}
    ,{   780,   780,   490,   770,   490}
    ,{   690,   190,  -140,   690,  -140}
    }
   ,{{  1140,   600,   310,  1140,   310}
    ,{  1140,   600,   310,  1140,   310}
    ,{   690,   150,  -140,   690,  -140}
    ,{  -580,  -580,  -640,  -890,  -640}
    ,{   690,   150,  -140,   690,  -140}
    }
   ,{{   770,   240,   -50,   770,   -50}
    ,{   770,   240,   -50,   770,   -50}
    ,{   470,   -60,  -360,   470,  -360}
    ,{   770,   240,   -50,   770,   -50}
    ,{   470,   -60,  -360,   470,  -360}
    }
   ,{{   780,   780,   490,   690,   490}
    ,{  -110,  -110,  -170,  -420,  -170}
    ,{   690,   150,  -140,   690,  -140}
    ,{   780,   780,   490,   470,   490}
    ,{   690,   150,  -140,   690,  -140}
    }
   ,{{   770,   240,   -50,   770,   -50}
    ,{   770,   240,   -50,   770,   -50}
    ,{   160,  -370,  -670,   160,  -670}
    ,{   770,   240,   -50,   770,   -50}
    ,{   190,   190,  -340,  -360,  -340}
    }
   }
  ,{{{   780,   780,   490,  -600,   490}
    ,{   600,   600,   310,  -600,   310}
    ,{   150,   150,  -140, -1030,  -140}
    ,{   780,   780,   490,  -970,   490}
    ,{   190,   190,  -140, -1030,  -140}
    }
   ,{{   600,   600,   310,  -600,   310}
    ,{   600,   600,   310,  -600,   310}
    ,{   150,   150,  -140, -1050,  -140}
    ,{  -580,  -580,  -880, -1790,  -880}
    ,{   150,   150,  -140, -1050,  -140}
    }
   ,{{   240,   240,   -50,  -970,   -50}
    ,{   240,   240,   -50,  -970,   -50}
    ,{   -60,   -60,  -360, -1030,  -360}
    ,{   240,   240,   -50,  -970,   -50}
    ,{   -60,   -60,  -360, -1030,  -360}
    }
   ,{{   780,   780,   490, -1050,   490}
    ,{  -110,  -110,  -400, -1320,  -400}
    ,{   150,   150,  -140, -1050,  -140}
    ,{   780,   780,   490, -1680,   490}
    ,{   150,   150,  -140, -1050,  -140}
    }
   ,{{   240,   240,   -50,  -970,   -50}
    ,{   240,   240,   -50,  -970,   -50}
    ,{  -370,  -370,  -670, -1340,  -670}
    ,{   240,   240,   -50,  -970,   -50}
    ,{   190,   190,  -340, -1260,  -340}
    }
   }
  ,{{{   480,   470,   480,   470,   480}
    ,{   300,   290,   300,   290,   300}
    ,{  -140,  -150,  -140,  -150,  -140}
    ,{   480,   470,   480,   470,   480}
    ,{  -140,  -150,  -140,  -150,  -140}
    }
   ,{{   300,   290,   300,   290,   300}
    ,{   300,   290,   300,   290,   300}
    ,{  -140,  -150,  -140,  -150,  -140}
    ,{  -640,  -890,  -640,  -890,  -640}
    ,{  -140,  -150,  -140,  -150,  -140}
    }
   ,{{   -60,   -70,   -60,   -70,   -60}
    ,{   -60,   -70,   -60,   -70,   -60}
    ,{  -360,  -370,  -360,  -370,  -360}
    ,{   -60,   -70,   -60,   -70,   -60}
    ,{  -360,  -370,  -360,  -370,  -360}
    }
   ,{{   480,   470,   480,   470,   480}
    ,{  -170,  -420,  -170,  -420,  -170}
    ,{  -140,  -150,  -140,  -150,  -140}
    ,{   480,   470,   480,   470,   480}
    ,{  -140,  -150,  -140,  -150,  -140}
    }
   ,{{   -60,   -70,   -60,   -70,   -60}
    ,{   -60,   -70,   -60,   -70,   -60}
    ,{  -670,  -680,  -670,  -680,  -670}
    ,{   -60,   -70,   -60,   -70,   -60}
    ,{  -350,  -360,  -350,  -360,  -350}
    }
   }
  ,{{{  1140,  -780,   490,  1140,   490}
    ,{  1140,  -780,   310,  1140,   310}
    ,{   690, -1210,  -140,   690,  -140}
    ,{   770, -1150,   490,   770,   490}
    ,{   690, -1210,  -140,   690,  -140}
    }
   ,{{  1140,  -780,   310,  1140,   310}
    ,{  1140,  -780,   310,  1140,   310}
    ,{   690, -1230,  -140,   690,  -140}
    ,{  -880, -1970,  -880, -1300,  -880}
    ,{   690, -1230,  -140,   690,  -140}
    }
   ,{{   770, -1150,   -50,   770,   -50}
    ,{   770, -1150,   -50,   770,   -50}
    ,{   470, -1210,  -360,   470,  -360}
    ,{   770, -1150,   -50,   770,   -50}
    ,{   470, -1210,  -360,   470,  -360}
    }
   ,{{   690, -1230,   490,   690,   490}
    ,{  -400, -1500,  -400,  -830,  -400}
    ,{   690, -1230,  -140,   690,  -140}
    ,{   490, -1860,   490, -1190,   490}
    ,{   690, -1230,  -140,   690,  -140}
    }
   ,{{   770, -1150,   -50,   770,   -50}
    ,{   770, -1150,   -50,   770,   -50}
    ,{   160, -1520,  -670,   160,  -670}
    ,{   770, -1150,   -50,   770,   -50}
    ,{  -340, -1440,  -340,  -770,  -340}
    }
   }
  ,{{{   480,   470,   480,   470,  -430}
    ,{   300,   290,   300,   290,  -430}
    ,{  -140,  -150,  -140,  -150, -1120}
    ,{   480,   470,   480,   470, -1040}
    ,{  -140,  -150,  -140,  -150, -1120}
    }
   ,{{   300,   290,   300,   290,  -430}
    ,{   300,   290,   300,   290,  -430}
    ,{  -140,  -150,  -140,  -150, -1120}
    ,{  -640,  -890,  -640,  -890, -1860}
    ,{  -140,  -150,  -140,  -150, -1120}
    }
   ,{{   -60,   -70,   -60,   -70, -1040}
    ,{   -60,   -70,   -60,   -70, -1040}
    ,{  -360,  -370,  -360,  -370, -1340}
    ,{   -60,   -70,   -60,   -70, -1040}
    ,{  -360,  -370,  -360,  -370, -1340}
    }
   ,{{   480,   470,   480,   470, -1120}
    ,{  -170,  -420,  -170,  -420, -1390}
    ,{  -140,  -150,  -140,  -150, -1120}
    ,{   480,   470,   480,   470, -1750}
    ,{  -140,  -150,  -140,  -150, -1120}
    }
   ,{{   -60,   -70,   -60,   -70, -1040}
    ,{   -60,   -70,   -60,   -70, -1040}
    ,{  -670,  -680,  -670,  -680, -1650}
    ,{   -60,   -70,   -60,   -70, -1040}
    ,{  -350,  -360,  -350,  -360, -1330}
    }
   }
  }
 ,{{{{   940,   940,   650,   630,   650}
    ,{   220,  -130,  -190,   220,  -190}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   940,   940,   650,   630,   650}
    ,{   220,   -70,  -600,   220,  -600}
    }
   ,{{   220,  -310,  -380,   220,  -380}
    ,{    40,  -490,  -780,    40,  -780}
    ,{   220,  -310,  -600,   220,  -600}
    ,{  -320,  -320,  -380,  -630,  -380}
    ,{   220,  -310,  -600,   220,  -600}
    }
   ,{{   220,  -310,  -600,   220,  -600}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   220,  -310,  -600,   220,  -600}
    }
   ,{{   940,   940,   650,   630,   650}
    ,{  -130,  -130,  -190,  -440,  -190}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   940,   940,   650,   630,   650}
    ,{   220,  -310,  -600,   220,  -600}
    }
   ,{{   220,   -70,  -600,   220,  -600}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   220,  -310,  -600,   220,  -600}
    ,{   -70,   -70,  -600,  -620,  -600}
    }
   }
  ,{{{   940,   940,   650, -1280,   650}
    ,{  -130,  -130,  -430, -1340,  -430}
    ,{  -310,  -310,  -600, -1280,  -600}
    ,{   940,   940,   650, -1520,   650}
    ,{   -70,   -70,  -600, -1280,  -600}
    }
   ,{{  -310,  -310,  -600, -1520,  -600}
    ,{  -490,  -490,  -780, -1700,  -780}
    ,{  -310,  -310,  -600, -1520,  -600}
    ,{  -320,  -320,  -620, -1530,  -620}
    ,{  -310,  -310,  -600, -1520,  -600}
    }
   ,{{  -310,  -310,  -600, -1280,  -600}
    ,{  -310,  -310,  -600, -1520,  -600}
    ,{  -310,  -310,  -600, -1280,  -600}
    ,{  -310,  -310,  -600, -1520,  -600}
    ,{  -310,  -310,  -600, -1280,  -600}
    }
   ,{{   940,   940,   650, -1340,   650}
    ,{  -130,  -130,  -430, -1340,  -430}
    ,{  -310,  -310,  -600, -1520,  -600}
    ,{   940,   940,   650, -1520,   650}
    ,{  -310,  -310,  -600, -1520,  -600}
    }
   ,{{   -70,   -70,  -600, -1280,  -600}
    ,{  -310,  -310,  -600, -1520,  -600}
    ,{  -310,  -310,  -600, -1280,  -600}
    ,{  -310,  -310,  -600, -1520,  -600}
    ,{   -70,   -70,  -600, -1520,  -600}
    }
   }
  ,{{{   640,   630,   640,   630,   640}
    ,{  -190,  -440,  -190,  -440,  -190}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{   640,   630,   640,   630,   640}
    ,{  -610,  -620,  -610,  -620,  -610}
    }
   ,{{  -380,  -620,  -380,  -620,  -380}
    ,{  -790,  -800,  -790,  -800,  -790}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{  -380,  -630,  -380,  -630,  -380}
    ,{  -610,  -620,  -610,  -620,  -610}
    }
   ,{{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    }
   ,{{   640,   630,   640,   630,   640}
    ,{  -190,  -440,  -190,  -440,  -190}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{   640,   630,   640,   630,   640}
    ,{  -610,  -620,  -610,  -620,  -610}
    }
   ,{{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    ,{  -610,  -620,  -610,  -620,  -610}
    }
   }
  ,{{{   650, -1460,   650,   220,   650}
    ,{   220, -1520,  -430,   220,  -430}
    ,{   220, -1460,  -600,   220,  -600}
    ,{   650, -1700,   650,   220,   650}
    ,{   220, -1460,  -600,   220,  -600}
    }
   ,{{   220, -1700,  -600,   220,  -600}
    ,{    40, -1880,  -780,    40,  -780}
    ,{   220, -1700,  -600,   220,  -600}
    ,{  -620, -1710,  -620, -1040,  -620}
    ,{   220, -1700,  -600,   220,  -600}
    }
   ,{{   220, -1460,  -600,   220,  -600}
    ,{   220, -1700,  -600,   220,  -600}
    ,{   220, -1460,  -600,   220,  -600}
    ,{   220, -1700,  -600,   220,  -600}
    ,{   220, -1460,  -600,   220,  -600}
    }
   ,{{   650, -1520,   650,   220,   650}
    ,{  -430, -1520,  -430,  -850,  -430}
    ,{   220, -1700,  -600,   220,  -600}
    ,{   650, -1700,   650, -1030,   650}
    ,{   220, -1700,  -600,   220,  -600}
    }
   ,{{   220, -1460,  -600,   220,  -600}
    ,{   220, -1700,  -600,   220,  -600}
    ,{   220, -1460,  -600,   220,  -600}
    ,{   220, -1700,  -600,   220,  -600}
    ,{  -600, -1700,  -600, -1030,  -600}
    }
   }
  ,{{{   640,   630,   640,   630, -1410}
    ,{  -190,  -440,  -190,  -440, -1410}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{   640,   630,   640,   630, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    }
   ,{{  -380,  -620,  -380,  -620, -1530}
    ,{  -790,  -800,  -790,  -800, -1530}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{  -380,  -630,  -380,  -630, -1600}
    ,{  -610,  -620,  -610,  -620, -1590}
    }
   ,{{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    }
   ,{{   640,   630,   640,   630, -1410}
    ,{  -190,  -440,  -190,  -440, -1410}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{   640,   630,   640,   630, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    }
   ,{{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    ,{  -610,  -620,  -610,  -620, -1590}
    }
   }
  }
 ,{{{{  1490,  1490,  1200,  1280,  1200}
    ,{  1280,   750,   460,  1280,   460}
    ,{   780,   240,   -50,   780,   -50}
    ,{  1490,  1490,  1200,  1190,  1200}
    ,{   780,   480,   -50,   780,   -50}
    }
   ,{{  1280,   750,   460,  1280,   460}
    ,{  1280,   750,   460,  1280,   460}
    ,{   780,   240,   -50,   780,   -50}
    ,{   -90,   -90,  -150,  -400,  -150}
    ,{   780,   240,   -50,   780,   -50}
    }
   ,{{   780,   240,   -50,   780,   -50}
    ,{   780,   240,   -50,   780,   -50}
    ,{   780,   240,   -50,   780,   -50}
    ,{   780,   240,   -50,   780,   -50}
    ,{   780,   240,   -50,   780,   -50}
    }
   ,{{  1490,  1490,  1200,  1190,  1200}
    ,{  -260,  -260,  -320,  -570,  -320}
    ,{   780,   240,   -50,   780,   -50}
    ,{  1490,  1490,  1200,  1190,  1200}
    ,{   780,   240,   -50,   780,   -50}
    }
   ,{{   780,   480,   -50,   780,   -50}
    ,{   780,   240,   -50,   780,   -50}
    ,{   780,   240,   -50,   780,   -50}
    ,{   780,   240,   -50,   780,   -50}
    ,{   480,   480,   -50,   -60,   -50}
    }
   }
  ,{{{  1490,  1490,  1200,  -450,  1200}
    ,{   750,   750,   460,  -450,   460}
    ,{   240,   240,   -50,  -720,   -50}
    ,{  1490,  1490,  1200,  -960,  1200}
    ,{   480,   480,   -50,  -720,   -50}
    }
   ,{{   750,   750,   460,  -450,   460}
    ,{   750,   750,   460,  -450,   460}
    ,{   240,   240,   -50,  -960,   -50}
    ,{   -90,   -90,  -390, -1300,  -390}
    ,{   240,   240,   -50,  -960,   -50}
    }
   ,{{   240,   240,   -50,  -720,   -50}
    ,{   240,   240,   -50,  -960,   -50}
    ,{   240,   240,   -50,  -720,   -50}
    ,{   240,   240,   -50,  -960,   -50}
    ,{   240,   240,   -50,  -720,   -50}
    }
   ,{{  1490,  1490,  1200,  -960,  1200}
    ,{  -260,  -260,  -560, -1470,  -560}
    ,{   240,   240,   -50,  -960,   -50}
    ,{  1490,  1490,  1200,  -960,  1200}
    ,{   240,   240,   -50,  -960,   -50}
    }
   ,{{   480,   480,   -50,  -720,   -50}
    ,{   240,   240,   -50,  -960,   -50}
    ,{   240,   240,   -50,  -720,   -50}
    ,{   240,   240,   -50,  -960,   -50}
    ,{   480,   480,   -50,  -960,   -50}
    }
   }
  ,{{{  1200,  1190,  1200,  1190,  1200}
    ,{   450,   440,   450,   440,   450}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{  1200,  1190,  1200,  1190,  1200}
    ,{   -50,   -60,   -50,   -60,   -50}
    }
   ,{{   450,   440,   450,   440,   450}
    ,{   450,   440,   450,   440,   450}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{  -150,  -400,  -150,  -400,  -150}
    ,{   -50,   -60,   -50,   -60,   -50}
    }
   ,{{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    }
   ,{{  1200,  1190,  1200,  1190,  1200}
    ,{  -320,  -570,  -320,  -570,  -320}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{  1200,  1190,  1200,  1190,  1200}
    ,{   -50,   -60,   -50,   -60,   -50}
    }
   ,{{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    ,{   -50,   -60,   -50,   -60,   -50}
    }
   }
  ,{{{  1280,  -630,  1200,  1280,  1200}
    ,{  1280,  -630,   460,  1280,   460}
    ,{   780,  -900,   -50,   780,   -50}
    ,{  1200, -1140,  1200,   780,  1200}
    ,{   780,  -900,   -50,   780,   -50}
    }
   ,{{  1280,  -630,   460,  1280,   460}
    ,{  1280,  -630,   460,  1280,   460}
    ,{   780, -1140,   -50,   780,   -50}
    ,{  -390, -1480,  -390,  -810,  -390}
    ,{   780, -1140,   -50,   780,   -50}
    }
   ,{{   780,  -900,   -50,   780,   -50}
    ,{   780, -1140,   -50,   780,   -50}
    ,{   780,  -900,   -50,   780,   -50}
    ,{   780, -1140,   -50,   780,   -50}
    ,{   780,  -900,   -50,   780,   -50}
    }
   ,{{  1200, -1140,  1200,   780,  1200}
    ,{  -560, -1650,  -560,  -980,  -560}
    ,{   780, -1140,   -50,   780,   -50}
    ,{  1200, -1140,  1200,  -470,  1200}
    ,{   780, -1140,   -50,   780,   -50}
    }
   ,{{   780,  -900,   -50,   780,   -50}
    ,{   780, -1140,   -50,   780,   -50}
    ,{   780,  -900,   -50,   780,   -50}
    ,{   780, -1140,   -50,   780,   -50}
    ,{   -50, -1140,   -50,  -470,   -50}
    }
   }
  ,{{{  1200,  1190,  1200,  1190,  -280}
    ,{   450,   440,   450,   440,  -280}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{  1200,  1190,  1200,  1190, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    }
   ,{{   450,   440,   450,   440,  -280}
    ,{   450,   440,   450,   440,  -280}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{  -150,  -400,  -150,  -400, -1370}
    ,{   -50,   -60,   -50,   -60, -1030}
    }
   ,{{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    }
   ,{{  1200,  1190,  1200,  1190, -1030}
    ,{  -320,  -570,  -320,  -570, -1540}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{  1200,  1190,  1200,  1190, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    }
   ,{{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    ,{   -50,   -60,   -50,   -60, -1030}
    }
   }
  }
 ,{{{{  1870,  1870,  1570,  1870,  1570}
    ,{  1870,  1340,  1040,  1870,  1040}
    ,{  1570,  1040,   740,  1570,   740}
    ,{  1870,  1870,  1570,  1570,  1570}
    ,{  1570,  1040,   740,  1570,   740}
    }
   ,{{  1870,  1340,  1040,  1870,  1040}
    ,{  1870,  1340,  1040,  1870,  1040}
    ,{  1560,  1030,   730,  1560,   730}
    ,{   -50,   -50,  -110,  -360,  -110}
    ,{  1560,  1030,   730,  1560,   730}
    }
   ,{{  1570,  1040,   750,  1570,   750}
    ,{  1570,  1040,   750,  1570,   750}
    ,{  1570,  1040,   740,  1570,   740}
    ,{  1570,  1040,   750,  1570,   750}
    ,{  1570,  1040,   740,  1570,   740}
    }
   ,{{  1870,  1870,  1570,  1560,  1570}
    ,{   130,   130,    70,  -180,    70}
    ,{  1560,  1030,   730,  1560,   730}
    ,{  1870,  1870,  1570,  1560,  1570}
    ,{  1560,  1030,   730,  1560,   730}
    }
   ,{{  1570,  1040,   750,  1570,   750}
    ,{  1570,  1040,   750,  1570,   750}
    ,{  1570,  1040,   740,  1570,   740}
    ,{  1570,  1040,   750,  1570,   750}
    ,{   300,   300,  -230,  -250,  -230}
    }
   }
  ,{{{  1870,  1870,  1570,   130,  1570}
    ,{  1340,  1340,  1040,   130,  1040}
    ,{  1040,  1040,   740,    70,   740}
    ,{  1870,  1870,  1570,  -160,  1570}
    ,{  1040,  1040,   740,    70,   740}
    }
   ,{{  1340,  1340,  1040,   130,  1040}
    ,{  1340,  1340,  1040,   130,  1040}
    ,{  1030,  1030,   730,  -180,   730}
    ,{   -50,   -50,  -340, -1260,  -340}
    ,{  1030,  1030,   730,  -180,   730}
    }
   ,{{  1040,  1040,   750,    70,   750}
    ,{  1040,  1040,   750,  -160,   750}
    ,{  1040,  1040,   740,    70,   740}
    ,{  1040,  1040,   750,  -160,   750}
    ,{  1040,  1040,   740,    70,   740}
    }
   ,{{  1870,  1870,  1570,  -180,  1570}
    ,{   130,   130,  -160, -1080,  -160}
    ,{  1030,  1030,   730,  -180,   730}
    ,{  1870,  1870,  1570,  -590,  1570}
    ,{  1030,  1030,   730,  -180,   730}
    }
   ,{{  1040,  1040,   750,    70,   750}
    ,{  1040,  1040,   750,  -160,   750}
    ,{  1040,  1040,   740,    70,   740}
    ,{  1040,  1040,   750,  -160,   750}
    ,{   300,   300,  -230, -1150,  -230}
    }
   }
  ,{{{  1570,  1560,  1570,  1560,  1570}
    ,{  1040,  1030,  1040,  1030,  1040}
    ,{   740,   730,   740,   730,   740}
    ,{  1570,  1560,  1570,  1560,  1570}
    ,{   740,   730,   740,   730,   740}
    }
   ,{{  1040,  1030,  1040,  1030,  1040}
    ,{  1040,  1030,  1040,  1030,  1040}
    ,{   730,   720,   730,   720,   730}
    ,{  -110,  -360,  -110,  -360,  -110}
    ,{   730,   720,   730,   720,   730}
    }
   ,{{   740,   730,   740,   730,   740}
    ,{   740,   730,   740,   730,   740}
    ,{   740,   730,   740,   730,   740}
    ,{   740,   730,   740,   730,   740}
    ,{   740,   730,   740,   730,   740}
    }
   ,{{  1570,  1560,  1570,  1560,  1570}
    ,{    70,  -180,    70,  -180,    70}
    ,{   730,   720,   730,   720,   730}
    ,{  1570,  1560,  1570,  1560,  1570}
    ,{   730,   720,   730,   720,   730}
    }
   ,{{   740,   730,   740,   730,   740}
    ,{   740,   730,   740,   730,   740}
    ,{   740,   730,   740,   730,   740}
    ,{   740,   730,   740,   730,   740}
    ,{  -240,  -250,  -240,  -250,  -240}
    }
   }
  ,{{{  1870,   -50,  1570,  1870,  1570}
    ,{  1870,   -50,  1040,  1870,  1040}
    ,{  1570,  -110,   740,  1570,   740}
    ,{  1570,  -340,  1570,  1570,  1570}
    ,{  1570,  -110,   740,  1570,   740}
    }
   ,{{  1870,   -50,  1040,  1870,  1040}
    ,{  1870,   -50,  1040,  1870,  1040}
    ,{  1560,  -360,   730,  1560,   730}
    ,{  -340, -1440,  -340,  -770,  -340}
    ,{  1560,  -360,   730,  1560,   730}
    }
   ,{{  1570,  -110,   750,  1570,   750}
    ,{  1570,  -340,   750,  1570,   750}
    ,{  1570,  -110,   740,  1570,   740}
    ,{  1570,  -340,   750,  1570,   750}
    ,{  1570,  -110,   740,  1570,   740}
    }
   ,{{  1570,  -360,  1570,  1560,  1570}
    ,{  -160, -1260,  -160,  -590,  -160}
    ,{  1560,  -360,   730,  1560,   730}
    ,{  1570,  -770,  1570,  -100,  1570}
    ,{  1560,  -360,   730,  1560,   730}
    }
   ,{{  1570,  -110,   750,  1570,   750}
    ,{  1570,  -340,   750,  1570,   750}
    ,{  1570,  -110,   740,  1570,   740}
    ,{  1570,  -340,   750,  1570,   750}
    ,{  -230, -1330,  -230,  -660,  -230}
    }
   }
  ,{{{  1570,  1560,  1570,  1560,   300}
    ,{  1040,  1030,  1040,  1030,   300}
    ,{   740,   730,   740,   730,  -240}
    ,{  1570,  1560,  1570,  1560,  -230}
    ,{   740,   730,   740,   730,  -240}
    }
   ,{{  1040,  1030,  1040,  1030,   300}
    ,{  1040,  1030,  1040,  1030,   300}
    ,{   730,   720,   730,   720,  -250}
    ,{  -110,  -360,  -110,  -360, -1330}
    ,{   730,   720,   730,   720,  -250}
    }
   ,{{   740,   730,   740,   730,  -230}
    ,{   740,   730,   740,   730,  -230}
    ,{   740,   730,   740,   730,  -240}
    ,{   740,   730,   740,   730,  -230}
    ,{   740,   730,   740,   730,  -240}
    }
   ,{{  1570,  1560,  1570,  1560,  -250}
    ,{    70,  -180,    70,  -180, -1150}
    ,{   730,   720,   730,   720,  -250}
    ,{  1570,  1560,  1570,  1560,  -660}
    ,{   730,   720,   730,   720,  -250}
    }
   ,{{   740,   730,   740,   730,  -230}
    ,{   740,   730,   740,   730,  -230}
    ,{   740,   730,   740,   730,  -240}
    ,{   740,   730,   740,   730,  -230}
    ,{  -240,  -250,  -240,  -250, -1220}
    }
   }
  }
 ,{{{{  2050,  2050,  1760,  1930,  1760}
    ,{  1930,  1400,  1110,  1930,  1110}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  2050,  2050,  1760,  1800,  1760}
    ,{  1670,  1140,   850,  1670,   850}
    }
   ,{{  1930,  1400,  1110,  1930,  1110}
    ,{  1930,  1400,  1110,  1930,  1110}
    ,{  1650,  1120,   830,  1650,   830}
    ,{     0,     0,   -60,  -310,   -60}
    ,{  1650,  1120,   830,  1650,   830}
    }
   ,{{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1670,  1140,   850,  1670,   850}
    }
   ,{{  2050,  2050,  1760,  1740,  1760}
    ,{  -300,  -300,  -360,  -610,  -360}
    ,{  1650,  1120,   830,  1650,   830}
    ,{  2050,  2050,  1760,  1740,  1760}
    ,{  1650,  1120,   830,  1650,   830}
    }
   ,{{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1360,   830,   540,  1360,   540}
    ,{  1800,  1270,   980,  1800,   980}
    ,{   570,   570,    40,    20,    40}
    }
   }
  ,{{{  2050,  2050,  1760,   300,  1760}
    ,{  1400,  1400,  1110,   190,  1110}
    ,{  1270,  1270,   980,   300,   980}
    ,{  2050,  2050,  1760,    60,  1760}
    ,{  1140,  1140,   850,   180,   850}
    }
   ,{{  1400,  1400,  1110,   190,  1110}
    ,{  1400,  1400,  1110,   190,  1110}
    ,{  1120,  1120,   830,   -80,   830}
    ,{     0,     0,  -290, -1210,  -290}
    ,{  1120,  1120,   830,   -80,   830}
    }
   ,{{  1270,  1270,   980,   300,   980}
    ,{  1270,  1270,   980,    60,   980}
    ,{  1270,  1270,   980,   300,   980}
    ,{  1270,  1270,   980,    60,   980}
    ,{  1140,  1140,   850,   180,   850}
    }
   ,{{  2050,  2050,  1760,   -80,  1760}
    ,{  -300,  -300,  -590, -1510,  -590}
    ,{  1120,  1120,   830,   -80,   830}
    ,{  2050,  2050,  1760,  -400,  1760}
    ,{  1120,  1120,   830,   -80,   830}
    }
   ,{{  1270,  1270,   980,    60,   980}
    ,{  1270,  1270,   980,    60,   980}
    ,{   830,   830,   540,  -130,   540}
    ,{  1270,  1270,   980,    60,   980}
    ,{   570,   570,    40,  -870,    40}
    }
   }
  ,{{{  1750,  1740,  1750,  1740,  1750}
    ,{  1100,  1090,  1100,  1090,  1100}
    ,{   970,   960,   970,   960,   970}
    ,{  1750,  1740,  1750,  1740,  1750}
    ,{   840,   830,   840,   830,   840}
    }
   ,{{  1100,  1090,  1100,  1090,  1100}
    ,{  1100,  1090,  1100,  1090,  1100}
    ,{   820,   810,   820,   810,   820}
    ,{   -60,  -310,   -60,  -310,   -60}
    ,{   820,   810,   820,   810,   820}
    }
   ,{{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   840,   830,   840,   830,   840}
    }
   ,{{  1750,  1740,  1750,  1740,  1750}
    ,{  -360,  -610,  -360,  -610,  -360}
    ,{   820,   810,   820,   810,   820}
    ,{  1750,  1740,  1750,  1740,  1750}
    ,{   820,   810,   820,   810,   820}
    }
   ,{{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   530,   520,   530,   520,   530}
    ,{   970,   960,   970,   960,   970}
    ,{    30,    20,    30,    20,    30}
    }
   }
  ,{{{  1930,   130,  1760,  1930,  1760}
    ,{  1930,    10,  1110,  1930,  1110}
    ,{  1800,   130,   980,  1800,   980}
    ,{  1800,  -110,  1760,  1800,  1760}
    ,{  1670,     0,   850,  1670,   850}
    }
   ,{{  1930,    10,  1110,  1930,  1110}
    ,{  1930,    10,  1110,  1930,  1110}
    ,{  1650,  -260,   830,  1650,   830}
    ,{  -290, -1390,  -290,  -720,  -290}
    ,{  1650,  -260,   830,  1650,   830}
    }
   ,{{  1800,   130,   980,  1800,   980}
    ,{  1800,  -110,   980,  1800,   980}
    ,{  1800,   130,   980,  1800,   980}
    ,{  1800,  -110,   980,  1800,   980}
    ,{  1670,     0,   850,  1670,   850}
    }
   ,{{  1760,  -260,  1760,  1650,  1760}
    ,{  -590, -1690,  -590, -1020,  -590}
    ,{  1650,  -260,   830,  1650,   830}
    ,{  1760,  -580,  1760,    80,  1760}
    ,{  1650,  -260,   830,  1650,   830}
    }
   ,{{  1800,  -110,   980,  1800,   980}
    ,{  1800,  -110,   980,  1800,   980}
    ,{  1360,  -310,   540,  1360,   540}
    ,{  1800,  -110,   980,  1800,   980}
    ,{    40, -1050,    40,  -380,    40}
    }
   }
  ,{{{  1750,  1740,  1750,  1740,   360}
    ,{  1100,  1090,  1100,  1090,   360}
    ,{   970,   960,   970,   960,     0}
    ,{  1750,  1740,  1750,  1740,     0}
    ,{   840,   830,   840,   830,  -130}
    }
   ,{{  1100,  1090,  1100,  1090,   360}
    ,{  1100,  1090,  1100,  1090,   360}
    ,{   820,   810,   820,   810,  -150}
    ,{   -60,  -310,   -60,  -310, -1280}
    ,{   820,   810,   820,   810,  -150}
    }
   ,{{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   840,   830,   840,   830,  -130}
    }
   ,{{  1750,  1740,  1750,  1740,  -150}
    ,{  -360,  -610,  -360,  -610, -1580}
    ,{   820,   810,   820,   810,  -150}
    ,{  1750,  1740,  1750,  1740,  -470}
    ,{   820,   810,   820,   810,  -150}
    }
   ,{{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   530,   520,   530,   520,  -440}
    ,{   970,   960,   970,   960,     0}
    ,{    30,    20,    30,    20,  -940}
    }
   }
  }
 ,{{{{  2050,  2050,  1760,  1930,  1760}
    ,{  1930,  1400,  1110,  1930,  1110}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  2050,  2050,  1760,  1800,  1760}
    ,{  1670,  1140,   850,  1670,   850}
    }
   ,{{  1930,  1400,  1110,  1930,  1110}
    ,{  1930,  1400,  1110,  1930,  1110}
    ,{  1650,  1120,   830,  1650,   830}
    ,{   220,   220,   170,   -80,   170}
    ,{  1650,  1120,   830,  1650,   830}
    }
   ,{{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1670,  1140,   850,  1670,   850}
    }
   ,{{  2050,  2050,  1760,  1740,  1760}
    ,{   130,   130,    70,  -180,    70}
    ,{  1650,  1120,   830,  1650,   830}
    ,{  2050,  2050,  1760,  1740,  1760}
    ,{  1650,  1120,   830,  1650,   830}
    }
   ,{{  1800,  1270,   980,  1800,   980}
    ,{  1800,  1270,   980,  1800,   980}
    ,{  1570,  1040,   740,  1570,   740}
    ,{  1800,  1270,   980,  1800,   980}
    ,{   570,   570,    40,    20,    40}
    }
   }
  ,{{{  2050,  2050,  1760,   300,  1760}
    ,{  1400,  1400,  1110,   190,  1110}
    ,{  1270,  1270,   980,   300,   980}
    ,{  2050,  2050,  1760,    60,  1760}
    ,{  1140,  1140,   850,   180,   850}
    }
   ,{{  1400,  1400,  1110,   190,  1110}
    ,{  1400,  1400,  1110,   190,  1110}
    ,{  1120,  1120,   830,   -80,   830}
    ,{   220,   220,   -70,  -980,   -70}
    ,{  1120,  1120,   830,   -80,   830}
    }
   ,{{  1270,  1270,   980,   300,   980}
    ,{  1270,  1270,   980,    60,   980}
    ,{  1270,  1270,   980,   300,   980}
    ,{  1270,  1270,   980,    60,   980}
    ,{  1140,  1140,   850,   180,   850}
    }
   ,{{  2050,  2050,  1760,   -80,  1760}
    ,{   130,   130,  -160, -1080,  -160}
    ,{  1120,  1120,   830,   -80,   830}
    ,{  2050,  2050,  1760,  -400,  1760}
    ,{  1120,  1120,   830,   -80,   830}
    }
   ,{{  1270,  1270,   980,    70,   980}
    ,{  1270,  1270,   980,    60,   980}
    ,{  1040,  1040,   740,    70,   740}
    ,{  1270,  1270,   980,    60,   980}
    ,{   570,   570,    40,  -870,    40}
    }
   }
  ,{{{  1750,  1740,  1750,  1740,  1750}
    ,{  1100,  1090,  1100,  1090,  1100}
    ,{   970,   960,   970,   960,   970}
    ,{  1750,  1740,  1750,  1740,  1750}
    ,{   840,   830,   840,   830,   840}
    }
   ,{{  1100,  1090,  1100,  1090,  1100}
    ,{  1100,  1090,  1100,  1090,  1100}
    ,{   820,   810,   820,   810,   820}
    ,{   170,   -80,   170,   -80,   170}
    ,{   820,   810,   820,   810,   820}
    }
   ,{{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   840,   830,   840,   830,   840}
    }
   ,{{  1750,  1740,  1750,  1740,  1750}
    ,{    70,  -180,    70,  -180,    70}
    ,{   820,   810,   820,   810,   820}
    ,{  1750,  1740,  1750,  1740,  1750}
    ,{   820,   810,   820,   810,   820}
    }
   ,{{   970,   960,   970,   960,   970}
    ,{   970,   960,   970,   960,   970}
    ,{   740,   730,   740,   730,   740}
    ,{   970,   960,   970,   960,   970}
    ,{    30,    20,    30,    20,    30}
    }
   }
  ,{{{  1930,   130,  1760,  1930,  1760}
    ,{  1930,    10,  1110,  1930,  1110}
    ,{  1800,   130,   980,  1800,   980}
    ,{  1800,  -110,  1760,  1800,  1760}
    ,{  1670,     0,   850,  1670,   850}
    }
   ,{{  1930,    10,  1110,  1930,  1110}
    ,{  1930,    10,  1110,  1930,  1110}
    ,{  1650,  -260,   830,  1650,   830}
    ,{   -70, -1160,   -70,  -490,   -70}
    ,{  1650,  -260,   830,  1650,   830}
    }
   ,{{  1800,   130,   980,  1800,   980}
    ,{  1800,  -110,   980,  1800,   980}
    ,{  1800,   130,   980,  1800,   980}
    ,{  1800,  -110,   980,  1800,   980}
    ,{  1670,     0,   850,  1670,   850}
    }
   ,{{  1760,  -260,  1760,  1650,  1760}
    ,{  -160, -1260,  -160,  -590,  -160}
    ,{  1650,  -260,   830,  1650,   830}
    ,{  1760,  -580,  1760,    80,  1760}
    ,{  1650,  -260,   830,  1650,   830}
    }
   ,{{  1800,  -110,   980,  1800,   980}
    ,{  1800,  -110,   980,  1800,   980}
    ,{  1570,  -110,   740,  1570,   740}
    ,{  1800,  -110,   980,  1800,   980}
    ,{    40, -1050,    40,  -380,    40}
    }
   }
  ,{{{  1750,  1740,  1750,  1740,   360}
    ,{  1100,  1090,  1100,  1090,   360}
    ,{   970,   960,   970,   960,     0}
    ,{  1750,  1740,  1750,  1740,     0}
    ,{   840,   830,   840,   830,  -130}
    }
   ,{{  1100,  1090,  1100,  1090,   360}
    ,{  1100,  1090,  1100,  1090,   360}
    ,{   820,   810,   820,   810,  -150}
    ,{   170,   -80,   170,   -80, -1050}
    ,{   820,   810,   820,   810,  -150}
    }
   ,{{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   840,   830,   840,   830,  -130}
    }
   ,{{  1750,  1740,  1750,  1740,  -150}
    ,{    70,  -180,    70,  -180, -1150}
    ,{   820,   810,   820,   810,  -150}
    ,{  1750,  1740,  1750,  1740,  -470}
    ,{   820,   810,   820,   810,  -150}
    }
   ,{{   970,   960,   970,   960,     0}
    ,{   970,   960,   970,   960,     0}
    ,{   740,   730,   740,   730,  -240}
    ,{   970,   960,   970,   960,     0}
    ,{    30,    20,    30,    20,  -940}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{  1350,   850,   720,  1350,   720}
    ,{  1300,   650,   520,  1300,   520}
    ,{  1350,   700,   570,  1350,   570}
    ,{  1300,   850,   720,  1300,   720}
    ,{  1250,   590,   460,  1250,   460}
    }
   ,{{  1160,   500,   400,  1160,   370}
    ,{  1160,   500,   370,  1160,   370}
    ,{   850,   190,    60,   850,    60}
    ,{   400,   290,   400,    10,   160}
    ,{   850,   190,    60,   850,    60}
    }
   ,{{  1300,   650,   520,  1300,   520}
    ,{  1300,   650,   520,  1300,   520}
    ,{  1290,   640,   510,  1290,   510}
    ,{  1300,   650,   520,  1300,   520}
    ,{  1250,   590,   460,  1250,   460}
    }
   ,{{   850,   850,   720,   850,   720}
    ,{   120,     0,   120,  -270,  -120}
    ,{   850,   190,    60,   850,    60}
    ,{   850,   850,   720,   570,   720}
    ,{   850,   190,    60,   850,    60}
    }
   ,{{  1350,   700,   570,  1350,   570}
    ,{  1300,   650,   520,  1300,   520}
    ,{  1350,   700,   570,  1350,   570}
    ,{  1300,   650,   520,  1300,   520}
    ,{   100,   100,  -270,  -420,  -270}
    }
   }
  ,{{{   850,   850,   720,  -760,   720}
    ,{   650,   650,   520, -1050,   520}
    ,{   700,   700,   570,  -760,   570}
    ,{   850,   850,   720, -1050,   720}
    ,{   590,   590,   460,  -870,   460}
    }
   ,{{   500,   500,   370, -1200,   370}
    ,{   500,   500,   370, -1200,   370}
    ,{   190,   190,    60, -1510,    60}
    ,{   290,   290,   160, -1410,   160}
    ,{   190,   190,    60, -1510,    60}
    }
   ,{{   650,   650,   520,  -820,   520}
    ,{   650,   650,   520, -1050,   520}
    ,{   640,   640,   510,  -820,   510}
    ,{   650,   650,   520, -1050,   520}
    ,{   590,   590,   460,  -870,   460}
    }
   ,{{   850,   850,   720, -1510,   720}
    ,{     0,     0,  -120, -1700,  -120}
    ,{   190,   190,    60, -1510,    60}
    ,{   850,   850,   720, -2110,   720}
    ,{   190,   190,    60, -1510,    60}
    }
   ,{{   700,   700,   570,  -760,   570}
    ,{   650,   650,   520, -1050,   520}
    ,{   700,   700,   570,  -760,   570}
    ,{   650,   650,   520, -1050,   520}
    ,{   100,   100,  -270, -1840,  -270}
    }
   }
  ,{{{   720,   570,   720,   570,   280}
    ,{   520,   370,   520,   370,    80}
    ,{   570,   420,   570,   420,   130}
    ,{   720,   570,   720,   570,   280}
    ,{   460,   310,   460,   310,    20}
    }
   ,{{   400,   220,   400,   220,   -40}
    ,{   370,   220,   370,   220,   -60}
    ,{    60,   -80,    60,   -80,  -370}
    ,{   400,    10,   400,    10,   -40}
    ,{    60,   -80,    60,   -80,  -370}
    }
   ,{{   520,   370,   520,   370,    80}
    ,{   520,   370,   520,   370,    80}
    ,{   510,   360,   510,   360,    70}
    ,{   520,   370,   520,   370,    80}
    ,{   460,   310,   460,   310,    20}
    }
   ,{{   720,   570,   720,   570,   280}
    ,{   120,  -270,   120,  -270,  -320}
    ,{    60,   -80,    60,   -80,  -370}
    ,{   720,   570,   720,   570,   280}
    ,{    60,   -80,    60,   -80,  -370}
    }
   ,{{   570,   420,   570,   420,   130}
    ,{   520,   370,   520,   370,    80}
    ,{   570,   420,   570,   420,   130}
    ,{   520,   370,   520,   370,    80}
    ,{  -270,  -420,  -270,  -420,  -710}
    }
   }
  ,{{{  1350,  -460,   720,  1350,   720}
    ,{  1300,  -750,   520,  1300,   520}
    ,{  1350,  -460,   570,  1350,   570}
    ,{  1300,  -750,   720,  1300,   720}
    ,{  1250,  -570,   460,  1250,   460}
    }
   ,{{  1160,  -900,   370,  1160,   370}
    ,{  1160,  -900,   370,  1160,   370}
    ,{   850, -1210,    60,   850,    60}
    ,{   160, -1110,   160,  -310,   160}
    ,{   850, -1210,    60,   850,    60}
    }
   ,{{  1300,  -520,   520,  1300,   520}
    ,{  1300,  -750,   520,  1300,   520}
    ,{  1290,  -520,   510,  1290,   510}
    ,{  1300,  -750,   520,  1300,   520}
    ,{  1250,  -570,   460,  1250,   460}
    }
   ,{{   850, -1210,   720,   850,   720}
    ,{  -120, -1400,  -120,  -590,  -120}
    ,{   850, -1210,    60,   850,    60}
    ,{   720, -1810,   720, -1000,   720}
    ,{   850, -1210,    60,   850,    60}
    }
   ,{{  1350,  -460,   570,  1350,   570}
    ,{  1300,  -750,   520,  1300,   520}
    ,{  1350,  -460,   570,  1350,   570}
    ,{  1300,  -750,   520,  1300,   520}
    ,{  -270, -1540,  -270,  -740,  -270}
    }
   }
  ,{{{   590,   570,   590,   570,  -320}
    ,{   390,   370,   390,   370,  -320}
    ,{   440,   420,   440,   420,  -360}
    ,{   590,   570,   590,   570,  -420}
    ,{   330,   310,   330,   310,  -470}
    }
   ,{{   270,   220,   270,   220,  -320}
    ,{   240,   220,   240,   220,  -320}
    ,{   -60,   -80,   -60,   -80,  -870}
    ,{   270,    10,   270,    10,  -780}
    ,{   -60,   -80,   -60,   -80,  -870}
    }
   ,{{   390,   370,   390,   370,  -420}
    ,{   390,   370,   390,   370,  -420}
    ,{   380,   360,   380,   360,  -420}
    ,{   390,   370,   390,   370,  -420}
    ,{   330,   310,   330,   310,  -470}
    }
   ,{{   590,   570,   590,   570,  -870}
    ,{   -10,  -270,   -10,  -270, -1060}
    ,{   -60,   -80,   -60,   -80,  -870}
    ,{   590,   570,   590,   570, -1470}
    ,{   -60,   -80,   -60,   -80,  -870}
    }
   ,{{   440,   420,   440,   420,  -360}
    ,{   390,   370,   390,   370,  -420}
    ,{   440,   420,   440,   420,  -360}
    ,{   390,   370,   390,   370,  -420}
    ,{  -400,  -420,  -400,  -420, -1210}
    }
   }
  }
 ,{{{{  1320,   850,   720,  1320,   720}
    ,{  1320,   670,   540,  1320,   540}
    ,{   870,   220,    90,   870,    90}
    ,{   960,   850,   720,   960,   720}
    ,{   870,   250,    90,   870,    90}
    }
   ,{{  1320,   670,   540,  1320,   540}
    ,{  1320,   670,   540,  1320,   540}
    ,{   870,   220,    90,   870,    90}
    ,{  -410,  -520,  -410,  -800,  -650}
    ,{   870,   220,    90,   870,    90}
    }
   ,{{   960,   300,   170,   960,   170}
    ,{   960,   300,   170,   960,   170}
    ,{   650,     0,  -130,   650,  -130}
    ,{   960,   300,   170,   960,   170}
    ,{   650,     0,  -130,   650,  -130}
    }
   ,{{   870,   850,   720,   870,   720}
    ,{    70,   -40,    70,  -320,  -170}
    ,{   870,   220,    90,   870,    90}
    ,{   850,   850,   720,   570,   720}
    ,{   870,   220,    90,   870,    90}
    }
   ,{{   960,   300,   170,   960,   170}
    ,{   960,   300,   170,   960,   170}
    ,{   340,  -310,  -440,   340,  -440}
    ,{   960,   300,   170,   960,   170}
    ,{   250,   250,  -110,  -260,  -110}
    }
   }
  ,{{{   850,   850,   720, -1030,   720}
    ,{   670,   670,   540, -1030,   540}
    ,{   220,   220,    90, -1460,    90}
    ,{   850,   850,   720, -1400,   720}
    ,{   250,   250,    90, -1460,    90}
    }
   ,{{   670,   670,   540, -1030,   540}
    ,{   670,   670,   540, -1030,   540}
    ,{   220,   220,    90, -1480,    90}
    ,{  -520,  -520,  -650, -2220,  -650}
    ,{   220,   220,    90, -1480,    90}
    }
   ,{{   300,   300,   170, -1400,   170}
    ,{   300,   300,   170, -1400,   170}
    ,{     0,     0,  -130, -1460,  -130}
    ,{   300,   300,   170, -1400,   170}
    ,{     0,     0,  -130, -1460,  -130}
    }
   ,{{   850,   850,   720, -1480,   720}
    ,{   -40,   -40,  -170, -1750,  -170}
    ,{   220,   220,    90, -1480,    90}
    ,{   850,   850,   720, -2110,   720}
    ,{   220,   220,    90, -1480,    90}
    }
   ,{{   300,   300,   170, -1400,   170}
    ,{   300,   300,   170, -1400,   170}
    ,{  -310,  -310,  -440, -1770,  -440}
    ,{   300,   300,   170, -1400,   170}
    ,{   250,   250,  -110, -1690,  -110}
    }
   }
  ,{{{   720,   570,   720,   570,   280}
    ,{   540,   390,   540,   390,   100}
    ,{    90,   -60,    90,   -60,  -350}
    ,{   720,   570,   720,   570,   280}
    ,{    90,   -60,    90,   -60,  -350}
    }
   ,{{   540,   390,   540,   390,   100}
    ,{   540,   390,   540,   390,   100}
    ,{    90,   -60,    90,   -60,  -350}
    ,{  -410,  -800,  -410,  -800,  -850}
    ,{    90,   -60,    90,   -60,  -350}
    }
   ,{{   170,    20,   170,    20,  -260}
    ,{   170,    20,   170,    20,  -260}
    ,{  -130,  -280,  -130,  -280,  -570}
    ,{   170,    20,   170,    20,  -260}
    ,{  -130,  -280,  -130,  -280,  -570}
    }
   ,{{   720,   570,   720,   570,   280}
    ,{    70,  -320,    70,  -320,  -370}
    ,{    90,   -60,    90,   -60,  -350}
    ,{   720,   570,   720,   570,   280}
    ,{    90,   -60,    90,   -60,  -350}
    }
   ,{{   170,    20,   170,    20,  -260}
    ,{   170,    20,   170,    20,  -260}
    ,{  -440,  -590,  -440,  -590,  -880}
    ,{   170,    20,   170,    20,  -260}
    ,{  -110,  -260,  -110,  -260,  -550}
    }
   }
  ,{{{  1320,  -730,   720,  1320,   720}
    ,{  1320,  -730,   540,  1320,   540}
    ,{   870, -1160,    90,   870,    90}
    ,{   960, -1100,   720,   960,   720}
    ,{   870, -1160,    90,   870,    90}
    }
   ,{{  1320,  -730,   540,  1320,   540}
    ,{  1320,  -730,   540,  1320,   540}
    ,{   870, -1180,    90,   870,    90}
    ,{  -650, -1920,  -650, -1120,  -650}
    ,{   870, -1180,    90,   870,    90}
    }
   ,{{   960, -1100,   170,   960,   170}
    ,{   960, -1100,   170,   960,   170}
    ,{   650, -1160,  -130,   650,  -130}
    ,{   960, -1100,   170,   960,   170}
    ,{   650, -1160,  -130,   650,  -130}
    }
   ,{{   870, -1180,   720,   870,   720}
    ,{  -170, -1450,  -170,  -640,  -170}
    ,{   870, -1180,    90,   870,    90}
    ,{   720, -1810,   720, -1000,   720}
    ,{   870, -1180,    90,   870,    90}
    }
   ,{{   960, -1100,   170,   960,   170}
    ,{   960, -1100,   170,   960,   170}
    ,{   340, -1470,  -440,   340,  -440}
    ,{   960, -1100,   170,   960,   170}
    ,{  -110, -1390,  -110,  -580,  -110}
    }
   }
  ,{{{   590,   570,   590,   570,  -160}
    ,{   410,   390,   410,   390,  -160}
    ,{   -40,   -60,   -40,   -60,  -850}
    ,{   590,   570,   590,   570,  -760}
    ,{   -40,   -60,   -40,   -60,  -850}
    }
   ,{{   410,   390,   410,   390,  -160}
    ,{   410,   390,   410,   390,  -160}
    ,{   -40,   -60,   -40,   -60,  -850}
    ,{  -540,  -800,  -540,  -800, -1590}
    ,{   -40,   -60,   -40,   -60,  -850}
    }
   ,{{    40,    20,    40,    20,  -760}
    ,{    40,    20,    40,    20,  -760}
    ,{  -260,  -280,  -260,  -280, -1070}
    ,{    40,    20,    40,    20,  -760}
    ,{  -260,  -280,  -260,  -280, -1070}
    }
   ,{{   590,   570,   590,   570,  -850}
    ,{   -60,  -320,   -60,  -320, -1110}
    ,{   -40,   -60,   -40,   -60,  -850}
    ,{   590,   570,   590,   570, -1470}
    ,{   -40,   -60,   -40,   -60,  -850}
    }
   ,{{    40,    20,    40,    20,  -760}
    ,{    40,    20,    40,    20,  -760}
    ,{  -570,  -590,  -570,  -590, -1380}
    ,{    40,    20,    40,    20,  -760}
    ,{  -240,  -260,  -240,  -260, -1050}
    }
   }
  }
 ,{{{{  1010,  1010,   880,   730,   880}
    ,{   410,   -70,    40,   410,  -200}
    ,{   410,  -240,  -370,   410,  -370}
    ,{  1010,  1010,   880,   730,   880}
    ,{   410,     0,  -370,   410,  -370}
    }
   ,{{   410,  -240,  -150,   410,  -370}
    ,{   230,  -420,  -550,   230,  -550}
    ,{   410,  -240,  -370,   410,  -370}
    ,{  -150,  -260,  -150,  -540,  -390}
    ,{   410,  -240,  -370,   410,  -370}
    }
   ,{{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    }
   ,{{  1010,  1010,   880,   730,   880}
    ,{    40,   -70,    40,  -350,  -200}
    ,{   410,  -240,  -370,   410,  -370}
    ,{  1010,  1010,   880,   730,   880}
    ,{   410,  -240,  -370,   410,  -370}
    }
   ,{{   410,     0,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{     0,     0,  -370,  -520,  -370}
    }
   }
  ,{{{  1010,  1010,   880, -1710,   880}
    ,{   -70,   -70,  -200, -1770,  -200}
    ,{  -240,  -240,  -370, -1710,  -370}
    ,{  1010,  1010,   880, -1950,   880}
    ,{     0,     0,  -370, -1710,  -370}
    }
   ,{{  -240,  -240,  -370, -1950,  -370}
    ,{  -420,  -420,  -550, -2130,  -550}
    ,{  -240,  -240,  -370, -1950,  -370}
    ,{  -260,  -260,  -390, -1960,  -390}
    ,{  -240,  -240,  -370, -1950,  -370}
    }
   ,{{  -240,  -240,  -370, -1710,  -370}
    ,{  -240,  -240,  -370, -1950,  -370}
    ,{  -240,  -240,  -370, -1710,  -370}
    ,{  -240,  -240,  -370, -1950,  -370}
    ,{  -240,  -240,  -370, -1710,  -370}
    }
   ,{{  1010,  1010,   880, -1770,   880}
    ,{   -70,   -70,  -200, -1770,  -200}
    ,{  -240,  -240,  -370, -1950,  -370}
    ,{  1010,  1010,   880, -1950,   880}
    ,{  -240,  -240,  -370, -1950,  -370}
    }
   ,{{     0,     0,  -370, -1710,  -370}
    ,{  -240,  -240,  -370, -1950,  -370}
    ,{  -240,  -240,  -370, -1710,  -370}
    ,{  -240,  -240,  -370, -1950,  -370}
    ,{     0,     0,  -370, -1950,  -370}
    }
   }
  ,{{{   880,   730,   880,   730,   440}
    ,{    40,  -350,    40,  -350,  -400}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{   880,   730,   880,   730,   440}
    ,{  -370,  -520,  -370,  -520,  -810}
    }
   ,{{  -150,  -520,  -150,  -520,  -590}
    ,{  -550,  -700,  -550,  -700,  -990}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{  -150,  -540,  -150,  -540,  -590}
    ,{  -370,  -520,  -370,  -520,  -810}
    }
   ,{{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    }
   ,{{   880,   730,   880,   730,   440}
    ,{    40,  -350,    40,  -350,  -400}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{   880,   730,   880,   730,   440}
    ,{  -370,  -520,  -370,  -520,  -810}
    }
   ,{{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    ,{  -370,  -520,  -370,  -520,  -810}
    }
   }
  ,{{{   880, -1410,   880,   410,   880}
    ,{   410, -1470,  -200,   410,  -200}
    ,{   410, -1410,  -370,   410,  -370}
    ,{   880, -1650,   880,   410,   880}
    ,{   410, -1410,  -370,   410,  -370}
    }
   ,{{   410, -1650,  -370,   410,  -370}
    ,{   230, -1830,  -550,   230,  -550}
    ,{   410, -1650,  -370,   410,  -370}
    ,{  -390, -1660,  -390,  -860,  -390}
    ,{   410, -1650,  -370,   410,  -370}
    }
   ,{{   410, -1410,  -370,   410,  -370}
    ,{   410, -1650,  -370,   410,  -370}
    ,{   410, -1410,  -370,   410,  -370}
    ,{   410, -1650,  -370,   410,  -370}
    ,{   410, -1410,  -370,   410,  -370}
    }
   ,{{   880, -1470,   880,   410,   880}
    ,{  -200, -1470,  -200,  -670,  -200}
    ,{   410, -1650,  -370,   410,  -370}
    ,{   880, -1650,   880,  -840,   880}
    ,{   410, -1650,  -370,   410,  -370}
    }
   ,{{   410, -1410,  -370,   410,  -370}
    ,{   410, -1650,  -370,   410,  -370}
    ,{   410, -1410,  -370,   410,  -370}
    ,{   410, -1650,  -370,   410,  -370}
    ,{  -370, -1650,  -370,  -840,  -370}
    }
   }
  ,{{{   750,   730,   750,   730, -1140}
    ,{   -90,  -350,   -90,  -350, -1140}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{   750,   730,   750,   730, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{  -280,  -520,  -280,  -520, -1250}
    ,{  -680,  -700,  -680,  -700, -1250}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -280,  -540,  -280,  -540, -1330}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{   750,   730,   750,   730, -1140}
    ,{   -90,  -350,   -90,  -350, -1140}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{   750,   730,   750,   730, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   }
  }
 ,{{{{  1560,  1560,  1430,  1470,  1430}
    ,{  1470,   820,   690,  1470,   690}
    ,{   960,   310,   180,   960,   180}
    ,{  1560,  1560,  1430,  1280,  1430}
    ,{   960,   550,   180,   960,   180}
    }
   ,{{  1470,   820,   690,  1470,   690}
    ,{  1470,   820,   690,  1470,   690}
    ,{   960,   310,   180,   960,   180}
    ,{    80,   -30,    80,  -310,  -160}
    ,{   960,   310,   180,   960,   180}
    }
   ,{{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    }
   ,{{  1560,  1560,  1430,  1280,  1430}
    ,{   -90,  -200,   -90,  -480,  -330}
    ,{   960,   310,   180,   960,   180}
    ,{  1560,  1560,  1430,  1280,  1430}
    ,{   960,   310,   180,   960,   180}
    }
   ,{{   960,   550,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   550,   550,   180,    30,   180}
    }
   }
  ,{{{  1560,  1560,  1430,  -880,  1430}
    ,{   820,   820,   690,  -880,   690}
    ,{   310,   310,   180, -1150,   180}
    ,{  1560,  1560,  1430, -1390,  1430}
    ,{   550,   550,   180, -1150,   180}
    }
   ,{{   820,   820,   690,  -880,   690}
    ,{   820,   820,   690,  -880,   690}
    ,{   310,   310,   180, -1390,   180}
    ,{   -30,   -30,  -160, -1730,  -160}
    ,{   310,   310,   180, -1390,   180}
    }
   ,{{   310,   310,   180, -1150,   180}
    ,{   310,   310,   180, -1390,   180}
    ,{   310,   310,   180, -1150,   180}
    ,{   310,   310,   180, -1390,   180}
    ,{   310,   310,   180, -1150,   180}
    }
   ,{{  1560,  1560,  1430, -1390,  1430}
    ,{  -200,  -200,  -330, -1900,  -330}
    ,{   310,   310,   180, -1390,   180}
    ,{  1560,  1560,  1430, -1390,  1430}
    ,{   310,   310,   180, -1390,   180}
    }
   ,{{   550,   550,   180, -1150,   180}
    ,{   310,   310,   180, -1390,   180}
    ,{   310,   310,   180, -1150,   180}
    ,{   310,   310,   180, -1390,   180}
    ,{   550,   550,   180, -1390,   180}
    }
   }
  ,{{{  1430,  1280,  1430,  1280,   990}
    ,{   690,   540,   690,   540,   250}
    ,{   180,    30,   180,    30,  -260}
    ,{  1430,  1280,  1430,  1280,   990}
    ,{   180,    30,   180,    30,  -260}
    }
   ,{{   690,   540,   690,   540,   250}
    ,{   690,   540,   690,   540,   250}
    ,{   180,    30,   180,    30,  -260}
    ,{    80,  -310,    80,  -310,  -360}
    ,{   180,    30,   180,    30,  -260}
    }
   ,{{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    }
   ,{{  1430,  1280,  1430,  1280,   990}
    ,{   -90,  -480,   -90,  -480,  -530}
    ,{   180,    30,   180,    30,  -260}
    ,{  1430,  1280,  1430,  1280,   990}
    ,{   180,    30,   180,    30,  -260}
    }
   ,{{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    ,{   180,    30,   180,    30,  -260}
    }
   }
  ,{{{  1470,  -580,  1430,  1470,  1430}
    ,{  1470,  -580,   690,  1470,   690}
    ,{   960,  -850,   180,   960,   180}
    ,{  1430, -1090,  1430,   960,  1430}
    ,{   960,  -850,   180,   960,   180}
    }
   ,{{  1470,  -580,   690,  1470,   690}
    ,{  1470,  -580,   690,  1470,   690}
    ,{   960, -1090,   180,   960,   180}
    ,{  -160, -1430,  -160,  -630,  -160}
    ,{   960, -1090,   180,   960,   180}
    }
   ,{{   960,  -850,   180,   960,   180}
    ,{   960, -1090,   180,   960,   180}
    ,{   960,  -850,   180,   960,   180}
    ,{   960, -1090,   180,   960,   180}
    ,{   960,  -850,   180,   960,   180}
    }
   ,{{  1430, -1090,  1430,   960,  1430}
    ,{  -330, -1600,  -330,  -800,  -330}
    ,{   960, -1090,   180,   960,   180}
    ,{  1430, -1090,  1430,  -290,  1430}
    ,{   960, -1090,   180,   960,   180}
    }
   ,{{   960,  -850,   180,   960,   180}
    ,{   960, -1090,   180,   960,   180}
    ,{   960,  -850,   180,   960,   180}
    ,{   960, -1090,   180,   960,   180}
    ,{   180, -1090,   180,  -290,   180}
    }
   }
  ,{{{  1300,  1280,  1300,  1280,   -10}
    ,{   560,   540,   560,   540,   -10}
    ,{    50,    30,    50,    30,  -760}
    ,{  1300,  1280,  1300,  1280,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{   560,   540,   560,   540,   -10}
    ,{   560,   540,   560,   540,   -10}
    ,{    50,    30,    50,    30,  -760}
    ,{   -50,  -310,   -50,  -310, -1100}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{  1300,  1280,  1300,  1280,  -760}
    ,{  -220,  -480,  -220,  -480, -1270}
    ,{    50,    30,    50,    30,  -760}
    ,{  1300,  1280,  1300,  1280,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   }
  }
 ,{{{{  2050,  1930,  1800,  2050,  1800}
    ,{  2050,  1400,  1270,  2050,  1270}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1930,  1930,  1800,  1760,  1800}
    ,{  1750,  1100,   970,  1750,   970}
    }
   ,{{  2050,  1400,  1270,  2050,  1270}
    ,{  2050,  1400,  1270,  2050,  1270}
    ,{  1740,  1090,   960,  1740,   960}
    ,{   130,    10,   130,  -260,  -110}
    ,{  1740,  1090,   960,  1740,   960}
    }
   ,{{  1760,  1110,   980,  1760,   980}
    ,{  1760,  1110,   980,  1760,   980}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1760,  1110,   980,  1760,   980}
    ,{  1750,  1100,   970,  1750,   970}
    }
   ,{{  1930,  1930,  1800,  1740,  1800}
    ,{   300,   190,   300,   -80,    60}
    ,{  1740,  1090,   960,  1740,   960}
    ,{  1930,  1930,  1800,  1650,  1800}
    ,{  1740,  1090,   960,  1740,   960}
    }
   ,{{  1760,  1110,   980,  1760,   980}
    ,{  1760,  1110,   980,  1760,   980}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1760,  1110,   980,  1760,   980}
    ,{   360,   360,     0,  -150,     0}
    }
   }
  ,{{{  1930,  1930,  1800,  -300,  1800}
    ,{  1400,  1400,  1270,  -300,  1270}
    ,{  1100,  1100,   970,  -360,   970}
    ,{  1930,  1930,  1800,  -590,  1800}
    ,{  1100,  1100,   970,  -360,   970}
    }
   ,{{  1400,  1400,  1270,  -300,  1270}
    ,{  1400,  1400,  1270,  -300,  1270}
    ,{  1090,  1090,   960,  -610,   960}
    ,{    10,    10,  -110, -1690,  -110}
    ,{  1090,  1090,   960,  -610,   960}
    }
   ,{{  1110,  1110,   980,  -360,   980}
    ,{  1110,  1110,   980,  -590,   980}
    ,{  1100,  1100,   970,  -360,   970}
    ,{  1110,  1110,   980,  -590,   980}
    ,{  1100,  1100,   970,  -360,   970}
    }
   ,{{  1930,  1930,  1800,  -610,  1800}
    ,{   190,   190,    60, -1510,    60}
    ,{  1090,  1090,   960,  -610,   960}
    ,{  1930,  1930,  1800, -1020,  1800}
    ,{  1090,  1090,   960,  -610,   960}
    }
   ,{{  1110,  1110,   980,  -360,   980}
    ,{  1110,  1110,   980,  -590,   980}
    ,{  1100,  1100,   970,  -360,   970}
    ,{  1110,  1110,   980,  -590,   980}
    ,{   360,   360,     0, -1580,     0}
    }
   }
  ,{{{  1800,  1650,  1800,  1650,  1360}
    ,{  1270,  1120,  1270,  1120,   830}
    ,{   970,   820,   970,   820,   530}
    ,{  1800,  1650,  1800,  1650,  1360}
    ,{   970,   820,   970,   820,   530}
    }
   ,{{  1270,  1120,  1270,  1120,   830}
    ,{  1270,  1120,  1270,  1120,   830}
    ,{   960,   810,   960,   810,   520}
    ,{   130,  -260,   130,  -260,  -310}
    ,{   960,   810,   960,   810,   520}
    }
   ,{{   980,   830,   980,   830,   540}
    ,{   980,   830,   980,   830,   540}
    ,{   970,   820,   970,   820,   530}
    ,{   980,   830,   980,   830,   540}
    ,{   970,   820,   970,   820,   530}
    }
   ,{{  1800,  1650,  1800,  1650,  1360}
    ,{   300,   -80,   300,   -80,  -130}
    ,{   960,   810,   960,   810,   520}
    ,{  1800,  1650,  1800,  1650,  1360}
    ,{   960,   810,   960,   810,   520}
    }
   ,{{   980,   830,   980,   830,   540}
    ,{   980,   830,   980,   830,   540}
    ,{   970,   820,   970,   820,   530}
    ,{   980,   830,   980,   830,   540}
    ,{     0,  -150,     0,  -150,  -440}
    }
   }
  ,{{{  2050,     0,  1800,  2050,  1800}
    ,{  2050,     0,  1270,  2050,  1270}
    ,{  1750,   -60,   970,  1750,   970}
    ,{  1800,  -290,  1800,  1760,  1800}
    ,{  1750,   -60,   970,  1750,   970}
    }
   ,{{  2050,     0,  1270,  2050,  1270}
    ,{  2050,     0,  1270,  2050,  1270}
    ,{  1740,  -310,   960,  1740,   960}
    ,{  -110, -1390,  -110,  -580,  -110}
    ,{  1740,  -310,   960,  1740,   960}
    }
   ,{{  1760,   -60,   980,  1760,   980}
    ,{  1760,  -290,   980,  1760,   980}
    ,{  1750,   -60,   970,  1750,   970}
    ,{  1760,  -290,   980,  1760,   980}
    ,{  1750,   -60,   970,  1750,   970}
    }
   ,{{  1800,  -310,  1800,  1740,  1800}
    ,{    60, -1210,    60,  -400,    60}
    ,{  1740,  -310,   960,  1740,   960}
    ,{  1800,  -720,  1800,    80,  1800}
    ,{  1740,  -310,   960,  1740,   960}
    }
   ,{{  1760,   -60,   980,  1760,   980}
    ,{  1760,  -290,   980,  1760,   980}
    ,{  1750,   -60,   970,  1750,   970}
    ,{  1760,  -290,   980,  1760,   980}
    ,{     0, -1280,     0,  -470,     0}
    }
   }
  ,{{{  1670,  1650,  1670,  1650,   570}
    ,{  1140,  1120,  1140,  1120,   570}
    ,{   840,   820,   840,   820,    30}
    ,{  1670,  1650,  1670,  1650,    40}
    ,{   840,   820,   840,   820,    30}
    }
   ,{{  1140,  1120,  1140,  1120,   570}
    ,{  1140,  1120,  1140,  1120,   570}
    ,{   830,   810,   830,   810,    20}
    ,{     0,  -260,     0,  -260, -1050}
    ,{   830,   810,   830,   810,    20}
    }
   ,{{   850,   830,   850,   830,    40}
    ,{   850,   830,   850,   830,    40}
    ,{   840,   820,   840,   820,    30}
    ,{   850,   830,   850,   830,    40}
    ,{   840,   820,   840,   820,    30}
    }
   ,{{  1670,  1650,  1670,  1650,    20}
    ,{   180,   -80,   180,   -80,  -870}
    ,{   830,   810,   830,   810,    20}
    ,{  1670,  1650,  1670,  1650,  -380}
    ,{   830,   810,   830,   810,    20}
    }
   ,{{   850,   830,   850,   830,    40}
    ,{   850,   830,   850,   830,    40}
    ,{   840,   820,   840,   820,    30}
    ,{   850,   830,   850,   830,    40}
    ,{  -130,  -150,  -130,  -150,  -940}
    }
   }
  }
 ,{{{{  2120,  2120,  1990,  2120,  1990}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  2120,  2120,  1990,  1990,  1990}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  1470,  1340,  2120,  1340}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{   180,    60,   180,  -210,   -60}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  2120,  1990,  1840,  1990}
    ,{  -120,  -230,  -120,  -510,  -360}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{  2120,  2120,  1990,  1840,  1990}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1550,   900,   770,  1550,   770}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{   640,   640,   270,   120,   270}
    }
   }
  ,{{{  2120,  2120,  1990,  -120,  1990}
    ,{  1470,  1470,  1340,  -230,  1340}
    ,{  1340,  1340,  1210,  -120,  1210}
    ,{  2120,  2120,  1990,  -360,  1990}
    ,{  1210,  1210,  1080,  -250,  1080}
    }
   ,{{  1470,  1470,  1340,  -230,  1340}
    ,{  1470,  1470,  1340,  -230,  1340}
    ,{  1190,  1190,  1060,  -510,  1060}
    ,{    60,    60,   -60, -1640,   -60}
    ,{  1190,  1190,  1060,  -510,  1060}
    }
   ,{{  1340,  1340,  1210,  -120,  1210}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{  1340,  1340,  1210,  -120,  1210}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{  1210,  1210,  1080,  -250,  1080}
    }
   ,{{  2120,  2120,  1990,  -510,  1990}
    ,{  -230,  -230,  -360, -1940,  -360}
    ,{  1190,  1190,  1060,  -510,  1060}
    ,{  2120,  2120,  1990,  -830,  1990}
    ,{  1190,  1190,  1060,  -510,  1060}
    }
   ,{{  1340,  1340,  1210,  -360,  1210}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{   900,   900,   770,  -560,   770}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{   640,   640,   270, -1300,   270}
    }
   }
  ,{{{  1990,  1840,  1990,  1840,  1550}
    ,{  1340,  1190,  1340,  1190,   900}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1990,  1840,  1990,  1840,  1550}
    ,{  1080,   930,  1080,   930,   640}
    }
   ,{{  1340,  1190,  1340,  1190,   900}
    ,{  1340,  1190,  1340,  1190,   900}
    ,{  1060,   910,  1060,   910,   620}
    ,{   180,  -210,   180,  -210,  -260}
    ,{  1060,   910,  1060,   910,   620}
    }
   ,{{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1080,   930,  1080,   930,   640}
    }
   ,{{  1990,  1840,  1990,  1840,  1550}
    ,{  -120,  -510,  -120,  -510,  -560}
    ,{  1060,   910,  1060,   910,   620}
    ,{  1990,  1840,  1990,  1840,  1550}
    ,{  1060,   910,  1060,   910,   620}
    }
   ,{{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{   770,   620,   770,   620,   330}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{   270,   120,   270,   120,  -170}
    }
   }
  ,{{{  2120,   180,  1990,  2120,  1990}
    ,{  2120,    60,  1340,  2120,  1340}
    ,{  1990,   180,  1210,  1990,  1210}
    ,{  1990,   -60,  1990,  1990,  1990}
    ,{  1860,    50,  1080,  1860,  1080}
    }
   ,{{  2120,    60,  1340,  2120,  1340}
    ,{  2120,    60,  1340,  2120,  1340}
    ,{  1840,  -210,  1060,  1840,  1060}
    ,{   -60, -1340,   -60,  -530,   -60}
    ,{  1840,  -210,  1060,  1840,  1060}
    }
   ,{{  1990,   180,  1210,  1990,  1210}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{  1990,   180,  1210,  1990,  1210}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{  1860,    50,  1080,  1860,  1080}
    }
   ,{{  1990,  -210,  1990,  1840,  1990}
    ,{  -360, -1640,  -360,  -830,  -360}
    ,{  1840,  -210,  1060,  1840,  1060}
    ,{  1990,  -530,  1990,   270,  1990}
    ,{  1840,  -210,  1060,  1840,  1060}
    }
   ,{{  1990,   -60,  1210,  1990,  1210}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{  1550,  -260,   770,  1550,   770}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{   270, -1000,   270,  -200,   270}
    }
   }
  ,{{{  1860,  1840,  1860,  1840,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1860,  1840,  1860,  1840,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1210,  1190,  1210,  1190,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{   930,   910,   930,   910,   120}
    ,{    50,  -210,    50,  -210, -1000}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1860,  1840,  1860,  1840,   120}
    ,{  -250,  -510,  -250,  -510, -1300}
    ,{   930,   910,   930,   910,   120}
    ,{  1860,  1840,  1860,  1840,  -200}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   640,   620,   640,   620,  -170}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   140,   120,   140,   120,  -670}
    }
   }
  }
 ,{{{{  2120,  2120,  1990,  2120,  1990}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  2120,  2120,  1990,  1990,  1990}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  1470,  1340,  2120,  1340}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{   400,   290,   400,    10,   160}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  2120,  1990,  1840,  1990}
    ,{   300,   190,   300,   -80,    60}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{  2120,  2120,  1990,  1840,  1990}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{   640,   640,   270,   120,   270}
    }
   }
  ,{{{  2120,  2120,  1990,  -120,  1990}
    ,{  1470,  1470,  1340,  -230,  1340}
    ,{  1340,  1340,  1210,  -120,  1210}
    ,{  2120,  2120,  1990,  -360,  1990}
    ,{  1210,  1210,  1080,  -250,  1080}
    }
   ,{{  1470,  1470,  1340,  -230,  1340}
    ,{  1470,  1470,  1340,  -230,  1340}
    ,{  1190,  1190,  1060,  -510,  1060}
    ,{   290,   290,   160, -1410,   160}
    ,{  1190,  1190,  1060,  -510,  1060}
    }
   ,{{  1340,  1340,  1210,  -120,  1210}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{  1340,  1340,  1210,  -120,  1210}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{  1210,  1210,  1080,  -250,  1080}
    }
   ,{{  2120,  2120,  1990,  -510,  1990}
    ,{   190,   190,    60, -1510,    60}
    ,{  1190,  1190,  1060,  -510,  1060}
    ,{  2120,  2120,  1990,  -830,  1990}
    ,{  1190,  1190,  1060,  -510,  1060}
    }
   ,{{  1340,  1340,  1210,  -360,  1210}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{  1100,  1100,   970,  -360,   970}
    ,{  1340,  1340,  1210,  -360,  1210}
    ,{   640,   640,   270, -1300,   270}
    }
   }
  ,{{{  1990,  1840,  1990,  1840,  1550}
    ,{  1340,  1190,  1340,  1190,   900}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1990,  1840,  1990,  1840,  1550}
    ,{  1080,   930,  1080,   930,   640}
    }
   ,{{  1340,  1190,  1340,  1190,   900}
    ,{  1340,  1190,  1340,  1190,   900}
    ,{  1060,   910,  1060,   910,   620}
    ,{   400,    10,   400,    10,   -40}
    ,{  1060,   910,  1060,   910,   620}
    }
   ,{{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{  1080,   930,  1080,   930,   640}
    }
   ,{{  1990,  1840,  1990,  1840,  1550}
    ,{   300,   -80,   300,   -80,  -130}
    ,{  1060,   910,  1060,   910,   620}
    ,{  1990,  1840,  1990,  1840,  1550}
    ,{  1060,   910,  1060,   910,   620}
    }
   ,{{  1210,  1060,  1210,  1060,   770}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{   970,   820,   970,   820,   530}
    ,{  1210,  1060,  1210,  1060,   770}
    ,{   270,   120,   270,   120,  -170}
    }
   }
  ,{{{  2120,   180,  1990,  2120,  1990}
    ,{  2120,    60,  1340,  2120,  1340}
    ,{  1990,   180,  1210,  1990,  1210}
    ,{  1990,   -60,  1990,  1990,  1990}
    ,{  1860,    50,  1080,  1860,  1080}
    }
   ,{{  2120,    60,  1340,  2120,  1340}
    ,{  2120,    60,  1340,  2120,  1340}
    ,{  1840,  -210,  1060,  1840,  1060}
    ,{   160, -1110,   160,  -310,   160}
    ,{  1840,  -210,  1060,  1840,  1060}
    }
   ,{{  1990,   180,  1210,  1990,  1210}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{  1990,   180,  1210,  1990,  1210}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{  1860,    50,  1080,  1860,  1080}
    }
   ,{{  1990,  -210,  1990,  1840,  1990}
    ,{    60, -1210,    60,  -400,    60}
    ,{  1840,  -210,  1060,  1840,  1060}
    ,{  1990,  -530,  1990,   270,  1990}
    ,{  1840,  -210,  1060,  1840,  1060}
    }
   ,{{  1990,   -60,  1210,  1990,  1210}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{  1750,   -60,   970,  1750,   970}
    ,{  1990,   -60,  1210,  1990,  1210}
    ,{   270, -1000,   270,  -200,   270}
    }
   }
  ,{{{  1860,  1840,  1860,  1840,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1860,  1840,  1860,  1840,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1210,  1190,  1210,  1190,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{   930,   910,   930,   910,   120}
    ,{   270,    10,   270,    10,  -780}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1860,  1840,  1860,  1840,   120}
    ,{   180,   -80,   180,   -80,  -870}
    ,{   930,   910,   930,   910,   120}
    ,{  1860,  1840,  1860,  1840,  -200}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   840,   820,   840,   820,    30}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   140,   120,   140,   120,  -670}
    }
   }
  }
 }
,{{{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  ,{{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   ,{{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    ,{   INF,   INF,   INF,   INF,   INF}
    }
   }
  }
 ,{{{{  1350,   850,   720,  1350,   720}
    ,{  1300,   650,   540,  1300,   520}
    ,{  1350,   700,   570,  1350,   570}
    ,{  1300,   850,   720,  1300,   720}
    ,{  1250,   590,   460,  1250,   460}
    }
   ,{{  1160,   500,   400,  1160,   370}
    ,{  1160,   500,   370,  1160,   370}
    ,{   850,   190,    60,   850,    60}
    ,{   400,   290,   400,    10,   170}
    ,{   850,   190,    60,   850,    60}
    }
   ,{{  1300,   650,   520,  1300,   520}
    ,{  1300,   650,   520,  1300,   520}
    ,{  1290,   640,   510,  1290,   510}
    ,{  1300,   650,   520,  1300,   520}
    ,{  1250,   590,   460,  1250,   460}
    }
   ,{{   850,   850,   720,   850,   720}
    ,{   540,     0,   540,  -270,  -120}
    ,{   850,   190,    60,   850,    60}
    ,{   850,   850,   720,   570,   720}
    ,{   850,   190,    60,   850,    60}
    }
   ,{{  1350,   700,   570,  1350,   570}
    ,{  1300,   650,   520,  1300,   520}
    ,{  1350,   700,   570,  1350,   570}
    ,{  1300,   650,   520,  1300,   520}
    ,{   100,   100,  -270,  -230,  -270}
    }
   }
  ,{{{   850,   850,   720,  -330,   720}
    ,{   650,   650,   520,  -620,   520}
    ,{   700,   700,   570,  -330,   570}
    ,{   850,   850,   720,  -620,   720}
    ,{   590,   590,   460,  -440,   460}
    }
   ,{{   500,   500,   370,  -770,   370}
    ,{   500,   500,   370,  -770,   370}
    ,{   190,   190,    60, -1070,    60}
    ,{   290,   290,   160,  -980,   160}
    ,{   190,   190,    60, -1080,    60}
    }
   ,{{   650,   650,   520,  -390,   520}
    ,{   650,   650,   520,  -620,   520}
    ,{   640,   640,   510,  -390,   510}
    ,{   650,   650,   520,  -620,   520}
    ,{   590,   590,   460,  -440,   460}
    }
   ,{{   850,   850,   720, -1080,   720}
    ,{    10,     0,    10, -1270,  -120}
    ,{   190,   190,    60, -1080,    60}
    ,{   850,   850,   720, -1080,   720}
    ,{   190,   190,    60, -1080,    60}
    }
   ,{{   700,   700,   570,  -330,   570}
    ,{   650,   650,   520,  -620,   520}
    ,{   700,   700,   570,  -330,   570}
    ,{   650,   650,   520,  -620,   520}
    ,{   100,   100,  -270, -1300,  -270}
    }
   }
  ,{{{   720,   570,   720,   570,   480}
    ,{   540,   370,   540,   370,   280}
    ,{   570,   420,   570,   420,   340}
    ,{   720,   570,   720,   570,   480}
    ,{   460,   310,   460,   310,   230}
    }
   ,{{   400,   220,   400,   220,   170}
    ,{   370,   220,   370,   220,   140}
    ,{    60,   -80,    60,   -80,  -170}
    ,{   400,    10,   400,    10,   170}
    ,{    60,   -80,    60,   -80,  -170}
    }
   ,{{   520,   370,   520,   370,   280}
    ,{   520,   370,   520,   370,   280}
    ,{   510,   360,   510,   360,   280}
    ,{   520,   370,   520,   370,   280}
    ,{   460,   310,   460,   310,   230}
    }
   ,{{   720,   570,   720,   570,   480}
    ,{   540,  -100,   540,  -270,  -120}
    ,{    60,   -80,    60,   -80,  -170}
    ,{   720,   570,   720,   570,   480}
    ,{    60,   -80,    60,   -80,  -170}
    }
   ,{{   570,   420,   570,   420,   340}
    ,{   520,   370,   520,   370,   280}
    ,{   570,   420,   570,   420,   340}
    ,{   520,   370,   520,   370,   280}
    ,{  -270,  -420,  -270,  -420,  -500}
    }
   }
  ,{{{  1350,  -230,   720,  1350,   720}
    ,{  1300,  -530,   520,  1300,   520}
    ,{  1350,  -230,   570,  1350,   570}
    ,{  1300,  -530,   720,  1300,   720}
    ,{  1250,  -340,   460,  1250,   460}
    }
   ,{{  1160,  -670,   370,  1160,   370}
    ,{  1160,  -670,   370,  1160,   370}
    ,{   850,  -980,    60,   850,    60}
    ,{   160,  -890,   160,  -310,   160}
    ,{   850,  -980,    60,   850,    60}
    }
   ,{{  1300,  -290,   520,  1300,   520}
    ,{  1300,  -530,   520,  1300,   520}
    ,{  1290,  -290,   510,  1290,   510}
    ,{  1300,  -530,   520,  1300,   520}
    ,{  1250,  -340,   460,  1250,   460}
    }
   ,{{   850,  -980,   720,   850,   720}
    ,{  -120, -1170,  -120,  -590,  -120}
    ,{   850,  -980,    60,   850,    60}
    ,{   720, -1580,   720, -1000,   720}
    ,{   850,  -980,    60,   850,    60}
    }
   ,{{  1350,  -230,   570,  1350,   570}
    ,{  1300,  -530,   520,  1300,   520}
    ,{  1350,  -230,   570,  1350,   570}
    ,{  1300,  -530,   520,  1300,   520}
    ,{  -230, -1320,  -270,  -230,  -270}
    }
   }
  ,{{{   590,   570,   590,   570,   -90}
    ,{   390,   370,   390,   370,   -90}
    ,{   440,   420,   440,   420,  -360}
    ,{   590,   570,   590,   570,  -420}
    ,{   330,   310,   330,   310,  -470}
    }
   ,{{   270,   220,   270,   220,  -320}
    ,{   240,   220,   240,   220,  -320}
    ,{   -60,   -80,   -60,   -80,  -830}
    ,{   270,    10,   270,    10,  -780}
    ,{   -60,   -80,   -60,   -80,  -870}
    }
   ,{{   390,   370,   390,   370,   -90}
    ,{   390,   370,   390,   370,   -90}
    ,{   380,   360,   380,   360,  -420}
    ,{   390,   370,   390,   370,  -420}
    ,{   330,   310,   330,   310,  -470}
    }
   ,{{   590,   570,   590,   570,  -810}
    ,{   -10,  -270,   -10,  -270,  -810}
    ,{   -60,   -80,   -60,   -80,  -870}
    ,{   590,   570,   590,   570, -1470}
    ,{   -60,   -80,   -60,   -80,  -870}
    }
   ,{{   440,   420,   440,   420,  -360}
    ,{   390,   370,   390,   370,  -420}
    ,{   440,   420,   440,   420,  -360}
    ,{   390,   370,   390,   370,  -420}
    ,{  -400,  -420,  -400,  -420, -1210}
    }
   }
  }
 ,{{{{  1320,   850,   720,  1320,   720}
    ,{  1320,   670,   540,  1320,   540}
    ,{   870,   220,    90,   870,    90}
    ,{   960,   850,   720,   960,   720}
    ,{   870,   250,    90,   870,    90}
    }
   ,{{  1320,   670,   540,  1320,   540}
    ,{  1320,   670,   540,  1320,   540}
    ,{   870,   220,    90,   870,    90}
    ,{  -410,  -520,  -410,  -800,  -640}
    ,{   870,   220,    90,   870,    90}
    }
   ,{{   960,   300,   170,   960,   170}
    ,{   960,   300,   170,   960,   170}
    ,{   650,     0,  -130,   650,  -130}
    ,{   960,   300,   170,   960,   170}
    ,{   650,     0,  -130,   650,  -130}
    }
   ,{{   870,   850,   720,   870,   720}
    ,{    70,   -40,    70,  -320,  -170}
    ,{   870,   220,    90,   870,    90}
    ,{   850,   850,   720,   570,   720}
    ,{   870,   220,    90,   870,    90}
    }
   ,{{   960,   300,   170,   960,   170}
    ,{   960,   300,   170,   960,   170}
    ,{   340,  -310,  -440,   340,  -440}
    ,{   960,   300,   170,   960,   170}
    ,{   250,   250,   -90,  -260,  -110}
    }
   }
  ,{{{   850,   850,   720,   540,   720}
    ,{   670,   670,   540,    10,   540}
    ,{   540,   220,    90,   540,    90}
    ,{   850,   850,   720,  -970,   720}
    ,{   250,   250,    90,  -810,    90}
    }
   ,{{   670,   670,   540,  -100,   540}
    ,{   670,   670,   540,  -600,   540}
    ,{   220,   220,    90,  -100,    90}
    ,{  -520,  -520,  -650, -1790,  -650}
    ,{   220,   220,    90, -1050,    90}
    }
   ,{{   540,   300,   170,   540,   170}
    ,{   300,   300,   170,    10,   170}
    ,{   540,     0,  -130,   540,  -130}
    ,{   300,   300,   170,  -970,   170}
    ,{     0,     0,  -130, -1030,  -130}
    }
   ,{{   850,   850,   720, -1050,   720}
    ,{   -40,   -40,  -170, -1320,  -170}
    ,{   220,   220,    90, -1050,    90}
    ,{   850,   850,   720, -1680,   720}
    ,{   220,   220,    90, -1050,    90}
    }
   ,{{   300,   300,   170,  -810,   170}
    ,{   300,   300,   170,  -970,   170}
    ,{  -310,  -310,  -440, -1340,  -440}
    ,{   300,   300,   170,  -970,   170}
    ,{   250,   250,   -90,  -810,  -110}
    }
   }
  ,{{{   720,   570,   720,   570,   480}
    ,{   540,   390,   540,   390,   300}
    ,{    90,   -60,    90,   -60,  -140}
    ,{   720,   570,   720,   570,   480}
    ,{    90,   -60,    90,   -60,  -140}
    }
   ,{{   540,   390,   540,   390,   300}
    ,{   540,   390,   540,   390,   300}
    ,{    90,   -60,    90,   -60,  -140}
    ,{  -410,  -800,  -410,  -800,  -640}
    ,{    90,   -60,    90,   -60,  -140}
    }
   ,{{   170,    20,   170,    20,   -60}
    ,{   170,    20,   170,    20,   -60}
    ,{  -130,  -280,  -130,  -280,  -360}
    ,{   170,    20,   170,    20,   -60}
    ,{  -130,  -280,  -130,  -280,  -360}
    }
   ,{{   720,   570,   720,   570,   480}
    ,{    70,  -320,    70,  -320,  -170}
    ,{    90,   -60,    90,   -60,  -140}
    ,{   720,   570,   720,   570,   480}
    ,{    90,   -60,    90,   -60,  -140}
    }
   ,{{   170,    20,   170,    20,   -60}
    ,{   170,    20,   170,    20,   -60}
    ,{  -440,  -590,  -440,  -590,  -670}
    ,{   170,    20,   170,    20,   -60}
    ,{  -110,  -260,  -110,  -260,  -350}
    }
   }
  ,{{{  1320,  -350,   720,  1320,   720}
    ,{  1320,  -730,   540,  1320,   540}
    ,{   870,  -350,    90,   870,    90}
    ,{   960,  -870,   720,   960,   720}
    ,{   870,  -940,    90,   870,    90}
    }
   ,{{  1320,  -350,   540,  1320,   540}
    ,{  1320,  -730,   540,  1320,   540}
    ,{   870,  -350,    90,   870,    90}
    ,{  -650, -1920,  -650, -1120,  -650}
    ,{   870,  -960,    90,   870,    90}
    }
   ,{{   960,  -870,   170,   960,   170}
    ,{   960, -1100,   170,   960,   170}
    ,{   650,  -940,  -130,   650,  -130}
    ,{   960,  -870,   170,   960,   170}
    ,{   650,  -940,  -130,   650,  -130}
    }
   ,{{   870,  -960,   720,   870,   720}
    ,{  -170, -1450,  -170,  -640,  -170}
    ,{   870,  -960,    90,   870,    90}
    ,{   720, -1370,   720, -1000,   720}
    ,{   870,  -960,    90,   870,    90}
    }
   ,{{   960,  -870,   170,   960,   170}
    ,{   960,  -870,   170,   960,   170}
    ,{   340, -1250,  -440,   340,  -440}
    ,{   960,  -870,   170,   960,   170}
    ,{  -110, -1360,  -110,  -580,  -110}
    }
   }
  ,{{{   590,   570,   590,   570,  -160}
    ,{   410,   390,   410,   390,  -160}
    ,{   -40,   -60,   -40,   -60,  -850}
    ,{   590,   570,   590,   570,  -230}
    ,{   -40,   -60,   -40,   -60,  -850}
    }
   ,{{   410,   390,   410,   390,  -160}
    ,{   410,   390,   410,   390,  -160}
    ,{   -40,   -60,   -40,   -60,  -850}
    ,{  -540,  -800,  -540,  -800, -1520}
    ,{   -40,   -60,   -40,   -60,  -850}
    }
   ,{{    40,    20,    40,    20,  -400}
    ,{    40,    20,    40,    20,  -400}
    ,{  -260,  -280,  -260,  -280, -1070}
    ,{    40,    20,    40,    20,  -760}
    ,{  -260,  -280,  -260,  -280, -1070}
    }
   ,{{   590,   570,   590,   570,  -230}
    ,{   -60,  -320,   -60,  -320, -1110}
    ,{   -40,   -60,   -40,   -60,  -850}
    ,{   590,   570,   590,   570,  -230}
    ,{   -40,   -60,   -40,   -60,  -850}
    }
   ,{{    40,    20,    40,    20,  -760}
    ,{    40,    20,    40,    20,  -760}
    ,{  -570,  -590,  -570,  -590, -1380}
    ,{    40,    20,    40,    20,  -760}
    ,{  -240,  -260,  -240,  -260, -1050}
    }
   }
  }
 ,{{{{  1010,  1010,   880,   730,   880}
    ,{   410,   -30,    40,   410,  -190}
    ,{   410,  -240,  -370,   410,  -370}
    ,{  1010,  1010,   880,   730,   880}
    ,{   410,     0,  -370,   410,  -370}
    }
   ,{{   410,   -70,  -150,   410,  -370}
    ,{   230,   -70,  -550,   230,  -550}
    ,{   410,  -240,  -370,   410,  -370}
    ,{  -150,  -260,  -150,  -540,  -380}
    ,{   410,  -240,  -370,   410,  -370}
    }
   ,{{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    }
   ,{{  1010,  1010,   880,   730,   880}
    ,{    40,   -30,    40,  -350,  -190}
    ,{   410,  -240,  -370,   410,  -370}
    ,{  1010,  1010,   880,   730,   880}
    ,{   410,  -240,  -370,   410,  -370}
    }
   ,{{   410,     0,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{   410,  -240,  -370,   410,  -370}
    ,{     0,     0,  -370,  -520,  -370}
    }
   }
  ,{{{  1010,  1010,   880, -1280,   880}
    ,{   -30,   -30,  -200, -1340,  -200}
    ,{  -240,  -240,  -370, -1280,  -370}
    ,{  1010,  1010,   880, -1520,   880}
    ,{     0,     0,  -370, -1280,  -370}
    }
   ,{{   -70,   -70,  -370, -1520,  -370}
    ,{   -70,   -70,  -550, -1700,  -550}
    ,{  -240,  -240,  -370, -1520,  -370}
    ,{  -260,  -260,  -390, -1530,  -390}
    ,{  -240,  -240,  -370, -1520,  -370}
    }
   ,{{  -240,  -240,  -370, -1280,  -370}
    ,{  -240,  -240,  -370, -1520,  -370}
    ,{  -240,  -240,  -370, -1280,  -370}
    ,{  -240,  -240,  -370, -1520,  -370}
    ,{  -240,  -240,  -370, -1280,  -370}
    }
   ,{{  1010,  1010,   880, -1340,   880}
    ,{   -30,   -30,  -200, -1340,  -200}
    ,{  -240,  -240,  -370, -1520,  -370}
    ,{  1010,  1010,   880, -1520,   880}
    ,{  -240,  -240,  -370, -1520,  -370}
    }
   ,{{     0,     0,  -370, -1280,  -370}
    ,{  -240,  -240,  -370, -1520,  -370}
    ,{  -240,  -240,  -370, -1280,  -370}
    ,{  -240,  -240,  -370, -1520,  -370}
    ,{     0,     0,  -370, -1520,  -370}
    }
   }
  ,{{{   880,   730,   880,   730,   640}
    ,{    40,  -350,    40,  -350,  -190}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{   880,   730,   880,   730,   640}
    ,{  -370,  -520,  -370,  -520,  -610}
    }
   ,{{  -150,  -520,  -150,  -520,  -380}
    ,{  -550,  -700,  -550,  -700,  -790}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{  -150,  -540,  -150,  -540,  -380}
    ,{  -370,  -520,  -370,  -520,  -610}
    }
   ,{{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    }
   ,{{   880,   730,   880,   730,   640}
    ,{    40,  -350,    40,  -350,  -190}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{   880,   730,   880,   730,   640}
    ,{  -370,  -520,  -370,  -520,  -610}
    }
   ,{{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    ,{  -370,  -520,  -370,  -520,  -610}
    }
   }
  ,{{{   880, -1180,   880,   410,   880}
    ,{   410, -1250,  -200,   410,  -200}
    ,{   410, -1180,  -370,   410,  -370}
    ,{   880, -1420,   880,   410,   880}
    ,{   410, -1180,  -370,   410,  -370}
    }
   ,{{   410, -1420,  -370,   410,  -370}
    ,{   230, -1600,  -550,   230,  -550}
    ,{   410, -1420,  -370,   410,  -370}
    ,{  -390, -1440,  -390,  -860,  -390}
    ,{   410, -1420,  -370,   410,  -370}
    }
   ,{{   410, -1180,  -370,   410,  -370}
    ,{   410, -1420,  -370,   410,  -370}
    ,{   410, -1180,  -370,   410,  -370}
    ,{   410, -1420,  -370,   410,  -370}
    ,{   410, -1180,  -370,   410,  -370}
    }
   ,{{   880, -1250,   880,   410,   880}
    ,{  -200, -1250,  -200,  -670,  -200}
    ,{   410, -1420,  -370,   410,  -370}
    ,{   880, -1420,   880,  -840,   880}
    ,{   410, -1420,  -370,   410,  -370}
    }
   ,{{   410, -1180,  -370,   410,  -370}
    ,{   410, -1420,  -370,   410,  -370}
    ,{   410, -1180,  -370,   410,  -370}
    ,{   410, -1420,  -370,   410,  -370}
    ,{  -370, -1420,  -370,  -840,  -370}
    }
   }
  ,{{{   750,   730,   750,   730, -1140}
    ,{   -90,  -350,   -90,  -350, -1140}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{   750,   730,   750,   730, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{  -280,  -520,  -280,  -520, -1250}
    ,{  -680,  -700,  -680,  -700, -1250}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -280,  -540,  -280,  -540, -1330}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{   750,   730,   750,   730, -1140}
    ,{   -90,  -350,   -90,  -350, -1140}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{   750,   730,   750,   730, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   ,{{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    ,{  -500,  -520,  -500,  -520, -1310}
    }
   }
  }
 ,{{{{  1560,  1560,  1430,  1470,  1430}
    ,{  1470,   820,   690,  1470,   690}
    ,{   960,   310,   180,   960,   180}
    ,{  1560,  1560,  1430,  1280,  1430}
    ,{   960,   550,   180,   960,   180}
    }
   ,{{  1470,   820,   690,  1470,   690}
    ,{  1470,   820,   690,  1470,   690}
    ,{   960,   310,   180,   960,   180}
    ,{    80,   -30,    80,  -310,  -150}
    ,{   960,   310,   180,   960,   180}
    }
   ,{{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    }
   ,{{  1560,  1560,  1430,  1280,  1430}
    ,{   -90,  -200,   -90,  -480,  -320}
    ,{   960,   310,   180,   960,   180}
    ,{  1560,  1560,  1430,  1280,  1430}
    ,{   960,   310,   180,   960,   180}
    }
   ,{{   960,   550,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   960,   310,   180,   960,   180}
    ,{   550,   550,   180,    30,   180}
    }
   }
  ,{{{  1560,  1560,  1430,   -30,  1430}
    ,{   820,   820,   690,   -30,   690}
    ,{   310,   310,   180,  -720,   180}
    ,{  1560,  1560,  1430,  -960,  1430}
    ,{   550,   550,   180,  -720,   180}
    }
   ,{{   820,   820,   690,   -30,   690}
    ,{   820,   820,   690,   -30,   690}
    ,{   310,   310,   180,  -960,   180}
    ,{   -30,   -30,  -160, -1300,  -160}
    ,{   310,   310,   180,  -960,   180}
    }
   ,{{   310,   310,   180,  -720,   180}
    ,{   310,   310,   180,  -960,   180}
    ,{   310,   310,   180,  -720,   180}
    ,{   310,   310,   180,  -960,   180}
    ,{   310,   310,   180,  -720,   180}
    }
   ,{{  1560,  1560,  1430,  -960,  1430}
    ,{  -200,  -200,  -330, -1470,  -330}
    ,{   310,   310,   180,  -960,   180}
    ,{  1560,  1560,  1430,  -960,  1430}
    ,{   310,   310,   180,  -960,   180}
    }
   ,{{   550,   550,   180,  -720,   180}
    ,{   310,   310,   180,  -960,   180}
    ,{   310,   310,   180,  -720,   180}
    ,{   310,   310,   180,  -960,   180}
    ,{   550,   550,   180,  -960,   180}
    }
   }
  ,{{{  1430,  1280,  1430,  1280,  1200}
    ,{   690,   540,   690,   540,   450}
    ,{   180,    30,   180,    30,   -50}
    ,{  1430,  1280,  1430,  1280,  1200}
    ,{   180,    30,   180,    30,   -50}
    }
   ,{{   690,   540,   690,   540,   450}
    ,{   690,   540,   690,   540,   450}
    ,{   180,    30,   180,    30,   -50}
    ,{    80,  -310,    80,  -310,  -150}
    ,{   180,    30,   180,    30,   -50}
    }
   ,{{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    }
   ,{{  1430,  1280,  1430,  1280,  1200}
    ,{   -90,  -480,   -90,  -480,  -320}
    ,{   180,    30,   180,    30,   -50}
    ,{  1430,  1280,  1430,  1280,  1200}
    ,{   180,    30,   180,    30,   -50}
    }
   ,{{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    ,{   180,    30,   180,    30,   -50}
    }
   }
  ,{{{  1470,  -360,  1430,  1470,  1430}
    ,{  1470,  -360,   690,  1470,   690}
    ,{   960,  -630,   180,   960,   180}
    ,{  1430,  -870,  1430,   960,  1430}
    ,{   960,  -630,   180,   960,   180}
    }
   ,{{  1470,  -360,   690,  1470,   690}
    ,{  1470,  -360,   690,  1470,   690}
    ,{   960,  -870,   180,   960,   180}
    ,{  -160, -1210,  -160,  -630,  -160}
    ,{   960,  -870,   180,   960,   180}
    }
   ,{{   960,  -630,   180,   960,   180}
    ,{   960,  -870,   180,   960,   180}
    ,{   960,  -630,   180,   960,   180}
    ,{   960,  -870,   180,   960,   180}
    ,{   960,  -630,   180,   960,   180}
    }
   ,{{  1430,  -870,  1430,   960,  1430}
    ,{  -330, -1380,  -330,  -800,  -330}
    ,{   960,  -870,   180,   960,   180}
    ,{  1430,  -870,  1430,  -290,  1430}
    ,{   960,  -870,   180,   960,   180}
    }
   ,{{   960,  -630,   180,   960,   180}
    ,{   960,  -870,   180,   960,   180}
    ,{   960,  -630,   180,   960,   180}
    ,{   960,  -870,   180,   960,   180}
    ,{   180,  -870,   180,  -290,   180}
    }
   }
  ,{{{  1300,  1280,  1300,  1280,   -10}
    ,{   560,   540,   560,   540,   -10}
    ,{    50,    30,    50,    30,  -760}
    ,{  1300,  1280,  1300,  1280,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{   560,   540,   560,   540,   -10}
    ,{   560,   540,   560,   540,   -10}
    ,{    50,    30,    50,    30,  -760}
    ,{   -50,  -310,   -50,  -310, -1100}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{  1300,  1280,  1300,  1280,  -760}
    ,{  -220,  -480,  -220,  -480, -1270}
    ,{    50,    30,    50,    30,  -760}
    ,{  1300,  1280,  1300,  1280,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   ,{{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    ,{    50,    30,    50,    30,  -760}
    }
   }
  }
 ,{{{{  2050,  1930,  1800,  2050,  1800}
    ,{  2050,  1400,  1270,  2050,  1270}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1930,  1930,  1800,  1760,  1800}
    ,{  1750,  1100,   970,  1750,   970}
    }
   ,{{  2050,  1400,  1270,  2050,  1270}
    ,{  2050,  1400,  1270,  2050,  1270}
    ,{  1740,  1090,   960,  1740,   960}
    ,{   130,    10,   130,  -260,  -110}
    ,{  1740,  1090,   960,  1740,   960}
    }
   ,{{  1760,  1110,   980,  1760,   980}
    ,{  1760,  1110,   980,  1760,   980}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1760,  1110,   980,  1760,   980}
    ,{  1750,  1100,   970,  1750,   970}
    }
   ,{{  1930,  1930,  1800,  1740,  1800}
    ,{   300,   190,   300,   -80,    70}
    ,{  1740,  1090,   960,  1740,   960}
    ,{  1930,  1930,  1800,  1650,  1800}
    ,{  1740,  1090,   960,  1740,   960}
    }
   ,{{  1760,  1110,   980,  1760,   980}
    ,{  1760,  1110,   980,  1760,   980}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1760,  1110,   980,  1760,   980}
    ,{   360,   360,     0,  -150,     0}
    }
   }
  ,{{{  1930,  1930,  1800,   130,  1800}
    ,{  1400,  1400,  1270,   130,  1270}
    ,{  1100,  1100,   970,    70,   970}
    ,{  1930,  1930,  1800,  -160,  1800}
    ,{  1100,  1100,   970,    70,   970}
    }
   ,{{  1400,  1400,  1270,   130,  1270}
    ,{  1400,  1400,  1270,   130,  1270}
    ,{  1090,  1090,   960,  -180,   960}
    ,{    10,    10,  -110, -1260,  -110}
    ,{  1090,  1090,   960,  -180,   960}
    }
   ,{{  1110,  1110,   980,    70,   980}
    ,{  1110,  1110,   980,  -160,   980}
    ,{  1100,  1100,   970,    70,   970}
    ,{  1110,  1110,   980,  -160,   980}
    ,{  1100,  1100,   970,    70,   970}
    }
   ,{{  1930,  1930,  1800,  -180,  1800}
    ,{   190,   190,    60, -1080,    60}
    ,{  1090,  1090,   960,  -180,   960}
    ,{  1930,  1930,  1800,  -590,  1800}
    ,{  1090,  1090,   960,  -180,   960}
    }
   ,{{  1110,  1110,   980,    70,   980}
    ,{  1110,  1110,   980,  -160,   980}
    ,{  1100,  1100,   970,    70,   970}
    ,{  1110,  1110,   980,  -160,   980}
    ,{   360,   360,     0, -1150,     0}
    }
   }
  ,{{{  1800,  1650,  1800,  1650,  1570}
    ,{  1270,  1120,  1270,  1120,  1040}
    ,{   970,   820,   970,   820,   740}
    ,{  1800,  1650,  1800,  1650,  1570}
    ,{   970,   820,   970,   820,   740}
    }
   ,{{  1270,  1120,  1270,  1120,  1040}
    ,{  1270,  1120,  1270,  1120,  1040}
    ,{   960,   810,   960,   810,   730}
    ,{   130,  -260,   130,  -260,  -110}
    ,{   960,   810,   960,   810,   730}
    }
   ,{{   980,   830,   980,   830,   740}
    ,{   980,   830,   980,   830,   740}
    ,{   970,   820,   970,   820,   740}
    ,{   980,   830,   980,   830,   740}
    ,{   970,   820,   970,   820,   740}
    }
   ,{{  1800,  1650,  1800,  1650,  1570}
    ,{   300,   -80,   300,   -80,    70}
    ,{   960,   810,   960,   810,   730}
    ,{  1800,  1650,  1800,  1650,  1570}
    ,{   960,   810,   960,   810,   730}
    }
   ,{{   980,   830,   980,   830,   740}
    ,{   980,   830,   980,   830,   740}
    ,{   970,   820,   970,   820,   740}
    ,{   980,   830,   980,   830,   740}
    ,{     0,  -150,     0,  -150,  -240}
    }
   }
  ,{{{  2050,   220,  1800,  2050,  1800}
    ,{  2050,   220,  1270,  2050,  1270}
    ,{  1750,   170,   970,  1750,   970}
    ,{  1800,   -70,  1800,  1760,  1800}
    ,{  1750,   170,   970,  1750,   970}
    }
   ,{{  2050,   220,  1270,  2050,  1270}
    ,{  2050,   220,  1270,  2050,  1270}
    ,{  1740,   -80,   960,  1740,   960}
    ,{  -110, -1160,  -110,  -580,  -110}
    ,{  1740,   -80,   960,  1740,   960}
    }
   ,{{  1760,   170,   980,  1760,   980}
    ,{  1760,   -70,   980,  1760,   980}
    ,{  1750,   170,   970,  1750,   970}
    ,{  1760,   -70,   980,  1760,   980}
    ,{  1750,   170,   970,  1750,   970}
    }
   ,{{  1800,   -80,  1800,  1740,  1800}
    ,{    60,  -980,    60,  -400,    60}
    ,{  1740,   -80,   960,  1740,   960}
    ,{  1800,  -490,  1800,    80,  1800}
    ,{  1740,   -80,   960,  1740,   960}
    }
   ,{{  1760,   170,   980,  1760,   980}
    ,{  1760,   -70,   980,  1760,   980}
    ,{  1750,   170,   970,  1750,   970}
    ,{  1760,   -70,   980,  1760,   980}
    ,{     0, -1050,     0,  -470,     0}
    }
   }
  ,{{{  1670,  1650,  1670,  1650,   570}
    ,{  1140,  1120,  1140,  1120,   570}
    ,{   840,   820,   840,   820,    30}
    ,{  1670,  1650,  1670,  1650,    40}
    ,{   840,   820,   840,   820,    30}
    }
   ,{{  1140,  1120,  1140,  1120,   570}
    ,{  1140,  1120,  1140,  1120,   570}
    ,{   830,   810,   830,   810,    20}
    ,{     0,  -260,     0,  -260, -1050}
    ,{   830,   810,   830,   810,    20}
    }
   ,{{   850,   830,   850,   830,    40}
    ,{   850,   830,   850,   830,    40}
    ,{   840,   820,   840,   820,    30}
    ,{   850,   830,   850,   830,    40}
    ,{   840,   820,   840,   820,    30}
    }
   ,{{  1670,  1650,  1670,  1650,    20}
    ,{   180,   -80,   180,   -80,  -870}
    ,{   830,   810,   830,   810,    20}
    ,{  1670,  1650,  1670,  1650,  -380}
    ,{   830,   810,   830,   810,    20}
    }
   ,{{   850,   830,   850,   830,    40}
    ,{   850,   830,   850,   830,    40}
    ,{   840,   820,   840,   820,    30}
    ,{   850,   830,   850,   830,    40}
    ,{  -130,  -150,  -130,  -150,  -940}
    }
   }
  }
 ,{{{{  2120,  2120,  1990,  2120,  1990}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  2120,  2120,  1990,  1990,  1990}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  1470,  1340,  2120,  1340}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{   180,    60,   180,  -210,   -60}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  2120,  1990,  1840,  1990}
    ,{  -120,  -230,  -120,  -510,  -360}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{  2120,  2120,  1990,  1840,  1990}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1550,   900,   770,  1550,   770}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{   640,   640,   270,   120,   270}
    }
   }
  ,{{{  2120,  2120,  1990,   300,  1990}
    ,{  1470,  1470,  1340,   190,  1340}
    ,{  1340,  1340,  1210,   300,  1210}
    ,{  2120,  2120,  1990,    60,  1990}
    ,{  1210,  1210,  1080,   180,  1080}
    }
   ,{{  1470,  1470,  1340,   190,  1340}
    ,{  1470,  1470,  1340,   190,  1340}
    ,{  1190,  1190,  1060,   -80,  1060}
    ,{    60,    60,   -60, -1210,   -60}
    ,{  1190,  1190,  1060,   -80,  1060}
    }
   ,{{  1340,  1340,  1210,   300,  1210}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{  1340,  1340,  1210,   300,  1210}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{  1210,  1210,  1080,   180,  1080}
    }
   ,{{  2120,  2120,  1990,   -80,  1990}
    ,{  -230,  -230,  -360, -1510,  -360}
    ,{  1190,  1190,  1060,   -80,  1060}
    ,{  2120,  2120,  1990,  -400,  1990}
    ,{  1190,  1190,  1060,   -80,  1060}
    }
   ,{{  1340,  1340,  1210,    60,  1210}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{   900,   900,   770,  -130,   770}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{   640,   640,   270,  -870,   270}
    }
   }
  ,{{{  1990,  1840,  1990,  1840,  1750}
    ,{  1340,  1190,  1340,  1190,  1100}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1990,  1840,  1990,  1840,  1750}
    ,{  1080,   930,  1080,   930,   840}
    }
   ,{{  1340,  1190,  1340,  1190,  1100}
    ,{  1340,  1190,  1340,  1190,  1100}
    ,{  1060,   910,  1060,   910,   820}
    ,{   180,  -210,   180,  -210,   -60}
    ,{  1060,   910,  1060,   910,   820}
    }
   ,{{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1080,   930,  1080,   930,   840}
    }
   ,{{  1990,  1840,  1990,  1840,  1750}
    ,{  -120,  -510,  -120,  -510,  -360}
    ,{  1060,   910,  1060,   910,   820}
    ,{  1990,  1840,  1990,  1840,  1750}
    ,{  1060,   910,  1060,   910,   820}
    }
   ,{{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{   770,   620,   770,   620,   530}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{   270,   120,   270,   120,    30}
    }
   }
  ,{{{  2120,   400,  1990,  2120,  1990}
    ,{  2120,   290,  1340,  2120,  1340}
    ,{  1990,   400,  1210,  1990,  1210}
    ,{  1990,   160,  1990,  1990,  1990}
    ,{  1860,   270,  1080,  1860,  1080}
    }
   ,{{  2120,   290,  1340,  2120,  1340}
    ,{  2120,   290,  1340,  2120,  1340}
    ,{  1840,    10,  1060,  1840,  1060}
    ,{   -60, -1110,   -60,  -530,   -60}
    ,{  1840,    10,  1060,  1840,  1060}
    }
   ,{{  1990,   400,  1210,  1990,  1210}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{  1990,   400,  1210,  1990,  1210}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{  1860,   270,  1080,  1860,  1080}
    }
   ,{{  1990,    10,  1990,  1840,  1990}
    ,{  -360, -1410,  -360,  -830,  -360}
    ,{  1840,    10,  1060,  1840,  1060}
    ,{  1990,  -310,  1990,   270,  1990}
    ,{  1840,    10,  1060,  1840,  1060}
    }
   ,{{  1990,   160,  1210,  1990,  1210}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{  1550,   -40,   770,  1550,   770}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{   270,  -780,   270,  -200,   270}
    }
   }
  ,{{{  1860,  1840,  1860,  1840,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1860,  1840,  1860,  1840,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1210,  1190,  1210,  1190,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{   930,   910,   930,   910,   120}
    ,{    50,  -210,    50,  -210, -1000}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1860,  1840,  1860,  1840,   120}
    ,{  -250,  -510,  -250,  -510, -1300}
    ,{   930,   910,   930,   910,   120}
    ,{  1860,  1840,  1860,  1840,  -200}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   640,   620,   640,   620,  -170}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   140,   120,   140,   120,  -670}
    }
   }
  }
 ,{{{{  2120,  2120,  1990,  2120,  1990}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  2120,  2120,  1990,  1990,  1990}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  1470,  1340,  2120,  1340}
    ,{  2120,  1470,  1340,  2120,  1340}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{   400,   290,   400,    10,   170}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1860,  1210,  1080,  1860,  1080}
    }
   ,{{  2120,  2120,  1990,  1840,  1990}
    ,{   540,   190,   540,   -80,    70}
    ,{  1840,  1190,  1060,  1840,  1060}
    ,{  2120,  2120,  1990,  1840,  1990}
    ,{  1840,  1190,  1060,  1840,  1060}
    }
   ,{{  1990,  1340,  1210,  1990,  1210}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{  1750,  1100,   970,  1750,   970}
    ,{  1990,  1340,  1210,  1990,  1210}
    ,{   640,   640,   270,   120,   270}
    }
   }
  ,{{{  2120,  2120,  1990,   540,  1990}
    ,{  1470,  1470,  1340,   190,  1340}
    ,{  1340,  1340,  1210,   540,  1210}
    ,{  2120,  2120,  1990,    60,  1990}
    ,{  1210,  1210,  1080,   180,  1080}
    }
   ,{{  1470,  1470,  1340,   190,  1340}
    ,{  1470,  1470,  1340,   190,  1340}
    ,{  1190,  1190,  1060,   -80,  1060}
    ,{   290,   290,   160,  -980,   160}
    ,{  1190,  1190,  1060,   -80,  1060}
    }
   ,{{  1340,  1340,  1210,   540,  1210}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{  1340,  1340,  1210,   540,  1210}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{  1210,  1210,  1080,   180,  1080}
    }
   ,{{  2120,  2120,  1990,   -80,  1990}
    ,{   190,   190,    60, -1080,    60}
    ,{  1190,  1190,  1060,   -80,  1060}
    ,{  2120,  2120,  1990,  -400,  1990}
    ,{  1190,  1190,  1060,   -80,  1060}
    }
   ,{{  1340,  1340,  1210,    70,  1210}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{  1100,  1100,   970,    70,   970}
    ,{  1340,  1340,  1210,    60,  1210}
    ,{   640,   640,   270,  -810,   270}
    }
   }
  ,{{{  1990,  1840,  1990,  1840,  1750}
    ,{  1340,  1190,  1340,  1190,  1100}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1990,  1840,  1990,  1840,  1750}
    ,{  1080,   930,  1080,   930,   840}
    }
   ,{{  1340,  1190,  1340,  1190,  1100}
    ,{  1340,  1190,  1340,  1190,  1100}
    ,{  1060,   910,  1060,   910,   820}
    ,{   400,    10,   400,    10,   170}
    ,{  1060,   910,  1060,   910,   820}
    }
   ,{{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{  1080,   930,  1080,   930,   840}
    }
   ,{{  1990,  1840,  1990,  1840,  1750}
    ,{   540,   -80,   540,   -80,    70}
    ,{  1060,   910,  1060,   910,   820}
    ,{  1990,  1840,  1990,  1840,  1750}
    ,{  1060,   910,  1060,   910,   820}
    }
   ,{{  1210,  1060,  1210,  1060,   970}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{   970,   820,   970,   820,   740}
    ,{  1210,  1060,  1210,  1060,   970}
    ,{   270,   120,   270,   120,    30}
    }
   }
  ,{{{  2120,   400,  1990,  2120,  1990}
    ,{  2120,   290,  1340,  2120,  1340}
    ,{  1990,   400,  1210,  1990,  1210}
    ,{  1990,   160,  1990,  1990,  1990}
    ,{  1860,   270,  1080,  1860,  1080}
    }
   ,{{  2120,   290,  1340,  2120,  1340}
    ,{  2120,   290,  1340,  2120,  1340}
    ,{  1840,    10,  1060,  1840,  1060}
    ,{   160,  -890,   160,  -310,   160}
    ,{  1840,    10,  1060,  1840,  1060}
    }
   ,{{  1990,   400,  1210,  1990,  1210}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{  1990,   400,  1210,  1990,  1210}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{  1860,   270,  1080,  1860,  1080}
    }
   ,{{  1990,    10,  1990,  1840,  1990}
    ,{    60,  -980,    60,  -400,    60}
    ,{  1840,    10,  1060,  1840,  1060}
    ,{  1990,  -310,  1990,   270,  1990}
    ,{  1840,    10,  1060,  1840,  1060}
    }
   ,{{  1990,   170,  1210,  1990,  1210}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{  1750,   170,   970,  1750,   970}
    ,{  1990,   160,  1210,  1990,  1210}
    ,{   270,  -780,   270,  -200,   270}
    }
   }
  ,{{{  1860,  1840,  1860,  1840,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1860,  1840,  1860,  1840,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1210,  1190,  1210,  1190,   640}
    ,{  1210,  1190,  1210,  1190,   640}
    ,{   930,   910,   930,   910,   120}
    ,{   270,    10,   270,    10,  -780}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   950,   930,   950,   930,   140}
    }
   ,{{  1860,  1840,  1860,  1840,   120}
    ,{   180,   -80,   180,   -80,  -810}
    ,{   930,   910,   930,   910,   120}
    ,{  1860,  1840,  1860,  1840,  -200}
    ,{   930,   910,   930,   910,   120}
    }
   ,{{  1080,  1060,  1080,  1060,   270}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   840,   820,   840,   820,    30}
    ,{  1080,  1060,  1080,  1060,   270}
    ,{   140,   120,   140,   120,  -670}
    }
   }
  }
 }};


PRIVATE vrna_param_t *
get_scaled_params(vrna_md_t *md)
{
  unsigned int  i, j, k, l;
  int id = 0;
  double        tempf;
  vrna_param_t  *params;
  double salt = md->salt;
  double saltStandard = VRNA_MODEL_DEFAULT_SALT;

  params = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));


  params->model_details = *md;  /* copy over the model details */
  params->temperature   = md->temperature;
  // tempf                 = ((md->temperature + K0) / Tmeasure);
  tempf = 1.;
  params->ninio[2]              = RESCALE_dG(ninio37, niniodH, tempf);
  params->lxc                   = lxc37 * tempf;
  params->TripleC               = RESCALE_dG(TripleC37, TripleCdH, tempf);
  params->MultipleCA            = RESCALE_dG(MultipleCA37, MultipleCAdH, tempf);
  params->MultipleCB            = RESCALE_dG(MultipleCB37, MultipleCBdH, tempf);
  params->TerminalAU            = RESCALE_dG(TerminalAU37, TerminalAUdH, tempf);
  params->DuplexInit            = RESCALE_dG(DuplexInit37, DuplexInitdH, tempf);
  params->MLbase                = RESCALE_dG(ML_BASE37, ML_BASEdH, tempf);
  params->MLclosing             = RESCALE_dG(ML_closing37, ML_closingdH, tempf);
  params->gquadLayerMismatch    = RESCALE_dG(GQuadLayerMismatch37, GQuadLayerMismatchH, tempf);
  params->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, tempf);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, tempf);
      params->gquad[i][j] = (int)GQuadAlpha_T * (i - 1) + (int)(((double)GQuadBeta_T) * log(j - 2));
    }

  for (i = 0; i < 31; i++)
    params->hairpin[i] = RESCALE_dG(hairpin37[i], hairpindH[i], tempf);

  for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
    params->bulge[i]          = RESCALE_dG(bulge37[i], bulgedH[i], tempf);
    params->internal_loop[i]  = RESCALE_dG(internal_loop37[i], internal_loopdH[i], tempf);
  }

  for (; i <= MAXLOOP; i++) {
    params->bulge[i] = params->bulge[30] +
                       (int)(params->lxc * log((double)(i) / 30.));
    params->internal_loop[i] = params->internal_loop[30] +
                               (int)(params->lxc * log((double)(i) / 30.));
    params->SaltLoopDbl[i]   = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]      = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (i = 0; i <= MIN2(31, MAXLOOP+1); i++) {
    params->SaltLoopDbl[i]    = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]       = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (; i <= MAXLOOP; i++) {
    params->SaltLoopDbl[i]   = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]      = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (i = 0; (i * 7) < strlen(Tetraloops); i++)
    params->Tetraloop_E[i] = RESCALE_dG(Tetraloop37[i], TetraloopdH[i], tempf);

  for (i = 0; (i * 5) < strlen(Triloops); i++)
    params->Triloop_E[i] = RESCALE_dG(Triloop37[i], TriloopdH[i], tempf);

  for (i = 0; (i * 9) < strlen(Hexaloops); i++)
    params->Hexaloop_E[i] = RESCALE_dG(Hexaloop37[i], HexaloopdH[i], tempf);

  for (i = 0; i <= NBPAIRS; i++)
    params->MLintern[i] = RESCALE_dG(ML_intern37, ML_interndH, tempf);

  /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      params->stack[i][j] = RESCALE_dG(stack37[i][j],
                                       stackdH[i][j],
                                       tempf);

  /* mismatches */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
      for (k = 0; k < 5; k++) {
        int mm;
        params->mismatchI[i][j][k] = RESCALE_dG(mismatchI37[i][j][k],
                                                mismatchIdH[i][j][k],
                                                tempf);
        params->mismatchH[i][j][k] = RESCALE_dG(mismatchH37[i][j][k],
                                                mismatchHdH[i][j][k],
                                                tempf);
        params->mismatch1nI[i][j][k] = RESCALE_dG(mismatch1nI37[i][j][k],
                                                  mismatch1nIdH[i][j][k],
                                                  tempf);
        params->mismatch23I[i][j][k] = RESCALE_dG(mismatch23I37[i][j][k],
                                                  mismatch23IdH[i][j][k],
                                                  tempf);
        if (md->dangles) {
          mm = RESCALE_dG(mismatchM37[i][j][k],
                          mismatchMdH[i][j][k],
                          tempf);
          params->mismatchM[i][j][k]  = (mm > 0) ? 0 : mm;
          mm                          = RESCALE_dG(mismatchExt37[i][j][k],
                                                   mismatchExtdH[i][j][k],
                                                   tempf);
          params->mismatchExt[i][j][k] = (mm > 0) ? 0 : mm;
        } else {
          params->mismatchM[i][j][k] = params->mismatchExt[i][j][k] = 0;
        }
      }

  /* dangles */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++) {
      int dd;
      dd = RESCALE_dG(dangle5_37[i][j],
                      dangle5_dH[i][j],
                      tempf);
      params->dangle5[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
      dd                    = RESCALE_dG(dangle3_37[i][j],
                                         dangle3_dH[i][j],
                                         tempf);
      params->dangle3[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
    }

  /* interior 1x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++)
          params->int11[i][j][k][l] = RESCALE_dG(int11_37[i][j][k][l],
                                                 int11_dH[i][j][k][l],
                                                 tempf);

  /* interior 2x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++)
            params->int21[i][j][k][l][m] = RESCALE_dG(int21_37[i][j][k][l][m],
                                                      int21_dH[i][j][k][l][m],
                                                      tempf);
        }

  /* interior 2x2 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++)
              params->int22[i][j][k][l][m][n] = RESCALE_dG(int22_37[i][j][k][l][m][n],
                                                           int22_dH[i][j][k][l][m][n],
                                                           tempf);
        }

  strncpy(params->Tetraloops, Tetraloops, 281);
  strncpy(params->Triloops, Triloops, 241);
  strncpy(params->Hexaloops, Hexaloops, 361);

  /* Salt correction for stack and multiloop */
  params->SaltStack = (salt==saltStandard) ? 0 : vrna_salt_stack(salt, saltT, md->helical_rise);
  if (salt == saltStandard)
    params->SaltMLbase = params->SaltMLclosing = 0;
  else
    vrna_salt_ml(params->SaltLoopDbl, md->saltMLLower, md->saltMLUpper, &params->SaltMLbase, &params->SaltMLclosing);
  
  params->MLclosing += params->SaltMLbase;
  params->MLclosing += params->SaltMLclosing;
  params->MLbase += params->SaltMLbase;
  for (i = 0; i <= NBPAIRS; i++)
    params->MLintern[i] += params->SaltMLbase;

  params->SaltDPXInit = 0;
  if (salt != saltStandard)
  {
    if (md->saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT)
      params->SaltDPXInit = md->saltDPXInit;
    else if (md->saltDPXInit)
      params->SaltDPXInit = vrna_salt_duplex_init(md);
  }
  params->DuplexInit += params->SaltDPXInit;

  params->id = ++id;
  return params;
}


PUBLIC vrna_param_t *
vrna_params(vrna_md_t *md)
{
  if (md) {
    return get_scaled_params(md);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_params(&md);
  }
}

#define NBASES 8
/*@notnull@*/

#ifndef INLINE
# ifdef __GNUC__
#  define INLINE inline
# else
#  define INLINE
# endif
#endif

// static const char Law_and_Order[]         = "_ACGUTXKI";
// static int        BP_pair[NBASES][NBASES] =
//   /* _  A  C  G  U  X  K  I */
// { { 0, 0, 0, 0, 0, 0, 0, 0 },
//   { 0, 0, 0, 0, 5, 0, 0, 5 },
//   { 0, 0, 0, 1, 0, 0, 0, 0 },
//   { 0, 0, 2, 0, 3, 0, 0, 0 },
//   { 0, 6, 0, 4, 0, 0, 0, 6 },
//   { 0, 0, 0, 0, 0, 0, 2, 0 },
//   { 0, 0, 0, 0, 0, 1, 0, 0 },
//   { 0, 6, 0, 0, 5, 0, 0, 0 } };

#define MAXALPHA 20       /* maximal length of alphabet */

static short  alias[MAXALPHA + 1];
static int    pair[MAXALPHA + 1][MAXALPHA + 1];
/* rtype[pair[i][j]]:=pair[j][i] */
// static int    rtype[8] = {
//   0, 2, 1, 4, 3, 6, 5, 7
// };

#ifdef _OPENMP
#pragma omp threadprivate(Law_and_Order, BP_pair, alias, pair, rtype)
#endif

/* for backward compatibility */
#define ENCODE(c) encode_char(c)

/*
将碱基字符（例如 A, C, G, U/T 等）转换为对应的数字编码
*/
static INLINE int
encode_char(char c)
{
  /* return numerical representation of base used e.g. in pair[][] */
  int code;

  c = toupper(c);
  /*
  如果 energy_set 大于 0，使用简单的字符编码方案。
  将字符变量c 转换为相对于字母 'A' 的偏移量，并加 1。
  例如，'A' -> 1, 'B' -> 2, ..., 'Z' -> 26。
  */
  if (energy_set > 0) {
    code = (int)(c - 'A') + 1;
  } else {
    /*
    如果 energy_set 小于或等于 0，使用自定义的编码方案。
    Law_and_Order[10] = "_ACGUTXKI"
    如果字符串变量c 未找到，编码为 0；否则计算字符在 Law_and_Order 中的位置。
    */
    const char *pos;
    // strchr(Law_and_Order, c) 查找字符 c 在 Law_and_Order 中的位置。
    // 如果找到，返回指向该位置的指针；否则返回 NULL。
    pos = strchr(Law_and_Order, c);
    if (pos == NULL)
      code = 0;
    else
      // pos - Law_and_Order:计算指针之间的距离，即 pos 相对于 Law_and_Order 起始位置的偏移量（索引）
      code = (int)(pos - Law_and_Order);
    // 如果编码大于 5，设为 0。
    if (code > 5)
      code = 0;
    // 如果编码大于 4，减少 1，以使 'T' 和 'U' 等效（即将 'T' 和 'U' 都编码为 4）。
    if (code > 4)
      code--;           /* make T and U equivalent */
  }

  return code;
}


/*@+boolint +charint@*/
/*@null@*/
extern char *nonstandards;


/*
  配对矩阵
*/
static INLINE void
make_pair_matrix(void)
{

  /* i, j：循环变量。
  alias[]：碱基的别名，用于处理非标准碱基。
  pair[][]：碱基对配对矩阵。
  rtype[]：反向配对矩阵。
  energy_set：指定的能量集。
  nonstandards：允许的非标准碱基对。
  NBASES：碱基的总数:8
  MAXALPHA：最大字母数量:20（根据能量集可能有所不同）。
  BP_pair[][]：默认的碱基对配对矩阵。
  */

  int i, j;
  // 默认能量集 (energy_set == 0)
    // 将前5个碱基设置为其本身。
    // 将非标准碱基（例如X和K）设置为标准碱基的别名。
  for (i = 0; i < 5; i++)
    alias[i] = (short)i;
  alias[5]  = 3;  /* X <-> G */
  alias[6]  = 2;  /* K <-> C */
  alias[7]  = 0;  /* I <-> default base '@' */
  for (i = 0; i < NBASES; i++)
    for (j = 0; j < NBASES; j++)
      pair[i][j] = BP_pair[i][j];
  // 如果noGU为真，禁止G-U配对。
  if (noGU)
    pair[3][4] = pair[4][3] = 0;
  // 如果nonstandards不为空，允许非标准碱基对
  if (nonstandards != NULL) {
    /* allow nonstandard bp's */
    for (i = 0; i < (int)strlen(nonstandards); i += 2)
      pair[encode_char(nonstandards[i])]
      [encode_char(nonstandards[i + 1])] = 7;
  }
  // 设置反向配对矩阵
  for (i = 0; i < NBASES; i++)
    for (j = 0; j < NBASES; j++)
      rtype[pair[i][j]] = pair[j][i];
}


PUBLIC int
vrna_E_ext_stem(unsigned int  type,
                int           n5d,
                int           n3d,
                vrna_param_t  *p)
{
  int energy = 0;
  // 两个相邻核苷酸都有效 (n5d >= 0 && n3d >= 0)
  // 累加终端错配的能量贡献
  if (n5d >= 0 && n3d >= 0)
    energy += p->mismatchExt[type][n5d][n3d];
  // 只有 5' 端相邻核苷酸有效 (n5d >= 0)
  // 累加 5' 悬挂末端的能量贡献
  else if (n5d >= 0)
    energy += p->dangle5[type][n5d];
  else if (n3d >= 0)
    energy += p->dangle3[type][n3d];
  // 如果碱基对类型大于 2（例如 A-U 碱基对），累加终端 A-U 的能量贡献
  if (type > 2)
    energy += p->TerminalAU;
  // 返回总能量贡献
  // printf("%d,",type);
  // printf("%d,",n5d);
  // printf("%d,",n3d);
  return energy;
}


/* compute energy of degree 2 loop (stack bulge or interior) */
PRIVATE INLINE int
E_IntLoop(int           n1,
          int           n2,
          int           type,
          int           type_2,
          int           si1,
          int           sj1,
          int           sp1,
          int           sq1,
          vrna_param_t  *P)
{
  
  int nl, ns, u, energy, salt_stack_correction, salt_loop_correction, backbones;

  salt_stack_correction = P->SaltStack;
  salt_loop_correction = 0;

  if (n1 > n2) {
    nl  = n1;
    ns  = n2;
  } else {
    nl  = n2;
    ns  = n1;
  }

  if (nl == 0) {
    return P->stack[type][type_2] + salt_stack_correction;  /* stack */
  }
  
  backbones = nl+ns+2;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    /* salt correction for loop */
    if (backbones <= MAXLOOP+1)
      salt_loop_correction = P->SaltLoop[backbones];
    else
      salt_loop_correction = vrna_salt_loop_int(backbones, P->model_details.salt, P->temperature+K0, P->model_details.backbone_length);
  }

  if (ns == 0) {
    /* bulge */
    energy = (nl <= MAXLOOP) ? P->bulge[nl] :
             (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
    if (nl == 1) {
      energy += P->stack[type][type_2];
    } else {
      if (type > 2)
        energy += P->TerminalAU;

      if (type_2 > 2)
        energy += P->TerminalAU;
    }
  } else {
    /* interior loop */
    if (ns == 1) {
      if (nl == 1)                    /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1] + salt_loop_correction;

      if (nl == 2) {
        /* 2x1 loop */
        if (n1 == 1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];

        return energy + salt_loop_correction;
      } else {
        /* 1xn loop */
        energy =
          (nl + 1 <=
           MAXLOOP) ? (P->internal_loop[nl + 1]) : (P->internal_loop[30] +
                                                    (int)(P->lxc * log((nl + 1) / 30.)));
        energy  += MIN2(MAX_NINIO, (nl - ns) * P->ninio[2]);
        energy  += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy + salt_loop_correction;
      }
    } else if (ns == 2) {
      if (nl == 2) {
        /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1] + salt_loop_correction;
      } else if (nl == 3) {
        /* 2x3 loop */
        energy  = P->internal_loop[5] + P->ninio[2];
        energy  += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy + salt_loop_correction;
      }
    }

    {
      /* generic interior loop (no else here!)*/
      u       = nl + ns;
      energy  =
        (u <=
         MAXLOOP) ? (P->internal_loop[u]) : (P->internal_loop[30] + (int)(P->lxc * log((u) / 30.)));

      energy += MIN2(MAX_NINIO, (nl - ns) * P->ninio[2]);

      energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];
    }
  }

  return energy + salt_loop_correction;
}


/*
  用 encode_char 函数将每个字符编码为数字，并将结果存储在数组 S 中,how为0时S[0]为序列长度，否则为序列末尾碱基对应的数值
*/
static INLINE short *
encode_sequence(const char  *sequence,
                short       how)
{ 
  /*
  i 和 l：用于循环和序列长度。
  S：分配内存以存储编码后的序列，长度为 l + 2。这里额外添加两个元素的目的是为了处理循环结构中的边界情况。
  */
  unsigned int  i, l = (unsigned int)strlen(sequence);
  short         *S = (short *)vrna_alloc(sizeof(short) * (l + 2));

  switch (how) {
    /* standard encoding as always used for S */
    // 遍历输入序列 sequence，调用 encode_char 函数将每个字符编码为数字，并将结果存储在数组 S 中。
    // S[l + 1] = S[1]：将序列首尾相连，用于处理循环结构。
    // S[0] = (short)l：存储序列的长度在 S[0] 中。
    case 0:
      for (i = 1; i <= l; i++)    /* make numerical encoding of sequence */
        S[i] = (short)encode_char(sequence[i - 1]);
      S[l + 1]  = S[1];
      S[0]      = (short)l;
      break;
    /* encoding for mismatches of nostandard bases (normally used for S1) */
    // 同样遍历输入序列 sequence，但是通过 alias 数组将每个字符的编码映射为新的编码，并将结果存储在数组 S 中。
    // S[l + 1] = S[1]：同样将序列首尾相连。
    // S[0] = S[l]：存储序列末尾元素的值作为 S[0]，用于处理特定情况。
    case 1:
      for (i = 1; i <= l; i++)
        S[i] = alias[(short)encode_char(sequence[i - 1])];
      S[l + 1]  = S[1];
      S[0]      = S[l];
      break;
  }

  return S;
}




PUBLIC duplexT
duplexfold(const char *s1,
           const char *s2)
{
  return duplexfold_cu(s1, s2, 1);
}

// generate result if no -e option  clean_up: 一个标志，指示是否在函数结束时清理分配的内存。
PRIVATE duplexT
duplexfold_cu(const char  *s1,
              const char  *s2,
              int         clean_up)
{
  //   i 和 j: 循环变量。
  // Emin: 最小能量值，初始化为无穷大（INF）。
  // i_min 和 j_min: 记录能量最小时的位置。
  // struc: 存储回溯得到的结构字符串。
  // mfe: 用于存储最小自由能的双链结构。
  // md: 结构体，存储模型参数。
  int       i, j, Emin = INF, i_min = 0, j_min = 0;
  char      *struc;
  duplexT   mfe;
  vrna_md_t md;
  // 获取序列长度
  n1  = (int)strlen(s1);
  n2  = (int)strlen(s2);
  // 初始化模型参数 md。
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature) > 1e-6)) {
    if (P)
      free(P);
    // 基于模型参数 md 生成新的参数集 P。
    P = vrna_params(&md);
    // 创建配对矩阵 /data/ntc/Repository/ViennaRNA-2.6.4/src/pair_mat.h
    make_pair_matrix();
  }
  // 内存分配 c: 分配一个二维数组，用于存储能量值 
  c = (int **)vrna_alloc(sizeof(int *) * (n1 + 1));
  for (i = 1; i <= n1; i++)
    c[i] = (int *)vrna_alloc(sizeof(int) * (n2 + 1));
  // 将序列转换为数字编码串 此处分配空间大于len(seq) + 1 后面的以0填充---------------------------
  // S* = len(seq) + encode(seq)
  S1  = encode_sequence(s1, 0);
  S2  = encode_sequence(s2, 0);
  // SS* = encode(seq[-1]) + encode(seq)
  SS1 = encode_sequence(s1, 1);
  SS2 = encode_sequence(s2, 1);
  // for(i = 0; i <= n1 + 11; i++){
  //   printf("%d ",SS1[i]);
  // }
  // printf("%d\n", SS1[0]);
  // printf("%d\n", SS1[1]);
  // printf("%d\n", SS1[-1]);
  // 计算能量
  int no_type_num = 0;
  for (i = 1; i <= n1; i++) {
    for (j = n2; j > 0; j--) {
      // P->DuplexInit: 初始化能量。
      // vrna_E_ext_stem: 计算外部茎的能量。
      // E_IntLoop: 计算内部环的能量。
      // MIN2: 返回两个值中的最小值。
      /* c: energy array, given that i-j pair */
      int type, type2, E, k, l;
      // type: 序列 s1[i] 和 s2[j] 的配对类型 例如，A-U，G-C。
      // pair no erro !!!!!!!!!!!!!!!!!!!!
      type    = pair[S1[i]][S2[j]];
      // printf("%d,",S1[i]);
      // printf("%d,",S2[j]);
      // printf("%d,",type);
      // print_pair();
      // 初始化能量。如果 type 有效，则设置初始值，否则设置为无穷大。
      c[i][j] = type ? P->DuplexInit : INF;
      // printf("%d,", E);
      // printf("%d,", type);
      // printf("%d,", c[i][j]);
      // 若配对的type 自由能为0 则跳过能量累加计算
      if (!type){
        continue;

      }
      // 计算外部茎的能量 使用 type 和相邻的碱基（如果存在）计算能量。
      int n5dd = (j < n2) ? SS2[j + 1] : -1;
      c[i][j] += vrna_E_ext_stem(type, (i > 1) ? SS1[i - 1] : -1, (j < n2) ? SS2[j + 1] : -1, P);
      // no_type_num += c[i][j];
      // continue;
      // printf("%d,", c[i][j]);
      // no_type_num += c[i][j];
      // printf("%d,", c[i][j]);
      for (k = i - 1; k > 0 && k > i - MAXLOOP - 2; k--) {
        for (l = j + 1; l <= n2; l++) {
          if (i - k + l - j - 2 > MAXLOOP)
            break;
          type2 = pair[S1[k]][S2[l]];
          if (!type2)
            continue;

          E = E_IntLoop(i - k - 1, l - j - 1, type2, rtype[type],
                        SS1[k + 1], SS2[l - 1], SS1[i - 1], SS2[j + 1], P);
          c[i][j] = MIN2(c[i][j], c[k][l] + E);
          no_type_num += l;
        }
      }
      // different c !!!!!!!!!!
      E = c[i][j];
      E += vrna_E_ext_stem(rtype[type], (j > 1) ? SS2[j - 1] : -1, (i < n1) ? SS1[i + 1] : -1, P);
      if (E < Emin) {
        Emin  = E;
        i_min = i;
        j_min = j;
      }
    }
  }
  // printf("%d",i_min);
  // printf("%d",j_min);
  // printf("%d",no_type_num);

  // 回溯生成结构
  struc = backtrack(i_min, j_min);
  if (i_min < n1)
    i_min++;

  if (j_min > 1)
    j_min--;
  // 填充 mfe 结构体
  mfe.i         = i_min;
  mfe.j         = j_min;
  mfe.energy    = (float)Emin / 100.;
  mfe.structure = struc;
  // 清理内存
  if (clean_up) {
    for (i = 1; i <= n1; i++)
      free(c[i]);
    free(c);
    free(S1);
    free(S2);
    free(SS1);
    free(SS2);
  }
  // 返回结果
  return mfe;
}

#define MAX2(A, B)      ((A) > (B) ? (A) : (B))
// 用于计算给定的两个RNA序列 s1 和 s2 的次优（suboptimal）二级结构。这些次优结构的自由能在最低自由能结构（MFE）的基础上增加一个给定的能量范围 delta 以内。函数返回一个包含这些次优结构的数组 w: 一个窗口参数，用于限制次优结构的数量
PUBLIC duplexT *

duplex_subopt(const char  *s1,
              const char  *s2,
              int         delta,
              int         w)
{
  /*   i, j: 循环变量。
  n1, n2: 分别表示序列 s1 和 s2 的长度。
  thresh: 能量阈值。
  E: 当前能量值。
  n_subopt: 记录找到的次优结构的数量。
  n_max: 次优结构数组的初始大小。
  struc: 保存回溯得到的次优结构。
  mfe: 最低自由能结构。
  subopt: 存储次优结构的数组。*/

  int     i, j, n1, n2, thresh, E, n_subopt = 0, n_max;
  char    *struc;
  duplexT mfe;
  duplexT *subopt;


  /*
  初始化 n_max 为 16。
  使用 vrna_alloc 分配初始大小为 16 的 duplexT 数组 subopt。
  计算并获取最低自由能结构 mfe。
  释放 mfe.structure 所占用的内存。
  */

  n_max   = 16;
  subopt  = (duplexT *)vrna_alloc(n_max * sizeof(duplexT));
  mfe     = duplexfold_cu(s1, s2, 0);
  free(mfe.structure);

  // 根据 mfe.energy 和 delta 计算能量阈值 threshold。
  thresh  = (int)mfe.energy * 100 + 0.1 + delta * 100;
  n1      = strlen(s1);
  n2      = strlen(s2);
  // print_pair();
  // 二层遍历查找次优结构
  for (i = n1; i > 0; i--) {
    for (j = 1; j <= n2; j++) {
      int type, ii, jj, Ed;
      type = pair[S2[j]][S1[i]];
      if (!type)
        continue;
      // 计算能量
      E   = Ed = c[i][j];
      // from /data/ntc/Repository/ViennaRNA-2.6.4/src/loops/external.c
      Ed  += vrna_E_ext_stem(type, (j > 1) ? SS2[j - 1] : -1, (i < n1) ? SS1[i + 1] : -1, P);
      printf("%d,",Ed);
      // 超过阈值则跳过
      if (Ed > thresh)
        continue;

      /* too keep output small, remove hits that are dominated by a
       * better one close (w) by. For simplicity we do test without
       * adding dangles, which is slightly inaccurate.
       */
      for (ii = MAX2(i - w, 1); (ii <= MIN2(i + w, n1)) && type; ii++) {
        for (jj = MAX2(j - w, 1); jj <= MIN2(j + w, n2); jj++)
          if (c[ii][jj] < E) {
            type = 0;
            break;
          }
      }
      // printf("%d,", thresh);
      if (!type)
        continue;

      struc = backtrack(i, j);
      
      // 如果 subopt 数组已满，使用 vrna_realloc 扩展数组大小
      if (n_subopt + 1 >= n_max) {
        n_max   *= 2;
        subopt  = (duplexT *)vrna_realloc(subopt, n_max * sizeof(duplexT));
      }

      subopt[n_subopt].i            = MIN2(i + 1, n1);
      subopt[n_subopt].j            = MAX2(j - 1, 1);
      subopt[n_subopt].energy       = Ed * 0.01;
      subopt[n_subopt++].structure  = struc;
    }
  }
  /* free all static globals */
  for (i = 1; i <= n1; i++)
    free(c[i]);
  free(c);
  free(S1);
  free(S2);
  free(SS1);
  free(SS2);

  subopt[n_subopt].i          = 0;
  subopt[n_subopt].j          = 0;
  subopt[n_subopt].structure  = NULL;
  return subopt;
}


PRIVATE char *
backtrack(int i,
          int j)
{
  /* backtrack structure going backwards from i, and forwards from j
   * return structure in bracket notation with & as separator */
  int   k, l, type, type2, E, traced, i0, j0;
  char  *st1, *st2, *struc;

  st1 = (char *)vrna_alloc(sizeof(char) * (n1 + 1));
  st2 = (char *)vrna_alloc(sizeof(char) * (n2 + 1));

  i0  = MIN2(i + 1, n1);
  j0  = MAX2(j - 1, 1);

  while (i > 0 && j <= n2) {
    E           = c[i][j];
    traced      = 0;
    st1[i - 1]  = '(';
    st2[j - 1]  = ')';
    type        = pair[S1[i]][S2[j]];
    if (!type)
      printf("backtrack failed in fold duplex");

    for (k = i - 1; k > 0 && k > i - MAXLOOP - 2; k--) {
      for (l = j + 1; l <= n2; l++) {
        int LE;
        if (i - k + l - j - 2 > MAXLOOP)
          break;

        type2 = pair[S1[k]][S2[l]];
        if (!type2)
          continue;

        LE = E_IntLoop(i - k - 1, l - j - 1, type2, rtype[type],
                       SS1[k + 1], SS2[l - 1], SS1[i - 1], SS2[j + 1], P);
        if (E == c[k][l] + LE) {
          traced  = 1;
          i       = k;
          j       = l;
          break;
        }
      }
      if (traced)
        break;
    }
    if (!traced) {
      E -= vrna_E_ext_stem(type, (i > 1) ? SS1[i - 1] : -1, (j < n2) ? SS2[j + 1] : -1, P);
      if (E != P->DuplexInit)
        printf("backtrack failed in fold duplex");
      else
        break;
    }
  }
  if (i > 1)
    i--;

  if (j < n2)
    j++;

  struc = (char *)vrna_alloc(i0 - i + 1 + j - j0 + 1 + 2);
  for (k = MAX2(i, 1); k <= i0; k++)
    if (!st1[k - 1])
      st1[k - 1] = '.';

  for (k = j0; k <= j; k++)
    if (!st2[k - 1])
      st2[k - 1] = '.';

  strcpy(struc, st1 + MAX2(i - 1, 0));
  strcat(struc, "&");
  strcat(struc, st2 + j0 - 1);

  /* printf("%s %3d,%-3d : %3d,%-3d\n", struc, i,i0,j0,j);  */
  free(st1);
  free(st2);

  return struc;
}


PRIVATE void
print_struc(duplexT const *dup)
{
  int   l1;

  l1 = strchr(dup->structure, '&') - dup->structure;
  // char  *msg = fprintf(" %3d,%-3d : %3d,%-3d (%5.2f)",
  //                                 dup->i + 1 - l1,
  //                                 dup->i,
  //                                 dup->j,
  //                                 dup->j + (int)strlen(dup->structure) - l1 - 2,
  //                                 dup->energy);
  printf( " %3d,",dup->i + 1 - l1);
  printf( "%-3d :",dup->i);
  printf( " %3d,",dup->j);
  printf( "%-3d",dup->j + (int)strlen(dup->structure) - l1 - 2);
  printf( " (%5.2f)\n", dup->energy);
  printf("%s\n",dup->structure);

  // free(msg);
}

void process(char *s1, char *s2){
  int                               i, sym, istty, delta, noconv;
    delta     = 10;
    duplexT mfe, *subopt;
    // char *s1 = "CTAGCATGCTACG";
    // char *s2 = "CGTAGCATGCTAG";
    if (delta >= 0) {
      duplexT *sub;
      subopt = duplex_subopt(s1, s2, delta, 5);
      for (sub = subopt; sub->i > 0; sub++) {
        print_struc(sub);
        free(sub->structure);
      }
      free(subopt);
    } else {
      mfe = duplexfold(s1, s2);
      print_struc(&mfe);
      free(mfe.structure);
    }
}


int main(){
    char *s1 = "CTTCCTCGGGTTCAAAGCTGGATT";
    char *s2 = "GTCCAGTTTTCCCAGGAAT";
    process(s1,s2);
    return 0;
}