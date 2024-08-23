#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define STACK_BULGE1 1 /* stacking energies for bulges of size 1 */
#define NEW_NINIO 1    /* new asymetry penalty */
#define MAXSECTORS 500 /* dimension for a backtrack array */
#define LOCALITY 0.    /* locality parameter for base-pairs */
#define UNIT 100
#define MINPSCORE -2 * UNIT
#define NONE -10000 /* score for forbidden pairs */
#define PUBLIC
#define PRIVATE static
#define MAXALPHA 20 /* maximal length of alphabet */

static short alias[MAXALPHA + 1];
static int pair[MAXALPHA + 1][MAXALPHA + 1];
/* rtype[pair[i][j]]:=pair[j][i] */
static int rtype[8] = {
    0, 2, 1, 4, 3, 6, 5, 7};
#define NBASES 8
static const char Law_and_Order[] = "_ACGUTXKI";
static int BP_pair[NBASES][NBASES] =
    /* _  A  C  G  U  X  K  I */
    {{0, 0, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 5, 0, 0, 5},
     {0, 0, 0, 1, 0, 0, 0, 0},
     {0, 0, 2, 0, 3, 0, 0, 0},
     {0, 6, 0, 4, 0, 0, 0, 6},
     {0, 0, 0, 0, 0, 0, 2, 0},
     {0, 0, 0, 0, 0, 1, 0, 0},
     {0, 6, 0, 0, 5, 0, 0, 0}};

/** The gas constant */
#define GASCONST 1.98717 /* in [cal/K] */
/** 0 deg Celsius in Kelvin */
#define K0 273.15
/** Infinity as used in minimization routines */
#define INF 10000000 /* (INT_MAX/10) */

#define EMAX (INF / 10)
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

#define VRNA_GQUAD_MAX_STACK_SIZE 7
#define VRNA_GQUAD_MIN_STACK_SIZE 2
#define VRNA_GQUAD_MAX_LINKER_LENGTH 15
#define VRNA_GQUAD_MIN_LINKER_LENGTH 1
#define VRNA_GQUAD_MIN_BOX_SIZE ((4 * VRNA_GQUAD_MIN_STACK_SIZE) + \
                                 (3 * VRNA_GQUAD_MIN_LINKER_LENGTH))
#define VRNA_GQUAD_MAX_BOX_SIZE ((4 * VRNA_GQUAD_MAX_STACK_SIZE) + \
                                 (3 * VRNA_GQUAD_MAX_LINKER_LENGTH))

#define VRNA_MODEL_DEFAULT_TEMPERATURE 37.0

#define VRNA_MODEL_DEFAULT_TEMPERATURE 37.0
#define VRNA_MODEL_DEFAULT_PF_SCALE -1

/**
 *  @brief  Default scaling factor for absolute thermodynamic temperature in Boltzmann factors
 *
 *  @see    #vrna_exp_param_t.alpha, #vrna_md_t.betaScale, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BETA_SCALE 1.

/**
 *  @brief  Default dangling end model
 *
 *  @see  #vrna_md_t.dangles, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_DANGLES 2
// #define VRNA_MODEL_DEFAULT_DANGLES 2

/**
 *  @brief  Default model behavior for lookup of special tri-, tetra-, and hexa-loops
 *
 *  @see    #vrna_md_t.special_hp, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_SPECIAL_HP 1

/**
 *  @brief  Default model behavior for so-called 'lonely pairs'
 *
 *  @see    #vrna_md_t.noLP, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_LP 0

/**
 *  @brief  Default model behavior for G-U base pairs
 *
 *  @see    #vrna_md_t.noGU, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU 0

/**
 *  @brief  Default model behavior for G-U base pairs closing a loop
 *
 *  @see    #vrna_md_t.noGUclosure, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_NO_GU_CLOSURE 1

/**
 *  @brief  Default model behavior to treat a molecule as a circular RNA (DNA)
 *  @see    #vrna_md_t.circ, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_CIRC 0

/**
 *  @brief  Default model behavior regarding the treatment of G-Quadruplexes
 *
 *  @see    #vrna_md_t.gquad, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_GQUAD 0

/**
 *  @brief  Default behavior of the model regarding unique multi-branch loop decomposition
 *
 *  @see    #vrna_md_t.uniq_ML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_UNIQ_ML 0

/**
 *  @brief  Default model behavior on which energy set to use
 *
 *  @see    #vrna_md_t.energy_set, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ENERGY_SET 0

/**
 *  @brief  Default model behavior with regards to backtracking of structures
 *  @see    #vrna_md_t.backtrack, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_BACKTRACK 1

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
#define VRNA_MODEL_DEFAULT_COMPUTE_BPP 1

/**
 *  @brief  Default model behavior for the allowed maximum base pair span
 *
 *  @see    #vrna_md_t.max_bp_span, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_MAX_BP_SPAN -1

/**
 *  @brief  Default model behavior for the sliding window approach
 *
 *  @see    #vrna_md_t.window_size, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_WINDOW_SIZE -1

/**
 *  @brief  Default model behavior on how to evaluate the energy contribution of multi-branch loops
 *
 *  @see    #vrna_md_t.logML, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_LOG_ML 0

/**
 *  @brief  Default model behavior for consensus structure energy evaluation
 *  @see    #vrna_md_t.oldAliEn, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_OLD_EN 0

/**
 *  @brief  Default model behavior for consensus structure co-variance contribution assessment
 *
 *  @see    #vrna_md_t.ribo, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_RIBO 0

/**
 *  @brief  Default model behavior for weighting the co-variance score in consensus structure prediction
 *
 *  @see    #vrna_md_t.cv_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_CV_FACT 1.

/** @brief  Default model behavior for weighting the nucleotide conservation? in consensus structure prediction
 *
 *  @see    #vrna_md_t.nc_fact, vrna_md_defaults_reset(), vrna_md_set_default()
 */
#define VRNA_MODEL_DEFAULT_ALI_NC_FACT 1.

#define VRNA_MODEL_DEFAULT_PF_SMOOTH 1

/**
 *  @brief  Default model salt concentration (M)
 */
#define VRNA_MODEL_DEFAULT_SALT 1.021

/**
 *  @brief  Default model lower bound of multiloop size for salt correction fiting
 */
#define VRNA_MODEL_DEFAULT_SALT_MLLOWER 6

/**
 *  @brief  Default model upper bound of multiloop size for salt correction fiting
 */
#define VRNA_MODEL_DEFAULT_SALT_MLUPPER 24

/**
 *  @brief  Default model value to turn off user-provided salt correction for duplex initializtion
 */
#define VRNA_MODEL_DEFAULT_SALT_DPXINIT 99999

#define VRNA_MODEL_SALT_DPXINIT_FACT_RNA -45.324
#define VRNA_MODEL_SALT_DPXINIT_FACT_DNA -58.389

#define VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT VRNA_MODEL_SALT_DPXINIT_FACT_RNA

/* Geometric parameters for RNA and DNA */

#define VRNA_MODEL_HELICAL_RISE_RNA 2.8
#define VRNA_MODEL_HELICAL_RISE_DNA 3.4
/**
 *  @brief  Default helical rise
 */
#define VRNA_MODEL_DEFAULT_HELICAL_RISE VRNA_MODEL_HELICAL_RISE_RNA

#define VRNA_MODEL_BACKBONE_LENGTH_RNA 6.0
#define VRNA_MODEL_BACKBONE_LENGTH_DNA 6.76
/**
 *  @brief  Default backbone length
 */
#define VRNA_MODEL_DEFAULT_BACKBONE_LENGTH VRNA_MODEL_BACKBONE_LENGTH_RNA

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#ifndef MAXALPHA
/**
 *  @brief Maximal length of alphabet
 */
#define MAXALPHA 20

#endif

#endif

#define BP_REV_DEFAULT {0, 2, 1, 4, 3, 6, 5, 7}

#define BP_ALIAS_DEFAULT {0, 1, 2, 3, 4, 3, 2, 0}

#define BP_ENCODING_DEFAULT       \
  /*  _  A  C  G  U  X  K  I */   \
  {                               \
    {0, 0, 0, 0, 0, 0, 0, 0},     \
        {0, 0, 0, 0, 5, 0, 0, 5}, \
        {0, 0, 0, 1, 0, 0, 0, 0}, \
        {0, 0, 2, 0, 3, 0, 0, 0}, \
        {0, 6, 0, 4, 0, 0, 0, 6}, \
        {0, 0, 0, 0, 0, 0, 2, 0}, \
        {0, 0, 0, 0, 0, 1, 0, 0}, \
    {                             \
      0, 6, 0, 0, 5, 0, 0, 0      \
    }                             \
  }

#define DM_DEFAULT                                              \
  {                                                             \
    {0, 0, 0, 0, 0, 0, 0}, /* hamming distance between pairs */ \
        {0, 0, 2, 2, 1, 2, 2} /* CG */,                         \
        {0, 2, 0, 1, 2, 2, 2} /* GC */,                         \
        {0, 2, 1, 0, 2, 1, 2} /* GU */,                         \
        {0, 1, 2, 2, 0, 2, 1} /* UG */,                         \
        {0, 2, 2, 1, 2, 0, 2} /* AU */,                         \
    {                                                           \
      0, 2, 2, 2, 1, 2, 0                                       \
    } /* UA */                                                  \
  }

typedef double FLT_OR_DBL;
typedef struct vrna_fc_s vrna_fold_compound_t;

double temperature = VRNA_MODEL_DEFAULT_TEMPERATURE;
double pf_scale = VRNA_MODEL_DEFAULT_PF_SCALE;
int dangles = VRNA_MODEL_DEFAULT_DANGLES;
int tetra_loop = VRNA_MODEL_DEFAULT_SPECIAL_HP;
int noLonelyPairs = VRNA_MODEL_DEFAULT_NO_LP;
int noGU = 1;
int no_closingGU = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
int circ = VRNA_MODEL_DEFAULT_CIRC;
int gquad = VRNA_MODEL_DEFAULT_GQUAD;
int uniq_ML = VRNA_MODEL_DEFAULT_UNIQ_ML;
int energy_set = VRNA_MODEL_DEFAULT_ENERGY_SET;
int do_backtrack = VRNA_MODEL_DEFAULT_COMPUTE_BPP;
char backtrack_type = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE;
char *nonstandards = NULL;
int max_bp_span = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
int oldAliEn = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
int ribo = VRNA_MODEL_DEFAULT_ALI_RIBO;
double cv_fact = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
double nc_fact = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
int logML = VRNA_MODEL_DEFAULT_LOG_ML;

/* below are some more deprecated global symbols we need to get rid off */

int james_rule = 1;       /* interior loops of size 2 get energy 0.8Kcal and
                           * no mismatches (no longer used) */
char *RibosumFile = NULL; /* TODO: compile ribosums into program
                           * Warning: this variable will vanish */
int csv = 0;              /*generate comma seperated output*/
// vrna_bp_stack_t *base_pair        = NULL;
FLT_OR_DBL *pr = NULL;    /* base pairing prob. matrix */
int *iindx = NULL;        /* pr[i,j] -> pr[iindx[i]-j] */
int fold_constrained = 0; /* fold with constraints */

double salt = VRNA_MODEL_DEFAULT_SALT;
int saltDPXInit = VRNA_MODEL_DEFAULT_SALT_DPXINIT;
float saltDPXInitFact = VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT;
float helical_rise = VRNA_MODEL_DEFAULT_HELICAL_RISE;
float backbone_length = VRNA_MODEL_DEFAULT_BACKBONE_LENGTH;



/* math:expn and kn definition*/
#define BIG 1.44115188075855872e+17
#define MAXLOG 7.09782712893383996843e+02
#define MACHEP 1.11022302462515654042e-16
#define MAXNUM 1.7976931348623158e+308
#define EUL 0.57721566490153286060
#define PI 3.14159265358979323846

#define K0 273.15
#define Tmeasure 37 + K0

// #define lxc37  107.9
// #define VRNA_MODEL_DEFAULT_SALT  1.021

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
PUBLIC double lxc37 = 107.856;

PUBLIC int ML_BASE37 = 20;
PUBLIC int ML_BASEdH = 0;
PUBLIC int ML_closing37 = 300;
PUBLIC int ML_closingdH = 900;
PUBLIC int ML_intern37 = 20;
PUBLIC int ML_interndH = 0;

PUBLIC int ninio37 = 40;
PUBLIC int niniodH = 0;
PUBLIC int MAX_NINIO = 300;

PUBLIC int DuplexInit37 = 100;
PUBLIC int DuplexInitdH = -720;
PUBLIC int TerminalAU37 = 0;
PUBLIC int TerminalAUdH = 320;

PUBLIC int TripleC37 = 100;
PUBLIC int TripleCdH = 1860;
PUBLIC int MultipleCA37 = 30;
PUBLIC int MultipleCAdH = 340;
PUBLIC int MultipleCB37 = 160;
PUBLIC int MultipleCBdH = 760;

PUBLIC int GQuadAlpha37 = -1800;
PUBLIC int GQuadAlphadH = -11934;
PUBLIC int GQuadBeta37 = 1200;
PUBLIC int GQuadBetadH = 0;
PUBLIC int GQuadLayerMismatch37 = 300;
PUBLIC int GQuadLayerMismatchH = 0;
PUBLIC int GQuadLayerMismatchMax = 1;


PUBLIC int stack37[NBPAIRS + 1][NBPAIRS + 1] =
    {{INF, INF, INF, INF, INF, INF, INF, INF}, {INF, -220, -180, -30, -50, -130, -150, -30}, {INF, -180, -220, -60, 10, -140, -130, 10}, {INF, -30, -60, 120, 70, 10, 30, 120}, {INF, -50, 10, 70, 50, 70, 40, 70}, {INF, -130, -140, 10, 70, -90, -100, 70}, {INF, -150, -130, 30, 40, -100, -60, 40}, {INF, -30, 10, 120, 70, 70, 40, 120}};
PUBLIC int stackdH[NBPAIRS + 1][NBPAIRS + 1] =
    {{INF, INF, INF, INF, INF, INF, INF, INF}, {INF, -980, -750, -240, -430, -580, -990, -240}, {INF, -750, -790, -240, 430, -780, -850, 430}, {INF, -240, -240, 500, 800, -170, -90, 800}, {INF, -430, 430, 800, -670, 340, -90, 800}, {INF, -580, -780, -170, 370, -530, -720, 370}, {INF, -990, -850, -90, -90, -720, -670, -90}, {INF, -240, 430, 800, 800, 340, -90, 800}};

PUBLIC int hairpin37[31] = {INF,   INF,   INF,   340,   340,   350,   420,   420,   420,   430, 440, 450, 460, 470, 480, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650};
PUBLIC int hairpindH[31] = {INF,   INF,   INF,   -60,  -230, -1210, -1670, -1360, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410};
PUBLIC int bulge37[31] = {INF,   290,   230,   250,   270,   300,   320,   340,   350,   360, 370, 390, 390, 400, 410, 420, 430, 430, 440, 440, 450, 450, 460, 460, 470, 470, 470, 480, 490, 490, 490};
PUBLIC int bulgedH[31] = {INF,  1890,   -60,  -230, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410, -1410,-1410, -1410, -1410, -1410};
PUBLIC int internal_loop37[31] = {INF,   INF,   INF,   INF,   310,   350,   390,   410,   420,   430,
     450, 460, 460, 470, 480, 490, 500, 500, 510, 510 ,
     520, 530, 530, 530, 540, 540, 550, 550, 560, 560 ,
     560};
PUBLIC int internal_loopdH[31] = {INF,   INF,   INF,   INF,     0,     0,     0,     0,     0,     0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0};

PUBLIC int mismatchI37[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, {{ 0, 0, 0, 0, 0 },
    { 0, -100, -80, -90, 0 },
    { 0, -80, -50, 0, -70 },
    { 0, -100, 0, -90, -100 },
    { 0, 0, -60, -90, -90 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -100, -70, -80, 0 },
    { 0, -100, -60, 0, -70 },
    { 0, -100, 0, -100, -80 },
    { 0, 0, -60, -90, -90 }},
    {{ -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -60 },
    { -20, -20, -20, -90, -20 },
    { -50, -50, -90, -50, -50 },
    { -20, -50, -20, -50, -20 }},
    {{ -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -50 },
    { -20, -20, -20, -80, -20 },
    { -50, -50, -100, -50, -50 },
    { -20, -50, -20, -50, -20 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -70, -30, -50, 0 },
    { 0, -60, -20, 0, -30 },
    { 0, -60, 0, -40, -50 },
    { 0, 0, -30, -50, -40 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -60, -40, -50, 0 },
    { 0, -50, -20, 0, -50 },
    { 0, -60, 0, -40, -50 },
    { 0, 0, -30, -60, -30 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -50, -20, -50, 0 },
    { 0, -20, -20, 0, -20 },
    { 0, -50, 0, -40, -50 },
    { 0, 0, -20, -50, -20 }}};
PUBLIC int mismatchIdH[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, {{ -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 }},
    {{ -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 }},
    {{ 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 }},
    {{ 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 }},
    {{ 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 }},
    {{ 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 }},
    {{ 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 }}};

PUBLIC int mismatchH37[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, 
    {{ 0, 0, 0, 0, 0 },
    { 0, -100, -80, -90, 0 },
    { 0, -80, -50, 0, -70 },
    { 0, -100, 0, -90, -100 },
    { 0, 0, -60, -90, -90 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -100, -70, -80, 0 },
    { 0, -100, -60, 0, -70 },
    { 0, -100, 0, -100, -80 },
    { 0, 0, -60, -90, -90 }},
    {{ -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -60 },
    { -20, -20, -20, -90, -20 },
    { -50, -50, -90, -50, -50 },
    { -20, -50, -20, -50, -20 }},
    {{ -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -50 },
    { -20, -20, -20, -80, -20 },
    { -50, -50, -100, -50, -50 },
    { -20, -50, -20, -50, -20 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -70, -30, -50, 0 },
    { 0, -60, -20, 0, -30 },
    { 0, -60, 0, -40, -50 },
    { 0, 0, -30, -50, -40 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -60, -40, -50, 0 },
    { 0, -50, -20, 0, -50 },
    { 0, -60, 0, -40, -50 },
    { 0, 0, -30, -60, -30 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -50, -20, -50, 0 },
    { 0, -20, -20, 0, -20 },
    { 0, -50, 0, -40, -50 },
    { 0, 0, -20, -50, -20 }}};
PUBLIC int mismatchHdH[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, 
    {{ -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 }},
    {{ -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 }},
    {{ 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 }},
    {{ 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 }},
    {{ 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 }},
    {{ 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 }},
    {{ 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 }}};

/* mismatch_multi */
PUBLIC int mismatchM37[NBPAIRS + 1][5][5] =
    {{/* NP.. */
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF}},
     {/* CG.. */
      { 0, 0, 0, 0, 0 },
    { 0, -100, -100, -100, 0 },
    { 0, -70, -60, 0, -60 },
    { 0, -80, 0, -100, -90 },
    { 0, 0, -70, -80, -90 }},
     {/* GC.. */
    { 0, 0, 0, 0, 0 },
    { 0, -100, -80, -100, 0 },
    { 0, -80, -50, 0, -60 },
    { 0, -90, 0, -90, -90 },
    { 0, 0, -70, -100, -90 }},
     {/* GU.. */
{ -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -50 },
    { -20, -20, -20, -100, -20 },
    { -50, -50, -80, -50, -50 },
    { -20, -50, -20, -50, -20 }},
     {/* UG.. */
    { -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -50 },
    { -20, -20, -20, -90, -20 },
    { -50, -50, -90, -50, -50 },
    { -20, -60, -20, -50, -20 }},
     {/* AU.. */
    { 0, 0, 0, 0, 0 },
    { 0, -60, -50, -60, 0 },
    { 0, -40, -20, 0, -30 },
    { 0, -50, 0, -40, -60 },
    { 0, 0, -50, -50, -30 }},
     {/* UA.. */
    { 0, 0, 0, 0, 0 },
    { 0, -70, -60, -60, 0 },
    { 0, -30, -20, 0, -30 },
    { 0, -50, 0, -40, -50 },
    { 0, 0, -30, -50, -40 }},
     {/* NN.. */
    { 0, 0, 0, 0, 0 },
    { 0, -50, -20, -50, 0 },
    { 0, -20, -20, 0, -20 },
    { 0, -50, 0, -40, -50 },
    { 0, 0, -20, -50, -20 }}};

/* mismatch_multi_enthalpies */
PUBLIC int mismatchMdH[NBPAIRS + 1][5][5] =
    {{/* NP.. */
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF}},
     {/* CG.. */
      { -90, -90, -90, -90, -90 },
    { -620, -620, -620, -620, -620 },
    { -330, -330, -330, -330, -330 },
    { -510, -510, -510, -510, -510 },
    { -90, -90, -90, -90, -90 }},
     {/* GC.. */
      { -420, -420, -420, -420, -420 },
    { -460, -460, -460, -460, -460 },
    { -420, -420, -420, -420, -420 },
    { -470, -470, -470, -470, -470 },
    { -500, -500, -500, -500, -500 }},
     {/* GU.. */
      { 230, 230, 230, 230, 230 },
    { 230, 230, 230, 230, 230 },
    { 170, 170, 170, 170, 170 },
    { 70, 70, 70, 70, 70 },
    { -580, -580, -580, -580, -580 }},
     {/* UG.. */
      { 230, 230, 230, 230, 230 },
    { -300, -300, -300, -300, -300 },
    { -10, -10, -10, -10, -10 },
    { -190, -190, -190, -190, -190 },
    { 230, 230, 230, 230, 230 }},
     {/* AU.. */
      { 230, 230, 230, 230, 230 },
    { 230, 230, 230, 230, 230 },
    { 170, 170, 170, 170, 170 },
    { 70, 70, 70, 70, 70 },
    { -580, -580, -580, -580, -580 }},
     {/* UA.. */
      { 520, 520, 520, 520, 520 },
    { 400, 400, 400, 400, 400 },
    { 520, 520, 520, 520, 520 },
    { 470, 470, 470, 470, 470 },
    { 430, 430, 430, 430, 430 }},
     {/* NN.. */
      { 520, 520, 520, 520, 520 },
    { 400, 400, 400, 400, 400 },
    { 520, 520, 520, 520, 520 },
    { 470, 470, 470, 470, 470 },
    { 430, 430, 430, 430, 430 }}};

PUBLIC int mismatch1nI37[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }}};
PUBLIC int mismatch1nIdH[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0 }},
    {{ 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 }},
    {{ 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 }},
    {{ 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 }},
    {{ 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 }},
    {{ 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 },
    { 320, 320, 320, 320, 320 }}};

PUBLIC int mismatch23I37[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, {{ 0, 0, 0, 0, 0 },
    { 0, -100, -80, -90, 0 },
    { 0, -80, -50, 0, -70 },
    { 0, -100, 0, -90, -100 },
    { 0, 0, -60, -90, -90 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -100, -70, -80, 0 },
    { 0, -100, -60, 0, -70 },
    { 0, -100, 0, -100, -80 },
    { 0, 0, -60, -90, -90 }},
    {{ -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -60 },
    { -20, -20, -20, -90, -20 },
    { -50, -50, -90, -50, -50 },
    { -20, -50, -20, -50, -20 }},
    {{ -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -50 },
    { -20, -20, -20, -80, -20 },
    { -50, -50, -100, -50, -50 },
    { -20, -50, -20, -50, -20 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -70, -30, -50, 0 },
    { 0, -60, -20, 0, -30 },
    { 0, -60, 0, -40, -50 },
    { 0, 0, -30, -50, -40 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -60, -40, -50, 0 },
    { 0, -50, -20, 0, -50 },
    { 0, -60, 0, -40, -50 },
    { 0, 0, -30, -60, -30 }},
    {{ 0, 0, 0, 0, 0 },
    { 0, -50, -20, -50, 0 },
    { 0, -20, -20, 0, -20 },
    { 0, -50, 0, -40, -50 },
    { 0, 0, -20, -50, -20 }}};
PUBLIC int mismatch23IdH[NBPAIRS + 1][5][5] =
    {{{INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}, {INF, INF, INF, INF, INF}}, {{ -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 }},
    {{ -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 },
    { -90, -620, -330, -510, -90 }},
    {{ 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 },
    { 230, -300, -10, -190, 230 }},
    {{ 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 }},
    {{ 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 }},
    {{ 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 },
    { 230, 230, 170, 70, -580 }},
    {{ 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 },
    { 520, 400, 520, 470, 430 }}};

/* mismatch_exterior */
PUBLIC int mismatchExt37[NBPAIRS + 1][5][5] =
    {{/* NP.. */
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF}},
     {/* CG.. */
      { 0, 0, 0, 0, 0 },
    { 0, -100, -100, -100, 0 },
    { 0, -70, -60, 0, -60 },
    { 0, -80, 0, -100, -90 },
    { 0, 0, -70, -80, -90 }},
     {/* GC.. */
      { 0, 0, 0, 0, 0 },
    { 0, -100, -80, -100, 0 },
    { 0, -80, -50, 0, -60 },
    { 0, -90, 0, -90, -90 },
    { 0, 0, -70, -100, -90 }},
     {/* GU.. */
      { -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -50 },
    { -20, -20, -20, -100, -20 },
    { -50, -50, -80, -50, -50 },
    { -20, -50, -20, -50, -20 }},
     {/* UG.. */
      { -20, -20, -20, -50, -20 },
    { -20, -50, -20, -50, -50 },
    { -20, -20, -20, -90, -20 },
    { -50, -50, -90, -50, -50 },
    { -20, -60, -20, -50, -20 }},
     {/* AU.. */
      { 0, 0, 0, 0, 0 },
    { 0, -60, -50, -60, 0 },
    { 0, -40, -20, 0, -30 },
    { 0, -50, 0, -40, -60 },
    { 0, 0, -50, -50, -30 }},
     {/* UA.. */
      { 0, 0, 0, 0, 0 },
    { 0, -70, -60, -60, 0 },
    { 0, -30, -20, 0, -30 },
    { 0, -50, 0, -40, -50 },
    { 0, 0, -30, -50, -40 }},
     {/* NN.. */
      { 0, 0, 0, 0, 0 },
    { 0, -50, -20, -50, 0 },
    { 0, -20, -20, 0, -20 },
    { 0, -50, 0, -40, -50 },
    { 0, 0, -20, -50, -20 }}};

/* mismatch_exterior_enthalpies */
PUBLIC int mismatchExtdH[NBPAIRS + 1][5][5] =
    {{/* NP.. */
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF},
      {INF, INF, INF, INF, INF}},
     {/* CG.. */
      { -90, -90, -90, -90, -90 },
    { -620, -620, -620, -620, -620 },
    { -330, -330, -330, -330, -330 },
    { -510, -510, -510, -510, -510 },
    { -90, -90, -90, -90, -90 }},
     {/* GC.. */
      { -420, -420, -420, -420, -420 },
    { -460, -460, -460, -460, -460 },
    { -420, -420, -420, -420, -420 },
    { -470, -470, -470, -470, -470 },
    { -500, -500, -500, -500, -500 }},
     {/* GU.. */
      { -60, -60, -60, -60, -60 },
    { -90, -90, -90, -90, -90 },
    { -60, -60, -60, -60, -60 },
    { -250, -250, -250, -250, -250 },
    { -900, -900, -900, -900, -900 }},
     {/* UG.. */
      { -90, -90, -90, -90, -90 },
    { -620, -620, -620, -620, -620 },
    { -330, -330, -330, -330, -330 },
    { -510, -510, -510, -510, -510 },
    { -90, -90, -90, -90, -90 }},
     {/* AU.. */
      { -60, -60, -60, -60, -60 },
    { -90, -90, -90, -90, -90 },
    { -60, -60, -60, -60, -60 },
    { -250, -250, -250, -250, -250 },
    { -900, -900, -900, -900, -900 }},
     {/* UA.. */
      { 200, 200, 200, 200, 200 },
    { 80, 80, 80, 80, 80 },
    { 200, 200, 200, 200, 200 },
    { 150, 150, 150, 150, 150 },
    { 110, 110, 110, 110, 110 }},
     {/* NN.. */
      { 200, 200, 200, 200, 200 },
    { 80, 80, 80, 80, 80 },
    { 200, 200, 200, 200, 200 },
    { 150, 150, 150, 150, 150 },
    { 110, 110, 110, 110, 110 }}};

/* dangle5 */
PUBLIC int dangle5_37[NBPAIRS + 1][5] =
    {/*           N      A      C      G      U */
     /* NP */ {INF, INF, INF, INF, INF},
     /* CG */ { -50, -90, -50, -70, -60 },
    { -50, -80, -50, -80, -70 },
    { -10, -10, -10, -10, -10 },
    { -10, -10, -10, -10, -10 },
    { -20, -50, -20, -60, -30 },
    { -20, -50, -20, -50, -30 },
    { -10, -10, -10, -10, -10 }};

/* dangle3 */
PUBLIC int dangle3_37[NBPAIRS + 1][5] =
    {/*           N      A      C      G      U */
     /* NP */ {INF, INF, INF, INF, INF},
     /* CG */ { -20, -40, -20, -40, -30 },
    { -40, -110, -40, -110, -80 },
    { -20, -20, -20, -20, -20 },
    { -20, -20, -20, -20, -20 },
    { -20, -20, -20, -20, -20 },
    { -20, -20, -20, -20, -20 },
    { -20, -20, -20, -20, -20 }};

/* dangle5_enthalpies */
PUBLIC int dangle5_dH[NBPAIRS + 1][5] =
    {/*           N      A      C      G      U */
     /* NP */ {INF, INF, INF, INF, INF},
     /* CG */ { -90, -620, -330, -510, -90 },
    { -420, -460, -420, -470, -500 },
    { -420, -460, -420, -470, -500 },
    { -90, -620, -330, -510, -90 },
    { -60, -90, -60, -250, -900 },
    { 200, 80, 200, 150, 110 },
    { 200, 80, 200, 150, 110 }};

/* dangle3_enthalpies */
PUBLIC int dangle3_dH[NBPAIRS + 1][5] =
    {/*           N      A      C      G      U */
     /* NP */ {INF, INF, INF, INF, INF},
     /* CG */ { -20, -190, -20, -400, -550 },
    { -250, -800, -250, -310, -390 },
    { -250, -800, -250, -310, -390 },
    { -20, -190, -20, -400, -550 },
    { 200, -530, 200, -480, 90 },
    { 630, 430, 630, 250, 260 },
    { 630, 430, 630, 250, 260 }};



// PUBLIC char Hexaloops[361] =
//     "ACAGUACU "
//     "ACAGUGAU "
//     "ACAGUGCU "
//     "ACAGUGUU ";
// PUBLIC int Hexaloop37[40] = {280, 360, 290, 180};
// PUBLIC int HexaloopdH[40] = {-1680, -1140, -1280, -1540};
PUBLIC char Hexaloops[361];
PUBLIC int Hexaloop37[40];
PUBLIC int HexaloopdH[40];

// PUBLIC char Tetraloops[281] =
//     "CAACGG "
//     "CCAAGG "
//     "CCACGG "
//     "CCCAGG "
//     "CCGAGG "
//     "CCGCGG "
//     "CCUAGG "
//     "CCUCGG "
//     "CUAAGG "
//     "CUACGG "
//     "CUCAGG "
//     "CUCCGG "
//     "CUGCGG "
//     "CUUAGG "
//     "CUUCGG "
//     "CUUUGG ";
// PUBLIC int Tetraloop37[40] = {550, 330, 370, 340, 350, 360, 370, 250, 360, 280, 370, 270, 280, 350, 370, 370};
// PUBLIC int TetraloopdH[40] = {690, -1030, -330, -890, -660, -750, -350, -1390, -760, -1070, -660, -1290, -1070, -620, -1530, -680};
PUBLIC char Tetraloops[281];
PUBLIC int Tetraloop37[40];
PUBLIC int TetraloopdH[40];

// PUBLIC char Triloops[241] =
//     "CAACG "
//     "GUUAC ";
// PUBLIC int Triloop37[40] = {680, 690};
// PUBLIC int TriloopdH[40] = {2370, 1080};
PUBLIC char Triloops[241];
PUBLIC int Triloop37[40] ;
PUBLIC int TriloopdH[40] ;


#define saltT md->temperature + K0
#define NBASES 8
/*@notnull@*/

#ifndef INLINE
#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
#endif
#endif



#define MAXALPHA 20 /* maximal length of alphabet */

static short alias[MAXALPHA + 1];
static int pair[MAXALPHA + 1][MAXALPHA + 1];
/* rtype[pair[i][j]]:=pair[j][i] */
// static int    rtype[8] = {
//   0, 2, 1, 4, 3, 6, 5, 7
// };

#ifdef _OPENMP
#pragma omp threadprivate(Law_and_Order, BP_pair, alias, pair, rtype)
#endif

#define Rods_dist 20.
/* for backward compatibility */
#define ENCODE(c) encode_char(c)
#define Eular_const 0.58
#define BUFFER_SIZE 1024
#define roundint(x) ((int)(x + 0.5 - (x < 0)))
#define MIN2(A, B) ((A) < (B) ? (A) : (B))
#define MAX2(A, B) ((A) > (B) ? (A) : (B))

#define EINVAL 22
#define ENOMEM 12

// constrain

#define STATE_UNINITIALIZED (unsigned char)4
#define VRNA_CONSTRAINT_CONTEXT_EXT_LOOP      (unsigned char)0x01

/**
 *  @brief  Hard constraints flag, base pair encloses hairpin loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_HP_LOOP       (unsigned char)0x02

/**
 *  @brief  Hard constraints flag, base pair encloses an interior loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_INT_LOOP      (unsigned char)0x04

/**
 *  @brief  Hard constraints flag, base pair encloses a multi branch loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC  (unsigned char)0x08

/**
 *  @brief  Hard constraints flag, base pair is enclosed in an interior loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_MB_LOOP       (unsigned char)0x10

/**
 *  @brief  Hard constraints flag, base pair is enclosed in a multi branch loop
 *
 *  @ingroup  hard_constraints
 *
 */
#define VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC   (unsigned char)0x20

/**
 *  @brief  Hard constraint flag to indicate enforcement of constraints
 */
#define VRNA_CONSTRAINT_CONTEXT_ENFORCE       (unsigned char)0x40

/**
 *  @brief  Hard constraint flag to indicate not to remove base pairs that conflict with a given constraint
 */
#define VRNA_CONSTRAINT_CONTEXT_NO_REMOVE     (unsigned char)0x80


/**
 *  @brief  Constraint context flag that forbids a nucleotide or base pair to appear in any loop
 */
#define VRNA_CONSTRAINT_CONTEXT_NONE          (unsigned char)0
#define VRNA_CONSTRAINT_CONTEXT_CLOSING_LOOPS (unsigned char)(VRNA_CONSTRAINT_CONTEXT_EXT_LOOP | VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_INT_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP)
#define VRNA_CONSTRAINT_CONTEXT_ENCLOSED_LOOPS  (unsigned char)(VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC | \
                                                                VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
#define VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS (unsigned char)(VRNA_CONSTRAINT_CONTEXT_CLOSING_LOOPS | VRNA_CONSTRAINT_CONTEXT_ENCLOSED_LOOPS)

#define VRNA_OPTION_DEFAULT 0U
#define VRNA_OPTION_HYBRID (1 << 2)

#define VRNA_OPTION_WINDOW (1 << 4)
#define VRNA_OPTION_EVAL_ONLY (1 << 3)

#define WITH_PTYPE 1L
#define WITH_PTYPE_COMPAT 2L
#define VRNA_OPTION_EVAL_ONLY (1 << 3)
#define SCALE 10
#define CLIP_NEGATIVE(X) ((X) < 0 ? 0 : (X))
#define TRUNC_MAYBE(X) ((!pf_smooth) ? (double)((int)(X)) : (X))

#define SMOOTH(X) ((!pf_smooth)               ?   CLIP_NEGATIVE(X) : \
                   ((X) / SCALE < -1.2283697) ?                 0  : \
                   ((X) / SCALE > 0.8660254) ?                (X) : \
                   SCALE *0.38490018   \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   )
/* Rescale Free energy contribution according to deviation of temperature from measurement conditions */
#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)
#define RESCALE_BF(dG, dH, dT, kT)          ( \
    exp( \
      -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
      * 10. \
      / kT \
      ) \
    )

#define RESCALE_BF_SMOOTH(dG, dH, dT, kT)   ( \
    exp(  \
      SMOOTH( \
        -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
        ) \
      * 10. \
      / kT \
      ) \
    )

