#define PUBLIC
#define PRIVATE static
#define INLINE inline
#define INF 10000000
#define MAXALPHA 20       /* maximal length of alphabet */
#define MAXLOOP 30
#define NBPAIRS 7
#define VRNA_GQUAD_MAX_STACK_SIZE 7
#define VRNA_GQUAD_MAX_LINKER_LENGTH 15

//model
#define VRNA_MODEL_DEFAULT_TEMPERATURE    37.0
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
#define TURN 3
#define BP_REV_DEFAULT { 0, 2, 1, 4, 3, 6, 5, 7 }
#define BP_ALIAS_DEFAULT { 0, 1, 2, 3, 4, 3, 2, 0 }
#define BP_ENCODING_DEFAULT { { 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 5, 0, 0, 5 }, { 0, 0, 0, 1, 0, 0, 0, 0 }, { 0, 0, 2, 0, 3, 0, 0, 0 }, { 0, 6, 0, 4, 0, 0, 0, 6 }, { 0, 0, 0, 0, 0, 0, 2, 0 }, { 0, 0, 0, 0, 0, 1, 0, 0 }, { 0, 6, 0, 0, 5, 0, 0, 0 } }
#define DM_DEFAULT { { 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 2, 2, 1, 2, 2 } , { 0, 2, 0, 1, 2, 2, 2 } , { 0, 2, 1, 0, 2, 1, 2 } , { 0, 1, 2, 2, 0, 2, 1 } , { 0, 2, 2, 1, 2, 0, 2 } , { 0, 2, 2, 2, 1, 2, 0 } }

#define NBASES 8
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
  /* compute energy of degree 2 loop (stack bulge or interior) */
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

/**
 *  @brief The datastructure that contains temperature scaled energy parameters.
 */
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

/**
 *  @brief Typename for the fold_compound data structure #vrna_fc_s
 */
typedef struct vrna_fc_s vrna_fold_compound_t;

/** @brief Typename for the free energy parameter data structure #vrna_params */
typedef struct  vrna_param_s vrna_param_t;

/** @brief Typename for the model details data structure #vrna_md_s */
typedef struct vrna_md_s vrna_md_t;

/** @brief Typename for the soft constraints data structure #vrna_sc_s
 *  @ingroup  soft_constraints
 */
typedef struct  vrna_sc_s vrna_sc_t;

/**
 * @brief Callback to evaluate whether or not a particular decomposition step is contributing to the solution space
 *
 * @ingroup hard_constraints
 *
 * This is the prototype for callback functions used by the folding recursions to evaluate generic
 * hard constraints. The first four parameters passed indicate the delimiting nucleotide positions
 * of the decomposition, and the parameter @p denotes the decomposition step. The last parameter
 * @p data is the auxiliary data structure associated to the hard constraints via vrna_hc_add_data(),
 * or NULL if no auxiliary data was added.
 *
 * @callback
 * @parblock
 * This callback enables one to over-rule default hard constraints in secondary structure
 * decompositions.
 * @endparblock
 *
 * @see #VRNA_DECOMP_PAIR_HP, #VRNA_DECOMP_PAIR_IL, #VRNA_DECOMP_PAIR_ML, #VRNA_DECOMP_ML_ML_ML,
 *      #VRNA_DECOMP_ML_STEM, #VRNA_DECOMP_ML_ML, #VRNA_DECOMP_ML_UP, #VRNA_DECOMP_ML_ML_STEM,
 *      #VRNA_DECOMP_ML_COAXIAL, #VRNA_DECOMP_EXT_EXT, #VRNA_DECOMP_EXT_UP, #VRNA_DECOMP_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_EXT, #VRNA_DECOMP_EXT_STEM_EXT, #VRNA_DECOMP_EXT_EXT_STEM,
 *      #VRNA_DECOMP_EXT_EXT_STEM1, vrna_hc_add_f(), vrna_hc_add_data()
 *
 * @param i         Left (5') delimiter position of substructure
 * @param j         Right (3') delimiter position of substructure
 * @param k         Left delimiter of decomposition
 * @param l         Right delimiter of decomposition
 * @param d         Decomposition step indicator
 * @param data      Auxiliary data
 * @return          A non-zero value if the decomposition is valid, 0 otherwise
 */
typedef unsigned char (*vrna_hc_eval_f)(int           i,
                                        int           j,
                                        int           k,
                                        int           l,
                                        unsigned char d,
                                        void          *data);

struct hc_ext_def_dat {
  unsigned int              n;
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};

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

#define NULL ((void *)0)

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
