#define PUBLIC
#define PRIVATE static
#define MAXALPHA 20       /* maximal length of alphabet */
#define VRNA_OPTION_DEFAULT 0U
#define VRNA_OPTION_HYBRID (1 << 2)

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

struct options {
  int             filename_full;
  char            *filename_delim;
  int             pf;
  int             doT;
  int             doC;
  int             noPS;
  int             noconv;
  int             centroid;
  int             MEA;
  double          MEAgamma;
  double          bppmThreshold;
  int             verbose;
  vrna_md_t       md;
  // vrna_cmd_t      commands;

//   dataset_id      id_control;

  char            *concentration_file;

  char            *constraint_file;
  int             constraint_batch;
  int             constraint_enforce;
  int             constraint_canonical;

  int             shape;
  char            *shape_file;
  char            *shape_method;
  char            *shape_conversion;

//   vrna_sc_mod_param_t *mod_params;

  int             csv_output;
  int             csv_header;
  char            csv_output_delim;

  int             jobs;
  int             keep_order;
  unsigned int    next_record_number;
//   vrna_ostream_t  output_queue;
};


struct vrna_dimer_pf_s {
  /* free energies for: */
  double  F0AB; /**< @brief Null model without DuplexInit */
  double  FAB;  /**< @brief all states with DuplexInit correction */
  double  FcAB; /**< @brief true hybrid states only */
  double  FA;   /**< @brief monomer A */
  double  FB;   /**< @brief monomer B */
};

typedef struct vrna_dimer_pf_s vrna_dimer_pf_t;


typedef enum {
  VRNA_FC_TYPE_SINGLE = 0,      /**< Type is suitable for single, and hybridizing sequences */
  VRNA_FC_TYPE_COMPARATIVE = 1 /**< Type is suitable for sequence alignments (consensus structure prediction) */
} vrna_fc_type_e;

typedef enum {
  VRNA_SEQ_UNKNOWN = 0,   /**< @brief Nucleotide sequence represents an Unkown type */
  VRNA_SEQ_RNA = 1,       /**< @brief Nucleotide sequence represents an RNA type */
  VRNA_SEQ_DNA = 2        /**< @brief Nucleotide sequence represents a DNA type */
} vrna_seq_type_e;


struct vrna_sequence_s {
  vrna_seq_type_e type;       /**< @brief The type of sequence */
  char            *name;
  char            *string;    /**< @brief The string representation of the sequence */
  short           *encoding;  /**< @brief The integer representation of the sequence */
  short           *encoding5;
  short           *encoding3;
  unsigned int    length;     /**< @brief The length of the sequence */
};


typedef struct vrna_sequence_s vrna_seq_t;


typedef enum {
  VRNA_HC_DEFAULT = 0,  /**<  @brief  Default Hard Constraints */
  VRNA_HC_WINDOW  = 1  /**<  @brief  Hard Constraints suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
} vrna_hc_type_e;

typedef unsigned char (*vrna_hc_eval_f)(int           i,
                                        int           j,
                                        int           k,
                                        int           l,
                                        unsigned char d,
                                        void          *data);

typedef void (*vrna_auxdata_free_f)(void *data);

struct hc_nuc {
  int           direction;
  unsigned char context;
  unsigned char nonspec;
};

struct hc_basepair {
  size_t        list_size;
  size_t        list_mem;
  unsigned int  *j;
  unsigned int  *strand_j;
  unsigned char *context;
};

struct vrna_hc_depot_s {
  unsigned int            strands;
  size_t                  *up_size;
  struct hc_nuc           **up;
  size_t                  *bp_size;
  struct hc_basepair      **bp;
};

typedef struct vrna_hc_depot_s vrna_hc_depot_t;

struct vrna_hc_s {
  vrna_hc_type_e  type;
  unsigned int    n;

  unsigned char   state;

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  unsigned char *mx;
#ifndef VRNA_DISABLE_C11_FEATURES
};
struct {
#endif
  unsigned char **matrix_local;
#ifndef VRNA_DISABLE_C11_FEATURES
};
};
#endif

  int                 *up_ext;    /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in an exterior loop
                                   */
  int                 *up_hp;     /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in a hairpin loop
                                   */
  int                 *up_int;    /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in an interior loop
                                   */
  int                 *up_ml;     /**<  @brief  A linear array that holds the number of allowed
                                   *            unpaired nucleotides in a multi branched loop
                                   */

  vrna_hc_eval_f      f;          /**<  @brief  A function pointer that returns whether or
                                   *            not a certain decomposition may be evaluated
                                   */

  void                *data;      /**<  @brief  A pointer to some structure where the user
                                   *            may store necessary data to evaluate its
                                   *            generic hard constraint function
                                   */

  vrna_auxdata_free_f free_data;  /**<  @brief  A pointer to a function to free memory
                                   *            occupied by auxiliary data
                                   *
                                   *    The function this pointer is pointing to will be
                                   *    called upon destruction of the #vrna_hc_s, and
                                   *    provided with the vrna_hc_s.data pointer that
                                   *    may hold auxiliary data. Hence, to avoid leaking
                                   *    memory, the user may use this pointer to free
                                   *    memory occupied by auxiliary data.
                                   */

  vrna_hc_depot_t *depot;
};


typedef struct  vrna_hc_s vrna_hc_t;

typedef enum {
  VRNA_MX_DEFAULT = 0,  /**<  @brief  Default DP matrices */
  VRNA_MX_WINDOW = 1,   /**<  @brief  DP matrices suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
  VRNA_MX_2DFOLD = 2    /**<  @brief  DP matrices suitable for distance class partitioned structure prediction
                     *    @see  vrna_mfe_TwoD(), vrna_pf_TwoD()
                     */
} vrna_mx_type_e;

struct vrna_mx_mfe_s {
  /** @name Common fields for MFE matrices
   *  @{
   */
  const vrna_mx_type_e  type;     /**< Type of the DP matrices */
  unsigned int          length;   /**<  @brief  Length of the sequence, therefore an indicator of the size of the DP matrices */
  unsigned int          strands;  /**< Number of strands */
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  /** @name Default DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_DEFAULT @endcode
   * @{
   */
  int *c;           /**<  @brief  Energy array, given that i-j pair */
  int *f5;          /**<  @brief  Energy of 5' end */
  int *f3;          /**<  @brief  Energy of 3' end */
  int **fms5;       /**<  @brief  Energy for connected interstrand configurations */
  int **fms3;       /**<  @brief  nergy for connected interstrand configurations */
  int *fML;         /**<  @brief  Multi-loop auxiliary energy array */
  int *fM1;         /**<  @brief  Second ML array, only for unique multibrnach loop decomposition */
  int *fM2;         /**<  @brief  Energy for a multibranch loop region with exactly two stems, extending to 3' end */
  int *ggg;         /**<  @brief  Energies of g-quadruplexes */
  int Fc;           /**<  @brief  Minimum Free Energy of entire circular RNA */
  int FcH;          /**<  @brief  Minimum Free Energy of hairpin loop cases in circular RNA */
  int FcI;          /**<  @brief  Minimum Free Energy of internal loop cases in circular RNA */
  int FcM;          /**<  @brief  Minimum Free Energy of multibranch loop cases in circular RNA */
  /**
   * @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif
  /** @name Local Folding DP matrices using window approach
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
   * @{
   */
  int **c_local;            /**<  @brief  Energy array, given that i-j pair */
  int *f3_local;            /**<  @brief  Energy of 5' end */
  int **fML_local;          /**<  @brief  Multi-loop auxiliary energy array */
  int **ggg_local;          /**<  @brief  Energies of g-quadruplexes */
  /**
   * @}
   */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Distance Class DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_2DFOLD @endcode
   * @{
   */
  int           ***E_F5;
  int           **l_min_F5;
  int           **l_max_F5;
  int           *k_min_F5;
  int           *k_max_F5;

  int           ***E_F3;
  int           **l_min_F3;
  int           **l_max_F3;
  int           *k_min_F3;
  int           *k_max_F3;

  int           ***E_C;
  int           **l_min_C;
  int           **l_max_C;
  int           *k_min_C;
  int           *k_max_C;

  int           ***E_M;
  int           **l_min_M;
  int           **l_max_M;
  int           *k_min_M;
  int           *k_max_M;

  int           ***E_M1;
  int           **l_min_M1;
  int           **l_max_M1;
  int           *k_min_M1;
  int           *k_max_M1;

  int           ***E_M2;
  int           **l_min_M2;
  int           **l_max_M2;
  int           *k_min_M2;
  int           *k_max_M2;

  int           **E_Fc;
  int           *l_min_Fc;
  int           *l_max_Fc;
  int           k_min_Fc;
  int           k_max_Fc;

  int           **E_FcH;
  int           *l_min_FcH;
  int           *l_max_FcH;
  int           k_min_FcH;
  int           k_max_FcH;

  int           **E_FcI;
  int           *l_min_FcI;
  int           *l_max_FcI;
  int           k_min_FcI;
  int           k_max_FcI;

  int           **E_FcM;
  int           *l_min_FcM;
  int           *l_max_FcM;
  int           k_min_FcM;
  int           k_max_FcM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  int           *E_F5_rem;
  int           *E_F3_rem;
  int           *E_C_rem;
  int           *E_M_rem;
  int           *E_M1_rem;
  int           *E_M2_rem;

  int           E_Fc_rem;
  int           E_FcH_rem;
  int           E_FcI_rem;
  int           E_FcM_rem;

#ifdef COUNT_STATES
  unsigned long ***N_F5;
  unsigned long ***N_C;
  unsigned long ***N_M;
  unsigned long ***N_M1;
#endif

  /**
   * @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif
};



typedef struct  vrna_mx_mfe_s vrna_mx_mfe_t;


typedef double FLT_OR_DBL;

struct vrna_mx_pf_s {
  /** @name Common fields for DP matrices
   *  @{
   */
  const vrna_mx_type_e  type;       /**< Type of the DP matrices */
  unsigned int          length;     /**< Size of the DP matrices (i.e. sequence length) */
  FLT_OR_DBL            *scale;     /**< Boltzmann factor scaling */
  FLT_OR_DBL            *expMLbase; /**< Boltzmann factors for unpaired bases in multibranch loop */

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif

  /** @name Default PF matrices
   *  @note These data fields are available if
   *        @code vrna_mx_pf_t.type == VRNA_MX_DEFAULT @endcode
   *  @{
   */
  FLT_OR_DBL *q;
  FLT_OR_DBL *qb;
  FLT_OR_DBL *qm;
  FLT_OR_DBL *qm1;
  FLT_OR_DBL *probs;
  FLT_OR_DBL *q1k;
  FLT_OR_DBL *qln;
  FLT_OR_DBL *G;

  FLT_OR_DBL qo;
  FLT_OR_DBL *qm2;
  FLT_OR_DBL qho;
  FLT_OR_DBL qio;
  FLT_OR_DBL qmo;

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Local Folding DP matrices using window approach
   *  @note These data fields are available if
   *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
   * @{
   */
  FLT_OR_DBL **q_local;
  FLT_OR_DBL **qb_local;
  FLT_OR_DBL **qm_local;
  FLT_OR_DBL **pR;
  FLT_OR_DBL **qm2_local;
  FLT_OR_DBL **QI5;
  FLT_OR_DBL **q2l;
  FLT_OR_DBL **qmb;
  FLT_OR_DBL **G_local;
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /** @name Distance Class DP matrices
   *  @note These data fields are available if
   *        @code vrna_mx_pf_t.type == VRNA_MX_2DFOLD @endcode
   *  @{
   */
  FLT_OR_DBL ***Q;
  int **l_min_Q;
  int **l_max_Q;
  int *k_min_Q;
  int *k_max_Q;


  FLT_OR_DBL ***Q_B;
  int **l_min_Q_B;
  int **l_max_Q_B;
  int *k_min_Q_B;
  int *k_max_Q_B;

  FLT_OR_DBL ***Q_M;
  int **l_min_Q_M;
  int **l_max_Q_M;
  int *k_min_Q_M;
  int *k_max_Q_M;

  FLT_OR_DBL ***Q_M1;
  int **l_min_Q_M1;
  int **l_max_Q_M1;
  int *k_min_Q_M1;
  int *k_max_Q_M1;

  FLT_OR_DBL ***Q_M2;
  int **l_min_Q_M2;
  int **l_max_Q_M2;
  int *k_min_Q_M2;
  int *k_max_Q_M2;

  FLT_OR_DBL **Q_c;
  int *l_min_Q_c;
  int *l_max_Q_c;
  int k_min_Q_c;
  int k_max_Q_c;

  FLT_OR_DBL **Q_cH;
  int *l_min_Q_cH;
  int *l_max_Q_cH;
  int k_min_Q_cH;
  int k_max_Q_cH;

  FLT_OR_DBL **Q_cI;
  int *l_min_Q_cI;
  int *l_max_Q_cI;
  int k_min_Q_cI;
  int k_max_Q_cI;

  FLT_OR_DBL **Q_cM;
  int *l_min_Q_cM;
  int *l_max_Q_cM;
  int k_min_Q_cM;
  int k_max_Q_cM;

  /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
  FLT_OR_DBL *Q_rem;
  FLT_OR_DBL *Q_B_rem;
  FLT_OR_DBL *Q_M_rem;
  FLT_OR_DBL *Q_M1_rem;
  FLT_OR_DBL *Q_M2_rem;

  FLT_OR_DBL Q_c_rem;
  FLT_OR_DBL Q_cH_rem;
  FLT_OR_DBL Q_cI_rem;
  FLT_OR_DBL Q_cM_rem;
  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif
};

typedef struct  vrna_mx_pf_s vrna_mx_pf_t;



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

#define NBPAIRS 7
#define MAXLOOP 30
#define VRNA_GQUAD_MAX_STACK_SIZE 7
#define VRNA_GQUAD_MAX_LINKER_LENGTH 15

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



struct vrna_exp_param_s {
  int     id;   /**<  @brief  An identifier for the data structure
                 *    @deprecated This attribute will be removed in version 3
                 */
  double  expstack[NBPAIRS + 1][NBPAIRS + 1];
  double  exphairpin[31];
  double  expbulge[MAXLOOP + 1];
  double  expinternal[MAXLOOP + 1];
  double  expmismatchExt[NBPAIRS + 1][5][5];
  double  expmismatchI[NBPAIRS + 1][5][5];
  double  expmismatch23I[NBPAIRS + 1][5][5];
  double  expmismatch1nI[NBPAIRS + 1][5][5];
  double  expmismatchH[NBPAIRS + 1][5][5];
  double  expmismatchM[NBPAIRS + 1][5][5];
  double  expdangle5[NBPAIRS + 1][5];
  double  expdangle3[NBPAIRS + 1][5];
  double  expint11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  double  expint21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  double  expint22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  double  expninio[5][MAXLOOP + 1];
  double  lxc;
  double  expMLbase;
  double  expMLintern[NBPAIRS + 1];
  double  expMLclosing;
  double  expTermAU;
  double  expDuplexInit;
  double  exptetra[40];
  double  exptri[40];
  double  exphex[40];
  char    Tetraloops[1401];
  double  expTriloop[40];
  char    Triloops[241];
  char    Hexaloops[1801];
  double  expTripleC;
  double  expMultipleCA;
  double  expMultipleCB;
  double  expgquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  double  expgquadLayerMismatch;
  int     gquadLayerMismatchMax;

  double  kT;
  double  pf_scale;           /**<  @brief    Scaling factor to avoid over-/underflows */

  double  temperature;        /**<  @brief    Temperature used for loop contribution scaling */
  double  alpha;              /**<  @brief    Scaling factor for the thermodynamic temperature
                               *    @details  This allows for temperature scaling in Boltzmann
                               *              factors independently from the energy contributions.
                               *              The resulting Boltzmann factors are then computed by
                               *              @f$ e^{-E/(\alpha \cdot K \cdot T)} @f$
                               */

  vrna_md_t model_details;    /**<  @brief  Model details to be used in the recursions */
  char      param_file[256];  /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */

  double    expSaltStack;
  double    expSaltLoop[MAXLOOP + 2];
  double    SaltLoopDbl[MAXLOOP + 2];
  int       SaltMLbase;
  int       SaltMLintern;
  int       SaltMLclosing;
  int       SaltDPXInit;
};


typedef struct  vrna_exp_param_s vrna_exp_param_t;

typedef void (*vrna_recursion_status_f)(unsigned char status,
                                              void          *data);


typedef void (*vrna_ud_production_f)(vrna_fold_compound_t *fc,
                                           void                 *data);
typedef void (*vrna_ud_exp_production_f)(vrna_fold_compound_t *fc,
                                               void                 *data);

typedef int (*vrna_ud_f)(vrna_fold_compound_t  *fc,
                                      int                   i,
                                      int                   j,
                                      unsigned int          loop_type,
                                      void                  *data);

typedef FLT_OR_DBL (*vrna_ud_exp_f)(vrna_fold_compound_t *fc,
                                                 int                  i,
                                                 int                  j,
                                                 unsigned int         loop_type,
                                                 void                 *data);

typedef void (*vrna_ud_add_probs_f)(vrna_fold_compound_t  *fc,
                                          int                   i,
                                          int                   j,
                                          unsigned int          loop_type,
                                          FLT_OR_DBL            exp_energy,
                                          void                  *data);


typedef FLT_OR_DBL (*vrna_ud_get_probs_f)(vrna_fold_compound_t  *fc,
                                                int                   i,
                                                int                   j,
                                                unsigned int          loop_type,
                                                int                   motif,
                                                void                  *data);


struct vrna_unstructured_domain_s {
  /*
   **********************************
   * Keep track of all motifs added
   **********************************
   */
  int           uniq_motif_count;                   /**<  @brief The unique number of motifs of different lengths */
  unsigned int  *uniq_motif_size;                   /**<  @brief An array storing a unique list of motif lengths */

  int           motif_count;                        /**<  @brief Total number of distinguished motifs */
  char          **motif;                            /**<  @brief Motif sequences */
  char          **motif_name;                       /**<  @brief Motif identifier/name */
  unsigned int  *motif_size;                        /**<  @brief Motif lengths */
  double        *motif_en;                          /**<  @brief Ligand binding free energy contribution */
  unsigned int  *motif_type;                        /**<  @brief Type of motif, i.e. loop type the ligand binds to */

  /*
   **********************************
   * Grammar extension for ligand
   * binding
   **********************************
   */
  vrna_ud_production_f     prod_cb;        /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP free energy matrices
                                                   *    @details This callback will be executed right before the actual secondary structure decompositions,
                                                   *    and, therefore, any implementation must not interleave with the regular DP matrices.
                                                   */
  vrna_ud_exp_production_f exp_prod_cb;    /**<  @brief Callback to ligand binding production rule, i.e. create/fill DP partition function matrices */
  vrna_ud_f         energy_cb;      /**<  @brief Callback to evaluate free energy of ligand binding to a particular unpaired stretch */
  vrna_ud_exp_f     exp_energy_cb;  /**<  @brief Callback to evaluate Boltzmann factor of ligand binding to a particular unpaired stretch */
  void                            *data;          /**<  @brief Auxiliary data structure passed to energy evaluation callbacks */
  vrna_auxdata_free_f      free_data;      /**<  @brief Callback to free auxiliary data structure */
  vrna_ud_add_probs_f      probs_add;      /**<  @brief Callback to store/add outside partition function */
  vrna_ud_get_probs_f      probs_get;      /**<  @brief Callback to retrieve outside partition function */
};


typedef struct vrna_unstructured_domain_s vrna_ud_t;



typedef void (*vrna_grammar_cond_f)(vrna_fold_compound_t *fc,
                                     unsigned char        stage,
                                     void                 *data);

typedef int (*vrna_grammar_rule_f)(vrna_fold_compound_t  *fc,
                                    int                   i,
                                    int                   j,
                                    void                  *data);


typedef void (*vrna_grammar_rule_f_aux)(vrna_fold_compound_t  *fc,
                                    int                   i,
                                    int                   j,
                                    void                  *data);


typedef FLT_OR_DBL (*vrna_grammar_rule_f_exp)(vrna_fold_compound_t *fc,
                                               int                  i,
                                               int                  j,
                                               void                 *data);


typedef void (*vrna_grammar_rule_f_aux_exp)(vrna_fold_compound_t *fc,
                                               int                  i,
                                               int                  j,
                                               void                 *data);


typedef void (*vrna_grammar_data_free_f)(void *data);



struct vrna_gr_aux_s {
  vrna_grammar_cond_f       cb_proc; /**< @brief A callback for pre- and post-processing of auxiliary grammar rules */

  vrna_grammar_rule_f       cb_aux_f;
  vrna_grammar_rule_f       cb_aux_c;
  vrna_grammar_rule_f       cb_aux_m;
  vrna_grammar_rule_f       cb_aux_m1;
  vrna_grammar_rule_f_aux       cb_aux;

  vrna_grammar_rule_f_exp   cb_aux_exp_f;
  vrna_grammar_rule_f_exp   cb_aux_exp_c;
  vrna_grammar_rule_f_exp     cb_aux_exp_m;
  vrna_grammar_rule_f_exp     cb_aux_exp_m1;
  vrna_grammar_rule_f_aux_exp   cb_aux_exp;

  void                        *data;
  vrna_grammar_data_free_f  free_data;
};


typedef struct vrna_gr_aux_s vrna_gr_aux_t;


typedef enum {
  VRNA_SC_DEFAULT = 0,  /**<  @brief  Default Soft Constraints */
  VRNA_SC_WINDOW = 1   /**<  @brief  Soft Constraints suitable for local structure prediction using
                     *    window approach.
                     *    @see    vrna_mfe_window(), vrna_mfe_window_zscore(), pfl_fold()
                     */
} vrna_sc_type_e;

typedef struct {
  unsigned int  interval_start;
  unsigned int  interval_end;
  int           e;
} vrna_sc_bp_storage_t;

typedef int (*vrna_sc_f)(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data);

struct vrna_basepair_s {
  int i;
  int j;
};

typedef struct vrna_basepair_s vrna_basepair_t;

typedef vrna_basepair_t *(*vrna_sc_bt_f)(int            i,
                                         int            j,
                                         int            k,
                                         int            l,
                                         unsigned char  d,
                                         void           *data);

typedef FLT_OR_DBL (*vrna_sc_exp_f)(int           i,
                                    int           j,
                                    int           k,
                                    int           l,
                                    unsigned char d,
                                    void          *data);


typedef int (*vrna_auxdata_prepare_f)(vrna_fold_compound_t  *fc,
                                      void                  *data,
                                      unsigned int          event,
                                      void                  *event_data);




struct vrna_sc_s {
  const vrna_sc_type_e  type;
  unsigned int          n;

  unsigned char         state;

  int                   **energy_up;      /**<  @brief Energy contribution for stretches of unpaired nucleotides */
  FLT_OR_DBL            **exp_energy_up;  /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */

  int                   *up_storage;      /**<  @brief  Storage container for energy contributions per unpaired nucleotide */
  vrna_sc_bp_storage_t  **bp_storage;     /**<  @brief  Storage container for energy contributions per base pair */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif
  int *energy_bp;                               /**<  @brief Energy contribution for base pairs */
  FLT_OR_DBL *exp_energy_bp;                    /**<  @brief Boltzmann Factors of the energy contribution for base pairs */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif
  int         **energy_bp_local;                        /**<  @brief Energy contribution for base pairs (sliding window approach) */
  FLT_OR_DBL  **exp_energy_bp_local;                    /**<  @brief Boltzmann Factors of the energy contribution for base pairs (sliding window approach) */
#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
};
#endif

  int           *energy_stack;                    /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
  FLT_OR_DBL    *exp_energy_stack;                /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

  /* generic soft contraints below */
  vrna_sc_f     f;            /**<  @brief  A function pointer used for pseudo
                               *            energy contribution in MFE calculations
                               *    @see    vrna_sc_add_f()
                               */

  vrna_sc_bt_f  bt;           /**<  @brief  A function pointer used to obtain backtraced
                               *            base pairs in loop regions that were altered
                               *            by soft constrained pseudo energy contributions
                               *    @see    vrna_sc_add_bt()
                               */

  vrna_sc_exp_f exp_f;        /**<  @brief  A function pointer used for pseudo energy
                               *            contribution boltzmann factors in PF
                               *            calculations
                               *    @see    vrna_sc_add_exp_f()
                               */

  void                *data;  /**<  @brief  A pointer to the data object provided for
                               *            for pseudo energy contribution functions of the
                               *            generic soft constraints feature
                               */

  vrna_auxdata_prepare_f  prepare_data;
  vrna_auxdata_free_f     free_data;
};


typedef struct  vrna_sc_s vrna_sc_t;

#ifdef VRNA_WARN_DEPRECATED
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif


struct vrna_fc_s {
  /**
   *  @name Common data fields
   *  @{
   */
  const vrna_fc_type_e type;        /**<  @brief  The type of the #vrna_fold_compound_t.
                                     * @details Currently possible values are #VRNA_FC_TYPE_SINGLE, and #VRNA_FC_TYPE_COMPARATIVE
                                     * @warning Do not edit this attribute, it will be automagically set by
                                     *      the corresponding get() methods for the #vrna_fold_compound_t.
                                     *      The value specified in this attribute dictates the set of other
                                     *      attributes to use within this data structure.
                                     */
  unsigned int      length;         /**<  @brief  The length of the sequence (or sequence alignment) */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
  DEPRECATED(int cutpoint,
             "Use strand_* members instead");
                                    /**<  @brief  The position of the (cofold) cutpoint within the provided sequence.
                                     * If there is no cutpoint, this field will be set to -1
                                     */
#endif
  unsigned int      *strand_number; /**<  @brief  The strand number a particular nucleotide is associated with */
  unsigned int      *strand_order;  /**<  @brief  The strand order, i.e. permutation of current concatenated sequence */
  unsigned int      *strand_order_uniq; /**<  @brief  The strand order array where identical sequences have the same ID */
  unsigned int      *strand_start;  /**<  @brief  The start position of a particular strand within the current concatenated sequence */
  unsigned int      *strand_end;    /**<  @brief  The end (last) position of a particular strand within the current concatenated sequence */

  unsigned int      strands;        /**<  @brief  Number of interacting strands */
  vrna_seq_t        *nucleotides;   /**<  @brief  Set of nucleotide sequences */

  vrna_hc_t         *hc;            /**<  @brief  The hard constraints data structure used for structure prediction */

  vrna_mx_mfe_t     *matrices;      /**<  @brief  The MFE DP matrices */
  vrna_mx_pf_t      *exp_matrices;  /**<  @brief  The PF DP matrices  */

  vrna_param_t      *params;        /**<  @brief  The precomputed free energy contributions for each type of loop */
  vrna_exp_param_t  *exp_params;    /**<  @brief  The precomputed free energy contributions as Boltzmann factors  */

  int               *iindx;         /**<  @brief  DP matrix accessor  */
  int               *jindx;         /**<  @brief  DP matrix accessor  */

  /**
   *  @}
   *
   *  @name User-defined data fields
   *  @{
   */
  vrna_recursion_status_f   stat_cb;       /**<  @brief  Recursion status callback (usually called just before, and
                                            *            after recursive computations in the library
                                            *    @see    vrna_recursion_status_f(), vrna_fold_compound_add_callback()
                                            */

  void                      *auxdata;      /**<  @brief  A pointer to auxiliary, user-defined data
                                            *    @see vrna_fold_compound_add_auxdata(), #vrna_fold_compound_t.free_auxdata
                                            */

  vrna_auxdata_free_f       free_auxdata;  /**<  @brief A callback to free auxiliary user data whenever the fold_compound itself is free'd
                                            *    @see  #vrna_fold_compound_t.auxdata, vrna_auxdata_free_f()
                                            */

  /**
   *  @}
   *
   *  @name Secondary Structure Decomposition (grammar) related data fields
   *  @{
   */

  /* data structure to adjust additional structural domains, such as G-quadruplexes */
//   vrna_sd_t     *domains_struc;             /**<  @brief  Additional structured domains */

  /* data structure to adjust additional contributions to unpaired stretches, e.g. due to protein binding */
  vrna_ud_t     *domains_up;                /**<  @brief  Additional unstructured domains */

  /* auxiliary (user-defined) extension to the folding grammar */
  vrna_gr_aux_t *aux_grammar;               /**<  @brief  Additional decomposition grammar rules */

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
  union {
    struct {
#endif

  /**
   *  @name Data fields available for single/hybrid structure prediction
   *  @{
   */
      char *sequence;                   /**<  @brief  The input sequence string
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */
      short *sequence_encoding;         /**<  @brief  Numerical encoding of the input sequence
                                         *    @see    vrna_sequence_encode()
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */
      short *encoding5;
      short *encoding3;

      short *sequence_encoding2;


      char *ptype;                      /**<  @brief  Pair type array
                                         *
                                         *    Contains the numerical encoding of the pair type for each pair (i,j) used
                                         *    in MFE, Partition function and Evaluation computations.
                                         *    @note This array is always indexed via jindx, in contrast to previously
                                         *    different indexing between mfe and pf variants!
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         *    @see    vrna_idx_col_wise(), vrna_ptypes()
                                         */
      char *ptype_pf_compat;            /**<  @brief  ptype array indexed via iindx
                                         *    @deprecated  This attribute will vanish in the future!
                                         *    It's meant for backward compatibility only!
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */
      vrna_sc_t *sc;                    /**<  @brief  The soft constraints for usage in structure prediction and evaluation
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
                                         */

  /**
   *  @}
   */

#ifndef VRNA_DISABLE_C11_FEATURES
  /* C11 support for unnamed unions/structs */
};
struct {
#endif

  /**
   *  @name Data fields for consensus structure prediction
   *  @{
   */
      char          **sequences;        /**<  @brief  The aligned sequences
                                         *    @note   The end of the alignment is indicated by a NULL pointer in the second dimension
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      unsigned int  n_seq;              /**<  @brief  The number of sequences in the alignment
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      char          *cons_seq;          /**<  @brief  The consensus sequence of the aligned sequences
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         *S_cons;            /**<  @brief  Numerical encoding of the consensus sequence
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S;                /**<  @brief  Numerical encoding of the sequences in the alignment
                                         *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S5;               /**<  @brief    S5[s][i] holds next base 5' of i in sequence s
                                         *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
      short         **S3;               /**<  @brief    Sl[s][i] holds next base 3' of i in sequence s
                                         *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                         */
  char          **Ss;
  unsigned int  **a2s;
      int           *pscore;              /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      int           **pscore_local;       /**<  @brief  Precomputed array of pair types expressed as pairing scores
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      short         *pscore_pf_compat;    /**<  @brief  Precomputed array of pair types expressed as pairing scores indexed via iindx
                                           *    @deprecated  This attribute will vanish in the future!
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
      vrna_sc_t     **scs;                /**<  @brief  A set of soft constraints (for each sequence in the alignment)
                                           *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
                                           */
  int           oldAliEn;

  /**
   *  @}
   */
#ifndef VRNA_DISABLE_C11_FEATURES
};
};
#endif

  /**
   *  @name Additional data fields for Distance Class Partitioning
   *
   *  These data fields are typically populated with meaningful data only if used in the context of Distance Class Partitioning
   *  @{
   */
  unsigned int  maxD1;            /**<  @brief  Maximum allowed base pair distance to first reference */
  unsigned int  maxD2;            /**<  @brief  Maximum allowed base pair distance to second reference */
  short         *reference_pt1;   /**<  @brief  A pairtable of the first reference structure */
  short         *reference_pt2;   /**<  @brief  A pairtable of the second reference structure */

  unsigned int  *referenceBPs1;   /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
  unsigned int  *referenceBPs2;   /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
  unsigned int  *bpdist;          /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

  unsigned int  *mm1;             /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
  unsigned int  *mm2;             /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

  /**
   *  @name Additional data fields for local folding
   *
   *  These data fields are typically populated with meaningful data only if used in the context of local folding
   *  @{
   */
  int   window_size;              /**<  @brief  window size for local folding sliding window approach */
  char  **ptype_local;            /**<  @brief  Pair type array (for local folding) */

};



typedef struct vrna_fc_s vrna_fold_compound_t;


typedef enum {
  VRNA_FC_TYPE_SINGLE = 0,      /**< Type is suitable for single, and hybridizing sequences */
  VRNA_FC_TYPE_COMPARATIVE = 1 /**< Type is suitable for sequence alignments (consensus structure prediction) */
} vrna_fc_type_e;

struct vrna_elem_prob_s {
  int   i;    /**<  @brief  Start position (usually 5' nucleotide that starts the element, e.g. base pair) */
  int   j;    /**<  @brief  End position (usually 3' nucleotide that ends the element, e.g. base pair) */
  float p;    /**<  @brief  Probability of the element */
  int   type; /**<  @brief  Type of the element */
};

typedef struct vrna_elem_prob_s vrna_ep_t;


typedef struct {
  FLT_OR_DBL  *prm_l;
  FLT_OR_DBL  *prm_l1;
  FLT_OR_DBL  *prml;

  int         ud_max_size;
  FLT_OR_DBL  **pmlu;
  FLT_OR_DBL  *prm_MLbu;
} helper_arrays;

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

struct hc_hp_def_dat {
  int                       n;
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};


struct hc_int_def_dat {
  unsigned char             *mx;
  unsigned char             **mx_local;
  unsigned int              *sn;
  unsigned int              n;
  int                       *up;

  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};

typedef unsigned char (*eval_hc)(int   i,
                                int   j,
                                int   k,
                                int   l,
                                void  *data);

struct hc_mb_def_dat {
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  unsigned int              n;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};



typedef FLT_OR_DBL (*sc_ext_exp_cb)(int                    i,
                                   int                    j,
                                   int                    k,
                                   int                    l,
                                   struct sc_ext_exp_dat  *data);

typedef FLT_OR_DBL (*sc_ext_exp_red_up)(int                    i,
                                       int                    j,
                                       struct sc_ext_exp_dat  *data);

typedef FLT_OR_DBL (*sc_ext_exp_split)(int                   i,
                                      int                   j,
                                      int                   k,
                                      struct sc_ext_exp_dat *data);

struct sc_ext_exp_dat {
  FLT_OR_DBL                  **up;

  sc_ext_exp_cb               red_ext;
  sc_ext_exp_cb               red_stem;
  sc_ext_exp_red_up           red_up;
  sc_ext_exp_split            split;

  vrna_sc_exp_f user_cb;
  void                        *user_data;

  /* below attributes are for comparative structure prediction */
  int                         n_seq;
  unsigned int                **a2s;
  FLT_OR_DBL                  ***up_comparative;

  vrna_sc_exp_f *user_cb_comparative;
  void                        **user_data_comparative;
};


typedef FLT_OR_DBL (*sc_hp_exp_cb)(int                   i,
                                  int                   j,
                                  struct sc_hp_exp_dat  *data);
struct sc_hp_exp_dat {
  int                         n;
  unsigned int                n_seq;
  unsigned int                **a2s;
  int                         *idx;

  FLT_OR_DBL                  **up;
  FLT_OR_DBL                  ***up_comparative;
  FLT_OR_DBL                  *bp;
  FLT_OR_DBL                  **bp_comparative;
  FLT_OR_DBL                  **bp_local;
  FLT_OR_DBL                  ***bp_local_comparative;

  vrna_sc_exp_f user_cb;
  void                        *user_data;

  vrna_sc_exp_f *user_cb_comparative;
  void                        **user_data_comparative;

  sc_hp_exp_cb                pair;
  sc_hp_exp_cb                pair_ext;
};



typedef FLT_OR_DBL (*sc_int_exp_cb)(int                    i,
                                   int                    j,
                                   int                    k,
                                   int                    l,
                                   struct sc_int_exp_dat  *data);

struct sc_int_exp_dat {
  unsigned int                n;
  int                         n_seq;
  unsigned int                **a2s;

  int                         *idx;
  FLT_OR_DBL                  **up;
  FLT_OR_DBL                  ***up_comparative;
  FLT_OR_DBL                  *bp;
  FLT_OR_DBL                  **bp_comparative;
  FLT_OR_DBL                  **bp_local;
  FLT_OR_DBL                  ***bp_local_comparative;
  FLT_OR_DBL                  *stack;
  FLT_OR_DBL                  **stack_comparative;

  vrna_sc_exp_f user_cb;
  void                        *user_data;

  vrna_sc_exp_f *user_cb_comparative;
  void                        **user_data_comparative;

  sc_int_exp_cb               pair;
  sc_int_exp_cb               pair_ext;
};

typedef FLT_OR_DBL (*sc_mb_exp_pair_cb)(int                  i,
                                       int                  j,
                                       struct sc_mb_exp_dat *data);


typedef FLT_OR_DBL (*sc_mb_exp_red_cb)(int                   i,
                                      int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_mb_exp_dat  *data);


struct sc_mb_exp_dat {
  unsigned int                n;
  unsigned int                n_seq;
  unsigned int                **a2s;

  int                         *idx;

  FLT_OR_DBL                  **up;
  FLT_OR_DBL                  ***up_comparative;
  FLT_OR_DBL                  *bp;
  FLT_OR_DBL                  **bp_comparative;
  FLT_OR_DBL                  **bp_local;
  FLT_OR_DBL                  ***bp_local_comparative;

  sc_mb_exp_pair_cb           pair;
  sc_mb_exp_pair_cb           pair_ext;
  sc_mb_exp_red_cb            red_stem;
  sc_mb_exp_red_cb            red_ml;
  sc_mb_exp_red_cb            decomp_ml;

  vrna_sc_exp_f user_cb;
  void                        *user_data;

  vrna_sc_exp_f *user_cb_comparative;
  void                        **user_data_comparative;
};

typedef struct {
  struct hc_ext_def_dat     hc_dat_ext;
  vrna_hc_eval_f hc_eval_ext;

  struct hc_hp_def_dat      hc_dat_hp;
  vrna_hc_eval_f hc_eval_hp;

  struct hc_int_def_dat     hc_dat_int;
  eval_hc                   hc_eval_int;

  struct hc_mb_def_dat      hc_dat_mb;
  vrna_hc_eval_f hc_eval_mb;

  struct sc_ext_exp_dat     sc_wrapper_ext;
  struct sc_hp_exp_dat      sc_wrapper_hp;
  struct sc_int_exp_dat     sc_wrapper_int;
  struct sc_mb_exp_dat      sc_wrapper_mb;
} constraints_helper;

