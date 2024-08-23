#include "constant.h" 
#include "intl11.h"
#include "intl11dH.h"
#include "intl21.h"
#include "intl21dH.h"
#include "intl22.h"
#include "intl22dH.h"


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif


struct vrna_md_s
{
  double temperature;    /**<  @brief  The temperature used to scale the thermodynamic parameters */
  double betaScale;      /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
  int pf_smooth;         /**<  @brief  A flat specifying whether energies in Boltzmann factors need to be smoothed */
  int dangles;           /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3) */
  int special_hp;        /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
  int noLP;              /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
  int noGU;              /**<  @brief  Do not allow GU pairs */
  int noGUclosure;       /**<  @brief  Do not allow loops to be closed by GU pair */
  int logML;             /**<  @brief  Use logarithmic scaling for multiloops */
  int circ;              /**<  @brief  Assume RNA to be circular instead of linear */
  int gquad;             /**<  @brief  Include G-quadruplexes in structure prediction */
  int uniq_ML;           /**<  @brief  Flag to ensure unique multi-branch loop decomposition during folding */
  int energy_set;        /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
  int backtrack;         /**<  @brief  Specifies whether or not secondary structures should be backtraced */
  char backtrack_type;   /**<  @brief  Specifies in which matrix to backtrack */
  int compute_bpp;       /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
  char nonstandards[64]; /**<  @brief  contains allowed non standard bases */
  int max_bp_span;       /**<  @brief  maximum allowed base pair span */

  int min_loop_size;                    /**<  @brief  Minimum size of hairpin loops
                                         *
                                         *    The default value for this field is #TURN, however, it may
                                         *    be 0 in cofolding context.
                                         */
  int window_size;                      /**<  @brief  Size of the sliding window for locally optimal structure prediction */
  int oldAliEn;                         /**<  @brief  Use old alifold energy model */
  int ribo;                             /**<  @brief  Use ribosum scoring table in alifold energy model */
  double cv_fact;                       /**<  @brief  Co-variance scaling factor for consensus structure prediction */
  double nc_fact;                       /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
  double sfact;                         /**<  @brief  Scaling factor for partition function scaling */
  int rtype[8];                         /**<  @brief  Reverse base pair type array */
  short alias[MAXALPHA + 1];            /**<  @brief  alias of an integer nucleotide representation */
  int pair[MAXALPHA + 1][MAXALPHA + 1]; /**<  @brief  Integer representation of a base pair */
  float pair_dist[7][7];                /**<  @brief  Base pair dissimilarity, a.k.a. distance matrix */
  double salt;                          /**<  @brief  Salt (monovalent) concentration (M) in buffer */
  int saltMLLower;                      /**<  @brief  Lower bound of multiloop size to use in loop salt correction linear fitting */
  int saltMLUpper;                      /**<  @brief  Upper bound of multiloop size to use in loop salt correction linear fitting */
  int saltDPXInit;                      /**<  @brief  User-provided salt correction for duplex initialization (in dcal/mol).
                                         *    If set to 99999 the default salt correction is used.
                                         *    If set to 0 there is no salt correction for duplex initialization.
                                         */
  float saltDPXInitFact;                /**<  @brief  */
  float helical_rise;                   /**<  @brief  */
  float backbone_length;                /**<  @brief  */
};

typedef struct vrna_md_s vrna_md_t;

struct vrna_param_s
{
  int id;
  int stack[NBPAIRS + 1][NBPAIRS + 1];
  int hairpin[31];
  int bulge[MAXLOOP + 1];
  int internal_loop[MAXLOOP + 1];
  int mismatchExt[NBPAIRS + 1][5][5];
  int mismatchI[NBPAIRS + 1][5][5];
  int mismatch1nI[NBPAIRS + 1][5][5];
  int mismatch23I[NBPAIRS + 1][5][5];
  int mismatchH[NBPAIRS + 1][5][5];
  int mismatchM[NBPAIRS + 1][5][5];
  int dangle5[NBPAIRS + 1][5];
  int dangle3[NBPAIRS + 1][5];
  int int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
  int int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
  int int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
  int ninio[5];
  double lxc;
  int MLbase;
  int MLintern[NBPAIRS + 1];
  int MLclosing;
  int TerminalAU;
  int DuplexInit;
  int Tetraloop_E[200];
  char Tetraloops[1401];
  int Triloop_E[40];
  char Triloops[241];
  int Hexaloop_E[40];
  char Hexaloops[1801];
  int TripleC;
  int MultipleCA;
  int MultipleCB;
  int gquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
  int gquadLayerMismatch;
  int gquadLayerMismatchMax;

  double temperature; /**<  @brief  Temperature used for loop contribution scaling */

  vrna_md_t model_details; /**<  @brief  Model details to be used in the recursions */
  char param_file[256];    /**<  @brief  The filename the parameters were derived from, or empty string if they represent the default */
  int SaltStack;
  int SaltLoop[MAXLOOP + 2];
  double SaltLoopDbl[MAXLOOP + 2];
  int SaltMLbase;
  int SaltMLintern;
  int SaltMLclosing;
  int SaltDPXInit;
};

typedef struct vrna_param_s vrna_param_t;
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
    {0},
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
    VRNA_MODEL_DEFAULT_BACKBONE_LENGTH};

PRIVATE vrna_param_t *P = NULL;
PRIVATE int **c = NULL; /* energy array, given that i-j pair */
PRIVATE short *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;
PRIVATE int n1, n2; /* sequence lengths */

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

typedef struct
{
  int i;
  int j;
  int end;
  char *structure;
  double energy;
  double energy_backtrack;
  double opening_backtrack_x;
  double opening_backtrack_y;
  int offset;
  double dG1;
  double dG2;
  double ddG;
  int tb;
  int te;
  int qb;
  int qe;
} duplexT;

PRIVATE duplexT
duplexfold_cu(const char *s1,
              const char *s2,
              int clean_up);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_nucleotide_encode(char c,
                       vrna_md_t *md)
{
  /* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */
  int code = -1;

  c = toupper(c);

  if (md)
  {
    if (md->energy_set > 0)
    {
      code = (int)(c - 'A') + 1;
    }
    else
    {
      const char *pos;
      pos = strchr(Law_and_Order, c);
      if (pos == NULL)
        code = 0;
      else
        code = (int)(pos - Law_and_Order);

      if (code > 5)
        code = 0;

      if (code > 4)
        code--; /* make T and U equivalent */
    }
  }
  return code;
}

PUBLIC char
vrna_nucleotide_decode(int enc,
                       vrna_md_t *md)
{
  if (md)
  {
    if (md->energy_set > 0)
      return (char)enc + 'A' - 1;
    else
      return (char)Law_and_Order[enc];
  }
  else
  {
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

  md->alias[5] = 3; /* X <-> G */
  md->alias[6] = 2; /* K <-> C */
  md->alias[7] = 0; /* I <-> default base '@' */

  for (i = 0; i < NBASES; i++)
    for (j = 0; j < NBASES; j++)
      md->pair[i][j] = BP_pair[i][j];

  if (md->noGU)
    md->pair[3][4] = md->pair[4][3] = 0;

  if (md->nonstandards[0] != '\0')
  {
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
  switch (md->energy_set)
  {
  case 0:
    prepare_default_pairs(md);
    break;

  case 1:
    for (i = 1; i < MAXALPHA;)
    {
      md->alias[i++] = 3; /* A <-> G */
      md->alias[i++] = 2; /* B <-> C */
    }
    for (i = 1; i < MAXALPHA; i++)
    {
      md->pair[i][i + 1] = 2; /* AB <-> GC */
      i++;
      md->pair[i][i - 1] = 1; /* BA <-> CG */
    }

    break;

  case 2:
    for (i = 1; i < MAXALPHA;)
    {
      md->alias[i++] = 1; /* A <-> A*/
      md->alias[i++] = 4; /* B <-> U */
    }
    for (i = 1; i < MAXALPHA; i++)
    {
      md->pair[i][i + 1] = 5; /* AB <-> AU */
      i++;
      md->pair[i][i - 1] = 6; /* BA <-> UA */
    }

    break;

  case 3:
    for (i = 1; i < MAXALPHA - 2;)
    {
      md->alias[i++] = 3; /* A <-> G */
      md->alias[i++] = 2; /* B <-> C */
      md->alias[i++] = 1; /* C <-> A */
      md->alias[i++] = 4; /* D <-> U */
    }
    for (i = 1; i < MAXALPHA - 2; i++)
    {
      md->pair[i][i + 1] = 2; /* AB <-> GC */
      i++;
      md->pair[i][i - 1] = 1; /* BA <-> CG */
      i++;
      md->pair[i][i + 1] = 5; /* CD <-> AU */
      i++;
      md->pair[i][i - 1] = 6; /* DC <-> UA */
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
  md->rtype[0] = 0;
  md->rtype[7] = 7;

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
  if (md)
  {
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
    md->dangles = dangles;
    md->special_hp = tetra_loop;
    md->noLP = noLonelyPairs;
    md->noGU = noGU;
    md->noGUclosure = no_closingGU;
    md->logML = logML;
    md->gquad = gquad;
    md->circ = circ;
    md->uniq_ML = uniq_ML;
    md->compute_bpp = do_backtrack;
    md->backtrack = VRNA_MODEL_DEFAULT_BACKTRACK;
    md->backtrack_type = backtrack_type;
    md->energy_set = energy_set;
    md->max_bp_span = max_bp_span;
    md->min_loop_size = TURN;
    md->window_size = VRNA_MODEL_DEFAULT_WINDOW_SIZE;
    md->oldAliEn = oldAliEn;
    md->ribo = ribo;
    md->cv_fact = cv_fact;
    md->nc_fact = nc_fact;
    md->temperature = temperature;
    md->betaScale = VRNA_MODEL_DEFAULT_BETA_SCALE;
    md->pf_smooth = VRNA_MODEL_DEFAULT_PF_SMOOTH;
    md->sfact = 1.07;
    md->salt = defaults.salt;
    md->saltMLLower = defaults.saltMLLower;
    md->saltMLUpper = defaults.saltMLUpper;
    md->saltDPXInit = defaults.saltDPXInit;
    md->saltDPXInitFact = defaults.saltDPXInitFact;
    md->helical_rise = defaults.helical_rise;
    md->backbone_length = defaults.backbone_length;

    fill_pair_matrices(md);
  }
}


double expn(int n, double x)
{
  if (n < 0 || x < 0)
  {
    fprintf(stderr, "Invalid input: n and x must be non-negative\n");
    exit(EXIT_FAILURE);
  }

  if (x > MAXLOG)
  {
    return 0.0;
  }

  if (x == 0.0)
  {
    if (n < 2)
    {
      fprintf(stderr, "Singularity: E_n(x) is undefined for n < 2\n");
      exit(EXIT_FAILURE);
    }
    else
    {
      return 1.0 / (n - 1.0);
    }
  }

  if (n == 0)
  {
    return exp(-x) / x;
  }

  if (n > 5000)
  {
    double xk = x + n;
    double yk = 1.0 / (xk * xk);
    double t = n;
    double ans = yk * t * (6.0 * x * x - 8.0 * t * x + t * t);
    ans = yk * (ans + t * (t - 2.0 * x));
    ans = yk * (ans + t);
    ans = (ans + 1.0) * exp(-x) / xk;
    return ans;
  }

  if (x > 1.0)
  {
    double psi = -EUL - log(x);
    for (int i = 1; i < n; ++i)
    {
      psi += 1.0 / i;
    }

    double z = -x;
    double xk = 0.0;
    double yk = 1.0;
    double pk = 1.0 - n;
    double ans = n != 1 ? 1.0 / pk : 0.0;

    while (1)
    {
      xk += 1.0;
      yk *= z / xk;
      pk += 1.0;
      if (pk != 0.0)
      {
        ans += yk / pk;
      }
      double t = fabs(yk / ans);
      if (t <= MACHEP)
      {
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

  while (1)
  {
    ++k;
    double yk = (k & 1) ? 1.0 : x;
    double xk = (k & 1) ? n + (k - 1) / 2 : k / 2;

    double pk = pkm1 * yk + pkm2 * xk;
    double qk = qkm1 * yk + qkm2 * xk;

    if (qk != 0)
    {
      double r = pk / qk;
      double t = fabs((ans - r) / r);
      ans = r;
      if (t <= MACHEP)
      {
        break;
      }
    }
    else
    {
      if (fabs(pk) > BIG)
      {
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

double kn(int nn, double x)
{
  int n = abs(nn);

  if (x <= 0.0)
  {
    if (x < 0.0)
    {
      fprintf(stderr, "kn: Domain error\n");
      exit(EXIT_FAILURE);
    }
    else
    {
      fprintf(stderr, "kn: Singularity\n");
      exit(EXIT_FAILURE);
    }
  }

  if (x > 9.55)
  {
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

    while (1)
    {
      double z = pn - pk * pk;
      t *= z / (fn * z0);
      double nk1f = fabs(t);
      if (i >= n && nk1f > nkf)
      {
        break;
      }
      nkf = nk1f;
      s += t;
      fn += 1.0;
      pk += 2.0;
      ++i;
      if (fabs(t / s) <= MACHEP)
      {
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

  if (n > 0)
  {
    pn = -EUL;
    double k = 1.0;
    for (int i = 1; i < n; ++i)
    {
      pn += 1.0 / k;
      k += 1.0;
      fn *= k;
    }

    zmn = tox;

    if (n == 1)
    {
      ans = 1.0 / x;
    }
    else
    {
      double nk1f = fn / n;
      double kf = 1.0;
      double s = nk1f;
      double z = -z0;
      double zn = 1.0;
      for (int i = 1; i < n; ++i)
      {
        nk1f = nk1f / (n - i);
        kf *= i;
        zn *= z;
        double t = nk1f * zn / kf;
        s += t;
        if ((MAXNUM - fabs(t)) < fabs(s))
        {
          fprintf(stderr, "kn: Overflow\n");
          exit(EXIT_FAILURE);
        }
        if (tox > 1.0 && (MAXNUM / tox) < zmn)
        {
          fprintf(stderr, "kn: Overflow\n");
          exit(EXIT_FAILURE);
        }
        zmn *= tox;
      }

      s *= 0.5;
      double t = fabs(s);
      if (zmn > 1.0 && (MAXNUM / zmn) < t)
      {
        fprintf(stderr, "kn: Overflow\n");
        exit(EXIT_FAILURE);
      }
      if (t > 1.0 && (MAXNUM / t) < zmn)
      {
        fprintf(stderr, "kn: Overflow\n");
        exit(EXIT_FAILURE);
      }
      ans = s * zmn;
    }
  }

  double tlg = 2.0 * log(0.5 * x);
  double pk = -EUL;
  double t;
  if (n == 0)
  {
    pn = pk;
    t = 1.0;
  }
  else
  {
    pn += 1.0 / n;
    t = 1.0 / fn;
  }

  double s = (pk + pn - tlg) * t;
  double k = 1.0;
  while (1)
  {
    t *= z0 / (k * (k + n));
    pk += 1.0 / k;
    pn += 1.0 / (k + n);
    s += (pk + pn - tlg) * t;
    k += 1.0;
    if (fabs(t / s) <= MACHEP)
    {
      break;
    }
  }

  s = 0.5 * s / zmn;
  if (n % 2 != 0)
  {
    s = -s;
  }
  ans += s;
  return ans;
}

// extern int *__errno_location(void) __THROW __attribute_const__;
// #define errno (*__errno_location())

#ifndef WITH_DMALLOC
/* include the following two functions only if not including <dmalloc.h> */
PUBLIC void *
vrna_alloc(unsigned size)
{
  void *pointer;

  if ((pointer = (void *)calloc(1, (size_t)size)) == NULL)
  {
    // #ifdef EINVAL
    //     if (errno == EINVAL)
    //     {
    //       fprintf(stderr, "vrna_alloc: requested size: %d\n", size);
    //       printf("Memory allocation failure -> EINVAL");
    //     }

    //     if (errno == ENOMEM)
    // #endif
    printf("Memory allocation failure -> no memory");
  }

  return pointer;
}

PUBLIC void *
vrna_realloc(void *p,
             unsigned size)
{
  if (p == NULL)
    return vrna_alloc(size);

  p = (void *)realloc(p, size);
  if (p == NULL)
  {
    // #ifdef EINVAL
    //     if (errno == EINVAL)
    //     {
    //       fprintf(stderr, "vrna_realloc: requested size: %d\n", size);
    //       printf("vrna_realloc allocation failure -> EINVAL");
    //     }

    //     if (errno == ENOMEM)
    // #endif
    printf("vrna_realloc allocation failure -> no memory");
  }

  return p;
}

#endif



PRIVATE INLINE double
epsilonr(double T)
{
  return 5321 / T + 233.76 - 0.9297 * T + 1.417 * T * T / 1000 - 0.8292 * T * T * T / 1000000;
}

PRIVATE INLINE double
bjerrum_length(double T)
{
  return 167100.052 / (T * epsilonr(T));
};

PRIVATE INLINE double
ionic_strength(double rho)
{
  return rho;
};

PRIVATE INLINE double
kappa(double rho, double T)
{
  return sqrt(bjerrum_length(T) * ionic_strength(rho)) / 8.1284;
};



PRIVATE double
approx_hyper(double y)
{
  double a, b, c;

  a = 1 / (pow(y, 6.) / pow(2 * PI, 6.) + 1);
  b = pow(y, 4.) / (36 * pow(PI, 4.)) - pow(y, 3.) / (24 * PI * PI) + y * y / (2 * PI * PI) - y / 2;
  c = log(2 * PI / y) - 1.96351;

  return a * b + (1 - a) * c;
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

  bjerrum_length_inv = 1 / bjerrum_length(T);
  return MIN2(1 / backbonelen, bjerrum_length_inv);
};


PRIVATE double
loop_salt_aux(double kmlss, int L, double T, double backbonelen)
{
  double a, b;

  a = (GASCONST / 1000.) * T * bjerrum_length(T) * L * backbonelen * tau_ss(T, backbonelen) * tau_ss(T, backbonelen);
  b = log(kmlss) - log(PI / 2) + Eular_const + approx_hyper(kmlss) + 1 / kmlss * (1 - exp(-kmlss) + kmlss * expn(1, kmlss));

  return a * b * 100;
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
  return 2 * (GASCONST / 1000.) * T * bjerrum_length(T) * Helical_Rise * tau_ds(T, Helical_Rise) * tau_ds(T, Helical_Rise);
};


PUBLIC int
vrna_salt_stack(double rho, double T, double hrise)
{
  double correction, kn_ref, kappa_ref;

  kappa_ref = kappa(VRNA_MODEL_DEFAULT_SALT, T);
  kn_ref = kn(0, Rods_dist * kappa_ref);
  correction = 100 * pairing_salt_const(T, hrise) * (kn(0, Rods_dist * kappa(rho, T)) - kn_ref);
  return roundint(correction);
}

PUBLIC void
vrna_salt_ml(double saltLoop[], int lower, int upper, int *m, int *b)
{
  int sumx, sumxx;
  double y, sumy, sumyy, sumxy, denom, dm, db;

  sumx = sumxx = 0;
  sumy = sumyy = sumxy = 0.;

  for (int i = lower; i <= upper; i++)
  {
    sumx += i;
    sumxx += i * i;

    y = saltLoop[i];

    sumxy += i * y;
    sumy += y;
    sumyy += y * y;
  }

  denom = (double)((upper - lower + 1) * sumxx - sumx * sumx);
  dm = ((upper - lower + 1) * sumxy - sumx * sumy) / denom;
  db = (sumy * sumxx - sumx * sumxy) / denom;

  *m = roundint(dm);
  *b = roundint(db);
}

PUBLIC vrna_md_t *
vrna_md_copy(vrna_md_t *md_to,
             const vrna_md_t *md_from)
{
  int i;
  vrna_md_t *md;

  md = NULL;

  /* only process if md_from is non-NULL */
  if (md_from)
  {
    if (!md_to)
      /* create container to be filled */
      md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    else
      /* or directly write to target */
      md = md_to;

    /* check if not the same object */
    if (md_to != md_from)
    {
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
  double a, x, penalty;
  vrna_md_t md;

  if (md_p == NULL)
  {
    vrna_md_set_default(&md);
    md_p = &md;
  }

  if (md_p->saltDPXInit != 99999)
  {
    return md_p->saltDPXInit;
  }
  else
  {
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



PRIVATE vrna_param_t *
get_scaled_params(vrna_md_t *md)
{
  unsigned int i, j, k, l;
  int id = 0;
  double tempf;
  vrna_param_t *params;
  double salt = md->salt;
  double saltStandard = VRNA_MODEL_DEFAULT_SALT;

  params = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));

  params->model_details = *md; /* copy over the model details */
  params->temperature = md->temperature;
  // tempf                 = ((md->temperature + K0) / Tmeasure);
  tempf = 1.;
  params->ninio[2] = RESCALE_dG(ninio37, niniodH, tempf);
  params->lxc = lxc37 * tempf;
  params->TripleC = RESCALE_dG(TripleC37, TripleCdH, tempf);
  params->MultipleCA = RESCALE_dG(MultipleCA37, MultipleCAdH, tempf);
  params->MultipleCB = RESCALE_dG(MultipleCB37, MultipleCBdH, tempf);
  params->TerminalAU = RESCALE_dG(TerminalAU37, TerminalAUdH, tempf);
  params->DuplexInit = RESCALE_dG(DuplexInit37, DuplexInitdH, tempf);
  params->MLbase = RESCALE_dG(ML_BASE37, ML_BASEdH, tempf);
  params->MLclosing = RESCALE_dG(ML_closing37, ML_closingdH, tempf);
  params->gquadLayerMismatch = RESCALE_dG(GQuadLayerMismatch37, GQuadLayerMismatchH, tempf);
  params->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++)
    {
      double GQuadAlpha_T = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, tempf);
      double GQuadBeta_T = RESCALE_dG(GQuadBeta37, GQuadBetadH, tempf);
      params->gquad[i][j] = (int)GQuadAlpha_T * (i - 1) + (int)(((double)GQuadBeta_T) * log(j - 2));
    }

  for (i = 0; i < 31; i++)
    params->hairpin[i] = RESCALE_dG(hairpin37[i], hairpindH[i], tempf);

  for (i = 0; i <= MIN2(30, MAXLOOP); i++)
  {
    params->bulge[i] = RESCALE_dG(bulge37[i], bulgedH[i], tempf);
    params->internal_loop[i] = RESCALE_dG(internal_loop37[i], internal_loopdH[i], tempf);
  }

  for (; i <= MAXLOOP; i++)
  {
    params->bulge[i] = params->bulge[30] +
                       (int)(params->lxc * log((double)(i) / 30.));
    params->internal_loop[i] = params->internal_loop[30] +
                               (int)(params->lxc * log((double)(i) / 30.));
    params->SaltLoopDbl[i] = (salt == saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i] = (int)(params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i] < 0));
  }

  for (i = 0; i <= MIN2(31, MAXLOOP + 1); i++)
  {
    params->SaltLoopDbl[i] = (salt == saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i] = (int)(params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i] < 0));
  }

  for (; i <= MAXLOOP; i++)
  {
    params->SaltLoopDbl[i] = (salt == saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i] = (int)(params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i] < 0));
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
      for (k = 0; k < 5; k++)
      {
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
        if (md->dangles)
        {
          mm = RESCALE_dG(mismatchM37[i][j][k],
                          mismatchMdH[i][j][k],
                          tempf);
          params->mismatchM[i][j][k] = (mm > 0) ? 0 : mm;
          mm = RESCALE_dG(mismatchExt37[i][j][k],
                          mismatchExtdH[i][j][k],
                          tempf);
          params->mismatchExt[i][j][k] = (mm > 0) ? 0 : mm;
        }
        else
        {
          params->mismatchM[i][j][k] = params->mismatchExt[i][j][k] = 0;
        }
      }

  /* dangles */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
    {
      int dd;
      dd = RESCALE_dG(dangle5_37[i][j],
                      dangle5_dH[i][j],
                      tempf);
      params->dangle5[i][j] = (dd > 0) ? 0 : dd; /* must be <= 0 */
      dd = RESCALE_dG(dangle3_37[i][j],
                      dangle3_dH[i][j],
                      tempf);
      params->dangle3[i][j] = (dd > 0) ? 0 : dd; /* must be <= 0 */
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
        for (l = 0; l < 5; l++)
        {
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
        for (l = 0; l < 5; l++)
        {
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
  params->SaltStack = (salt == saltStandard) ? 0 : vrna_salt_stack(salt, saltT, md->helical_rise);
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
  if (md)
  {
    return get_scaled_params(md);
  }
  else
  {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_params(&md);
  }
}




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
  if (energy_set > 0)
  {
    code = (int)(c - 'A') + 1;
  }
  else
  {
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
      code--; /* make T and U equivalent */
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
  alias[5] = 3; /* X <-> G */
  alias[6] = 2; /* K <-> C */
  alias[7] = 0; /* I <-> default base '@' */
  for (i = 0; i < NBASES; i++)
    for (j = 0; j < NBASES; j++)
      pair[i][j] = BP_pair[i][j];
  // 如果noGU为真，禁止G-U配对。
  if (noGU)
    pair[3][4] = pair[4][3] = 0;
  // 如果nonstandards不为空，允许非标准碱基对
  if (nonstandards != NULL)
  {
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
vrna_E_ext_stem(unsigned int type,
                int n5d,
                int n3d,
                vrna_param_t *p)
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
E_IntLoop(int n1,
          int n2,
          int type,
          int type_2,
          int si1,
          int sj1,
          int sp1,
          int sq1,
          vrna_param_t *P)
{

  int nl, ns, u, energy, salt_stack_correction, salt_loop_correction, backbones;

  salt_stack_correction = P->SaltStack;
  salt_loop_correction = 0;

  if (n1 > n2)
  {
    nl = n1;
    ns = n2;
  }
  else
  {
    nl = n2;
    ns = n1;
  }

  if (nl == 0)
  {
    return P->stack[type][type_2] + salt_stack_correction; /* stack */
  }

  backbones = nl + ns + 2;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT)
  {
    /* salt correction for loop */
    if (backbones <= MAXLOOP + 1)
      salt_loop_correction = P->SaltLoop[backbones];
    else
      salt_loop_correction = vrna_salt_loop_int(backbones, P->model_details.salt, P->temperature + K0, P->model_details.backbone_length);
  }

  if (ns == 0)
  {
    /* bulge */
    energy = (nl <= MAXLOOP) ? P->bulge[nl] : (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
    if (nl == 1)
    {
      energy += P->stack[type][type_2];
    }
    else
    {
      if (type > 2)
        energy += P->TerminalAU;

      if (type_2 > 2)
        energy += P->TerminalAU;
    }
  }
  else
  {
    /* interior loop */
    if (ns == 1)
    {
      if (nl == 1) /* 1x1 loop */
        return P->int11[type][type_2][si1][sj1] + salt_loop_correction;

      if (nl == 2)
      {
        /* 2x1 loop */
        if (n1 == 1)
          energy = P->int21[type][type_2][si1][sq1][sj1];
        else
          energy = P->int21[type_2][type][sq1][si1][sp1];

        return energy + salt_loop_correction;
      }
      else
      {
        /* 1xn loop */
        energy =
            (nl + 1 <=
             MAXLOOP)
                ? (P->internal_loop[nl + 1])
                : (P->internal_loop[30] +
                   (int)(P->lxc * log((nl + 1) / 30.)));
        energy += MIN2(MAX_NINIO, (nl - ns) * P->ninio[2]);
        energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
        return energy + salt_loop_correction;
      }
    }
    else if (ns == 2)
    {
      if (nl == 2)
      {
        /* 2x2 loop */
        return P->int22[type][type_2][si1][sp1][sq1][sj1] + salt_loop_correction;
      }
      else if (nl == 3)
      {
        /* 2x3 loop */
        energy = P->internal_loop[5] + P->ninio[2];
        energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
        return energy + salt_loop_correction;
      }
    }

    {
      /* generic interior loop (no else here!)*/
      u = nl + ns;
      energy =
          (u <=
           MAXLOOP)
              ? (P->internal_loop[u])
              : (P->internal_loop[30] + (int)(P->lxc * log((u) / 30.)));

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
encode_sequence(const char *sequence,
                short how)
{
  /*
  i 和 l：用于循环和序列长度。
  S：分配内存以存储编码后的序列，长度为 l + 2。这里额外添加两个元素的目的是为了处理循环结构中的边界情况。
  */
  unsigned int i, l = (unsigned int)strlen(sequence);
  short *S = (short *)vrna_alloc(sizeof(short) * (l + 2));

  switch (how)
  {
  /* standard encoding as always used for S */
  // 遍历输入序列 sequence，调用 encode_char 函数将每个字符编码为数字，并将结果存储在数组 S 中。
  // S[l + 1] = S[1]：将序列首尾相连，用于处理循环结构。
  // S[0] = (short)l：存储序列的长度在 S[0] 中。
  case 0:
    for (i = 1; i <= l; i++) /* make numerical encoding of sequence */
      S[i] = (short)encode_char(sequence[i - 1]);
    S[l + 1] = S[1];
    S[0] = (short)l;
    break;
  /* encoding for mismatches of nostandard bases (normally used for S1) */
  // 同样遍历输入序列 sequence，但是通过 alias 数组将每个字符的编码映射为新的编码，并将结果存储在数组 S 中。
  // S[l + 1] = S[1]：同样将序列首尾相连。
  // S[0] = S[l]：存储序列末尾元素的值作为 S[0]，用于处理特定情况。
  case 1:
    for (i = 1; i <= l; i++)
      S[i] = alias[(short)encode_char(sequence[i - 1])];
    S[l + 1] = S[1];
    S[0] = S[l];
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
duplexfold_cu(const char *s1,
              const char *s2,
              int clean_up)
{
  //   i 和 j: 循环变量。
  // Emin: 最小能量值，初始化为无穷大（INF）。
  // i_min 和 j_min: 记录能量最小时的位置。
  // struc: 存储回溯得到的结构字符串。
  // mfe: 用于存储最小自由能的双链结构。
  // md: 结构体，存储模型参数。
  int i, j, Emin = INF, i_min = 0, j_min = 0;
  char *struc;
  duplexT mfe;
  vrna_md_t md;
  // 获取序列长度
  n1 = (int)strlen(s1);
  n2 = (int)strlen(s2);
  // 初始化模型参数 md。
  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature) > 1e-6))
  {
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
  S1 = encode_sequence(s1, 0);
  S2 = encode_sequence(s2, 0);
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
  for (i = 1; i <= n1; i++)
  {
    for (j = n2; j > 0; j--)
    {
      // P->DuplexInit: 初始化能量。
      // vrna_E_ext_stem: 计算外部茎的能量。
      // E_IntLoop: 计算内部环的能量。
      // MIN2: 返回两个值中的最小值。
      /* c: energy array, given that i-j pair */
      int type, type2, E, k, l;
      // type: 序列 s1[i] 和 s2[j] 的配对类型 例如，A-U，G-C。
      // pair no erro !!!!!!!!!!!!!!!!!!!!
      type = pair[S1[i]][S2[j]];
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
      if (!type)
      {
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
      for (k = i - 1; k > 0 && k > i - MAXLOOP - 2; k--)
      {
        for (l = j + 1; l <= n2; l++)
        {
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
      if (E < Emin)
      {
        Emin = E;
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
  mfe.i = i_min;
  mfe.j = j_min;
  mfe.energy = (float)Emin / 100.;
  mfe.structure = struc;
  // 清理内存
  if (clean_up)
  {
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

// 用于计算给定的两个RNA序列 s1 和 s2 的次优（suboptimal）二级结构。这些次优结构的自由能在最低自由能结构（MFE）的基础上增加一个给定的能量范围 delta 以内。函数返回一个包含这些次优结构的数组 w: 一个窗口参数，用于限制次优结构的数量
PUBLIC duplexT *

duplex_subopt(const char *s1,
              const char *s2,
              int delta,
              int w)
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

  int i, j, n1, n2, thresh, E, n_subopt = 0, n_max;
  char *struc;
  duplexT mfe;
  duplexT *subopt;

  /*
  初始化 n_max 为 16。
  使用 vrna_alloc 分配初始大小为 16 的 duplexT 数组 subopt。
  计算并获取最低自由能结构 mfe。
  释放 mfe.structure 所占用的内存。
  */

  n_max = 16;
  subopt = (duplexT *)vrna_alloc(n_max * sizeof(duplexT));
  mfe = duplexfold_cu(s1, s2, 0);
  free(mfe.structure);

  // 根据 mfe.energy 和 delta 计算能量阈值 threshold。
  // thresh = (int)mfe.energy * 100 + 0.1 + delta * 100;
  thresh = -500;
  n1 = strlen(s1);
  n2 = strlen(s2);
  // print_pair();
  // 二层遍历查找次优结构
  for (i = n1; i > 0; i--)
  {
    for (j = 1; j <= n2; j++)
    {
      int type, ii, jj, Ed;
      type = pair[S2[j]][S1[i]];
      if (!type)
        continue;
      // 计算能量
      E = Ed = c[i][j];
      // from /data/ntc/Repository/ViennaRNA-2.6.4/src/loops/external.c
      Ed += vrna_E_ext_stem(type, (j > 1) ? SS2[j - 1] : -1, (i < n1) ? SS1[i + 1] : -1, P);
      // 超过阈值则跳过
      if (Ed > thresh)
        continue;

      /* too keep output small, remove hits that are dominated by a
       * better one close (w) by. For simplicity we do test without
       * adding dangles, which is slightly inaccurate.
       */
      for (ii = MAX2(i - w, 1); (ii <= MIN2(i + w, n1)) && type; ii++)
      {
        for (jj = MAX2(j - w, 1); jj <= MIN2(j + w, n2); jj++)
          if (c[ii][jj] < E)
          {
            type = 0;
            break;
          }
      }
      // printf("%d,", thresh);
      if (!type)
        continue;

      struc = backtrack(i, j);

      // 如果 subopt 数组已满，使用 vrna_realloc 扩展数组大小
      if (n_subopt + 1 >= n_max)
      {
        n_max *= 2;
        subopt = (duplexT *)vrna_realloc(subopt, n_max * sizeof(duplexT));
      }

      subopt[n_subopt].i = MIN2(i + 1, n1);
      subopt[n_subopt].j = MAX2(j - 1, 1);
      subopt[n_subopt].energy = Ed * 0.01;
      subopt[n_subopt++].structure = struc;
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

  subopt[n_subopt].i = 0;
  subopt[n_subopt].j = 0;
  subopt[n_subopt].structure = NULL;
  return subopt;
}

PRIVATE char *
backtrack(int i,
          int j)
{
  /* backtrack structure going backwards from i, and forwards from j
   * return structure in bracket notation with & as separator */
  int k, l, type, type2, E, traced, i0, j0;
  char *st1, *st2, *struc;

  st1 = (char *)vrna_alloc(sizeof(char) * (n1 + 1));
  st2 = (char *)vrna_alloc(sizeof(char) * (n2 + 1));

  i0 = MIN2(i + 1, n1);
  j0 = MAX2(j - 1, 1);

  while (i > 0 && j <= n2)
  {
    E = c[i][j];
    traced = 0;
    st1[i - 1] = '(';
    st2[j - 1] = ')';
    type = pair[S1[i]][S2[j]];
    if (!type)
      printf("backtrack failed in fold duplex");

    for (k = i - 1; k > 0 && k > i - MAXLOOP - 2; k--)
    {
      for (l = j + 1; l <= n2; l++)
      {
        int LE;
        if (i - k + l - j - 2 > MAXLOOP)
          break;

        type2 = pair[S1[k]][S2[l]];
        if (!type2)
          continue;

        LE = E_IntLoop(i - k - 1, l - j - 1, type2, rtype[type],
                       SS1[k + 1], SS2[l - 1], SS1[i - 1], SS2[j + 1], P);
        if (E == c[k][l] + LE)
        {
          traced = 1;
          i = k;
          j = l;
          break;
        }
      }
      if (traced)
        break;
    }
    if (!traced)
    {
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


char *formatDuplexOutput(const duplexT *dup)
{
  int l1;
  l1 = strchr(dup->structure, '&') - dup->structure;
  char buffer[BUFFER_SIZE];
  buffer[0] = '\0'; // 清空缓冲区

  char temp[50];

  // 拼接各个部分
  sprintf(temp, "%d,", dup->i + 1 - l1);
  strcat(buffer, temp);

  sprintf(temp, "%d|", dup->i);
  strcat(buffer, temp);

  sprintf(temp, "%d,", dup->j);
  strcat(buffer, temp);

  sprintf(temp, "%d|", dup->j + (int)strlen(dup->structure) - l1 - 2);
  strcat(buffer, temp);

  sprintf(temp, "%5.2f|", dup->energy);
  strcat(buffer, temp);

  strcat(buffer, dup->structure);
  strcat(buffer, "\n");

  // 为结果字符串分配内存并复制缓冲区内容
  char *result = (char *)malloc(strlen(buffer) + 1);
  if (result != NULL)
  {
    strcpy(result, buffer);
  }

  return strdup(result);
}

char *process(char *s1, char *s2)
{
  int i, sym, istty, delta, noconv;
  char *result = (char *)malloc(1000000 * sizeof(char)); // 动态分配内存
  result[0] = '\0';
  delta = 10;
  duplexT mfe, *subopt;
  if (delta >= 0)
  {
    duplexT *sub;
    subopt = duplex_subopt(s1, s2, delta, 5);
    for (sub = subopt; sub->i > 0; sub++)
    {
      char *tmp_str = formatDuplexOutput(sub);
      strcat(result, tmp_str);
      free(tmp_str);
      free(sub->structure);
    }
    free(subopt);
  }
  // printf("%s",result);
  return result;
}

int main_dup(){
    char *s1 = "ACTGCCAAGTAGGAAAGTCCCATAAGGTCAT";
    char *s2 = "TGAACTTATGGGACTTTCCTACTTGGCAG";
    char *res;
    res = process(s1,s2);
    printf("%s", res);
    free(res);
    return 0;
}

// cofold start ----------------------------------------
// header include start--------------------------------


// structure definition start -------------------------------

typedef enum {
  VRNA_FC_TYPE_SINGLE,      /**< Type is suitable for single, and hybridizing sequences */
  VRNA_FC_TYPE_COMPARATIVE  /**< Type is suitable for sequence alignments (consensus structure prediction) */
} vrna_fc_type_e;

struct vrna_elem_prob_s {
  int   i;    /**<  @brief  Start position (usually 5' nucleotide that starts the element, e.g. base pair) */
  int   j;    /**<  @brief  End position (usually 3' nucleotide that ends the element, e.g. base pair) */
  float p;    /**<  @brief  Probability of the element */
  int   type; /**<  @brief  Type of the element */
};

typedef struct vrna_elem_prob_s vrna_ep_t;


struct vrna_cstr_s {
  char          *string;
  size_t        size;
  FILE          *output;
  unsigned char istty;
};

typedef struct vrna_cstr_s *vrna_cstr_t;
struct output_stream {
  vrna_cstr_t data;
  vrna_cstr_t err;
};

struct id_data {
  char      *name;
  int       auto_id;
  char      *prefix;
  char      *delimiter;
  int       digits;
  long int  number;
};
typedef struct id_data *dataset_id;


typedef void (*vrna_stream_output_f)(void        *auxdata,
                                            unsigned int i,
                                            void         *data);

struct vrna_ordered_stream_s {
  unsigned int                start;      /* first element index in queue, i.e. start of queue */
  unsigned int                end;        /* last element index in queue */
  unsigned int                size;       /* available memory size for 'data' and 'provided' */
  unsigned int                shift;      /* pointer offset for 'data' and 'provided' */

  vrna_stream_output_f output;    /* callback to execute if consecutive elements from head are available */
  void                        **data;     /* actual data passed to the callback */
  unsigned char               *provided;  /* for simplicity we use unsigned char instead of single bits per element */
  void                        *auxdata;   /* auxiliary data passed to the callback */
#if VRNA_WITH_PTHREADS
  pthread_mutex_t             mtx;        /* semaphore to provide concurrent access */
#endif
};

typedef struct vrna_ordered_stream_s *vrna_ostream_t;

#define MAX_ALPHABET (6)
#define MAX_PAIRS (NBPAIRS + 1 + 25)

struct vrna_sc_mod_param_s {
  unsigned int  available;

  char          *name;
  char          one_letter_code;
  char          unmodified;
  char          fallback;
  char          pairing_partners[7];
  unsigned int  pairing_partners_encoding[7];
  unsigned int  unmodified_encoding;
  unsigned int  fallback_encoding;

  size_t        num_ptypes;
  size_t        ptypes[MAX_ALPHABET][MAX_ALPHABET];

  int           stack_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int           stack_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int           dangle5_dG[MAX_PAIRS][MAX_ALPHABET];
  int           dangle5_dH[MAX_PAIRS][MAX_ALPHABET];
  int           dangle3_dG[MAX_PAIRS][MAX_ALPHABET];
  int           dangle3_dH[MAX_PAIRS][MAX_ALPHABET];

  int           mismatch_dG[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];
  int           mismatch_dH[MAX_PAIRS][MAX_ALPHABET][MAX_ALPHABET];

  int           terminal_dG[MAX_PAIRS];
  int           terminal_dH[MAX_PAIRS];
};

typedef struct vrna_sc_mod_param_s *vrna_sc_mod_param_t;
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

  dataset_id      id_control;

  char            *concentration_file;

  char            *constraint_file;
  int             constraint_batch;
  int             constraint_enforce;
  int             constraint_canonical;

  int             shape;
  char            *shape_file;
  char            *shape_method;
  char            *shape_conversion;

  vrna_sc_mod_param_t *mod_params;

  int             csv_output;
  int             csv_header;
  char            csv_output_delim;

  int             jobs;
  int             keep_order;
  unsigned int    next_record_number;
  vrna_ostream_t  output_queue;
};


struct record_data {
  unsigned int    number;
  char            *id;
  char            *sequence;
  char            *SEQ_ID;
  char            *input_filename;
  int             multiline_input;
  struct options  *options;
  int             tty;
};

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

typedef enum {
  VRNA_SEQ_UNKNOWN = 0,   /**< @brief Nucleotide sequence represents an Unkown type */
  VRNA_SEQ_RNA = 1,       /**< @brief Nucleotide sequence represents an RNA type */
  VRNA_SEQ = 2        /**< @brief Nucleotide sequence represents a DNA type */
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


typedef void (*vrna_recursion_status_f)(unsigned char status,
                                              void          *data);



struct vrna_structured_domains_s {
  char __placeholder; /* dummy placeholder to not leave this struct empty */
};
typedef struct vrna_structured_domains_s vrna_sd_t;





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
typedef struct  vrna_sc_s vrna_sc_t;


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
  vrna_sd_t     *domains_struc;             /**<  @brief  Additional structured domains */

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




struct vrna_dimer_pf_s {
  /* free energies for: */
  double  F0AB; /**< @brief Null model without DuplexInit */
  double  FAB;  /**< @brief all states with DuplexInit correction */
  double  FcAB; /**< @brief true hybrid states only */
  double  FA;   /**< @brief monomer A */
  double  FB;   /**< @brief monomer B */
};

typedef struct vrna_dimer_pf_s vrna_dimer_pf_t;

struct vrna_dimer_conc_s {
  double  Ac_start;   /**< @brief start concentration A */
  double  Bc_start;   /**< @brief start concentration B */
  double  ABc;        /**< @brief End concentration AB */
  double  AAc;
  double  BBc;
  double  Ac;
  double  Bc;
};



typedef struct vrna_dimer_conc_s vrna_dimer_conc_t;


#define VRNA_OPTION_PF (1 << 1)



PUBLIC int
vrna_gr_set_aux_exp_c(vrna_fold_compound_t      *fc,
                      vrna_grammar_rule_f_exp cb);

typedef struct {
  FLT_OR_DBL  *prm_l;
  FLT_OR_DBL  *prm_l1;
  FLT_OR_DBL  *prml;

  int         ud_max_size;
  FLT_OR_DBL  **pmlu;
  FLT_OR_DBL  *prm_MLbu;
} helper_arrays;

typedef unsigned char (*eval_hc)(int   i,
                                int   j,
                                int   k,
                                int   l,
                                void  *data);



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

struct hc_mb_def_dat {
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  unsigned int              n;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};

struct sc_ext_exp_dat;  // 前向声明


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
struct sc_hp_exp_dat;  // 前向声明


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

struct sc_int_exp_dat;  // 前向声明


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



struct sc_mb_exp_dat;

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



//function definition start --------------------------------



void
init_default_options(struct options *opt)
{
  opt->filename_full  = 0;
  opt->filename_delim = NULL;
  opt->pf             = 1;
  opt->doT            = 0; /* compute dimer free energies etc. */
  opt->noPS           = 1;
  opt->noconv         = 1;
  opt->centroid       = 0;  /* off by default due to historical reasons */
  opt->MEA            = 0;
  opt->MEAgamma       = 1.;
  opt->bppmThreshold  = 1e-5;
  opt->verbose        = 0;
  // opt->commands       = NULL;
  opt->id_control     = NULL;
  set_model_details(&(opt->md));

  opt->doC                = 0; /* toggle to compute concentrations */
  opt->concentration_file = NULL;

  opt->constraint_file      = NULL;
  opt->constraint_batch     = 0;
  opt->constraint_enforce   = 0;
  opt->constraint_canonical = 0;

  opt->shape            = 0;
  opt->shape_file       = NULL;
  opt->shape_method     = NULL;
  opt->shape_conversion = NULL;

  opt->mod_params         = NULL;

  opt->csv_output       = 0;    /* flag indicating whether we produce one-line outputs, a.k.a. CSV */
  opt->csv_header       = 1;    /* print header for one-line output */
  opt->csv_output_delim = ',';  /* delimiting character for one-line output */

  opt->jobs               = 1;
  opt->keep_order         = 1;
  opt->next_record_number = 0;
  opt->output_queue       = NULL;
}


PRIVATE INLINE void
nullify(vrna_fold_compound_t *fc)
{
  if (fc) {
    fc->length            = 0;
    fc->strands           = 0;
    fc->cutpoint          = -1;
    fc->strand_number     = NULL;
    fc->strand_order      = NULL;
    fc->strand_order_uniq = NULL;
    fc->strand_start      = NULL;
    fc->strand_end        = NULL;
    fc->nucleotides       = NULL;
    // fc->alignment         = NULL;

    fc->hc            = NULL;
    fc->matrices      = NULL;
    fc->exp_matrices  = NULL;
    fc->params        = NULL;
    fc->exp_params    = NULL;
    fc->iindx         = NULL;
    fc->jindx         = NULL;

    fc->stat_cb       = NULL;
    fc->auxdata       = NULL;
    fc->free_auxdata  = NULL;

    fc->domains_struc = NULL;
    fc->domains_up    = NULL;
    fc->aux_grammar   = NULL;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        fc->sequence            = NULL;
        fc->sequence_encoding   = NULL;
        fc->encoding5           = NULL;
        fc->encoding3           = NULL;
        fc->sequence_encoding2  = NULL;
        fc->ptype               = NULL;
        fc->ptype_pf_compat     = NULL;
        fc->sc                  = NULL;

        break;

    }

    fc->maxD1         = 0;
    fc->maxD2         = 0;
    fc->reference_pt1 = NULL;
    fc->reference_pt2 = NULL;
    fc->referenceBPs1 = NULL;
    fc->referenceBPs2 = NULL;
    fc->bpdist        = NULL;
    fc->mm1           = NULL;
    fc->mm2           = NULL;

    fc->window_size = -1;
    fc->ptype_local = NULL;
#ifdef VRNA_WITH_SVM
    fc->zscore_data = NULL;
#endif
  }
}



PRIVATE vrna_fold_compound_t *
init_fc_single(void)
{
  vrna_fold_compound_t  init = {
    .type = VRNA_FC_TYPE_SINGLE
  };
  vrna_fold_compound_t  *fc = vrna_alloc(sizeof(vrna_fold_compound_t));

  if (fc) {
    // printf("into nullify\n");
    memcpy(fc, &init, sizeof(vrna_fold_compound_t));
    nullify(fc);
  }

  return fc;
}

PRIVATE vrna_exp_param_t *
get_scaled_exp_params(vrna_md_t *md,
                      double    pfs)
{
  unsigned int      i, j, k, l;
  int               pf_smooth;
  double            kT, TT;
  double            GT;
  double            salt, saltStandard;
  vrna_exp_param_t  *pf;

  pf = (vrna_exp_param_t *)vrna_alloc(sizeof(vrna_exp_param_t));

  memset(pf->param_file, '\0', 256);
  // if (last_parameter_file() != NULL) /* not into here */
  //   strncpy(pf->param_file, last_parameter_file(), 255);

  pf->model_details = *md;
  pf->temperature   = md->temperature;
  pf->alpha         = md->betaScale;
  pf->kT            = kT = md->betaScale * (md->temperature + K0) * GASCONST; /* kT in cal/mol  */
  pf->pf_scale      = pfs;
  pf_smooth         = md->pf_smooth;
  TT                = (md->temperature + K0) / (Tmeasure);
  salt = md->salt;
  saltStandard = VRNA_MODEL_DEFAULT_SALT;

  pf->lxc                   = lxc37 * TT;
  pf->expDuplexInit         = RESCALE_BF(DuplexInit37, DuplexInitdH, TT, kT);
  pf->expTermAU             = RESCALE_BF(TerminalAU37, TerminalAUdH, TT, kT);
  pf->expMLbase             = RESCALE_BF(ML_BASE37, ML_BASEdH, TT, kT);
  pf->expMLclosing          = RESCALE_BF(ML_closing37, ML_closingdH, TT, kT);
  pf->expgquadLayerMismatch = RESCALE_BF(GQuadLayerMismatch37, GQuadLayerMismatchH, TT, kT);
  pf->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, TT);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, TT);
      GT = ((double)GQuadAlpha_T) * ((double)(i - 1)) + ((double)GQuadBeta_T) *
           log(((double)j) - 2.);
      pf->expgquad[i][j] = exp(-TRUNC_MAYBE(GT) * 10. / kT);
    }

  /* loop energies: hairpins, bulges, interior, mulit-loops */
  for (i = 0; i < 31; i++)
    pf->exphairpin[i] = RESCALE_BF(hairpin37[i], hairpindH[i], TT, kT);

  for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
    pf->expbulge[i]     = RESCALE_BF(bulge37[i], bulgedH[i], TT, kT);
    pf->expinternal[i]  = RESCALE_BF(internal_loop37[i], internal_loopdH[i], TT, kT);
  }

  if (salt==saltStandard) { /* into here */
    for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
      pf->SaltLoopDbl[i]   = 0.;
      pf->expSaltLoop[i]   = 1.;
    }

    for (i = 31; i <= MAXLOOP; i++) {
      pf->SaltLoopDbl[i] = 0.;
      pf->expSaltLoop[i] = 1.;
    }
  } else {
    for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
      pf->SaltLoopDbl[i]   = vrna_salt_loop(i, salt, saltT, md->backbone_length);
      int saltLoop = (int) (pf->SaltLoopDbl[i] + 0.5 - (pf->SaltLoopDbl[i]<0));
      pf->expSaltLoop[i]   = exp(-saltLoop * 10. / kT);
    }

    for (i = 31; i <= MAXLOOP; i++) {
      pf->SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md->backbone_length);
      int saltLoop = (int) (pf->SaltLoopDbl[i] + 0.5 - (pf->SaltLoopDbl[i]<0));
      pf->expSaltLoop[i] = exp(-saltLoop * 10. / kT);
    }
  }

  /* special case of size 2 interior loops (single mismatch) */
  if (james_rule) /* into */
    pf->expinternal[2] = exp(-80 * 10. / kT);

  GT = RESCALE_dG(bulge37[30],
                  bulgedH[30],
                  TT);
  for (i = 31; i <= MAXLOOP; i++)
    pf->expbulge[i] = exp(-TRUNC_MAYBE(GT + (pf->lxc * log(i / 30.))) * 10. / kT);

  GT = RESCALE_dG(internal_loop37[30],
                  internal_loopdH[30],
                  TT);
  for (i = 31; i <= MAXLOOP; i++)
    pf->expinternal[i] = exp(-TRUNC_MAYBE(GT + (pf->lxc * log(i / 30.))) * 10. / kT);

  GT = RESCALE_dG(ninio37, niniodH, TT);
  for (j = 0; j <= MAXLOOP; j++)
    pf->expninio[2][j] = exp(-MIN2(MAX_NINIO, j * TRUNC_MAYBE(GT)) * 10. / kT);

  for (i = 0; (i * 7) < strlen(Tetraloops); i++)
    pf->exptetra[i] = RESCALE_BF(Tetraloop37[i], TetraloopdH[i], TT, kT);

  for (i = 0; (i * 5) < strlen(Triloops); i++)
    pf->exptri[i] = RESCALE_BF(Triloop37[i], TriloopdH[i], TT, kT);

  for (i = 0; (i * 9) < strlen(Hexaloops); i++)
    pf->exphex[i] = RESCALE_BF(Hexaloop37[i], HexaloopdH[i], TT, kT);

  for (i = 0; i <= NBPAIRS; i++)
    pf->expMLintern[i] = RESCALE_BF(ML_intern37, ML_interndH, TT, kT);

  /* if dangles==0 just set their energy to 0,
   * don't let dangle energies become > 0 (at large temps),
   * but make sure go smoothly to 0                        */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= 4; j++) {
      if (md->dangles) {
        pf->expdangle5[i][j]  = RESCALE_BF_SMOOTH(dangle5_37[i][j], dangle5_dH[i][j], TT, kT);
        pf->expdangle3[i][j]  = RESCALE_BF_SMOOTH(dangle3_37[i][j], dangle3_dH[i][j], TT, kT);
      } else {
        pf->expdangle3[i][j] = pf->expdangle5[i][j] = 1;
      }
    }

  /* stacking energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      pf->expstack[i][j] = RESCALE_BF(stack37[i][j], stackdH[i][j], TT, kT);

  /* mismatch energies */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
      for (k = 0; k < 5; k++) {
        pf->expmismatchI[i][j][k] = RESCALE_BF(mismatchI37[i][j][k],
                                               mismatchIdH[i][j][k],
                                               TT,
                                               kT);
        pf->expmismatch1nI[i][j][k] = RESCALE_BF(mismatch1nI37[i][j][k],
                                                 mismatch1nIdH[i][j][k],
                                                 TT,
                                                 kT);
        pf->expmismatchH[i][j][k] = RESCALE_BF(mismatchH37[i][j][k],
                                               mismatchHdH[i][j][k],
                                               TT,
                                               kT);
        pf->expmismatch23I[i][j][k] = RESCALE_BF(mismatch23I37[i][j][k],
                                                 mismatch23IdH[i][j][k],
                                                 TT,
                                                 kT);

        if (md->dangles) {
          pf->expmismatchM[i][j][k] = RESCALE_BF_SMOOTH(mismatchM37[i][j][k],
                                                        mismatchMdH[i][j][k],
                                                        TT,
                                                        kT);
          pf->expmismatchExt[i][j][k] = RESCALE_BF_SMOOTH(mismatchExt37[i][j][k],
                                                          mismatchExtdH[i][j][k],
                                                          TT,
                                                          kT);
        } else {
          pf->expmismatchM[i][j][k] = pf->expmismatchExt[i][j][k] = 1.;
        }
      }

  /* interior lops of length 2 */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          pf->expint11[i][j][k][l] = RESCALE_BF(int11_37[i][j][k][l],
                                                int11_dH[i][j][k][l],
                                                TT,
                                                kT);
        }

  /* interior 2x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++) {
            pf->expint21[i][j][k][l][m] = RESCALE_BF(int21_37[i][j][k][l][m],
                                                     int21_dH[i][j][k][l][m],
                                                     TT,
                                                     kT);
          }
        }

  /* interior 2x2 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++) {
              pf->expint22[i][j][k][l][m][n] = RESCALE_BF(int22_37[i][j][k][l][m][n],
                                                          int22_dH[i][j][k][l][m][n],
                                                          TT,
                                                          kT);
            }
        }

  strncpy(pf->Tetraloops, Tetraloops, 281);
  strncpy(pf->Triloops, Triloops, 241);
  strncpy(pf->Hexaloops, Hexaloops, 361);

  /* Salt correction for stack and multiloop */
  pf->SaltMLbase = pf->SaltMLclosing = pf->SaltDPXInit = 0.;

  if (salt==saltStandard) { /* into */
    pf->expSaltStack = 1.;
  } else {
    pf->expSaltStack = exp(- vrna_salt_stack(salt, saltT, md->helical_rise) * 10. / kT);
    vrna_salt_ml(pf->SaltLoopDbl, md->saltMLLower, md->saltMLUpper, &pf->SaltMLbase, &pf->SaltMLclosing);

    if (md->saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT)
      pf->SaltDPXInit = md->saltDPXInit;
    else if (md->saltDPXInit) /* into */
      pf->SaltDPXInit = vrna_salt_duplex_init(md);

    pf->expMLclosing *= exp(- pf->SaltMLbase * 10. / kT);
    pf->expMLclosing *= exp(- pf->SaltMLclosing * 10. / kT);
    pf->expMLbase *= exp(- pf->SaltMLbase * 10. / kT);
    for (i = 0; i <= NBPAIRS; i++)
      pf->expMLintern[i] *= exp(- pf->SaltMLbase * 10. / kT);

    pf->expDuplexInit *= exp(- pf->SaltDPXInit*10. / kT);
  }

  return pf;
}




PUBLIC vrna_exp_param_t *
vrna_exp_params(vrna_md_t *md)
{
  if (md) {
    return get_scaled_exp_params(md, -1.);
  } else {
    vrna_md_t md;
    vrna_md_set_default(&md);
    return get_scaled_exp_params(&md, -1.);
  }
}

PUBLIC void
vrna_params_prepare(vrna_fold_compound_t  *fc,
                    unsigned int          options)
{
  if (fc) {
    vrna_md_t *md_p;

    /*
     *  every vrna_fold_compound_t must have a vrna_paramt_t structure attached
     *  to it that holds the current model details. So we just use this here as
     *  the reference model
     */
    md_p = &(fc->params->model_details);

    if (options & VRNA_OPTION_PF) {
      /* remove previous parameters if present and they differ from reference model */
      if (fc->exp_params) {
        if (memcmp(md_p, &(fc->exp_params->model_details), sizeof(vrna_md_t)) != 0) {
          free(fc->exp_params);
          fc->exp_params = NULL;
        }
      }

      if (!fc->exp_params)
        // printf("into exp\n");
        fc->exp_params = vrna_exp_params(md_p);
    }
  }
}


PRIVATE void
add_params(vrna_fold_compound_t *fc,
           vrna_md_t            *md_p,
           unsigned int         options)
{
  /*
   * ALWAYS provide regular energy parameters
   * remove previous parameters if present and they differ from current model
   */
  if (fc->params) {
    if (memcmp(md_p, &(fc->params->model_details), sizeof(vrna_md_t)) != 0) {
      free(fc->params);
      fc->params = NULL;
    }
  }

  if (!fc->params)
    fc->params = vrna_params(md_p);

  vrna_params_prepare(fc, options);
}


PRIVATE void
sanitize_bp_span(vrna_fold_compound_t *fc,
                 unsigned int         options)
{
  vrna_md_t *md;

  md = &(fc->params->model_details);

  /* make sure that min_loop_size, max_bp_span, and window_size are sane */
  if (options & VRNA_OPTION_WINDOW) {
    if (md->window_size <= 0)
      md->window_size = (int)fc->length;
    else if (md->window_size > (int)fc->length)
      md->window_size = (int)fc->length;

    fc->window_size = md->window_size;
  } else {
    /* non-local fold mode */
    md->window_size = (int)fc->length;
  }

  if ((md->max_bp_span <= 0) || (md->max_bp_span > md->window_size))
    md->max_bp_span = md->window_size;
}



PUBLIC int *
vrna_idx_col_wise(unsigned int length)
{
  unsigned int  i;
  int           *idx = (int *)vrna_alloc(sizeof(int) * (length + 1));

  for (i = 1; i <= length; i++)
    idx[i] = (i * (i - 1)) / 2;
  return idx;
}


PUBLIC char *
vrna_ptypes(const short *S,
            vrna_md_t   *md)
{
  char  *ptype;
  int   n, i, j, k, l, *idx;
  int   min_loop_size = md->min_loop_size;

  n = S[0];

  ptype = (char *)vrna_alloc(sizeof(char) * ((n * (n + 1)) / 2 + 2));
  idx   = vrna_idx_col_wise(n);

  for (k = 1; k < n - min_loop_size; k++)
    for (l = 1; l <= 2; l++) {
      int type, ntype = 0, otype = 0;
      i = k;
      j = i + min_loop_size + l;
      if (j > n)
        continue;

      type = md->pair[S[i]][S[j]];
      while ((i >= 1) && (j <= n)) {
        if ((i > 1) && (j < n))
          ntype = md->pair[S[i - 1]][S[j + 1]];

        if (md->noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */

        ptype[idx[j] + i] = (char)type;
        otype             = type;
        type              = ntype;
        i--;
        j++;
      }
    }
  free(idx);
  return ptype;
}


PUBLIC int *
vrna_idx_row_wise(unsigned int length)
{ 
  // printf("(%d)\n", length);
  int i;
  int *idx = (int *)vrna_alloc(sizeof(int) * (length + 1));

  for (i = 1; i <= length; i++)
    idx[i] = (((length + 1 - i) * (length - i)) / 2) + length + 1;
    // printf("%d-", idx[i]);
  return idx;
}


PRIVATE char *
wrap_get_ptypes(const short *S,
                vrna_md_t   *md)
{
  char  *ptype;
  int   n, i, j, k, l, *idx;

  n     = S[0];
  ptype = (char *)vrna_alloc(sizeof(char) * ((n * (n + 1)) / 2 + 2));
  idx   = vrna_idx_row_wise(n);
  int   min_loop_size = md->min_loop_size;

  for (k = 1; k < n - min_loop_size; k++)
    for (l = 1; l <= 2; l++) {
      int type, ntype = 0, otype = 0;
      i = k;
      j = i + min_loop_size + l;
      if (j > n)
        continue;

      type = md->pair[S[i]][S[j]];
      while ((i >= 1) && (j <= n)) {
        if ((i > 1) && (j < n))
          ntype = md->pair[S[i - 1]][S[j + 1]];

        if (md->noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */

        ptype[idx[i] - j] = (char)type;
        otype             = type;
        type              = ntype;
        i--;
        j++;
      }
    }
  free(idx);
  return ptype;
}

PUBLIC char *
get_ptypes(const short  *S,
           vrna_md_t    *md,
           unsigned int idx_type)
{
  if (S) {
    if ((unsigned int)S[0] > 1000000) {
      printf("get_ptypes@alphabet.c: sequence length of %d exceeds addressable range",
                           (int)S[0]);
      return NULL;
    }

    if (idx_type)
      return wrap_get_ptypes(S, md);
    else
      return vrna_ptypes(S, md);
  } else {
    return NULL;
  }
}

# if defined __USE_POSIX || defined __USE_MISC
#  define strtok_r(s, sep, nextp) __strtok_r (s, sep, nextp)
# endif


#include <stdarg.h>
#define INT_MAX __INT_MAX__
PUBLIC char *
vrna_strdup_vprintf(const char  *format,
                    va_list     argp)
{
  char    *result;
  int     r;

  result = NULL;

#ifndef HAVE_VASPRINTF
  int     count;
  va_list copy;
  va_copy(copy, argp);

  r = -1;

  /* retrieve the number of characters that the string requires */
#ifdef _WIN32
  /*
   * vsnprintf() in Windows is not ANSI compliant, although it's
   * "...included for compliance to the ANSI standard"
   * Thus, we use _vscprintf() that explicitly counts characters
   */
  count = _vscprintf(format, argp);
#else
  count = vsnprintf(NULL, 0, format, argp);
#endif

  if ((count >= 0) && (count < INT_MAX)) {
    char *buf = (char *)vrna_alloc(sizeof(char) * (count + 1));
    if (buf == NULL)
      r = -1;
    else if ((r = vsnprintf(buf, count + 1, format, copy)) < 0)
      free(buf);
    else
      result = buf;
  }

  va_end(copy);  /* Each va_start() or va_copy() needs a va_end() */
#else
  /* the default is to use vasprintf() if available */
  r = vasprintf(&result, format, argp);
#endif

  /* check for any memory allocation error indicated by r == -1 */
  if (r == -1) {
    printf("vrna_strdup_printf: memory allocation failure!");
    result = NULL;
  }

  return result;
}

PUBLIC char *
vrna_strdup_printf(const char *format,
                   ...)
{
  char    *result;
  va_list argp;

  va_start(argp, format);
  result = vrna_strdup_vprintf(format, argp);
  va_end(argp); /* Each va_start() or va_copy() needs a va_end() */

  return result;
}

PUBLIC char **
vrna_strsplit(const char  *string,
              const char  *delimiter)
{
  char          delim[2], *ptr, *ptr2, *token, *save, **split;
  unsigned int  n;

  split = NULL;
  n     = 0;

  if (string) {
    if ((delimiter) && (*delimiter))
      delim[0] = *delimiter;
    else
      delim[0] = '&';

    delim[1] = '\0';

    /* copy string such that we can alter it via strtok() */
    ptr2 = strdup(string);

    /* count how many elements we'll extract */
    ptr = ptr2;

    while (*ptr++)
      if (*ptr == *delim)
        n++;

    /*
     * allocate (n + 1) + 1 elements in split list
     * n + 1 elements plus 1 additional element to indicate
     * the last element in split
     */
    split = (char **)vrna_alloc(sizeof(char *) * (n + 2));

    n     = 0;
#ifdef _WIN32
# ifndef __MINGW32__
    token = strtok_s(ptr2, delim, &save);
# else
    token = strtok_r(ptr2, delim, &save);
# endif
#else
    token = strtok_r(ptr2, delim, &save);
#endif

    while (token != NULL) {
      split[n++]  = vrna_strdup_printf("%s", token);
#ifdef _WIN32
# ifndef __MINGW32__
      token = strtok_s(NULL, delim, &save);
# else
      token = strtok_r(NULL, delim, &save);
# endif
#else
      token = strtok_r(NULL, delim, &save);
#endif
    }

    split[n] = NULL;

    free(ptr2);
  }

  return split;
}

PUBLIC void
vrna_seq_toupper(char *sequence)
{
  unsigned int i;

  if (sequence)
    for (i = 0; sequence[i]; i++)
      sequence[i] = toupper(sequence[i]);
}



PUBLIC short *
vrna_seq_encode_simple(const char *sequence,
                       vrna_md_t  *md)
{
  unsigned int  i, l;
  short         *S = NULL;

  if (sequence && md) {
    l = (unsigned int)strlen(sequence);
    S = (short *)vrna_alloc(sizeof(short) * (l + 2));


    for (i = 1; i <= l; i++) /* make numerical encoding of sequence */
      S[i] = (short)vrna_nucleotide_encode(sequence[i - 1], md);
    S[l + 1]  = S[1];
    S[0]      = (short)l;
  }
  return S;
}

PUBLIC short *
vrna_seq_encode(const char  *sequence,
                vrna_md_t   *md)
{
  unsigned int  i, l;
  short         *S = NULL;

  if (sequence && md) {
    S = vrna_seq_encode_simple(sequence, md);
    l = (unsigned int)strlen(sequence);
    // printf("len is %d\n",l);
    for (i = 1; i <= l; i++)
      S[i] = md->alias[S[i]];

    S[l + 1]  = S[1];
    S[0]      = S[l];
    // for (i = 0; i <= l + 1; i++)
    //   printf("%d-", S[i]);
    // for (i = 0; i <= l + 1; i++)
    //   printf("%d-", S[i]);
  }

  return S;
}



PRIVATE void
set_sequence(vrna_seq_t   *obj,
             const char   *string,
             const char   *name,
             vrna_md_t    *md,
             unsigned int options)
{
  obj->name   = name ? strdup(name) : NULL;
  obj->string = strdup(string);
  vrna_seq_toupper(obj->string);
  obj->length = strlen(obj->string);

  switch (options) {
    default:
      obj->type = VRNA_SEQ_RNA;
  }

  obj->encoding   = vrna_seq_encode(obj->string, md);
  obj->encoding5  = (short *)vrna_alloc(sizeof(short) * (obj->length + 1));
  obj->encoding3  = (short *)vrna_alloc(sizeof(short) * (obj->length + 1));

  if (md->circ) {
    for (size_t i = obj->length; i > 0; i--) {
      if (obj->encoding[i] == 0) /* no nucleotide, i.e. gap */
        continue;

      obj->encoding5[1] = obj->encoding[i];
      break;
    }
    for (size_t i = 1; i <= obj->length; i++) {
      if (obj->encoding[i] == 0) /* no nucleotide, i.e. gap */
        continue;

      obj->encoding3[obj->length] = obj->encoding[i];
      break;
    }
  } else {
    obj->encoding5[1] = obj->encoding3[obj->length] = 0;
  }

  for (size_t i = 1; i < obj->length; i++) {
    if (obj->encoding[i] == 0)
      obj->encoding5[i + 1] = obj->encoding5[i];
    else
      obj->encoding5[i + 1] = obj->encoding[i];
  }

  for (size_t i = obj->length; i > 1; i--) {
    if (obj->encoding[i] == 0)
      obj->encoding3[i - 1] = obj->encoding3[i];
    else
      obj->encoding3[i - 1] = obj->encoding[i];
  }
}


PUBLIC int
vrna_sequence_add(vrna_fold_compound_t  *vc,
                  const char            *string,
                  unsigned int          options)
{
  unsigned int  add_length;
  int           ret = 0;

  if ((vc) &&
      (vc->type == VRNA_FC_TYPE_SINGLE) &&
      (string)) {
    add_length = strlen(string);

    /* add the sequence to the nucleotides container */
    vc->nucleotides = (vrna_seq_t *)vrna_realloc(vc->nucleotides,
                                                 sizeof(vrna_seq_t) *
                                                 (vc->strands + 1));
    set_sequence(&(vc->nucleotides[vc->strands]),
                 string,
                 NULL,
                 &(vc->params->model_details),
                 options);
    /* increase strands counter */
    vc->strands++;

    /* add new sequence to initial order of all strands */
    vc->sequence = (char *)vrna_realloc(vc->sequence,
                                        sizeof(char) *
                                        (vc->length + add_length + 1));
    memcpy(vc->sequence + vc->length,
           vc->nucleotides[vc->strands - 1].string,
           add_length * sizeof(char));
    vc->sequence[vc->length + add_length] = '\0';

    /* add encoding for new strand */
    vc->sequence_encoding = (short *)vrna_realloc(vc->sequence_encoding,
                                                  sizeof(short) *
                                                  (vc->length + add_length + 2));

    memcpy(vc->sequence_encoding + vc->length + 1,
           vc->nucleotides[vc->strands - 1].encoding + 1,
           add_length * sizeof(short));

    /* restore circular encoding */
    vc->sequence_encoding[vc->length + add_length + 1]  = vc->sequence_encoding[1];
    vc->sequence_encoding[0]                            =
      vc->sequence_encoding[vc->length + add_length];

    /* add encoding2 (simple encoding) for new strand */
    vc->sequence_encoding2 = (short *)vrna_realloc(vc->sequence_encoding2,
                                                   sizeof(short) *
                                                   (vc->length + add_length + 2));
    short *enc = vrna_seq_encode_simple(vc->nucleotides[vc->strands - 1].string,
                                        &(vc->params->model_details));
    memcpy(vc->sequence_encoding2 + vc->length + 1,
           enc + 1,
           add_length * sizeof(short));
    free(enc);
    vc->sequence_encoding2[vc->length + add_length + 1] = vc->sequence_encoding2[1];
    vc->sequence_encoding2[0]                           = (short)(vc->length + add_length);

    /* finally, increase length property of the fold compound */
    vc->length = vc->length + add_length;
    
    ret = 1;
  }

  return ret;
}


PUBLIC void
vrna_sequence_prepare(vrna_fold_compound_t *fc)
{
  unsigned int cnt, i;

  if (fc) {
    free(fc->strand_number);
    free(fc->strand_order);
    free(fc->strand_order_uniq);
    free(fc->strand_start);
    free(fc->strand_end);

    fc->strand_order      = NULL;
    fc->strand_order_uniq = NULL;
    fc->strand_start      = NULL;
    fc->strand_end        = NULL;

    fc->strand_number = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->length + 2));

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        /* 1. store initial strand order */
        fc->strand_order_uniq =
          (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));
        fc->strand_order =
          (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));
        for (cnt = 0; cnt < fc->strands; cnt++)
          fc->strand_order[cnt] = cnt;

        /* 2. mark start and end positions of sequences */
        fc->strand_start  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));
        fc->strand_end    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (fc->strands + 1));

        fc->strand_start[0] = 1;
        fc->strand_end[0]   = fc->strand_start[0] + fc->nucleotides[0].length - 1;

        for (cnt = 1; cnt < fc->strands; cnt++) {
          fc->strand_start[cnt] = fc->strand_end[cnt - 1] + 1;
          fc->strand_end[cnt]   = fc->strand_start[cnt] + fc->nucleotides[cnt].length - 1;
          for (i = fc->strand_start[cnt]; i <= fc->strand_end[cnt]; i++)
            fc->strand_number[i] = cnt;
        }
        /* this sets pos. 0 and n + 1 as well for convenience reasons */
        fc->strand_number[0]              = fc->strand_number[1];
        fc->strand_number[fc->length + 1] = fc->strand_number[fc->length];

        break;
    }
  }
}




PRIVATE void
set_fold_compound(vrna_fold_compound_t  *fc,
                  unsigned int          options,
                  unsigned int          aux)
{
  char          *sequence, **sequences, **ptr;
  unsigned int  length, s;
  vrna_md_t     *md_p;

  sequence  = NULL;
  sequences = NULL;

  md_p = &(fc->params->model_details);

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      // printf("VRNA_FC_TYPE_SINGLE");
      sequence = fc->sequence;

      fc->sequence  = NULL;
      fc->length    = 0;

      /* split input sequences at default delimiter '&' */
      sequences = vrna_strsplit(sequence, NULL);

      /* add individual sequences to fold compound */
      for (ptr = sequences; *ptr; ptr++) {
        vrna_sequence_add(fc, *ptr, 1U);
        free(*ptr);
      }

      free(sequences);
      free(sequence);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
      if (fc->strands > 1) {
        fc->cutpoint = fc->nucleotides[0].length + 1;

#if 0
        if (md_p->min_loop_size == TURN)
          md_p->min_loop_size = 0;                              /* is it safe to set this here? */
#endif
      }

#endif

      if (!(options & VRNA_OPTION_EVAL_ONLY)) {
        /* temporary hack for multi-strand case */
        if (fc->strands > 1) {
          int min_loop_size = md_p->min_loop_size;
          md_p->min_loop_size = 0;
          fc->ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc->sequence_encoding2, md_p) : NULL;
          md_p->min_loop_size = min_loop_size;
        } else {
          fc->ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc->sequence_encoding2, md_p) : NULL;
        }
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
        /* backward compatibility ptypes */
        fc->ptype_pf_compat =
          (aux & WITH_PTYPE_COMPAT) ? get_ptypes(fc->sequence_encoding2, md_p, 1) : NULL;
#endif
      }

      break;

    default:                      /* do nothing ? */
      break;
  }

  vrna_sequence_prepare(fc);

  if (!(options & VRNA_OPTION_WINDOW) && (fc->length <= 1000000)) {
    fc->iindx = vrna_idx_row_wise(fc->length);
    fc->jindx = vrna_idx_col_wise(fc->length);
  }
}

PRIVATE void
hc_depot_free(vrna_hc_t *hc)
{
  unsigned int    s, i;
  vrna_hc_depot_t *depot = hc->depot;

  if (depot) {
    if (depot->up) {
      for (s = 0; s < depot->strands; s++)
        free(depot->up[s]);

      free(depot->up);
    }

    if (depot->bp) {
      for (s = 0; s < depot->strands; s++) {
        for (i = 1; i <= depot->bp_size[s]; i++) {
          free(depot->bp[s][i].j);
          free(depot->bp[s][i].strand_j);
          free(depot->bp[s][i].context);
        }
        free(depot->bp[s]);
      }

      free(depot->bp);
    }

    free(depot->bp_size);
    free(depot->up_size);
    free(depot);
  }
  
  hc->depot = NULL;
}


PUBLIC void
vrna_hc_free(vrna_hc_t *hc)
{
  if (hc) {
    if (hc->type == VRNA_HC_DEFAULT)
      free(hc->mx);
    else if (hc->type == VRNA_HC_WINDOW)
      free(hc->matrix_local);

    hc_depot_free(hc);

    free(hc->up_ext);
    free(hc->up_hp);
    free(hc->up_int);
    free(hc->up_ml);

    if (hc->free_data)
      hc->free_data(hc->data);

    free(hc);
  }
}
PRIVATE unsigned char
default_pair_constraint(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j)
{
  unsigned char constraint, can_stack;
  short         *S;
  unsigned int  *sn;
  int           type;
  vrna_md_t     *md;

  sn          = fc->strand_number;
  md          = &(fc->params->model_details);
  constraint  = VRNA_CONSTRAINT_CONTEXT_NONE;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S = fc->sequence_encoding2;
      if (((j - i) < md->max_bp_span) &&
          ((sn[i] != sn[j]) || ((j - i) > md->min_loop_size))) {
        type = md->pair[S[i]][S[j]];
        switch (type) {
          case 0:
            break;
          case 3:
          /* fallthrough */
          case 4:
            if (md->noGU) {
              break;
            } else if (md->noGUclosure) {
              constraint  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              constraint  &= ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
              break;
            }

          /* else fallthrough */
          default:
            constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            break;
        }

        if (md->noLP) {
          // printf("into noLP\n");
          /* check, whether this nucleotide can actually stack with anything or only forms isolated pairs */
          can_stack = VRNA_CONSTRAINT_CONTEXT_NONE;

          /* can it be enclosed by another base pair? */
          if ((i > 1) &&
              (j < fc->length) &&
              (((j - i + 2) < md->max_bp_span) || (sn[i - 1] != sn[j + 1])) &&
              (md->pair[S[i - 1]][S[j + 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          /* can it enclose another base pair? */
          if ((i + 2 < j) &&
              (((j - i - 2) > md->min_loop_size) || (sn[i + 1] != sn[j - 1])) &&
              (md->pair[S[i + 1]][S[j - 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          constraint &= can_stack;
        }
      }
      break;
  }

  return constraint;
}

PRIVATE void
default_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned int  i, n;
  vrna_hc_t     *hc;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n = fc->length;

    for (i = 1; i <= n; i++)
      hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
  }
}

PRIVATE void
default_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned int  i, j, n;
  vrna_hc_t     *hc;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n = fc->length;
    for (j = n; j > 1; j--) {
      for (i = 1; i < j; i++) {
        hc->mx[n * i + j] = default_pair_constraint(fc, i, j);
        hc->mx[n * j + i] = hc->mx[n * i + j];
      }
    }
  }
}


PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc)
{
  vrna_hc_t *hc = vc->hc;

  /*
   * #########################
   * fill with default values
   * #########################
   */

  default_hc_up(vc, 0);

  default_hc_bp(vc, 0);

  /* should we reset the generalized hard constraint feature here? */
  if (hc->f || hc->data) {
    if (hc->free_data)
      hc->free_data(hc->data);

    hc->f         = NULL;
    hc->data      = NULL;
    hc->free_data = NULL;
  }
}


PRIVATE void
hc_update_up(vrna_fold_compound_t *vc)
{
  unsigned int  i, n;
  vrna_hc_t     *hc;

  n   = vc->length;
  hc  = vc->hc;

  if (hc->type == VRNA_HC_WINDOW) {
    /* do nothing for now! */
  } else {
    for (hc->up_ext[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
      hc->up_ext[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 +
                      hc->up_ext[i + 1] : 0;

    for (hc->up_hp[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
      hc->up_hp[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 +
                     hc->up_hp[i + 1] : 0;

    for (hc->up_int[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
      hc->up_int[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 +
                      hc->up_int[i + 1] : 0;

    for (hc->up_ml[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
      hc->up_ml[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 +
                     hc->up_ml[i + 1] : 0;

    /*
     *  loop arround once more until we find a nucleotide that mustn't
     *  be unpaired (needed for circular folding)
     *  Note, circular fold is only possible for single strand predictions
     */
    if (vc->strands < 2) {
      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        hc->up_ext[n + 1] = hc->up_ext[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
            hc->up_ext[i] = MIN2(n, 1 + hc->up_ext[i + 1]);
          else
            break;
        }
      }

      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
        hc->up_hp[n + 1] = hc->up_hp[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP)
            hc->up_hp[i] = MIN2(n, 1 + hc->up_hp[i + 1]);
          else
            break;
        }
      }

      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
        hc->up_int[n + 1] = hc->up_int[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
            hc->up_int[i] = MIN2(n, 1 + hc->up_int[i + 1]);
          else
            break;
        }
      }

      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        hc->up_ml[n + 1] = hc->up_ml[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP)
            hc->up_ml[i] = MIN2(n, 1 + hc->up_ml[i + 1]);
          else
            break;
        }
      }
    }
  }
}

PUBLIC void
vrna_hc_init(vrna_fold_compound_t *vc)
{
  unsigned int  n;
  vrna_hc_t     *hc;

  n = vc->length;

  /* free previous hard constraints */
  vrna_hc_free(vc->hc);

  /* allocate memory new hard constraints data structure */
  hc          = (vrna_hc_t *)vrna_alloc(sizeof(vrna_hc_t));
  hc->type    = VRNA_HC_DEFAULT;
  hc->n       = n;
  hc->mx      = (unsigned char *)vrna_alloc(sizeof(unsigned char) * ((n + 1) * (n + 1) + 1));
  hc->up_ext  = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_hp   = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_int  = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->up_ml   = (int *)vrna_alloc(sizeof(int) * (n + 2));
  hc->depot   = NULL;
  hc->state   = STATE_UNINITIALIZED;

  /* set new hard constraints */
  vc->hc = hc;

  /* prefill default values  */
  hc_reset_to_default(vc);

  /* add null pointers for the generalized hard constraint feature */
  hc->f         = NULL;
  hc->data      = NULL;
  hc->free_data = NULL;

  /* update */
  hc_update_up(vc);
}




#define ALLOC_NOTHING     0
#define ALLOC_F           1
#define ALLOC_F5          2
#define ALLOC_F3          4
#define ALLOC_FC          8
#define ALLOC_C           16
#define ALLOC_FML         32
#define ALLOC_PROBS       256
#define ALLOC_AUX         512
#define VRNA_OPTION_MFE (1 << 0)
#define ALLOC_UNIQ 4096
#define ALLOC_CIRC 1024
#define ALLOC_MULTISTRAND 2048
#define VRNA_OPTION_HYBRID (1 << 2)
#define ALLOC_MFE_LOCAL (ALLOC_F3 | ALLOC_C | ALLOC_FML)
#define ALLOC_MFE_DEFAULT (ALLOC_F5 | ALLOC_C | ALLOC_FML)
#define ALLOC_PF_WO_PROBS (ALLOC_F | ALLOC_C | ALLOC_FML)
#define ALLOC_PF_DEFAULT (ALLOC_PF_WO_PROBS | ALLOC_PROBS | ALLOC_AUX)


PUBLIC int
vrna_mx_add(vrna_fold_compound_t  *vc,
            vrna_mx_type_e        mx_type,
            unsigned int          options)
{
  int ret;

  ret = 1;

  // if (options & VRNA_OPTION_MFE)
  //   printf("into vrna_mx_mfe_add \n");
  //   ret &= vrna_mx_mfe_add(vc, mx_type, options);

  // if (options & VRNA_OPTION_PF)
  //   printf("into vrna_mx_pf_add \n");
  //   ret &= vrna_mx_pf_add(vc, mx_type, options);

  return ret;
}





PUBLIC vrna_fold_compound_t *
vrna_fold_compound(const char       *sequence,
                   const vrna_md_t  *md_p,
                   unsigned int     options)
{
  unsigned int          length, aux_options;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;
  // printf("vrna_fold_compound called, seq is %s", sequence);
  if (sequence == NULL)
    return NULL;

  /* sanity check */
  length = strlen(sequence);

  fc = init_fc_single();

  fc->length    = length;
  fc->sequence  = strdup(sequence);

  aux_options = 0L;


  /* get a copy of the model details */
  if (md_p)
    /* into  here*/
    md = *md_p;
  else
    printf("into vrna_md_set_default\n");
    vrna_md_set_default(&md);

  /* now for the energy parameters */
  add_params(fc, &md, options);

  sanitize_bp_span(fc, options);
  /* only into else*/
  if (options & VRNA_OPTION_WINDOW) {
    set_fold_compound(fc, options, aux_options);
    printf("into (options & VRNA_OPTION_EVAL_ONLY)");

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add minimal hard constraint data structure */
      // vrna_hc_init_window(fc);

      /* add DP matrices */
      // vrna_mx_add(fc, VRNA_MX_WINDOW, options);
    }
  } else {
    /* regular global structure prediction */
    aux_options |= WITH_PTYPE;

    if (options & VRNA_OPTION_PF)
      aux_options |= WITH_PTYPE_COMPAT;

    set_fold_compound(fc, options, aux_options);

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      // printf("into !(options & VRNA_OPTION_EVAL_ONLY)\n");
      /* add default hard constraints */
      vrna_hc_init(fc);

      /* add DP matrices (if required) */
      // vrna_mx_add(fc, VRNA_MX_DEFAULT, options);
    }
  }

  return fc;
}


PRIVATE void
rescale_params(vrna_fold_compound_t *vc)
{
  int               i;
  vrna_exp_param_t  *pf = vc->exp_params;
  vrna_mx_pf_t      *m  = vc->exp_matrices;

  if (m && pf) {
    m->scale[0]     = 1.;
    m->scale[1]     = (FLT_OR_DBL)(1. / pf->pf_scale);
    m->expMLbase[0] = 1;
    m->expMLbase[1] = (FLT_OR_DBL)(pf->expMLbase / pf->pf_scale);
    for (i = 2; i <= vc->length; i++) {
      m->scale[i]     = m->scale[i / 2] * m->scale[i - (i / 2)];
      m->expMLbase[i] = (FLT_OR_DBL)pow(pf->expMLbase, (double)i) * m->scale[i];
    }
  }
}


PUBLIC void
vrna_exp_params_rescale(vrna_fold_compound_t  *vc,
                        double                *mfe)
{
  vrna_exp_param_t  *pf;
  double            e_per_nt, kT;
  vrna_md_t         *md;
  // printf("\n%f\n",*mfe);
  if (vc) {
    if (!vc->exp_params) {
      vc->exp_params = vrna_exp_params(&(vc->params->model_details));
    } else if (memcmp(&(vc->params->model_details),
                      &(vc->exp_params->model_details),
                      sizeof(vrna_md_t)) != 0) {
      /* make sure that model details are matching */
      (void)vrna_md_copy(&(vc->exp_params->model_details), &(vc->params->model_details));
      /* we probably need some mechanism to check whether DP matrices still match the new model settings! */
    }

    pf = vc->exp_params;
    if (pf) {
      kT  = pf->kT;
      md  = &(pf->model_details);

      /* re-compute scaling factor if necessary */
      if ((mfe) || (pf->pf_scale < 1.)) {
        if (mfe)  /* use largest known Boltzmann factor for scaling */
          e_per_nt = *mfe * 1000. / vc->length;
        else      /* use mean energy for random sequences: 184.3*length cal for scaling */
          e_per_nt = -185 + (pf->temperature - 37.) * 7.27;

        /* apply user-defined scaling factor to allow scaling for unusually stable/unstable structure enembles */
        pf->pf_scale = exp(-(md->sfact * e_per_nt) / kT);
      }

      if (pf->pf_scale < 1.)
        pf->pf_scale = 1.;

      rescale_params(vc);
    }
  }
}


#define VRNA_STATUS_PF_PRE (unsigned char)3

PRIVATE void
add_aux_grammar(vrna_fold_compound_t *fc)
{
  fc->aux_grammar = (struct vrna_gr_aux_s *)vrna_alloc(sizeof(struct vrna_gr_aux_s));

  fc->aux_grammar->cb_proc = NULL;

  fc->aux_grammar->cb_aux     = NULL;
  fc->aux_grammar->cb_aux_f   = NULL;
  fc->aux_grammar->cb_aux_c   = NULL;
  fc->aux_grammar->cb_aux_m   = NULL;
  fc->aux_grammar->cb_aux_m1  = NULL;

  fc->aux_grammar->cb_aux_exp     = NULL;
  fc->aux_grammar->cb_aux_exp_f   = NULL;
  fc->aux_grammar->cb_aux_exp_c   = NULL;
  fc->aux_grammar->cb_aux_exp_m   = NULL;
  fc->aux_grammar->cb_aux_exp_m1  = NULL;

  fc->aux_grammar->data       = NULL;
  fc->aux_grammar->free_data  = NULL;
}



PUBLIC int
vrna_gr_set_aux_exp_c(vrna_fold_compound_t      *fc,
                      vrna_grammar_rule_f_exp cb)
{
  int ret = 0;

  if (fc) {
    if (!fc->aux_grammar)
      add_aux_grammar(fc);

    fc->aux_grammar->cb_aux_exp_c = cb;

    ret = 1;
  }

  return ret;
}



#define VRNA_DECOMP_EXT_STEM (unsigned char)14

#define VRNA_CONSTRAINT_FILE      0

/**
 *  @brief  Indicate generation of constraints for MFE folding
 *
 *  @deprecated   This flag has no meaning anymore, since constraints are now always stored! (since v2.2.6)
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_MFE  0

/**
 *  @brief  Indicate generation of constraints for partition function computation
 *  @deprecated   Use #VRNA_OPTION_PF instead!
 *  @ingroup  constraints
 *
 */
#define VRNA_CONSTRAINT_SOFT_PF   VRNA_OPTION_PF

/**
 *  @brief  Flag passed to generic softt constraints callback to indicate hairpin loop decomposition step
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a hairpin loop enclosed by the base pair @f$(i,j)@f$.
 *
 *  @image xml   decomp_hp.svg
 */
#define VRNA_DECOMP_PAIR_HP     (unsigned char)1

/**
 *  @brief  Indicator for interior loop decomposition step
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an interior loop enclosed by the base pair @f$(i,j)@f$,
 *  and enclosing the base pair @f$(k,l)@f$.
 *
 *  @image xml   decomp_il.svg
 */
#define VRNA_DECOMP_PAIR_IL     (unsigned char)2

/**
 *  @brief  Indicator for multibranch loop decomposition step
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop enclosed by the base pair @f$(i,j)@f$,
 *  and consisting of some enclosed multi loop content from k to l.
 *
 *  @image xml   decomp_ml.svg
 *
 */
#define VRNA_DECOMP_PAIR_ML     (unsigned char)3
#define VRNA_DECOMP_PAIR_ML_EXT     (unsigned char)23

#define VRNA_DECOMP_PAIR_ML_OUTSIDE     (unsigned char)4
/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into two multibranch loop parts @f$[i:k]@f$, and @f$[l:j]@f$.
 *
 *  @image xml   decomp_ml_ml_ml.svg
 */
#define VRNA_DECOMP_ML_ML_ML    (unsigned char)5

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be considered a single stem branching off with base pair @f$(k,l)@f$.
 *
 *  @image xml   decomp_ml_stem.svg
 */
#define VRNA_DECOMP_ML_STEM     (unsigned char)6

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into a (usually) smaller multibranch loop part @f$[k:l]@f$.
 *
 *  @image xml   decomp_ml_ml.svg
 */
#define VRNA_DECOMP_ML_ML       (unsigned char)7

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will be considered a multibranch loop part that only consists of unpaired
 *  nucleotides.
 *
 *  @image xml   decomp_ml_up.svg
 */
#define VRNA_DECOMP_ML_UP       (unsigned char)8

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  which will decomposed into a multibranch loop part @f$[i:k]@f$, and a stem with
 *  enclosing base pair @f$(l,j)@f$.
 *
 *  @image xml   decomp_ml_ml_stem.svg
 */
#define VRNA_DECOMP_ML_ML_STEM (unsigned char)9

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  where two stems with enclosing pairs @f$(i,k)@f$ and @f$(l,j)@f$ are coaxially stacking
 *  onto each other.
 *
 *  @image xml   decomp_ml_coaxial.svg
 */
#define VRNA_DECOMP_ML_COAXIAL  (unsigned char)10

/**
 *  @brief  Indicator for decomposition of multibranch loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates a multibranch loop part in the interval @f$[i:j]@f$,
 *  where two stems with enclosing pairs @f$(i,k)@f$ and @f$(l,j)@f$ are coaxially stacking
 *  onto each other.
 *
 *  @image xml   decomp_ml_coaxial.svg
 */
#define VRNA_DECOMP_ML_COAXIAL_ENC  (unsigned char)11

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @def VRNA_DECOMP_EXT_EXT
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into a (usually) smaller exterior loop part @f$[k:l]@f$.
 *
 *  @image xml   decomp_ext_ext.svg
 */
#define VRNA_DECOMP_EXT_EXT     (unsigned char)12

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be considered as an exterior loop component consisting of only unpaired
 *  nucleotides.
 *
 *  @image xml   decomp_ext_up.svg
 */
#define VRNA_DECOMP_EXT_UP      (unsigned char)13

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be considered a stem with enclosing pair @f$(k,l)@f$.
 *
 *  @image xml   decomp_ext_stem.svg
 */
#define VRNA_DECOMP_EXT_STEM (unsigned char)14

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into two exterior loop parts @f$[i:k]@f$ and @f$[l:j]@f$.
 *
 *  @image xml   decomp_ext_ext_ext.svg
 */
#define VRNA_DECOMP_EXT_EXT_EXT (unsigned char)15

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into a stem branching off with base pair @f$(i,k)@f$, and
 *  an exterior loop part @f$[l:j]@f$.
 *
 *  @image xml   decomp_ext_stem_ext.svg
 */
#define VRNA_DECOMP_EXT_STEM_EXT (unsigned char)16

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 */
#define VRNA_DECOMP_EXT_STEM_OUTSIDE (unsigned char)17

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into an exterior loop part @f$[i:k]@f$, and a stem
 *  branching off with base pair @f$(l,j)@f$.
 *
 *  @image xml   decomp_ext_ext_stem.svg
 */
#define VRNA_DECOMP_EXT_EXT_STEM (unsigned char)18

/**
 *  @brief  Indicator for decomposition of exterior loop part
 *
 *  @ingroup  constraints
 *
 *  @def VRNA_DECOMP_EXT_EXT_STEM1
 *  @details This flag notifies the soft or hard constraint callback function that the current
 *  decomposition step evaluates an exterior loop part in the interval @f$[i:j]@f$,
 *  which will be decomposed into an exterior loop part @f$[i:k]@f$, and a stem
 *  branching off with base pair @f$(l,j-1)@f$.
 *
 *  @image xml   decomp_ext_ext_stem1.svg
 */
#define VRNA_DECOMP_EXT_EXT_STEM1 (unsigned char)19

#define VRNA_DECOMP_EXT_STEM_EXT1 (unsigned char)20

#define VRNA_DECOMP_EXT_L         (unsigned char)21
#define VRNA_DECOMP_EXT_EXT_L     (unsigned char)22

/*
 * currently we do not allow for more than 31 different decomposition types
 * This must be changed as soon as the above macros turn to values above 32
 */
#define VRNA_DECOMP_TYPES_MAX     32

PRIVATE unsigned char
hc_ext_cb_def(int           i,
              int           j,
              int           k,
              int           l,
              unsigned char d,
              void          *data)
{
  int                   di, dj;
  unsigned char         eval;
  unsigned int          n;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  n     = dat->n;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      if (dat->mx[n * j + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      if (dat->mx[n * k + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (dat->mx[n * (j - 1) + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (dat->hc_up[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (dat->mx[n * k + i + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;

          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      eval  = (unsigned char)1;
      di    = l - k - 1;
      if ((di != 0) && (dat->hc_up[k + 1] < di))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_STEM:
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
        eval = (unsigned char)1;

      break;

    default:
      printf("hc_cb@exterior_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_sn(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             void           *data)
{
  unsigned int          *sn;
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  sn    = dat->sn;
  eval  = (unsigned char)0;

  /* for now with the 'old' cofold implementation, we allow for any decomposition */
  //return (unsigned char)1;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (sn[j - 1] != sn[j])
        break;

      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (sn[i] != sn[i + 1])
        break;

      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_EXT_STEM:
    /* fall through */
    case VRNA_DECOMP_EXT_STEM_EXT:
    /* fall through */
    case VRNA_DECOMP_EXT_EXT_EXT:
      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_STEM:
    /* fall through */
    case VRNA_DECOMP_EXT_EXT:
      if ((sn[i] == sn[k]) && (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_UP:
      if (sn[i] == sn[j])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      if (((k <= i) || sn[k - 1] == sn[k]) &&
          ((l >= j) || sn[l + 1] == sn[l]))
        eval = (unsigned char)1;

      break;

    default:
      printf("hc_cb@exterior_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_sn(int            i,
                 int            j,
                 int            k,
                 int            l,
                 unsigned char  d,
                 void           *data)
{
  unsigned char eval;

  eval  = hc_ext_cb_def(i, j, k, l, d, data);
  eval  = hc_ext_cb_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   void           *data)
{
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = hc_ext_cb_def(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_sn_user(int           i,
                      int           j,
                      int           k,
                      int           l,
                      unsigned char d,
                      void          *data)
{
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = hc_ext_cb_def(i, j, k, l, d, data);
  eval  = hc_ext_cb_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;
  eval  = dat->hc_f(i, j, k, l, d, dat->hc_dat) ? eval : (unsigned char)0;

  return eval;
}



PRIVATE INLINE vrna_hc_eval_f
prepare_hc_ext_def(vrna_fold_compound_t   *fc,
                   struct hc_ext_def_dat  *dat)
{
  dat->mx     = fc->hc->mx;
  dat->n      = fc->length;
  dat->hc_up  = fc->hc->up_ext;
  dat->sn     = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return (fc->strands == 1) ? &hc_ext_cb_def_user : &hc_ext_cb_def_sn_user;
  }

  return (fc->strands == 1) ? &hc_ext_cb_def : &hc_ext_cb_def_sn;
}


PUBLIC unsigned int
vrna_get_ptype_md(int       i,
                  int       j,
                  vrna_md_t *md)
{
  unsigned int tt = (unsigned int)md->pair[i][j];

  return (tt == 0) ? 7 : tt;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ext_stem(unsigned int      type,
                    int               n5d,
                    int               n3d,
                    vrna_exp_param_t  *p)
{
  double energy = 1.0;

  if (n5d >= 0 && n3d >= 0)
    energy = p->expmismatchExt[type][n5d][n3d];
  else if (n5d >= 0)
    energy = p->expdangle5[type][n5d];
  else if (n3d >= 0)
    energy = p->expdangle3[type][n3d];

  if (type > 2)
    energy *= p->expTermAU;

  return (FLT_OR_DBL)energy;
}




PRIVATE unsigned char
hc_ext_cb_def_window(int            i,
                     int            j,
                     int            k,
                     int            l,
                     unsigned char  d,
                     void           *data)
{
  int                   di, dj;
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      if (dat->mx_window[l][j - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      if (dat->mx_window[i][k - i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (j != k) {
          /* otherwise, stem spans from i to j */
          dj = l - k - 1;
          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (dat->mx_window[l][j - 1 - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (dat->mx_window[i + 1][k - (i + 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;

          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      eval  = (unsigned char)1;
      di    = l - k - 1;
      if ((di != 0) && (dat->hc_up[k + 1] < di))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    default:
      printf("hc_cb@exterior_loops.c: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_def_user_window(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          void          *data)
{
  unsigned char         eval;
  struct hc_ext_def_dat *dat = (struct hc_ext_def_dat *)data;

  eval  = hc_ext_cb_def_window(i, j, k, l, d, data);
  eval  = dat->hc_f(i, j, k, l, d, dat->hc_dat) ? eval : (unsigned char)0;

  return eval;
}

PRIVATE INLINE vrna_hc_eval_f
prepare_hc_ext_def_window(vrna_fold_compound_t  *fc,
                          struct hc_ext_def_dat *dat)
{
  dat->mx_window  = fc->hc->matrix_local;
  dat->hc_up      = fc->hc->up_ext;
  dat->sn         = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_ext_cb_def_user_window;
  }

  return &hc_ext_cb_def_window;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_red(int                   i,
                  int                   j,
                  int                   k,
                  int                   l,
                  struct sc_ext_exp_dat *data)
{
  unsigned int  start_2, length_1, length_2;
  FLT_OR_DBL    q_sc, **sc_up;

  sc_up = data->up;

  q_sc = 1.;

  length_1  = k - i;
  start_2   = l + 1;
  length_2  = j - l;

  if (length_1 != 0)
    q_sc *= sc_up[i][length_1];

  if (length_2 != 0)
    q_sc *= sc_up[start_2][length_2];

  return q_sc;
}
PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_red_user_to_ext(int                   i,
                              int                   j,
                              int                   k,
                              int                   l,
                              struct sc_ext_exp_dat *data)
{
  return data->user_cb(i, j, k, l, VRNA_DECOMP_EXT_EXT, data->user_data);
}
PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_red_user_def_to_ext(int                   i,
                                  int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_ext_exp_dat *data)
{
  return sc_ext_exp_cb_red(i, j, k, l, data) *
         sc_ext_exp_cb_red_user_to_ext(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_red_user_to_stem(int                    i,
                               int                    j,
                               int                    k,
                               int                    l,
                               struct sc_ext_exp_dat  *data)
{
  return data->user_cb(i, j, k, l, VRNA_DECOMP_EXT_STEM, data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_red_user_def_to_stem(int                    i,
                                   int                    j,
                                   int                    k,
                                   int                    l,
                                   struct sc_ext_exp_dat  *data)
{
  return sc_ext_exp_cb_red(i, j, k, l, data) *
         sc_ext_exp_cb_red_user_to_stem(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_up(int                    i,
                 int                    j,
                 struct sc_ext_exp_dat  *data)
{
  unsigned int  length;
  FLT_OR_DBL    q_sc, **sc_up;

  sc_up   = data->up;
  length  = j - i + 1;
  q_sc    = 1.;

  if (length != 0)
    q_sc *= sc_up[i][length];

  return q_sc;
}

PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_up_user(int                   i,
                      int                   j,
                      struct sc_ext_exp_dat *data)
{
  return data->user_cb(i, j, i, j, VRNA_DECOMP_EXT_UP, data->user_data);
}

PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_up_user_def(int                   i,
                          int                   j,
                          struct sc_ext_exp_dat *data)
{
  return sc_ext_exp_cb_up(i, j, data) *
         sc_ext_exp_cb_up_user(i, j, data);
}

PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_cb_split_user(int                    i,
                         int                    j,
                         int                    k,
                         struct sc_ext_exp_dat  *data)
{
  return data->user_cb(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, data->user_data);
}


PRIVATE INLINE void
init_sc_ext_exp(vrna_fold_compound_t  *fc,
                struct sc_ext_exp_dat *sc_wrapper)
{
  vrna_sc_t *sc, **scs;

  sc_wrapper->up                    = NULL;
  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->n_seq                 = 1;
  sc_wrapper->a2s                   = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  /* no soft constraints by default */
  sc_wrapper->red_ext   = NULL;
  sc_wrapper->red_stem  = NULL;
  sc_wrapper->red_up    = NULL;
  sc_wrapper->split     = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        sc_wrapper->up        = sc->exp_energy_up;
        sc_wrapper->user_cb   = sc->exp_f;
        sc_wrapper->user_data = sc->data;

        /* bind correct wrapper functions */
        if (sc->exp_energy_up) {
          if (sc->exp_f) {
            sc_wrapper->red_ext   = &sc_ext_exp_cb_red_user_def_to_ext;
            sc_wrapper->red_stem  = &sc_ext_exp_cb_red_user_def_to_stem;
            sc_wrapper->red_up    = &sc_ext_exp_cb_up_user_def;
            sc_wrapper->split     = &sc_ext_exp_cb_split_user;
          } else {
            sc_wrapper->red_ext   = &sc_ext_exp_cb_red;
            sc_wrapper->red_stem  = &sc_ext_exp_cb_red;
            sc_wrapper->red_up    = &sc_ext_exp_cb_up;
          }
        } else if (sc->exp_f) {
          sc_wrapper->red_ext   = &sc_ext_exp_cb_red_user_to_ext;
          sc_wrapper->red_stem  = &sc_ext_exp_cb_red_user_to_stem;
          sc_wrapper->red_up    = &sc_ext_exp_cb_up_user;
          sc_wrapper->split     = &sc_ext_exp_cb_split_user;
        }
      }

      break;
  }
}


struct vrna_mx_pf_aux_el_s {
  FLT_OR_DBL  *qq;
  FLT_OR_DBL  *qq1;

  int         qqu_size;
  FLT_OR_DBL  **qqu;
};
typedef struct vrna_mx_pf_aux_el_s *vrna_mx_pf_aux_el_t;
#define VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP 1U

PRIVATE INLINE FLT_OR_DBL
reduce_ext_up_fast(vrna_fold_compound_t       *fc,
                   int                        i,
                   int                        j,
                   struct vrna_mx_pf_aux_el_s *aux_mx,
                   vrna_hc_eval_f  evaluate,
                   struct hc_ext_def_dat      *hc_dat_local,
                   struct sc_ext_exp_dat      *sc_wrapper)
{
  int               u;
  FLT_OR_DBL        qbt, q_temp, *scale;
  vrna_ud_t         *domains_up;
  sc_ext_exp_red_up sc_red_up;

  sc_red_up = sc_wrapper->red_up;

  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;
  qbt         = 0.;

  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_UP, hc_dat_local)) {
    u       = j - i + 1;
    q_temp  = scale[u];

    if (sc_red_up)
      q_temp *= sc_red_up(i, j, sc_wrapper);

    qbt += q_temp;

    if ((domains_up) && (domains_up->exp_energy_cb)) {
      qbt += q_temp *
             domains_up->exp_energy_cb(fc,
                                       i, j,
                                       VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                       domains_up->data);
    }
  }

  return qbt;
}

PUBLIC struct vrna_mx_pf_aux_el_s *
vrna_exp_E_ext_fast_init(vrna_fold_compound_t *fc)
{
  struct vrna_mx_pf_aux_el_s *aux_mx = NULL;

  if (fc) {
    unsigned int              u;
    int                       i, j, max_j, d, n, turn, ij, *iidx, with_ud;
    FLT_OR_DBL                *q, **q_local;
    vrna_hc_eval_f evaluate;
    struct hc_ext_def_dat     hc_dat_local;
    struct sc_ext_exp_dat     sc_wrapper;
    vrna_ud_t                 *domains_up;

    n           = (int)fc->length;
    iidx        = fc->iindx;
    turn        = fc->exp_params->model_details.min_loop_size;
    domains_up  = fc->domains_up;
    with_ud     = (domains_up && domains_up->exp_energy_cb);

    if (fc->hc->type == VRNA_HC_WINDOW)
      evaluate = prepare_hc_ext_def_window(fc, &hc_dat_local);
    else
      evaluate = prepare_hc_ext_def(fc, &hc_dat_local);

    init_sc_ext_exp(fc, &sc_wrapper);

    /* allocate memory for helper arrays */
    aux_mx =
      (struct vrna_mx_pf_aux_el_s *)vrna_alloc(sizeof(struct vrna_mx_pf_aux_el_s));
    aux_mx->qq        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qq1       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqu_size  = 0;
    aux_mx->qqu       = NULL;

    /* pre-processing ligand binding production rule(s) and auxiliary memory */
    if (with_ud) {
      int ud_max_size = 0;
      for (u = 0; u < domains_up->uniq_motif_count; u++)
        if (ud_max_size < domains_up->uniq_motif_size[u])
          ud_max_size = domains_up->uniq_motif_size[u];

      aux_mx->qqu_size  = ud_max_size;
      aux_mx->qqu       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));

      for (u = 0; u <= ud_max_size; u++)
        aux_mx->qqu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    }

    if (fc->hc->type == VRNA_HC_WINDOW) {
      q_local = fc->exp_matrices->q_local;
      max_j   = MIN2(turn + 1, fc->window_size);
      max_j   = MIN2(max_j, n);
      for (j = 1; j <= max_j; j++)
        for (i = 1; i <= j; i++)
          q_local[i][j] =
            reduce_ext_up_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);
    } else {
      q = fc->exp_matrices->q;
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;

          q[ij] = reduce_ext_up_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);
        }

      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_f)) {
        for (d = 0; d <= turn; d++)
          for (i = 1; i <= n - d; i++) {
            j   = i + d;
            ij  = iidx[i] - j;

            q[ij] += fc->aux_grammar->cb_aux_exp_f(fc, i, j, fc->aux_grammar->data);
          }
      }
    }
  }

  return aux_mx;
}


PRIVATE FLT_OR_DBL
mf_rule_pair(vrna_fold_compound_t *fc,
             int                  i,
             int                  j,
             void                 *data)
{
  short                     *S1, *S2, s5, s3;
  unsigned int              *sn, *ends, type, nick;
  int                       *my_iindx;
  FLT_OR_DBL                contribution, *q, *scale, qbase, tmp, tmp2;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;
  vrna_sc_t                 *sc;

  contribution  = 0;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  sn            = fc->strand_number;
  ends          = fc->strand_end;
  q             = fc->exp_matrices->q;
  scale         = fc->exp_matrices->scale;
  my_iindx      = fc->iindx;
  sc            = fc->sc;
  evaluate      = prepare_hc_ext_def(fc, &hc_dat_local);

  if ((sn[i] != sn[j]) &&
      (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local))) {
    /* most obious strand nick is at end of sn[i] and start of sn[j] */
    type  = vrna_get_ptype_md(S2[j], S2[i], md);
    s5    = (sn[j] == sn[j - 1]) ? S1[j - 1] : -1;
    s3    = (sn[i] == sn[i + 1]) ? S1[i + 1] : -1;
    qbase = vrna_exp_E_ext_stem(type, s5, s3, pf_params) *
            scale[2];

    if (sc) {
      if (sc->exp_f)
        qbase *= sc->exp_f(j, i, j, i, VRNA_DECOMP_EXT_STEM, sc->data);
    }

    tmp = 0.;

    /*
     *  if (evaluate(i + 1,
     *               j - 1,
     *               ends[sn[i]],
     *               ends[sn[i]] + 1,
     *               VRNA_DECOMP_EXT_EXT_EXT,
     *               &hc_dat_local))
     */
    if (sn[i] != sn[i + 1]) {
      if ((sn[j - 1] != sn[j]) &&
          (i + 1 == j))
        tmp = 1.;
      else if (sn[j - 1] == sn[j])
        tmp = q[my_iindx[i + 1] - j + 1];
    } else if (sn[j - 1] != sn[j]) {
      tmp = q[my_iindx[i + 1] - j + 1];
    } else {
      tmp = q[my_iindx[i + 1] - ends[sn[i]]] *
            q[my_iindx[ends[sn[i]] + 1] - j + 1];

      /* check whether we find more strand nicks between i and j */
      nick = ends[sn[i]] + 1;
      while (sn[nick] != sn[j]) {
        /*
         *      if (evaluate(i + 1,
         *                   j - 1,
         *                   ends[sn[nick]],
         *                   ends[sn[nick]] + 1,
         *                   VRNA_DECOMP_EXT_EXT_EXT,
         *                   &hc_dat_local))
         */
        tmp2 = 1.;
        if (i + 1 <= ends[sn[nick]])
          tmp2 *= q[my_iindx[i + 1] - ends[sn[nick]]];

        if (ends[sn[nick]] + 1 <= j - 1)
          tmp2 *= q[my_iindx[ends[sn[nick]] + 1] - j + 1];

        tmp += tmp2;


        nick = ends[sn[nick]] + 1;
      }
    }

    contribution = qbase *
                   tmp;
  }

  return contribution;
}



PUBLIC int
vrna_pf_multifold_prepare(vrna_fold_compound_t *fc)
{
  if (fc)
    return vrna_gr_set_aux_exp_c(fc, &mf_rule_pair);

  return 0;
}


struct vrna_mx_pf_aux_ml_s {
  FLT_OR_DBL  *qqm;
  FLT_OR_DBL  *qqm1;

  int         qqmu_size;
  FLT_OR_DBL  **qqmu;
};

typedef struct vrna_mx_pf_aux_ml_s *vrna_mx_pf_aux_ml_t;
#define FLT_MAX __FLT_MAX__
#define DBL_MAX __DBL_MAX__

PUBLIC struct vrna_mx_pf_aux_ml_s *
vrna_exp_E_ml_fast_init(vrna_fold_compound_t *fc)
{
  struct vrna_mx_pf_aux_ml_s *aux_mx = NULL;

  if (fc) {
    int         i, j, d, n, u, turn, ij, *iidx;
    FLT_OR_DBL  *qm;

    n     = (int)fc->length;
    iidx  = fc->iindx;
    turn  = fc->exp_params->model_details.min_loop_size;
    qm    = fc->exp_matrices->qm;

    /* allocate memory for helper arrays */
    aux_mx =
      (struct vrna_mx_pf_aux_ml_s *)vrna_alloc(sizeof(struct vrna_mx_pf_aux_ml_s));
    aux_mx->qqm       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqm1      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqmu_size = 0;
    aux_mx->qqmu      = NULL;

    if (fc->type == VRNA_FC_TYPE_SINGLE) {
      vrna_ud_t *domains_up = fc->domains_up;
      int       with_ud     = (domains_up && domains_up->exp_energy_cb);
      int       ud_max_size = 0;

      /* pre-processing ligand binding production rule(s) and auxiliary memory */
      if (with_ud) {
        for (u = 0; u < domains_up->uniq_motif_count; u++)
          if (ud_max_size < domains_up->uniq_motif_size[u])
            ud_max_size = domains_up->uniq_motif_size[u];

        aux_mx->qqmu_size = ud_max_size;
        aux_mx->qqmu      = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));
        for (u = 0; u <= ud_max_size; u++)
          aux_mx->qqmu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
      }
    }

    if (fc->hc->type == VRNA_HC_WINDOW) {
    } else {
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;

          if (j > n)
            continue;

          qm[ij] = 0.;
        }

      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m)) {
        for (d = 0; d <= turn; d++)
          for (i = 1; i <= n - d; i++) {
            j   = i + d;
            ij  = iidx[i] - j;

            if (j > n)
              continue;

            qm[ij] += fc->aux_grammar->cb_aux_exp_m(fc, i, j, fc->aux_grammar->data);
          }
      }
    }
  }

  return aux_mx;
}





PRIVATE unsigned char
hc_hp_cb_def_window(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    void          *data)
{
  int                   u;
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

  eval = (unsigned char)0;

  u = j - i - 1;

  if (dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    eval = (unsigned char)1;
    if (dat->hc_up[i + 1] < u)
      eval = (unsigned char)0;
  }

  return eval;
}



PRIVATE unsigned char
hc_hp_cb_def_user_window(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data)
{
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

  eval  = hc_hp_cb_def_window(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}

PRIVATE INLINE vrna_hc_eval_f
prepare_hc_hp_def_window(vrna_fold_compound_t *fc,
                         struct hc_hp_def_dat *dat)
{
  dat->mx_window  = fc->hc->matrix_local;
  dat->hc_up      = fc->hc->up_hp;
  dat->n          = fc->length;
  dat->sn         = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_hp_cb_def_user_window;
  }

  return &hc_hp_cb_def_window;
}


PRIVATE unsigned char
hc_hp_cb_def(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             void           *data)
{
  int                   u, p, q;
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

  eval = (char)0;

  /* no strand nicks are allowed in hairpin loops */
  if (dat->sn[i] != dat->sn[j])
    return eval;

  if (j > i) {
    /* linear case */
    p = i;
    q = j;
    u = q - p - 1;
  } else {
    /* circular case */
    p = j;
    q = i;
    u = dat->n - q + p - 1;
  }

  if (dat->mx[dat->n * p + q] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    eval = (unsigned char)1;
    if (dat->hc_up[i + 1] < u)
      eval = (unsigned char)0;
  }

  return eval;
}




PRIVATE unsigned char
hc_hp_cb_def_user(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data)
{
  unsigned char         eval;
  struct hc_hp_def_dat  *dat = (struct hc_hp_def_dat *)data;

  eval  = hc_hp_cb_def(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_hp_def(vrna_fold_compound_t  *fc,
                  struct hc_hp_def_dat  *dat)
{
  dat->mx     = fc->hc->mx;
  dat->hc_up  = fc->hc->up_hp;
  dat->n      = fc->length;
  dat->sn     = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_hp_cb_def_user;
  }

  return &hc_hp_cb_def;
}




PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_up(int                   i,
                int                   j,
                struct sc_hp_exp_dat  *data)
{
  return data->up[i + 1][j - i - 1];
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_bp(int                   i,
                int                   j,
                struct sc_hp_exp_dat  *data)
{
  return data->bp[data->idx[j] + i];
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_bp_local(int                   i,
                      int                   j,
                      struct sc_hp_exp_dat  *data)
{
  return data->bp_local[i][j - i];
}


PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_user(int                   i,
                  int                   j,
                  struct sc_hp_exp_dat  *data)
{
  return data->user_cb(i, j, i, j,
                       VRNA_DECOMP_PAIR_HP,
                       data->user_data);
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_up_bp(int                  i,
                   int                  j,
                   struct sc_hp_exp_dat *data)
{
  return sc_hp_exp_cb_up(i, j, data) *
         sc_hp_exp_cb_bp(i, j, data);
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_up_bp_local(int                  i,
                         int                  j,
                         struct sc_hp_exp_dat *data)
{
  return sc_hp_exp_cb_up(i, j, data) *
         sc_hp_exp_cb_bp_local(i, j, data);
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_up_user(int                  i,
                     int                  j,
                     struct sc_hp_exp_dat *data)
{
  return sc_hp_exp_cb_up(i, j, data) *
         sc_hp_exp_cb_user(i, j, data);
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_bp_user(int                  i,
                     int                  j,
                     struct sc_hp_exp_dat *data)
{
  return sc_hp_exp_cb_bp(i, j, data) *
         sc_hp_exp_cb_user(i, j, data);
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_bp_local_user(int                  i,
                           int                  j,
                           struct sc_hp_exp_dat *data)
{
  return sc_hp_exp_cb_bp_local(i, j, data) *
         sc_hp_exp_cb_user(i, j, data);
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_up_bp_user(int                   i,
                        int                   j,
                        struct sc_hp_exp_dat  *data)
{
  return sc_hp_exp_cb_up(i, j, data) *
         sc_hp_exp_cb_bp(i, j, data) *
         sc_hp_exp_cb_user(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_up_bp_local_user(int                   i,
                              int                   j,
                              struct sc_hp_exp_dat  *data)
{
  return sc_hp_exp_cb_up(i, j, data) *
         sc_hp_exp_cb_bp_local(i, j, data) *
         sc_hp_exp_cb_user(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_ext_up(int                   i,
                    int                   j,
                    struct sc_hp_exp_dat  *data)
{
  FLT_OR_DBL  sc;
  int         u1, u2;

  u1  = data->n - j;
  u2  = i - 1;
  sc  = 1.;

  if (u1 > 0)
    sc *= data->up[j + 1][u1];

  if (u2 > 0)
    sc *= data->up[1][u2];

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_ext_user(int                   i,
                      int                   j,
                      struct sc_hp_exp_dat  *data)
{
  return data->user_cb(j, i, j, i,
                       VRNA_DECOMP_PAIR_HP,
                       data->user_data);
}

PRIVATE INLINE FLT_OR_DBL
sc_hp_exp_cb_ext_up_user(int                  i,
                         int                  j,
                         struct sc_hp_exp_dat *data)
{
  return sc_hp_exp_cb_ext_up(i, j, data) *
         sc_hp_exp_cb_ext_user(i, j, data);
}



PRIVATE INLINE void
init_sc_hp_exp(vrna_fold_compound_t *fc,
               struct sc_hp_exp_dat *sc_wrapper)
{
  unsigned char sliding_window;
  unsigned int  s;
  vrna_sc_t     *sc, **scs;

  if (fc->exp_matrices)
    sliding_window = (fc->exp_matrices->type == VRNA_MX_WINDOW) ? 1 : 0;
  else if ((fc->type == VRNA_FC_TYPE_SINGLE) && (fc->sc))
    sliding_window = (fc->sc->type == VRNA_SC_WINDOW) ? 1 : 0;
  else if (fc->hc)
    sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  else
    sliding_window = 0;

  sc_wrapper->n     = (int)fc->length;
  sc_wrapper->idx   = fc->jindx;
  sc_wrapper->n_seq = 1;
  sc_wrapper->a2s   = NULL;

  sc_wrapper->up                    = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->bp                    = NULL;
  sc_wrapper->bp_comparative        = NULL;
  sc_wrapper->bp_local              = NULL;
  sc_wrapper->bp_local_comparative  = NULL;

  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  sc_wrapper->pair      = NULL;
  sc_wrapper->pair_ext  = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        unsigned int provides_sc_up, provides_sc_bp, provides_sc_user;

        provides_sc_up    = 0;
        provides_sc_bp    = 0;
        provides_sc_user  = 0;

        sc_wrapper->up        = sc->exp_energy_up;
        sc_wrapper->bp        = (sliding_window) ? NULL : sc->exp_energy_bp;
        sc_wrapper->bp_local  = (sliding_window) ? sc->exp_energy_bp_local : NULL;
        sc_wrapper->user_cb   = sc->exp_f;
        sc_wrapper->user_data = sc->data;

        if (sc->exp_energy_up)
          provides_sc_up = 1;

        if (sliding_window) {
          if (sc->exp_energy_bp_local)
            provides_sc_bp = 1;
        } else if (sc->exp_energy_bp) {
          provides_sc_bp = 1;
        }

        if (sc->exp_f)
          provides_sc_user = 1;

        if (provides_sc_user) {
          sc_wrapper->pair_ext = &sc_hp_exp_cb_ext_user;
          if (provides_sc_up) {
            sc_wrapper->pair_ext = &sc_hp_exp_cb_ext_up_user;

            if (provides_sc_bp)
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_hp_exp_cb_up_bp_local_user :
                                 &sc_hp_exp_cb_up_bp_user;
            else
              sc_wrapper->pair = &sc_hp_exp_cb_up_user;
          } else if (provides_sc_bp) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_hp_exp_cb_bp_local_user :
                               &sc_hp_exp_cb_bp_user;
          } else {
            sc_wrapper->pair = &sc_hp_exp_cb_user;
          }
        } else {
          if (provides_sc_up) {
            sc_wrapper->pair_ext = &sc_hp_exp_cb_ext_up;
            if (provides_sc_bp) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_hp_exp_cb_up_bp_local :
                                 &sc_hp_exp_cb_up_bp;
            } else {
              sc_wrapper->pair = &sc_hp_exp_cb_up;
            }
          } else if (provides_sc_bp) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_hp_exp_cb_bp_local :
                               &sc_hp_exp_cb_bp;
          }
        }
      }

      break;
 }
}



PRIVATE INLINE FLT_OR_DBL
exp_E_Hairpin(int               u,
              int               type,
              short             si1,
              short             sj1,
              const char        *string,
              vrna_exp_param_t  *P)
{
  double q, kT, salt_correction;

  kT = P->kT;   /* kT in cal/mol  */
  salt_correction = 1.;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    if (u<=MAXLOOP)
      salt_correction = P->expSaltLoop[u+1];
    else
      salt_correction = exp(-vrna_salt_loop_int(u+1, P->model_details.salt, P->temperature+K0, P->model_details.backbone_length) * 10. / kT);
  }

  if (u <= 30)
    q = P->exphairpin[u];
  else
    q = P->exphairpin[30] * exp(-(P->lxc * log(u / 30.)) * 10. / kT);

  q *= salt_correction;

  if (u < 3)
    return (FLT_OR_DBL)q;         /* should only be the case when folding alignments */

  if ((string) && (P->model_details.special_hp)) {
    if (u == 4) {
      char tl[7] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 6);
      tl[6] = '\0';
      if ((ts = strstr(P->Tetraloops, tl))) {
        if (type != 7)
          return (FLT_OR_DBL)(P->exptetra[(ts - P->Tetraloops) / 7] * salt_correction);
        else
          q *= P->exptetra[(ts - P->Tetraloops) / 7];
      }
    } else if (u == 6) {
      char tl[9] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 8);
      tl[8] = '\0';
      if ((ts = strstr(P->Hexaloops, tl)))
        return (FLT_OR_DBL)(P->exphex[(ts - P->Hexaloops) / 9] * salt_correction);
    } else if (u == 3) {
      char tl[6] = {
        0
      }, *ts;
      memcpy(tl, string, sizeof(char) * 5);
      tl[5] = '\0';
      if ((ts = strstr(P->Triloops, tl)))
        return (FLT_OR_DBL)(P->exptri[(ts - P->Triloops) / 6] * salt_correction);

      if (type > 2)
        return (FLT_OR_DBL)(q * P->expTermAU);
      else
        return (FLT_OR_DBL)q;
    }
  }

  q *= P->expmismatchH[type][si1][sj1];

  return (FLT_OR_DBL)q;
}


PRIVATE INLINE void
free_sc_hp_exp(struct sc_hp_exp_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->bp_comparative);
  free(sc_wrapper->bp_local_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}



#define VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP 2U
PRIVATE FLT_OR_DBL
exp_eval_hp_loop(vrna_fold_compound_t *fc,
                 int                  i,
                 int                  j)
{
  char                  **Ss;
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  unsigned int          *sn;
  int                   u, type, n_seq, s;
  FLT_OR_DBL            q, qbt1, *scale;
  vrna_exp_param_t      *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct sc_hp_exp_dat  sc_wrapper;

  P           = fc->exp_params;
  md          = &(P->model_details);
  sn          = fc->strand_number;
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;

  init_sc_hp_exp(fc, &sc_wrapper);
  q = 0.;

  if (sn[j] != sn[i])
    return 0; //exp_eval_hp_loop_fake(fc, i, j);

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      u     = j - i - 1;
      type  = vrna_get_ptype_md(S2[i], S2[j], md);

      if (sn[j] == sn[i]) {
        /* regular hairpin loop */
        q = exp_E_Hairpin(u, type, S[i + 1], S[j - 1], fc->sequence + i - 1, P);
      } else {
        /*
         * hairpin-like exterior loop (for cofolding)
         * this is currently handle somewhere else
         */
      }

      break;
  }

  /* add soft constraints */
  if (sc_wrapper.pair)
    q *= sc_wrapper.pair(i, j, &sc_wrapper);

  if (domains_up && domains_up->exp_energy_cb) {
    /* we always consider both, bound and unbound state */
    q += q * domains_up->exp_energy_cb(fc,
                                       i + 1, j - 1,
                                       VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                       domains_up->data);
  }

  q *= scale[j - i + 1];

  free_sc_hp_exp(&sc_wrapper);

  return q;
}


PRIVATE FLT_OR_DBL
exp_eval_ext_hp_loop(vrna_fold_compound_t *fc,
                     int                  i,
                     int                  j)
{
  char                  **Ss, *sequence, loopseq[10] = {
    0
  };
  unsigned int          **a2s;
  short                 *S, *S2, **SS, **S5, **S3;
  int                   u1, u2, n, type, n_seq, s, noGUclosure;
  FLT_OR_DBL            q, qbt1, *scale;
  vrna_exp_param_t      *P;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  struct sc_hp_exp_dat  sc_wrapper;

  n           = fc->length;
  P           = fc->exp_params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;

  init_sc_hp_exp(fc, &sc_wrapper);

  q   = 0.;
  u1  = n - j;
  u2  = i - 1;

  if ((u1 + u2) < 3)
    return q;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sequence  = fc->sequence;
      S         = fc->sequence_encoding;
      S2        = fc->sequence_encoding2;
      type      = vrna_get_ptype_md(S2[j], S2[i], md);

      if (((type == 3) || (type == 4)) && noGUclosure)
        return q;

      /* get the loop sequence */
      if ((u1 + u2) < 7) {
        memcpy(loopseq, sequence + j - 1, sizeof(char) * (u1 + 1));
        memcpy(loopseq + u1 + 1, sequence, sizeof(char) * (u2 + 1));
        loopseq[u1 + u2 + 2] = '\0';
      }

      q = exp_E_Hairpin(u1 + u2, type, S[j + 1], S[i - 1], loopseq, P);

      break;

    default:
      break;
  }

  /* add soft constraints */
  if (sc_wrapper.pair_ext)
    q *= sc_wrapper.pair_ext(i, j, &sc_wrapper);

  if (domains_up && domains_up->exp_energy_cb) {
    /* we always consider both, bound and unbound state */
    q += q * domains_up->exp_energy_cb(fc,
                                       j + 1, i - 1,
                                       VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                       domains_up->data);
  }

  q *= scale[u1 + u2];

  free_sc_hp_exp(&sc_wrapper);

  return q;
}



PUBLIC FLT_OR_DBL
vrna_exp_E_hp_loop(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j)
{
  vrna_hc_eval_f evaluate;
  struct hc_hp_def_dat      hc_dat_local;

  if (fc->hc->type == VRNA_HC_WINDOW)
    evaluate = prepare_hc_hp_def_window(fc, &hc_dat_local);
  else
    evaluate = prepare_hc_hp_def(fc, &hc_dat_local);

  if ((i > 0) && (j > 0)) {
    if (evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, &hc_dat_local)) {
      if (j > i)  /* linear case */
        return exp_eval_hp_loop(fc, i, j);
      else        /* circular case */
        return exp_eval_ext_hp_loop(fc, j, i);
    }
  }

  return 0.;
}







PRIVATE unsigned char
hc_int_cb_def(int   i,
              int   j,
              int   k,
              int   l,
              void  *data)
{
  unsigned char pij, pkl;

  struct hc_int_def_dat *dat = (struct hc_int_def_dat *)data;

  if ((dat->sn[i] != dat->sn[k]) ||
      (dat->sn[l] != dat->sn[j]))
    return (unsigned char)0;

  if (dat->mx) {
    pij = dat->mx[dat->n * i + j];
    pkl = dat->mx[dat->n * k + l];
  } else {
    pij = dat->mx_local[i][j - i];
    pkl = dat->mx_local[k][l - k];
  }

  if ((pij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
      (pkl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
    return (unsigned char)1;

  return (unsigned char)0;
}


PRIVATE unsigned char
hc_int_cb_def_user(int  i,
                   int  j,
                   int  k,
                   int  l,
                   void *data)
{
  struct hc_int_def_dat *dat = (struct hc_int_def_dat *)data;

  unsigned char eval = hc_int_cb_def(i, j, k, l, data);

  return (dat->hc_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, dat->hc_dat)) ? eval : (unsigned char)0;
}


PRIVATE INLINE eval_hc
prepare_hc_int_def(vrna_fold_compound_t   *fc,
                   struct hc_int_def_dat  *dat)
{
  dat->mx       = (fc->hc->type == VRNA_HC_WINDOW) ? NULL : fc->hc->mx;
  dat->mx_local = (fc->hc->type == VRNA_HC_WINDOW) ? fc->hc->matrix_local : NULL;
  dat->n        = fc->length;
  dat->up       = fc->hc->up_int;
  dat->sn       = fc->strand_number;
  dat->hc_f     = NULL;
  dat->hc_dat   = NULL;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return &hc_int_cb_def_user;
  }

  return &hc_int_cb_def;
}

PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up(int                    i,
                 int                    j,
                 int                    k,
                 int                    l,
                 struct sc_int_exp_dat  *data)
{
  int         u1, u2;
  FLT_OR_DBL  sc;

  u1  = k - i - 1;
  u2  = j - l - 1;

  sc = 1.;

  if (u1 > 0)
    sc *= data->up[i + 1][u1];

  if (u2 > 0)
    sc *= data->up[l + 1][u2];

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp(int                    i,
                 int                    j,
                 int                    k,
                 int                    l,
                 struct sc_int_exp_dat  *data)
{
  return data->bp[data->idx[j] + i];
}

PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local(int                    i,
                       int                    j,
                       int                    k,
                       int                    l,
                       struct sc_int_exp_dat  *data)
{
  return data->bp_local[i][j - i];
}

PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_stack(int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    struct sc_int_exp_dat *data)
{
  FLT_OR_DBL sc;

  sc = 1.;

  if ((i + 1 == k) && (l + 1 == j)) {
    sc *= data->stack[i] *
          data->stack[k] *
          data->stack[l] *
          data->stack[j];
  }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_user(int                    i,
                   int                    j,
                   int                    k,
                   int                    l,
                   struct sc_int_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_PAIR_IL,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp(int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local(int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_stack(int                    i,
                       int                    j,
                       int                    k,
                       int                    l,
                       struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_user(int                   i,
                      int                   j,
                      int                   k,
                      int                   l,
                      struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_stack(int                    i,
                       int                    j,
                       int                    k,
                       int                    l,
                       struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_stack(int                    i,
                             int                    j,
                             int                    k,
                             int                    l,
                             struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_user(int                   i,
                      int                   j,
                      int                   k,
                      int                   l,
                      struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_user(int                   i,
                            int                   j,
                            int                   k,
                            int                   l,
                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_stack_user(int                    i,
                         int                    j,
                         int                    k,
                         int                    l,
                         struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_stack(int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_stack(int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_user(int                    i,
                         int                    j,
                         int                    k,
                         int                    l,
                         struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_user(int                    i,
                               int                    j,
                               int                    k,
                               int                    l,
                               struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_stack_user(int                   i,
                            int                   j,
                            int                   k,
                            int                   l,
                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_stack_user(int                   i,
                            int                   j,
                            int                   k,
                            int                   l,
                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_stack_user(int                   i,
                                  int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_stack_user(int                    i,
                               int                    j,
                               int                    k,
                               int                    l,
                               struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_stack_user(int                    i,
                                     int                    j,
                                     int                    k,
                                     int                    l,
                                     struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up(int                    i,
                     int                    j,
                     int                    k,
                     int                    l,
                     struct sc_int_exp_dat  *data)
{
  int         u1, u2, u3;
  FLT_OR_DBL  sc;

  sc  = 1.;
  u1  = i - 1;
  u2  = k - j - 1;
  u3  = data->n - l;

  if (u1 > 0)
    sc *= data->up[1][u1];

  if (u2 > 0)
    sc *= data->up[j + 1][u2];

  if (u3 > 0)
    sc *= data->up[l + 1][u3];

  return sc;
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_stack(int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        struct sc_int_exp_dat *data)
{
  FLT_OR_DBL sc;

  sc = 1.;

  if ((i == 1) && (j + 1 == k) && (l == data->n)) {
    sc *= data->stack[i] *
          data->stack[k] *
          data->stack[l] *
          data->stack[j];
  }

  return sc;
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_user(int                    i,
                       int                    j,
                       int                    k,
                       int                    l,
                       struct sc_int_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_PAIR_IL,
                       data->user_data);
}




PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_stack(int                    i,
                           int                    j,
                           int                    k,
                           int                    l,
                           struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_ext_up(i, j, k, l, data) *
         sc_int_exp_cb_ext_stack(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_user(int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_ext_up(i, j, k, l, data) *
         sc_int_exp_cb_ext_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_stack_user(int                    i,
                             int                    j,
                             int                    k,
                             int                    l,
                             struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_ext_stack(i, j, k, l, data) *
         sc_int_exp_cb_ext_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_stack_user(int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_ext_up(i, j, k, l, data) *
         sc_int_exp_cb_ext_stack(i, j, k, l, data) *
         sc_int_exp_cb_ext_user(i, j, k, l, data);
}





PRIVATE INLINE void
init_sc_int_exp(vrna_fold_compound_t  *fc,
                struct sc_int_exp_dat *sc_wrapper)
{
  unsigned char sliding_window;
  unsigned int  s, provides_sc_up, provides_sc_bp, provides_sc_stack, provides_sc_user;
  vrna_sc_t     *sc, **scs;

  if (fc->exp_matrices)
    sliding_window = (fc->exp_matrices->type == VRNA_MX_WINDOW) ? 1 : 0;
  else if ((fc->type == VRNA_FC_TYPE_SINGLE) && (fc->sc))
    sliding_window = (fc->sc->type == VRNA_SC_WINDOW) ? 1 : 0;
  else if (fc->hc)
    sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  else
    sliding_window = 0;

  provides_sc_up    = 0;
  provides_sc_bp    = 0;
  provides_sc_stack = 0;
  provides_sc_user  = 0;

  sc_wrapper->n                     = fc->length;
  sc_wrapper->n_seq                 = 1;
  sc_wrapper->a2s                   = NULL;
  sc_wrapper->idx                   = fc->jindx;
  sc_wrapper->up                    = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->bp                    = NULL;
  sc_wrapper->bp_comparative        = NULL;
  sc_wrapper->bp_local              = NULL;
  sc_wrapper->bp_local_comparative  = NULL;
  sc_wrapper->stack                 = NULL;
  sc_wrapper->stack_comparative     = NULL;
  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->user_data_comparative = NULL;

  sc_wrapper->pair      = NULL;
  sc_wrapper->pair_ext  = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        sc_wrapper->up        = sc->exp_energy_up;
        sc_wrapper->bp        = (sliding_window) ? NULL : sc->exp_energy_bp;
        sc_wrapper->bp_local  = (sliding_window) ? sc->exp_energy_bp_local : NULL;
        sc_wrapper->stack     = sc->exp_energy_stack;
        sc_wrapper->user_cb   = sc->exp_f;
        sc_wrapper->user_data = sc->data;

        if (sc->exp_energy_up)
          provides_sc_up = 1;

        if (sliding_window) {
          if (sc->exp_energy_bp_local)
            provides_sc_bp = 1;
        } else if (sc->exp_energy_bp) {
          provides_sc_bp = 1;
        }

        if (sc->exp_energy_stack)
          provides_sc_stack = 1;

        if (sc->exp_f)
          provides_sc_user = 1;

        if (provides_sc_user) {
          if (provides_sc_up) {
            if (provides_sc_bp) {
              if (provides_sc_stack) {
                sc_wrapper->pair = (sliding_window) ?
                                   &sc_int_exp_cb_up_bp_local_stack_user :
                                   &sc_int_exp_cb_up_bp_stack_user;
                sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_stack_user;
              } else {
                sc_wrapper->pair = (sliding_window) ?
                                   &sc_int_exp_cb_up_bp_local_user :
                                   &sc_int_exp_cb_up_bp_user;
                sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_user;
              }
            } else if (provides_sc_stack) {
              sc_wrapper->pair      = &sc_int_exp_cb_up_stack_user;
              sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_stack_user;
            } else {
              sc_wrapper->pair      = &sc_int_exp_cb_up_user;
              sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_user;
            }
          } else if (provides_sc_bp) {
            if (provides_sc_stack) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_bp_local_stack_user :
                                 &sc_int_exp_cb_bp_stack_user;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_stack_user;
            } else {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_bp_local_user :
                                 &sc_int_exp_cb_bp_user;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_user;
            }
          } else if (provides_sc_stack) {
            sc_wrapper->pair      = &sc_int_exp_cb_stack_user;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_stack_user;
          } else {
            sc_wrapper->pair      = &sc_int_exp_cb_user;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_user;
          }
        } else if (provides_sc_bp) {
          if (provides_sc_up) {
            if (provides_sc_stack) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_up_bp_local_stack :
                                 &sc_int_exp_cb_up_bp_stack;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_stack;
            } else {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_up_bp_local :
                                 &sc_int_exp_cb_up_bp;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up;
            }
          } else if (provides_sc_stack) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_int_exp_cb_bp_local_stack :
                               &sc_int_exp_cb_bp_stack;
            sc_wrapper->pair_ext = &sc_int_exp_cb_ext_stack;
          } else {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_int_exp_cb_bp_local :
                               &sc_int_exp_cb_bp;
          }
        } else if (provides_sc_up) {
          if (provides_sc_stack) {
            sc_wrapper->pair      = &sc_int_exp_cb_up_stack;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_stack;
          } else {
            sc_wrapper->pair      = &sc_int_exp_cb_up;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up;
          }
        } else if (provides_sc_stack) {
          sc_wrapper->pair      = &sc_int_exp_cb_stack;
          sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_stack;
        }
      }

      break;
  }
}




PRIVATE INLINE FLT_OR_DBL
exp_E_IntLoop(int               u1,
              int               u2,
              int               type,
              int               type2,
              short             si1,
              short             sj1,
              short             sp1,
              short             sq1,
              vrna_exp_param_t  *P)
{
  int     ul, us, no_close = 0;
  double  z           = 0.;
  int     noGUclosure = P->model_details.noGUclosure;
  int     backbones;
  double  salt_stack_correction = P->expSaltStack;
  double  salt_loop_correction = 1.;

  if ((noGUclosure) && ((type2 == 3) || (type2 == 4) || (type == 3) || (type == 4)))
    no_close = 1;

  if (u1 > u2) {
    ul  = u1;
    us  = u2;
  } else {
    ul  = u2;
    us  = u1;
  }

  /* salt correction for loop */
  backbones = ul+us+2;

  if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
    if (backbones <= MAXLOOP+1)
      salt_loop_correction = P->expSaltLoop[backbones];
    else
      salt_loop_correction = exp(-vrna_salt_loop_int(backbones, P->model_details.salt, P->temperature+K0, P->model_details.backbone_length) * 10. / P->kT);
  }

  if (ul == 0) {
    /* stack */
    z = P->expstack[type][type2] * salt_stack_correction;
  } else if (!no_close) {
    if (us == 0) {
      /* bulge */
      z = P->expbulge[ul];
      if (ul == 1) {
        z *= P->expstack[type][type2];
      } else {
        if (type > 2)
          z *= P->expTermAU;

        if (type2 > 2)
          z *= P->expTermAU;
      }

      return (FLT_OR_DBL)(z * salt_loop_correction);
    } else if (us == 1) {
      if (ul == 1)                     /* 1x1 loop */
        return (FLT_OR_DBL)(P->expint11[type][type2][si1][sj1] * salt_loop_correction);

      if (ul == 2) {
        /* 2x1 loop */
        if (u1 == 1)
          return (FLT_OR_DBL)(P->expint21[type][type2][si1][sq1][sj1] * salt_loop_correction);
        else
          return (FLT_OR_DBL)(P->expint21[type2][type][sq1][si1][sp1] * salt_loop_correction);
      } else {
        /* 1xn loop */
        z = P->expinternal[ul + us] * P->expmismatch1nI[type][si1][sj1] *
            P->expmismatch1nI[type2][sq1][sp1];
        return (FLT_OR_DBL)(z * P->expninio[2][ul - us] * salt_loop_correction);
      }
    } else if (us == 2) {
      if (ul == 2) {
        /* 2x2 loop */
        return (FLT_OR_DBL)(P->expint22[type][type2][si1][sp1][sq1][sj1] * salt_loop_correction);
      } else if (ul == 3) {
        /* 2x3 loop */
        z = P->expinternal[5] * P->expmismatch23I[type][si1][sj1] *
            P->expmismatch23I[type2][sq1][sp1];
        return (FLT_OR_DBL)(z * P->expninio[2][1] * salt_loop_correction);
      }
    }

    /* generic interior loop (no else here!)*/
    z = P->expinternal[ul + us] * P->expmismatchI[type][si1][sj1] *
        P->expmismatchI[type2][sq1][sp1];
    return (FLT_OR_DBL)(z * P->expninio[2][ul - us] * salt_loop_correction);
  }

  return (FLT_OR_DBL)z;
}


PRIVATE INLINE void
free_sc_int_exp(struct sc_int_exp_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->bp_comparative);
  free(sc_wrapper->bp_local_comparative);
  free(sc_wrapper->stack_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}

#define VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP 4U

PRIVATE FLT_OR_DBL
exp_E_ext_int_loop(vrna_fold_compound_t *fc,
                   int                  i,
                   int                  j)
{
  unsigned char         *hc_mx, eval_loop;
  short                 *S, *S2, **SS, **S5, **S3;
  unsigned int          *tt, n_seq, s, **a2s, type, type2;
  int                   k, l, u1, u2, u3, qmin, with_ud,
                        n, *my_iindx, *hc_up,
                        u1_local, u2_local, u3_local;
  FLT_OR_DBL            q, q_temp, *qb, *scale;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  eval_hc               evaluate;
  struct hc_int_def_dat hc_dat_local;
  struct sc_int_exp_dat sc_wrapper;

  n           = fc->length;
  n_seq       = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  S           = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  S2          = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  my_iindx    = fc->iindx;
  qb          = fc->exp_matrices->qb;
  scale       = fc->exp_matrices->scale;
  hc_mx       = fc->hc->mx;
  hc_up       = fc->hc->up_int;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  type        = 0;
  tt          = NULL;
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->exp_energy_cb)) ? 1 : 0;

  q = 0.;

  evaluate = prepare_hc_int_def(fc, &hc_dat_local);

  init_sc_int_exp(fc, &sc_wrapper);

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    /* prepare necessary variables */
    if (fc->type == VRNA_FC_TYPE_SINGLE) {
      type = vrna_get_ptype_md(S2[j], S2[i], md);
    } else {
      tt = (unsigned int *)vrna_alloc(sizeof(unsigned int) * n_seq);

      for (s = 0; s < n_seq; s++)
        tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
    }

    for (k = j + 1; k < n; k++) {
      u2 = k - j - 1;
      if (u2 + i - 1 > MAXLOOP)
        break;

      if (hc_up[j + 1] < u2)
        break;

      qmin = u2 + i - 1 + n - MAXLOOP;
      if (qmin < k + 1)
        qmin = k + 1;

      for (l = n; l >= qmin; l--) {
        u1  = i - 1;
        u3  = n - l;
        if (hc_up[l + 1] < (u1 + u3))
          break;

        if (u1 + u2 + u3 > MAXLOOP)
          continue;

        eval_loop = hc_mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP;

        if (eval_loop && evaluate(i, j, k, l, &hc_dat_local)) {
          q_temp = qb[my_iindx[k] - l];

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              type2 = vrna_get_ptype_md(S2[l], S2[k], md);

              /* regular interior loop */
              q_temp *=
                exp_E_IntLoop(u2,
                              u1 + u3,
                              type,
                              type2,
                              S[j + 1],
                              S[i - 1],
                              S[k - 1],
                              S[l + 1],
                              pf_params);
              break;

          }

          if (sc_wrapper.pair_ext)
            q_temp *= sc_wrapper.pair_ext(i, j, k, l, &sc_wrapper);

          q += q_temp *
               scale[u1 + u2 + u3];

          if (with_ud) {
            FLT_OR_DBL q5, q3;

            q5  = q3 = 0.;
            u1  = i - 1;
            u2  = k - j - 1;
            u3  = n - l;

            if (u2 > 0) {
              q5 = domains_up->exp_energy_cb(fc,
                                             j + 1, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
            }

            if (u1 + u3 > 0) {
              q3 = domains_up->exp_energy_cb(fc,
                                             l + 1, i - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
            }

            q += q_temp *
                 q5 *
                 scale[u1 + u2 + u3];
            q += q_temp *
                 q3 *
                 scale[u1 + u2 + u3];
            q += q_temp *
                 q5 *
                 q3 *
                 scale[u1 + u2 + u3];
          }
        }
      }
    }
  }

  free(tt);
  free_sc_int_exp(&sc_wrapper);

  return q;
}


PUBLIC unsigned int
vrna_get_ptype(int  ij,
               char *ptype)
{
  unsigned int tt = (unsigned int)ptype[ij];

  return (tt == 0) ? 7 : tt;
}

PUBLIC unsigned int
vrna_get_ptype_window(int   i,
                      int   j,
                      char  **ptype)
{
  unsigned int tt = (unsigned int)ptype[i][j - i];

  return (tt == 0) ? 7 : tt;
}


PRIVATE FLT_OR_DBL
exp_E_int_loop(vrna_fold_compound_t *fc,
               int                  i,
               int                  j)
{
  unsigned char         sliding_window, hc_decompose_ij, hc_decompose_kl;
  char                  *ptype, **ptype_local;
  unsigned char         *hc_mx, **hc_mx_local;
  short                 *S1, **SS, **S5, **S3;
  unsigned int          *sn, *se, *ss, n_seq, s, **a2s, n;
  int                   *rtype, noclose, *my_iindx, *jindx, *hc_up, ij,
                        with_gquad, with_ud;
  FLT_OR_DBL            qbt1, q_temp, *qb, **qb_local, *G, *scale;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  eval_hc               evaluate;
  struct hc_int_def_dat hc_dat_local;
  struct sc_int_exp_dat sc_wrapper;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n               = fc->length;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  sn              = fc->strand_number;
  se              = fc->strand_end;
  ss              = fc->strand_start;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S1          = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  qb          = (sliding_window) ? NULL : fc->exp_matrices->qb;
  G           = (sliding_window) ? NULL : fc->exp_matrices->G;
  qb_local    = (sliding_window) ? fc->exp_matrices->qb_local : NULL;
  scale       = fc->exp_matrices->scale;
  my_iindx    = fc->iindx;
  jindx       = fc->jindx;
  hc_mx       = (sliding_window) ? NULL : fc->hc->mx;
  hc_mx_local = (sliding_window) ? fc->hc->matrix_local : NULL;
  hc_up       = fc->hc->up_int;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  domains_up  = fc->domains_up;
  with_ud     = ((domains_up) && (domains_up->exp_energy_cb)) ? 1 : 0;
  rtype       = &(md->rtype[0]);
  qbt1        = 0.;
  evaluate    = prepare_hc_int_def(fc, &hc_dat_local);

  init_sc_int_exp(fc, &sc_wrapper);

  ij = (sliding_window) ? 0 : jindx[j] + i;

  hc_decompose_ij = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_decompose_ij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    unsigned int  type, type2, *tt;
    int           k, l, kl, last_k, first_l, u1, u2, noGUclosure;

    noGUclosure = md->noGUclosure;
    tt          = NULL;
    type        = 0;

    if (fc->type == VRNA_FC_TYPE_SINGLE)
      type = sliding_window ?
             vrna_get_ptype_window(i, j + i, ptype_local) :
             vrna_get_ptype(ij, ptype);

    noclose = ((noGUclosure) && (type == 3 || type == 4)) ? 1 : 0;


    /* handle stacks separately */
    k = i + 1;
    l = j - 1;
    if ((k < l) && (sn[i] == sn[k]) && (sn[l] == sn[j])) {
      kl              = (sliding_window) ? 0 : jindx[l] + k;
      hc_decompose_kl = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[n * k + l];

      if ((hc_decompose_kl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
          (evaluate(i, j, k, l, &hc_dat_local))) {
        q_temp = (sliding_window) ? qb_local[k][l] : qb[my_iindx[k] - l];

        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            type2 = sliding_window ?
                    rtype[vrna_get_ptype_window(k, l + k, ptype_local)] :
                    rtype[vrna_get_ptype(kl, ptype)];

            q_temp *= exp_E_IntLoop(0,
                                    0,
                                    type,
                                    type2,
                                    S1[i + 1],
                                    S1[j - 1],
                                    S1[k - 1],
                                    S1[l + 1],
                                    pf_params);

            break;

          
        }

        if (sc_wrapper.pair)
          q_temp *= sc_wrapper.pair(i, j, k, l, &sc_wrapper);

        qbt1 += q_temp *
                scale[2];
      }
    }

    if (!noclose) {
      /* only proceed if the enclosing pair is allowed */

      /* handle bulges in 5' side */
      l = j - 1;
      if ((l > i + 2) && (sn[j] == sn[l])) {
        last_k = l - 1;

        if (last_k > i + 1 + MAXLOOP)
          last_k = i + 1 + MAXLOOP;

        if (last_k > i + 1 + hc_up[i + 1])
          last_k = i + 1 + hc_up[i + 1];

        if (last_k > se[sn[i]])
          last_k = se[sn[i]];

        u1 = 1;

        k     = i + 2;
        kl    = (sliding_window) ? 0 : jindx[l] + k;
        hc_mx += n * l;

        for (; k <= last_k; k++, u1++, kl++) {
          hc_decompose_kl = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[k];

          if ((hc_decompose_kl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (evaluate(i, j, k, l, &hc_dat_local))) {
            q_temp = (sliding_window) ? qb_local[k][l] : qb[my_iindx[k] - l];

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type2 = sliding_window ?
                        rtype[vrna_get_ptype_window(k, l + k, ptype_local)] :
                        rtype[vrna_get_ptype(kl, ptype)];

                if ((noGUclosure) && (type2 == 3 || type2 == 4))
                  continue;

                q_temp *= exp_E_IntLoop(u1,
                                        0,
                                        type,
                                        type2,
                                        S1[i + 1],
                                        S1[j - 1],
                                        S1[k - 1],
                                        S1[l + 1],
                                        pf_params);

                break;

              
            }

            if (sc_wrapper.pair)
              q_temp *= sc_wrapper.pair(i, j, k, l, &sc_wrapper);

            qbt1 += q_temp *
                    scale[u1 + 2];

            if (with_ud) {
              q_temp *= domains_up->exp_energy_cb(fc,
                                                  i + 1, k - 1,
                                                  VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                  domains_up->data);
              qbt1 += q_temp *
                      scale[u1 + 2];
            }
          }
        }

        hc_mx -= n * l;
      }

      /* handle bulges in 3' side */
      k = i + 1;
      if ((k < j - 2) && (sn[i] == sn[k])) {
        first_l = k + 1;
        if (first_l < j - 1 - MAXLOOP)
          first_l = j - 1 - MAXLOOP;

        if (first_l < ss[sn[j]])
          first_l = ss[sn[j]];

        u2    = 1;
        hc_mx += n * k;

        for (l = j - 2; l >= first_l; l--, u2++) {
          if (u2 > hc_up[l + 1])
            break;

          kl              = (sliding_window) ? 0 : jindx[l] + k;
          hc_decompose_kl = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[l];

          if ((hc_decompose_kl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (evaluate(i, j, k, l, &hc_dat_local))) {
            q_temp = (sliding_window) ? qb_local[k][l] : qb[my_iindx[k] - l];

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type2 = sliding_window ?
                        rtype[vrna_get_ptype_window(k, l + k, ptype_local)] :
                        rtype[vrna_get_ptype(kl, ptype)];

                if ((noGUclosure) && (type2 == 3 || type2 == 4))
                  continue;

                q_temp *= exp_E_IntLoop(0,
                                        u2,
                                        type,
                                        type2,
                                        S1[i + 1],
                                        S1[j - 1],
                                        S1[k - 1],
                                        S1[l + 1],
                                        pf_params);

                break;

              
            }

            if (sc_wrapper.pair)
              q_temp *= sc_wrapper.pair(i, j, k, l, &sc_wrapper);

            qbt1 += q_temp *
                    scale[u2 + 2];

            if (with_ud) {
              q_temp *= domains_up->exp_energy_cb(fc,
                                                  l + 1, j - 1,
                                                  VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                  domains_up->data);
              qbt1 += q_temp *
                      scale[u2 + 2];
            }
          }
        }

        hc_mx -= n * k;
      }

      /* last but not least, all other internal loops */
      last_k = j - 3;

      if (last_k > i + MAXLOOP + 1)
        last_k = i + MAXLOOP + 1;

      if (last_k > i + 1 + hc_up[i + 1])
        last_k = i + 1 + hc_up[i + 1];

      if (last_k > se[sn[i]])
        last_k = se[sn[i]];

      u1 = 1;

      for (k = i + 2; k <= last_k; k++, u1++) {
        first_l = k + 1;

        if (first_l < j - 1 - MAXLOOP + u1)
          first_l = j - 1 - MAXLOOP + u1;

        if (first_l < ss[sn[j]])
          first_l = ss[sn[j]];

        u2 = 1;

        hc_mx += n * k;

        for (l = j - 2; l >= first_l; l--, u2++) {
          if (hc_up[l + 1] < u2)
            break;

          kl              = (sliding_window) ? 0 : jindx[l] + k;
          hc_decompose_kl = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[l];

          if ((hc_decompose_kl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (evaluate(i, j, k, l, &hc_dat_local))) {
            q_temp = (sliding_window) ? qb_local[k][l] : qb[my_iindx[k] - l];

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type2 = sliding_window ?
                        rtype[vrna_get_ptype_window(k, l + k, ptype_local)] :
                        rtype[vrna_get_ptype(kl, ptype)];

                if ((noGUclosure) && (type2 == 3 || type2 == 4))
                  continue;

                q_temp *= exp_E_IntLoop(u1,
                                        u2,
                                        type,
                                        type2,
                                        S1[i + 1],
                                        S1[j - 1],
                                        S1[k - 1],
                                        S1[l + 1],
                                        pf_params);

                break;

            }

            if (sc_wrapper.pair)
              q_temp *= sc_wrapper.pair(i, j, k, l, &sc_wrapper);

            qbt1 += q_temp *
                    scale[u1 + u2 + 2];

            if (with_ud) {
              FLT_OR_DBL q5, q3;

              q5 = domains_up->exp_energy_cb(fc,
                                             i + 1, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);
              q3 = domains_up->exp_energy_cb(fc,
                                             l + 1, j - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                             domains_up->data);

              qbt1 += q_temp *
                      q5 *
                      scale[u1 + u2 + 2];
              qbt1 += q_temp *
                      q3 *
                      scale[u1 + u2 + 2];
              qbt1 += q_temp *
                      q5 *
                      q3 *
                      scale[u1 + u2 + 2];
            }
          }
        }

        hc_mx -= n * k;
      }

    }

    free(tt);
  }

  free_sc_int_exp(&sc_wrapper);

  return qbt1;
}



PUBLIC FLT_OR_DBL
vrna_exp_E_int_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j)
{
  FLT_OR_DBL q = 0.;

  if ((fc) && (i > 0) && (j > 0)) {
    if (j < i) {
      /* Note: j < i indicates that we want to evaluate exterior int loop (for circular RNAs)! */
      if (fc->hc->type == VRNA_HC_WINDOW) {
        printf(
          "vrna_exp_E_int_loop: invalid sequence positions for pair (i,j) = (%d,%d)!",
          i,
          j);
      } else {
        q = exp_E_ext_int_loop(fc, j, i);
      }
    } else {
      q = exp_E_int_loop(fc, i, j);
    }
  }

  return q;
}





PRIVATE unsigned char
hc_mb_cb_def(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             void           *data)
{
  unsigned char         eval;
  unsigned int          n;
  int                   di, dj, u;
  struct hc_mb_def_dat  *dat = (struct hc_mb_def_dat *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  n     = dat->n;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) &&
          (dat->hc_up[k + 1] < u))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;

      if ((di != 0) &&
          (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) &&
          (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) &&
            (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) &&
            (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      if (dat->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) &&
            (dat->hc_up[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) &&
            (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML_EXT:
      if (dat->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di++;
        dj++;
        if ((di != 0) &&
            (dat->hc_up[k + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) &&
            (dat->hc_up[j + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_ML_STEM:
      u     = l - k - 1;
      if (dat->mx[n * j + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      if ((u != 0) && (dat->hc_up[k + 1] < u))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if ((dat->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (dat->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      printf("hc_mb_cb_def@multibranch_hc.inc: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}




PRIVATE unsigned char
hc_sn(int           i,
      int           j,
      int           k,
      int           l,
      unsigned char d,
      void          *data)
{
  unsigned int          *sn;
  unsigned char         eval;
  struct hc_mb_def_dat  *dat = (struct hc_mb_def_dat *)data;

  sn    = dat->sn;
  eval  = (unsigned char)0;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      /* fall through */
    case VRNA_DECOMP_ML_ML_STEM:
      if (sn[k] == sn[l])
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_STEM:
    /* fall through */

    case VRNA_DECOMP_ML_ML:
      if ((sn[i] == sn[k]) &&
          (sn[l] == sn[j]) &&
          (sn[i - 1] == sn[i]) &&
          (sn[j + 1] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_PAIR_ML_EXT:
      /* fall through */
    case VRNA_DECOMP_PAIR_ML:
      if ((sn[i] == sn[k]) &&
          (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if ((i == k - 1) &&
          (sn[i] == sn[k]))
        eval = (unsigned char)1;
      else if ((l + 1 == j) &&
               (sn[l] == sn[j]))
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if (sn[j] == sn[k])
        eval = (unsigned char)1;

      break;

    default:
      printf("hc_sn@multibranch_hc.inc: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE unsigned char
hc_mb_cb_def_window(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    void          *data)
{
  int                   di, dj, u;
  unsigned char         eval;
  struct hc_mb_def_dat  *dat = (struct hc_mb_def_dat *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) && (dat->hc_up[k + 1] < u))
        eval = (unsigned char)0;

      if (dat->sn[k] != dat->sn[l])
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;
      if ((di != 0) && ((dat->hc_up[i] < di) || (dat->sn[i] != dat->sn[k])))
        eval = (unsigned char)0;

      if ((dj != 0) && ((dat->hc_up[l + 1] < dj) || (dat->sn[l] != dat->sn[j])))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      if (dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) && (dat->hc_up[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if ((dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      printf("hc_mb_cb_def_window@multibranch_hc.inc: "
                           "Unrecognized decomposition %d",
                           d);
  }

  return eval;
}


PRIVATE INLINE unsigned char
hc_mb_cb_def_sn(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data)
{
  unsigned char eval;

  eval  = hc_mb_cb_def(i, j, k, l, d, data);
  eval  = hc_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_mb_cb_def_user(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data)
{
  unsigned char         eval;
  struct hc_mb_def_dat  *dat = (struct hc_mb_def_dat *)data;

  eval  = hc_mb_cb_def(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_mb_cb_def_sn_user(int            i,
                     int            j,
                     int            k,
                     int            l,
                     unsigned char  d,
                     void           *data)
{
  unsigned char         eval;
  struct hc_mb_def_dat  *dat = (struct hc_mb_def_dat *)data;

  eval  = hc_mb_cb_def(i, j, k, l, d, data);
  eval  = hc_sn(i, j, k, l, d, data) ? eval : (unsigned char)0;
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_mb_cb_def_user_window(int            i,
                         int            j,
                         int            k,
                         int            l,
                         unsigned char  d,
                         void           *data)
{
  unsigned char         eval;
  struct hc_mb_def_dat  *dat = (struct hc_mb_def_dat *)data;

  eval  = hc_mb_cb_def_window(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE INLINE vrna_hc_eval_f
prepare_hc_mb_def(vrna_fold_compound_t  *fc,
                  struct hc_mb_def_dat  *dat)
{
  dat->mx         = fc->hc->mx;
  dat->n          = fc->hc->n;
  dat->mx_window  = fc->hc->matrix_local;
  dat->hc_up      = fc->hc->up_ml;
  dat->sn         = fc->strand_number;

  if (fc->hc->f) {
    dat->hc_f   = fc->hc->f;
    dat->hc_dat = fc->hc->data;
    return (fc->hc->type == VRNA_HC_WINDOW) ?
           &hc_mb_cb_def_user_window :
           ((fc->strands == 1) ?
            &hc_mb_cb_def_user :
            &hc_mb_cb_def_sn_user);
  }

  return (fc->hc->type == VRNA_HC_WINDOW) ?
         &hc_mb_cb_def_window :
         ((fc->strands == 1) ?
          &hc_mb_cb_def :
          hc_mb_cb_def_sn);
}




PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp(int                  i,
                     int                  j,
                     struct sc_mb_exp_dat *data)
{
  return data->bp[data->idx[j] + i];
}



PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_local(int                  i,
                           int                  j,
                           struct sc_mb_exp_dat *data)
{
  return data->bp_local[i][j - i];
}




PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_user(int                  i,
                       int                  j,
                       struct sc_mb_exp_dat *data)
{
  return data->user_cb(i, j, i + 1, j - 1,
                       VRNA_DECOMP_PAIR_ML,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_ext_cb_user(int                  i,
                           int                  j,
                           struct sc_mb_exp_dat *data)
{
  return data->user_cb(i, j, i - 1, j + 1,
                       VRNA_DECOMP_PAIR_ML,
                       data->user_data);
}





PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_user(int                   i,
                          int                   j,
                          struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_pair_cb_bp(i, j, data) *
         sc_mb_exp_pair_cb_user(i, j, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_local_user(int                   i,
                                int                   j,
                                struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_pair_cb_bp_local(i, j, data) *
         sc_mb_exp_pair_cb_user(i, j, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_up(int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    struct sc_mb_exp_dat  *data)
{
  int         l1  = k - i;
  int         l2  = j - l;
  FLT_OR_DBL  sc  = 1.;

  if (l1 > 0)
    sc *= data->up[i][l1];

  if (l2 > 0)
    sc *= data->up[l + 1][l2];

  return sc;
}




PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_user(int                   i,
                      int                   j,
                      int                   k,
                      int                   l,
                      struct sc_mb_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_ML,
                       data->user_data);
}




PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_up_user(int                  i,
                         int                  j,
                         int                  k,
                         int                  l,
                         struct sc_mb_exp_dat *data)
{
  return sc_mb_exp_red_cb_up(i, j, k, l, data) *
         sc_mb_exp_red_cb_user(i, j, k, l, data);
}




PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_stem_user(int                  i,
                           int                  j,
                           int                  k,
                           int                  l,
                           struct sc_mb_exp_dat *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_STEM,
                       data->user_data);
}




PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_stem_up_user(int                   i,
                              int                   j,
                              int                   k,
                              int                   l,
                              struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_red_cb_up(i, j, k, l, data) *
         sc_mb_exp_red_cb_stem_user(i, j, k, l, data);
}



PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_split_cb_user(int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        struct sc_mb_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_ML_ML,
                       data->user_data);
}



PRIVATE INLINE void
init_sc_mb_exp(vrna_fold_compound_t *fc,
               struct sc_mb_exp_dat *sc_wrapper)
{
  unsigned char sliding_window;
  vrna_sc_t     *sc, **scs;

  sc_wrapper->n     = fc->length;
  sc_wrapper->n_seq = 1;
  sc_wrapper->idx   = fc->jindx;
  sc_wrapper->a2s   = NULL;

  sc_wrapper->up                    = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->bp                    = NULL;
  sc_wrapper->bp_comparative        = NULL;
  sc_wrapper->bp_local              = NULL;
  sc_wrapper->bp_local_comparative  = NULL;

  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  sc_wrapper->pair      = NULL;
  sc_wrapper->pair_ext  = NULL;
  sc_wrapper->red_stem  = NULL;
  sc_wrapper->red_ml    = NULL;
  sc_wrapper->decomp_ml = NULL;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        unsigned int provides_sc_up, provides_sc_bp, provides_sc_user;

        provides_sc_up    = 0;
        provides_sc_bp    = 0;
        provides_sc_user  = 0;

        sc_wrapper->up        = sc->exp_energy_up;
        sc_wrapper->user_cb   = sc->exp_f;
        sc_wrapper->user_data = sc->data;

        if (sliding_window)
          sc_wrapper->bp_local = sc->exp_energy_bp_local;
        else
          sc_wrapper->bp = sc->exp_energy_bp;

        if (sc->exp_energy_up)
          provides_sc_up = 1;

        if (sliding_window) {
          if (sc->exp_energy_bp_local)
            provides_sc_bp = 1;
        } else if (sc->exp_energy_bp) {
          provides_sc_bp = 1;
        }

        if (sc->exp_f)
          provides_sc_user = 1;

        /* done initializing, now assign function pointers */
        if (provides_sc_user) {
          sc_wrapper->decomp_ml = &sc_mb_exp_split_cb_user;
          sc_wrapper->red_stem  = &sc_mb_exp_red_cb_stem_user;
          sc_wrapper->red_ml    = &sc_mb_exp_red_cb_user;
          sc_wrapper->pair      = &sc_mb_exp_pair_cb_user;
          if (!sliding_window)
            sc_wrapper->pair_ext  = &sc_mb_exp_pair_ext_cb_user;

          if (provides_sc_bp) {
            if (sliding_window) {
              sc_wrapper->pair = &sc_mb_exp_pair_cb_bp_local_user;
            } else {
              sc_wrapper->pair      = &sc_mb_exp_pair_cb_bp_user;
              sc_wrapper->pair_ext  = &sc_mb_exp_pair_ext_cb_user;
            }
          }

          if (provides_sc_up) {
            sc_wrapper->red_stem  = &sc_mb_exp_red_cb_stem_up_user;
            sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up_user;
          }
        } else if (provides_sc_bp) {
          if (sliding_window) {
            sc_wrapper->pair = &sc_mb_exp_pair_cb_bp_local;
          } else {
            sc_wrapper->pair      = &sc_mb_exp_pair_cb_bp;
          }

          if (provides_sc_up) {
            sc_wrapper->red_stem  = &sc_mb_exp_red_cb_up;
            sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up;
          }
        } else if (provides_sc_up) {
          sc_wrapper->red_stem  = &sc_mb_exp_red_cb_up;
          sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up;
        }
      }

      break;
  }
}



PRIVATE INLINE FLT_OR_DBL
exp_E_MLstem(int              type,
             int              si1,
             int              sj1,
             vrna_exp_param_t *P)
{
  double energy = 1.0;

  if (si1 >= 0 && sj1 >= 0)
    energy = P->expmismatchM[type][si1][sj1];
  else if (si1 >= 0)
    energy = P->expdangle5[type][si1];
  else if (sj1 >= 0)
    energy = P->expdangle3[type][sj1];

  if (type > 2)
    energy *= P->expTermAU;

  energy *= P->expMLintern[type];
  return (FLT_OR_DBL)energy;
}



PRIVATE INLINE void
free_sc_mb_exp(struct sc_mb_exp_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->bp_comparative);
  free(sc_wrapper->bp_local_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}



PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast(vrna_fold_compound_t       *fc,
                   int                        i,
                   int                        j,
                   struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  unsigned char             sliding_window;
  char                      *ptype, **ptype_local;
  short                     *S1, **SS, **S5, **S3;
  unsigned int              *sn, n_seq, s, *se;
  int                       ij, k, kl, *my_iindx, *jindx, *rtype, tt;
  FLT_OR_DBL                qbt1, temp, qqqmmm, *qm, **qm_local, *scale, expMLclosing, *qqm1;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_hc_eval_f evaluate;
  struct hc_mb_def_dat      hc_dat_local;
  struct sc_mb_exp_dat      sc_wrapper;

  qqm1            = aux_mx->qqm1;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  se              = fc->strand_end;
  my_iindx        = (sliding_window) ? NULL : fc->iindx;
  jindx           = (sliding_window) ? NULL : fc->jindx;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     = (sliding_window) ? fc->ptype_local : NULL;
  S1              = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  qm              = (sliding_window) ? NULL : fc->exp_matrices->qm;
  qm_local        = (sliding_window) ? fc->exp_matrices->qm_local : NULL;
  scale           = fc->exp_matrices->scale;
  pf_params       = fc->exp_params;
  md              = &(pf_params->model_details);
  ij              = (sliding_window) ? 0 : jindx[j] + i;
  sn              = fc->strand_number;
  hc              = fc->hc;
  expMLclosing    = pf_params->expMLclosing;
  qbt1            = 0.;
  rtype           = &(md->rtype[0]);
  evaluate        = prepare_hc_mb_def(fc, &hc_dat_local);

  init_sc_mb_exp(fc, &sc_wrapper);

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    qqqmmm = pow(expMLclosing, (double)n_seq) *
             scale[2];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        tt = (sliding_window) ?
             rtype[vrna_get_ptype_window(i, j + i, ptype_local)] :
             rtype[vrna_get_ptype(ij, ptype)];
        qqqmmm *= exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params);

        break;
    }

    if (sc_wrapper.pair)
      qqqmmm *= sc_wrapper.pair(i, j, &sc_wrapper);

    FLT_OR_DBL *qqm1_tmp = qqm1;

    if (hc->f) {
      qqm1_tmp  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm1_tmp  -= i;

      for (k = i + 2; k <= j - 1; k++) {
        qqm1_tmp[k] = qqm1[k];
        if (!evaluate(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, &hc_dat_local))
          qqm1_tmp[k] = 0.;
      }
    }

    if (sc_wrapper.decomp_ml) {
      if (qqm1_tmp == qqm1) {
        qqm1_tmp  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
        qqm1_tmp  -= i;

        for (k = i + 2; k <= j - 1; k++)
          qqm1_tmp[k] = qqm1[k];
      }

      for (k = i + 2; k <= j - 1; k++)
        qqm1_tmp[k] *= sc_wrapper.decomp_ml(i + 1, j - 1, k - 1, k, &sc_wrapper);
    }

    temp = 0.0;

    /* set initial decomposition split point */
    k = i + 2;

    if (sliding_window) {
      for (; k <= j - 1; k++, kl--)
        temp += qm_local[i + 1][k - 1] *
                qqm1_tmp[k];
    } else {
      kl = my_iindx[i + 1] - (i + 1);
      /*
       *  loop over entire range but skip decompositions with in-between strand nick,
       *  this should be faster than evaluating hard constraints callback for each
       *  decomposition
       */
      while (1) {
        /* limit for-loop to last nucleotide of 5' part strand */
        int stop = MIN2(j - 1, se[sn[k - 1]]);

        for (; k <= stop; k++, kl--)
          temp += qm[kl] *
                  qqm1_tmp[k];

        k++;
        kl--;

        if (stop == j - 1)
          break;
      }
    }

    if (qqm1_tmp != qqm1) {
      qqm1_tmp += i;
      free(qqm1_tmp);
    }

    qbt1 += temp *
            qqqmmm;
  }

  free_sc_mb_exp(&sc_wrapper);

  return qbt1;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_mb_loop_fast(vrna_fold_compound_t        *fc,
                        int                         i,
                        int                         j,
                        struct vrna_mx_pf_aux_ml_s  *aux_mx)
{
  FLT_OR_DBL q = 0.;

  if ((fc) && (aux_mx))
    q = exp_E_mb_loop_fast(fc, i, j, aux_mx);

  return q;
}



PRIVATE FLT_OR_DBL
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               vrna_mx_pf_aux_ml_t  aux_mx_ml)
{
  unsigned int  n;
  int           *jindx, *pscore;
  FLT_OR_DBL    contribution;
  double        kTn;
  vrna_hc_t     *hc;

  contribution  = 0.;
  n             = fc->length;
  hc            = fc->hc;

  if (hc->mx[j * n + i]) {
    /* process hairpin loop(s) */
    contribution += vrna_exp_E_hp_loop(fc, i, j);
    /* process interior loop(s) */
    contribution += vrna_exp_E_int_loop(fc, i, j);
    /* process multibranch loop(s) */
    contribution += vrna_exp_E_mb_loop_fast(fc, i, j, aux_mx_ml);

    if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_c))
      contribution += fc->aux_grammar->cb_aux_exp_c(fc, i, j, fc->aux_grammar->data);
  }

  return contribution;
}


#define VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP 8U
#define VRNA_UNSTRUCTURED_DOMAIN_MOTIF 16U

PRIVATE FLT_OR_DBL
exp_E_ml_fast(vrna_fold_compound_t        *fc,
              int                         i,
              int                         j,
              struct vrna_mx_pf_aux_ml_s  *aux_mx)
{
  unsigned char             sliding_window;
  short                     *S1, *S2, **SS, **S5, **S3;
  unsigned int              *sn, *ss, *se, n_seq, s;
  int                       n, *iidx, k, ij, kl, maxk, ii, with_ud, u, circular, with_gquad,
                            *hc_up_ml, type;
  FLT_OR_DBL                qbt1, temp, *qm, *qb, *qqm, *qqm1, **qqmu, q_temp, q_temp2, *G,
                            *expMLbase, **qb_local, **qm_local, **G_local;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_hc_eval_f evaluate;
  struct hc_mb_def_dat      hc_dat_local;
  struct sc_mb_exp_dat      sc_wrapper;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n               = (int)fc->length;
  sn              = fc->strand_number;
  ss              = fc->strand_start;
  se              = fc->strand_end;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  iidx            = (sliding_window) ? NULL : fc->iindx;
  ij              = (sliding_window) ? 0 : iidx[i] - j;
  qqm             = aux_mx->qqm;
  qqm1            = aux_mx->qqm1;
  qqmu            = aux_mx->qqmu;
  qm              = (sliding_window) ? NULL : fc->exp_matrices->qm;
  qb              = (sliding_window) ? NULL : fc->exp_matrices->qb;
  G               = (sliding_window) ? NULL : fc->exp_matrices->G;
  qm_local        = (sliding_window) ? fc->exp_matrices->qm_local : NULL;
  qb_local        = (sliding_window) ? fc->exp_matrices->qb_local : NULL;
  G_local         = (sliding_window) ? fc->exp_matrices->G_local : NULL;
  expMLbase       = fc->exp_matrices->expMLbase;
  pf_params       = fc->exp_params;
  md              = &(pf_params->model_details);
  hc              = fc->hc;
  domains_up      = fc->domains_up;
  circular        = md->circ;
  with_gquad      = md->gquad;
  with_ud         = (domains_up && domains_up->exp_energy_cb);
  hc_up_ml        = hc->up_ml;
  evaluate        = prepare_hc_mb_def(fc, &hc_dat_local);

  init_sc_mb_exp(fc, &sc_wrapper);

  qbt1    = 0;
  q_temp  = 0.;

  qqm[i] = 0.;

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (sc_wrapper.red_ml)
      q_temp *= sc_wrapper.red_ml(i, j, i, j - 1, &sc_wrapper);

    qqm[i] += q_temp;
  }

  if (with_ud) {
    q_temp = 0.;

    int cnt;
    for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
      u = domains_up->uniq_motif_size[cnt];
      if (j - u >= i) {
        if (evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          q_temp2 = qqmu[u][i] *
                    domains_up->exp_energy_cb(fc,
                                              j - u + 1,
                                              j,
                                              VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP |
                                              VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                              domains_up->data) *
                    expMLbase[u];

          if (sc_wrapper.red_ml)
            q_temp2 *= sc_wrapper.red_ml(i, j, i, j - u, &sc_wrapper);

          q_temp += q_temp2;
        }
      }
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    qbt1 = (sliding_window) ? qb_local[i][j] : qb[ij];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S1    = fc->sequence_encoding;
        S2    = fc->sequence_encoding2;
        type  = vrna_get_ptype_md(S2[i], S2[j], md);

        qbt1 *= exp_E_MLstem(type,
                             ((i > 1) || circular) ? S1[i - 1] : -1,
                             ((j < n) || circular) ? S1[j + 1] : -1,
                             pf_params);

        break;
    }

    if (sc_wrapper.red_stem)
      qbt1 *= sc_wrapper.red_stem(i, j, i, j, &sc_wrapper);

    qqm[i] += qbt1;
  }

  if (with_gquad) {
    q_temp  = (sliding_window) ? G_local[i][j] : G[ij];
    qqm[i]  += q_temp *
               pow(exp_E_MLstem(0, -1, -1, pf_params), (double)n_seq);
  }

  if (with_ud)
    qqmu[0][i] = qqm[i];

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  FLT_OR_DBL *qqm_tmp = qqm;

  /* apply hard constraints if necessary */
  if (hc->f) {
    qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
    qqm_tmp -= i;

    for (k = j; k > i; k--) {
      qqm_tmp[k] = qqm[k];
      if (!evaluate(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, &hc_dat_local))
        qqm_tmp[k] = 0.;
    }
  }

  /* apply soft constraints if necessary */
  if (sc_wrapper.decomp_ml) {
    if (qqm_tmp == qqm) {
      qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm_tmp -= i;

      for (k = j; k > i; k--)
        qqm_tmp[k] = qqm[k];
    }

    for (k = j; k > i; k--)
      qqm_tmp[k] *= sc_wrapper.decomp_ml(i, j, k - 1, k, &sc_wrapper);
  }

  /* finally, decompose segment */
  temp  = 0.0;
  k     = j;

  if (sliding_window) {
    for (; k > i; k--)
      temp += qm_local[i][k - 1] *
              qqm_tmp[k];
  } else {
    kl = iidx[i] - j + 1; /* ii-k=[i,k-1] */

    while (1) {
      /* limit for-loop to first nucleotide of 3' part strand */
      int stop = MAX2(i, ss[sn[k]]);
      for (; k > stop; k--, kl++)
        temp += qm[kl] *
                qqm_tmp[k];

      k--;
      kl++;

      if (stop == i)
        break;
    }
  }

  /* determine last nucleotide position for unpaired stretch */
  maxk = j;

  /* obey hard constraints (those that we can simply look-up) */
  if (maxk > i + hc_up_ml[i])
    maxk = i + hc_up_ml[i];

  /* obey connected components constraint (for multi-RNA folding) */
  if (maxk > se[sn[i]])
    maxk = se[sn[i]];

  if (qqm_tmp != qqm) {
    /* we've been using this helper array, so it's likely we use it again... */
    for (k = maxk; k > i; k--)
      qqm_tmp[k] = qqm[k];
  }

  /* apply hard constraints if necessary */
  if (hc->f) {
    if (qqm_tmp == qqm) {
      qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm_tmp -= i;

      for (k = maxk; k > i; k--)
        qqm_tmp[k] = qqm[k];
    }

    for (k = maxk; k > i; k--)
      if (!evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, &hc_dat_local))
        qqm_tmp[k] = 0.;
  }

  /* apply soft constraints if necessary */
  if (sc_wrapper.red_ml) {
    if (qqm_tmp == qqm) {
      qqm_tmp = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 2));
      qqm_tmp -= i;

      for (k = maxk; k > i; k--)
        qqm_tmp[k] = qqm[k];
    }

    for (k = maxk; k > i; k--)
      qqm_tmp[k] *= sc_wrapper.red_ml(i, j, k, j, &sc_wrapper);
  }

  ii = maxk - i; /* length of unpaired stretch */

  /* finally, decompose segment */
  for (k = maxk; k > i; k--, ii--)
    temp += expMLbase[ii] *
            qqm_tmp[k];

  if (with_ud) {
    ii = maxk - i; /* length of unpaired stretch */

    for (k = maxk; k > i; k--, ii--)
      temp += expMLbase[ii] *
              qqm_tmp[k] *
              domains_up->exp_energy_cb(fc,
                                        i, k - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                        domains_up->data);
  }

  if (qqm_tmp != qqm) {
    qqm_tmp += i;
    free(qqm_tmp);
  }

  /* apply auxiliary grammar rule for multibranch loop case */
  if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m))
    temp += fc->aux_grammar->cb_aux_exp_m(fc, i, j, fc->aux_grammar->data);

  free_sc_mb_exp(&sc_wrapper);

  return temp + qqm[i];
}




PUBLIC FLT_OR_DBL
vrna_exp_E_ml_fast(vrna_fold_compound_t       *fc,
                   int                        i,
                   int                        j,
                   struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  FLT_OR_DBL q = 0.;

  if ((fc) && (aux_mx))
    q = exp_E_ml_fast(fc, i, j, aux_mx);

  return q;
}

PUBLIC const FLT_OR_DBL *
vrna_exp_E_ml_fast_qqm(struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  if (aux_mx)
    return (const FLT_OR_DBL *)aux_mx->qqm;

  return NULL;
}




PRIVATE INLINE FLT_OR_DBL
reduce_ext_ext_fast(vrna_fold_compound_t        *fc,
                    int                         i,
                    int                         j,
                    struct vrna_mx_pf_aux_el_s  *aux_mx,
                    vrna_hc_eval_f   evaluate,
                    struct hc_ext_def_dat       *hc_dat_local,
                    struct sc_ext_exp_dat       *sc_wrapper)
{
  int           u;
  FLT_OR_DBL    q_temp, q_temp2, q, *qq1, **qqu, *scale;
  vrna_ud_t     *domains_up;
  sc_ext_exp_cb sc_red_ext;

  domains_up  = fc->domains_up;
  qq1         = aux_mx->qq1;
  qqu         = aux_mx->qqu;
  scale       = fc->exp_matrices->scale;
  sc_red_ext  = sc_wrapper->red_ext;

  q = 0.;

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
    q_temp = qq1[i] * scale[1];

    if (sc_red_ext)
      q_temp *= sc_red_ext(i, j, i, j - 1, sc_wrapper);

    if ((domains_up) && (domains_up->exp_energy_cb)) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_EXT_EXT, hc_dat_local)) {
            q_temp2 = qqu[u][i] *
                      domains_up->exp_energy_cb(fc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      scale[u];

            if (sc_red_ext)
              q_temp2 *= sc_red_ext(i, j, i, j - u, sc_wrapper);

            q_temp += q_temp2;
          }
        }
      }
    }

    q = q_temp;
  }

  return q;
}


PRIVATE INLINE FLT_OR_DBL
reduce_ext_stem_fast(vrna_fold_compound_t       *fc,
                     int                        i,
                     int                        j,
                     struct vrna_mx_pf_aux_el_s *aux_mx,
                     vrna_hc_eval_f  evaluate,
                     struct hc_ext_def_dat      *hc_dat_local,
                     struct sc_ext_exp_dat      *sc_wrapper)
{
  short             **S, **S5, **S3, *S1, *S2, s5, s3;
  unsigned int      type, *sn, n, s, n_seq, **a2s;
  int               *idx, circular;
  FLT_OR_DBL        qbt, q_temp, qb;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;
  sc_ext_exp_cb     sc_red_stem;

  sc_red_stem = sc_wrapper->red_stem;
  n           = fc->length;
  sn          = fc->strand_number;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  circular    = md->circ;
  idx         = fc->iindx;
  qb          = (fc->hc->type == VRNA_HC_WINDOW) ?
                fc->exp_matrices->qb_local[i][j] :
                fc->exp_matrices->qb[idx[i] - j];
  qbt = 0.;

  /* exterior loop part with stem (i, j) */
  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc_dat_local)) {
    q_temp = qb;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S1      = fc->sequence_encoding;
        S2      = fc->sequence_encoding2;
        type    = vrna_get_ptype_md(S2[i], S2[j], md);
        s5      = (((i > 1) || circular) && (sn[i] == sn[i - 1])) ? S1[i - 1] : -1;
        s3      = (((j < n) || circular) && (sn[j + 1] == sn[j])) ? S1[j + 1] : -1;
        q_temp  *= vrna_exp_E_ext_stem(type, s5, s3, pf_params);
        break;
    }

    if (sc_red_stem)
      q_temp *= sc_red_stem(i, j, i, j, sc_wrapper);

    qbt += q_temp;
  }

  return qbt;
}



PRIVATE INLINE FLT_OR_DBL
split_ext_fast(vrna_fold_compound_t       *fc,
               int                        i,
               int                        j,
               struct vrna_mx_pf_aux_el_s *aux_mx,
               vrna_hc_eval_f  evaluate,
               struct hc_ext_def_dat      *hc_dat_local,
               struct sc_ext_exp_dat      *sc_wrapper)
{
  int               *idx, k, ij1;
  FLT_OR_DBL        qbt, *q, *qq, *qqq;
  sc_ext_exp_split  sc_split;

  sc_split = sc_wrapper->split;

  idx = fc->iindx;
  q   = (fc->hc->type == VRNA_HC_WINDOW) ?
        fc->exp_matrices->q_local[i] :
        fc->exp_matrices->q + idx[i];
  qq  = aux_mx->qq;
  qbt = 0.;

  /*
   *  use an auxiliary array qqq that already contains soft constraint
   *  contributions when we split the exterior loops
   */
  if (sc_split) {
    /* pre-process qq array */
    qqq = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (j - i + 1));
    qqq -= i;

    for (k = j; k > i; k--)
      qqq[k] = qq[k] * sc_split(i, j, k, sc_wrapper);
  } else {
    /* pre-process qq array */
    qqq = qq;
  }

  /*
   *  the index for access in q array now is:
   *
   *  (a) q[- (j - 1)] for global, and
   *  (b) q[j - 1] for local structure prediction
   *
   *  Thus, we may use a single variable to address both cases:
   *  k = j:
   *    ij1 = (j - 1) * factor = j * factor - factor
   *  k = j - 1:
   *    ij1 = (j - 2) * factor = j * factor - 2 * factor = j * factor - factor - factor
   *  ...
   */
  int factor = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : -1;

  ij1 = factor * (j - 1);

  /* do actual decomposition (skip hard constraint checks if we use default settings) */
#if SPEEDUP_HC
  /*
   *  checking whether we actually are provided with hard constraints and
   *  otherwise not evaluating the default ones within the loop drastically
   *  increases speed. However, once we check for the split point between
   *  strands in hard constraints, we have to think of something else...
   */
  if ((evaluate == &hc_ext_cb_def) || (evaluate == &hc_ext_cb_def_window)) {
    for (k = j; k > i; k--) {
      qbt += q[ij1] *
             qqq[k];
      ij1 -= factor;
    }
  } else {
    for (k = j; k > i; k--) {
      if (evaluate(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, hc_dat_local))
        qbt += q[ij1] *
               qqq[k];

      ij1 -= factor;
    }
  }

#else
  for (k = j; k > i; k--) {
    if (evaluate(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, hc_dat_local))
      qbt += q[ij1] *
             qqq[k];

    ij1 -= factor;
  }
#endif

  if (qqq != qq) {
    qqq += i;
    free(qqq);
  }

  return qbt;
}

PRIVATE INLINE void
free_sc_ext_exp(struct sc_ext_exp_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}

PRIVATE FLT_OR_DBL
exp_E_ext_fast(vrna_fold_compound_t       *fc,
               int                        i,
               int                        j,
               struct vrna_mx_pf_aux_el_s *aux_mx)
{
  int                       *iidx, ij, with_ud, with_gquad;
  FLT_OR_DBL                qbt1, *qq, **qqu, *G, **G_local;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;
  struct sc_ext_exp_dat     sc_wrapper;

  qq          = aux_mx->qq;
  qqu         = aux_mx->qqu;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  domains_up  = fc->domains_up;
  with_gquad  = md->gquad;
  with_ud     = (domains_up && domains_up->exp_energy_cb);

  if (fc->hc->type == VRNA_HC_WINDOW) // not into
    evaluate = prepare_hc_ext_def_window(fc, &hc_dat_local);
  else
    evaluate = prepare_hc_ext_def(fc, &hc_dat_local);

  init_sc_ext_exp(fc, &sc_wrapper);

  qbt1 = 0.;

  /* all exterior loop parts [i, j] with exactly one stem (i, u) i < u < j */
  qbt1 += reduce_ext_ext_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);
  /* exterior loop part with stem (i, j) */
  qbt1 += reduce_ext_stem_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);

  if (with_gquad) {
    if (fc->hc->type == VRNA_HC_WINDOW) {
      G_local = fc->exp_matrices->G_local;
      qbt1    += G_local[i][j];
    } else {
      G     = fc->exp_matrices->G;
      iidx  = fc->iindx;
      ij    = iidx[i] - j;
      qbt1  += G[ij];
    }
  }

  qq[i] = qbt1;

  if (with_ud)
    qqu[0][i] = qbt1;

  /* the entire stretch [i,j] is unpaired */
  qbt1 += reduce_ext_up_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);

  qbt1 += split_ext_fast(fc, i, j, aux_mx, evaluate, &hc_dat_local, &sc_wrapper);

  /* apply auxiliary grammar rule for exterior loop case */
  if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_f))
    qbt1 += fc->aux_grammar->cb_aux_exp_f(fc, i, j, fc->aux_grammar->data);

  free_sc_ext_exp(&sc_wrapper);

  return qbt1;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ext_fast(vrna_fold_compound_t        *fc,
                    int                         i,
                    int                         j,
                    struct vrna_mx_pf_aux_el_s  *aux_mx)
{
  if (fc) {
    if (j < i) {
      int t = j;
      printf(
        "vrna_exp_E_ext_fast: i (%d) larger than j (%d)! Swapping coordinates...",
        i,
        j);
      j = i;
      i = t;
    } else if ((j < 1) || (i < 1)) {
      printf(
        "vrna_exp_E_ext_fast: Indices too small [i = %d, j = %d]! "
        "Refusing to compute anything...",
        i,
        j);
      return 0.;
    } else if (j > fc->length) {
      printf(
        "vrna_exp_E_ext_fast: Indices exceed sequence length (%d) [i = %d, j = %d]! "
        "Refusing to compute anything...",
        fc->length,
        i,
        j);
      return 0.;
    }

    return exp_E_ext_fast(fc, i, j, aux_mx);
  }

  return 0.;
}




PUBLIC void
vrna_exp_E_ext_fast_rotate(struct vrna_mx_pf_aux_el_s *aux_mx)
{
  if (aux_mx) {
    int         u;
    FLT_OR_DBL  *tmp;

    tmp         = aux_mx->qq1;
    aux_mx->qq1 = aux_mx->qq;
    aux_mx->qq  = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqu) {
      tmp = aux_mx->qqu[aux_mx->qqu_size];
      for (u = aux_mx->qqu_size; u > 0; u--)
        aux_mx->qqu[u] = aux_mx->qqu[u - 1];
      aux_mx->qqu[0] = tmp;
    }
  }
}



PUBLIC void
vrna_exp_E_ml_fast_rotate(struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  if (aux_mx) {
    int         u;
    FLT_OR_DBL  *tmp;

    tmp           = aux_mx->qqm1;
    aux_mx->qqm1  = aux_mx->qqm;
    aux_mx->qqm   = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqmu) {
      tmp = aux_mx->qqmu[aux_mx->qqmu_size];
      for (u = aux_mx->qqmu_size; u > 0; u--)
        aux_mx->qqmu[u] = aux_mx->qqmu[u - 1];
      aux_mx->qqmu[0] = tmp;
    }
  }
}




PUBLIC void
vrna_exp_E_ml_fast_free(struct vrna_mx_pf_aux_ml_s *aux_mx)
{
  if (aux_mx) {
    int u;

    free(aux_mx->qqm);
    free(aux_mx->qqm1);

    if (aux_mx->qqmu) {
      for (u = 0; u <= aux_mx->qqmu_size; u++)
        free(aux_mx->qqmu[u]);

      free(aux_mx->qqmu);
    }

    free(aux_mx);
  }
}


PUBLIC void
vrna_exp_E_ext_fast_free(struct vrna_mx_pf_aux_el_s *aux_mx)
{
  if (aux_mx) {
    int u;

    free(aux_mx->qq);
    free(aux_mx->qq1);

    if (aux_mx->qqu) {
      for (u = 0; u <= aux_mx->qqu_size; u++)
        free(aux_mx->qqu[u]);

      free(aux_mx->qqu);
    }

    free(aux_mx);
  }
}



PRIVATE int
fill_arrays(vrna_fold_compound_t *fc)
{ 
    
  int                 n, i, j, k, ij, *my_iindx, *jindx, with_gquad, with_ud;
  FLT_OR_DBL          temp, Qmax, *q, *qb, *qm, *qm1, *q1k, *qln;
  double              max_real;
  vrna_ud_t           *domains_up;
  vrna_md_t           *md;
  vrna_mx_pf_t        *matrices;
  vrna_mx_pf_aux_el_t aux_mx_el;
  vrna_mx_pf_aux_ml_t aux_mx_ml;
  vrna_exp_param_t    *pf_params;

  n           = fc->length;
  my_iindx    = fc->iindx;
  jindx       = fc->jindx;
  matrices    = fc->exp_matrices;
  pf_params   = fc->exp_params;
  domains_up  = fc->domains_up;
  q           = matrices->q;
  qb          = matrices->qb;
  qm          = matrices->qm;
  qm1         = matrices->qm1;
  q1k         = matrices->q1k;
  qln         = matrices->qln;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  with_ud = (domains_up && domains_up->exp_energy_cb && (!(fc->type == VRNA_FC_TYPE_COMPARATIVE)));
  Qmax    = 0;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (with_ud && domains_up->exp_prod_cb) // not into
    domains_up->exp_prod_cb(fc, domains_up->data);

  /* init auxiliary arrays for fast exterior/multibranch loops */
  aux_mx_el = vrna_exp_E_ext_fast_init(fc);
  aux_mx_ml = vrna_exp_E_ml_fast_init(fc);

  /*array initialization ; qb,qm,q
   * qb,qm,q (i,j) are stored as ((n+1-i)*(n-i) div 2 + n+1-j */
  for (i = 1; i <= n; i++) {
    ij      = my_iindx[i] - i;
    qb[ij]  = 0.0;
  }

  for (j = 2; j <= n; j++) {
    for (i = j - 1; i >= 1; i--) {
      ij = my_iindx[i] - j;

      qb[ij] = decompose_pair(fc, i, j, aux_mx_ml);

      /* Multibranch loop */
      qm[ij] = vrna_exp_E_ml_fast(fc, i, j, aux_mx_ml);

      if (qm1) { // not into
        temp = vrna_exp_E_ml_fast_qqm(aux_mx_ml)[i]; /* for stochastic backtracking and circfold */

        /* apply auxiliary grammar rule for multibranch loop (M1) case */
        if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp_m1))
          temp += fc->aux_grammar->cb_aux_exp_m1(fc, i, j, fc->aux_grammar->data);

        qm1[jindx[j] + i] = temp;
      }

      /* Exterior loop */
      q[ij] = vrna_exp_E_ext_fast(fc, i, j, aux_mx_el);

      /* apply auxiliary grammar rule (storage takes place in user-defined data structure */
      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_exp)) // not into
        fc->aux_grammar->cb_aux_exp(fc, i, j, fc->aux_grammar->data);

      if (q[ij] > Qmax) {
        Qmax = q[ij];
        if (Qmax > max_real / 10.) // not into
          printf("Q close to overflow: %d %d %g", i, j, q[ij]);
      }

      if (q[ij] >= max_real) { // not into
        printf("overflow while computing partition function for segment q[%d,%d]\n"
                             "use larger pf_scale", i, j);

        vrna_exp_E_ml_fast_free(aux_mx_ml);
        vrna_exp_E_ext_fast_free(aux_mx_el);

        return 0; /* failure */
      }
    }

    /* rotate auxiliary arrays */
    vrna_exp_E_ext_fast_rotate(aux_mx_el);
    vrna_exp_E_ml_fast_rotate(aux_mx_ml);
  }

  /* prefill linear qln, q1k arrays */
  if (q1k && qln) {
    for (k = 1; k <= n; k++) {
      q1k[k]  = q[my_iindx[1] - k];
      qln[k]  = q[my_iindx[k] - n];
    }
    q1k[0]      = 1.0;
    qln[n + 1]  = 1.0;
  }

  /* free memory occupied by auxiliary arrays for fast exterior/multibranch loops */
  vrna_exp_E_ml_fast_free(aux_mx_ml);
  vrna_exp_E_ext_fast_free(aux_mx_el);

  return 1;
}


/**/
PRIVATE void
postprocess_circular(vrna_fold_compound_t *fc)
{
  unsigned int      **a2s;
  int               u, p, q, k, turn, n, *my_iindx, *jindx, s;
  FLT_OR_DBL        *scale, *qb, *qm, *qm1, *qm2, qo, qho, qio, qmo,
                    qbt1, qot, expMLclosing, n_seq;
  unsigned char     eval;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc, **scs;

  n             = fc->length;
  n_seq         = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  matrices      = fc->exp_matrices;
  my_iindx      = fc->iindx;
  jindx         = fc->jindx;
  pf_params     = fc->exp_params;
  hc            = fc->hc;
  qb            = matrices->qb;
  qm            = matrices->qm;
  qm1           = matrices->qm1;
  qm2           = matrices->qm2;
  scale         = matrices->scale;
  expMLclosing  = pf_params->expMLclosing;
  turn          = pf_params->model_details.min_loop_size;
  hc            = fc->hc;
  sc            = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sc : NULL;
  scs           = NULL;
  a2s           = NULL;
  qo            = qho = qio = qmo = 0.;

  for (p = 1; p < n; p++) {
    for (q = p + turn + 1; q <= n; q++) {
      u = n - q + p - 1;
      if (u < turn)
        continue;

      /* 1. get exterior hairpin contribution  */
      qho += qb[my_iindx[p] - q] *
             vrna_exp_E_hp_loop(fc, q, p);

      /* 2. get exterior interior loop contribution */
      qio += qb[my_iindx[p] - q] *
             vrna_exp_E_int_loop(fc, q, p);
    }
  } /* end of pq double loop */

  /* 3. Multiloops  */

  /* construct qm2 matrix for exterior multibranch loop computation */
  if (hc->f) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              if (hc->f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
                qot += qm1[jindx[u] + k] *
                       qm1[jindx[n] + (u + 1)] *
                       sc->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              }

            qm2[k] = qot;
          }
        } else {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              if (hc->f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
                qot += qm1[jindx[u] + k] *
                       qm1[jindx[n] + (u + 1)];

            qm2[k] = qot;
          }
        }

        break;

    }
  } else {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              qot += qm1[jindx[u] + k] *
                     qm1[jindx[n] + (u + 1)] *
                     sc->exp_f(k, n, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            qm2[k] = qot;
          }
        } else {
          for (k = 1; k < n - turn - 1; k++) {
            qot = 0.;

            for (u = k + turn + 1; u < n - turn - 1; u++)
              qot += qm1[jindx[u] + k] *
                     qm1[jindx[n] + (u + 1)];

            qm2[k] = qot;
          }
        }

        break;

    }
  }

  qbt1 = 0.;
  /* go through exterior multibranch loop configurations */
  if (hc->f) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            if (hc->f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
              qbt1 += qm[my_iindx[1] - k] *
                      qm2[k + 1] *
                      sc->exp_f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        } else {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            if (hc->f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data))
              qbt1 += qm[my_iindx[1] - k] *
                      qm2[k + 1];
        }

        qbt1 *= expMLclosing;
        break;
    }
  } else {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if ((sc) && (sc->exp_f)) {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            qbt1 += qm[my_iindx[1] - k] *
                    qm2[k + 1] *
                    sc->exp_f(1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
        } else {
          for (k = turn + 2; k < n - 2 * turn - 3; k++)
            qbt1 += qm[my_iindx[1] - k] *
                    qm2[k + 1];
        }

        qbt1 *= expMLclosing;
        break;

    }
  }

  qmo += qbt1;

  /* add an additional pf of 1.0 to take the open chain into account too */
  eval = (hc->up_ext[1] >= n) ? 1 : 0;
  if (hc->f)
    eval = (hc->f(1, n, 1, n, VRNA_DECOMP_EXT_UP, hc->data)) ? eval : 0;

  if (eval) {
    qbt1 = scale[n];

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (sc) {
          if (sc->exp_energy_up)
            qbt1 *= sc->exp_energy_up[1][n];

          if (sc->exp_f)
            qbt1 *= sc->exp_f(1, n, 1, n, VRNA_DECOMP_EXT_UP, sc->data);
        }

        break;

    }
    qo += qbt1;
  }

  qo += qho + qio + qmo;

  matrices->qo  = qo;
  matrices->qho = qho;
  matrices->qio = qio;
  matrices->qmo = qmo;
}


PUBLIC int
vrna_gr_reset(vrna_fold_compound_t *fc)
{
  int ret = 0;

  if ((fc) && (fc->aux_grammar)) {
    if (fc->aux_grammar->free_data)
      fc->aux_grammar->free_data(fc->aux_grammar->data);

    free(fc->aux_grammar);
    fc->aux_grammar = NULL;
  }

  return ret;
}




PRIVATE size_t *
get_BM_BCT(const char *needle,
           size_t     needle_size)
{
  size_t *table, i;

  table = vrna_alloc(sizeof(size_t) * (127 + 2));

  /* store maximum element value at position 0 */
  table[0] = 127;

  /* use remainder of array for actual bad character table */
  for (i = 1; i <= 127 + 1; i++)
    table[i] = needle_size;

  for (i = 0; i < needle_size - 1; i++)
    table[needle[i] + 1] = needle_size - i - 1;

  return table;
}



#define BMH { \
    hit     = NULL; \
    shift   = start; \
    margin  = (cyclic) ? 0 : needle_size; \
    max     = bad_chars[0]; \
    /* pop first element since the Bad Character Table starts at position 1 */ \
    bad_chars++; \
    /* main loop - go through haystack */ \
    while (shift + margin < haystack_size) { \
      /* \
       *  matching loop, note that we allow for possibly wrapping the \
       *  pattern around the haystack \
       */\
      for (i = needle_size - 1; \
           haystack[(shift + i) % haystack_size] == needle[i]; \
           i--) { \
        if (i == 0) \
        return haystack + shift; \
      } \
      val = haystack[(shift + needle_size - 1) % haystack_size]; \
      shift += bad_chars[(size_t)val]; \
    } \
}



PRIVATE const char *
BoyerMooreHorspool(const char     *needle,
                   size_t         needle_size,
                   const char     *haystack,
                   size_t         haystack_size,
                   size_t         start,
                   size_t         *bad_chars,
                   unsigned char  cyclic)
{
  const char  *hit;
  char        val, max;
  size_t      i, shift, margin;

  /* empty pattern matches element in haystack */
  if (!needle)
    return haystack;

  /* empty pattern matches element in haystack */
  if (needle_size == 0)
    return haystack;

  /* empty haystack can't contain any pattern */
  if (haystack_size == 0)
    return NULL;

  /* haystack mustn't be shorter than needle */
  if (haystack_size < needle_size)
    return NULL;

  /* begin actual algorithm */
  BMH;

  return hit;
}




PUBLIC const char *
vrna_search_BMH(const char    *needle,
                size_t        needle_size,
                const char    *haystack,
                size_t        haystack_size,
                size_t        start,
                size_t        *badchars,
                unsigned char cyclic)
{
  const char  *hit;
  size_t      *bc;

  if ((!needle) || (!haystack) || (start > haystack_size))
    return NULL;

  bc = badchars;

  /* create bad character table in case none is supplied */
  if (!bc)
    bc = get_BM_BCT(needle, needle_size);

  /* perform actual search */
  hit = BoyerMooreHorspool(needle,
                           needle_size,
                           haystack,
                           haystack_size,
                           start,
                           bc,
                           cyclic);

  if (bc != badchars)
    free(bc);

  return hit;
}


PUBLIC size_t *
vrna_search_BM_BCT(const char *pattern)
{
  if (!pattern)
    return NULL;

  size_t pattern_size = strlen(pattern);

  return get_BM_BCT(pattern, pattern_size);
}


PUBLIC unsigned int
vrna_rotational_symmetry_pos(const char   *string,
                             unsigned int **positions)
{
  const char    *ptr;
  unsigned int  matches, shifts_size;
  size_t        *badchars, shift, i, string_length;

  if (!string) {
    if (positions)
      *positions = NULL;

    return 0;
  }

  string_length = strlen(string);

  if (string_length == 0) {
    if (positions)
      *positions = NULL;

    return 0;
  }

  /* any string is at least order 1 */
  matches = 1;

  if (positions) {
    shifts_size = 10; /* initial guess for the order of rotational symmetry */
    *positions  = vrna_alloc(sizeof(unsigned int) * shifts_size);

    /* store trivial symmetry */
    (*positions)[matches - 1] = 0;
  }

  /* strings of length 1 are order 1 */
  if (string_length == 1) {
    /* resize positions array to actual length */
    if (positions)
      *positions = vrna_realloc(*positions, sizeof(unsigned int) * matches);

    return matches;
  }

  /* determine largest number/character in string */
  badchars = vrna_search_BM_BCT(string);

  shift = 1; /* skip trivial symmetry */

  /* detect order of rotational symmetry */

  /*
   *  Note, that finding the smallest shift s of the string that
   *  results in an identity mapping of the string to itself
   *  already determines the order of rotational symmetry R, i.e.
   *  R = n / r where n is the length of the string
   */
  ptr = vrna_search_BMH(string,
                        string_length,
                        string,
                        string_length,
                        shift,
                        badchars,
                        1);

  if (ptr) {
    shift   = ptr - string;
    matches = string_length / shift;
    if (positions) {
      *positions = vrna_realloc(*positions, sizeof(unsigned int) * matches);
      for (i = 0; i < matches; i++)
        (*positions)[i] = i * shift;
    }
  }

  free(badchars);

  return matches;
}



PUBLIC unsigned int
vrna_rotational_symmetry(const char *string)
{
  return vrna_rotational_symmetry_pos(string,
                                      NULL);
}





PRIVATE helper_arrays *
get_ml_helper_arrays(vrna_fold_compound_t *fc)
{
  unsigned int  n, u;
  int           with_ud;
  vrna_ud_t     *domains_up;
  helper_arrays *ml_helpers;

  n           = fc->length;
  domains_up  = fc->domains_up;
  with_ud     = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;

  ml_helpers = (helper_arrays *)vrna_alloc(sizeof(helper_arrays));

  ml_helpers->prm_l   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  ml_helpers->prm_l1  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
  ml_helpers->prml    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

  ml_helpers->ud_max_size = 0;
  ml_helpers->pmlu        = NULL;
  ml_helpers->prm_MLbu    = NULL;

  if (with_ud) {
    /* find out maximum size of any unstructured domain */
    for (u = 0; u < domains_up->uniq_motif_count; u++)
      if (ml_helpers->ud_max_size < domains_up->uniq_motif_size[u])
        ml_helpers->ud_max_size = domains_up->uniq_motif_size[u];

    ml_helpers->pmlu = (FLT_OR_DBL **)vrna_alloc(
      sizeof(FLT_OR_DBL *) * (ml_helpers->ud_max_size + 1)); /* maximum motif size */

    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      ml_helpers->pmlu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));

    ml_helpers->prm_MLbu =
      (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (ml_helpers->ud_max_size + 1));

    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      ml_helpers->prm_MLbu[u] = 0.;
  }

  return ml_helpers;
}




PRIVATE constraints_helper *
get_constraints_helper(vrna_fold_compound_t *fc)
{
  constraints_helper *helpers;

  helpers = (constraints_helper *)vrna_alloc(sizeof(constraints_helper));

  helpers->hc_eval_ext  = prepare_hc_ext_def(fc, &(helpers->hc_dat_ext));
  helpers->hc_eval_hp   = prepare_hc_hp_def(fc, &(helpers->hc_dat_hp));
  helpers->hc_eval_int  = prepare_hc_int_def(fc, &(helpers->hc_dat_int));
  helpers->hc_eval_mb   = prepare_hc_mb_def(fc, &(helpers->hc_dat_mb));

  init_sc_ext_exp(fc, &(helpers->sc_wrapper_ext));
  init_sc_hp_exp(fc, &(helpers->sc_wrapper_hp));
  init_sc_int_exp(fc, &(helpers->sc_wrapper_int));
  init_sc_mb_exp(fc, &(helpers->sc_wrapper_mb));

  return helpers;
}




PRIVATE void
compute_bpp_internal(vrna_fold_compound_t *fc,
                     int                  l,
                     vrna_ep_t            **bp_correction,
                     int                  *corr_cnt,
                     int                  *corr_size,
                     FLT_OR_DBL           *Qmax,
                     int                  *ov,
                     constraints_helper   *constraints)
{
  unsigned char         type, type_2;
  char                  *ptype;
  short                 *S1;
  int                   i, j, k, n, ij, kl, u1, u2, *my_iindx, *jindx, *rtype,
                        with_ud, *hc_up_int;
  FLT_OR_DBL            temp, tmp2, *qb, *probs, *scale;
  double                max_real;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  vrna_sc_t             *sc;
  vrna_ud_t             *domains_up;
  struct hc_int_def_dat *hc_dat_local;
  eval_hc               hc_eval;
  struct sc_int_exp_dat *sc_wrapper_int;

  hc_dat_local    = &(constraints->hc_dat_int);
  hc_eval         = constraints->hc_eval_int;
  sc_wrapper_int  = &(constraints->sc_wrapper_int);

  n           = (int)fc->length;
  ptype       = fc->ptype;
  S1          = fc->sequence_encoding;
  my_iindx    = fc->iindx;
  jindx       = fc->jindx;
  pf_params   = fc->exp_params;
  md          = &(pf_params->model_details);
  rtype       = &(md->rtype[0]);
  hc          = fc->hc;
  sc          = fc->sc;
  domains_up  = fc->domains_up;
  with_ud     = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  hc_up_int   = hc->up_int;

  qb    = fc->exp_matrices->qb;
  probs = fc->exp_matrices->probs;
  scale = fc->exp_matrices->scale;

  max_real = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  /* 2. bonding k,l as substem of 2:loop enclosed by i,j */
  for (k = 1; k < l; k++) {
    kl = my_iindx[k] - l;

    if (qb[kl] == 0.)
      continue;

    if (hc->mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
      type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)];

      for (i = MAX2(1, k - MAXLOOP - 1); i <= k - 1; i++) {
        u1 = k - i - 1;
        if (hc_up_int[i + 1] < u1)
          continue;

        int max_j = l + 1 + MAXLOOP - u1;

        if (max_j > n)
          max_j = n;

        if (max_j > l + 1 + hc_up_int[l + 1])
          max_j = l + 1 + hc_up_int[l + 1];

        u2 = 0;

        for (j = l + 1; j <= max_j; j++, u2++) {
          ij = my_iindx[i] - j;

          if (probs[ij] == 0.)
            continue;

          if (hc_eval(i, j, k, l, hc_dat_local)) {
            int jij = jindx[j] + i;
            type  = vrna_get_ptype(jij, ptype);
            tmp2  = probs[ij] *
                    exp_E_IntLoop(u1,
                                  u2,
                                  type,
                                  type_2,
                                  S1[i + 1],
                                  S1[j - 1],
                                  S1[k - 1],
                                  S1[l + 1],
                                  pf_params) *
                    scale[u1 + u2 + 2];

            if (sc_wrapper_int->pair)
              tmp2 *= sc_wrapper_int->pair(i, j, k, l, sc_wrapper_int);

            if (with_ud) {
              FLT_OR_DBL qql, qqr;

              qql = qqr = 0.;

              if (u1 > 0) {
                qql = domains_up->exp_energy_cb(fc,
                                                i + 1, k - 1,
                                                VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                domains_up->data);
              }

              if (u2 > 0) {
                qqr = domains_up->exp_energy_cb(fc,
                                                l + 1, j - 1,
                                                VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                domains_up->data);
              }

              temp  = tmp2;
              tmp2  += temp * qql;
              tmp2  += temp * qqr;
              tmp2  += temp * qql * qqr;
            }

            if ((sc) &&
                (sc->exp_f) &&
                (sc->bt)) {
              /* store probability correction for auxiliary pairs in interior loop motif */
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = sc->bt(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                (*bp_correction)[*corr_cnt].i     = ptr->i;
                (*bp_correction)[*corr_cnt].j     = ptr->j;
                (*bp_correction)[(*corr_cnt)++].p = tmp2 * qb[kl];
                if ((*corr_cnt) == (*corr_size)) {
                  (*corr_size)      += 5;
                  (*bp_correction)  = vrna_realloc(*bp_correction,
                                                   sizeof(vrna_ep_t) * (*corr_size));
                }
              }
              free(aux_bps);
            }

            probs[kl] += tmp2;
          }
        }
      }
    }

    if (probs[kl] > (*Qmax)) {
      (*Qmax) = probs[kl];
      if ((*Qmax) > max_real / 10.)
        printf("P close to overflow: %d %d %g %g\n",
                             k, l, probs[kl], qb[kl]);
    }

    if (probs[kl] >= max_real) {
      (*ov)++;
      probs[kl] = FLT_MAX;
    }
  }
}

PRIVATE void
rotate_ml_helper_arrays_inner(helper_arrays *ml_helpers)
{
  unsigned int u;

  /* rotate prm_MLbu entries required for unstructured domain feature */
  if (ml_helpers->prm_MLbu)
    for (u = ml_helpers->ud_max_size; u > 0; u--)
      ml_helpers->prm_MLbu[u] = ml_helpers->prm_MLbu[u - 1];
}


PRIVATE void
rotate_ml_helper_arrays_outer(helper_arrays *ml_helpers)
{
  unsigned int  u;
  FLT_OR_DBL    *tmp;

  /* rotate prm_l and prm_l1 arrays */
  tmp                 = ml_helpers->prm_l1;
  ml_helpers->prm_l1  = ml_helpers->prm_l;
  ml_helpers->prm_l   = tmp;

  /* rotate pmlu entries required for unstructured domain feature */
  if (ml_helpers->pmlu) {
    tmp = ml_helpers->pmlu[ml_helpers->ud_max_size];

    for (u = ml_helpers->ud_max_size; u > 0; u--)
      ml_helpers->pmlu[u] = ml_helpers->pmlu[u - 1];

    ml_helpers->pmlu[0] = tmp;

    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      ml_helpers->prm_MLbu[u] = 0.;
  }
}



PRIVATE void
compute_bpp_multibranch(vrna_fold_compound_t  *fc,
                        int                   l,
                        helper_arrays         *ml_helpers,
                        FLT_OR_DBL            *Qmax,
                        int                   *ov,
                        constraints_helper    *constraints)
{
  unsigned char             tt;
  char                      *ptype;
  short                     *S, *S1, s5, s3;
  unsigned int              *sn;
  int                       cnt, i, j, k, n, u, ii, ij, kl, lj, *my_iindx, *jindx,
                            *rtype, with_gquad, with_ud;
  FLT_OR_DBL                temp, ppp, prm_MLb, prmt, prmt1, *qb, *probs, *qm, *G, *scale,
                            *expMLbase, expMLclosing, expMLstem;
  double                    max_real;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  struct hc_mb_def_dat      *hc_dat;
  vrna_hc_eval_f hc_eval;
  struct sc_mb_exp_dat      *sc_wrapper;


  n             = (int)fc->length;
  sn            = fc->strand_number;
  S             = fc->sequence_encoding2;
  S1            = fc->sequence_encoding;
  my_iindx      = fc->iindx;
  jindx         = fc->jindx;
  pf_params     = fc->exp_params;
  md            = &(pf_params->model_details);
  rtype         = &(md->rtype[0]);
  ptype         = fc->ptype;
  qb            = fc->exp_matrices->qb;
  qm            = fc->exp_matrices->qm;
  G             = fc->exp_matrices->G;
  probs         = fc->exp_matrices->probs;
  scale         = fc->exp_matrices->scale;
  expMLbase     = fc->exp_matrices->expMLbase;
  expMLclosing  = pf_params->expMLclosing;
  domains_up    = fc->domains_up;
  with_ud       = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  with_gquad    = md->gquad;
  expMLstem     = (with_gquad) ? exp_E_MLstem(0, -1, -1, pf_params) : 0;

  hc_dat      = &(constraints->hc_dat_mb);
  hc_eval     = constraints->hc_eval_mb;
  sc_wrapper  = &(constraints->sc_wrapper_mb);

  prm_MLb   = 0.;
  max_real  = (sizeof(FLT_OR_DBL) == sizeof(float)) ? FLT_MAX : DBL_MAX;

  if (sn[l + 1] != sn[l]) {
    /* set prm_l to 0 to get prm_l1 in the next round to be 0 */
    for (i = 0; i <= n; i++)
      ml_helpers->prm_l[i] = 0;
  } else {
    for (k = 2; k < l; k++) {
      kl    = my_iindx[k] - l;
      i     = k - 1;
      prmt  = prmt1 = 0.0;

      ij  = my_iindx[i] - (l + 2);
      lj  = my_iindx[l + 1] - (l + 1);
      s3  = S1[i + 1];
      if (sn[k] == sn[i]) {
        for (j = l + 2; j <= n; j++, ij--, lj--) {
          if (hc_eval(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc_dat)) {
            tt = vrna_get_ptype_md(S[j], S[i], md);

            /* which decomposition is covered here? =>
             * i + 1 = k < l < j:
             * (i,j)       -> enclosing pair
             * (k, l)      -> enclosed pair
             * (l+1, j-1)  -> multiloop part with at least one stem
             * a.k.a. (k,l) is left-most stem in multiloop closed by (k-1, j)
             */
            ppp = probs[ij] *
                  exp_E_MLstem(tt, S1[j - 1], s3, pf_params) *
                  qm[lj];

            if (sc_wrapper->pair)
              ppp *= sc_wrapper->pair(i, j, sc_wrapper);

            prmt += ppp;
          }
        }

        ii  = my_iindx[i];  /* ii-j=[i,j]     */
        tt  = vrna_get_ptype(jindx[l + 1] + i, ptype);
        tt  = rtype[tt];
        if (hc_eval(i, l + 1, i + 1, l, VRNA_DECOMP_PAIR_ML, hc_dat)) {
          prmt1 = probs[ii - (l + 1)] *
                  exp_E_MLstem(tt,
                               S1[l],
                               S1[i + 1],
                               pf_params) *
                  expMLclosing;

          if (sc_wrapper->pair)
            prmt1 *= sc_wrapper->pair(i, l + 1, sc_wrapper);
        }
      }

      prmt *= expMLclosing;

      ml_helpers->prml[i] = prmt;

      /* l+1 is unpaired */
      if (hc_eval(k, l + 1, k, l, VRNA_DECOMP_ML_ML, hc_dat)) {
        ppp = ml_helpers->prm_l1[i] *
              expMLbase[1];

        if (sc_wrapper->red_ml)
          ppp *= sc_wrapper->red_ml(k, l + 1, k, l, sc_wrapper);

        /* add contributions of MB loops where any unstructured domain starts at l+1 */
        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (l + u < n) {
              if (hc_eval(k, l + u, k, l, VRNA_DECOMP_ML_ML, hc_dat)) {
                temp = domains_up->exp_energy_cb(fc,
                                                 l + 1,
                                                 l + u,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data) *
                       ml_helpers->pmlu[u][i] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(k, l + u, k, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          ml_helpers->pmlu[0][i] = ppp + prmt1;
        }

        ml_helpers->prm_l[i] = ppp + prmt1;
      } else {
        /* skip configuration where l+1 is unpaired */
        ml_helpers->prm_l[i] = prmt1;

        if (with_ud)
          ml_helpers->pmlu[0][i] = prmt1;
      }

      if (hc_eval(i, l, i + 1, l, VRNA_DECOMP_ML_ML, hc_dat)) {
        ppp = prm_MLb *
              expMLbase[1];

        if (sc_wrapper->red_ml)
          ppp *= sc_wrapper->red_ml(i, l, i + 1, l, sc_wrapper);

        if (with_ud) {
          for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
            u = domains_up->uniq_motif_size[cnt];
            if (1 + u <= i) {
              if (hc_eval(i - u + 1, l, i + 1, l, VRNA_DECOMP_ML_ML, hc_dat)) {
                temp = domains_up->exp_energy_cb(fc,
                                                 i - u + 1,
                                                 i,
                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 domains_up->data) *
                       ml_helpers->prm_MLbu[u] *
                       expMLbase[u];

                if (sc_wrapper->red_ml)
                  temp *= sc_wrapper->red_ml(i - u + 1, l, i + 1, l, sc_wrapper);

                ppp += temp;
              }
            }
          }
          ml_helpers->prm_MLbu[0] = ppp + ml_helpers->prml[i];
        }

        prm_MLb = ppp + ml_helpers->prml[i];
        /* same as:    prm_MLb = 0;
         * for (i=1; i<=k-1; i++) prm_MLb += prml[i]*expMLbase[k-i-1]; */
      } else {
        /* skip all configurations where i is unpaired */
        prm_MLb = ml_helpers->prml[i];

        if (with_ud)
          ml_helpers->prm_MLbu[0] = ml_helpers->prml[i];
      }

      ml_helpers->prml[i] = ml_helpers->prml[i] + ml_helpers->prm_l[i];

      tt = ptype[jindx[l] + k];

      if (with_gquad) {
        if ((!tt) &&
            (G[kl] == 0.))
          continue;
      } else {
        if (qb[kl] == 0.)
          continue;
      }

      temp = prm_MLb;

      if (sn[k] == sn[k - 1]) {
        if (sc_wrapper->decomp_ml) {
          for (i = 1; i <= k - 2; i++)
            temp += ml_helpers->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)] *
                    sc_wrapper->decomp_ml(i + 1, l, k - 1, k, sc_wrapper);
        } else {
          for (i = 1; i <= k - 2; i++)
            temp += ml_helpers->prml[i] *
                    qm[my_iindx[i + 1] - (k - 1)];
        }
      }

      s5  = ((k > 1) && (sn[k] == sn[k - 1])) ? S1[k - 1] : -1;
      s3  = ((l < n) && (sn[l + 1] == sn[l])) ? S1[l + 1] : -1;

      if ((with_gquad) &&
          (qb[kl] == 0.)) {
        temp *= G[kl] *
                expMLstem;
      } else if (hc_eval(k, l, k, l, VRNA_DECOMP_ML_STEM, hc_dat)) {
        if (tt == 0)
          tt = 7;

        temp *= exp_E_MLstem(tt, s5, s3, pf_params);
      }

      if (sc_wrapper->red_stem)
        temp *= sc_wrapper->red_stem(k, l, k, l, sc_wrapper);

      probs[kl] += temp *
                   scale[2];

      if (probs[kl] > (*Qmax)) {
        (*Qmax) = probs[kl];
        if ((*Qmax) > max_real / 10.)
          printf("P close to overflow: %d %d %g %g\n",
                               k, l, probs[kl], qb[kl]);
      }

      if (probs[kl] >= max_real) {
        (*ov)++;
        probs[kl] = FLT_MAX;
      }

      /* rotate prm_MLbu entries required for unstructured domain feature */
      rotate_ml_helper_arrays_inner(ml_helpers);
    } /* end for (k=..) */
  }

  rotate_ml_helper_arrays_outer(ml_helpers);
}


PRIVATE FLT_OR_DBL
numerator_single(vrna_fold_compound_t *vc,
                 int                  i,
                 int                  j)
{
  return 1.;
}
/* calculate base pairing probs */
PRIVATE INLINE void
bppm_circ(vrna_fold_compound_t  *fc,
          constraints_helper    *constraints)
{
  char                      *ptype;
  unsigned char             *hard_constraints, eval;
  short                     *S, *S1, **SS, **S5, **S3;
  unsigned int              s, n_seq, type, rt, *tt, **a2s;
  int                       n, i, j, k, l, ij, *rtype, *my_iindx, *jindx;
  FLT_OR_DBL                tmp, tmp2, expMLclosing, *qb, *qm, *qm1, *probs, *scale, *expMLbase, qo;
  vrna_hc_t                 *hc;
  vrna_exp_param_t          *pf_params;
  vrna_mx_pf_t              *matrices;
  vrna_md_t                 *md;
  struct hc_mb_def_dat      *hc_dat_mb;
  vrna_hc_eval_f hc_eval_mb;
  struct sc_int_exp_dat     *sc_dat_int;
  struct sc_mb_exp_dat      *sc_dat_mb;

  FLT_OR_DBL                (*numerator_f)(vrna_fold_compound_t *fc,
                                           int                  i,
                                           int                  j);

  n                 = (int)fc->length;
  n_seq             = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  SS                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3                = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s               = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  pf_params         = fc->exp_params;
  md                = &(pf_params->model_details);
  S                 = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
  S1                = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  my_iindx          = fc->iindx;
  jindx             = fc->jindx;
  ptype             = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->ptype : NULL;
  type              = 0;
  rt                = 0;
  tt                = NULL;
  hc                = fc->hc;
  matrices          = fc->exp_matrices;
  qb                = matrices->qb;
  qm                = matrices->qm;
  qm1               = matrices->qm1;
  probs             = matrices->probs;
  scale             = matrices->scale;
  expMLbase         = matrices->expMLbase;
  qo                = matrices->qo;
  hard_constraints  = hc->mx;
  hc_dat_mb         = &(constraints->hc_dat_mb);
  hc_eval_mb        = constraints->hc_eval_mb;
  sc_dat_int        = &(constraints->sc_wrapper_int);
  sc_dat_mb         = &(constraints->sc_wrapper_mb);

  expMLclosing  = pf_params->expMLclosing;
  rtype         = &(pf_params->model_details.rtype[0]);

  switch (fc->type) {
    case  VRNA_FC_TYPE_SINGLE:
      numerator_f = numerator_single;
      tt          = NULL;
      break;
    default:
      numerator_f = NULL;
      break;
  }

  /* 1. exterior pair i,j */
  for (i = 1; i <= n; i++) {
    probs[my_iindx[i] - i] = 0;

    for (j = i + 1; j <= n; j++) {
      ij = my_iindx[i] - j;
      if (qb[ij] > 0.) {
        probs[ij] = numerator_f(fc, i, j) / qo;

        if (fc->type == VRNA_FC_TYPE_SINGLE) {
          type  = vrna_get_ptype_md(S[i], S[j], md);
          rt    = vrna_get_ptype_md(S[j], S[i], md);
        } else {
          for (s = 0; s < n_seq; s++)
            tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
        }

        /* 1.1. Exterior Hairpin Contribution */
        tmp2 = vrna_exp_E_hp_loop(fc, j, i);

        if (hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
          /*
           * 1.2. Exterior Interior Loop Contribution
           * 1.2.1. i,j  delimtis the "left" part of the interior loop
           * (j,i) is "outer pair"
           */
          for (k = 1; k < i - 1; k++) {
            int ln1, ln3, lstart;
            ln1 = k - 1;
            ln3 = n - j;

            if (hc->up_int[j + 1] < (ln1 + ln3))
              break;

            if ((ln1 + ln3) > MAXLOOP)
              break;

            lstart = (ln1 + ln3) + i - 1 - MAXLOOP;
            if (lstart < k + 1)
              lstart = k + 1;

            for (l = lstart; l < i; l++) {
              int ln2, ln2a, ln1a, type_2;
              ln2 = i - l - 1;

              if ((ln1 + ln2 + ln3) > MAXLOOP)
                continue;

              if (hc->up_int[l + 1] < ln2)
                continue;

              if (qb[my_iindx[k] - l] == 0.)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(k, l, i, j, VRNA_DECOMP_PAIR_IL, hc->data);

              if (eval) {
                tmp = qb[my_iindx[k] - l];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2 = vrna_get_ptype(jindx[l] + k, ptype);

                  tmp *= exp_E_IntLoop(ln1 + ln3,
                                       ln2,
                                       rt,
                                       rtype[type_2],
                                       S1[j + 1],
                                       S1[i - 1],
                                       S1[k - 1],
                                       S1[l + 1],
                                       pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    ln2a    = a2s[s][i - 1];
                    ln2a    -= a2s[s][l];
                    ln1a    = a2s[s][n] - a2s[s][j];
                    ln1a    += a2s[s][k - 1];
                    type_2  = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    tmp     *= exp_E_IntLoop(ln1a, ln2a, tt[s], type_2,
                                             S3[s][j],
                                             S5[s][i],
                                             S5[s][k],
                                             S3[s][l],
                                             pf_params);
                  }
                }

                if (sc_dat_int->pair_ext)
                  tmp *= sc_dat_int->pair_ext(k, l, i, j, sc_dat_int);

                tmp2 += tmp *
                        scale[ln1 + ln2 + ln3];
              }
            }
          }

          /* 1.2.2. i,j  delimtis the "right" part of the interior loop  */
          for (k = j + 1; k < n; k++) {
            int ln1, lstart;
            ln1 = k - j - 1;

            if (hc->up_int[j + 1] < ln1)
              break;

            if ((ln1 + i - 1) > MAXLOOP)
              break;

            lstart = ln1 + i - 1 + n - MAXLOOP;
            if (lstart < k + 1)
              lstart = k + 1;

            for (l = lstart; l <= n; l++) {
              int ln2, ln3, ln2a, ln1a, type_2;
              ln2 = i - 1;
              ln3 = n - l;

              if ((ln1 + ln2 + ln3) > MAXLOOP)
                continue;

              if (hc->up_int[l + 1] < (ln2 + ln3))
                continue;

              if (qb[my_iindx[k] - l] == 0.)
                continue;

              eval = (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 : 0;
              if (hc->f)
                eval = hc->f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc->data) ? eval : 0;

              if (eval) {
                tmp = qb[my_iindx[k] - l];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  type_2  = vrna_get_ptype(jindx[l] + k, ptype);
                  tmp     *= exp_E_IntLoop(ln2 + ln3,
                                           ln1,
                                           rtype[type_2],
                                           rt,
                                           S1[l + 1],
                                           S1[k - 1],
                                           S1[i - 1],
                                           S1[j + 1],
                                           pf_params);
                } else {
                  for (s = 0; s < n_seq; s++) {
                    ln1a    = a2s[s][k] - a2s[s][j + 1];
                    ln2a    = a2s[s][i - 1] + a2s[s][n] - a2s[s][l];
                    type_2  = vrna_get_ptype_md(SS[s][l], SS[s][k], md);
                    tmp     *= exp_E_IntLoop(ln2a, ln1a, type_2, tt[s],
                                             S3[s][l],
                                             S5[s][k],
                                             S5[s][i],
                                             S3[s][j],
                                             pf_params);
                  }
                }

                if (sc_dat_int->pair_ext)
                  tmp *= sc_dat_int->pair_ext(i, j, k, l, sc_dat_int);

                tmp2 += tmp *
                        scale[ln1 + ln2 + ln3];
              }
            }
          }
        }

        /* 1.3 Exterior multiloop decomposition */
        if (hc_eval_mb(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML_EXT, hc_dat_mb)) {
          FLT_OR_DBL sc_contrib = 1.;

          if (sc_dat_mb->pair_ext)
            sc_contrib = sc_dat_mb->pair_ext(i, j, sc_dat_mb);

          /* 1.3.1 Middle part                    */
          if ((i > 2) &&
              (j < n - 1)) {
            tmp = qm[my_iindx[1] - i + 1] *
                  qm[my_iindx[j + 1] - n];

            if (fc->type == VRNA_FC_TYPE_SINGLE) {
              tmp *= exp_E_MLstem(type,
                                  S1[i - 1],
                                  S1[j + 1],
                                  pf_params) *
                     expMLclosing;
            } else {
              for (s = 0; s < n_seq; s++)
                tmp *= exp_E_MLstem(rtype[tt[s]],
                                    S5[s][i],
                                    S3[s][j],
                                    pf_params);

              tmp *= pow(expMLclosing, n_seq);
            }

            tmp2 += tmp *
                    sc_contrib;
          }

          /* 1.3.2 Left part  */
          if (hc->up_ml[j + 1] >= (n - j)) {
            for (k = 2; k < i - 2; k++) {
              if ((hc_eval_mb(i, n, i, j, VRNA_DECOMP_ML_ML, hc_dat_mb)) &&
                  (hc_eval_mb(1, i - 1, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc_dat_mb))) {
                tmp = qm[my_iindx[1] - k] *
                      qm1[jindx[i - 1] + k + 1] *
                      expMLbase[n - j];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= exp_E_MLstem(type,
                                      S1[i - 1],
                                      S1[j + 1],
                                      pf_params) *
                         expMLclosing;
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= exp_E_MLstem(rtype[tt[s]],
                                        S5[s][i],
                                        S3[s][j],
                                        pf_params);

                  tmp *= pow(expMLclosing, n_seq);
                }

                if (sc_dat_mb->red_ml)
                  tmp *= sc_dat_mb->red_ml(i, n, i, j, sc_dat_mb);

                if (sc_dat_mb->decomp_ml)
                  tmp *= sc_dat_mb->decomp_ml(1, i - 1, k, k + 1, sc_dat_mb);

                tmp2 += tmp *
                        sc_contrib;
              }
            }
          }

          /* 1.3.3 Right part */
          if (hc->up_ml[1] >= (i - 1)) {
            for (k = j + 2; k < n - 1; k++) {
              if ((hc_eval_mb(1, j, i, j, VRNA_DECOMP_ML_ML, hc_dat_mb)) &&
                  (hc_eval_mb(j + 1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc_dat_mb))) {
                tmp = qm[my_iindx[j + 1] - k] *
                      qm1[jindx[n] + k + 1] *
                      expMLbase[i - 1];

                if (fc->type == VRNA_FC_TYPE_SINGLE) {
                  tmp *= exp_E_MLstem(type,
                                      S1[i - 1],
                                      S1[j + 1],
                                      pf_params) *
                         expMLclosing;
                } else {
                  for (s = 0; s < n_seq; s++)
                    tmp *= exp_E_MLstem(rtype[tt[s]],
                                        S5[s][i],
                                        S3[s][j],
                                        pf_params);

                  tmp *= pow(expMLclosing, n_seq);
                }

                if (sc_dat_mb->red_ml)
                  tmp *= sc_dat_mb->red_ml(1, j, i, j, sc_dat_mb);

                if (sc_dat_mb->decomp_ml)
                  tmp *= sc_dat_mb->decomp_ml(j + 1, n, k, k + 1, sc_dat_mb);

                tmp2 += tmp *
                        sc_contrib;
              }
            }
          }

          /* all exterior loop decompositions for pair i,j done  */
        }

        probs[ij] *= tmp2;
      } else {
        probs[ij] = 0;
      }
    }
  }

  free(tt);
}



PRIVATE FLT_OR_DBL
contrib_ext_pair(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 constraints_helper   *constraints)
{
  unsigned char     type;
  char              *ptype;
  short             *S1, s5, s3;
  unsigned int      *sn, n;
  int               *jindx;
  FLT_OR_DBL        contribution;
  vrna_exp_param_t  *pf_params;
  vrna_sc_t         *sc;

  n         = fc->length;
  pf_params = fc->exp_params;
  S1        = fc->sequence_encoding;
  sn        = fc->strand_number;
  ptype     = fc->ptype;
  jindx     = fc->jindx;
  sc        = fc->sc;

  type  = vrna_get_ptype(jindx[j] + i, ptype);
  s5    = ((i > 1) && (sn[i] == sn[i - 1])) ? S1[i - 1] : -1;
  s3    = ((j < n) && (sn[j + 1] == sn[j])) ? S1[j + 1] : -1;

  contribution = vrna_exp_E_ext_stem(type, s5, s3, pf_params);

  if ((sc) &&
      (sc->exp_f))
    contribution *= sc->exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, sc->data);

  return contribution;
}


PRIVATE void
compute_bpp_external(vrna_fold_compound_t *fc,
                     constraints_helper   *constraints)
{
  unsigned int              i, j, n;
  int                       *my_iindx, ij;
  FLT_OR_DBL                *probs, *q1k, *qln, *qb;
  vrna_mx_pf_t              *matrices;
  struct hc_ext_def_dat     *hc_dat;
  vrna_hc_eval_f evaluate;

  FLT_OR_DBL                (*contrib_f)(vrna_fold_compound_t *,
                                         unsigned int,
                                         unsigned int,
                                         constraints_helper *);

  n         = fc->length;
  my_iindx  = fc->iindx;
  matrices  = fc->exp_matrices;
  qb        = matrices->qb;
  probs     = matrices->probs;
  q1k       = matrices->q1k;
  qln       = matrices->qln;
  hc_dat    = &(constraints->hc_dat_ext);
  evaluate  = constraints->hc_eval_ext;

  contrib_f =  &contrib_ext_pair;

  for (i = 1; i <= n; i++) {
    for (j = i + 1; j <= n; j++) {
      ij        = my_iindx[i] - j;
      probs[ij] = 0.;

      if ((evaluate(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, hc_dat)) &&
          (qb[ij] > 0.)) {
        probs[ij] = q1k[i - 1] *
                    qln[j + 1] /
                    q1k[n];
        probs[ij] *= contrib_f(fc, i, j, constraints);
      }
    }
  }
}



PRIVATE void
multistrand_update_Y5(vrna_fold_compound_t  *fc,
                      int                   l,
                      FLT_OR_DBL            *Y5,
                      FLT_OR_DBL            **Y5p,
                      constraints_helper    *constraints)
{
  short             *S, *S1;
  unsigned int      i, s, *se, *sn, type, n, end;
  int               *my_iindx;
  FLT_OR_DBL        *probs, *q, *scale, qtmp;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;

  n         = fc->length;
  sn        = fc->strand_number;
  se        = fc->strand_end;
  my_iindx  = fc->iindx;
  q         = fc->exp_matrices->q;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  sc_wrapper  = &(constraints->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;

  /* compute Y5 for all strands */
  for (s = 0; s < fc->strands; s++) {
    unsigned int j;

    Y5[s] = 0;

    if ((se[s] < l) &&
        (sn[l] == sn[l + 1])) {
      /* pre-compute newly available Y5p[s][j] with j == l + 1 */
      end = se[s];
      j   = l + 1;

      Y5p[s][j] = 0.;

      if (probs[my_iindx[end] - j] > 0) {
        type  = vrna_get_ptype_md(S[j], S[end], md);
        qtmp  = probs[my_iindx[end] - j] *
                vrna_exp_E_ext_stem(type,
                                    S1[j - 1],
                                    -1,
                                    pf_params) *
                scale[2];

        if (sc_red_stem)
          qtmp *= sc_red_stem(j, end, j, end, sc_wrapper);

        Y5p[s][j] += qtmp;
      }

      for (i = 1; i < end; i++) {
        if ((probs[my_iindx[i] - j] > 0) &&
            (sn[i] == sn[i + 1])) {
          type  = vrna_get_ptype_md(S[j], S[i], md);
          qtmp  = probs[my_iindx[i] - j] *
                  vrna_exp_E_ext_stem(type,
                  S1[j - 1],
                  S1[i + 1],
                  pf_params) *
                  q[my_iindx[i + 1] - end] *
                  scale[2];

          if (sc_red_stem)
            qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);
          if (sc_split)
            qtmp *= sc_split(i, end, i + 1, sc_wrapper);

          Y5p[s][j] += qtmp;
        }
      }

      if ((probs[my_iindx[i] - j] > 0) &&
          (sn[i] == sn[i + 1])) {
        type  = vrna_get_ptype_md(S[j], S[i], md);
        qtmp  = probs[my_iindx[i] - j] *
                vrna_exp_E_ext_stem(type,
                                    S1[j - 1],
                                    S1[i + 1],
                                    pf_params) *
                scale[2];

        if (sc_red_stem)
          qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);

        Y5p[s][j] += qtmp;
      }

      /* recompute Y5[s] */
      Y5[s] += Y5p[s][l + 1];
      for (j = l + 2; j <= n; j++) {
        qtmp = q[my_iindx[l + 1] - (j - 1)] *
               Y5p[s][j];

        if (sc_split)
          qtmp *= sc_split(l + 1, j, j, sc_wrapper);

        Y5[s] += qtmp;
      }
    }
  }
}


PRIVATE void
multistrand_update_Y3(vrna_fold_compound_t  *fc,
                      int                   l,
                      FLT_OR_DBL            **Y3,
                      FLT_OR_DBL            **Y3p,
                      constraints_helper    *constraints)
{
  short             *S, *S1;
  unsigned int      i, j, k, s, n, start, type, *ss, *sn;
  int               *my_iindx;
  FLT_OR_DBL        *q, *probs, *scale, qtmp;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;

  n         = fc->length;
  sn        = fc->strand_number;
  ss        = fc->strand_start;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  my_iindx  = fc->iindx;
  q         = fc->exp_matrices->q;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  sc_wrapper  = &(constraints->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;

  for (s = 0; s < fc->strands; s++) {
    start = ss[s];
    if (start == l + 1) {
      for (i = 1; i < start; i++) {
        Y3p[s][i] = 0;

        if (sn[i] == sn[i + 1]) {
          if (probs[my_iindx[i] - start] > 0) {
            type = vrna_get_ptype_md(S[start], S[i], md);
            qtmp = probs[my_iindx[i] - start] *
                         vrna_exp_E_ext_stem(type,
                                             -1,
                                             S1[i + 1],
                                             pf_params) *
                         scale[2];

            if (sc_red_stem)
              qtmp *= sc_red_stem(start, i, start, i, sc_wrapper);

            Y3p[s][i] += qtmp;
          }

          for (j = start + 1; j <= n; j++) {
            if ((probs[my_iindx[i] - j] > 0) &&
                (sn[j - 1] == sn[j])) {
              type = vrna_get_ptype_md(S[j], S[i], md);
              qtmp = probs[my_iindx[i] - j] *
                           vrna_exp_E_ext_stem(type,
                                               S1[j - 1],
                                               S1[i + 1],
                                               pf_params) *
                           q[my_iindx[start] - (j - 1)] *
                           scale[2];

              if (sc_red_stem)
                qtmp *= sc_red_stem(j, i, j, i, sc_wrapper);
              if (sc_split)
                qtmp *= sc_split(start, j, j, sc_wrapper);

              Y3p[s][i] += qtmp;
            }
          }
        }
      }

      for (k = 1; k < start; k++) {
        Y3[s][k] = 0.;

        if (sn[k - 1] == sn[k]) {
          for (i = 1; i < k - 1; i++) {
            if (sn[i] == sn[i + 1]) {
              qtmp = q[my_iindx[i + 1] - (k - 1)] *
                     Y3p[s][i];

              if (sc_split)
                qtmp *= sc_split(i, k - 1, i + 1, sc_wrapper);

              Y3[s][k] += qtmp;
            }
          }

          Y3[s][k] += Y3p[s][k - 1];
        }
      }
    }
  }
}


PRIVATE void
multistrand_contrib(vrna_fold_compound_t  *fc,
                    int                   l,
                    FLT_OR_DBL            *Y5,
                    FLT_OR_DBL            **Y3,
                    constraints_helper    *constraints,
                    FLT_OR_DBL            *Qmax,
                    int                   *ov)
{
  short             *S, *S1, s5, s3;
  unsigned int      k, s, *sn, *se, *ss, end, start, type;
  int               *my_iindx, kl;
  FLT_OR_DBL        *q, *qb, *probs, tmp, qtmp;
  vrna_md_t         *md;
  vrna_exp_param_t  *pf_params;
  struct sc_ext_exp_dat   *sc_wrapper;
  sc_ext_exp_cb           sc_red_stem;
  sc_ext_exp_split        sc_split;

  sn        = fc->strand_number;
  ss        = fc->strand_start;
  se        = fc->strand_end;
  S         = fc->sequence_encoding2;
  S1        = fc->sequence_encoding;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  my_iindx  = fc->iindx;
  q         = fc->exp_matrices->q;
  qb        = fc->exp_matrices->qb;
  probs     = fc->exp_matrices->probs;
  sc_wrapper  = &(constraints->sc_wrapper_ext);
  sc_red_stem = sc_wrapper->red_stem;
  sc_split    = sc_wrapper->split;

  for (k = l - 1; k > 1; k--) {
    kl = my_iindx[k] - l;
    if (qb[kl] > 0) {
      tmp = 0.;
      for (s = 0; s < fc->strands; s++) {
        end   = se[s];
        start = ss[s];
        if (end == k - 1)
          tmp += Y5[s];
        else if ((end < k - 1) &&
                 (sn[k - 1] == sn[k]))
          tmp += Y5[s] *
                 q[my_iindx[end + 1] - (k - 1)];
        else if (start == l + 1)
          tmp += Y3[s][k];
        else if ((start > l + 1) &&
                 (sn[l] == sn[l + 1]))
          tmp += Y3[s][k] *
                 q[my_iindx[l + 1] - (start - 1)];
      }

      type      = vrna_get_ptype_md(S[k], S[l], md);
      s5        = (sn[k - 1] == sn[k]) ? S1[k - 1] : -1;
      s3        = (sn[l] == sn[l + 1]) ? S1[l + 1] : -1;
      qtmp      = vrna_exp_E_ext_stem(type,
                                      s5,
                                      s3,
                                      pf_params);

      if (sc_red_stem)
        qtmp *= sc_red_stem(k, l, k, l, sc_wrapper);

      probs[kl] += tmp * qtmp;
    }
  }
}


PRIVATE FLT_OR_DBL
exp_E_interior_loop(vrna_fold_compound_t  *fc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l)
{
  unsigned char         sliding_window, type, type2;
  char                  *ptype, **ptype_local;
  unsigned char         *hc_mx, **hc_mx_local, eval_loop, hc_decompose_ij, hc_decompose_kl;
  short                 *S1, **SS, **S5, **S3;
  unsigned int          n, *sn, n_seq, s, **a2s;
  int                   u1, u2, *rtype, *jindx, *hc_up;
  FLT_OR_DBL            qbt1, q_temp, *scale;
  vrna_exp_param_t      *pf_params;
  vrna_md_t             *md;
  vrna_ud_t             *domains_up;
  eval_hc               evaluate;
  struct hc_int_def_dat hc_dat_local;
  struct sc_int_exp_dat sc_wrapper;

  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  n               = fc->length;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? NULL : fc->ptype) : NULL;
  ptype_local     =
    (fc->type == VRNA_FC_TYPE_SINGLE) ? (sliding_window ? fc->ptype_local : NULL) : NULL;
  S1          = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : NULL;
  SS          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  S5          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
  S3          = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
  a2s         = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
  jindx       = fc->jindx;
  hc_mx       = (sliding_window) ? NULL : fc->hc->mx;
  hc_mx_local = (sliding_window) ? fc->hc->matrix_local : NULL;
  hc_up       = fc->hc->up_int;
  pf_params   = fc->exp_params;
  sn          = fc->strand_number;
  md          = &(pf_params->model_details);
  scale       = fc->exp_matrices->scale;
  domains_up  = fc->domains_up;
  rtype       = &(md->rtype[0]);
  qbt1        = 0.;
  u1          = k - i - 1;
  u2          = j - l - 1;

  if ((sn[k] != sn[i]) || (sn[j] != sn[l]))
    return qbt1;

  if (hc_up[l + 1] < u2)
    return qbt1;

  if (hc_up[i + 1] < u1)
    return qbt1;

  evaluate = prepare_hc_int_def(fc, &hc_dat_local);

  init_sc_int_exp(fc, &sc_wrapper);

  hc_decompose_ij = (sliding_window) ? hc_mx_local[i][j - i] : hc_mx[n * i + j];
  hc_decompose_kl = (sliding_window) ? hc_mx_local[k][l - k] : hc_mx[n * k + l];
  eval_loop       = ((hc_decompose_ij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
                     (hc_decompose_kl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) ?
                    1 : 0;

  /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
  if (eval_loop && evaluate(i, j, k, l, &hc_dat_local)) {
    q_temp = 0;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = (sliding_window) ?
               vrna_get_ptype_window(i, j, ptype_local) :
               vrna_get_ptype(jindx[j] + i, ptype);
        type2 = (sliding_window) ?
                rtype[vrna_get_ptype_window(k, l, ptype_local)] :
                rtype[vrna_get_ptype(jindx[l] + k, ptype)];

        q_temp = exp_E_IntLoop(u1,
                               u2,
                               type,
                               type2,
                               S1[i + 1],
                               S1[j - 1],
                               S1[k - 1],
                               S1[l + 1],
                               pf_params);

        break;
    }

    /* soft constraints */
    if (sc_wrapper.pair)
      q_temp *= sc_wrapper.pair(i, j, k, l, &sc_wrapper);

    qbt1 += q_temp *
            scale[u1 + u2 + 2];

    /* unstructured domains */
    if (domains_up && domains_up->exp_energy_cb) {
      FLT_OR_DBL qq5, qq3;

      qq5 = qq3 = 0.;

      if (u1 > 0) {
        qq5 = domains_up->exp_energy_cb(fc,
                                        i + 1, k - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      if (u2 > 0) {
        qq3 = domains_up->exp_energy_cb(fc,
                                        l + 1, j - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      qbt1 += q_temp *
              qq5 *
              scale[u1 + u2 + 2];      /* only motifs in 5' part */
      qbt1 += q_temp *
              qq3 *
              scale[u1 + u2 + 2];      /* only motifs in 3' part */
      qbt1 += q_temp *
              qq5 *
              qq3 *
              scale[u1 + u2 + 2]; /* motifs in both parts */
    }
  }

  free_sc_int_exp(&sc_wrapper);

  return qbt1;
}




PUBLIC FLT_OR_DBL
vrna_exp_E_interior_loop(vrna_fold_compound_t *fc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l)
{
  if (fc)
    return exp_E_interior_loop(fc, i, j, k, l);

  return 0.;
}



PUBLIC char
vrna_bpp_symbol(const float *x)
{
  /*  if( ((x[1]-x[2])*(x[1]-x[2]))<0.1&&x[0]<=0.677) return '|'; */
  if (x[0] > 0.667)
    return '.';

  if (x[1] > 0.667)
    return '(';

  if (x[2] > 0.667)
    return ')';

  if ((x[1] + x[2]) > x[0]) {
    if ((x[1] / (x[1] + x[2])) > 0.667)
      return '{';

    if ((x[2] / (x[1] + x[2])) > 0.667)
      return '}';
    else
      return '|';
  }

  if (x[0] > (x[1] + x[2]))
    return ',';

  return ':';
}

PUBLIC char *
vrna_db_from_probs(const FLT_OR_DBL *p,
                   unsigned int     length)
{
  int   i, j, *index;
  float P[3];    /* P[][0] unpaired, P[][1] upstream p, P[][2] downstream p */
  char  *s;

  s = NULL;

  if (p) {
    index = vrna_idx_row_wise(length);
    s     = (char *)vrna_alloc(sizeof(char) * (length + 1));

    for (j = 1; j <= length; j++) {
      P[0]  = 1.0;
      P[1]  = P[2] = 0.0;
      for (i = 1; i < j; i++) {
        P[2]  += (float)p[index[i] - j];  /* j is paired downstream */
        P[0]  -= (float)p[index[i] - j];  /* j is unpaired */
      }
      for (i = j + 1; i <= length; i++) {
        P[1]  += (float)p[index[j] - i];  /* j is paired upstream */
        P[0]  -= (float)p[index[j] - i];  /* j is unpaired */
      }
      s[j - 1] = vrna_bpp_symbol(P);
    }
    s[length] = '\0';
    free(index);
  }

  return s;
}


PRIVATE void
free_ml_helper_arrays(helper_arrays *ml_helpers)
{
  unsigned int u;

  free(ml_helpers->prm_l);
  free(ml_helpers->prm_l1);
  free(ml_helpers->prml);

  if (ml_helpers->pmlu) {
    for (u = 0; u <= ml_helpers->ud_max_size; u++)
      free(ml_helpers->pmlu[u]);
    free(ml_helpers->pmlu);
  }

  free(ml_helpers->prm_MLbu);
  free(ml_helpers);
}

PRIVATE void
free_constraints_helper(constraints_helper *helper)
{
  free_sc_ext_exp(&(helper->sc_wrapper_ext));
  free_sc_hp_exp(&(helper->sc_wrapper_hp));
  free_sc_int_exp(&(helper->sc_wrapper_int));
  free_sc_mb_exp(&(helper->sc_wrapper_mb));

  free(helper);
}


/* calculate base pairing probs  changed structure -------start--------------*/
PRIVATE int
pf_create_bppm(vrna_fold_compound_t *vc,
               char                 *structure)
{
  unsigned int      s;
  int               n, i, j, l, ij, *pscore, *jindx, ov = 0;
  FLT_OR_DBL        Qmax = 0;
  FLT_OR_DBL        *qb, *G, *probs;
  FLT_OR_DBL        *q1k, *qln;

  int               with_gquad;
  vrna_hc_t         *hc;
  vrna_sc_t         *sc;
  int               *my_iindx;
  int               circular, with_ud, with_ud_outside;
  vrna_exp_param_t  *pf_params;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;

  n           = vc->length;
  pscore      = NULL;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  circular    = md->circ;
  with_gquad  = md->gquad;

  hc  = vc->hc;
  sc  = vc->sc;

  domains_up  = vc->domains_up;
  matrices    = vc->exp_matrices;

  qb    = matrices->qb;
  G     = matrices->G;
  probs = matrices->probs;
  /* probs is 0.0000 now */
  q1k   = matrices->q1k;
  qln   = matrices->qln;

  with_ud         = (domains_up && domains_up->exp_energy_cb) ? 1 : 0;
  with_ud_outside = (with_ud && domains_up->probs_add) ? 1 : 0;

  /*
   * the following is a crude check whether the partition function forward recursion
   * has already been taken place
   */
  if ((qb) &&
      (probs) &&
      (circular ? matrices->qm2 != NULL : (q1k != NULL && qln != NULL))) {
    with_gquad = pf_params->model_details.gquad;
    /* probs is 0.0000 now */
    // printf("probs is :%f", probs);

    double      kTn = pf_params->kT / 10.;               /* kT in cal/mol  */
    int         corr_size = 5;
    int         corr_cnt = 0;
    vrna_ep_t   *bp_correction = vrna_alloc(sizeof(vrna_ep_t) * corr_size);
    FLT_OR_DBL  *Y5, **Y5p, **Y3, **Y3p;

    Y5  = NULL;
    Y5p = NULL;
    Y3  = NULL;
    Y3p = NULL;

    helper_arrays       *ml_helpers;
    constraints_helper  *constraints;

    ml_helpers  = get_ml_helper_arrays(vc);
    constraints = get_constraints_helper(vc);

    /* function pointer here igrone it */
    void                (*compute_bpp_int)(vrna_fold_compound_t *fc,
                                           int                  l,
                                           vrna_ep_t            **bp_correction,
                                           int                  *corr_cnt,
                                           int                  *corr_size,
                                           FLT_OR_DBL           *Qmax,
                                           int                  *ov,
                                           constraints_helper   *constraints);

    void (*compute_bpp_mul)(vrna_fold_compound_t  *fc,
                            int                   l,
                            helper_arrays         *ml_helpers,
                            FLT_OR_DBL            *Qmax,
                            int                   *ov,
                            constraints_helper    *constraints);

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      compute_bpp_int = &compute_bpp_internal;
      compute_bpp_mul = &compute_bpp_multibranch;
    }

    Qmax = 0;

    if (vc->strands > 1) {
      Y5  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * vc->strands);
      Y5p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * vc->strands);
      for (s = 0; s < vc->strands; s++)
        Y5p[s] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));

      Y3  = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * vc->strands);
      Y3p = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * vc->strands);
      for (s = 0; s < vc->strands; s++) {
        Y3[s]   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
        Y3p[s]  = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
      }
    }

    /* init diagonal entries unable to pair in pr matrix */
    for (i = 1; i <= n; i++)
      probs[my_iindx[i] - i] = 0.;
      // printf("\nprobs is :%f\n", probs);
    /* 1. external loop pairs, i.e. pairs not enclosed by any other pair (or external loop for circular RNAs) */
    // printf("into cir %d\n", circular);
    if (circular)
      bppm_circ(vc, constraints);
    else
      compute_bpp_external(vc, constraints);

    /* 2. all cases where base pair (k,l) is enclosed by another pair (i,j) */
    l = n;
    compute_bpp_int(vc,
                    l,
                    &bp_correction,
                    &corr_cnt,
                    &corr_size,
                    &Qmax,
                    &ov,
                    constraints);

    for (l = n - 1; l > 1; l--) {
      compute_bpp_int(vc,
                      l,
                      &bp_correction,
                      &corr_cnt,
                      &corr_size,
                      &Qmax,
                      &ov,
                      constraints);

      compute_bpp_mul(vc,
                      l,
                      ml_helpers,
                      &Qmax,
                      &ov,
                      constraints);

      if (vc->strands > 1) {
        multistrand_update_Y5(vc, l, Y5, Y5p, constraints);
        multistrand_update_Y3(vc, l, Y3, Y3p, constraints);
        multistrand_contrib(vc,
                            l,
                            Y5,
                            Y3,
                            constraints,
                            &Qmax,
                            &ov);
      }
    }
    // linked vc 
    if (vc->type == VRNA_FC_TYPE_SINGLE) {


      if ((sc) &&
          (sc->f) &&
          (sc->bt)) {
        for (i = 1; i <= n; i++)
          for (j = i + 1; j <= n; j++) {
            ij = my_iindx[i] - j;
            /*  search for possible auxiliary base pairs in hairpin loop motifs to store
             *  the corresponding probability corrections
             */
            if (hc->mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = sc->bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
              if (aux_bps) {
                FLT_OR_DBL qhp = vrna_exp_E_hp_loop(vc, i, j);
                /*  kind  different*/
                for (ptr = aux_bps; (ptr) && (ptr->i != 0); ptr++) {
                  bp_correction[corr_cnt].i   = ptr->i;
                  bp_correction[corr_cnt].j   = ptr->j;
                  bp_correction[corr_cnt++].p = probs[ij] * qhp;
                  if (corr_cnt == corr_size) {
                    corr_size     += 5;
                    bp_correction = vrna_realloc(bp_correction, sizeof(vrna_ep_t) * corr_size);
                  }
                }
              }

              free(aux_bps);
            }
          }

        /*  correct pairing probabilities for auxiliary base pairs from hairpin-, or interior loop motifs
         *  as augmented by the generalized soft constraints feature
         */
        // linked bp_correction  corr_cnt corr_size qhp aux_bps ptr sc->bt sc->data  hc->mx sc->f sc->bt
        for (i = 0; i < corr_cnt; i++) {
          ij = my_iindx[bp_correction[i].i] - bp_correction[i].j;
          /* printf("correcting pair %d, %d by %f\n", bp_correction[i].i, bp_correction[i].j, bp_correction[i].p); */
          probs[ij] += bp_correction[i].p / qb[ij];
        }
      }
    }

    for (i = 1; i <= n; i++)
      for (j = i + 1; j <= n; j++) {
        ij = my_iindx[i] - j;

        if (with_gquad) {
          if (qb[ij] > 0.) {
            probs[ij] *= qb[ij];
          } else if (G[ij] > 0.) {
            probs[ij] += q1k[i - 1] *
                         G[ij] *
                         qln[j + 1] /
                         q1k[n];
          }
        } else {
          if (qb[ij] > 0.) {
            probs[ij] *= qb[ij];
          }
        }
        // printf("%f,", probs[ij]);
      }
    if (structure != NULL) {
      /* s generate here -------------------*/
      // printf("---%d\n",n); 
      char *s = vrna_db_from_probs(probs, (unsigned int)n);
      /* strcuture here  generate by s  ---------------- */
      memcpy(structure, s, n);
      structure[n] = '\0';
      // printf("%s\n\n",structure);
      free(s);
    }

    if (ov > 0)
      printf("%d overflows occurred while backtracking;\n"
                           "you might try a smaller pf_scale than %g\n",
                           ov, pf_params->pf_scale);

    /* clean up */
    free_ml_helper_arrays(ml_helpers);

    free_constraints_helper(constraints);

    free(bp_correction);
    free(Y5);
    if (Y5p)
      for (unsigned int s = 0; s < vc->strands; s++)
        free(Y5p[s]);

    free(Y5p);

    if (Y3)
      for (unsigned int s = 0; s < vc->strands; s++)
        free(Y3[s]);

    free(Y3);

    if (Y3p)
      for (unsigned int s = 0; s < vc->strands; s++)
        free(Y3p[s]);

    free(Y3p);
  } /* end if 'check for forward recursion' */
  else {
    printf("bppm calculations have to be done after calling forward recursion");
    return 0;
  }

  return 1;
}




PUBLIC int
vrna_pairing_probs(vrna_fold_compound_t *vc,
                   char                 *structure)
{
  if (vc)
    return pf_create_bppm(vc, structure);

  return 0;
}

#define VRNA_STATUS_PF_POST (unsigned char)4
#define FLT_MIN __FLT_MIN__




PUBLIC void
vrna_ptypes_prepare(vrna_fold_compound_t  *fc,
                    unsigned int          options)
{
  if (!fc)
    return;

  if (options & VRNA_OPTION_MFE) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (options & VRNA_OPTION_WINDOW) {
          fc->ptype_local =
            (char **)vrna_realloc(fc->ptype_local, sizeof(char *) * (fc->length + 1));
        } else {
          if (!fc->ptype) {
            /* temporary hack for multi-strand case */
            if (fc->strands > 1) {
              int min_loop_size = fc->params->model_details.min_loop_size;
              fc->params->model_details.min_loop_size = 0;
              fc->ptype = vrna_ptypes(fc->sequence_encoding2,
                                      &(fc->params->model_details));
              fc->params->model_details.min_loop_size = min_loop_size;
            } else {
              fc->ptype = vrna_ptypes(fc->sequence_encoding2,
                                      &(fc->params->model_details));
            }
          }
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        break;

      default:
        break;
    }
  }

  if (options & VRNA_OPTION_PF) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (options & VRNA_OPTION_WINDOW) {
          fc->ptype_local =
            (char **)vrna_realloc(fc->ptype_local, sizeof(char *) * (fc->length + 1));
        } else {
          if (!fc->ptype) {
            /* temporary hack for multi-strand case */
            if (fc->strands > 1) {
              int min_loop_size = fc->exp_params->model_details.min_loop_size;
              fc->exp_params->model_details.min_loop_size = 0;
              fc->ptype = vrna_ptypes(fc->sequence_encoding2,
                                      &(fc->exp_params->model_details));
              fc->exp_params->model_details.min_loop_size = min_loop_size;
            } else {
              fc->ptype = vrna_ptypes(fc->sequence_encoding2, &(fc->exp_params->model_details));
            }
          }
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
          /* backward compatibility ptypes */
          if (!fc->ptype_pf_compat)
            fc->ptype_pf_compat = get_ptypes(fc->sequence_encoding2,
                                             &(fc->exp_params->model_details),
                                             1);

#endif
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        break;

      default:
        break;
    }
  }
}



PUBLIC void
vrna_hc_init_window(vrna_fold_compound_t *vc)
{
  unsigned int  n;
  vrna_hc_t     *hc;

  n = vc->length;

  /* free previous hard constraints */
  vrna_hc_free(vc->hc);

  /* allocate memory new hard constraints data structure */
  hc                = (vrna_hc_t *)vrna_alloc(sizeof(vrna_hc_t));
  hc->type          = VRNA_HC_WINDOW;
  hc->n             = n;
  hc->matrix_local  = (unsigned char **)vrna_alloc(sizeof(unsigned char *) * (n + 2));
  hc->up_ext        = NULL;
  hc->up_hp         = NULL;
  hc->up_int        = NULL;
  hc->up_ml         = NULL;
  hc->depot         = NULL;
  hc->state         = STATE_UNINITIALIZED;

  /* set new hard constraints */
  vc->hc = hc;

  /* add null pointers for the generalized hard constraint feature */
  hc->f         = NULL;
  hc->data      = NULL;
  hc->free_data = NULL;
}


PRIVATE void
prepare_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned char   option, type, t1, t2;
  unsigned int    i, j, k, n, s, *ss;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc = fc->hc;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n     = fc->length;
    ss    = fc->strand_start;
    depot = hc->depot;

    /* 2. apply constraints as stored in depot */
    if ((depot) && (depot->up)) {
      for (s = 0; s < depot->strands; s++) {
        for (k = 1; k <= depot->up_size[s]; k++) {
          /* process nucleotide-specific constraint */
          option  = depot->up[s][k].context;
          i       = ss[s] + k - 1; /* constraint position in current strand order */

          if (depot->up[s][k].nonspec) {
            /* this is actually a must-pair constraint */
            type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            /* acknowledge pairing direction */
            t1  = (depot->up[s][k].direction <= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;
            t2  = (depot->up[s][k].direction >= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;

            if (option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
              /* only allow for possibly non-canonical pairs, do not enforce them */
              for (j = 1; j < i; j++) {
                hc->mx[n * i + j] |= t1;
                hc->mx[n * j + i] |= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->mx[n * i + j] |= t2;
                hc->mx[n * j + i] |= t2;
              }
            } else {
              /* force pairing direction */
              for (j = 1; j < i; j++) {
                hc->mx[n * i + j] &= t1;
                hc->mx[n * j + i] &= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->mx[n * i + j] &= t2;
                hc->mx[n * j + i] &= t2;
              }
            }

            /* nucleotide mustn't be unpaired */
            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE)
              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
          } else {
            /* 'regular' nucleotide-specific constraint */
            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /*
               * force nucleotide to appear unpaired within a certain type of loop
               * do not allow i to be paired with any other nucleotide
               */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }

              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              hc->mx[n * i + i] = type;
            } else {
              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              /* do not allow i to be paired with any other nucleotide (in context type) */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
              }

              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            }
          }
        }
      }
    }
  }
}


PRIVATE void
prepare_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned char   option;
  unsigned int    i, j, k, p, q, n, actual_i, actual_j, strand_j, s, *ss;
  int             *idx;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc    = fc->hc;
  depot = hc->depot;
  ss    = fc->strand_start;

  if ((!depot) || (!depot->bp))
    return;

  if (options & VRNA_OPTION_WINDOW) {
  } else {
    n   = fc->length;
    idx = fc->jindx;

    /* 2. apply constraints stored in depot */
    for (s = 0; s < depot->strands; s++) {
      for (actual_i = 1; actual_i <= depot->bp_size[s]; actual_i++) {
        for (k = 0; k < depot->bp[s][actual_i].list_size; k++) {
          option    = depot->bp[s][actual_i].context[k];
          strand_j  = depot->bp[s][actual_i].strand_j[k];
          actual_j  = depot->bp[s][actual_i].j[k];
          i         = ss[s] + actual_i - 1;         /* constraint position in current strand order */
          j         = ss[strand_j] + actual_j - 1;  /* constraint position in current strand order */

          if (i < j) {
            /* apply the constraint */
            hc->mx[n * i + j] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            hc->mx[n * j + i] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

            /* is the ptype reset actually required??? */
            if ((fc->type == VRNA_FC_TYPE_SINGLE) &&
                (option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)) {
              /* reset ptype in case (i,j) is a non-canonical pair */
              if (fc->ptype[idx[j] + i] == 0)
                fc->ptype[idx[j] + i] = 7;
            }

            if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
              /*
               * remove all conflicting base pairs, i.e. do not allow i,j to pair
               * with any other nucleotide k
               */
              for (p = 1; p < i; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

                for (q = i + 1; q < j; q++) {
                  hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
              for (p = i + 1; p < j; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

                for (q = j + 1; q <= n; q++) {
                  hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
              for (p = j + 1; p <= n; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
              }
            }

            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /* do not allow i,j to be unpaired */
              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * j + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          }
        }
      }
    }
  }
}

#define STATE_CLEAN         (unsigned char)0
#define STATE_DIRTY_UP      (unsigned char)1
#define STATE_DIRTY_BP      (unsigned char)2
#define STATE_UNINITIALIZED (unsigned char)4

PUBLIC int
vrna_hc_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options)
{
  int ret = 0;


  if (fc) {
    if (options & VRNA_OPTION_WINDOW) {
      /* check for minimal hard constraints structure */
      if ((!fc->hc) || (fc->hc->type != VRNA_HC_WINDOW) || (!fc->hc->matrix_local))
        vrna_hc_init_window(fc);
    } else {
      if (fc->hc->state & STATE_UNINITIALIZED) {
        default_hc_up(fc, options);
        default_hc_bp(fc, options);
      }

      if (fc->hc->state & STATE_DIRTY_UP)
        prepare_hc_up(fc, options);

      if (fc->hc->state & STATE_DIRTY_BP)
        prepare_hc_bp(fc, options);

      if (fc->hc->state & ~STATE_CLEAN)
        hc_update_up(fc);
    }

    fc->hc->state = STATE_CLEAN;
    ret           = 1;
  }

  return ret;
}


PRIVATE unsigned int
get_mx_alloc_vector(vrna_fold_compound_t  *fc,
                    vrna_mx_type_e        mx_type,
                    unsigned int          options)
{
  unsigned int  v;
  vrna_md_t     *md_p;

  md_p = &(fc->params->model_details);

  v = ALLOC_NOTHING;

  /* default MFE matrices ? */
  if (options & VRNA_OPTION_MFE)
    v |= (mx_type == VRNA_MX_WINDOW) ? ALLOC_MFE_LOCAL : ALLOC_MFE_DEFAULT;

  /* default PF matrices ? */
  if (options & VRNA_OPTION_PF)
    v |= (md_p->compute_bpp) ? ALLOC_PF_DEFAULT : ALLOC_PF_WO_PROBS;

  if ((fc->strands > 1) || (options & VRNA_OPTION_HYBRID))
    v |= ALLOC_MULTISTRAND;

  /* matrices for circular folding ? */
  if (md_p->circ) {
    md_p->uniq_ML = 1; /* we need unique ML arrays for circular folding */
    v             |= ALLOC_CIRC;
  }

  /* unique ML decomposition ? */
  if (md_p->uniq_ML)
    v |= ALLOC_UNIQ;

  return v;
}

PRIVATE unsigned int
get_mx_mfe_alloc_vector_current(vrna_mx_mfe_t   *mx,
                                vrna_mx_type_e  mx_type)
{
  unsigned int mx_alloc_vector = ALLOC_NOTHING;

  if (mx) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        if (mx->f5)
          mx_alloc_vector |= ALLOC_F5;

        if (mx->f3)
          mx_alloc_vector |= ALLOC_F3;

        if ((mx->fms5) ||
            (mx->fms3))
          mx_alloc_vector |= ALLOC_MULTISTRAND;

        if (mx->c)
          mx_alloc_vector |= ALLOC_C;

        if (mx->fML)
          mx_alloc_vector |= ALLOC_FML;

        if (mx->fM1)
          mx_alloc_vector |= ALLOC_UNIQ;

        if (mx->fM2)
          mx_alloc_vector |= ALLOC_CIRC;

        break;

      default:
        break;
    }
  }

  return mx_alloc_vector;
}


PRIVATE void
mfe_matrices_free_default(vrna_mx_mfe_t *self)
{
  free(self->f5);
  free(self->f3);

  if (self->fms5)
    for (unsigned int s = 0; s < self->strands; s++)
      free(self->fms5[s]);

  free(self->fms5);

  if (self->fms3)
    for (unsigned int s = 0; s < self->strands; s++)
      free(self->fms3[s]);

  free(self->fms3);

  free(self->c);
  free(self->fML);
  free(self->fM1);
  free(self->fM2);
  free(self->ggg);
}

PRIVATE void
mfe_matrices_free_window(vrna_mx_mfe_t  *self,
                         unsigned int   length,
                         unsigned int   window_size)
{
  free(self->c_local);
  free(self->fML_local);
  free(self->ggg_local);
  free(self->f3_local);
}


PRIVATE void
mfe_matrices_free_2Dfold(vrna_mx_mfe_t  *self,
                         unsigned int   length,
                         int            turn,
                         int            *indx)
{
  unsigned int  i, j, ij;
  int           cnt1;

  /* This will be some fun... */
#ifdef COUNT_STATES
  if (self->N_F5 != NULL) {
    for (i = 1; i <= length; i++) {
      if (!self->N_F5[i])
        continue;

      for (cnt1 = self->k_min_F5[i]; cnt1 <= vars->k_max_F5[i]; cnt1++)
        if (vars->l_min_F5[i][cnt1] < INF) {
          vars->N_F5[i][cnt1] += vars->l_min_F5[i][cnt1] / 2;
          free(vars->N_F5[i][cnt1]);
        }

      if (vars->k_min_F5[i] < INF) {
        vars->N_F5[i] += vars->k_min_F5[i];
        free(vars->N_F5[i]);
      }
    }
    free(vars->N_F5);
  }

#endif

  if (self->E_F5 != NULL) {
    for (i = 1; i <= length; i++) {
      if (!self->E_F5[i])
        continue;

      for (cnt1 = self->k_min_F5[i]; cnt1 <= self->k_max_F5[i]; cnt1++)
        if (self->l_min_F5[i][cnt1] < INF) {
          self->E_F5[i][cnt1] += self->l_min_F5[i][cnt1] / 2;
          free(self->E_F5[i][cnt1]);
        }

      if (self->k_min_F5[i] < INF) {
        self->E_F5[i] += self->k_min_F5[i];
        free(self->E_F5[i]);
        self->l_min_F5[i] += self->k_min_F5[i];
        self->l_max_F5[i] += self->k_min_F5[i];
        free(self->l_min_F5[i]);
        free(self->l_max_F5[i]);
      }
    }
    free(self->E_F5);
    free(self->l_min_F5);
    free(self->l_max_F5);
    free(self->k_min_F5);
    free(self->k_max_F5);
  }

  if (self->E_F3 != NULL) {
    for (i = 1; i <= length; i++) {
      if (!self->E_F3[i])
        continue;

      for (cnt1 = self->k_min_F3[i]; cnt1 <= self->k_max_F3[i]; cnt1++)
        if (self->l_min_F3[i][cnt1] < INF) {
          self->E_F3[i][cnt1] += self->l_min_F3[i][cnt1] / 2;
          free(self->E_F3[i][cnt1]);
        }

      if (self->k_min_F3[i] < INF) {
        self->E_F3[i] += self->k_min_F3[i];
        free(self->E_F3[i]);
        self->l_min_F3[i] += self->k_min_F3[i];
        self->l_max_F3[i] += self->k_min_F3[i];
        free(self->l_min_F3[i]);
        free(self->l_max_F3[i]);
      }
    }
    free(self->E_F3);
    free(self->l_min_F3);
    free(self->l_max_F3);
    free(self->k_min_F3);
    free(self->k_max_F3);
  }

#ifdef COUNT_STATES
  if (self->N_C != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->N_C[ij])
          continue;

        for (cnt1 = self->k_min_C[ij]; cnt1 <= self->k_max_C[ij]; cnt1++)
          if (self->l_min_C[ij][cnt1] < INF) {
            self->N_C[ij][cnt1] += self->l_min_C[ij][cnt1] / 2;
            free(self->N_C[ij][cnt1]);
          }

        if (self->k_min_C[ij] < INF) {
          self->N_C[ij] += self->k_min_C[ij];
          free(self->N_C[ij]);
        }
      }
    }
    free(self->N_C);
  }

#endif

  if (self->E_C != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->E_C[ij])
          continue;

        for (cnt1 = self->k_min_C[ij]; cnt1 <= self->k_max_C[ij]; cnt1++)
          if (self->l_min_C[ij][cnt1] < INF) {
            self->E_C[ij][cnt1] += self->l_min_C[ij][cnt1] / 2;
            free(self->E_C[ij][cnt1]);
          }

        if (self->k_min_C[ij] < INF) {
          self->E_C[ij] += self->k_min_C[ij];
          free(self->E_C[ij]);
          self->l_min_C[ij] += self->k_min_C[ij];
          self->l_max_C[ij] += self->k_min_C[ij];
          free(self->l_min_C[ij]);
          free(self->l_max_C[ij]);
        }
      }
    }
    free(self->E_C);
    free(self->l_min_C);
    free(self->l_max_C);
    free(self->k_min_C);
    free(self->k_max_C);
  }

#ifdef COUNT_STATES
  if (self->N_M != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->N_M[ij])
          continue;

        for (cnt1 = self->k_min_M[ij]; cnt1 <= self->k_max_M[ij]; cnt1++)
          if (self->l_min_M[ij][cnt1] < INF) {
            self->N_M[ij][cnt1] += self->l_min_M[ij][cnt1] / 2;
            free(self->N_M[ij][cnt1]);
          }

        if (self->k_min_M[ij] < INF) {
          self->N_M[ij] += self->k_min_M[ij];
          free(self->N_M[ij]);
        }
      }
    }
    free(self->N_M);
  }

#endif

  if (self->E_M != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->E_M[ij])
          continue;

        for (cnt1 = self->k_min_M[ij]; cnt1 <= self->k_max_M[ij]; cnt1++)
          if (self->l_min_M[ij][cnt1] < INF) {
            self->E_M[ij][cnt1] += self->l_min_M[ij][cnt1] / 2;
            free(self->E_M[ij][cnt1]);
          }

        if (self->k_min_M[ij] < INF) {
          self->E_M[ij] += self->k_min_M[ij];
          free(self->E_M[ij]);
          self->l_min_M[ij] += self->k_min_M[ij];
          self->l_max_M[ij] += self->k_min_M[ij];
          free(self->l_min_M[ij]);
          free(self->l_max_M[ij]);
        }
      }
    }
    free(self->E_M);
    free(self->l_min_M);
    free(self->l_max_M);
    free(self->k_min_M);
    free(self->k_max_M);
  }

#ifdef COUNT_STATES
  if (self->N_M1 != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->N_M1[ij])
          continue;

        for (cnt1 = self->k_min_M1[ij]; cnt1 <= self->k_max_M1[ij]; cnt1++)
          if (self->l_min_M1[ij][cnt1] < INF) {
            self->N_M1[ij][cnt1] += self->l_min_M1[ij][cnt1] / 2;
            free(self->N_M1[ij][cnt1]);
          }

        if (self->k_min_M1[ij] < INF) {
          self->N_M1[ij] += self->k_min_M1[ij];
          free(self->N_M1[ij]);
        }
      }
    }
    free(self->N_M1);
  }

#endif

  if (self->E_M1 != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->E_M1[ij])
          continue;

        for (cnt1 = self->k_min_M1[ij]; cnt1 <= self->k_max_M1[ij]; cnt1++)
          if (self->l_min_M1[ij][cnt1] < INF) {
            self->E_M1[ij][cnt1] += self->l_min_M1[ij][cnt1] / 2;
            free(self->E_M1[ij][cnt1]);
          }

        if (self->k_min_M1[ij] < INF) {
          self->E_M1[ij] += self->k_min_M1[ij];
          free(self->E_M1[ij]);
          self->l_min_M1[ij]  += self->k_min_M1[ij];
          self->l_max_M1[ij]  += self->k_min_M1[ij];
          free(self->l_min_M1[ij]);
          free(self->l_max_M1[ij]);
        }
      }
    }
    free(self->E_M1);
    free(self->l_min_M1);
    free(self->l_max_M1);
    free(self->k_min_M1);
    free(self->k_max_M1);
  }

  if (self->E_M2 != NULL) {
    for (i = 1; i < length - turn - 1; i++) {
      if (!self->E_M2[i])
        continue;

      for (cnt1 = self->k_min_M2[i]; cnt1 <= self->k_max_M2[i]; cnt1++)
        if (self->l_min_M2[i][cnt1] < INF) {
          self->E_M2[i][cnt1] += self->l_min_M2[i][cnt1] / 2;
          free(self->E_M2[i][cnt1]);
        }

      if (self->k_min_M2[i] < INF) {
        self->E_M2[i] += self->k_min_M2[i];
        free(self->E_M2[i]);
        self->l_min_M2[i] += self->k_min_M2[i];
        self->l_max_M2[i] += self->k_min_M2[i];
        free(self->l_min_M2[i]);
        free(self->l_max_M2[i]);
      }
    }
    free(self->E_M2);
    free(self->l_min_M2);
    free(self->l_max_M2);
    free(self->k_min_M2);
    free(self->k_max_M2);
  }

  if (self->E_Fc != NULL) {
    for (cnt1 = self->k_min_Fc; cnt1 <= self->k_max_Fc; cnt1++)
      if (self->l_min_Fc[cnt1] < INF) {
        self->E_Fc[cnt1] += self->l_min_Fc[cnt1] / 2;
        free(self->E_Fc[cnt1]);
      }

    if (self->k_min_Fc < INF) {
      self->E_Fc += self->k_min_Fc;
      free(self->E_Fc);
      self->l_min_Fc  += self->k_min_Fc;
      self->l_max_Fc  += self->k_min_Fc;
      free(self->l_min_Fc);
      free(self->l_max_Fc);
    }
  }

  if (self->E_FcI != NULL) {
    for (cnt1 = self->k_min_FcI; cnt1 <= self->k_max_FcI; cnt1++)
      if (self->l_min_FcI[cnt1] < INF) {
        self->E_FcI[cnt1] += self->l_min_FcI[cnt1] / 2;
        free(self->E_FcI[cnt1]);
      }

    if (self->k_min_FcI < INF) {
      self->E_FcI += self->k_min_FcI;
      free(self->E_FcI);
      self->l_min_FcI += self->k_min_FcI;
      self->l_max_FcI += self->k_min_FcI;
      free(self->l_min_FcI);
      free(self->l_max_FcI);
    }
  }

  if (self->E_FcH != NULL) {
    for (cnt1 = self->k_min_FcH; cnt1 <= self->k_max_FcH; cnt1++)
      if (self->l_min_FcH[cnt1] < INF) {
        self->E_FcH[cnt1] += self->l_min_FcH[cnt1] / 2;
        free(self->E_FcH[cnt1]);
      }

    if (self->k_min_FcH < INF) {
      self->E_FcH += self->k_min_FcH;
      free(self->E_FcH);
      self->l_min_FcH += self->k_min_FcH;
      self->l_max_FcH += self->k_min_FcH;
      free(self->l_min_FcH);
      free(self->l_max_FcH);
    }
  }

  if (self->E_FcM != NULL) {
    for (cnt1 = self->k_min_FcM; cnt1 <= self->k_max_FcM; cnt1++)
      if (self->l_min_FcM[cnt1] < INF) {
        self->E_FcM[cnt1] += self->l_min_FcM[cnt1] / 2;
        free(self->E_FcM[cnt1]);
      }

    if (self->k_min_FcM < INF) {
      self->E_FcM += self->k_min_FcM;
      free(self->E_FcM);
      self->l_min_FcM += self->k_min_FcM;
      self->l_max_FcM += self->k_min_FcM;
      free(self->l_min_FcM);
      free(self->l_max_FcM);
    }
  }

  free(self->E_F5_rem);
  free(self->E_F3_rem);
  free(self->E_C_rem);
  free(self->E_M_rem);
  free(self->E_M1_rem);
  free(self->E_M2_rem);
}


PUBLIC void
vrna_mx_mfe_free(vrna_fold_compound_t *vc)
{
  if (vc) {
    vrna_mx_mfe_t *self = vc->matrices;
    if (self) {
      switch (self->type) {
        case VRNA_MX_DEFAULT:
          mfe_matrices_free_default(self);
          break;

        case VRNA_MX_WINDOW:
          mfe_matrices_free_window(self, vc->length, vc->window_size);
          break;

        case VRNA_MX_2DFOLD:
          mfe_matrices_free_2Dfold(self,
                                   vc->length,
                                   vc->params->model_details.min_loop_size,
                                   vc->iindx);
          break;

        default:                /* do nothing */
          break;
      }
      free(self);
      vc->matrices = NULL;
    }
  }
}



PRIVATE INLINE void
nullify_mfe(vrna_mx_mfe_t *mx)
{
  if (mx) {
    mx->length  = 0;
    mx->strands = 0;

    switch (mx->type) {
      case VRNA_MX_DEFAULT:
        mx->c     = NULL;
        mx->f5    = NULL;
        mx->f3    = NULL;
        mx->fms5  = NULL;
        mx->fms3  = NULL;
        mx->fML   = NULL;
        mx->fM1   = NULL;
        mx->fM2   = NULL;
        mx->ggg   = NULL;
        mx->Fc    = INF;
        mx->FcH   = INF;
        mx->FcI   = INF;
        mx->FcM   = INF;
        break;

      case VRNA_MX_WINDOW:
        mx->c_local   = NULL;
        mx->f3_local  = NULL;
        mx->fML_local = NULL;
        mx->ggg_local = NULL;
        break;

      case VRNA_MX_2DFOLD:
        mx->E_F5      = NULL;
        mx->l_min_F5  = NULL;
        mx->l_max_F5  = NULL;
        mx->k_min_F5  = NULL;
        mx->k_max_F5  = NULL;
        mx->E_F5_rem  = NULL;

        mx->E_F3      = NULL;
        mx->l_min_F3  = NULL;
        mx->l_max_F3  = NULL;
        mx->k_min_F3  = NULL;
        mx->k_max_F3  = NULL;
        mx->E_F3_rem  = NULL;

        mx->E_C     = NULL;
        mx->l_min_C = NULL;
        mx->l_max_C = NULL;
        mx->k_min_C = NULL;
        mx->k_max_C = NULL;
        mx->E_C_rem = NULL;

        mx->E_M     = NULL;
        mx->l_min_M = NULL;
        mx->l_max_M = NULL;
        mx->k_min_M = NULL;
        mx->k_max_M = NULL;
        mx->E_M_rem = NULL;

        mx->E_M1      = NULL;
        mx->l_min_M1  = NULL;
        mx->l_max_M1  = NULL;
        mx->k_min_M1  = NULL;
        mx->k_max_M1  = NULL;
        mx->E_M1_rem  = NULL;

        mx->E_M2      = NULL;
        mx->l_min_M2  = NULL;
        mx->l_max_M2  = NULL;
        mx->k_min_M2  = NULL;
        mx->k_max_M2  = NULL;
        mx->E_M2_rem  = NULL;

        mx->E_Fc      = NULL;
        mx->l_min_Fc  = NULL;
        mx->l_max_Fc  = NULL;
        mx->k_min_Fc  = 0;
        mx->k_max_Fc  = 0;
        mx->E_Fc_rem  = INF;

        mx->E_FcH     = NULL;
        mx->l_min_FcH = NULL;
        mx->l_max_FcH = NULL;
        mx->k_min_FcH = 0;
        mx->k_max_FcH = 0;
        mx->E_FcH_rem = INF;

        mx->E_FcI     = NULL;
        mx->l_min_FcI = NULL;
        mx->l_max_FcI = NULL;
        mx->k_min_FcI = 0;
        mx->k_max_FcI = 0;
        mx->E_FcI_rem = INF;

        mx->E_FcM     = NULL;
        mx->l_min_FcM = NULL;
        mx->l_max_FcM = NULL;
        mx->k_min_FcM = 0;
        mx->k_max_FcM = 0;
        mx->E_FcM_rem = INF;

#ifdef COUNT_STATES
        mx->N_F5  = NULL;
        mx->N_C   = NULL;
        mx->N_M   = NULL;
        mx->N_M1  = NULL;
#endif

        break;
    }
  }
}



PRIVATE vrna_mx_mfe_t *
init_mx_mfe_default(vrna_fold_compound_t  *fc,
                    unsigned int          alloc_vector)
{
  unsigned int  n, size, lin_size, s, strands;
  vrna_mx_mfe_t *mx;
  vrna_mx_mfe_t init = {
    .type = VRNA_MX_DEFAULT
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    printf("init_mx_mfe_default(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_mfe_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_mfe_t));
    nullify_mfe(mx);

    strands     = fc->strands;
    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;
    mx->strands = strands;

    if (alloc_vector & ALLOC_F5)
      mx->f5 = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_F3)
      mx->f3 = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_MULTISTRAND) {
      mx->fms5  = (int **)vrna_alloc(sizeof(int *) * strands);
      mx->fms3  = (int **)vrna_alloc(sizeof(int *) * strands);

      for (s = 0; s < strands; s++) {
        mx->fms5[s] = (int *)vrna_alloc(sizeof(int) * (n + 1));
        mx->fms3[s] = (int *)vrna_alloc(sizeof(int) * (n + 1));
      }
    }

    if (alloc_vector & ALLOC_C)
      mx->c = (int *)vrna_alloc(sizeof(int) * size);

    if (alloc_vector & ALLOC_FML)
      mx->fML = (int *)vrna_alloc(sizeof(int) * size);

    if (alloc_vector & ALLOC_UNIQ)
      mx->fM1 = (int *)vrna_alloc(sizeof(int) * size);

    if (alloc_vector & ALLOC_CIRC)
      mx->fM2 = (int *)vrna_alloc(sizeof(int) * lin_size);
  }

  return mx;
}


PRIVATE vrna_mx_mfe_t *
init_mx_mfe_window(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector)
{
  unsigned int  n, m, lin_size;
  vrna_mx_mfe_t *mx;
  vrna_mx_mfe_t init = {
    .type = VRNA_MX_WINDOW
  };

  n = fc->length;
  m = fc->window_size;

  if ((int)(n * m) >= (int)INT_MAX) {
    printf("init_mx_mfe_window(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_mfe_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_mfe_t));
    nullify_mfe(mx);

    lin_size    = n + 2;
    mx->length  = n;
    mx->strands = fc->strands;

    if (alloc_vector & ALLOC_F3)
      mx->f3_local = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_C)
      mx->c_local = (int **)vrna_alloc(sizeof(int *) * lin_size);

    if (alloc_vector & ALLOC_FML)
      mx->fML_local = (int **)vrna_alloc(sizeof(int *) * lin_size);
  }

  return mx;
}




PRIVATE vrna_mx_mfe_t *
init_mx_mfe_2Dfold(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector)
{
  unsigned int  n, i, size, lin_size;
  vrna_mx_mfe_t *mx;
  vrna_mx_mfe_t init = {
    .type = VRNA_MX_2DFOLD
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    printf("init_mx_mfe_2Dfold(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_mfe_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_mfe_t));
    nullify_mfe(mx);

    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;
    mx->strands = fc->strands;

    if (alloc_vector & ALLOC_F5) {
      mx->E_F5      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
      mx->l_min_F5  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_F5  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_F5  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_F5  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->E_F5_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
      for (i = 0; i <= n; i++)
        mx->E_F5_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_F3) {
      mx->E_F3      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
      mx->l_min_F3  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_F3  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_F3  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_F3  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->E_F3_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
      for (i = 0; i <= n; i++)
        mx->E_F3_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_C) {
      mx->E_C     = (int ***)vrna_alloc(sizeof(int **) * size);
      mx->l_min_C = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_C = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_C = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_C = (int *)vrna_alloc(sizeof(int) * size);
      mx->E_C_rem = (int *)vrna_alloc(sizeof(int) * size);
      for (i = 0; i < size; i++)
        mx->E_C_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_FML) {
      mx->E_M     = (int ***)vrna_alloc(sizeof(int **) * size);
      mx->l_min_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->E_M_rem = (int *)vrna_alloc(sizeof(int) * size);
      for (i = 0; i < size; i++)
        mx->E_M_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_UNIQ) {
      mx->E_M1      = (int ***)vrna_alloc(sizeof(int **) * size);
      mx->l_min_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->E_M1_rem  = (int *)vrna_alloc(sizeof(int) * size);
      for (i = 0; i < size; i++)
        mx->E_M1_rem[i] = INF;
    }

    if (alloc_vector & ALLOC_CIRC) {
      mx->E_M2      = (int ***)vrna_alloc(sizeof(int **) * lin_size);
      mx->l_min_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->E_M2_rem  = (int *)vrna_alloc(sizeof(int) * lin_size);
      for (i = 0; i <= n; i++)
        mx->E_M2_rem[i] = INF;
    }

#ifdef COUNT_STATES
    mx->N_C   = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
    mx->N_F5  = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * lin_size);
    mx->N_M   = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
    mx->N_M1  = (unsigned long ***)vrna_alloc(sizeof(unsigned long **) * size);
#endif
  }

  return mx;
}



PRIVATE int
add_mfe_matrices(vrna_fold_compound_t *vc,
                 vrna_mx_type_e       mx_type,
                 unsigned int         alloc_vector)
{
  if (vc) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        vc->matrices = init_mx_mfe_default(vc, alloc_vector);
        break;

      case VRNA_MX_WINDOW:
        vc->matrices = init_mx_mfe_window(vc, alloc_vector);
        break;

      case VRNA_MX_2DFOLD:
        vc->matrices = init_mx_mfe_2Dfold(vc, alloc_vector);
        break;

      default:
        return 0;
    }

    if (!vc->matrices)
      return 0;
  }

  return 1;
}




PUBLIC int
vrna_mx_mfe_add(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        mx_type,
                unsigned int          options)
{
  unsigned int mx_alloc_vector;

  if (vc->params) {
    options |= VRNA_OPTION_MFE;

    mx_alloc_vector = get_mx_alloc_vector(vc,
                                          mx_type,
                                          options);
    vrna_mx_mfe_free(vc);
    return add_mfe_matrices(vc, mx_type, mx_alloc_vector);
  }

  return 0;
}




PRIVATE unsigned int
get_mx_pf_alloc_vector_current(vrna_mx_pf_t   *mx,
                               vrna_mx_type_e mx_type)
{
  unsigned int mx_alloc_vector = ALLOC_NOTHING;

  if (mx) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        if (mx->q)
          mx_alloc_vector |= ALLOC_F;

        if (mx->qb)
          mx_alloc_vector |= ALLOC_C;

        if (mx->qm)
          mx_alloc_vector |= ALLOC_FML;

        if (mx->qm1)
          mx_alloc_vector |= ALLOC_UNIQ;

        if (mx->qm2)
          mx_alloc_vector |= ALLOC_CIRC;

        if (mx->probs)
          mx_alloc_vector |= ALLOC_PROBS;

        if (mx->q1k && mx->qln)
          mx_alloc_vector |= ALLOC_AUX;

        break;

      default:
        break;
    }
  }

  return mx_alloc_vector;
}



PRIVATE void
pf_matrices_free_default(vrna_mx_pf_t *self)
{
  free(self->q);
  free(self->qb);
  free(self->qm);
  free(self->qm1);
  free(self->qm2);
  free(self->probs);
  free(self->G);
  free(self->q1k);
  free(self->qln);
}


PRIVATE void
pf_matrices_free_window(vrna_mx_pf_t  *self,
                        unsigned int  length,
                        unsigned int  window_size)
{
  free(self->q_local);
  free(self->qb_local);
  free(self->qm_local);
  free(self->qm2_local);
  free(self->pR);
  free(self->QI5);
  free(self->q2l);
  free(self->qmb);
  free(self->G_local);
}


PRIVATE void
pf_matrices_free_2Dfold(vrna_mx_pf_t  *self,
                        unsigned int  length,
                        int           turn,
                        int           *indx,
                        int           *jindx)
{
  unsigned int  i, j, ij;
  int           cnt1;

  /* This will be some fun... */
  if (self->Q != NULL) {
    for (i = 1; i <= length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->Q[ij])
          continue;

        for (cnt1 = self->k_min_Q[ij]; cnt1 <= self->k_max_Q[ij]; cnt1++)
          if (self->l_min_Q[ij][cnt1] < INF) {
            self->Q[ij][cnt1] += self->l_min_Q[ij][cnt1] / 2;
            free(self->Q[ij][cnt1]);
          }

        if (self->k_min_Q[ij] < INF) {
          self->Q[ij] += self->k_min_Q[ij];
          free(self->Q[ij]);
          self->l_min_Q[ij] += self->k_min_Q[ij];
          self->l_max_Q[ij] += self->k_min_Q[ij];
          free(self->l_min_Q[ij]);
          free(self->l_max_Q[ij]);
        }
      }
    }
  }

  free(self->Q);
  free(self->l_min_Q);
  free(self->l_max_Q);
  free(self->k_min_Q);
  free(self->k_max_Q);

  if (self->Q_B != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->Q_B[ij])
          continue;

        for (cnt1 = self->k_min_Q_B[ij]; cnt1 <= self->k_max_Q_B[ij]; cnt1++)
          if (self->l_min_Q_B[ij][cnt1] < INF) {
            self->Q_B[ij][cnt1] += self->l_min_Q_B[ij][cnt1] / 2;
            free(self->Q_B[ij][cnt1]);
          }

        if (self->k_min_Q_B[ij] < INF) {
          self->Q_B[ij] += self->k_min_Q_B[ij];
          free(self->Q_B[ij]);
          self->l_min_Q_B[ij] += self->k_min_Q_B[ij];
          self->l_max_Q_B[ij] += self->k_min_Q_B[ij];
          free(self->l_min_Q_B[ij]);
          free(self->l_max_Q_B[ij]);
        }
      }
    }
  }

  free(self->Q_B);
  free(self->l_min_Q_B);
  free(self->l_max_Q_B);
  free(self->k_min_Q_B);
  free(self->k_max_Q_B);

  if (self->Q_M != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = indx[i] - j;
        if (!self->Q_M[ij])
          continue;

        for (cnt1 = self->k_min_Q_M[ij]; cnt1 <= self->k_max_Q_M[ij]; cnt1++)
          if (self->l_min_Q_M[ij][cnt1] < INF) {
            self->Q_M[ij][cnt1] += self->l_min_Q_M[ij][cnt1] / 2;
            free(self->Q_M[ij][cnt1]);
          }

        if (self->k_min_Q_M[ij] < INF) {
          self->Q_M[ij] += self->k_min_Q_M[ij];
          free(self->Q_M[ij]);
          self->l_min_Q_M[ij] += self->k_min_Q_M[ij];
          self->l_max_Q_M[ij] += self->k_min_Q_M[ij];
          free(self->l_min_Q_M[ij]);
          free(self->l_max_Q_M[ij]);
        }
      }
    }
  }

  free(self->Q_M);
  free(self->l_min_Q_M);
  free(self->l_max_Q_M);
  free(self->k_min_Q_M);
  free(self->k_max_Q_M);

  if (self->Q_M1 != NULL) {
    for (i = 1; i < length; i++) {
      for (j = i; j <= length; j++) {
        ij = jindx[j] + i;
        if (!self->Q_M1[ij])
          continue;

        for (cnt1 = self->k_min_Q_M1[ij]; cnt1 <= self->k_max_Q_M1[ij]; cnt1++)
          if (self->l_min_Q_M1[ij][cnt1] < INF) {
            self->Q_M1[ij][cnt1] += self->l_min_Q_M1[ij][cnt1] / 2;
            free(self->Q_M1[ij][cnt1]);
          }

        if (self->k_min_Q_M1[ij] < INF) {
          self->Q_M1[ij] += self->k_min_Q_M1[ij];
          free(self->Q_M1[ij]);
          self->l_min_Q_M1[ij]  += self->k_min_Q_M1[ij];
          self->l_max_Q_M1[ij]  += self->k_min_Q_M1[ij];
          free(self->l_min_Q_M1[ij]);
          free(self->l_max_Q_M1[ij]);
        }
      }
    }
  }

  free(self->Q_M1);
  free(self->l_min_Q_M1);
  free(self->l_max_Q_M1);
  free(self->k_min_Q_M1);
  free(self->k_max_Q_M1);

  if (self->Q_M2 != NULL) {
    for (i = 1; i < length - turn - 1; i++) {
      if (!self->Q_M2[i])
        continue;

      for (cnt1 = self->k_min_Q_M2[i]; cnt1 <= self->k_max_Q_M2[i]; cnt1++)
        if (self->l_min_Q_M2[i][cnt1] < INF) {
          self->Q_M2[i][cnt1] += self->l_min_Q_M2[i][cnt1] / 2;
          free(self->Q_M2[i][cnt1]);
        }

      if (self->k_min_Q_M2[i] < INF) {
        self->Q_M2[i] += self->k_min_Q_M2[i];
        free(self->Q_M2[i]);
        self->l_min_Q_M2[i] += self->k_min_Q_M2[i];
        self->l_max_Q_M2[i] += self->k_min_Q_M2[i];
        free(self->l_min_Q_M2[i]);
        free(self->l_max_Q_M2[i]);
      }
    }
  }

  free(self->Q_M2);
  free(self->l_min_Q_M2);
  free(self->l_max_Q_M2);
  free(self->k_min_Q_M2);
  free(self->k_max_Q_M2);

  if (self->Q_c != NULL) {
    for (cnt1 = self->k_min_Q_c; cnt1 <= self->k_max_Q_c; cnt1++)
      if (self->l_min_Q_c[cnt1] < INF) {
        self->Q_c[cnt1] += self->l_min_Q_c[cnt1] / 2;
        free(self->Q_c[cnt1]);
      }

    if (self->k_min_Q_c < INF) {
      self->Q_c += self->k_min_Q_c;
      free(self->Q_c);
      self->l_min_Q_c += self->k_min_Q_c;
      self->l_max_Q_c += self->k_min_Q_c;
      free(self->l_min_Q_c);
      free(self->l_max_Q_c);
    }
  }

  if (self->Q_cI != NULL) {
    for (cnt1 = self->k_min_Q_cI; cnt1 <= self->k_max_Q_cI; cnt1++)
      if (self->l_min_Q_cI[cnt1] < INF) {
        self->Q_cI[cnt1] += self->l_min_Q_cI[cnt1] / 2;
        free(self->Q_cI[cnt1]);
      }

    if (self->k_min_Q_cI < INF) {
      self->Q_cI += self->k_min_Q_cI;
      free(self->Q_cI);
      self->l_min_Q_cI  += self->k_min_Q_cI;
      self->l_max_Q_cI  += self->k_min_Q_cI;
      free(self->l_min_Q_cI);
      free(self->l_max_Q_cI);
    }
  }

  if (self->Q_cH != NULL) {
    for (cnt1 = self->k_min_Q_cH; cnt1 <= self->k_max_Q_cH; cnt1++)
      if (self->l_min_Q_cH[cnt1] < INF) {
        self->Q_cH[cnt1] += self->l_min_Q_cH[cnt1] / 2;
        free(self->Q_cH[cnt1]);
      }

    if (self->k_min_Q_cH < INF) {
      self->Q_cH += self->k_min_Q_cH;
      free(self->Q_cH);
      self->l_min_Q_cH  += self->k_min_Q_cH;
      self->l_max_Q_cH  += self->k_min_Q_cH;
      free(self->l_min_Q_cH);
      free(self->l_max_Q_cH);
    }
  }

  if (self->Q_cM != NULL) {
    for (cnt1 = self->k_min_Q_cM; cnt1 <= self->k_max_Q_cM; cnt1++)
      if (self->l_min_Q_cM[cnt1] < INF) {
        self->Q_cM[cnt1] += self->l_min_Q_cM[cnt1] / 2;
        free(self->Q_cM[cnt1]);
      }

    if (self->k_min_Q_cM < INF) {
      self->Q_cM += self->k_min_Q_cM;
      free(self->Q_cM);
      self->l_min_Q_cM  += self->k_min_Q_cM;
      self->l_max_Q_cM  += self->k_min_Q_cM;
      free(self->l_min_Q_cM);
      free(self->l_max_Q_cM);
    }
  }

  free(self->Q_rem);
  free(self->Q_B_rem);
  free(self->Q_M_rem);
  free(self->Q_M1_rem);
  free(self->Q_M2_rem);
}



PUBLIC void
vrna_mx_pf_free(vrna_fold_compound_t *vc)
{
  if (vc) {
    vrna_mx_pf_t *self = vc->exp_matrices;
    if (self) {
      switch (self->type) {
        case VRNA_MX_DEFAULT:
          pf_matrices_free_default(self);
          break;

        case VRNA_MX_WINDOW:
          pf_matrices_free_window(self, vc->length, vc->window_size);
          break;

        case VRNA_MX_2DFOLD:
          pf_matrices_free_2Dfold(self,
                                  vc->length,
                                  vc->exp_params->model_details.min_loop_size,
                                  vc->iindx,
                                  vc->jindx);
          break;

        default:                /* do nothing */
          break;
      }

      free(self->expMLbase);
      free(self->scale);

      free(self);
      vc->exp_matrices = NULL;
    }
  }
}



PRIVATE INLINE void
nullify_pf(vrna_mx_pf_t *mx)
{
  if (mx) {
    mx->length    = 0;
    mx->scale     = NULL;
    mx->expMLbase = NULL;

    switch (mx->type) {
      case VRNA_MX_DEFAULT:
        mx->q     = NULL;
        mx->qb    = NULL;
        mx->qm    = NULL;
        mx->qm1   = NULL;
        mx->qm2   = NULL;
        mx->probs = NULL;
        mx->q1k   = NULL;
        mx->qln   = NULL;
        break;

      case VRNA_MX_WINDOW:
        mx->q_local   = NULL;
        mx->qb_local  = NULL;
        mx->qm_local  = NULL;
        mx->qm2_local = NULL;
        mx->pR        = NULL;
        mx->QI5       = NULL;
        mx->q2l       = NULL;
        mx->qmb       = NULL;
        mx->G_local   = NULL;
        break;

      case VRNA_MX_2DFOLD:
        mx->Q       = NULL;
        mx->l_min_Q = NULL;
        mx->l_max_Q = NULL;
        mx->k_min_Q = NULL;
        mx->k_max_Q = NULL;
        mx->Q_rem   = NULL;

        mx->Q_B       = NULL;
        mx->l_min_Q_B = NULL;
        mx->l_max_Q_B = NULL;
        mx->k_min_Q_B = NULL;
        mx->k_max_Q_B = NULL;
        mx->Q_B_rem   = NULL;

        mx->Q_M       = NULL;
        mx->l_min_Q_M = NULL;
        mx->l_max_Q_M = NULL;
        mx->k_min_Q_M = NULL;
        mx->k_max_Q_M = NULL;
        mx->Q_M_rem   = NULL;

        mx->Q_M1        = NULL;
        mx->l_min_Q_M1  = NULL;
        mx->l_max_Q_M1  = NULL;
        mx->k_min_Q_M1  = NULL;
        mx->k_max_Q_M1  = NULL;
        mx->Q_M1_rem    = NULL;

        mx->Q_M2        = NULL;
        mx->l_min_Q_M2  = NULL;
        mx->l_max_Q_M2  = NULL;
        mx->k_min_Q_M2  = NULL;
        mx->k_max_Q_M2  = NULL;
        mx->Q_M2_rem    = NULL;

        mx->Q_c       = NULL;
        mx->l_min_Q_c = NULL;
        mx->l_max_Q_c = NULL;
        mx->k_min_Q_c = 0;
        mx->k_max_Q_c = 0;
        mx->Q_c_rem   = 0.;

        mx->Q_cH        = NULL;
        mx->l_min_Q_cH  = NULL;
        mx->l_max_Q_cH  = NULL;
        mx->k_min_Q_cH  = 0;
        mx->k_max_Q_cH  = 0;
        mx->Q_cH_rem    = 0.;

        mx->Q_cI        = NULL;
        mx->l_min_Q_cI  = NULL;
        mx->l_max_Q_cI  = NULL;
        mx->k_min_Q_cI  = 0;
        mx->k_max_Q_cI  = 0;
        mx->Q_cI_rem    = 0.;

        mx->Q_cM        = NULL;
        mx->l_min_Q_cM  = NULL;
        mx->l_max_Q_cM  = NULL;
        mx->k_min_Q_cM  = 0;
        mx->k_max_Q_cM  = 0;
        mx->Q_cM_rem    = 0.;

        break;
    }
  }
}




PRIVATE vrna_mx_pf_t *
init_mx_pf_default(vrna_fold_compound_t *fc,
                   unsigned int         alloc_vector)
{
  unsigned int  n, size, lin_size;
  vrna_mx_pf_t  *mx;
  vrna_mx_pf_t  init = {
    .type = VRNA_MX_DEFAULT
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    printf("init_mx_pf_default(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_pf_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_pf_t));
    nullify_pf(mx);

    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;

    if (alloc_vector & ALLOC_F)
      mx->q = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_C)
      mx->qb = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_FML)
      mx->qm = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_UNIQ)
      mx->qm1 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_CIRC)
      mx->qm2 = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);

    if (alloc_vector & ALLOC_PROBS)
      mx->probs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);

    if (alloc_vector & ALLOC_AUX) {
      mx->q1k = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
      mx->qln = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    }

    /*
     *  always alloc the helper arrays for unpaired nucleotides in multi-
     *  branch loops and scaling
     */
    mx->scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    mx->expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }

  return mx;
}


PRIVATE vrna_mx_pf_t *
init_mx_pf_window(vrna_fold_compound_t  *fc,
                  unsigned int          alloc_vector)
{
  unsigned int  n, m, lin_size;
  vrna_mx_pf_t  *mx;
  vrna_mx_pf_t  init = {
    .type = VRNA_MX_WINDOW
  };

  n = fc->length;
  m = fc->window_size;

  if ((int)(n * m) >= (int)INT_MAX) {
    printf("init_mx_pf_window(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_pf_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_pf_t));
    nullify_pf(mx);

    lin_size    = n + 2;
    mx->length  = n;

    if (alloc_vector & ALLOC_F)
      mx->q_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    if (alloc_vector & ALLOC_C)
      mx->qb_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    if (alloc_vector & ALLOC_FML)
      mx->qm_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    mx->pR = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);

    if (alloc_vector & ALLOC_PROBS) {
      mx->QI5       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
      mx->qmb       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
      mx->qm2_local = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
      mx->q2l       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * lin_size);
    }

    /*
     *  always alloc the helper arrays for unpaired nucleotides in multi-
     *  branch loops and scaling
     */
    mx->scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    mx->expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }

  return mx;
}


PRIVATE vrna_mx_pf_t *
init_mx_pf_2Dfold(vrna_fold_compound_t  *fc,
                  unsigned int          alloc_vector)
{
  unsigned int  n, size, lin_size;
  vrna_mx_pf_t  *mx;
  vrna_mx_pf_t  init = {
    .type = VRNA_MX_2DFOLD
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    printf("init_mx_pf_2Dfold(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_pf_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_pf_t));
    nullify_pf(mx);

    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;

    if (alloc_vector & ALLOC_F) {
      mx->Q       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_C) {
      mx->Q_B       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q_B = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q_B = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q_B = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q_B = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_B_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_FML) {
      mx->Q_M       = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q_M = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q_M = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_M_rem   = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_UNIQ) {
      mx->Q_M1        = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * size);
      mx->l_min_Q_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->l_max_Q_M1  = (int **)vrna_alloc(sizeof(int *) * size);
      mx->k_min_Q_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->k_max_Q_M1  = (int *)vrna_alloc(sizeof(int) * size);
      mx->Q_M1_rem    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    }

    if (alloc_vector & ALLOC_CIRC) {
      mx->Q_M2        = (FLT_OR_DBL ***)vrna_alloc(sizeof(FLT_OR_DBL * *) * lin_size);
      mx->l_min_Q_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->l_max_Q_M2  = (int **)vrna_alloc(sizeof(int *) * lin_size);
      mx->k_min_Q_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->k_max_Q_M2  = (int *)vrna_alloc(sizeof(int) * lin_size);
      mx->Q_M2_rem    = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    }

    /*
     *  always alloc the helper arrays for unpaired nucleotides in multi-
     *  branch loops and scaling
     */
    mx->scale     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
    mx->expMLbase = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * lin_size);
  }

  return mx;
}




PRIVATE int
add_pf_matrices(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        mx_type,
                unsigned int          alloc_vector)
{
  if (vc) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        vc->exp_matrices = init_mx_pf_default(vc, alloc_vector);
        break;

      case VRNA_MX_WINDOW:
        vc->exp_matrices = init_mx_pf_window(vc, alloc_vector);
        break;

      case VRNA_MX_2DFOLD:
        vc->exp_matrices = init_mx_pf_2Dfold(vc, alloc_vector);
        break;

      default:                /* do nothing */
        return 0;
    }

    if (!vc->exp_matrices)
      return 0;

    if (vc->exp_params->model_details.gquad) {
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          vc->exp_matrices->G = NULL;
          /* can't do that here, since scale[] is not filled yet :(
           * vc->exp_matrices->G = get_gquad_pf_matrix(vc->sequence_encoding2, vc->exp_matrices->scale, vc->exp_params);
           */
          break;
        default:                    /* do nothing */
          break;
      }
    }

    vrna_exp_params_rescale(vc, NULL);
  }

  return 1;
}



PUBLIC int
vrna_mx_pf_add(vrna_fold_compound_t *vc,
               vrna_mx_type_e       mx_type,
               unsigned int         options)
{
  unsigned int mx_alloc_vector;

  if (vc->exp_params) {
    mx_alloc_vector = get_mx_alloc_vector(vc,
                                          mx_type,
                                          options | VRNA_OPTION_PF);
    vrna_mx_pf_free(vc);
    return add_pf_matrices(vc, mx_type, mx_alloc_vector);
  }

  return 0;
}

PUBLIC int
vrna_mx_prepare(vrna_fold_compound_t  *vc,
                unsigned int          options)
{
  int             ret, realloc;
  unsigned int    mx_alloc_vector, mx_alloc_vector_current;
  vrna_mx_type_e  mx_type;

  ret = 1;

  if (vc) {
    /*  check whether we have the correct DP matrices attached, and if there is
     *  enough memory allocated
     */
    if (options & VRNA_OPTION_MFE) {
      /* prepare for MFE computation */
      if (options & VRNA_OPTION_WINDOW) /* Windowing approach, a.k.a. locally optimal */
        mx_type = VRNA_MX_WINDOW;
      else                              /* default is regular MFE */
        mx_type = VRNA_MX_DEFAULT;

      if (vc->strands > 1)
        options |= VRNA_OPTION_HYBRID;

      realloc = 0;

      if (!vc->matrices || (vc->matrices->type != mx_type) || (vc->matrices->length < vc->length)) {
        realloc = 1;
      } else {
        mx_alloc_vector =
          get_mx_alloc_vector(vc, mx_type, options);
        mx_alloc_vector_current = get_mx_mfe_alloc_vector_current(vc->matrices, mx_type);
        if ((mx_alloc_vector & mx_alloc_vector_current) != mx_alloc_vector)
          realloc = 1;
      }

      if (realloc) /* Add DP matrices, if not they are not present */
        ret &= vrna_mx_mfe_add(vc, mx_type, options);
    }

    if (options & VRNA_OPTION_PF) {
      /* prepare for partition function computations */
      if (!vc->exp_params) /* return failure if exp_params data is not present */
        return 0;

      if (options & VRNA_OPTION_WINDOW) /* Windowing approach, a.k.a. locally optimal */
        mx_type = VRNA_MX_WINDOW;
      else                              /* default is regular MFE */
        mx_type = VRNA_MX_DEFAULT;

      if (vc->strands > 1)
        options |= VRNA_OPTION_HYBRID;

      realloc = 0;

      /*  Add DP matrices, if not they are not present */
      if (!vc->exp_matrices || (vc->exp_matrices->type != mx_type) ||
          (vc->exp_matrices->length < vc->length)) {
        realloc = 1;
      } else {
        mx_alloc_vector = get_mx_alloc_vector(vc,
                                              mx_type,
                                              options);
        mx_alloc_vector_current = get_mx_pf_alloc_vector_current(vc->exp_matrices, mx_type);
        if ((mx_alloc_vector & mx_alloc_vector_current) != mx_alloc_vector)
          realloc = 1;
      }

      if (realloc) /* Add DP matrices, if not they are not present */
        ret &= vrna_mx_pf_add(vc, mx_type, options);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
      else   /* re-compute pf_scale and MLbase contributions (for RNAup)*/
        vrna_exp_params_rescale(vc, NULL);

#endif
    }
  } else {
    ret = 0;
  }

  return ret;
}



PUBLIC int
vrna_fold_compound_prepare(vrna_fold_compound_t *fc,
                           unsigned int         options)
{
  int ret = 1; /* success */

  /* check maximum sequence length restrictions */
  // if (fc->length > 1000000) {
  //   printf(
  //     "vrna_fold_compound_prepare@data_structures.c: sequence length of %d exceeds addressable range",
  //     fc->length);
  //   return 0;
  // }

  /* make sure to always provide sane bp-span settings */
  // sanitize_bp_span(fc, options);

  /* prepare Boltzmann factors if required */
  vrna_params_prepare(fc, options);

  /* prepare ptype array(s) */
  // vrna_ptypes_prepare(fc, options);

  if (options & VRNA_OPTION_PF) {
    /* prepare for partition function computation */

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:     /* get pre-computed Boltzmann factors if not present*/
        if (fc->domains_up)         /* turn on unique ML decomposition with qm1 array */
          fc->exp_params->model_details.uniq_ML = 1;

        break;

      default:
        /* not doing anything here... */
        break;
    }
  }

  /* prepare hard constraints */
  // vrna_hc_prepare(fc, options);

  /* prepare soft constraints data structure, if required */
//   vrna_sc_prepare(fc, options);

  /* Add DP matrices, if not they are not present or do not fit current settings */
  vrna_mx_prepare(fc, options);

  return ret;
}



PUBLIC FLT_OR_DBL
vrna_pf(vrna_fold_compound_t  *fc,
        char                  *structure)
{
  int               n;
  FLT_OR_DBL        Q, dG;
  vrna_md_t         *md;
  vrna_exp_param_t  *params;
  vrna_mx_pf_t      *matrices;

  dG = (FLT_OR_DBL)(INF / 100.);

  if (fc) {
    /* make sure, everything is set up properly to start partition function computations */
    if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_PF)) {
      printf("vrna_pf@part_func.c: Failed to prepare vrna_fold_compound");
      return dG;
    }

    n         = fc->length;
    params    = fc->exp_params;
    matrices  = fc->exp_matrices;
    md        = &(params->model_details);
    // md->helical_rise = VRNA_MODEL_HELICAL_RISE;
    // md->backbone_length = VRNA_MODEL_BACKBONE_LENGTH;
    // md->saltDPXInitFact = VRNA_MODEL_SALT_DPXINIT_FACT;
    /* call user-defined recursion status callback function */
    if (fc->stat_cb) // not into
      fc->stat_cb(VRNA_STATUS_PF_PRE, fc->auxdata);

    /* for now, multi-strand folding is implemented as additional grammar rule */
    if (fc->strands > 1)
      vrna_pf_multifold_prepare(fc);

    if (!fill_arrays(fc)) {
      return dG;
    }


    if (fc->strands > 1)
      vrna_gr_reset(fc);

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(VRNA_STATUS_PF_POST, fc->auxdata);

    switch (md->backtrack_type) {
      case 'C':
        Q = matrices->qb[fc->iindx[1] - n];
        break;

      case 'M':
        Q = matrices->qm[fc->iindx[1] - n];
        break;

      default:
        Q = (md->circ) ? matrices->qo : matrices->q[fc->iindx[1] - n];
        break;
    }

    /* ensemble free energy in Kcal/mol              */
    if (Q <= FLT_MIN)
      printf("pf_scale too large");

    if (fc->strands > 1) {
      /* check for rotational symmetry correction */
      unsigned int sym = vrna_rotational_symmetry(fc->sequence);
      Q /= (FLT_OR_DBL)sym;

      /* add interaction penalty */
      Q *= pow(params->expDuplexInit, (FLT_OR_DBL)(fc->strands - 1));
    }

    dG = (FLT_OR_DBL)((-log(Q) - n * log(params->pf_scale)) *
                      params->kT /
                      1000.0);

    /* calculate base pairing probability matrix (bppm)  compute in here !!!!!!!*/
    if (md->compute_bpp) {
      vrna_pairing_probs(fc, structure);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

      /*
       *  Backward compatibility:
       *  This block may be removed if deprecated functions
       *  relying on the global variable "pr" vanish from within the package!
       */
      pr = matrices->probs;

#endif
    }

  }
  printf("%f",dG);
  return dG;
}





PRIVATE void
extract_dimer_props(vrna_fold_compound_t  *fc,
                    double                *F0AB,  /**< @brief Null model without DuplexInit */
                    double                *FAB,   /**< @brief all states with DuplexInit correction */
                    double                *FcAB,  /**< @brief true hybrid states only */
                    double                *FA,    /**< @brief monomer A */
                    double                *FB     /**< @brief monomer B */
                    )
{
  unsigned int      n, sym, *ss, *so, *se;
  double            kT, QAB, QToT, Qzero;
  vrna_mx_pf_t      *matrices;
  vrna_exp_param_t  *params;

  n         = fc->length;
  ss        = fc->strand_start;
  se        = fc->strand_end;
  so        = fc->strand_order;
  params    = fc->exp_params;
  matrices  = fc->exp_matrices;

  if (fc->strands > 1) {
    kT  = params->kT / 1000.0;
    QAB = matrices->q[fc->iindx[1] - n];

    /* check for rotational symmetry correction */
    sym = vrna_rotational_symmetry(fc->sequence);
    QAB /= (FLT_OR_DBL)sym;

    /* add interaction penalty */
    QAB *= pow(params->expDuplexInit, (FLT_OR_DBL)(fc->strands - 1));

    Qzero = matrices->q[fc->iindx[1] - n] +
            matrices->q[fc->iindx[1] - se[so[0]]] *
            matrices->q[fc->iindx[ss[so[1]]] - n];

    QToT = matrices->q[fc->iindx[1] - se[so[0]]] *
           matrices->q[fc->iindx[ss[so[1]]] - n] +
           QAB;

    *FAB  = -kT * (log(QToT) + n * log(params->pf_scale));
    *F0AB = -kT * (log(Qzero) + n * log(params->pf_scale));
    *FcAB = (QAB > 1e-17) ? -kT * (log(QAB) + n * log(params->pf_scale)) : 999;
    *FA   = -kT *
            (log(matrices->q[fc->iindx[1] - se[so[0]]]) + (se[so[0]]) *
             log(params->pf_scale));
    *FB = -kT *
          (log(matrices->q[fc->iindx[ss[so[1]]] - n]) + (n - ss[so[1]] + 1) *
           log(params->pf_scale));
  } else {
    *FA = *FB = *FAB = *F0AB = (-log(matrices->q[fc->iindx[1] - n]) - n * log(params->pf_scale)) *
                               params->kT / 1000.0;
    *FcAB = 0;
  }
}



PUBLIC vrna_dimer_pf_t
vrna_pf_dimer(vrna_fold_compound_t  *fc,
              char                  *structure)
{
  vrna_dimer_pf_t X;

  X.F0AB = X.FAB = X.FcAB = X.FA = X.FB = 0.;
  /* structure changed by vrna_pf */
  if (fc) {
    (void)vrna_pf(fc, structure);

    /* backward compatibility partition function and ensemble energy computation */
    extract_dimer_props(fc,
                        &(X.F0AB),
                        &(X.FAB),
                        &(X.FcAB),
                        &(X.FA),
                        &(X.FB));
  }

  return X;
}


PUBLIC char *
vrna_cut_point_insert(const char  *string,
                      int         cp)
{
  char  *ctmp;
  int   len;

  if (cp > 0) {
    len   = strlen(string);
    ctmp  = (char *)vrna_alloc((len + 2) * sizeof(char));
    /* first sequence */
    (void)strncpy(ctmp, string, cp - 1);
    /* spacer */
    ctmp[cp - 1] = '&';
    /* second sequence */
    (void)strcat(ctmp, string + cp - 1);
  } else {
    ctmp = strdup(string);
  }

  return ctmp;
}

char *process_record(char *sequence){
  char                  *mfe_structure, **rec_rest;
  unsigned int          n, i;
  double                min_en, kT, *concentrations;
  vrna_ep_t             *prAB, *prAA, *prBB, *prA, *prB, *mfAB, *mfAA, *mfBB, *mfA, *mfB;
  struct options        opt;
  struct output_stream  *o_stream;
  size_t                **mod_positions;
  size_t                mod_param_sets;
  mfAB            = mfAA = mfBB = mfA = mfB = NULL;
  prAB            = prAA = prBB = prA = prB = NULL;
  concentrations  = NULL;
  init_default_options(&opt);
  o_stream        = (struct output_stream *)vrna_alloc(sizeof(struct output_stream));
  vrna_fold_compound_t *vc = vrna_fold_compound(sequence,
                                                &(opt.md),
                                                VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
  n = vc->length;
  // min_en  = vrna_mfe_dimer(vc, mfe_structure);
  min_en = -49.7;
  // pf calc

  char              *Astring, *Bstring, *orig_Astring, *orig_Bstring, *pairing_propensity;
  int               Blength, Alength;
  vrna_dimer_pf_t   AB, AA, BB;
  vrna_dimer_conc_t *conc_result;

  conc_result = NULL;
  prAB        = NULL;
  prAA        = NULL;
  prBB        = NULL;
  prA         = NULL;
  prB         = NULL;

  Astring = Bstring = orig_Astring = orig_Bstring = NULL;
  Alength = Blength = 0;

  pairing_propensity = (char *)vrna_alloc(sizeof(char) * (n + 1));
  vrna_exp_params_rescale(vc, &min_en); // expmismatchExt(same global),1nI(same global),<int11,int21,int22>,tetra,tri,hex,pf_scale,bug in here!!!!!!!!
  kT = vc->exp_params->kT / 1000.;
  // pairing_propensity is the finnal structure
  AB = AA = BB = vrna_pf_dimer(vc, pairing_propensity); /* exp_matrices changed*/
  char *costruc;
  // prAB = vrna_plist_from_probs(vc, opt->bppmThreshold);
  costruc = vrna_cut_point_insert(pairing_propensity, vc->cutpoint);
  // printf("%s", costruc);
  return costruc;
}




int main_cofold(){
  char *s1 = "ACTGCCAAGTAGGAAAGTCCCATAAGGTCAT&TGAACTTATGGGACTTTCCTACTTGGCAG";
  char *res;
  res = process_record(s1);
  printf("%s",res);
  return 0;
}

int main(){
  main_dup();
  main_cofold();
}




