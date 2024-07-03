
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <myall.h>


PUBLIC int    subopt_sorted = 0;      /* output sorted by energy */
static short  alias[MAXALPHA + 1];
static int    pair[MAXALPHA + 1][MAXALPHA + 1];
static int    rtype[8] = {
  0, 2, 1, 4, 3, 6, 5, 7
};

/* Standard streams.  */
extern FILE *stdin;		/* Standard input stream.  */
extern FILE *stdout;		/* Standard output stream.  */
extern FILE *stderr;		/* Standard error output stream.  */
/* C89/C99 say they're macros.  Make them happy.  */
#define stdin stdin
#define stdout stdout
#define stderr stderr

PRIVATE vrna_param_t  *P = NULL;
PRIVATE int           **c = NULL;     /* energy array, given that i-j pair */
PRIVATE short         *S1 = NULL, *SS1 = NULL, *S2 = NULL, *SS2 = NULL;
PRIVATE int           n1, n2;         /* sequence lengths */

struct vrna_fc_s
{
  /* data */
  int i,j;
};


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
      vrna_message_error("backtrack failed in fold duplex");

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
        vrna_message_error("backtrack failed in fold duplex");
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

PUBLIC void *
vrna_alloc(unsigned size)
{
  void *pointer;

  if ((pointer = (void *)calloc(1, (size_t)size)) == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "vrna_alloc: requested size: %d\n", size);
      vrna_message_error("Memory allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    vrna_message_error("Memory allocation failure -> no memory");
  }

  return pointer;
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

PRIVATE int
compare(const void  *sub1,
        const void  *sub2)
{
  int d;

  if (((duplexT *)sub1)->energy > ((duplexT *)sub2)->energy)
    return 1;

  if (((duplexT *)sub1)->energy < ((duplexT *)sub2)->energy)
    return -1;

  d = ((duplexT *)sub1)->i - ((duplexT *)sub2)->i;
  if (d != 0)
    return d;

  return ((duplexT *)sub1)->j - ((duplexT *)sub2)->j;
}

PUBLIC int
vrna_eval_ext_stem(vrna_fold_compound_t  *fc,
                   int                   i,
                   int                   j)
{
  char                      *ptype;
  short                     *S;
  unsigned int              type;
  int                       ij, en, e, *idx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc;
  vrna_hc_eval_f evaluate;
  struct hc_ext_def_dat     hc_dat_local;

  S         = fc->sequence_encoding;
  idx       = fc->jindx;
  ptype     = fc->ptype;
  P         = fc->params;
  md        = &(P->model_details);
  sc        = fc->sc;
  evaluate  = prepare_hc_ext_def(fc, &hc_dat_local);

  e     = INF;
  ij    = idx[j] + i;
  type  = vrna_get_ptype(ij, ptype);

  if (evaluate(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
    switch (md->dangles) {
      case 2:
        e = vrna_E_ext_stem(type, S[i - 1], S[j + 1], P);
        break;

      case 0:
      /* fall through */

      default:
        e = vrna_E_ext_stem(type, -1, -1, P);
        break;
    }
    if (sc)
      if (sc->f)
        e += sc->f(i, j, i, j, VRNA_DECOMP_EXT_STEM, sc->data);
  }

  if (md->dangles % 2) {
    ij = idx[j - 1] + i;
    if (evaluate(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
      type = vrna_get_ptype(ij, ptype);

      en = vrna_E_ext_stem(type, -1, S[j], P);

      if (sc)
        if (sc->f)
          en += sc->f(i, j, i, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);

      e = MIN2(e, en);
    }

    ij = idx[j] + i + 1;
    if (evaluate(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
      type = vrna_get_ptype(ij, ptype);

      en = vrna_E_ext_stem(type, S[i], -1, P);

      if (sc)
        if (sc->f)
          en += sc->f(i, j, i + 1, j, VRNA_DECOMP_EXT_STEM, sc->data);

      e = MIN2(e, en);
    }

    ij = idx[j - 1] + i + 1;
    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
      type = vrna_get_ptype(ij, ptype);

      en = vrna_E_ext_stem(type, S[i], S[j], P);

      if (sc)
        if (sc->f)
          en += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_EXT_STEM, sc->data);

      e = MIN2(e, en);
    }
  }

  return e;
}
// model
double          temperature     = VRNA_MODEL_DEFAULT_TEMPERATURE;
double          pf_scale        = VRNA_MODEL_DEFAULT_PF_SCALE;
int             dangles         = VRNA_MODEL_DEFAULT_DANGLES;
int             tetra_loop      = VRNA_MODEL_DEFAULT_SPECIAL_HP;
int             noLonelyPairs   = VRNA_MODEL_DEFAULT_NO_LP;
int             noGU            = VRNA_MODEL_DEFAULT_NO_GU;
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


PUBLIC void
set_model_details(vrna_md_t *md)
{
  if (md) {
    /* make sure there are no uninitialized data fields */
    memset(md, 0, sizeof(vrna_md_t));

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

    if (nonstandards)
      copy_nonstandards(md, nonstandards);

    /* set default values for the pair/rtype[pair] stuff */
    vrna_md_update(md);
  }
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
static INLINE void
make_pair_matrix(void)
{
  int i, j;

  if (energy_set == 0) {
    for (i = 0; i < 5; i++)
      alias[i] = (short)i;
    alias[5]  = 3;  /* X <-> G */
    alias[6]  = 2;  /* K <-> C */
    alias[7]  = 0;  /* I <-> default base '@' */
    for (i = 0; i < NBASES; i++)
      for (j = 0; j < NBASES; j++)
        pair[i][j] = BP_pair[i][j];
    if (noGU)
      pair[3][4] = pair[4][3] = 0;

    if (nonstandards != NULL) {
      /* allow nonstandard bp's */
      for (i = 0; i < (int)strlen(nonstandards); i += 2)
        pair[encode_char(nonstandards[i])]
        [encode_char(nonstandards[i + 1])] = 7;
    }

    for (i = 0; i < NBASES; i++)
      for (j = 0; j < NBASES; j++)
        rtype[pair[i][j]] = pair[j][i];
  } else {
    for (i = 0; i <= MAXALPHA; i++)
      for (j = 0; j <= MAXALPHA; j++)
        pair[i][j] = 0;
    if (energy_set == 1) {
      for (i = 1; i < MAXALPHA; ) {
        alias[i++]  = 3;      /* A <-> G */
        alias[i++]  = 2;      /* B <-> C */
      }
      for (i = 1; i < MAXALPHA; i++) {
        pair[i][i + 1] = 2;       /* AB <-> GC */
        i++;
        pair[i][i - 1] = 1;       /* BA <-> CG */
      }
    } else if (energy_set == 2) {
      for (i = 1; i < MAXALPHA; ) {
        alias[i++]  = 1;      /* A <-> A*/
        alias[i++]  = 4;      /* B <-> U */
      }
      for (i = 1; i < MAXALPHA; i++) {
        pair[i][i + 1] = 5;       /* AB <-> AU */
        i++;
        pair[i][i - 1] = 6;       /* BA <-> UA */
      }
    } else if (energy_set == 3) {
      for (i = 1; i < MAXALPHA - 2; ) {
        alias[i++]  = 3;    /* A <-> G */
        alias[i++]  = 2;    /* B <-> C */
        alias[i++]  = 1;    /* C <-> A */
        alias[i++]  = 4;    /* D <-> U */
      }
      for (i = 1; i < MAXALPHA - 2; i++) {
        pair[i][i + 1] = 2;     /* AB <-> GC */
        i++;
        pair[i][i - 1] = 1;     /* BA <-> CG */
        i++;
        pair[i][i + 1] = 5;     /* CD <-> AU */
        i++;
        pair[i][i - 1] = 6;     /* DC <-> UA */
      }
    } else {
      vrna_message_error("What energy_set are YOU using??");
    }

    for (i = 0; i <= MAXALPHA; i++)
      for (j = 0; j <= MAXALPHA; j++)
        rtype[pair[i][j]] = pair[j][i];
  }
}
static INLINE short *
encode_sequence(const char  *sequence,
                short       how)
{
  unsigned int  i, l = (unsigned int)strlen(sequence);
  short         *S = (short *)vrna_alloc(sizeof(short) * (l + 2));

  switch (how) {
    /* standard encoding as always used for S */
    case 0:
      for (i = 1; i <= l; i++)    /* make numerical encoding of sequence */
        S[i] = (short)encode_char(sequence[i - 1]);
      S[l + 1]  = S[1];
      S[0]      = (short)l;
      break;
    /* encoding for mismatches of nostandard bases (normally used for S1) */
    case 1:
      for (i = 1; i <= l; i++)
        S[i] = alias[(short)encode_char(sequence[i - 1])];
      S[l + 1]  = S[1];
      S[0]      = S[l];
      break;
  }

  return S;
}
// generate result if no -e option 
PRIVATE duplexT
duplexfold_cu(const char  *s1,
              const char  *s2,
              int         clean_up)
{
  int       i, j, Emin = INF, i_min = 0, j_min = 0;
  char      *struc;
  duplexT   mfe;
  vrna_md_t md;

  n1  = (int)strlen(s1);
  n2  = (int)strlen(s2);

  set_model_details(&md);
  if ((!P) || (fabs(P->temperature - temperature) > 1e-6)) {
    if (P)
      free(P);

    P = vrna_params(&md);
    make_pair_matrix();
  }

  c = (int **)vrna_alloc(sizeof(int *) * (n1 + 1));
  for (i = 1; i <= n1; i++)
    c[i] = (int *)vrna_alloc(sizeof(int) * (n2 + 1));

  S1  = encode_sequence(s1, 0);
  S2  = encode_sequence(s2, 0);
  SS1 = encode_sequence(s1, 1);
  SS2 = encode_sequence(s2, 1);

  for (i = 1; i <= n1; i++) {
    for (j = n2; j > 0; j--) {
      int type, type2, E, k, l;
      type    = pair[S1[i]][S2[j]];
      c[i][j] = type ? P->DuplexInit : INF;
      if (!type)
        continue;

      c[i][j] += vrna_E_ext_stem(type, (i > 1) ? SS1[i - 1] : -1, (j < n2) ? SS2[j + 1] : -1, P);
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
        }
      }
      E = c[i][j];
      E += vrna_E_ext_stem(rtype[type], (j > 1) ? SS2[j - 1] : -1, (i < n1) ? SS1[i + 1] : -1, P);
      if (E < Emin) {
        Emin  = E;
        i_min = i;
        j_min = j;
      }
    }
  }

  struc = backtrack(i_min, j_min);
  if (i_min < n1)
    i_min++;

  if (j_min > 1)
    j_min--;

  mfe.i         = i_min;
  mfe.j         = j_min;
  mfe.energy    = (float)Emin / 100.;
  mfe.structure = struc;
  if (clean_up) {
    for (i = 1; i <= n1; i++)
      free(c[i]);
    free(c);
    free(S1);
    free(S2);
    free(SS1);
    free(SS2);
  }

  return mfe;
}

PUBLIC int
vrna_E_ext_stem(unsigned int  type,
                int           n5d,
                int           n3d,
                vrna_param_t  *p)
{
  int energy = 0;

  if (n5d >= 0 && n3d >= 0)
    energy += p->mismatchExt[type][n5d][n3d];
  else if (n5d >= 0)
    energy += p->dangle5[type][n5d];
  else if (n3d >= 0)
    energy += p->dangle3[type][n3d];

  if (type > 2)
    energy += p->TerminalAU;

  return energy;
}

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
  thresh  = (int)mfe.energy * 100 + 0.1 + delta;
  n1      = strlen(s1);
  n2      = strlen(s2);

  // 二层遍历查找次优结构
  for (i = n1; i > 0; i--) {
    for (j = 1; j <= n2; j++) {
      int type, ii, jj, Ed;
      type = pair[S2[j]][S1[i]];
      if (!type)
        continue;
      // 计算能量
      E   = Ed = c[i][j];
      // from /data/ntc/Repository/ViennaRNA-2.6.4/src/ViennaRNA/loops/external.c
      Ed  += vrna_E_ext_stem(type, (j > 1) ? SS2[j - 1] : -1, (i < n1) ? SS1[i + 1] : -1, P);
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
      if (!type)
        continue;

      struc = backtrack(i, j);
      vrna_message_info(stderr, "%d %d %d", i, j, E);
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

  if (subopt_sorted)
    qsort(subopt, n_subopt, sizeof(duplexT), compare);

  subopt[n_subopt].i          = 0;
  subopt[n_subopt].j          = 0;
  subopt[n_subopt].structure  = NULL;
  return subopt;
}