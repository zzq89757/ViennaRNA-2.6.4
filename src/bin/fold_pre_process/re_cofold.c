#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "re_cofold.h"
#include "matrices/all_mat.h"
#include "dup.h"


void
init_default_options(struct options *opt)
{
  opt->filename_full  = 0;
  opt->filename_delim = NULL;
  opt->pf             = 1;
  opt->doT            = 0; /* compute dimer free energies etc. */
  opt->noPS           = 1;
  opt->noconv         = 0;
  opt->centroid       = 0;  /* off by default due to historical reasons */
  opt->MEA            = 0;
  opt->MEAgamma       = 1.;
  opt->bppmThreshold  = 1e-5;
  opt->verbose        = 0;
  // opt->commands       = NULL;
  // opt->id_control     = NULL;
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

  // opt->mod_params         = NULL;

  opt->csv_output       = 0;    /* flag indicating whether we produce one-line outputs, a.k.a. CSV */
  opt->csv_header       = 1;    /* print header for one-line output */
  opt->csv_output_delim = ',';  /* delimiting character for one-line output */

  opt->jobs               = 1;
  opt->keep_order         = 1;
  opt->next_record_number = 0;
  // opt->output_queue       = NULL;
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

    // fc->domains_struc = NULL;
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

      case VRNA_FC_TYPE_COMPARATIVE:
        /*
         *  for now, comparative structure prediction does not allow for RNA-RNA interactions,
         *  so we pretend working on a single strand
         */
        fc->nucleotides = (vrna_seq_t *)vrna_realloc(fc->nucleotides,
                                                     sizeof(vrna_seq_t) * (fc->strands + 1));
        fc->nucleotides[0].string = NULL;
        fc->nucleotides[0].type   = VRNA_SEQ_RNA;
        fc->nucleotides[0].length = fc->length;

        /* 1. store initial strand order */
        fc->strand_order_uniq = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_order      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);

        /* 2. mark start and end positions of sequences */
        fc->strand_start    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_end      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
        fc->strand_start[0] = 1;
        fc->strand_end[0]   = fc->strand_start[0] + fc->length - 1;

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

      if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
        kT /= vc->n_seq;

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
  pscore      = (vc->type == VRNA_FC_TYPE_COMPARATIVE) ? vc->pscore : NULL;
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

    compute_bpp_int = &compute_bpp_internal;
    compute_bpp_mul = &compute_bpp_multibranch;

    Qmax = 0;
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

    /* init diagonal entries unable to pair in pr matrix */
    for (i = 1; i <= n; i++)
      probs[my_iindx[i] - i] = 0.;
    // printf("probs is :%f", probs);
    /* 1. external loop pairs, i.e. pairs not enclosed by any other pair (or external loop for circular RNAs) */
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
    // linked vc 
    for (i = 1; i <= n; i++)
      for (j = i + 1; j <= n; j++) {
        ij = my_iindx[i] - j;

        if (qb[ij] > 0.) {
          probs[ij] *= qb[ij];
          }
        // printf("%d-", probs[ij]);
      }
    char *s = vrna_db_from_probs(probs, (unsigned int)n);
    /* strcuture here  generate by s  ---------------- */
    memcpy(structure, s, n);
    structure[n] = '\0';
    free(s);

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
    // if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_PF)) {
    //   printf("vrna_pf@part_func.c: Failed to prepare vrna_fold_compound");
    //   return dG;
    // }

    n         = fc->length;
    params    = fc->exp_params;
    matrices  = fc->exp_matrices;
    md        = &(params->model_details);

// #ifdef _OPENMP
//     /* Explicitly turn off dynamic threads */
//     omp_set_dynamic(0);
// #endif

// #ifdef SUN4
//     nonstandard_arithmetic();
// #elif defined(HP9)
//     fpsetfastmode(1);
// #endif

    /* call user-defined recursion status callback function */
    if (fc->stat_cb) // not into
      fc->stat_cb(VRNA_STATUS_PF_PRE, fc->auxdata);

    /* for now, multi-strand folding is implemented as additional grammar rule */
    if (fc->strands > 1)
      vrna_pf_multifold_prepare(fc);

    /* call user-defined grammar pre-condition callback function */
    // if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc)) // not into
    //   fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_PRE, fc->aux_grammar->data);

    if (!fill_arrays(fc)) {
// #ifdef SUN4
//       standard_arithmetic();
// #elif defined(HP9)
//       fpsetfastmode(0);
// #endif
      return dG;
    }

    if (md->circ)
      /* do post processing step for circular RNAs */
      postprocess_circular(fc);

    /* call user-defined grammar post-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc))
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_POST, fc->aux_grammar->data);

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

// #ifdef SUN4
//     standard_arithmetic();
// #elif defined(HP9)
//     fpsetfastmode(0);
// #endif
  }

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
  // printf("%s*",structure);

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




int main(){
  struct  options             *opt;
  char                        *sequence, *pairing_propensity;
  int                         n;
  double                      min_en = -20.;
  double                      kT;
  vrna_dimer_pf_t             AB;
  init_default_options(&opt);
  vrna_fold_compound_t *vc = vrna_fold_compound(sequence,
                                                &(opt->md),
                                                VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
  n = vc->length;
  vrna_exp_params_rescale(vc, &min_en);
  kT = vc->exp_params->kT / 1000.;
  AB = vrna_pf_dimer(vc, pairing_propensity);
  return 1;
}