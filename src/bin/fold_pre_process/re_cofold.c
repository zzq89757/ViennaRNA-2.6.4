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


PUBLIC unsigned int
vrna_get_ptype(int  ij,
               char *ptype)
{
  unsigned int tt = (unsigned int)ptype[ij];

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


PUBLIC unsigned int
vrna_get_ptype_md(int       i,
                  int       j,
                  vrna_md_t *md)
{
  unsigned int tt = (unsigned int)md->pair[i][j];

  return (tt == 0) ? 7 : tt;
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