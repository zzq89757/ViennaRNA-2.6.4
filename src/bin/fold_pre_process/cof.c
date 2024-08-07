// header include start--------------------------------
#include "cof.h"
#include "matrices/all_mat.h"
#include "dup.h"


// structure definition start -------------------------------


//function definition start --------------------------------


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

      case VRNA_FC_TYPE_COMPARATIVE:
        fc->sequences         = NULL;
        fc->n_seq             = 0;
        fc->cons_seq          = NULL;
        fc->S_cons            = NULL;
        fc->S                 = NULL;
        fc->S5                = NULL;
        fc->S3                = NULL;
        fc->Ss                = NULL;
        fc->a2s               = NULL;
        fc->pscore            = NULL;
        fc->pscore_local      = NULL;
        fc->pscore_pf_compat  = NULL;
        fc->scs               = NULL;
        fc->oldAliEn          = 0;

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

// her
// PUBLIC int
// vrna_mx_mfe_add(vrna_fold_compound_t  *vc,
//                 vrna_mx_type_e        mx_type,
//                 unsigned int          options)
// {
//   unsigned int mx_alloc_vector;

//   if (vc->params) {
//     options |= VRNA_OPTION_MFE;

//     mx_alloc_vector = get_mx_alloc_vector(vc,
//                                           mx_type,
//                                           options);
//     vrna_mx_mfe_free(vc);
//     return add_mfe_matrices(vc, mx_type, mx_alloc_vector);
//   }

//   return 0;
// }

PUBLIC int
vrna_mx_add(vrna_fold_compound_t  *vc,
            vrna_mx_type_e        mx_type,
            unsigned int          options)
{
  int ret;

  ret = 1;

  if (options & VRNA_OPTION_MFE)
    printf("into vrna_mx_mfe_add \n");
    ret &= vrna_mx_mfe_add(vc, mx_type, options);

  if (options & VRNA_OPTION_PF)
    printf("into vrna_mx_pf_add \n");
    ret &= vrna_mx_pf_add(vc, mx_type, options);

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


struct hc_ext_def_dat {
  unsigned int              n;
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};


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
    //   vrna_message_warning("vrna_pf@part_func.c: Failed to prepare vrna_fold_compound");
    //   return dG;
    // }

    n         = fc->length;
    params    = fc->exp_params;
    matrices  = fc->exp_matrices;
    md        = &(params->model_details);

#ifdef _OPENMP
    /* Explicitly turn off dynamic threads */
    omp_set_dynamic(0);
#endif

#ifdef SUN4
    nonstandard_arithmetic();
#elif defined(HP9)
    fpsetfastmode(1);
#endif

    /* call user-defined recursion status callback function */
    if (fc->stat_cb) // not into
      fc->stat_cb(VRNA_STATUS_PF_PRE, fc->auxdata);

    /* for now, multi-strand folding is implemented as additional grammar rule */
    if (fc->strands > 1)
      vrna_pf_multifold_prepare(fc);

    /* call user-defined grammar pre-condition callback function */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc)) // not into
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_PRE, fc->aux_grammar->data);

    if (!fill_arrays(fc)) {
#ifdef SUN4
      standard_arithmetic();
#elif defined(HP9)
      fpsetfastmode(0);
#endif
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

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      dG /= fc->n_seq;

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

#ifdef SUN4
    standard_arithmetic();
#elif defined(HP9)
    fpsetfastmode(0);
#endif
  }

  return dG;
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



static void 
process_record(struct record_data *record){
  char                  *mfe_structure, *sequence, **rec_rest;
  unsigned int          n, i;
  double                min_en, kT, *concentrations;
  vrna_ep_t             *prAB, *prAA, *prBB, *prA, *prB, *mfAB, *mfAA, *mfBB, *mfA, *mfB;
  struct options        *opt;
  struct output_stream  *o_stream;
  size_t                **mod_positions;
  size_t                mod_param_sets;

  mfAB            = mfAA = mfBB = mfA = mfB = NULL;
  prAB            = prAA = prBB = prA = prB = NULL;
  concentrations  = NULL;
  init_default_options(&opt);
  o_stream        = (struct output_stream *)vrna_alloc(sizeof(struct output_stream));
  sequence        = strdup(record->sequence);
  // rec_rest        = record->rest;
  vrna_fold_compound_t *vc = vrna_fold_compound(sequence,
                                                &(opt->md),
                                                VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
  n = vc->length;
  // min_en  = vrna_mfe_dimer(vc, mfe_structure);

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
  vrna_exp_params_rescale(vc, &min_en);
  kT = vc->exp_params->kT / 1000.;
  AB = AA = BB = vrna_pf_dimer(vc, pairing_propensity); /* exp_matrices changed*/


}




int main(){
  return 0;
}