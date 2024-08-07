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

struct vrna_mx_pf_aux_el_s {
  FLT_OR_DBL  *qq;
  FLT_OR_DBL  *qq1;

  int         qqu_size;
  FLT_OR_DBL  **qqu;
};


typedef struct vrna_mx_pf_aux_el_s *vrna_mx_pf_aux_el_t;

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


struct hc_hp_def_dat {
  int                       n;
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};


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


typedef unsigned char (*eval_hc)(int   i,
                                int   j,
                                int   k,
                                int   l,
                                void  *data);

struct hc_int_def_dat {
  unsigned char             *mx;
  unsigned char             **mx_local;
  unsigned int              *sn;
  unsigned int              n;
  int                       *up;

  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
};


struct sc_int_exp_dat;

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



PRIVATE INLINE
FLT_OR_DBL
exp_E_GQuad_IntLoop(int               i,
                    int               j,
                    int               type,
                    short             *S,
                    FLT_OR_DBL        *G,
                    FLT_OR_DBL        *scale,
                    int               *index,
                    vrna_exp_param_t  *pf)
{
  int         k, l, minl, maxl, u, r;
  FLT_OR_DBL  q, qe;
  double      *expintern;
  short       si, sj;

  q         = 0;
  si        = S[i + 1];
  sj        = S[j - 1];
  qe        = (FLT_OR_DBL)pf->expmismatchI[type][si][sj];
  expintern = &(pf->expinternal[0]);

  if (type > 2)
    qe *= (FLT_OR_DBL)pf->expTermAU;

  k = i + 1;
  if (S[k] == 3) {
    if (k < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minl  = j - MAXLOOP - 1;
      u     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minl  = MAX2(u, minl);
      u     = j - 3;
      maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxl  = MIN2(u, maxl);
      for (l = minl; l < maxl; l++) {
        if (S[l] != 3)
          continue;

        if (G[index[k] - l] == 0.)
          continue;

        q += qe
             * G[index[k] - l]
             * (FLT_OR_DBL)expintern[j - l - 1]
             * scale[j - l + 1];
      }
    }
  }

  for (k = i + 2;
       k <= j - VRNA_GQUAD_MIN_BOX_SIZE;
       k++) {
    u = k - i - 1;
    if (u > MAXLOOP)
      break;

    if (S[k] != 3)
      continue;

    minl  = j - i + k - MAXLOOP - 2;
    r     = k + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minl  = MAX2(r, minl);
    maxl  = k + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    r     = j - 1;
    maxl  = MIN2(r, maxl);
    for (l = minl; l < maxl; l++) {
      if (S[l] != 3)
        continue;

      if (G[index[k] - l] == 0.)
        continue;

      q += qe
           * G[index[k] - l]
           * (FLT_OR_DBL)expintern[u + j - l - 1]
           * scale[u + j - l + 1];
    }
  }

  l = j - 1;
  if (S[l] == 3)
    for (k = i + 4; k <= j - VRNA_GQUAD_MIN_BOX_SIZE; k++) {
      u = k - i - 1;
      if (u > MAXLOOP)
        break;

      if (S[k] != 3)
        continue;

      if (G[index[k] - l] == 0.)
        continue;

      q += qe
           * G[index[k] - l]
           * (FLT_OR_DBL)expintern[u]
           * scale[u + 2];
    }

  return q;
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

      if ((with_gquad) && (!noclose)) {
        switch (fc->type) {
          case VRNA_FC_TYPE_SINGLE:
            if (sliding_window) {
              /* no G-Quadruplex support for sliding window partition function yet! */
            } else if (sn[j] == sn[i]) {
              qbt1 += exp_E_GQuad_IntLoop(i, j, type, S1, G, scale, my_iindx, pf_params);
            }

            break;

        }
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


struct hc_mb_def_dat {
  unsigned char             *mx;
  unsigned char             **mx_window;
  unsigned int              *sn;
  unsigned int              n;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_hc_eval_f hc_f;
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


      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          tt      = vrna_get_ptype_md(SS[s][j], SS[s][i], md);
          qqqmmm  *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params);
        }
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

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE) {
      jindx         = fc->jindx;
      pscore        = fc->pscore;
      kTn           = fc->exp_params->kT / 10.;  /* kT in cal/mol */
      contribution  *= exp(pscore[jindx[j] + i] / kTn);
    }
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

      case VRNA_FC_TYPE_COMPARATIVE:
        q_temp = 1.;
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          q_temp  *= exp_E_MLstem(type,
                                  ((i > 1) || circular) ? S5[s][i] : -1,
                                  ((j < n) || circular) ? S3[s][j] : -1,
                                  pf_params);
        }
        qbt1 *= q_temp;
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

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = fc->n_seq;
        S     = fc->S;
        S5    = fc->S5;
        S3    = fc->S3;
        a2s   = fc->a2s;
        for (s = 0; s < n_seq; s++) {
          type    = vrna_get_ptype_md(S[s][i], S[s][j], md);
          q_temp  *= vrna_exp_E_ext_stem(type,
                                         ((a2s[s][i] > 1) || circular) ? S5[s][i] : -1,
                                         ((a2s[s][j] < a2s[s][n]) || circular) ? S3[s][j] : -1,
                                         pf_params);
        }
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

  /* no G-Quadruplexes for comparative partition function (yet) */
  // if (with_gquad) { // not into
  //   free(fc->exp_matrices->G);
  //   fc->exp_matrices->G = NULL;

  //   switch (fc->type) {
  //     case VRNA_FC_TYPE_SINGLE:
  //       fc->exp_matrices->G = get_gquad_pf_matrix(fc->sequence_encoding2,
  //                                                 fc->exp_matrices->scale,
  //                                                 fc->exp_params);
  //       break;
  //   }
  // }

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

      // if (q[ij] >= max_real) { // not into
      //   printf("overflow while computing partition function for segment q[%d,%d]\n"
      //                        "use larger pf_scale", i, j);

      //   vrna_exp_E_ml_fast_free(aux_mx_ml);
      //   vrna_exp_E_ext_fast_free(aux_mx_el);

      //   return 0; /* failure */
      // }
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
  scs           = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->scs : NULL;
  a2s           = (fc->type == VRNA_FC_TYPE_COMPARATIVE) ? fc->a2s : NULL;
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
      if (val > max) { \
        printf("vrna_search_BMH: " \
                             "haystack value %d at hit %d " \
                             "out of bad character table range [%d : %d]\n" \
                             "Aborting search...", \
                             (shift + needle_size - 1) % haystack_size, \
                             val, \
                             0, \
                             max); \
        return NULL; \
      } \
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

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      compute_bpp_int = &compute_bpp_internal;
      compute_bpp_mul = &compute_bpp_multibranch;
    } else {
      compute_bpp_int = &compute_bpp_internal_comparative;
      compute_bpp_mul = &compute_bpp_multibranch_comparative;
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
    // printf("probs is :%f", probs);
    /* 1. external loop pairs, i.e. pairs not enclosed by any other pair (or external loop for circular RNAs) */
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
      if (with_ud_outside) {
        /*
         *  The above recursions only deal with base pairs, and how they might be
         *  enclosed by other pairs. However, for unstructrued domains, we have
         *  unpaired stretches, and require information about how these are enclosed
         *  by base pairs.
         */

        /* 1. Exterior loops */
        ud_outside_ext_loops(vc);

        /* 2. Hairpin loops */
        ud_outside_hp_loops(vc);

        /* 3. Interior loops */
        ud_outside_int_loops(vc);

        /* 4. Multi branch loops */
        ud_outside_mb_loops(vc);
      }

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
            if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
              probs[ij] *= exp(-pscore[jindx[j] + i] / kTn);
          } else if (G[ij] > 0.) {
            probs[ij] += q1k[i - 1] *
                         G[ij] *
                         qln[j + 1] /
                         q1k[n];
          }
        } else {
          if (qb[ij] > 0.) {
            probs[ij] *= qb[ij];

            if (vc->type == VRNA_FC_TYPE_COMPARATIVE)
              probs[ij] *= exp(-pscore[jindx[j] + i] / kTn);
          }
        }
        // printf("%d-", probs[ij]);
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
      vrna_message_warning("%d overflows occurred while backtracking;\n"
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
    vrna_message_warning("bppm calculations have to be done after calling forward recursion");
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
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_proc)) // not into
      fc->aux_grammar->cb_proc(fc, VRNA_STATUS_PF_PRE, fc->aux_grammar->data);

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