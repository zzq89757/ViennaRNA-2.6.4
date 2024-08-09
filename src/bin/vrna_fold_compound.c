
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
