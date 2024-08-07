# define INLINE inline

PUBLIC void
set_model_details(vrna_md_t *md);

PUBLIC void *
vrna_alloc(unsigned size);

PUBLIC void
vrna_md_set_default(vrna_md_t *md);

PUBLIC vrna_md_t *
vrna_md_copy(vrna_md_t        *md_to,
             const vrna_md_t  *md_from);


PUBLIC vrna_param_t *
vrna_params(vrna_md_t *md);

PUBLIC int
vrna_nucleotide_encode(char       c,
                       vrna_md_t  *md);
