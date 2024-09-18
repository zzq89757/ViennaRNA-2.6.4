# define INLINE inline
# define PUBLIC

typedef struct vrna_md_s vrna_md_t;
typedef struct  vrna_param_s vrna_param_t;
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


PUBLIC int
vrna_salt_loop_int(int L, double rho, double T, double backbonelen);

#define MIN2(A,B) ((A) < (B) ? (A) : (B))
#define MAX2(A,B) ((A) > (B) ? (A) : (B))
#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)

PUBLIC double
vrna_salt_loop(int L, double rho, double T, double backbonelen);

PUBLIC int
vrna_salt_stack(double rho, double T, double hrise);


PUBLIC void
vrna_salt_ml(double saltLoop[], int lower, int upper, int* m, int* b);
PUBLIC int
vrna_salt_duplex_init(vrna_md_t *md_p);