// header include start--------------------------------
#include "cof.h"
#include "matrices/all_mat.h"



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


PUBLIC void
vrna_exp_params_rescale(vrna_fold_compound_t  *vc,
                        double                *mfe)
{
  vrna_exp_param_t  *pf;
  double            e_per_nt, kT;
  vrna_md_t         *md;

  if (vc) {
    if (!vc->exp_params) {
      switch (vc->type) {
        case VRNA_FC_TYPE_SINGLE:
          vc->exp_params = vrna_exp_params(&(vc->params->model_details));
          break;
        case VRNA_FC_TYPE_COMPARATIVE:
          vc->exp_params = vrna_exp_params_comparative(vc->n_seq, &(vc->params->model_details));
          break;
      }
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