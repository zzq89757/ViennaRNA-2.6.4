/**********************************************/
/* BEGIN interface for model details          */
/**********************************************/


%{
#include <sstream>
%}

/* scripting language access through 'md' instead of 'vrna_md_t' */
%rename (md) vrna_md_t;

%nodefaultctor vrna_md_t;
%nodefaultdtor vrna_md_t;

/* hide all attributes, except for some trivial ones */
typedef struct {
  double  temperature;
  double  betaScale;
  int     pf_smooth;
  int     dangles;
  int     special_hp;
  int     noLP;
  int     noGU;
  int     noGUclosure;
  int     logML;
  int     circ;
  int     gquad;
  int     uniq_ML;
  int     energy_set;
  int     backtrack;
  char    backtrack_type;
  int     compute_bpp;
  char    nonstandards[64];
  int     max_bp_span;
  int     min_loop_size;
  int     window_size;
  int     oldAliEn;
  int     ribo;
  double  cv_fact;
  double  nc_fact;
  double  sfact;
  const int     rtype[8];
  const short   alias[MAXALPHA+1];
  const int const pair[MAXALPHA+1][MAXALPHA+1];
  double  salt;
  int saltMLLower;
  int saltMLUpper;
  int saltDPXInit;
  float   saltDPXInitFact;
  float   helical_rise;
  float   backbone_length;
} vrna_md_t;


/* make a nice object oriented interface to vrna_md_t */
%extend vrna_md_t {

#ifdef SWIGPYTHON
%feature("autodoc")vrna_md_t::vrna_md_t;
%feature("kwargs")vrna_md_t::vrna_md_t;
#endif

  /*  Default constructor */
  vrna_md_t(
    const double  temperature     = vrna_md_defaults_temperature_get(),
    const double  betaScale       = vrna_md_defaults_betaScale_get(),
    const int     pf_smooth       = vrna_md_defaults_pf_smooth_get(),
    const int     dangles         = vrna_md_defaults_dangles_get(),
    const int     special_hp      = vrna_md_defaults_special_hp_get(),
    const int     noLP            = vrna_md_defaults_noLP_get(),
    const int     noGU            = vrna_md_defaults_noGU_get(),
    const int     noGUclosure     = vrna_md_defaults_noGUclosure_get(),
    const int     logML           = vrna_md_defaults_logML_get(),
    const int     circ            = vrna_md_defaults_circ_get(),
    const int     gquad           = vrna_md_defaults_gquad_get(),
    const int     uniq_ML         = vrna_md_defaults_uniq_ML_get(),
    const int     energy_set      = vrna_md_defaults_energy_set_get(),
    const int     backtrack       = vrna_md_defaults_backtrack_get(),
    const char    backtrack_type  = vrna_md_defaults_backtrack_type_get(),
    const int     compute_bpp     = vrna_md_defaults_compute_bpp_get(),
    const int     max_bp_span     = vrna_md_defaults_max_bp_span_get(),
    const int     min_loop_size   = vrna_md_defaults_min_loop_size_get(),
    const int     window_size     = vrna_md_defaults_window_size_get(),
    const int     oldAliEn        = vrna_md_defaults_oldAliEn_get(),
    const int     ribo            = vrna_md_defaults_ribo_get(),
    const double  cv_fact         = vrna_md_defaults_cv_fact_get(),
    const double  nc_fact         = vrna_md_defaults_nc_fact_get(),
    const double  sfact           = vrna_md_defaults_sfact_get(),
    const double  salt            = vrna_md_defaults_salt_get(),
    const int     saltMLLower     = vrna_md_defaults_saltMLLower_get(),
    const int     saltMLUpper     = vrna_md_defaults_saltMLUpper_get(),
    const int     saltDPXInit     = vrna_md_defaults_saltDPXInit_get(),
    const float   saltDPXInitFact = vrna_md_defaults_saltDPXInitFact_get(),
    const float   helical_rise    = vrna_md_defaults_helical_rise_get(),
    const float   backbone_length = vrna_md_defaults_backbone_length_get())
  {
    vrna_md_t *md       = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    md->temperature     = temperature;
    md->betaScale       = betaScale;
    md->pf_smooth       = pf_smooth;
    md->dangles         = dangles;
    md->special_hp      = special_hp;
    md->noLP            = noLP;
    md->noGU            = noGU;
    md->noGUclosure     = noGUclosure;
    md->logML           = logML;
    md->circ            = circ;
    md->gquad           = gquad;
    md->uniq_ML         = uniq_ML;
    md->energy_set      = energy_set;
    md->backtrack       = backtrack;
    md->backtrack_type  = backtrack_type;
    md->compute_bpp     = compute_bpp;
    md->max_bp_span     = max_bp_span;
    md->min_loop_size   = min_loop_size;
    md->window_size     = window_size;
    md->oldAliEn        = oldAliEn;
    md->ribo            = ribo;
    md->cv_fact         = cv_fact;
    md->nc_fact         = nc_fact;
    md->sfact           = sfact;
    md->salt            = salt;
    md->saltMLLower     = saltMLLower;
    md->saltMLUpper     = saltMLUpper;
    md->saltDPXInit     = saltDPXInit;
    md->saltDPXInitFact = saltDPXInitFact;
    md->helical_rise    = helical_rise;
    md->backbone_length = backbone_length;

    vrna_md_update(md);

    return md;
  }

  ~vrna_md_t()
  {
    free($self);
  }

  void
  reset()
  {
    vrna_md_set_default($self);
  }

  void
  set_from_globals()
  {
    set_model_details($self);
  }

  char *
  option_string()
  {
    return vrna_md_option_string($self);
  }

#ifdef SWIGPYTHON
  std::string
  __str__()
  {
    std::ostringstream out;
    out << "{ temperature: " << $self->temperature ;
    out << ", dangles: " << $self->dangles;
    out << ", betaScale: " << $self->betaScale ;
    out << ", pf_smooth: " << $self->pf_smooth ;
    out << ", special_hp: " << $self->special_hp ;
    out << ", noLP: " << $self->noLP ;
    out << ", noGU: " << $self->noGU ;
    out << ", noGUclosure: " << $self->noGUclosure ;
    out << ", logML: " << $self->logML ;
    out << ", circ: " << $self->circ ;
    out << ", gquad: " << $self->gquad ;
    out << ", uniq_ML: " << $self->uniq_ML ;
    out << ", energy_set: " << $self->energy_set  ;
    out << ", backtrack: " << $self->backtrack ;
    out << ", backtrack_type: " << $self->backtrack_type ;
    out << ", compute_bpp: " << $self->compute_bpp ;
    out << ", max_bp_span: " << $self->max_bp_span ;
    out << ", min_loop_size: " << $self->min_loop_size;
    out << ", window_size: " << $self->window_size ;
    out << ", oldAliEn: " << $self->oldAliEn;
    out << ", ribo: " << $self->ribo;
    out << ", cv_fact: " << $self->cv_fact ;
    out << ", nc_fact: " << $self->nc_fact ;
    out << ", sfact: " << $self->sfact ;
    out << ", salt: " << $self->salt ;
    out << ", saltMLLower: " << $self->saltMLLower ;
    out << ", saltMLUpper: " << $self->saltMLUpper ;
    out << ", saltDPXInit: " << $self->saltDPXInit ;
    out << ", saltDPXInitFact: " << $self->saltDPXInitFact ;
    out << ", helical_rise: " << $self->helical_rise ;
    out << ", backbone_length: " << $self->backbone_length ;
    out << " }";

    return std::string(out.str());
  }

%pythoncode %{
def __repr__(self):
    # reformat string representation (self.__str__()) to something
    # that looks like a constructor argument list
    strthis = self.__str__().replace(": ", "=").replace("{ ", "").replace(" }", "")
    return  "%s.%s(%s)" % (self.__class__.__module__, self.__class__.__name__, strthis) 
%}
#endif

}


/*
 * Hide the entire interface
 */
%ignore set_model_details;
%ignore option_string;

/*
 * Include typemaps for augmenting reading and assignment of
 * global variables with their corresponding getter/setter
 * functions in the library
 */
%include globals-md.i

/* these two global variables will be handled by getter/setter methods */
%ignore helical_rise;
%ignore backbone_length;

/*
 * The list of available global variables in the scripting
 * interface. Access to these variables is defined in the
 * typemap definition file above
 */
extern double temperature;
extern int    dangles;
extern double betaScale;
extern int    pf_smooth;
extern int    tetra_loop;     /* this is an alias of special_hp */
extern int    special_hp;
extern int    noLonelyPairs;  /* this is an alias of noLP */
extern int    noLP;
extern int    noGU;
extern int    no_closingGU;    /* this is an alias of noGUclosure */
extern int    noGUclosure;
extern int    logML;
extern int    circ;
extern int    gquad;
extern int    uniq_ML;
extern int    energy_set;
extern int    backtrack;
extern char   backtrack_type; /* still needs implementation */
extern int    do_backtrack;   /* this is an alias of compute_bpp */
extern int    compute_bpp;
extern int    max_bp_span;
extern int    min_loop_size;
extern int    window_size;
extern int    oldAliEn;
extern int    ribo;
extern double cv_fact;
extern double nc_fact;
extern double sfact;
extern double salt;
extern int    saltDPXInit;
extern double saltDPXInitFact;
extern double helical_rise;
extern double backbone_length;


%constant double  MODEL_DEFAULT_TEMPERATURE         = VRNA_MODEL_DEFAULT_TEMPERATURE;
%constant double  MODEL_DEFAULT_PF_SCALE          = VRNA_MODEL_DEFAULT_PF_SCALE;
%constant double  MODEL_DEFAULT_BETA_SCALE        = VRNA_MODEL_DEFAULT_BETA_SCALE;
%constant int     MODEL_DEFAULT_DANGLES           = VRNA_MODEL_DEFAULT_DANGLES;
%constant int     MODEL_DEFAULT_SPECIAL_HP        = VRNA_MODEL_DEFAULT_SPECIAL_HP;
%constant int     MODEL_DEFAULT_NO_LP             = VRNA_MODEL_DEFAULT_NO_LP;
%constant int     MODEL_DEFAULT_NO_GU             = VRNA_MODEL_DEFAULT_NO_GU;
%constant int     MODEL_DEFAULT_NO_GU_CLOSURE     = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
%constant int     MODEL_DEFAULT_CIRC              = VRNA_MODEL_DEFAULT_CIRC;
%constant int     MODEL_DEFAULT_GQUAD             = VRNA_MODEL_DEFAULT_GQUAD;
%constant int     MODEL_DEFAULT_UNIQ_ML           = VRNA_MODEL_DEFAULT_UNIQ_ML;
%constant int     MODEL_DEFAULT_ENERGY_SET        = VRNA_MODEL_DEFAULT_ENERGY_SET;
%constant int     MODEL_DEFAULT_BACKTRACK         = VRNA_MODEL_DEFAULT_BACKTRACK;
%constant char    MODEL_DEFAULT_BACKTRACK_TYPE    = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE;
%constant int     MODEL_DEFAULT_COMPUTE_BPP       = VRNA_MODEL_DEFAULT_COMPUTE_BPP;
%constant int     MODEL_DEFAULT_MAX_BP_SPAN       = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
%constant int     MODEL_DEFAULT_WINDOW_SIZE       = VRNA_MODEL_DEFAULT_WINDOW_SIZE;
%constant int     MODEL_DEFAULT_LOG_ML            = VRNA_MODEL_DEFAULT_LOG_ML;
%constant int     MODEL_DEFAULT_ALI_OLD_EN        = VRNA_MODEL_DEFAULT_ALI_OLD_EN;
%constant int     MODEL_DEFAULT_ALI_RIBO          = VRNA_MODEL_DEFAULT_ALI_RIBO;
%constant double  MODEL_DEFAULT_ALI_CV_FACT       = VRNA_MODEL_DEFAULT_ALI_CV_FACT;
%constant double  MODEL_DEFAULT_ALI_NC_FACT       = VRNA_MODEL_DEFAULT_ALI_NC_FACT;
%constant int     MODEL_DEFAULT_PF_SMOOTH         = VRNA_MODEL_DEFAULT_PF_SMOOTH;
%constant double  MODEL_DEFAULT_SALT              = VRNA_MODEL_DEFAULT_SALT;
%constant int     MODEL_DEFAULT_SALT_MLLOWER      = VRNA_MODEL_DEFAULT_SALT_MLLOWER;
%constant int     MODEL_DEFAULT_SALT_MLUPPER      = VRNA_MODEL_DEFAULT_SALT_MLUPPER;
%constant int     MODEL_DEFAULT_SALT_DPXINIT      = VRNA_MODEL_DEFAULT_SALT_DPXINIT;
%constant double  MODEL_DEFAULT_SALT_DPXINIT_FACT = VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT;
%constant double  MODEL_DEFAULT_HELICAL_RISE      = VRNA_MODEL_DEFAULT_HELICAL_RISE;
%constant double  MODEL_DEFAULT_BACKBONE_LENGTH   = VRNA_MODEL_DEFAULT_BACKBONE_LENGTH;
%constant double  MODEL_SALT_DPXINIT_FACT_RNA     = VRNA_MODEL_SALT_DPXINIT_FACT_RNA;
%constant double  MODEL_SALT_DPXINIT_FACT_DNA     = VRNA_MODEL_SALT_DPXINIT_FACT_DNA;
%constant double  MODEL_HELICAL_RISE_RNA          = VRNA_MODEL_HELICAL_RISE_RNA;
%constant double  MODEL_HELICAL_RISE_DNA          = VRNA_MODEL_HELICAL_RISE_DNA;
%constant double  MODEL_BACKBONE_LENGTH_RNA       = VRNA_MODEL_BACKBONE_LENGTH_RNA;
%constant double  MODEL_BACKBONE_LENGTH_DNA       = VRNA_MODEL_BACKBONE_LENGTH_DNA;

%include <ViennaRNA/model.h>

