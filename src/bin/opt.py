from remove import vrna_md_t, Duplex


class options:
    def __init__(self) -> None:
        self.filename_full = self.doT = self.doC = 0
        self.md = vrna_md_t()

def init_default_options(opt:options) -> None:
    opt.filename_full  = 0
    opt.filename_delim = None
    opt.pf             = 1
    opt.doT            = 0 # compute dimer free energies etc. */
    opt.noPS           = 1
    opt.noconv         = 0
    opt.centroid       = 0  # off by default due to historical reasons */
    opt.MEA            = 0
    opt.MEAgamma       = 1.
    opt.bppmThreshold  = 1e-5
    opt.verbose        = 0
    opt.commands       = None
    opt.id_control     = None
    Duplex.set_model_details(opt.md)

    opt.doC                = 0 # toggle to compute concentrations */
    opt.concentration_file = None

    opt.constraint_file      = None
    opt.constraint_batch     = 0
    opt.constraint_enforce   = 0
    opt.constraint_canonical = 0

    opt.shape            = 0
    opt.shape_file       = None
    opt.shape_method     = None
    opt.shape_conversion = None

    opt.mod_params         = None

    opt.csv_output       = 0    # flag indicating whether we produce one-line outputs, a.k.a. CSV */
    opt.csv_header       = 1    # print header for one-line output */
    opt.csv_output_delim = ','  # delimiting character for one-line output */

    opt.jobs               = 1
    opt.keep_order         = 1
    opt.next_record_number = 0
    opt.output_queue       = None