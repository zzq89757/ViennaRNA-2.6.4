



class vrna_fc_type_e:
    def __init__(self) -> None:
        pass


class vrna_fold_compound_t:
    def __init__(self) -> None:
        pass
    
    def init_fc_single(self):
        self.type = VRNA_FC_TYPE_SINGLE
        
    def nullify(self):
        self.length = 0
        self.length            = 0
        self.strands           = 0
        self.cutpoint          = -1
        self.strand_number     = None
        self.strand_order      = None
        self.strand_order_uniq = None
        self.strand_start      = None
        self.strand_end        = None
        self.nucleotides       = None
        self.alignment         = None

        self.hc            = None
        self.matrices      = None
        self.exp_matrices  = None
        self.params        = None
        self.exp_params    = None
        self.iindx         = None
        self.jindx         = None

        self.stat_cb       = None
        self.auxdata       = None
        self.free_auxdata  = None

        self.domains_struc = None
        self.domains_up    = None
        self.aux_grammar   = None
        
        if self.type == VRNA_FC_TYPE_SINGLE:
            self.sequence            = None
            self.sequence_encoding   = None
            self.encoding5           = None
            self.encoding3           = None
            self.sequence_encoding2  = None
            self.ptype               = None
            self.ptype_pf_compat     = None
            self.sc                  = None
        elif self.type == VRNA_FC_TYPE_COMPARATIVE:
            self.sequences         = None
            self.n_seq             = 0
            self.cons_seq          = None
            self.S_cons            = None
            self.S                 = None
            self.S5                = None
            self.S3                = None
            self.Ss                = None
            self.a2s               = None
            self.pscore            = None
            self.pscore_local      = None
            self.pscore_pf_compat  = None
            self.scs               = None
            self.oldAliEn          = 0
        self.maxD1         = 0
        self.maxD2         = 0
        self.reference_pt1 = None
        self.reference_pt2 = None
        self.referenceBPs1 = None
        self.referenceBPs2 = None
        self.bpdist        = None
        self.mm1           = None
        self.mm2           = None

        self.window_size = -1
        self.ptype_local = None
        
        if VRNA_WITH_SVM:self.zscore_data = None
    
    def vrna_fold_compound(self, sequence, md_p, options):
        if sequence is None:
            return None

        # sanity check
        length = len(sequence)
        if length == 0:
            print("vrna_fold_compound@data_structures.c: sequence length must be greater than 0")
            return None

        if length > vrna_sequence_length_max(options):
            print(f"vrna_fold_compound@data_structures.c: sequence length of {length} exceeds addressable range")
            return None

        fc = init_fc_single()

        fc.length = length
        fc.sequence = sequence  # In Python, no need for strdup

        aux_options = 0

        # get a copy of the model details
        if md_p:
            md = md_p
        else:
            md = vrna_md_t()
            vrna_md_set_default(md)

        # now for the energy parameters
        add_params(fc, md, options)

        sanitize_bp_span(fc, options)

        if options & VRNA_OPTION_WINDOW:
            set_fold_compound(fc, options, aux_options)

            if not (options & VRNA_OPTION_EVAL_ONLY):
                # add minimal hard constraint data structure
                vrna_hc_init_window(fc)

                # add DP matrices
                vrna_mx_add(fc, VRNA_MX_WINDOW, options)
        else:
            # regular global structure prediction
            aux_options |= WITH_PTYPE

            if options & VRNA_OPTION_PF:
                aux_options |= WITH_PTYPE_COMPAT

            set_fold_compound(fc, options, aux_options)

            if not (options & VRNA_OPTION_EVAL_ONLY):
                # add default hard constraints
                vrna_hc_init(fc)

                # add DP matrices (if required)
                vrna_mx_add(fc, VRNA_MX_DEFAULT, options)

        return fc
    
    
    def do_partfunc(self, mystring, length, Switch, tpr, mfpl, kT, opt):
        




if __name__ == "__main__":
    VRNA_OPTION_MFE = 1 << 0
    VRNA_OPTION_PF = 1 << 1
    VRNA_OPTION_HYBRID = 1 << 2
    vc = vrna_fold_compound(Newstring, md, VRNA_OPTION_MFE | VRNA_OPTION_PF | VRNA_OPTION_HYBRID)