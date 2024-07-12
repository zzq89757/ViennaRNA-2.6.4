



class vrna_fc_type_e:
    def __init__(self) -> None:
        pass


class vrna_fold_compound_t:
    def __init__(self) -> None:
        pass
    
    def init_fc_single(self):
        self.type = VRNA_FC_TYPE_SINGLE
        
    
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




if __name__ == "__main__":
    pass