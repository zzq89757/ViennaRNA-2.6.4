import numpy as np

# constant define
VRNA_FC_TYPE_COMPARATIVE = 1
VRNA_CONSTRAINT_CONTEXT_HP_LOOP = 0x02
VRNA_DECOMP_PAIR_HP = 1
VRNA_FC_TYPE_SINGLE = 0
VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP = 4
VRNA_DECOMP_PAIR_ML = 3
VRNA_DECOMP_ML_ML = 7
VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP = 8
VRNA_UNSTRUCTURED_DOMAIN_MOTIF = 16
VRNA_DECOMP_ML_STEM = 6
VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC = 0x08
MAXLOOP = 30
VRNA_DECOMP_PAIR_IL = 2
VRNA_CONSTRAINT_CONTEXT_MB_LOOP = 0x10
VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC = 0x20
VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP = 1
VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP = 2
VRNA_HC_WINDOW = 1
VRNA_DECOMP_EXT_STEM_OUTSIDE = 17

# def generate_probs(vc, matrices):
#     my_iindx = vc.iindx
#     n = vc.length
#     probs = matrices.probs
    
    
#     for i in range(1, n + 1):
#         probs[my_iindx[i] - i] = 0

class vrna_fc_type_e:
    def __init__(self) -> None:
        ...


class vrna_fold_compound_t:
    def __init__(self) -> None:
        # Common data fields
        self.type = vrna_fc_type_e  # vrna_fc_type_e
        self.length = 0  # unsigned int
        self.strand_number = []  # unsigned int array
        self.strand_order = []  # unsigned int array
        self.strand_order_uniq = []  # unsigned int array
        self.strand_start = []  # unsigned int array
        self.strand_end = []  # unsigned int array
        self.strands = 0  # unsigned int
        self.nucleotides = vrna_seq_t  # vrna_seq_t array
        self.alignment = vrna_msa_t  # vrna_msa_t array
        self.hc = vrna_hc_t  # vrna_hc_t
        self.matrices = vrna_mx_mfe_t  # vrna_mx_mfe_t
        self.exp_matrices = vrna_mx_pf_t  # vrna_mx_pf_t
        self.params = vrna_param_t  # vrna_param_t
        self.exp_params = vrna_exp_param_t  # vrna_exp_param_t
        self.iindx = 0  # int array
        self.jindx = 0  # int array

        # User-defined data fields
        self.stat_cb = vrna_recursion_status_f  # vrna_recursion_status_f
        self.auxdata = 0  # void pointer
        self.free_auxdata = vrna_auxdata_free_f  # vrna_auxdata_free_f

        # Secondary Structure Decomposition (grammar) related data fields
        self.domains_struc = vrna_sd_t  # vrna_sd_t
        self.domains_up = vrna_ud_t  # vrna_ud_t
        self.aux_grammar = vrna_gr_aux_t  # vrna_gr_aux_t

        # Data fields available for single/hybrid structure prediction
        self.sequence = ''  # char array
        self.sequence_encoding = 0  # short array
        self.encoding5 = 0  # short array
        self.encoding3 = 0  # short array
        self.sequence_encoding2 = 0  # short array
        self.ptype = 'ptype'  # char array
        self.ptype_pf_compat = 'ptype_pf_compat'  # char array
        self.sc = vrna_sc_t  # vrna_sc_t

        # Data fields for consensus structure prediction
        self.n_seq = 0  # unsigned int
        self.cons_seq = 'cons_seq'  # char array
        self.S_cons = 0  # short array
        self.S = 0  # short array of arrays
        self.S5 = 0  # short array of arrays
        self.S3 = 0  # short array of arrays
        self.Ss = 'Ss'  # char array of arrays
        self.a2s = 0  # unsigned int array of arrays
        self.pscore = 0  # int array
        self.pscore_local = 0  # int array of arrays
        self.pscore_pf_compat = 0  # short array
        self.scs = vrna_sc_t  # vrna_sc_t array of arrays
        self.oldAliEn = 0  # int

        # Additional data fields for Distance Class Partitioning
        self.maxD1 = 0  # unsigned int
        self.maxD2 = 0  # unsigned int
        self.reference_pt1 = 0  # short array
        self.reference_pt2 = 0  # short array
        self.referenceBPs1 = 0  # unsigned int array
        self.referenceBPs2 = 0  # unsigned int array
        self.bpdist = 0  # unsigned int array
        self.mm1 = 0  # unsigned int array
        self.mm2 = 0  # unsigned int array

        # Additional data fields for local folding
        self.window_size = 0  # int
        self.ptype_local = 'ptype_local'  # char array of arrays
        self.zscore_data = vrna_zsc_dat_t  # vrna_zsc_dat_t









    


def vrna_idx_row_wise(length):
    idx = []
    for i in range(1, length + 1):
        idx.append((((length + 1 - i) * (length - i)) // 2) + length + 1)
    return idx

def vrna_bpp_symbol(x:list[float]):
    if x[0] > 0.667:
        return '.'

    if x[1] > 0.667:
        return '('

    if x[2] > 0.667:
        return ')'

    if (x[1] + x[2]) > x[0]:
        if (x[1] / (x[1] + x[2])) > 0.667:
            return '{'
        if (x[2] / (x[1] + x[2])) > 0.667:
            return '}'
        else:
            return '|'

    if x[0] > (x[1] + x[2]):
        return ','

    return ':'


def vrna_db_from_probs(p, length):
    s = [ 0 ] * (length + 1)
    if p is not None:
        index = vrna_idx_row_wise(length)
        for j in range(1, length + 1):
            _p = [1.0, 0.0, 0.0]
            for i in range(1, j):
                _p[2] += float(p[index[i] - j])  # j is paired downstream
                _p[0] -= float(p[index[i] - j])  # j is unpaired
            for i in range(j + 1, length + 1):
                _p[1] += float(p[index[j] - i])  # j is paired upstream
                _p[0] -= float(p[index[j] - i])  # j is unpaired
            s[j - 1] = vrna_bpp_symbol(_p)
        s[length] = '\0'
    return ''.join(s)
    


class HelperArrays:
    def __init__(self):
        self.prm_l = None
        self.prm_l1 = None
        self.prml = None
        self.ud_max_size = 0
        self.pmlu = []
        self.prm_MLbu = []


def get_ml_helper_arrays(fc:vrna_fold_compound_t):
    n = fc.length
    domains_up = fc.domains_up
    with_ud = 1 if (domains_up and domains_up.exp_energy_cb) else 0

    ml_helpers = HelperArrays()

    ml_helpers.prm_l = [0.0] * (n + 2)
    ml_helpers.prm_l1 = [0.0] * (n + 2)
    ml_helpers.prml = [0.0] * (n + 2)

    if with_ud:
        # Find out maximum size of any unstructured domain
        for u in range(domains_up.uniq_motif_count):
            if ml_helpers.ud_max_size < domains_up.uniq_motif_size[u]:
                ml_helpers.ud_max_size = domains_up.uniq_motif_size[u]

        ml_helpers.pmlu = [[0.0] * (n + 2) for _ in range(ml_helpers.ud_max_size + 1)]
        ml_helpers.prm_MLbu = [0.0] * (ml_helpers.ud_max_size + 1)

    return ml_helpers


class ConstraintsHelper:
    def __init__(self):
        self.hc_dat_ext = None
        self.hc_eval_ext = None

        self.hc_dat_hp = None
        self.hc_eval_hp = None

        self.hc_dat_int = None
        self.hc_eval_int = None

        self.hc_dat_mb = None
        self.hc_eval_mb = None

        self.sc_wrapper_ext = None
        self.sc_wrapper_hp = None
        self.sc_wrapper_int = None
        self.sc_wrapper_mb = None


def get_constraints_helper(fc):
    helpers = ConstraintsHelper()

    helpers.hc_eval_ext = prepare_hc_ext_def(fc, helpers.hc_dat_ext)
    helpers.hc_eval_hp = prepare_hc_hp_def(fc, helpers.hc_dat_hp)
    helpers.hc_eval_int = prepare_hc_int_def(fc, helpers.hc_dat_int)
    helpers.hc_eval_mb = prepare_hc_mb_def(fc, helpers.hc_dat_mb)

    init_sc_ext_exp(fc, helpers.sc_wrapper_ext)
    init_sc_hp_exp(fc, helpers.sc_wrapper_hp)
    init_sc_int_exp(fc, helpers.sc_wrapper_int)
    init_sc_mb_exp(fc, helpers.sc_wrapper_mb)

    return helpers


def compute_bpp_internal(fc,
                         l,
                         bp_correction,
                         corr_cnt,
                         corr_size,
                         Qmax,
                         ov,
                         constraints):
    import math
    from sys import float_info

    max_real = float_info.max

    hc_dat_local = constraints.hc_dat_int
    hc_eval = constraints.hc_eval_int
    sc_wrapper_int = constraints.sc_wrapper_int

    n = fc.length
    ptype = fc.ptype
    S1 = fc.sequence_encoding
    my_iindx = fc.iindx
    jindx = fc.jindx
    pf_params = fc.exp_params
    md = pf_params.model_details
    rtype = md.rtype
    hc = fc.hc
    sc = fc.sc
    domains_up = fc.domains_up
    with_ud = 1 if domains_up and domains_up.exp_energy_cb else 0
    hc_up_int = hc.up_int

    qb = fc.exp_matrices.qb
    probs = fc.exp_matrices.probs
    scale = fc.exp_matrices.scale

    for k in range(1, l):
        kl = my_iindx[k] - l

        if qb[kl] == 0.0:
            continue

        if hc.mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC:
            type_2 = rtype[vrna_get_ptype(jindx[l] + k, ptype)]

            for i in range(max(1, k - MAXLOOP - 1), k):
                u1 = k - i - 1
                if hc_up_int[i + 1] < u1:
                    continue

                max_j = l + 1 + MAXLOOP - u1

                if max_j > n:
                    max_j = n

                if max_j > l + 1 + hc_up_int[l + 1]:
                    max_j = l + 1 + hc_up_int[l + 1]

                u2 = 0

                for j in range(l + 1, max_j + 1):
                    ij = my_iindx[i] - j

                    if probs[ij] == 0.0:
                        continue

                    if hc_eval(i, j, k, l, hc_dat_local):
                        jij = jindx[j] + i
                        type = vrna_get_ptype(jij, ptype)
                        tmp2 = (probs[ij] *
                                exp_E_IntLoop(u1, u2, type, type_2,
                                              S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1], pf_params) *
                                scale[u1 + u2 + 2])

                        if sc_wrapper_int.pair:
                            tmp2 *= sc_wrapper_int.pair(i, j, k, l, sc_wrapper_int)

                        if with_ud:
                            qql = qqr = 0.0

                            if u1 > 0:
                                qql = domains_up.exp_energy_cb(fc, i + 1, k - 1,
                                                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP, domains_up.data)

                            if u2 > 0:
                                qqr = domains_up.exp_energy_cb(fc, l + 1, j - 1,
                                                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP, domains_up.data)

                            temp = tmp2
                            tmp2 += temp * qql
                            tmp2 += temp * qqr
                            tmp2 += temp * qql * qqr

                        if sc and sc.exp_f and sc.bt:
                            aux_bps = sc.bt(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc.data)
                            for ptr in aux_bps:
                                if ptr.i != 0:
                                    bp_correction[corr_cnt].i = ptr.i
                                    bp_correction[corr_cnt].j = ptr.j
                                    bp_correction[corr_cnt].p = tmp2 * qb[kl]
                                    corr_cnt += 1
                                    if corr_cnt == corr_size:
                                        corr_size += 5
                                        bp_correction.extend([None] * 5)

                        probs[kl] += tmp2

        if probs[kl] > Qmax[0]:
            Qmax[0] = probs[kl]
            if Qmax[0] > max_real / 10.0:
                print(f"P close to overflow: {k} {l} {probs[kl]} {qb[kl]}")

        if probs[kl] >= max_real:
            ov[0] += 1
            probs[kl] = float_info.max

    if md.gquad:
        compute_gquad_prob_internal(fc, l)



def compute_bpp_multibranch(fc,
                            l,
                            ml_helpers,
                            Qmax,
                            ov,
                            constraints):
    import math
    from sys import float_info

    max_real = float_info.max

    n = fc.length
    sn = fc.strand_number
    S = fc.sequence_encoding2
    S1 = fc.sequence_encoding
    my_iindx = fc.iindx
    jindx = fc.jindx
    pf_params = fc.exp_params
    md = pf_params.model_details
    rtype = md.rtype
    ptype = fc.ptype
    qb = fc.exp_matrices.qb
    qm = fc.exp_matrices.qm
    G = fc.exp_matrices.G
    probs = fc.exp_matrices.probs
    scale = fc.exp_matrices.scale
    expMLbase = fc.exp_matrices.expMLbase
    expMLclosing = pf_params.expMLclosing
    domains_up = fc.domains_up
    with_ud = 1 if domains_up and domains_up.exp_energy_cb else 0
    with_gquad = md.gquad
    expMLstem = exp_E_MLstem(0, -1, -1, pf_params) if with_gquad else 0

    hc_dat = constraints.hc_dat_mb
    hc_eval = constraints.hc_eval_mb
    sc_wrapper = constraints.sc_wrapper_mb

    prm_MLb = 0.0

    if sn[l + 1] != sn[l]:
        # set prm_l to 0 to get prm_l1 in the next round to be 0
        for i in range(n + 1):
            ml_helpers.prm_l[i] = 0.0
    else:
        for k in range(2, l):
            kl = my_iindx[k] - l
            i = k - 1
            prmt = prmt1 = 0.0

            ij = my_iindx[i] - (l + 2)
            lj = my_iindx[l + 1] - (l + 1)
            s3 = S1[i + 1]

            if sn[k] == sn[i]:
                for j in range(l + 2, n + 1):
                    if hc_eval(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, hc_dat):
                        tt = vrna_get_ptype_md(S[j], S[i], md)
                        ppp = (probs[ij] *
                               exp_E_MLstem(tt, S1[j - 1], s3, pf_params) *
                               qm[lj])

                        if sc_wrapper.pair:
                            ppp *= sc_wrapper.pair(i, j, sc_wrapper)

                        prmt += ppp

                ii = my_iindx[i]
                tt = vrna_get_ptype(jindx[l + 1] + i, ptype)
                tt = rtype[tt]

                if hc_eval(i, l + 1, i + 1, l, VRNA_DECOMP_PAIR_ML, hc_dat):
                    prmt1 = (probs[ii - (l + 1)] *
                             exp_E_MLstem(tt, S1[l], S1[i + 1], pf_params) *
                             expMLclosing)

                    if sc_wrapper.pair:
                        prmt1 *= sc_wrapper.pair(i, l + 1, sc_wrapper)

            prmt *= expMLclosing
            ml_helpers.prml[i] = prmt

            if hc_eval(k, l + 1, k, l, VRNA_DECOMP_ML_ML, hc_dat):
                ppp = (ml_helpers.prm_l1[i] *
                       expMLbase[1])

                if sc_wrapper.red_ml:
                    ppp *= sc_wrapper.red_ml(k, l + 1, k, l, sc_wrapper)

                if with_ud:
                    for cnt in range(domains_up.uniq_motif_count):
                        u = domains_up.uniq_motif_size[cnt]
                        if l + u < n and hc_eval(k, l + u, k, l, VRNA_DECOMP_ML_ML, hc_dat):
                            temp = (domains_up.exp_energy_cb(fc, l + 1, l + u,
                                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                             domains_up.data) *
                                    ml_helpers.pmlu[u][i] *
                                    expMLbase[u])

                            if sc_wrapper.red_ml:
                                temp *= sc_wrapper.red_ml(k, l + u, k, l, sc_wrapper)

                            ppp += temp

                    ml_helpers.pmlu[0][i] = ppp + prmt1

                ml_helpers.prm_l[i] = ppp + prmt1
            else:
                ml_helpers.prm_l[i] = prmt1

                if with_ud:
                    ml_helpers.pmlu[0][i] = prmt1

            if hc_eval(i, l, i + 1, l, VRNA_DECOMP_ML_ML, hc_dat):
                ppp = (prm_MLb *
                       expMLbase[1])

                if sc_wrapper.red_ml:
                    ppp *= sc_wrapper.red_ml(i, l, i + 1, l, sc_wrapper)

                if with_ud:
                    for cnt in range(domains_up.uniq_motif_count):
                        u = domains_up.uniq_motif_size[cnt]
                        if 1 + u <= i and hc_eval(i - u + 1, l, i + 1, l, VRNA_DECOMP_ML_ML, hc_dat):
                            temp = (domains_up.exp_energy_cb(fc, i - u + 1, i,
                                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                             domains_up.data) *
                                    ml_helpers.prm_MLbu[u] *
                                    expMLbase[u])

                            if sc_wrapper.red_ml:
                                temp *= sc_wrapper.red_ml(i - u + 1, l, i + 1, l, sc_wrapper)

                            ppp += temp

                    ml_helpers.prm_MLbu[0] = ppp + ml_helpers.prml[i]

                prm_MLb = ppp + ml_helpers.prml[i]

            else:
                prm_MLb = ml_helpers.prml[i]

                if with_ud:
                    ml_helpers.prm_MLbu[0] = ml_helpers.prml[i]

            ml_helpers.prml[i] += ml_helpers.prm_l[i]

            tt = ptype[jindx[l] + k]

            if with_gquad:
                if not tt and G[kl] == 0.0:
                    continue
            else:
                if qb[kl] == 0.0:
                    continue

            temp = prm_MLb

            if sn[k] == sn[k - 1]:
                if sc_wrapper.decomp_ml:
                    for i in range(1, k - 1):
                        temp += (ml_helpers.prml[i] *
                                 qm[my_iindx[i + 1] - (k - 1)] *
                                 sc_wrapper.decomp_ml(i + 1, l, k - 1, k, sc_wrapper))
                else:
                    for i in range(1, k - 1):
                        temp += ml_helpers.prml[i] * qm[my_iindx[i + 1] - (k - 1)]

            s5 = S1[k - 1] if k > 1 and sn[k] == sn[k - 1] else -1
            s3 = S1[l + 1] if l < n and sn[l + 1] == sn[l] else -1

            if with_gquad and qb[kl] == 0.0:
                temp *= G[kl] * expMLstem
            elif hc_eval(k, l, k, l, VRNA_DECOMP_ML_STEM, hc_dat):
                tt = vrna_get_ptype_md(S[j], S[i], md) if tt == 0 else tt
                temp *= exp_E_MLstem(tt, s5, s3, pf_params)

            if sc_wrapper.red_stem:
                temp *= sc_wrapper.red_stem(k, l, k, l, sc_wrapper)

            probs[kl] += temp * scale[2]

            if probs[kl] > Qmax[0]:
                Qmax[0] = probs[kl]
                if Qmax[0] > max_real / 10.0:
                    print(f"P close to overflow: {k} {l} {probs[kl]} {qb[kl]}")

            if probs[kl] >= max_real:
                ov[0] += 1
                probs[kl] = float_info.max

            # rotate prm_MLbu entries required for unstructured domain feature
            rotate_ml_helper_arrays_inner(ml_helpers)

    rotate_ml_helper_arrays_outer(ml_helpers)



def compute_bpp_internal_comparative(fc,
                                     l,
                                     bp_correction,
                                     corr_cnt,
                                     corr_size,
                                     Qmax,
                                     ov,
                                     constraints):
    import math
    from sys import float_info

    # Constants and initializations
    max_real = float_info.max
    kTn = fc.exp_params.kT / 10.0  # kT in cal/mol
    n = fc.length
    n_seq = fc.n_seq
    pscore = fc.pscore
    SS = fc.S
    S5 = fc.S5
    S3 = fc.S3
    a2s = fc.a2s
    my_iindx = fc.iindx
    jindx = fc.jindx
    pf_params = fc.exp_params
    md = pf_params.model_details
    hc = fc.hc
    hc_up_int = hc.up_int

    qb = fc.exp_matrices.qb
    probs = fc.exp_matrices.probs
    scale = fc.exp_matrices.scale

    hc_eval = constraints.hc_eval_int
    hc_dat_local = constraints.hc_dat_int
    sc_wrapper_int = constraints.sc_wrapper_int

    # Allocate memory for type array
    tt = [0] * n_seq

    # Loop over k for internal loops
    for k in range(1, l):
        kl = my_iindx[k] - l

        if qb[kl] == 0.0:
            continue

        # Check if (k,l) forms a substem of an internal loop
        if hc.mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC:
            psc_exp = math.exp(pscore[jindx[l] + k] / kTn)

            # Compute type array for current k,l pair
            for s in range(n_seq):
                tt[s] = vrna_get_ptype_md(SS[s][l], SS[s][k], md)

            # Loop over i and j to compute probabilities
            for i in range(max(1, k - md.max_nesting_pair - 1), k):
                u1 = k - i - 1
                if hc_up_int[i + 1] < u1:
                    continue

                for j in range(l + 1, min(l + md.max_nesting_pair - k + i + 2, n + 1)):
                    ij = my_iindx[i] - j

                    if probs[ij] == 0.0:
                        continue

                    u2 = j - l - 1

                    if hc_up_int[l + 1] < u2:
                        break

                    if hc_eval(i, j, k, l, hc_dat_local):
                        tmp2 = (probs[ij] *
                                scale[u1 + u2 + 2] *
                                psc_exp)

                        for s in range(n_seq):
                            u1_loc = a2s[s][k - 1] - a2s[s][i]
                            u2_loc = a2s[s][j - 1] - a2s[s][l]
                            type = vrna_get_ptype_md(SS[s][i], SS[s][j], md)
                            tmp2 *= exp_E_IntLoop(u1_loc,
                                                  u2_loc,
                                                  type,
                                                  tt[s],
                                                  S3[s][i],
                                                  S5[s][j],
                                                  S5[s][k],
                                                  S3[s][l],
                                                  pf_params)

                        if sc_wrapper_int.pair:
                            tmp2 *= sc_wrapper_int.pair(i, j, k, l, sc_wrapper_int)

                        probs[kl] += tmp2

        # Check and handle overflow conditions
        if probs[kl] > Qmax[0]:
            Qmax[0] = probs[kl]
            if Qmax[0] > max_real / 10.0:
                print(f"P close to overflow: {k} {l} {probs[kl]} {qb[kl]}")

        if probs[kl] >= max_real:
            ov[0] += 1
            probs[kl] = float_info.max

    # Free allocated memory
    del tt

    # Call function to compute probabilities for G-quadruplexes if enabled
    if md.gquad:
        compute_gquad_prob_internal_comparative(fc, l)



def compute_bpp_multibranch_comparative(fc,
                                        l,
                                        ml_helpers,
                                        Qmax,
                                        ov,
                                        constraints):
    import math
    from sys import float_info

    # Constants and initializations
    max_real = float_info.max
    kTn = fc.exp_params.kT / 10.0  # kT in cal/mol
    n = fc.length
    n_seq = fc.n_seq
    S = fc.S
    S5 = fc.S5
    S3 = fc.S3
    a2s = fc.a2s
    sn = fc.strand_number
    pscore = fc.pscore
    my_iindx = fc.iindx
    jindx = fc.jindx
    pf_params = fc.exp_params
    md = pf_params.model_details
    qb = fc.exp_matrices.qb
    qm = fc.exp_matrices.qm
    G = fc.exp_matrices.G
    probs = fc.exp_matrices.probs
    scale = fc.exp_matrices.scale
    expMLbase = fc.exp_matrices.expMLbase
    expMLclosing = pf_params.expMLclosing
    with_gquad = md.gquad
    hc = fc.hc
    scs = fc.scs
    expMLstem = pow(exp_E_MLstem(0, -1, -1, pf_params), n_seq) if with_gquad else 0.0

    prm_MLb = 0.0
    temp = 0.0

    if sn[l + 1] != sn[l]:
        for i in range(n + 1):
            ml_helpers.prm_l[i] = 0
    else:
        for k in range(2, l):
            kl = my_iindx[k] - l
            i = k - 1
            prmt = 0.0
            prmt1 = 0.0

            ij = my_iindx[i] - (l + 2)
            lj = my_iindx[l + 1] - (l + 1)

            if sn[k] == sn[i]:
                for j in range(l + 2, n + 1):
                    if (hc.mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) and (sn[j] == sn[j - 1]):
                        ppp = probs[ij] * qm[lj]

                        for s in range(n_seq):
                            tt = vrna_get_ptype_md(S[s][j], S[s][i], md)
                            ppp *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params)

                        if scs:
                            for s in range(n_seq):
                                if scs[s]:
                                    if scs[s].exp_energy_bp:
                                        ppp *= scs[s].exp_energy_bp[jindx[j] + i]

                        prmt += ppp

                ii = my_iindx[i]

                if hc.mx[(l + 1) * n + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP:
                    prmt1 = probs[ii - (l + 1)] * pow(expMLclosing, n_seq)

                    for s in range(n_seq):
                        tt = vrna_get_ptype_md(S[s][l + 1], S[s][i], md)
                        prmt1 *= exp_E_MLstem(tt, S5[s][l + 1], S3[s][i], pf_params)

                    if scs:
                        for s in range(n_seq):
                            if scs[s]:
                                if scs[s].exp_energy_bp:
                                    prmt1 *= scs[s].exp_energy_bp[jindx[l + 1] + i]

                prmt *= pow(expMLclosing, n_seq)

                ml_helpers.prml[i] = prmt

                if hc.up_ml[l + 1]:
                    ppp = ml_helpers.prm_l1[i] * expMLbase[1]

                    if scs:
                        for s in range(n_seq):
                            if scs[s]:
                                if scs[s].exp_energy_up:
                                    ppp *= scs[s].exp_energy_up[a2s[s][l + 1]][1]

                    ml_helpers.prm_l[i] = ppp + prmt1
                else:
                    ml_helpers.prm_l[i] = prmt1

                if hc.up_ml[i]:
                    ppp = prm_MLb * expMLbase[1]

                    if scs:
                        for s in range(n_seq):
                            if scs[s]:
                                if scs[s].exp_energy_up:
                                    ppp *= scs[s].exp_energy_up[a2s[s][i]][1]

                    prm_MLb = ppp + ml_helpers.prml[i]
                else:
                    prm_MLb = ml_helpers.prml[i]

                ml_helpers.prml[i] += ml_helpers.prm_l[i]

                if with_gquad:
                    if qb[kl] == 0 and G[kl] == 0:
                        continue

                temp = prm_MLb

                if sn[k] == sn[k - 1]:
                    for i in range(1, k - 1):
                        if sn[i + 1] == sn[i]:
                            temp += ml_helpers.prml[i] * qm[my_iindx[i + 1] - (k - 1)]

                if with_gquad and qb[kl] == 0:
                    temp *= G[kl] * expMLstem
                elif hc.mx[l * n + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC:
                    for s in range(n_seq):
                        tt = vrna_get_ptype_md(S[s][k], S[s][l], md)
                        temp *= exp_E_MLstem(tt, S5[s][k], S3[s][l], pf_params)

                probs[kl] += temp * scale[2] * math.exp(pscore[jindx[l] + k] / kTn)

                if probs[kl] > Qmax[0]:
                    Qmax[0] = probs[kl]
                    if Qmax[0] > max_real / 10.0:
                        print(f"P close to overflow: {k} {l} {probs[kl]} {qb[kl]}")

                if probs[kl] >= max_real:
                    ov[0] += 1
                    probs[kl] = float_info.max

                # Rotate ml_helpers arrays
                rotate_ml_helper_arrays_inner(ml_helpers)

    # Rotate ml_helpers arrays
    rotate_ml_helper_arrays_outer(ml_helpers)






def compute_bpp_external(fc, constraints):
    n = fc.length
    my_iindx = fc.iindx
    matrices = fc.exp_matrices
    qb = matrices.qb
    probs = matrices.probs
    q1k = matrices.q1k
    qln = matrices.qln
    hc_dat = constraints.hc_dat_ext
    evaluate = constraints.hc_eval_ext

    if fc.type == VRNA_FC_TYPE_SINGLE:
        contrib_f = contrib_ext_pair
    else:
        contrib_f = contrib_ext_pair_comparative

    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            ij = my_iindx[i] - j
            probs[ij] = 0.0

            if evaluate(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, hc_dat) and qb[ij] > 0.0:
                probs[ij] = q1k[i - 1] * qln[j + 1] / q1k[n]
                probs[ij] *= contrib_f(fc, i, j, constraints)


def multistrand_update_Y5(fc, l, Y5, Y5p, constraints):
    n = fc.length
    sn = fc.strand_number
    se = fc.strand_end
    my_iindx = fc.iindx
    q = fc.exp_matrices.q
    probs = fc.exp_matrices.probs
    scale = fc.exp_matrices.scale
    pf_params = fc.exp_params
    md = pf_params.model_details
    S = fc.sequence_encoding2
    S1 = fc.sequence_encoding
    sc_wrapper = constraints.sc_wrapper_ext
    sc_red_stem = sc_wrapper.red_stem
    sc_split = sc_wrapper.split

    # Compute Y5 for all strands
    for s in range(fc.strands):
        Y5[s] = 0.0

        if se[s] < l and sn[l] == sn[l + 1]:
            # Pre-compute newly available Y5p[s][j] with j == l + 1
            end = se[s]
            j = l + 1
            Y5p[s][j] = 0.0

            if probs[my_iindx[end] - j] > 0:
                type = vrna_get_ptype_md(S[j], S[end], md)
                qtmp = probs[my_iindx[end] - j] * \
                       vrna_exp_E_ext_stem(type, S1[j - 1], -1, pf_params) * \
                       scale[2]

                if sc_red_stem:
                    qtmp *= sc_red_stem(j, end, j, end, sc_wrapper)

                Y5p[s][j] += qtmp

            for i in range(1, end):
                if probs[my_iindx[i] - j] > 0 and sn[i] == sn[i + 1]:
                    type = vrna_get_ptype_md(S[j], S[i], md)
                    qtmp = probs[my_iindx[i] - j] * \
                           vrna_exp_E_ext_stem(type, S1[j - 1], S1[i + 1], pf_params) * \
                           q[my_iindx[i + 1] - end] * \
                           scale[2]

                    if sc_red_stem:
                        qtmp *= sc_red_stem(j, i, j, i, sc_wrapper)
                    if sc_split:
                        qtmp *= sc_split(i, end, i + 1, sc_wrapper)

                    Y5p[s][j] += qtmp

            if probs[my_iindx[i] - j] > 0 and sn[i] == sn[i + 1]:
                type = vrna_get_ptype_md(S[j], S[i], md)
                qtmp = probs[my_iindx[i] - j] * \
                       vrna_exp_E_ext_stem(type, S1[j - 1], S1[i + 1], pf_params) * \
                       scale[2]

                if sc_red_stem:
                    qtmp *= sc_red_stem(j, i, j, i, sc_wrapper)

                Y5p[s][j] += qtmp

            # Recompute Y5[s]
            Y5[s] += Y5p[s][l + 1]
            for j in range(l + 2, n + 1):
                qtmp = q[my_iindx[l + 1] - (j - 1)] * Y5p[s][j]

                if sc_split:
                    qtmp *= sc_split(l + 1, j, j, sc_wrapper)

                Y5[s] += qtmp


def multistrand_update_Y3(fc, l, Y3, Y3p, constraints):
    n = fc.length
    sn = fc.strand_number
    ss = fc.strand_start
    S = fc.sequence_encoding2
    S1 = fc.sequence_encoding
    my_iindx = fc.iindx
    q = fc.exp_matrices.q
    probs = fc.exp_matrices.probs
    scale = fc.exp_matrices.scale
    pf_params = fc.exp_params
    md = pf_params.model_details
    sc_wrapper = constraints.sc_wrapper_ext
    sc_red_stem = sc_wrapper.red_stem
    sc_split = sc_wrapper.split

    for s in range(fc.strands):
        start = ss[s]
        if start == l + 1:
            for i in range(1, start):
                Y3p[s][i] = 0.0

                if sn[i] == sn[i + 1]:
                    if probs[my_iindx[i] - start] > 0:
                        type = vrna_get_ptype_md(S[start], S[i], md)
                        qtmp = probs[my_iindx[i] - start] * \
                               vrna_exp_E_ext_stem(type, -1, S1[i + 1], pf_params) * \
                               scale[2]

                        if sc_red_stem:
                            qtmp *= sc_red_stem(start, i, start, i, sc_wrapper)

                        Y3p[s][i] += qtmp

                    for j in range(start + 1, n + 1):
                        if probs[my_iindx[i] - j] > 0 and sn[j - 1] == sn[j]:
                            type = vrna_get_ptype_md(S[j], S[i], md)
                            qtmp = probs[my_iindx[i] - j] * \
                                   vrna_exp_E_ext_stem(type, S1[j - 1], S1[i + 1], pf_params) * \
                                   q[my_iindx[start] - (j - 1)] * \
                                   scale[2]

                            if sc_red_stem:
                                qtmp *= sc_red_stem(j, i, j, i, sc_wrapper)
                            if sc_split:
                                qtmp *= sc_split(start, j, j, sc_wrapper)

                            Y3p[s][i] += qtmp

            for k in range(1, start):
                Y3[s][k] = 0.0

                if sn[k - 1] == sn[k]:
                    for i in range(1, k):
                        if sn[i] == sn[i + 1]:
                            qtmp = q[my_iindx[i + 1] - (k - 1)] * Y3p[s][i]

                            if sc_split:
                                qtmp *= sc_split(i, k - 1, i + 1, sc_wrapper)

                            Y3[s][k] += qtmp

                    Y3[s][k] += Y3p[s][k - 1]


def multistrand_contrib(fc, l, Y5, Y3, constraints, Qmax, ov):
    sn = fc.strand_number
    ss = fc.strand_start
    se = fc.strand_end
    S = fc.sequence_encoding2
    S1 = fc.sequence_encoding
    pf_params = fc.exp_params
    md = pf_params.model_details
    my_iindx = fc.iindx
    q = fc.exp_matrices.q
    qb = fc.exp_matrices.qb
    probs = fc.exp_matrices.probs
    sc_wrapper = constraints.sc_wrapper_ext
    sc_red_stem = sc_wrapper.red_stem

    for k in range(l - 1, 1, -1):
        kl = my_iindx[k] - l
        if qb[kl] > 0:
            tmp = 0.0
            for s in range(fc.strands):
                end = se[s]
                start = ss[s]
                if end == k - 1:
                    tmp += Y5[s]
                elif end < k - 1 and sn[k - 1] == sn[k]:
                    tmp += Y5[s] * q[my_iindx[end + 1] - (k - 1)]
                elif start == l + 1:
                    tmp += Y3[s][k]
                elif start > l + 1 and sn[l] == sn[l + 1]:
                    tmp += Y3[s][k] * q[my_iindx[l + 1] - (start - 1)]

            type = vrna_get_ptype_md(S[k], S[l], md)
            s5 = S1[k - 1] if sn[k - 1] == sn[k] else -1
            s3 = S1[l + 1] if sn[l] == sn[l + 1] else -1
            qtmp = vrna_exp_E_ext_stem(type, s5, s3, pf_params)

            if sc_red_stem:
                qtmp *= sc_red_stem(k, l, k, l, sc_wrapper)

            probs[kl] += tmp * qtmp


def ud_outside_ext_loops(vc):
    n = vc.length
    q1k = vc.exp_matrices.q1k
    qln = vc.exp_matrices.qln
    scale = vc.exp_matrices.scale
    hc_up = vc.hc.up_ext
    domains_up = vc.domains_up
    sc = vc.sc

    for i in range(1, n + 1):
        motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP)

        # 1. Exterior loops
        if motif_list:
            cnt = 0
            while motif_list[cnt] != -1:
                u = motif_list[cnt]
                j = i + u - 1
                if j <= n:
                    if hc_up[i] >= u:
                        temp = q1k[i - 1] * qln[j + 1] / q1k[n]
                        temp *= domains_up.exp_energy_cb(vc,
                                                         i,
                                                         j,
                                                         VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                         domains_up.data)

                        if sc and sc.exp_energy_up:
                            temp *= sc.exp_energy_up[i][u]

                        temp *= scale[u]

                        if temp > 0:
                            domains_up.probs_add(vc,
                                                 i,
                                                 j,
                                                 VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                 temp,
                                                 domains_up.data)
                cnt += 1

            
def ud_outside_hp_loops(vc):
    n = vc.length
    my_iindx = vc.iindx
    probs = vc.exp_matrices.probs
    hc_up = vc.hc.up_hp
    domains_up = vc.domains_up

    for i in range(1, n + 1):
        motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP)

        # 2. Hairpin loops
        if motif_list:
            cnt = 0
            while motif_list[cnt] != -1:
                u = motif_list[cnt]
                outside = 0.
                j = i + u - 1
                if j < n:
                    if hc_up[i] >= u:
                        exp_motif_en = domains_up.exp_energy_cb(vc,
                                                                 i,
                                                                 j,
                                                                 VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                                 domains_up.data)

                        # Compute the contribution of all hairpins with bound motif
                        for k in range(1, i):
                            for l in range(j + 1, n + 1):
                                kl = my_iindx[k] - l
                                if probs[kl] > 0:
                                    ud_bak = vc.domains_up
                                    vc.domains_up = None
                                    temp = vrna_exp_E_hp_loop(vc, k, l)
                                    vc.domains_up = ud_bak

                                    # Add contribution of motif
                                    if temp > 0:
                                        temp *= exp_motif_en * probs[kl]

                                        q1 = q2 = 0.
                                        # Add contributions of other motifs in remaining unpaired segments
                                        if (i - k - 1) > 0:
                                            q1 = domains_up.exp_energy_cb(vc,
                                                                          k + 1, i - 1,
                                                                          VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                                                          domains_up.data)

                                        if (l - j - 1) > 0:
                                            q2 = domains_up.exp_energy_cb(vc,
                                                                          j + 1, l - 1,
                                                                          VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                                                          domains_up.data)

                                        outside += temp
                                        outside += temp * q1
                                        outside += temp * q1 * q2
                                        outside += temp * q2

                    if outside > 0:
                        domains_up.probs_add(vc,
                                             i, j,
                                             VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                             outside,
                                             domains_up.data)

                cnt += 1


def ud_outside_int_loops(vc):
    n = vc.length
    my_iindx = vc.iindx
    qb = vc.exp_matrices.qb
    probs = vc.exp_matrices.probs
    hc_up = vc.hc.up_int
    domains_up = vc.domains_up

    for i in range(2, n + 1):
        motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP)

        # 3. Interior loops
        if motif_list:
            cnt = 0
            while motif_list[cnt] != -1:
                u = motif_list[cnt]
                outside = 0.
                j = i + u - 1

                if j < n:
                    if hc_up[i] >= u:
                        exp_motif_en = domains_up.exp_energy_cb(vc,
                                                                 i,
                                                                 j,
                                                                 VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                                 domains_up.data)

                        # 3.1 motif is within 5' loop
                        kmin = j - MAXLOOP - 1
                        kmin = max(kmin, 1)
                        for k in range(kmin, i):
                            pmax = k + MAXLOOP + 1
                            pmax = min(pmax, n)
                            for p in range(j + 1, n):
                                for q in range(p + 1, n):
                                    pq = my_iindx[p] - q
                                    if qb[pq] == 0:
                                        continue

                                    lmax = k + MAXLOOP + q - p + 2
                                    lmax = min(lmax, n)
                                    for l in range(q + 1, lmax + 1):
                                        kl = my_iindx[k] - l
                                        if probs[kl] > 0:
                                            ud_bak = vc.domains_up
                                            vc.domains_up = None
                                            temp = vrna_exp_E_interior_loop(vc, k, l, p, q)
                                            vc.domains_up = ud_bak

                                            if temp > 0:
                                                temp *= probs[kl] * qb[pq] * exp_motif_en

                                                q1 = q2 = q3 = 0.
                                                if (l - q - 1) > 0:
                                                    q1 = domains_up.exp_energy_cb(vc,
                                                                                   q + 1, l - 1,
                                                                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                                                   domains_up.data)

                                                if (i - k - 1) > 0:
                                                    q2 = domains_up.exp_energy_cb(vc,
                                                                                   k + 1, i - 1,
                                                                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                                                   domains_up.data)

                                                if (p - j - 1) > 0:
                                                    q3 = domains_up.exp_energy_cb(vc,
                                                                                   j + 1, p - 1,
                                                                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                                                   domains_up.data)

                                                outside += temp
                                                outside += temp * q1
                                                outside += temp * q1 * q2
                                                outside += temp * q1 * q2 * q3
                                                outside += temp * q2
                                                outside += temp * q2 * q3
                                                outside += temp * q3

                        # 3.2 motif is within 3' loop
                        for k in range(1, i - 1):
                            pmax = k + i + MAXLOOP - j
                            pmax = min(pmax, n)
                            for p in range(k + 1, pmax + 1):
                                qmin = p + j - k - MAXLOOP - 1
                                qmin = max(qmin, p + 1)
                                for q in range(i - 1, qmin - 1, -1):
                                    pq = my_iindx[p] - q
                                    if qb[pq] == 0:
                                        continue

                                    lmax = k + q - p + MAXLOOP + 2
                                    lmax = min(lmax, n)
                                    for l in range(j + 1, lmax):
                                        kl = my_iindx[k] - l
                                        if probs[kl] > 0:
                                            ud_bak = vc.domains_up
                                            vc.domains_up = None
                                            temp = vrna_exp_E_interior_loop(vc, k, l, p, q)
                                            vc.domains_up = ud_bak

                                            if temp > 0:
                                                q1 = q2 = q3 = 0.
                                                temp *= probs[kl] * qb[pq] * exp_motif_en

                                                if (l - j - 1) > 0:
                                                    q1 = domains_up.exp_energy_cb(vc,
                                                                                   j + 1, l - 1,
                                                                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                                                   domains_up.data)

                                                if (i - q - 1) > 0:
                                                    q2 = domains_up.exp_energy_cb(vc,
                                                                                   q + 1, i - 1,
                                                                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                                                   domains_up.data)

                                                if (p - k - 1) > 0:
                                                    q3 = domains_up.exp_energy_cb(vc,
                                                                                   k + 1, p - 1,
                                                                                   VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                                                                   domains_up.data)

                                                outside += temp
                                                outside += temp * q1
                                                outside += temp * q1 * q2
                                                outside += temp * q1 * q2 * q3
                                                outside += temp * q2
                                                outside += temp * q2 * q3
                                                outside += temp * q3

                    if outside > 0:
                        domains_up.probs_add(vc,
                                             i, j,
                                             VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                             outside,
                                             domains_up.data)

                cnt += 1


def ud_outside_mb_loops(vc):
    n = vc.length
    S = vc.sequence_encoding
    my_iindx = vc.iindx
    jindx = vc.jindx
    ptype = vc.ptype
    pf_params = vc.exp_params
    md = pf_params.model_details
    qb = vc.exp_matrices.qb
    qm = vc.exp_matrices.qm
    probs = vc.exp_matrices.probs
    scale = vc.exp_matrices.scale
    hc_up = vc.hc.up_ml
    hc = vc.hc.mx
    domains_up = vc.domains_up
    sc = vc.sc
    rtype = md.rtype
    expMLbase = vc.exp_matrices.expMLbase
    expMLclosing = pf_params.expMLclosing

    ud_max_size = max(domains_up.uniq_motif_size) if domains_up.uniq_motif_count > 0 else 0

    for i in range(1, n + 1):
        motif_list = vrna_ud_get_motif_size_at(vc, i, VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP)

        # 4. Multibranch loops
        if motif_list:
            cnt = 0
            while motif_list[cnt] != -1:
                u = motif_list[cnt]
                outside = 0.
                j = i + u - 1

                if j < n:
                    if hc_up[i] >= u:
                        exp_motif_en = domains_up.exp_energy_cb(vc,
                                                                 i,
                                                                 j,
                                                                 VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                                 domains_up.data)

                        exp_motif_en *= expMLbase[u]

                        if sc and sc.exp_energy_up:
                            exp_motif_en *= sc.exp_energy_up[i][u]

                        temp = 0.

                        # 4.1 Motif [i:j] is somewhere in between branching stems
                        for l in range(j + 1, n + 1):
                            for k in range(i - 1, 0, -1):
                                kl = my_iindx[k] - l
                                if probs[kl] > 0:
                                    if hc[n * l + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP:
                                        jkl = jindx[l] + k
                                        tt = rtype[vrna_get_ptype(jkl, ptype)]
                                        qqq = (probs[kl] *
                                               qm[my_iindx[k + 1] - (i - 1)] *
                                               qm[my_iindx[j + 1] - (l - 1)] *
                                               exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params) *
                                               expMLclosing *
                                               scale[2])

                                        if sc and sc.exp_energy_bp:
                                            qqq *= sc.exp_energy_bp[jkl]

                                        temp += qqq

                        outside += temp * exp_motif_en

                        # 4.2 Motif is in left-most unpaired stretch of multiloop
                        qm1ui = [None] * (ud_max_size + 1)
                        for l in range(ud_max_size + 1):
                            qm1ui[l] = [0.0] * (n + 2)

                        exp_motif_ml_left = 0.0
                        for l in range(j + 1, n + 1):
                            lqq = 0.0
                            rqq = 0.0
                            for k in range(i - 1, 0, -1):
                                up = i - k - 1
                                kl = my_iindx[k] - l
                                if (hc[n * l + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP and
                                        probs[kl] > 0 and
                                        hc_up[k + 1] >= up):
                                    jkl = jindx[l] + k
                                    tt = rtype[vrna_get_ptype(jkl, ptype)]
                                    temp = (probs[kl] *
                                            expMLbase[up] *
                                            exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params) *
                                            expMLclosing *
                                            scale[2])

                                    if sc:
                                        if sc.exp_energy_up:
                                            temp *= sc.exp_energy_up[k + 1][up]

                                    lqq += temp
                                    lqq += temp * domains_up.exp_energy_cb(vc,
                                                                            k + 1, i - 1,
                                                                            VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                                            domains_up.data)

                            for u in range(j + 1, l):
                                if hc_up[l - 1]:
                                    temp = qm1ui[1][u] * expMLbase[1]
                                    if sc and sc.exp_energy_up:
                                        temp *= sc.exp_energy_up[l - 1][1]

                                    qm1ui[0][u] = temp
                                else:
                                    qm1ui[0][u] = 0.0

                                for cnt in range(domains_up.uniq_motif_count):
                                    size = domains_up.uniq_motif_size[cnt]
                                    if (u < l - size and
                                            hc_up[l - size] >= size):
                                        temp = (qm1ui[size][u] *
                                                expMLbase[size] *
                                                domains_up.exp_energy_cb(vc,
                                                                        l - size,
                                                                        l - 1,
                                                                        VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                                        domains_up.data))

                                        if sc and sc.exp_energy_up:
                                            temp *= sc.exp_energy_up[l - size][size]

                                        qm1ui[0][u] += temp

                                if hc[n * (l - 1) + u] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC:
                                    tt = vrna_get_ptype(jindx[l - 1] + u, ptype)
                                    temp = (qb[my_iindx[u] - (l - 1)] *
                                            exp_E_MLstem(tt, S[u - 1], S[l], pf_params))

                                    qm1ui[0][u] += temp

                                rqq += (qm[my_iindx[j + 1] - (u - 1)] *
                                        qm1ui[0][u])

                            exp_motif_ml_left += (lqq * rqq)

                            tmp = qm1ui[ud_max_size]
                            for cnt in range(ud_max_size, 0, -1):
                                qm1ui[cnt] = qm1ui[cnt - 1]
                            qm1ui[0] = tmp

                        outside += (exp_motif_ml_left * exp_motif_en)

                        # Cleanup memory
                        for l in range(ud_max_size + 1):
                            del qm1ui[l]
                        del qm1ui

                        # 4.3 Motif is in right-most unpaired stretch of multiloop
                        qmli = [0.0] * n
                        exp_motif_ml_right = 0.0
                        for k in range(i - 1, 0, -1):
                            lqq = 0.0
                            rqq = 0.0

                            for u in range(k + 1, i):
                                if (hc[n * u + k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC and
                                        sc and sc.exp_energy_up and
                                        sc.exp_energy_up[u + 1][i - u - 1] >= u + 1):
                                    temp = (qb[my_iindx[k] - u] *
                                            expMLbase[i - 1 - (u + 1) + 1])

                                    if sc and sc.exp_energy_up:
                                        temp *= sc.exp_energy_up[u + 1][i - u - 1]

                                    qmli[k] += temp

                                    qmli[k] += (temp *
                                                domains_up.exp_energy_cb(vc,
                                                                        u + 1, i - 1,
                                                                        VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                                        domains_up.data))

                            for u in range(k, i):
                                lqq += (qm[my_iindx[k + 1] - (u - 1)] *
                                        qmli[u])

                            for l in range(j + 1, n + 1):
                                kl = my_iindx[k] - l
                                if (hc[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP and
                                        sc and
                                        hc_up[j + 1] >= l - j - 1):
                                    jkl = jindx[l] + k
                                    tt = rtype[vrna_get_ptype(jkl, ptype)]
                                    temp = (probs[kl] *
                                            exp_E_MLstem(tt, S[l - 1], S[k + 1], pf_params) *
                                            expMLclosing *
                                            scale[2] *
                                            expMLbase[l - j - 1])

                                    if sc and sc.exp_energy_bp:
                                        temp *= sc.exp_energy_bp[jkl]

                                    rqq += temp

                                    rqq += (temp *
                                            domains_up.exp_energy_cb(vc,
                                                                    j + 1, l - 1,
                                                                    VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                                                    domains_up.data))

                            exp_motif_ml_right += (rqq * lqq)

                        outside += (exp_motif_ml_right * exp_motif_en)

                    if outside > 0:
                        domains_up.probs_add(vc,
                                             i, j,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                             outside,
                                             domains_up.data)

                cnt += 1


def vrna_exp_E_hp_loop(fc, i, j):
    evaluate = None
    hc_dat_local = None

    if fc.hc.type == VRNA_HC_WINDOW:
        evaluate, hc_dat_local = prepare_hc_hp_def_window(fc, hc_dat_local)
    else:
        evaluate, hc_dat_local = prepare_hc_hp_def(fc, hc_dat_local)

    if i > 0 and j > 0:
        if evaluate(i, j, i, j, VRNA_DECOMP_PAIR_HP, hc_dat_local):
            if j > i:  # linear case
                return exp_eval_hp_loop(fc, i, j)
            else:  # circular case
                return exp_eval_ext_hp_loop(fc, j, i)

    return 0.0



            






def pf_create_bppm(vc, structure):
    ov = 0
    n = vc.length
    pscore = vc.pscore if vc.type == VRNA_FC_TYPE_COMPARATIVE else None
    pf_params = vc.exp_params
    md = pf_params.model_details
    my_iindx = vc.iindx
    jindx = vc.jindx
    circular = md.circ
    with_gquad = md.gquad

    hc = vc.hc
    sc = vc.sc

    domains_up = vc.domains_up
    matrices = vc.exp_matrices

    qb = matrices.qb
    G = matrices.G
    probs = matrices.probs
    q1k = matrices.q1k
    qln = matrices.qln

    with_ud = 1 if (domains_up and domains_up.exp_energy_cb) else 0
    with_ud_outside = 1 if (with_ud and domains_up.probs_add) else 0

    if qb is not None and probs is not None and (circular or (q1k is not None and qln is not None)):
        Qmax = 0
        corr_size = 5
        corr_cnt = 0
        bp_correction = np.zeros(corr_size, dtype=[('i', int), ('j', int), ('p', float)])
        ml_helpers  = get_ml_helper_arrays(vc)
        constraints = get_constraints_helper(vc)
        
        if vc.type == VRNA_FC_TYPE_SINGLE:
            compute_bpp_int = compute_bpp_internal
            compute_bpp_mul = compute_bpp_multibranch
        else:
            compute_bpp_int = compute_bpp_internal_comparative
            compute_bpp_mul = compute_bpp_multibranch_comparative
            
        
        if vc.strands > 1:
            Y5 = np.zeros(vc.strands)
            Y5p = [np.zeros(n + 1) for _ in range(vc.strands)]
            Y3 = [np.zeros(n + 1) for _ in range(vc.strands)]
            Y3p = [np.zeros(n + 1) for _ in range(vc.strands)]

        probs[my_iindx[1] - 1:my_iindx[n] - n + 1] = 0

        if circular:
            bppm_circ(vc, constraints)
        else:
            compute_bpp_external(vc, constraints)

        l = n
        compute_bpp_int(vc, l, bp_correction, corr_cnt, corr_size, Qmax, ov, constraints=None)

        for l in range(n - 1, 1, -1):
            compute_bpp_int(vc, l, bp_correction, corr_cnt, corr_size, Qmax, ov, constraints=None)
            compute_bpp_mul(vc, l, ml_helpers=None, Qmax=Qmax, ov=ov, constraints=None)

            if vc.strands > 1:
                multistrand_update_Y5(vc, l, Y5, Y5p, constraints=None)
                multistrand_update_Y3(vc, l, Y3, Y3p, constraints=None)
                multistrand_contrib(vc, l, Y5, Y3, constraints=None, Qmax=Qmax, ov=ov)

        if vc.type == VRNA_FC_TYPE_SINGLE:
            if with_ud_outside:
                ud_outside_ext_loops(vc)
                ud_outside_hp_loops(vc)
                ud_outside_int_loops(vc)
                ud_outside_mb_loops(vc)

            if sc and sc.f and sc.bt:
                for i in range(1, n + 1):
                    for j in range(i + 1, n + 1):
                        ij = my_iindx[i] - j
                        if hc.mx[i * n + j] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP:
                            aux_bps = sc.bt(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc.data)
                            if aux_bps:
                                qhp = vrna_exp_E_hp_loop(vc, i, j)
                                for ptr in aux_bps:
                                    if ptr.i != 0:
                                        if corr_cnt == corr_size:
                                            corr_size += 5
                                            bp_correction = np.resize(bp_correction, corr_size)
                                        bp_correction[corr_cnt] = (ptr.i, ptr.j, probs[ij] * qhp)
                                        corr_cnt += 1

                for i in range(corr_cnt):
                    ij = my_iindx[bp_correction[i]['i']] - bp_correction[i]['j']
                    probs[ij] += bp_correction[i]['p'] / qb[ij]

        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                ij = my_iindx[i] - j
                if with_gquad:
                    if qb[ij] > 0.:
                        probs[ij] *= qb[ij]
                        if vc.type == VRNA_FC_TYPE_COMPARATIVE:
                            probs[ij] *= np.exp(-pscore[jindx[j] + i] / (pf_params.kT / 10))
                    elif G[ij] > 0.:
                        probs[ij] += q1k[i - 1] * G[ij] * qln[j + 1] / q1k[n]
                else:
                    if qb[ij] > 0.:
                        probs[ij] *= qb[ij]
                        if vc.type == VRNA_FC_TYPE_COMPARATIVE:
                            probs[ij] *= np.exp(-pscore[jindx[j] + i] / (pf_params.kT / 10))

        if structure is not None:
            s = vrna_db_from_probs(probs, n)
            structure[:] = s[:]
            structure[n] = '\0'

        if ov > 0:
            print(f"{ov} overflows occurred while backtracking; you might try a smaller pf_scale than {pf_params.pf_scale}")

    else:
        print("bppm calculations have to be done after calling forward recursion")
        return 0

    return 1


def generate_struc(probs, n):
    struc = "(((((((((((((((((((((()))))))))))...)))))))))))"
    seq = "aactgcccactcagtacatcaaTTGATGTACTGCCAAGTGGGCAGTT"
    # The length of the sequence (or sequence alignment) vrna_fold_compound_t
    # n = vc.length
    n = len(struc)
    # vrna_mx_pf_t *matrices
    
    # probs = matrices.probs
    probs = [0] * n**2

    return vrna_db_from_probs(probs, n)

if __name__ == "__main__":
    print(generate_struc('s',1))

