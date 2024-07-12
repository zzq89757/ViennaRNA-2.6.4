import numpy as np

def generate_probs(vc, matrices):
    my_iindx = vc.iindx
    n = vc.length
    probs = matrices.probs
    
    
    for i in range(1, n + 1):
        probs[my_iindx[i] - i] = 0
    
    




def pf_create_bppm(vc, structure):
    n = vc.length
    pscore = vc.pscore if vc.type == "VRNA_FC_TYPE_COMPARATIVE" else None
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

        if vc.type == "VRNA_FC_TYPE_SINGLE":
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

    


def vrna_idx_row_wise(length):
    idx = []
    for i in range(1, length + 1):
        idx.append((((length + 1 - i) * (length - i)) // 2) + length + 1)
    print(idx)
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

