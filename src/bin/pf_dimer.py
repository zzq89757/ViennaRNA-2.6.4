from math import exp, log, sin
from dimer_cofold import vrna_fold_compound_t, vrna_md_set_default, vrna_md_t, vrna_exp_param_t, GASCONST, K0, VRNA_MODEL_DEFAULT_SALT, vrna_salt_loop
from remove import MAX_NINIO, NBPAIRS, VRNA_GQUAD_MAX_LINKER_LENGTH, VRNA_GQUAD_MAX_STACK_SIZE, VRNA_GQUAD_MIN_LINKER_LENGTH, VRNA_GQUAD_MIN_STACK_SIZE, Duplex, MAXLOOP, DuplexInit37, DuplexInitdH, GQuadAlpha37, GQuadAlphadH, GQuadBeta37, GQuadBetadH, GQuadLayerMismatch37, GQuadLayerMismatchH, GQuadLayerMismatchMax, Hexaloop37, HexaloopdH, Hexaloops, ML_BASEdH, ML_closing37, ML_closingdH, ML_intern37, ML_interndH, RESCALE_dG, TerminalAU37, ML_BASE37, TerminalAUdH, Tetraloop37, TetraloopdH, Tetraloops, Triloop37, TriloopdH, Triloops

## constant define 
Tmeasure = 310.15
lxc37 = 107.9
SCALE = 10
CLIP_NEGATIVE = lambda X:0 if X < 0 else X
def SMOOTH(X):
    if not pf_smooth:
        return CLIP_NEGATIVE(X)
    elif (X / SCALE) < -1.2283697:
        return 0
    elif (X / SCALE) > 0.8660254:
        return X
    else:
        factor = (sin((X / SCALE) - 0.34242663) + 1)
        return SCALE * 0.38490018 * factor * factor
TRUNC_MAYBE = lambda X:X if pf_smooth else float(int(X))
RESCALE_BF = lambda dG,dH,dT,kT: exp(-TRUNC_MAYBE(float(RESCALE_dG((dG), (dH), (dT))) * 10. / kT ))
RESCALE_BF_SMOOTH = lambda dG, dH, dT, kT:exp(SMOOTH(-TRUNC_MAYBE(float(RESCALE_dG(dG,dH,dT)))) * 10. / kT)


def get_scaled_exp_params(md:vrna_md_t, pfs:float) -> vrna_exp_param_t:
    pf = vrna_exp_param_t()
    pf.model_details = md
    pf.temperature   = md.temperature
    pf.alpha         = md.betaScale
    pf.kT            = kT = md.betaScale * (md.temperature + K0) * GASCONST # kT in cal/mol
    pf.pf_scale      = pfs
    pf_smooth         = md.pf_smooth
    TT                = (md.temperature + K0) / (Tmeasure)
    salt = md.salt
    saltStandard = VRNA_MODEL_DEFAULT_SALT

    pf.lxc                   = lxc37 * TT
    pf.expDuplexInit         = RESCALE_BF(DuplexInit37, DuplexInitdH, TT, kT)
    pf.expTermAU             = RESCALE_BF(TerminalAU37, TerminalAUdH, TT, kT)
    pf.expMLbase             = RESCALE_BF(ML_BASE37, ML_BASEdH, TT, kT)
    pf.expMLclosing          = RESCALE_BF(ML_closing37, ML_closingdH, TT, kT)
    pf.expgquadLayerMismatch = RESCALE_BF(GQuadLayerMismatch37, GQuadLayerMismatchH, TT, kT)
    pf.gquadLayerMismatchMax = GQuadLayerMismatchMax
    
    for i in range(VRNA_GQUAD_MIN_STACK_SIZE, VRNA_GQUAD_MAX_STACK_SIZE + 1):
        for j in range(3 * VRNA_GQUAD_MIN_LINKER_LENGTH, 3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1):
            GQuadAlpha_T:float = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, TT)
            GQuadBeta_T:float = RESCALE_dG(GQuadBeta37, GQuadBetadH, TT)
            GT = float(GQuadAlpha_T) * float(i - 1) + float(GQuadBeta_T) * log((float(j)) - 2.)
            pf.expgquad[i][j] = exp(-TRUNC_MAYBE(GT) * 10. / kT)
    
    for i in range(31):
        pf.exphairpin[i] = RESCALE_BF(hairpin37[i], hairpindH[i], TT, kT)
        
    for i in range(min(30, MAXLOOP) + 1):
        pf.expbulge[i]     = RESCALE_BF(bulge37[i], bulgedH[i], TT, kT)
        pf.expinternal[i]  = RESCALE_BF(internal_loop37[i], internal_loopdH[i], TT, kT)
        
    if salt == saltStandard:
        for i in range(min(30, MAXLOOP) + 1):
            pf.SaltLoopDbl[i] = 0.
            pf.expSaltLoop[i] = 1.

        for i in range(31, MAXLOOP + 1):
            pf.SaltLoopDbl[i] = 0.
            pf.expSaltLoop[i] = 1.
    else:
        for i in range(min(30, MAXLOOP) + 1):
            pf.SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md.backbone_length)
            saltLoop = int(pf.SaltLoopDbl[i] + 0.5 - (pf.SaltLoopDbl[i] < 0))
            pf.expSaltLoop[i] = exp(-saltLoop * 10. / kT)

        for i in range(31, MAXLOOP + 1):
            pf.SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md.backbone_length)
            saltLoop = int(pf.SaltLoopDbl[i] + 0.5 - (pf.SaltLoopDbl[i] < 0))
            pf.expSaltLoop[i] = exp(-saltLoop * 10. / kT)

    if james_rule:
        pf.expinternal[2] = exp(-80 * 10. / kT)

    GT = RESCALE_dG(bulge37[30], bulgedH[30], TT)
    for i in range(31, MAXLOOP + 1):
        pf.expbulge[i] = exp(-TRUNC_MAYBE(GT + (pf.lxc * log(i / 30.))) * 10. / kT)

    GT = RESCALE_dG(internal_loop37[30], internal_loopdH[30], TT)
    for i in range(31, MAXLOOP + 1):
        pf.expinternal[i] = exp(-TRUNC_MAYBE(GT + (pf.lxc * log(i / 30.))) * 10. / kT)

    GT = RESCALE_dG(ninio37, niniodH, TT)
    for j in range(MAXLOOP + 1):
        pf.expninio[2][j] = exp(-min(MAX_NINIO, j * TRUNC_MAYBE(GT)) * 10. / kT)

    for i in range(len(Tetraloops) // 7):
        pf.exptetra[i] = RESCALE_BF(Tetraloop37[i], TetraloopdH[i], TT, kT)

    for i in range(len(Triloops) // 5):
        pf.exptri[i] = RESCALE_BF(Triloop37[i], TriloopdH[i], TT, kT)

    for i in range(len(Hexaloops) // 9):
        pf.exphex[i] = RESCALE_BF(Hexaloop37[i], HexaloopdH[i], TT, kT)

    for i in range(NBPAIRS + 1):
        pf.expMLintern[i] = RESCALE_BF(ML_intern37, ML_interndH, TT, kT)

    for i in range(NBPAIRS + 1):
        for j in range(5):
            if md.dangles:
                pf.expdangle5[i][j] = RESCALE_BF_SMOOTH(dangle5_37[i][j], dangle5_dH[i][j], TT, kT)
                pf.expdangle3[i][j] = RESCALE_BF_SMOOTH(dangle3_37[i][j], dangle3_dH[i][j], TT, kT)
            else:
                pf.expdangle3[i][j] = pf.expdangle5[i][j] = 1

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            pf.expstack[i][j] = RESCALE_BF(stack37[i][j], stackdH[i][j], TT, kT)

    for i in range(NBPAIRS + 1):
        for j in range(5):
            for k in range(5):
                pf.expmismatchI[i][j][k] = RESCALE_BF(mismatchI37[i][j][k], mismatchIdH[i][j][k], TT, kT)
                pf.expmismatch1nI[i][j][k] = RESCALE_BF(mismatch1nI37[i][j][k], mismatch1nIdH[i][j][k], TT, kT)
                pf.expmismatchH[i][j][k] = RESCALE_BF(mismatchH37[i][j][k], mismatchHdH[i][j][k], TT, kT)
                pf.expmismatch23I[i][j][k] = RESCALE_BF(mismatch23I37[i][j][k], mismatch23IdH[i][j][k], TT, kT)

                if md.dangles:
                    pf.expmismatchM[i][j][k] = RESCALE_BF_SMOOTH(mismatchM37[i][j][k], mismatchMdH[i][j][k], TT, kT)
                    pf.expmismatchExt[i][j][k] = RESCALE_BF_SMOOTH(mismatchExt37[i][j][k], mismatchExtdH[i][j][k], TT, kT)
                else:
                    pf.expmismatchM[i][j][k] = pf.expmismatchExt[i][j][k] = 1

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            for k in range(5):
                for l in range(5):
                    pf.expint11[i][j][k][l] = RESCALE_BF(int11_37[i][j][k][l], int11_dH[i][j][k][l], TT, kT)

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            for k in range(5):
                for l in range(5):
                    for m in range(5):
                        pf.expint21[i][j][k][l][m] = RESCALE_BF(int21_37[i][j][k][l][m], int21_dH[i][j][k][l][m], TT, kT)

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            for k in range(5):
                for l in range(5):
                    for m in range(5):
                        for n in range(5):
                            pf.expint22[i][j][k][l][m][n] = RESCALE_BF(int22_37[i][j][k][l][m][n], int22_dH[i][j][k][l][m][n], TT, kT)

    pf.Tetraloops = Tetraloops[:281]
    pf.Triloops = Triloops[:241]
    pf.Hexaloops = Hexaloops[:361]
    pf.SaltMLbase = pf.SaltMLclosing = pf.SaltDPXInit = 0.
    if salt != saltStandard:
        pf.expSaltStack = exp(-vrna_salt_stack(salt, saltT, md.helical_rise) * 10. / kT)
        vrna_salt_ml(pf.SaltLoopDbl, md.saltMLLower, md.saltMLUpper, pf.SaltMLbase, pf.SaltMLclosing)

        if md.saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT:
            pf.SaltDPXInit = md.saltDPXInit
        elif md.saltDPXInit:
            pf.SaltDPXInit = vrna_salt_duplex_init(md)

        pf.expMLclosing *= exp(-pf.SaltMLbase * 10. / kT)
        pf.expMLclosing *= exp(-pf.SaltMLclosing * 10. / kT)
        pf.expMLbase *= exp(-pf.SaltMLbase * 10. / kT)

        for i in range(NBPAIRS + 1):
            pf.expMLintern[i] *= exp(-pf.SaltMLbase * 10. / kT)

        pf.expDuplexInit *= exp(-pf.SaltDPXInit * 10. / kT)

    return pf




def vrna_exp_params(md:vrna_md_t):
    if md:
        return get_scaled_exp_params(md, -1.)
    else:
        md = vrna_md_t
        vrna_md_set_default(md)
        return get_scaled_exp_params(md, -1.)


def rescale_params(vc:vrna_fold_compound_t):
    pf = vc.exp_params
    m = vc.exp_matrices
    
    if pf and m :
        m.scale[0] = 1.
        m.scale[1] = float(1. / pf.pf_scale)
        m.expMLbase[0] = 1
        m.expMLbase[1] = float(pf.expMLbase / pf.pf_scale)
        
        for i in range(2, vc.length + 1):
            m.scale[i] = m.scale[i / 2] * m.scale[i - (i / 2)]
            m.expMLbase = float(pow(pf.expMLbase, float(i)) * m.scale[i])

def vrna_exp_params_rescale(vc:vrna_fold_compound_t, mfe:float) -> None:
    if not vc:return
    if vc.exp_params:
        vc.exp_params = vrna_exp_params(vc.params.model_details)
    else:
        vc.exp_params.model_details = vc.params.model_details
    pf = vc.exp_params
    
    if pf:
        kT = pf.kT
        md = pf.model_details
        
        if mfe or pf.pf_scale < 1:
            if mfe:
                e_per_nt = mfe * 1000. / vc.length
            else:
                e_per_nt = -185 + (pf.temperature - 37.) * 7.27
            pf.pf_scale = exp(-(md.sfact * e_per_nt) / kT)
        if pf.pf_scale < 1.:
            pf.pf_scale = 1.
        
        rescale_params(vc)


# def decompose_pair(fc:vrna_fold_compound_t, i:int, j:int, aux_mx_ml:vrna_mx_pf_aux_ml_t) -> float:
def decompose_pair() -> float:
    return 0.


def exp_E_ml_fast(fc:vrna_fold_compound_t, i:int, j:int, aux_mx:vrna_mx_pf_aux_ml_s) -> float:
    from struc_utils import init_sc_mb_exp, prepare_hc_mb_def, vrna_get_ptype_md, exp_E_MLstem
    sliding_window = (fc.hc.type == VRNA_HC_WINDOW)
    n = fc.length
    sn = fc.strand_number
    ss = fc.strand_start
    se = fc.strand_end
    n_seq = 1 if (fc.type == VRNA_FC_TYPE_SINGLE) else fc.n_seq
    SS = None if (fc.type == VRNA_FC_TYPE_SINGLE) else fc.S
    S5 = None if (fc.type == VRNA_FC_TYPE_SINGLE) else fc.S5
    S3 = None if (fc.type == VRNA_FC_TYPE_SINGLE) else fc.S3
    iidx = None if sliding_window else fc.iindx
    ij = 0 if sliding_window else iidx[i] - j
    qqm = aux_mx.qqm
    qqm1 = aux_mx.qqm1
    qqmu = aux_mx.qqmu
    qm = None if sliding_window else fc.exp_matrices.qm
    qb = None if sliding_window else fc.exp_matrices.qb
    G = None if sliding_window else fc.exp_matrices.G
    qm_local = fc.exp_matrices.qm_local if sliding_window else None
    qb_local = fc.exp_matrices.qb_local if sliding_window else None
    G_local = fc.exp_matrices.G_local if sliding_window else None
    expMLbase = fc.exp_matrices.expMLbase
    pf_params = fc.exp_params
    md = pf_params.model_details
    hc = fc.hc
    domains_up = fc.domains_up
    circular = md.circ
    with_gquad = md.gquad
    with_ud = (domains_up and domains_up.exp_energy_cb)
    hc_up_ml = hc.up_ml
    evaluate = prepare_hc_mb_def(fc)

    sc_wrapper = init_sc_mb_exp(fc)

    qbt1 = 0
    q_temp = 0.0

    qqm[i] = 0.0

    if evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML):
        q_temp = qqm1[i] * expMLbase[1]

        if sc_wrapper.red_ml:
            q_temp *= sc_wrapper.red_ml(i, j, i, j - 1)

        qqm[i] += q_temp

    if with_ud:
        q_temp = 0.0

        for cnt in range(domains_up.uniq_motif_count):
            u = domains_up.uniq_motif_size[cnt]
            if j - u >= i:
                if evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML):
                    q_temp2 = (
                        qqmu[u][i]
                        * domains_up.exp_energy_cb(fc, j - u + 1, j, VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF, domains_up.data)
                        * expMLbase[u]
                    )

                    if sc_wrapper.red_ml:
                        q_temp2 *= sc_wrapper.red_ml(i, j, i, j - u)

                    q_temp += q_temp2

        qqm[i] += q_temp

    if evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM):
        qbt1 = qb_local[i][j] if sliding_window else qb[ij]

        if fc.type == VRNA_FC_TYPE_SINGLE:
            S1 = fc.sequence_encoding
            S2 = fc.sequence_encoding2
            type = vrna_get_ptype_md(S2[i], S2[j], md)

            qbt1 *= exp_E_MLstem(type, (S1[i - 1] if (i > 1 or circular) else -1), (S1[j + 1] if (j < n or circular) else -1), pf_params)

        elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
            q_temp = 1.0
            for s in range(n_seq):
                type = vrna_get_ptype_md(SS[s][i], SS[s][j], md)
                q_temp *= exp_E_MLstem(type, (S5[s][i] if (i > 1 or circular) else -1), (S3[s][j] if (j < n) or circular else -1), pf_params)

            qbt1 *= q_temp

        if sc_wrapper.red_stem:
            qbt1 *= sc_wrapper.red_stem(i, j, i, j)

        qqm[i] += qbt1

    if with_gquad:
        q_temp = G_local[i][j] if sliding_window else G[ij]
        qqm[i] += q_temp * pow(exp_E_MLstem(0, -1, -1, pf_params), n_seq)

    if with_ud:
        qqmu[0][i] = qqm[i]

    # Construction of qm matrix
    qqm_tmp = qqm

    if hc.f:
        qqm_tmp = [0.0] * (j - i + 2)
        qqm_tmp = [None] * i + qqm_tmp  # Adjust the index offset

        for k in range(j, i, -1):
            qqm_tmp[k] = qqm[k]
            if not evaluate(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML):
                qqm_tmp[k] = 0.0

    if sc_wrapper.decomp_ml:
        if qqm_tmp == qqm:
            qqm_tmp = [0.0] * (j - i + 2)
            qqm_tmp = [None] * i + qqm_tmp  # Adjust the index offset

            for k in range(j, i, -1):
                qqm_tmp[k] = qqm[k]

        for k in range(j, i, -1):
            qqm_tmp[k] *= sc_wrapper.decomp_ml(i, j, k - 1, k)

    temp = 0.0
    k = j

    if sliding_window:
        while k > i:
            temp += qm_local[i][k - 1] * qqm_tmp[k]
            k -= 1
    else:
        kl = iidx[i] - j + 1

        while True:
            stop = max(i, ss[sn[k]])
            while k > stop:
                temp += qm[kl] * qqm_tmp[k]
                k -= 1
                kl += 1

            if stop == i:
                break

            k -= 1
            kl += 1

    maxk = j

    if maxk > i + hc_up_ml[i]:
        maxk = i + hc_up_ml[i]

    if maxk > se[sn[i]]:
        maxk = se[sn[i]]

    if qqm_tmp != qqm:
        for k in range(maxk, i, -1):
            qqm_tmp[k] = qqm[k]

    if hc.f:
        if qqm_tmp == qqm:
            qqm_tmp = [0.0] * (j - i + 2)
            qqm_tmp = [None] * i + qqm_tmp  # Adjust the index offset

            for k in range(maxk, i, -1):
                qqm_tmp[k] = qqm[k]

        for k in range(maxk, i, -1):
            if not evaluate(i, j, k, j, VRNA_DECOMP_ML_ML):
                qqm_tmp[k] = 0.0

    if sc_wrapper.red_ml:
        if qqm_tmp == qqm:
            qqm_tmp = [0.0] * (j - i + 2)
            qqm_tmp = [None] * i + qqm_tmp  # Adjust the index offset

            for k in range(maxk, i, -1):
                qqm_tmp[k] = qqm[k]

        for k in range(maxk, i, -1):
            qqm_tmp[k] *= sc_wrapper.red_ml(i, j, k, j)

    ii = maxk - i

    for k in range(maxk, i, -1):
        temp += expMLbase[ii] * qqm_tmp[k]
        ii -= 1

    if with_ud:
        ii = maxk - i

        for k in range(maxk, i, -1):
            temp += expMLbase[ii] * qqm_tmp[k] * domains_up.exp_energy_cb(fc, i, k - 1, VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP, domains_up.data)
            ii -= 1

    if qqm_tmp != qqm:
        del qqm_tmp[i:]

    if fc.aux_grammar and fc.aux_grammar.cb_aux_exp_m:
        temp += fc.aux_grammar.cb_aux_exp_m(fc, i, j, fc.aux_grammar.data)

    # free_sc_mb_exp(sc_wrapper)

    return temp + qqm[i]
    


def vrna_exp_E_ml_fast(fc:vrna_fold_compound_t, i:int, j:int, aux_mx:vrna_mx_pf_aux_ml_s) -> float:
    q = 0.
    if fc and aux_mx:
        q = exp_E_ml_fast(fc, i, j, aux_mx)
    return q


def vrna_exp_E_ext_fast_init(fc:vrna_fold_compound_t) -> vrna_mx_pf_aux_el_s:
    from struc_utils import init_sc_ext_exp, prepare_hc_ext_def, sc_ext_exp_dat, hc_ext_def_dat
    hc_dat_local:hc_ext_def_dat
    sc_wrapper:sc_ext_exp_dat
    if fc is None:
        return None

    n = int(fc.length)
    iidx = fc.iindx
    turn = fc.exp_params.model_details.min_loop_size
    domains_up = fc.domains_up
    with_ud = domains_up and domains_up.exp_energy_cb

    if fc.hc.type == VRNA_HC_WINDOW:
        evaluate = prepare_hc_ext_def_window(fc, hc_dat_local)
    else:
        evaluate = prepare_hc_ext_def(fc, hc_dat_local)

    init_sc_ext_exp(fc, sc_wrapper)

    # Allocate memory for helper arrays
    aux_mx = vrna_mx_pf_aux_el_s()

    # Pre-processing ligand binding production rule(s) and auxiliary memory
    if with_ud:
        ud_max_size = max(domains_up.uniq_motif_size)
        aux_mx.qqu_size = ud_max_size
        aux_mx.qqu = [[0.0] * (n + 2) for _ in range(ud_max_size + 1)]

    if fc.hc.type == VRNA_HC_WINDOW:
        q_local = fc.exp_matrices.q_local
        max_j = min(turn + 1, fc.window_size, n)
        for j in range(1, max_j + 1):
            for i in range(1, j + 1):
                q_local[i][j] = reduce_ext_up_fast(fc, i, j, aux_mx, evaluate, hc_dat_local, sc_wrapper)
    else:
        q = fc.exp_matrices.q
        for d in range(turn + 1):
            for i in range(1, n - d + 1):
                j = i + d
                ij = iidx[i] - j
                q[ij] = reduce_ext_up_fast(fc, i, j, aux_mx, evaluate, hc_dat_local, sc_wrapper)

        if fc.aux_grammar and fc.aux_grammar.cb_aux_exp_f:
            for d in range(turn + 1):
                for i in range(1, n - d + 1):
                    j = i + d
                    ij = iidx[i] - j
                    q[ij] += fc.aux_grammar.cb_aux_exp_f(fc, i, j, fc.aux_grammar.data)

    return aux_mx

def vrna_exp_E_ml_fast_init(fc:vrna_fold_compound_t) -> vrna_mx_pf_aux_ml_s:
    if fc is None:
        return None

    n = int(fc.length)
    iidx = fc.iindx
    turn = fc.exp_params.model_details.min_loop_size
    qm = fc.exp_matrices.qm

    # Allocate memory for helper arrays
    aux_mx = vrna_mx_pf_aux_ml_s()

    if fc.type == VRNA_FC_TYPE_SINGLE:
        domains_up = fc.domains_up
        with_ud = domains_up and domains_up.exp_energy_cb
        ud_max_size = 0

        # Pre-processing ligand binding production rule(s) and auxiliary memory
        if with_ud:
            for u in domains_up.uniq_motif_size:
                if ud_max_size < u:
                    ud_max_size = u

            aux_mx.qqmu_size = ud_max_size
            aux_mx.qqmu = [[0.0] * (n + 2) for _ in range(ud_max_size + 1)]

    if fc.hc.type != VRNA_HC_WINDOW:
        for d in range(turn + 1):
            for i in range(1, n - d + 1):
                j = i + d
                ij = iidx[i] - j

                if j > n:
                    continue

                qm[ij] = 0.0

        if fc.aux_grammar and fc.aux_grammar.cb_aux_exp_m:
            for d in range(turn + 1):
                for i in range(1, n - d + 1):
                    j = i + d
                    ij = iidx[i] - j

                    if j > n:
                        continue

                    qm[ij] += fc.aux_grammar.cb_aux_exp_m(fc, i, j, fc.aux_grammar.data)

    return aux_mx
  

def fill_arrays(fc:vrna_fold_compound_t) -> int:
    n           = fc.length
    my_iindx    = fc.iindx
    jindx       = fc.jindx
    matrices    = fc.exp_matrices
    pf_params   = fc.exp_params
    domains_up  = fc.domains_up
    q           = matrices.q
    qb          = matrices.qb
    qm          = matrices.qm
    qm1         = matrices.qm1
    q1k         = matrices.q1k
    qln         = matrices.qln
    md          = pf_params.model_details
    with_gquad  = md.gquad
    
    with_ud = (domains_up and domains_up.exp_energy_cb)
    Qmax    = 0
    
    aux_mx_el = vrna_exp_E_ext_fast_init(fc)
    aux_mx_ml = vrna_exp_E_ml_fast_init(fc)
    
    for i in range(1, n + 1):
        ij      = my_iindx[i] - i
        qb[ij]  = 0.0
    
    for j in range(2, n + 1):
        for i in range(j - 1, 0, -1):
            ij = my_iindx[i] - j
            qb[ij] = decompose_pair(fc, i, j, aux_mx_ml)
            qm[ij] = vrna_exp_E_ml_fast(fc, i, j, aux_mx_ml)
            q[ij] = vrna_exp_E_ext_fast(fc, i, j, aux_mx_el) 
            if q[ij] > Qmax:Qmax = q[ij]
        vrna_exp_E_ext_fast_rotate(aux_mx_el)
        vrna_exp_E_ml_fast_rotate(aux_mx_ml)
    
    if q1k and qln:
        for k in range(1,n + 1):
            q1k[k]  = q[my_iindx[1] - k]
            qln[k]  = q[my_iindx[k] - n]
        q1k[0]      = 1.0
        qln[n + 1]  = 1.0
        
    # vrna_exp_E_ml_fast_free(aux_mx_ml)
    # vrna_exp_E_ext_fast_free(aux_mx_el)
    
    return 1
    
    

def vrna_pf(fc:vrna_fold_compound_t, structure:str) -> float:
    dG = float(10000000/100.)
    if fc:
        n         = fc.length
        params    = fc.exp_params
        matrices  = fc.exp_matrices
        md        = params.model_details

        # Explicitly turn off dynamic threads
        # if 'OMP' in globals():
        omp_set_dynamic(0)

        # Set appropriate arithmetic mode
        # if 'SUN4' in globals():
        #     nonstandard_arithmetic()
        # elif 'HP9' in globals():
        #     fpsetfastmode(1)

        # Call user-defined recursion status callback function
        if fc.stat_cb: # not into
            fc.stat_cb(VRNA_STATUS_PF_PRE, fc.auxdata)

        # Prepare multi-strand folding
        if fc.strands > 1: # aux_grammar changed
            vrna_pf_multifold_prepare(fc)

        # Call user-defined grammar pre-condition callback function
        if fc.aux_grammar and fc.aux_grammar.cb_proc: # not into
            fc.aux_grammar.cb_proc(fc, VRNA_STATUS_PF_PRE, fc.aux_grammar.data)

        # Fill arrays
        if not fill_arrays(fc): # not into
            # if 'SUN4' in globals():
            #     standard_arithmetic()
            # elif 'HP9' in globals():
            #     fpsetfastmode(0)
            return dG

        if md.circ: # not into
            # Post-process step for circular RNAs
            postprocess_circular(fc)

        # Call user-defined grammar post-condition callback function
        if fc.aux_grammar and fc.aux_grammar.cb_proc: # not into 
            fc.aux_grammar.cb_proc(fc, VRNA_STATUS_PF_POST, fc.aux_grammar.data)

        if fc.strands > 1:
            vrna_gr_reset(fc)

        # Call user-defined recursion status callback function
        if fc.stat_cb: # not into
            fc.stat_cb(VRNA_STATUS_PF_POST, fc.auxdata)

        # Calculate Q based on backtrack type
        if md.backtrack_type == 'C':
            Q = matrices.qb[fc.iindx[1] - n]
        elif md.backtrack_type == 'M':
            Q = matrices.qm[fc.iindx[1] - n]
        else:
            Q = matrices.qo if md.circ else matrices.q[fc.iindx[1] - n]

        # Ensemble free energy in Kcal/mol
        if Q <= 1.17549435082228750796873653722224568e-38:
            print("pf_scale too large")

        if fc.strands > 1:
            # Check for rotational symmetry correction
            sym = vrna_rotational_symmetry(fc.sequence)
            Q /= sym

            # Add interaction penalty
            Q *= power(params.expDuplexInit, fc.strands - 1)

        dG = (-log(Q) - n * log(params.pf_scale)) * params.kT / 1000.0

        # Calculate base pairing probability matrix (bppm)
        if md.compute_bpp: # into
            vrna_pairing_probs(fc, structure)

            # Backward compatibility block
            if 'VRNA_DISABLE_BACKWARD_COMPATIBILITY' not in globals(): # into
                pr = matrices.probs
    return dG

def vrna_rotational_symmetry_pos(string:str, positions:int) -> int:
    shift_size:int
    if not string:
        if positions:
            positions = None
        return 0
    string_length = len(string)
    if string_length == 0:
        if positions:
            positions = None
        return 0
    matches = 1
    if positions:
        shift_size = 10
    


def vrna_rotational_symmetry(string:str):
    return vrna_rotational_symmetry_pos(string, None)



def extract_dimer_props(fc:vrna_fold_compound_t, F0AB, FAB, FcAB, FA, FB):
    """
    Extract dimer properties for a given fold compound.

    :param fc: The fold compound object
    :param F0AB: Null model without DuplexInit (output)
    :param FAB: All states with DuplexInit correction (output)
    :param FcAB: True hybrid states only (output)
    :param FA: Monomer A (output)
    :param FB: Monomer B (output)
    """
    n = fc.length
    sym = None
    kT = fc.exp_params.kT / 1000.0
    QAB = fc.exp_matrices.q[fc.iindx[1] - n]

    if fc.strands > 1:
        # Check for rotational symmetry correction
        sym = vrna_rotational_symmetry(fc.sequence)
        QAB /= float(sym)

        # Add interaction penalty
        QAB *= pow(fc.exp_params.expDuplexInit, float(fc.strands - 1))

        Qzero = (fc.exp_matrices.q[fc.iindx[1] - n] +
                 fc.exp_matrices.q[fc.iindx[1] - fc.strand_end[fc.strand_order[0]]] *
                 fc.exp_matrices.q[fc.iindx[fc.strand_start[fc.strand_order[1]]] - n])

        QToT = (fc.exp_matrices.q[fc.iindx[1] - fc.strand_end[fc.strand_order[0]]] *
                fc.exp_matrices.q[fc.iindx[fc.strand_start[fc.strand_order[1]]] - n] +
                QAB)

        FAB = -kT * (log(QToT) + n * log(fc.exp_params.pf_scale))
        F0AB = -kT * (log(Qzero) + n * log(fc.exp_params.pf_scale))
        FcAB = -kT * (log(QAB) + n * log(fc.exp_params.pf_scale)) if QAB > 1e-17 else 999
        FA = -kT * (log(fc.exp_matrices.q[fc.iindx[1] - fc.strand_end[fc.strand_order[0]]]) +
                        fc.strand_end[fc.strand_order[0]] * log(fc.exp_params.pf_scale))
        FB = -kT * (log(fc.exp_matrices.q[fc.iindx[fc.strand_start[fc.strand_order[1]]] - n]) +
                        (n - fc.strand_start[fc.strand_order[1]] + 1) * log(fc.exp_params.pf_scale))
    else:
        FA = FB = FAB = F0AB = (-log(fc.exp_matrices.q[fc.iindx[1] - n]) - n * log(fc.exp_params.pf_scale)) * fc.exp_params.kT / 1000.0
        FcAB = 0



def vrna_pf_dimer(fc:vrna_fold_compound_t|None, structure:str) -> vrna_dimer_pf_t:
    X:vrna_dimer_pf_t
    X.F0AB = X.FAB = X.FcAB = X.FA = X.FB = 0.
    if fc:
        vrna_pf(fc, structure)
        extract_dimer_props(fc,
                            X.F0AB,
                            X.FAB,
                            X.FcAB,
                            X.FA,
                            X.FB)
    return X
        
    
    

def main():
    AB:vrna_dimer_pf_t
    pairing_propensity:str
    mfe_structure:str
    # vc init by vrna_fold_compound
    
    
    # vc.matrices changed by vrna_mfe_dimer(no erro after remove vc.matrices)
    # min_en  = vrna_mfe_dimer(vc, mfe_structure)
    min_en = 
    # init vc  exp params
    vrna_exp_params_rescale(vc, min_en)
    # obtain pairing_propensity  ----exp_matrices changed
    AB = vrna_pf_dimer(vc, pairing_propensity)
    # compute bpp ---vc not change
    # prAB = vrna_plist_from_probs(vc, opt.bppmThreshold)
    # insert & character
    costruc = vrna_cut_point_insert(pairing_propensity, vc.cutpoint)
    print(costruc)

if __name__ == "__main__":
    main()
    