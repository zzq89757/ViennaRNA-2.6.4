from math import exp
from dimer_cofold import vrna_fold_compound_t, vrna_md_set_default, vrna_md_t




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


def main():
    AB:vrna_dimer_pf_t
    pairing_propensity:str
    # init exp params
    vrna_exp_params_rescale(vc:vrna_fold_compound_t, min_en:float)
    # obtain pairing_propensity  ----exp_matrices changed
    AB = vrna_pf_dimer(vc, pairing_propensity)
    # compute bpp ---vc not change
    prAB = vrna_plist_from_probs(vc, opt.bppmThreshold)
    # insert & character
    costruc = vrna_cut_point_insert(pairing_propensity, vc.cutpoint)
    print(costruc)

if __name__ == "__main__":
    main()
    