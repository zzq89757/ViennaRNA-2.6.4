
def generate_probs(vc, matrices):
    my_iindx = vc.iindx
    n = vc.length
    probs = matrices.probs
    
    
    for i in range(1, n + 1):
        probs[my_iindx[i] - i] = 0
    
    
    


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

