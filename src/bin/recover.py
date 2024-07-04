import math

########### global variable ################
#0 = BP; 1=any with GC; 2=any with AU-parameter If set to 1 or 2: fold sequences from an artificial alphabet ABCD..., where A pairs B, C pairs D, etc. using either GC (1) or AU parameters (2); default is 0, you probably don't want to change it.
energy_set = 0 
Law_and_Order = "_ACGUTXKI"
MIN2 = lambda a, b:(a + b)/2 - abs((a - b)/2)
NBASES = 8
NBPAIRS = 7
MAXALPHA = 20
MAXLOOP = 30
INF = 10000000
VRNA_GQUAD_MAX_STACK_SIZE = 7
VRNA_GQUAD_MAX_LINKER_LENGTH = 15
BP_pair = [
   # _  A  C  G  U  X  K  I
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 5, 0, 0, 5],
    [0, 0, 0, 1, 0, 0, 0, 0],
    [0, 0, 2, 0, 3, 0, 0, 0],
    [0, 6, 0, 4, 0, 0, 0, 6],
    [0, 0, 0, 0, 0, 0, 2, 0],
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 6, 0, 0, 5, 0, 0, 0]
]

rtype = [0, 2, 1, 4, 3, 6, 5, 7]
noGU = 1
nonstandards = None

# salt
VRNA_MODEL_DEFAULT_SALT = 1.021
K0 = 273.15
MAX_NINIO = 300
GASCONST = 1.98717
Eular_const = 0.58
PI = 3.141592653589793

######### model default params
VRNA_MODEL_DEFAULT_TEMPERATURE = 37.0
VRNA_MODEL_DEFAULT_PF_SCALE = -1
VRNA_MODEL_DEFAULT_BETA_SCALE = 1.0
VRNA_MODEL_DEFAULT_DANGLES = 2
VRNA_MODEL_DEFAULT_SPECIAL_HP = 1
VRNA_MODEL_DEFAULT_NO_LP = 0
VRNA_MODEL_DEFAULT_NO_GU = 0
VRNA_MODEL_DEFAULT_NO_GU_CLOSURE = 0
VRNA_MODEL_DEFAULT_CIRC = 0
VRNA_MODEL_DEFAULT_GQUAD = 0
VRNA_MODEL_DEFAULT_UNIQ_ML = 0
VRNA_MODEL_DEFAULT_ENERGY_SET = 0
VRNA_MODEL_DEFAULT_BACKTRACK = 1
VRNA_MODEL_DEFAULT_BACKTRACK_TYPE = 'F'
VRNA_MODEL_DEFAULT_COMPUTE_BPP = 1
VRNA_MODEL_DEFAULT_MAX_BP_SPAN = -1
VRNA_MODEL_DEFAULT_WINDOW_SIZE = -1
VRNA_MODEL_DEFAULT_LOG_ML = 0
VRNA_MODEL_DEFAULT_ALI_OLD_EN = 0
VRNA_MODEL_DEFAULT_ALI_RIBO = 0
VRNA_MODEL_DEFAULT_ALI_CV_FACT = 1.0
VRNA_MODEL_DEFAULT_ALI_NC_FACT = 1.0
VRNA_MODEL_DEFAULT_PF_SMOOTH = 1
VRNA_MODEL_DEFAULT_SALT = 1.021
VRNA_MODEL_DEFAULT_SALT_MLLOWER = 6
VRNA_MODEL_DEFAULT_SALT_MLUPPER = 24
VRNA_MODEL_DEFAULT_SALT_DPXINIT = 99999
VRNA_MODEL_SALT_DPXINIT_FACT_RNA = -45.324
VRNA_MODEL_SALT_DPXINIT_FACT_DNA = -58.389
VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT = VRNA_MODEL_SALT_DPXINIT_FACT_RNA
VRNA_MODEL_HELICAL_RISE_RNA = 2.8
VRNA_MODEL_HELICAL_RISE_DNA = 3.4
VRNA_MODEL_DEFAULT_HELICAL_RISE = VRNA_MODEL_HELICAL_RISE_RNA
VRNA_MODEL_BACKBONE_LENGTH_RNA = 6.0
VRNA_MODEL_BACKBONE_LENGTH_DNA = 6.76
TURN = 3
DM_DEFAULT = [
    [0, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 2, 1, 2, 2],
    [0, 2, 0, 1, 2, 2, 2],
    [0, 2, 1, 0, 2, 1, 2],
    [0, 1, 2, 2, 0, 2, 1],
    [0, 2, 2, 1, 2, 0, 2],
    [0, 2, 2, 2, 1, 2, 0]
]



class MFE:
    def __init__(self, i, j, energy, structure):
        self.i = i
        self.j = j
        self.energy = energy
        self.structure = structure
        
        
        
class vrna_md_t():
    def __init__(self):
        self.temperature = 37.0
        self.betaScale = 1.0
        self.pf_smooth = 0
        self.dangles = 2
        self.special_hp = 1
        self.noLP = 0
        self.noGU = 0
        self.noGUclosure = 0
        self.logML = 0
        self.circ = 0
        self.gquad = 0
        self.uniq_ML = 0
        self.energy_set = 0
        self.backtrack = 1
        self.backtrack_type = 'F'
        self.compute_bpp = 1
        self.nonstandards = [''] * 64
        self.max_bp_span = 30
        self.min_loop_size = 3
        self.window_size = -1
        self.oldAliEn = 0
        self.ribo = 0
        self.cv_fact = 1.0
        self.nc_fact = 1.0
        self.sfact = 1.07
        self.rtype = [0] * 8
        self.alias = [0] * (MAXALPHA + 1)
        self.pair = [[0] * (MAXALPHA + 1) for _ in range(MAXALPHA + 1)]
        self.pair_dist = [[0.0] * 7 for _ in range(7)]
        self.salt = 0.1
        self.saltMLLower = 3
        self.saltMLUpper = 30
        self.saltDPXInit = 0
        self.saltDPXInitFact = 0.0
        self.helical_rise = 2.8
        self.backbone_length = 6.0

        

class vrna_param_s():
    def __init__(self):
        self.id = 0
        self.stack = [[0] * (NBPAIRS + 1) for _ in range(NBPAIRS + 1)]
        self.hairpin = [0] * 31
        self.bulge = [0] * (MAXLOOP + 1)
        self.internal_loop = [0] * (MAXLOOP + 1)
        self.mismatchExt = [[[0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.mismatchI = [[[0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.mismatch1nI = [[[0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.mismatch23I = [[[0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.mismatchH = [[[0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.mismatchM = [[[0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.dangle5 = [[0] * 5 for _ in range(NBPAIRS + 1)]
        self.dangle3 = [[0] * 5 for _ in range(NBPAIRS + 1)]
        self.int11 = [[[[0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)] for _ in range(NBPAIRS + 1)]
        self.int21 = [[[[[0] * 5 for _ in range(5)] for _ in range(5)] for _ in range(NBPAIRS + 1)] for _ in range(NBPAIRS + 1)]
        self.int22 = [[[[[[0] * 5 for _ in range(5)] for _ in range(5)] for _ in range(5)] for _ in range(NBPAIRS + 1)] for _ in range(NBPAIRS + 1)]
        self.ninio = [0] * 5
        self.lxc = 0.0
        self.MLbase = 0
        self.MLintern = [0] * (NBPAIRS + 1)
        self.MLclosing = 0
        self.TerminalAU = 0
        self.DuplexInit = 0
        self.Tetraloop_E = [0] * 200
        self.Tetraloops = [''] * 1401
        self.Triloop_E = [0] * 40
        self.Triloops = [''] * 241
        self.Hexaloop_E = [0] * 40
        self.Hexaloops = [''] * 1801
        self.TripleC = 0
        self.MultipleCA = 0
        self.MultipleCB = 0
        self.gquad = [[0] * (3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1) for _ in range(VRNA_GQUAD_MAX_STACK_SIZE + 1)]
        self.gquadLayerMismatch = 0
        self.gquadLayerMismatchMax = 0
        self.temperature = 37.0
        self.model_details = vrna_md_t()
        self.param_file = ""
        self.SaltStack = 0
        self.SaltLoop = [0] * (MAXLOOP + 2)
        self.SaltLoopDbl = [0.0] * (MAXLOOP + 2)
        self.SaltMLbase = 0
        self.SaltMLintern = 0
        self.SaltMLclosing = 0
        self.SaltDPXInit = 0





class Duplex():
    def __init__(self, s1, s2) -> None:
        self.s1, self.s2 = s1, s2
        # 默认配对矩阵
        self.pair = [[BP_pair[i][j] for j in range(NBASES)] for i in range(NBASES)]
        self.n1 = len(s1)
        self.n2 = len(s2)

    ### encode 
    def encode_base(self, base:str):
        code = 0
        base = base.upper()
        if energy_set > 0:
            code =  ord(base) - ord('A') + 1
        else:
            idx = Law_and_Order.find(base)
            if idx > -1:code = idx
            if code > 5:code = 0
            if code > 4:code -= 1
        return code
        
    
    
    def encode_sequence(self, seq, how):
        s = [0]
        l = len(seq)
        s[0] = self.encode_base(seq[-1]) if how else l
        for i in range(1, l+1):
            s.append(self.encode_base(seq[i-1]))
        s.append(s[1])
        return s 
    
    #### salt
    def tau_ss(self, T, backbonelen):
        bjerrum_length_inv = 1 / self.bjerrum_length(T)
        return MIN2(1 / backbonelen, bjerrum_length_inv)
    
    def approx_hyper(self, y):
        a = 1 / (pow(y, 6.)/pow(2 * PI, 6.) + 1)
        b = pow(y, 4.)/(36 * pow(PI, 4.)) - pow(y, 3.) / (24 * PI * PI) + y * y / (2 * PI * PI) - y / 2
        c = math.log(2 * PI / y) - 1.96351
        return a * b + (1 - a) * c
    
    def loop_salt_aux(self, kmlss, L, T, backbonelen):
        expn = lambda n, k:math.exp(n * k) # unsure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        a = (GASCONST / 1000.) * T * self.bjerrum_length(T) * L * backbonelen * self.tau_ss(T, backbonelen) * self.tau_ss(T, backbonelen)
        b = math.log(kmlss) - math.log(PI / 2) + Eular_const + self.approx_hyper(kmlss) + 1 / kmlss * (1 - math.exp(-kmlss) + kmlss * expn(1, kmlss))
        return a * b * 100
    
    def vrna_salt_loop(self, L, rho, T, backbonelen):
        ionic_strength = lambda rho:rho
        epsilonr = lambda T:5321 / T + 233.76 - 0.9297 * T + 1.417 * T * T / 1000 - 0.8292 * T * T * T / 1000000
        self.bjerrum_length = lambda T:167100.052/(T * epsilonr(T))
        kappa = lambda rho, T:math.sqrt(self.bjerrum_length(T) * ionic_strength(rho)) / 8.1284
        if L == 0:return 0
        kmlss_ref = kappa(VRNA_MODEL_DEFAULT_SALT, T) * L * backbonelen
        kmlss = kappa(rho, T) * L * backbonelen
        correction = self.loop_salt_aux(kmlss, L, T, backbonelen) - self.loop_salt_aux(kmlss_ref, L, T, backbonelen)
        return correction
    
    
    def vrna_salt_loop_int(self, L, rho, T, backbonelen):
        correction = self.vrna_salt_loop(L, rho, T, backbonelen)
        return int(correction + 0.5 - (correction<0))
    
    def set_model_details(self, md:vrna_md_t):
        temperature     = VRNA_MODEL_DEFAULT_TEMPERATURE
        pf_scale        = VRNA_MODEL_DEFAULT_PF_SCALE
        dangles         = VRNA_MODEL_DEFAULT_DANGLES
        tetra_loop      = VRNA_MODEL_DEFAULT_SPECIAL_HP
        noLonelyPairs   = VRNA_MODEL_DEFAULT_NO_LP
        noGU            = VRNA_MODEL_DEFAULT_NO_GU
        no_closingGU    = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE
        circ            = VRNA_MODEL_DEFAULT_CIRC
        gquad           = VRNA_MODEL_DEFAULT_GQUAD
        uniq_ML         = VRNA_MODEL_DEFAULT_UNIQ_ML
        energy_set      = VRNA_MODEL_DEFAULT_ENERGY_SET
        do_backtrack    = VRNA_MODEL_DEFAULT_COMPUTE_BPP
        backtrack_type  = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE
        nonstandards    = None
        max_bp_span     = VRNA_MODEL_DEFAULT_MAX_BP_SPAN
        oldAliEn        = VRNA_MODEL_DEFAULT_ALI_OLD_EN
        ribo            = VRNA_MODEL_DEFAULT_ALI_RIBO
        cv_fact         = VRNA_MODEL_DEFAULT_ALI_CV_FACT
        nc_fact         = VRNA_MODEL_DEFAULT_ALI_NC_FACT
        logML           = VRNA_MODEL_DEFAULT_LOG_ML
        if md:
            md.dangles         = dangles
            md.special_hp      = tetra_loop
            md.noLP            = noLonelyPairs
            md.noGU            = noGU
            md.noGUclosure     = no_closingGU
            md.logML           = logML
            md.gquad           = gquad
            md.circ            = circ
            md.uniq_ML         = uniq_ML
            md.compute_bpp     = do_backtrack
            md.backtrack       = VRNA_MODEL_DEFAULT_BACKTRACK
            md.backtrack_type  = backtrack_type
            md.energy_set      = energy_set
            md.max_bp_span     = max_bp_span
            md.min_loop_size   = TURN
            md.window_size     = VRNA_MODEL_DEFAULT_WINDOW_SIZE
            md.oldAliEn        = oldAliEn
            md.ribo            = ribo
            md.cv_fact         = cv_fact
            md.nc_fact         = nc_fact
            md.temperature     = temperature
            md.betaScale       = VRNA_MODEL_DEFAULT_BETA_SCALE
            md.pf_smooth       = VRNA_MODEL_DEFAULT_PF_SMOOTH
            md.sfact           = 1.07
            md.salt            = self.P.model_details.salt
            md.saltMLLower     = self.P.model_details.saltMLLower
            md.saltMLUpper     = self.P.model_details.saltMLUpper
            md.saltDPXInit     = self.P.model_details.saltDPXInit
            md.saltDPXInitFact = self.P.model_details.saltDPXInitFact
            md.helical_rise    = self.P.model_details.helical_rise
            md.backbone_length = self.P.model_details.backbone_length 
            self.fill_pair_matrices(md)
      
            
    def fill_pair_matrices(self, md:vrna_md_t):
        dm_default = DM_DEFAULT
        def vrna_nucleotide_encode(c, md:vrna_md_t):
            code = -1

            c = c.upper()

            if md:
                if md.energy_set > 0:
                    code = ord(c) - ord('A') + 1
                else:
                    law_and_order = "ACGUT"
                    pos = law_and_order.find(c)
                    if pos == -1:
                        code = 0
                    else:
                        code = pos

                    if code > 4:
                        code -= 1  # make T and U equivalent

            return code

        def prepare_default_pairs(md:vrna_md_t):
            for i in range(5):
                md.alias[i] = i

            md.alias[5] = 3  # X <-> G
            md.alias[6] = 2  # K <-> C
            md.alias[7] = 0  # I <-> default base '@'

            for i in range(NBASES):
                for j in range(NBASES):
                    md.pair[i][j] = BP_pair[i][j]

            if md.noGU:
                md.pair[3][4] = md.pair[4][3] = 0

            if md.nonstandards:
                for i in range(0, len(md.nonstandards), 2):
                    md.pair[vrna_nucleotide_encode(md.nonstandards[i], md)][vrna_nucleotide_encode(md.nonstandards[i + 1], md)] = 7
        # nullify everything
        for i in range(MAXALPHA + 1):
            md.pair[i] = [0] * (MAXALPHA + 1)

        md.alias = [0] * (MAXALPHA + 1)

        # start setting actual base pair type encodings
        energy_set = md.energy_set
        if energy_set == 0:
            prepare_default_pairs(md)
        elif energy_set == 1:
            i = 1
            while i < MAXALPHA:
                md.alias[i] = 3  # A <-> G
                i += 1
                md.alias[i] = 2  # B <-> C
                i += 1

            i = 1
            while i < MAXALPHA:
                md.pair[i][i + 1] = 2  # AB <-> GC
                i += 1
                md.pair[i][i - 1] = 1  # BA <-> CG
                i += 1

        elif energy_set == 2:
            i = 1
            while i < MAXALPHA:
                md.alias[i] = 1  # A <-> A
                i += 1
                md.alias[i] = 4  # B <-> U
                i += 1

            i = 1
            while i < MAXALPHA:
                md.pair[i][i + 1] = 5  # AB <-> AU
                i += 1
                md.pair[i][i - 1] = 6  # BA <-> UA
                i += 1

        elif energy_set == 3:
            i = 1
            while i < MAXALPHA - 2:
                md.alias[i] = 3  # A <-> G
                i += 1
                md.alias[i] = 2  # B <-> C
                i += 1
                md.alias[i] = 1  # C <-> A
                i += 1
                md.alias[i] = 4  # D <-> U
                i += 1

            i = 1
            while i < MAXALPHA - 2:
                md.pair[i][i + 1] = 2  # AB <-> GC
                i += 1
                md.pair[i][i - 1] = 1  # BA <-> CG
                i += 1
                md.pair[i][i + 1] = 5  # CD <-> AU
                i += 1
                md.pair[i][i - 1] = 6  # DC <-> UA
                i += 1

        else:
            print("vrna_md_update: Unknown energy_set = {}. Using defaults!".format(md.energy_set))
            md.energy_set = 0
            prepare_default_pairs(md)

        # set the reverse base pair types
        for i in range(MAXALPHA + 1):
            for j in range(MAXALPHA + 1):
                md.rtype[md.pair[i][j]] = md.pair[j][i]

        # handle special cases separately
        md.rtype[0] = 0
        md.rtype[7] = 7

        for i in range(7):
            for j in range(7):
                md.pair_dist[i][j] = dm_default[i][j]
        
    
    def make_pair_matrix(self):
        # 默认能量集 (energy_set == 0)
        if energy_set == 0:
            # 将前5个碱基设置为其本身
            alias = list(range(5))
            alias += [3, 2, 0]  # X <-> G, K <-> C, I <-> default base '@'

            # 如果noGU为真，禁止G-U配对
            if noGU:
                self.pair[3][4] = self.pair[4][3] = 0

            # 如果nonstandards不为空，允许非标准碱基对
            if nonstandards is not None:
                for i in range(0, len(nonstandards), 2):
                    self.pair[self.encode_base(nonstandards[i])][self.encode_base(nonstandards[i + 1])] = 7

            # 设置反向配对矩阵
            # rtype = [[0 for _ in range(NBASES)] for _ in range(NBASES)]
            for i in range(NBASES):
                for j in range(NBASES):
                    rtype[self.pair[i][j]] = self.pair[j][i]

        else:
            self.pair = [[0 for _ in range(MAXALPHA + 1)] for _ in range(MAXALPHA + 1)]
            if energy_set == 1:
                alias = [0] + [3, 2] * ((MAXALPHA - 1) // 2)
                for i in range(1, MAXALPHA, 2):
                    self.pair[i][i + 1] = 2  # AB <-> GC
                    self.pair[i + 1][i] = 1  # BA <-> CG
            elif energy_set == 2:
                alias = [0] + [1, 4] * ((MAXALPHA - 1) // 2)
                for i in range(1, MAXALPHA, 2):
                    self.pair[i][i + 1] = 5  # AB <-> AU
                    self.pair[i + 1][i] = 6  # BA <-> UA
            elif energy_set == 3:
                alias = [0]
                for i in range(1, MAXALPHA - 2, 4):
                    alias += [3, 2, 1, 4]  # A <-> G, B <-> C, C <-> A, D <-> U
                for i in range(1, MAXALPHA - 2, 4):
                    self.pair[i][i + 1] = 2  # AB <-> GC
                    self.pair[i + 1][i] = 1  # BA <-> CG
                    self.pair[i + 2][i + 3] = 5  # CD <-> AU
                    self.pair[i + 3][i + 2] = 6  # DC <-> UA
            else:
                print("What energy_set are YOU using??")

            # 设置反向配对矩阵
            # rtype = [[0 for _ in range(MAXALPHA + 1)] for _ in range(MAXALPHA + 1)]
            for i in range(MAXALPHA + 1):
                for j in range(MAXALPHA + 1):
                    rtype[self.pair[i][j]] = self.pair[j][i]
    
        
    
    def vrna_E_ext_stem(self, type, n5d, n3d, P:vrna_param_s):
        energy = 0
        if n5d >= 0 and n3d >= 0:
            energy += P.mismatchExt[type][n5d][n3d]
        elif n5d >= 0:
            energy += P.dangle5[type][n5d]
        elif n3d >= 0:
            energy += P.dangle3[type][n3d]
        if type > 2:
            energy += P.TerminalAU
        return energy
    
    def E_IntLoop(self, n1, n2, type, type_2, si1, sj1, sp1, sq1, P:vrna_param_s):
        salt_stack_correction = P.SaltStack
        salt_loop_correction = 0

        if n1 > n2:
            nl = n1
            ns = n2
        else:
            nl = n2
            ns = n1

        if nl == 0:
            return P.stack[type][type_2] + salt_stack_correction  # stack

        backbones = nl + ns + 2

        if P.model_details.salt != VRNA_MODEL_DEFAULT_SALT:
            # salt correction for loop
            if backbones <= MAXLOOP + 1:
                salt_loop_correction = P.SaltLoop[backbones]
            else:
                salt_loop_correction = self.vrna_salt_loop_int(backbones, P.model_details.salt, P.temperature + K0, P.model_details.backbone_length)

        if ns == 0:
            # bulge
            energy = P.bulge[nl] if nl <= MAXLOOP else P.bulge[30] + int(P.lxc * math.log(nl / 30.))
            if nl == 1:
                energy += P.stack[type][type_2]
            else:
                if type > 2:
                    energy += P.TerminalAU
                if type_2 > 2:
                    energy += P.TerminalAU
        else:
            # interior loop
            if ns == 1:
                if nl == 1:  # 1x1 loop
                    return P.int11[type][type_2][si1][sj1] + salt_loop_correction
                if nl == 2:
                    # 2x1 loop
                    if n1 == 1:
                        energy = P.int21[type][type_2][si1][sq1][sj1]
                    else:
                        energy = P.int21[type_2][type][sq1][si1][sp1]
                    return energy + salt_loop_correction
                else:
                    # 1xn loop
                    energy = P.internal_loop[nl + 1] if nl + 1 <= MAXLOOP else P.internal_loop[30] + int(P.lxc * math.log((nl + 1) / 30.))
                    energy += min(MAX_NINIO, (nl - ns) * P.ninio[2])
                    energy += P.mismatch1nI[type][si1][sj1] + P.mismatch1nI[type_2][sq1][sp1]
                    return energy + salt_loop_correction
            elif ns == 2:
                if nl == 2:
                    # 2x2 loop
                    return P.int22[type][type_2][si1][sp1][sq1][sj1] + salt_loop_correction
                elif nl == 3:
                    # 2x3 loop
                    energy = P.internal_loop[5] + P.ninio[2]
                    energy += P.mismatch23I[type][si1][sj1] + P.mismatch23I[type_2][sq1][sp1]
                    return energy + salt_loop_correction

            # generic interior loop (no else here!)
            u = nl + ns
            energy = P.internal_loop[u] if u <= MAXLOOP else P.internal_loop[30] + int(P.lxc * math.log(u / 30.))
            energy += min(MAX_NINIO, (nl - ns) * P.ninio[2])
            energy += P.mismatchI[type][si1][sj1] + P.mismatchI[type_2][sq1][sp1]

        return energy + salt_loop_correction
    
    def backtrack(self, i, j):
        k = l = type = type2 = Energy = traced = i0 = j0 = 0
        
        st1 = ['.'] * (self.n1)
        st2 = ['.'] * (self.n2)

        i0 = min(i + 1, self.n1)
        j0 = max(j - 1, 1)

        while i > 0 and j <= self.n2:
            Energy = self.c[i][j]
            traced = False
            st1[i - 1] = '('
            st2[j - 1] = ')'
            type = self.pair[self.S1[i]][self.S2[j]]
            if not type:
                raise ValueError("backtrack failed in fold duplex")

            for k in range(i - 1, max(0, i - MAXLOOP - 2), -1):
                for l in range(j + 1, self.n2 + 1):
                    if i - k + l - j - 2 > MAXLOOP:
                        break

                    type2 = self.pair[self.S1[k]][self.S2[l]]
                    if not type2:
                        continue

                    LE = self.E_IntLoop(i - k - 1, l - j - 1, type2, rtype[type],
                                self.SS1[k + 1], self.SS2[l - 1], self.SS1[i - 1], self.SS2[j + 1], self.P)
                    if Energy == self.c[k][l] + LE:
                        traced = True
                        i = k
                        j = l
                        break
                if traced:
                    break

            if not traced:
                Energy -= self.vrna_E_ext_stem(type, self.SS1[i - 1] if i > 1 else -1, self.SS2[j + 1] if j < self.n2 else -1, self.P)
                if Energy != self.P.DuplexInit:
                    raise ValueError("backtrack failed in fold duplex")
                else:
                    break

        if i > 1:
            i -= 1

        if j < self.n2:
            j += 1

        struc = ''.join(st1[max(i - 1, 0):i0]) + '&' + ''.join(st2[j0 - 1:j])

        return struc
    
    def duplexfold_cu(self, clean_up):
        Emin = INF
        self.P = vrna_param_s()
        n1 = len(self.s1)
        n2 = len(self.s2)
        md = vrna_md_t()
        self.set_model_details(md)
        self.make_pair_matrix()
        self.c = [[0 for _ in range(n2 + 1)] for _ in range(n1 + 1)]
        self.S1 = self.encode_sequence(self.s1, 0)
        self.S2 = self.encode_sequence(self.s2, 0)
        self.SS1 = self.encode_sequence(self.s1, 1)
        self.SS2 = self.encode_sequence(self.s2, 1)
        for i in range(1, n1 + 1):
            for j in range(n2, 0 , -1):
                type = self.pair[self.S1[i]][self.S2[i]]
                self.c[i][j] = self.P.DuplexInit if type else INF
                if not type:
                    continue
                self.c[i][j] += self.vrna_E_ext_stem(type, self.SS1[i - 1] if i > 1 else -1, self.SS2[j + 1] if j < n2 else -1, self.P)
                for k in range(i - 1, 0, -1):
                    if i - k + n2 - j - 2 > MAXLOOP:
                        break
                    for l in range(j + 1, n2 + 1):
                        if i - k + l - j - 2 > MAXLOOP:
                            break
                        type2 = self.pair[self.S1[k]][self.S2[l]]
                        if not type2:
                            continue
                        Energy = self.E_IntLoop(i - k - 1, l - j - 1, type2, rtype[type], self.SS1[k + 1], self.SS2[l - 1], self.SS1[i - 1], self.SS2[j + 1], self.P)
                        self.c[i][j] = MIN2(self.c[i][j], self.c[k][l] + Energy)
                Energy = self.c[i][j]
                Energy += self.vrna_E_ext_stem(rtype[type], self.SS2[j - 1] if j > 1 else -1, self.SS1[i + 1] if i < n1 else -1, self.P)
                if Energy < Emin:
                    Emin = Energy
                    i_min = i
                    j_min = j
        struc = self.backtrack(i_min, j_min)
        if i_min < n1:
            i_min += 1

        if j_min > 1:
            j_min -= 1
        # 填充 MFE 结构体
        mfe = MFE(i_min, j_min, float(Emin) / 100., struc)

        # 清理内存
        if clean_up:
            self.c.clear()
            self.S1.clear()
            self.S2.clear()
            self.SS1.clear()
            self.SS2.clear()

        # 返回结果
        return mfe    
                

if __name__ == "__main__":
    d = Duplex("CTTCCTCGGGTTCAAAGCTGGATT","GTCCAGTTTTCCCAGGAAT")
    print(d.duplexfold_cu(0))
        