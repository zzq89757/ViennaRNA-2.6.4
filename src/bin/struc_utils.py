import math
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
VRNA_CONSTRAINT_CONTEXT_INT_LOOP = 0x04
VRNA_DECOMP_PAIR_ML_EXT = 23
VRNA_DECOMP_ML_ML_ML = 5
VRNA_DECOMP_ML_ML_STEM = 9
VRNA_GQUAD_MAX_STACK_SIZE = 7
VRNA_GQUAD_MIN_STACK_SIZE = 2
VRNA_GQUAD_MAX_LINKER_LENGTH = 15
VRNA_GQUAD_MIN_LINKER_LENGTH = 1
NBPAIRS = 7
VRNA_GQUAD_MIN_BOX_SIZE = ((4 * VRNA_GQUAD_MIN_STACK_SIZE) + (3 * VRNA_GQUAD_MIN_LINKER_LENGTH))
VRNA_GQUAD_MAX_BOX_SIZE = ((4 * VRNA_GQUAD_MAX_STACK_SIZE) + (3 * VRNA_GQUAD_MAX_LINKER_LENGTH))


VRNA_DECOMP_EXT_EXT_STEM = 18
VRNA_CONSTRAINT_CONTEXT_EXT_LOOP = 0x01
VRNA_DECOMP_EXT_STEM_EXT = 16
VRNA_DECOMP_EXT_EXT_STEM1 = 19
VRNA_DECOMP_EXT_STEM_EXT1 = 20
VRNA_DECOMP_EXT_EXT_EXT = 15
VRNA_DECOMP_EXT_STEM = 14
VRNA_DECOMP_EXT_EXT = 12
VRNA_DECOMP_EXT_UP = 13
VRNA_DECOMP_EXT_STEM_OUTSIDE = 17
VRNA_DECOMP_ML_COAXIAL_ENC  = 11
VRNA_DECOMP_ML_COAXIAL = 10

# data structure/class start ######################################
from enum import Enum
class vrna_fc_type_e(Enum):
    VRNA_FC_TYPE_SINGLE = 0
    VRNA_FC_TYPE_COMPARATIVE = 1

from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Callable

class vrna_seq_t:
    def __init__(self,
                 seq_type: int,               # Assuming vrna_seq_type_e is an integer or enum
                 name: Optional[str] = None,
                 string: Optional[str] = None,
                 encoding: Optional[List[int]] = None,
                 encoding5: Optional[List[int]] = None,
                 encoding3: Optional[List[int]] = None,
                 length: int = 0):
        self.type = seq_type
        self.name = name
        self.string = string
        self.encoding = encoding if encoding is not None else []
        self.encoding5 = encoding5 if encoding5 is not None else []
        self.encoding3 = encoding3 if encoding3 is not None else []
        self.length = length


class vrna_hc_t:
    def __init__(self) -> None:
        self.type = vrna_fc_type_e.VRNA_FC_TYPE_COMPARATIVE.value
        self.mx = ''
        self.up_ext = self.up_hp = self.up_int = self.up_ml = 0
        self.f = vrna_hc_eval_f()
        self.data = None
     
     
class vrna_sc_type_e(Enum):
    VRNA_SC_DEFAULT = 0
    VRNA_SC_WINDOW  = 1


    


class vrna_sc_t:
    def __init__(self) -> None:
        self.type = vrna_sc_type_e.VRNA_SC_DEFAULT.value
        self.n = 0
        self.state = 0  # unsigned char state

        self.energy_up = None  # int **energy_up
        self.exp_energy_up = None  # FLT_OR_DBL **exp_energy_up

        self.up_storage = None  # int *up_storage
        # self.bp_storage = vrna_sc_bp_storage_t  # vrna_sc_bp_storage_t **bp_storage

        self.energy_bp = None  # int *energy_bp
        self.exp_energy_bp = None  # FLT_OR_DBL *exp_energy_bp

        self.energy_bp_local = None  # int **energy_bp_local
        self.exp_energy_bp_local = None  # FLT_OR_DBL **exp_energy_bp_local

        self.energy_stack = None  # int *energy_stack
        self.exp_energy_stack = None  # FLT_OR_DBL *exp_energy_stack

        self.f = vrna_sc_f  # vrna_sc_f f
        self.bt = vrna_sc_bt_f  # vrna_sc_bt_f bt
        self.exp_f = vrna_sc_exp_f  # vrna_sc_exp_f exp_f

        self.data = None  # void *data
        # self.prepare_data = vrna_auxdata_prepare_f  # vrna_auxdata_prepare_f prepare_data
        # self.free_data = vrna_auxdata_free_f  # vrna_auxdata_free_f free_data
class vrna_basepair_t:
    def __init__(self) -> None:
        self.i = self.j = 0

vrna_sc_f = Callable[[int, int, int, int, str, None], int]
vrna_sc_bt_f = Callable[[int, int, int, int, str, None], vrna_basepair_t]


class vrna_mx_mfe_t:
    def __init__(self):
        self.type: vrna_mx_type_e.VRNA_MX_DEFAULT.value
        self.length: int
        self.strands: int

        # Default DP matrices
        self.c: Optional[List[int]] = None
        self.f5: Optional[List[int]] = None
        self.f3: Optional[List[int]] = None
        self.fms5: Optional[List[List[int]]] = None
        self.fms3: Optional[List[List[int]]] = None
        self.fML: Optional[List[int]] = None
        self.fM1: Optional[List[int]] = None
        self.fM2: Optional[List[int]] = None
        self.ggg: Optional[List[int]] = None
        self.Fc: Optional[int] = None
        self.FcH: Optional[int] = None
        self.FcI: Optional[int] = None
        self.FcM: Optional[int] = None

        # Local Folding DP matrices using window approach
        self.c_local: Optional[List[List[int]]] = None
        self.f3_local: Optional[List[int]] = None
        self.fML_local: Optional[List[List[int]]] = None
        self.ggg_local: Optional[List[List[int]]] = None

        # Distance Class DP matrices
        self.E_F5: Optional[List[List[List[int]]]] = None
        self.l_min_F5: Optional[List[List[int]]] = None
        self.l_max_F5: Optional[List[List[int]]] = None
        self.k_min_F5: Optional[List[int]] = None
        self.k_max_F5: Optional[List[int]] = None

        self.E_F3: Optional[List[List[List[int]]]] = None
        self.l_min_F3: Optional[List[List[int]]] = None
        self.l_max_F3: Optional[List[List[int]]] = None
        self.k_min_F3: Optional[List[int]] = None
        self.k_max_F3: Optional[List[int]] = None

        self.E_C: Optional[List[List[List[int]]]] = None
        self.l_min_C: Optional[List[List[int]]] = None
        self.l_max_C: Optional[List[List[int]]] = None
        self.k_min_C: Optional[List[int]] = None
        self.k_max_C: Optional[List[int]] = None

        self.E_M: Optional[List[List[List[int]]]] = None
        self.l_min_M: Optional[List[List[int]]] = None
        self.l_max_M: Optional[List[List[int]]] = None
        self.k_min_M: Optional[List[int]] = None
        self.k_max_M: Optional[List[int]] = None

        self.E_M1: Optional[List[List[List[int]]]] = None
        self.l_min_M1: Optional[List[List[int]]] = None
        self.l_max_M1: Optional[List[List[int]]] = None
        self.k_min_M1: Optional[List[int]] = None
        self.k_max_M1: Optional[List[int]] = None

        self.E_M2: Optional[List[List[List[int]]]] = None
        self.l_min_M2: Optional[List[List[int]]] = None
        self.l_max_M2: Optional[List[List[int]]] = None
        self.k_min_M2: Optional[List[int]] = None
        self.k_max_M2: Optional[List[int]] = None

        self.E_Fc: Optional[List[List[int]]] = None
        self.l_min_Fc: Optional[List[int]] = None
        self.l_max_Fc: Optional[List[int]] = None
        self.k_min_Fc: Optional[int] = None
        self.k_max_Fc: Optional[int] = None

        self.E_FcH: Optional[List[List[int]]] = None
        self.l_min_FcH: Optional[List[int]] = None
        self.l_max_FcH: Optional[List[int]] = None
        self.k_min_FcH: Optional[int] = None
        self.k_max_FcH: Optional[int] = None

        self.E_FcI: Optional[List[List[int]]] = None
        self.l_min_FcI: Optional[List[int]] = None
        self.l_max_FcI: Optional[List[int]] = None
        self.k_min_FcI: Optional[int] = None
        self.k_max_FcI: Optional[int] = None

        self.E_FcM: Optional[List[List[int]]] = None
        self.l_min_FcM: Optional[List[int]] = None
        self.l_max_FcM: Optional[List[int]] = None
        self.k_min_FcM: Optional[int] = None
        self.k_max_FcM: Optional[int] = None

        # Auxiliary arrays for remaining set of coarse graining (k,l) > (k_max, l_max)
        self.E_F5_rem: Optional[List[int]] = None
        self.E_F3_rem: Optional[List[int]] = None
        self.E_C_rem: Optional[List[int]] = None
        self.E_M_rem: Optional[List[int]] = None
        self.E_M1_rem: Optional[List[int]] = None
        self.E_M2_rem: Optional[List[int]] = None

        self.E_Fc_rem: Optional[int] = None
        self.E_FcH_rem: Optional[int] = None
        self.E_FcI_rem: Optional[int] = None
        self.E_FcM_rem: Optional[int] = None

        # Additional fields
        # Uncomment the following if COUNT_STATES is defined
        # self.N_F5: Optional[List[List[List[int]]]] = None
        # self.N_C: Optional[List[List[List[int]]]] = None
        # self.N_M: Optional[List[List[List[int]]]] = None
        # self.N_M1: Optional[List[List[List[int]]]] = None

class vrna_mx_type_e(Enum):
    VRNA_MX_DEFAULT = 0,
    VRNA_MX_WINDOW = 1,
    VRNA_MX_2DFOLD = 2

class vrna_mx_pf_t:
    def __init__(self):
        self.type = vrna_mx_type_e.VRNA_MX_DEFAULT.value
        self.length = 0
        self.scale = 0.
        self.expMLbase = 0.
        
        # Default PF matrices
        self.q = None
        self.qb = None
        self.qm = None
        self.qm1 = None
        self.probs = None
        self.q1k = None
        self.qln = None
        self.G = None
        
        self.qo = None
        self.qm2 = None
        self.qho = None
        self.qio = None
        self.qmo = None
        
        # Local Folding DP matrices using window approach
        self.q_local = None
        self.qb_local = None
        self.qm_local = None
        self.pR = None
        self.qm2_local = None
        self.QI5 = None
        self.q2l = None
        self.qmb = None
        self.G_local = None
        
        # Distance Class DP matrices
        self.Q = None
        self.l_min_Q = None
        self.l_max_Q = None
        self.k_min_Q = None
        self.k_max_Q = None
        
        self.Q_B = None
        self.l_min_Q_B = None
        self.l_max_Q_B = None
        self.k_min_Q_B = None
        self.k_max_Q_B = None
        
        self.Q_M = None
        self.l_min_Q_M = None
        self.l_max_Q_M = None
        self.k_min_Q_M = None
        self.k_max_Q_M = None
        
        self.Q_M1 = None
        self.l_min_Q_M1 = None
        self.l_max_Q_M1 = None
        self.k_min_Q_M1 = None
        self.k_max_Q_M1 = None
        
        self.Q_M2 = None
        self.l_min_Q_M2 = None
        self.l_max_Q_M2 = None
        self.k_min_Q_M2 = None
        self.k_max_Q_M2 = None
        
        self.Q_c = None
        self.l_min_Q_c = None
        self.l_max_Q_c = None
        self.k_min_Q_c = None
        self.k_max_Q_c = None
        
        self.Q_cH = None
        self.l_min_Q_cH = None
        self.l_max_Q_cH = None
        self.k_min_Q_cH = None
        self.k_max_Q_cH = None
        
        self.Q_cI = None
        self.l_min_Q_cI = None
        self.l_max_Q_cI = None
        self.k_min_Q_cI = None
        self.k_max_Q_cI = None
        
        self.Q_cM = None
        self.l_min_Q_cM = None
        self.l_max_Q_cM = None
        self.k_min_Q_cM = None
        self.k_max_Q_cM = None
        
        self.Q_rem = None
        self.Q_B_rem = None
        self.Q_M_rem = None
        self.Q_M1_rem = None
        self.Q_M2_rem = None
        
        self.Q_c_rem = None
        self.Q_cH_rem = None
        self.Q_cI_rem = None
        self.Q_cM_rem = None

        
        
class vrna_exp_param_t:
    def __init__(self) -> None:
        self.id: int = 0  # An identifier for the data structure (deprecated)
        self.expstack: List[List[float]] = [[0.0] * (NBPAIRS + 1) for _ in range(NBPAIRS + 1)]
        self.exphairpin: List[float] = [0.0] * 31
        self.expbulge: List[float] = [0.0] * (MAXLOOP + 1)
        self.expinternal: List[float] = [0.0] * (MAXLOOP + 1)
        self.expmismatchExt: List[List[List[float]]] = [[[0.0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.expmismatchI: List[List[List[float]]] = [[[0.0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.expmismatch23I: List[List[List[float]]] = [[[0.0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.expmismatch1nI: List[List[List[float]]] = [[[0.0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.expmismatchH: List[List[List[float]]] = [[[0.0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.expmismatchM: List[List[List[float]]] = [[[0.0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.expdangle5: List[List[float]] = [[0.0] * 5 for _ in range(NBPAIRS + 1)]
        self.expdangle3: List[List[float]] = [[0.0] * 5 for _ in range(NBPAIRS + 1)]
        self.expint11: List[List[List[List[float]]]] = [[[[0.0] * 5 for _ in range(5)] for _ in range(NBPAIRS + 1)] for _ in range(NBPAIRS + 1)]
        self.expint21: List[List[List[List[List[float]]]]] = [[[[[0.0] * 5 for _ in range(5)] for _ in range(5)] for _ in range(NBPAIRS + 1)] for _ in range(NBPAIRS + 1)]
        self.expint22: List[List[List[List[List[float]]]]] = [[[[[0.0] * 5 for _ in range(5)] for _ in range(5)] for _ in range(5)] for _ in range(NBPAIRS + 1)]
        self.expninio: List[List[float]] = [[0.0] * (MAXLOOP + 1) for _ in range(5)]
        self.lxc: float = 0.0
        self.expMLbase: float = 0.0
        self.expMLintern: List[float] = [0.0] * (NBPAIRS + 1)
        self.expMLclosing: float = 0.0
        self.expTermAU: float = 0.0
        self.expDuplexInit: float = 0.0
        self.exptetra: List[float] = [0.0] * 40
        self.exptri: List[float] = [0.0] * 40
        self.exphex: List[float] = [0.0] * 40
        self.Tetraloops: str = ""  # Assuming it's a string
        self.expTriloop: List[float] = [0.0] * 40
        self.Triloops: str = ""  # Assuming it's a string
        self.Hexaloops: str = ""  # Assuming it's a string
        self.expTripleC: float = 0.0
        self.expMultipleCA: float = 0.0
        self.expMultipleCB: float = 0.0
        self.expgquad: List[List[float]] = [[0.0] * (3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1) for _ in range(VRNA_GQUAD_MAX_STACK_SIZE + 1)]
        self.expgquadLayerMismatch: float = 0.0
        self.gquadLayerMismatchMax: int = 0

        self.kT: float = 0.0
        self.pf_scale: float = 0.0  # Scaling factor to avoid over-/underflows

        self.temperature: float = 0.0  # Temperature used for loop contribution scaling
        self.alpha: float = 0.0  # Scaling factor for the thermodynamic temperature

        self.model_details: vrna_md_t = vrna_md_t()  # Model details to be used in the recursions
        self.param_file: str = ""  # The filename the parameters were derived from

        self.expSaltStack: float = 0.0
        self.expSaltLoop: List[float] = [0.0] * (MAXLOOP + 2)
        self.SaltLoopDbl: List[float] = [0.0] * (MAXLOOP + 2)
        self.SaltMLbase: int = 0
        self.SaltMLintern: int = 0
        self.SaltMLclosing: int = 0
        self.SaltDPXInit: int = 0




    

 
class vrna_ud_t:
    def __init__(self) -> None:
        self.exp_energy_cb = vrna_ud_exp_f
        self.uniq_motif_count = self.uniq_motif_size = self.motif_count = self.motif = self.motif_name = self.motif_size = self.motif_en = self.motif_type = 0
        self.probs_add = vrna_ud_add_probs_f
     
  
vrna_ud_add_probs_f = Callable[[int, int, int, float, None], None]
vrna_ud_exp_f = Callable[[int, int, int, None], float]

   



class vrna_fold_compound_t:
    type:vrna_fc_type_e.VRNA_FC_TYPE_COMPARATIVE.value
    length: int
    strand_number: Optional[List[int]] = None
    strand_order: Optional[List[int]] = None
    strand_order_uniq: Optional[List[int]] = None
    strand_start: Optional[List[int]] = None
    strand_end: Optional[List[int]] = None
    strands: int = 0
    cutpoint: int = 0
    nucleotides: Optional[List[str]] = vrna_seq_t
    # alignment: Optional[List[str]] = vrna_msa_t
    hc: Optional[str] = vrna_hc_t  # This should be updated to the correct type
    matrices: Optional[str] = vrna_mx_mfe_t  # This should be updated to the correct type
    exp_matrices: Optional[str] = vrna_mx_pf_t  # This should be updated to the correct type
    params: Optional[str] = vrna_param_t  # This should be updated to the correct type
    exp_params: Optional[str] = vrna_exp_param_t  # This should be updated to the correct type
    iindx: Optional[List[int]] = None
    jindx: Optional[List[int]] = None
    stat_cb: Optional[Callable] = vrna_recursion_status_f
    auxdata: Optional[object] = None
    free_auxdata: Optional[Callable] = vrna_auxdata_free_f
    domains_struc: Optional[str] = vrna_sd_t  # This should be updated to the correct type
    domains_up: Optional[str] = vrna_ud_t  # This should be updated to the correct type
    aux_grammar: Optional[str] = vrna_gr_aux_t  # This should be updated to the correct type

    # Fields for single/hybrid structure prediction
    sequence: Optional[str] = None
    sequence_encoding: Optional[List[int]] = None
    encoding5: Optional[List[int]] = None
    encoding3: Optional[List[int]] = None
    sequence_encoding2: Optional[List[int]] = None
    ptype: Optional[str] = None
    ptype_pf_compat: Optional[str] = None
    sc: Optional[str] = vrna_sc_t  # This should be updated to the correct type

    # Fields for consensus structure prediction
    sequences: Optional[List[str]] = None
    n_seq: int = 0
    cons_seq: Optional[str] = None
    S_cons: Optional[List[int]] = None
    S: Optional[List[List[int]]] = None
    S5: Optional[List[List[int]]] = None
    S3: Optional[List[List[int]]] = None
    Ss: Optional[List[str]] = None
    a2s: Optional[List[List[int]]] = None
    pscore: Optional[List[int]] = None
    pscore_local: Optional[List[List[int]]] = None
    pscore_pf_compat: Optional[List[int]] = None
    scs: Optional[List[str]] = vrna_sc_t  # This should be updated to the correct type
    oldAliEn: int = 0

    # Fields for Distance Class Partitioning
    maxD1: int = 0
    maxD2: int = 0
    reference_pt1: Optional[List[int]] = None
    reference_pt2: Optional[List[int]] = None
    referenceBPs1: Optional[List[int]] = None
    referenceBPs2: Optional[List[int]] = None
    bpdist: Optional[List[int]] = None
    mm1: Optional[List[int]] = None
    mm2: Optional[List[int]] = None

    # Fields for local folding
    window_size: int = 0
    ptype_local: Optional[List[str]] = None
    zscore_data: Optional[str] = None  # This should be updated to the correct type


  

class vrna_fold_compound_t_:
    def __init__(self) -> None:
        # Common data fields
        self.type = vrna_fc_type_e.VRNA_FC_TYPE_COMPARATIVE.value # vrna_fc_type_e
        self.length = 0  # unsigned int
        self.strand_number = []  # unsigned int array
        self.strand_order = []  # unsigned int array
        self.strand_order_uniq = []  # unsigned int array
        self.strand_start = []  # unsigned int array
        self.strand_end = []  # unsigned int array
        self.strands = 0  # unsigned int
        # self.nucleotides = vrna_seq_t  # vrna_seq_t array
        # self.alignment = vrna_msa_t  # vrna_msa_t array
        self.hc = vrna_hc_t()  # vrna_hc_t
        self.matrices = vrna_mx_mfe_t()  # vrna_mx_mfe_t
        self.exp_matrices = vrna_mx_pf_t()  # vrna_mx_pf_t
        # self.params = vrna_param_t  # vrna_param_t
        self.exp_params = vrna_exp_param_t()  # vrna_exp_param_t
        self.iindx = 0  # int array
        self.jindx = 0  # int array

        # User-defined data fields
        # self.stat_cb = vrna_recursion_status_f  # vrna_recursion_status_f
        self.auxdata = 0  # void pointer
        # self.free_auxdata = vrna_auxdata_free_f  # vrna_auxdata_free_f

        # Secondary Structure Decomposition (grammar) related data fields
        # self.domains_struc = vrna_sd_t  # vrna_sd_t
        self.domains_up = vrna_ud_t()  # vrna_ud_t
        # self.aux_grammar = vrna_gr_aux_t  # vrna_gr_aux_t

        # Data fields available for single/hybrid structure prediction
        self.sequence = ''  # char array
        self.sequence_encoding = 0  # short array
        self.encoding5 = 0  # short array
        self.encoding3 = 0  # short array
        self.sequence_encoding2 = 0  # short array
        self.ptype = 'ptype'  # char array
        self.ptype_pf_compat = 'ptype_pf_compat'  # char array
        self.sc = vrna_sc_t()  # vrna_sc_t

        # Data fields for consensus structure prediction
        self.n_seq = 0  # unsigned int
        self.cons_seq = 'cons_seq'  # char array
        self.S_cons = 0  # short array
        self.S = 0  # short array of arrays
        self.S5 = 0  # short array of arrays
        self.S3 = 0  # short array of arrays
        self.Ss = 'Ss'  # char array of arrays
        self.a2s = 0  # unsigned int array of arrays
        self.pscore = list[int]  # int array
        self.pscore_local = 0  # int array of arrays
        self.pscore_pf_compat = 0  # short array
        self.scs = vrna_sc_t()  # vrna_sc_t array of arrays
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
        # self.zscore_data = vrna_zsc_dat_t  # vrna_zsc_dat_t

# data structure/class end ######################################

# function/method start ######################################





    


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


class constraints_helper:
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


class hc_ext_def_dat:
    def __init__(self) -> None:
        self.n = self.mx = self.mx_window =self.sn = self.hc_up = self.hc_dat = None
        self.hc_f = vrna_hc_eval_f()


def hc_ext_cb_def(i, j, k, l, d, data:hc_ext_def_dat):
    dat = data  # Assume data is an instance of HC_EXT_Def_Dat
    
    eval = 0
    di = k - i
    dj = j - l
    n = dat.n
    # Handle different cases based on d
    if d == VRNA_DECOMP_EXT_EXT_STEM:
        if dat.mx[dat.n * j + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP:
            eval = 1
            if i != l:
                di = l - k - 1
                if di != 0 and dat.hc_up[k + 1] < di:
                    eval = 0

    elif d == VRNA_DECOMP_EXT_STEM_EXT:
        if dat.mx[dat.n * k + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP:
            eval = 1
            if i != l:
                di = l - k - 1
                if di != 0 and dat.hc_up[k + 1] < di:
                    eval = 0

    elif d == VRNA_DECOMP_EXT_EXT_STEM1:
        if dat.mx[dat.n * (j - 1) + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP:
            eval = 1
            if dat.hc_up[j] == 0:
                eval = 0
            if i != l:
                di = l - k - 1
                if di != 0 and dat.hc_up[k + 1] < di:
                    eval = 0

    elif d == VRNA_DECOMP_EXT_STEM_EXT1:
        if dat.mx[dat.n * k + i + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP:
            eval = 1
            if dat.hc_up[i] == 0:
                eval = 0
            if j != k:
                dj = l - k - 1
                if dj != 0 and dat.hc_up[k + 1] < dj:
                    eval = 0

    elif d == VRNA_DECOMP_EXT_EXT_EXT:
        eval = 1
        di = l - k - 1
        if di != 0 and dat.hc_up[k + 1] < di:
            eval = 0

    elif d == VRNA_DECOMP_EXT_STEM:
        if dat.mx[dat.n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP:
            eval = 1
            if di != 0 and dat.hc_up[i] < di:
                eval = 0
            if dj != 0 and dat.hc_up[l + 1] < dj:
                eval = 0

    elif d == VRNA_DECOMP_EXT_EXT:
        eval = 1
        if di != 0 and dat.hc_up[i] < di:
            eval = 0
        if dj != 0 and dat.hc_up[l + 1] < dj:
            eval = 0

    elif d == VRNA_DECOMP_EXT_UP:
        di = j - i + 1
        eval = 1 if dat.hc_up[i] >= di else 0

    elif d == VRNA_DECOMP_EXT_STEM_OUTSIDE:
        if dat.mx[dat.n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP:
            eval = 1

    else:
        print(f"Warning: Unrecognized decomposition {d}")

    return eval


def hc_ext_cb_sn(i, j, k, l, d, data:hc_ext_def_dat):
    dat = data  # Assume data is an instance of HC_EXT_Def_Dat
    
    sn = dat.sn
    eval = 0

    # Handle different cases based on d
    if d == VRNA_DECOMP_EXT_EXT_STEM1:
        if sn[j - 1] != sn[j]:
            eval = 0
        elif sn[k] == sn[l]:
            eval = 1

    elif d == VRNA_DECOMP_EXT_STEM_EXT1:
        if sn[i] != sn[i + 1]:
            eval = 0
        elif sn[k] == sn[l]:
            eval = 1

    elif d == VRNA_DECOMP_EXT_EXT_STEM or d == VRNA_DECOMP_EXT_STEM_EXT or d == VRNA_DECOMP_EXT_EXT_EXT:
        if sn[k] == sn[l]:
            eval = 1

    elif d == VRNA_DECOMP_EXT_STEM or d == VRNA_DECOMP_EXT_EXT:
        if sn[i] == sn[k] and sn[l] == sn[j]:
            eval = 1

    elif d == VRNA_DECOMP_EXT_UP:
        if sn[i] == sn[j]:
            eval = 1

    elif d == VRNA_DECOMP_EXT_STEM_OUTSIDE:
        if (k <= i or sn[k - 1] == sn[k]) and (l >= j or sn[l + 1] == sn[l]):
            eval = 1

    else:
        print(f"Warning: Unrecognized decomposition {d}")

    return eval



def hc_ext_cb_def_sn(i, j, k, l, d, data:hc_ext_def_dat):
    eval = hc_ext_cb_def(i, j, k, l, d, data) if hc_ext_cb_sn(i, j, k, l, d, data) else 0
    return eval

def hc_ext_cb_def_user(i, j, k, l, d, data:hc_ext_def_dat):
    dat = data
    eval = hc_ext_cb_def(i, j, k, l, d, data) if dat.hc_f(i, j, k, l, d, dat.hc_dat) else 0
    return eval


def hc_ext_cb_def_sn_user(i, j, k, l, d, data:hc_ext_def_dat):
    dat = data
    eval = hc_ext_cb_def(i, j, k, l, d, data) if hc_ext_cb_sn(i, j, k, l, d, data) else 0
    eval = hc_ext_cb_def(i, j, k, l, d, data) if dat.hc_f(i, j, k, l, d, dat.hc_dat) else 0
    return eval

def prepare_hc_ext_def(fc: vrna_fold_compound_t, dat: hc_ext_def_dat) -> Callable:
    dat.mx = fc.hc.mx
    dat.n = fc.length
    dat.hc_up = fc.hc.up_ext
    dat.sn = fc.strand_number

    if fc.hc.f:
        dat.hc_f = fc.hc.f
        dat.hc_dat = fc.hc.data
        return hc_ext_cb_def_user if fc.strands == 1 else hc_ext_cb_def_sn_user

    return hc_ext_cb_def if fc.strands == 1 else hc_ext_cb_def_sn

class hc_hp_def_dat:
    def __init__(self) -> None:
        self.n = self.mx = self.mx_window = self.sn = self.hc_up = self.hc_dat = None
        self.hc_f = vrna_hc_eval_f()


def hc_hp_cb_def(i, j, k, l, d, data:hc_hp_def_dat):
    dat = data  # Assume data is an instance of HC_HP_Def_Dat
    
    eval = 0

    # No strand nicks are allowed in hairpin loops
    if dat.sn[i] != dat.sn[j]:
        return eval

    if j > i:
        # Linear case
        p = i
        q = j
        u = q - p - 1
    else:
        # Circular case
        p = j
        q = i
        u = dat.n - q + p - 1

    if dat.mx[dat.n * p + q] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP:
        eval = 1
        if dat.hc_up[i + 1] < u:
            eval = 0

    return eval

def hc_hp_cb_def_user(i, j, k, l, d, data:hc_hp_def_dat):
    dat = data
    eval = hc_hp_cb_def(i, j, k, l, d, data) if (dat.hc_f(i, j, k, l, d, dat.hc_dat)) else 0
    return eval


def prepare_hc_hp_def(fc: vrna_fold_compound_t, dat: hc_hp_def_dat) -> Callable:
    dat.mx = fc.hc.mx
    dat.hc_up = fc.hc.up_hp
    dat.n = fc.length
    dat.sn = fc.strand_number

    if fc.hc.f:
        dat.hc_f = fc.hc.f
        dat.hc_dat = fc.hc.data
        return hc_hp_cb_def_user

    return hc_hp_cb_def


# hc_int_def_dat

class hc_int_def_dat:
    def __init__(self):
        self.mx = None
        self.mx_local = None
        self.n = 0
        self.up = 0
        self.sn = 0
        self.hc_f = None
        self.hc_dat = None
        

def hc_int_cb_def(i, j, k, l, data:hc_int_def_dat):
    dat = data  # Assume data is an instance of HC_INT_Def_Dat
    
    pij = 0
    pkl = 0

    # Check strand conditions
    if dat.sn[i] != dat.sn[k] or dat.sn[l] != dat.sn[j]:
        return 0

    # Retrieve values based on mx or mx_local
    if dat.mx:
        pij = dat.mx[dat.n * i + j]
        pkl = dat.mx[dat.n * k + l]
    else:
        pij = dat.mx_local[i][j - i]
        pkl = dat.mx_local[k][l - k]

    # Check conditions for internal loops
    if (pij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) and (pkl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC):
        return 1

    return 0


def hc_int_cb_def_user(i, j, k, l, data:hc_int_def_dat):
    dat = data
    eval = hc_int_cb_def(i, j, k, l, data) if (dat.hc_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, dat.hc_dat)) else 0
    return eval

def prepare_hc_int_def(fc:vrna_fold_compound_t, dat:hc_int_def_dat) -> Callable:
    if fc.hc.type == VRNA_HC_WINDOW:
        dat.mx = None
        dat.mx_local = fc.hc.matrix_local
    else:
        dat.mx = fc.hc.mx
        dat.mx_local = None
    
    dat.n = fc.length
    dat.up = fc.hc.up_int
    dat.sn = fc.strand_number
    
    if fc.hc.f:
        dat.hc_f = fc.hc.f
        dat.hc_dat = fc.hc.data
        return hc_int_cb_def_user
    
    return hc_int_cb_def

class vrna_hc_eval_f:
    def __init__(self) -> None:
        self.i = self.j = self.k = self.l = self.d = self.data = None

class hc_mb_def_dat:
    def __init__(self) -> None:
        self.mx = None
        self.n = None
        self.mx_window = None
        self.hc_up = None
        self.sn = None
        self.hc_f = vrna_hc_eval_f()
        self.hc_dat = None


def hc_mb_cb_def_window(i, j, k, l, d, data:hc_mb_def_dat):
    dat = data  # Assume data is an instance of HC_MB_Def_Dat
    
    eval = 0
    di = k - i
    dj = j - l

    # Switch based on decomposition type d
    if d == VRNA_DECOMP_ML_ML_ML:
        u = l - k - 1
        eval = 1
        if u != 0 and dat.hc_up[k + 1] < u:
            eval = 0
        if dat.sn[k] != dat.sn[l]:
            eval = 0

    elif d == VRNA_DECOMP_ML_ML:
        eval = 1
        if di != 0 and (dat.hc_up[i] < di or dat.sn[i] != dat.sn[k]):
            eval = 0
        if dj != 0 and (dat.hc_up[l + 1] < dj or dat.sn[l] != dat.sn[j]):
            eval = 0

    elif d == VRNA_DECOMP_ML_STEM:
        if dat.mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC:
            eval = 1
            if di != 0 and dat.hc_up[i] < di:
                eval = 0
            if dj != 0 and dat.hc_up[l + 1] < dj:
                eval = 0

    elif d == VRNA_DECOMP_PAIR_ML:
        if dat.mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP:
            eval = 1
            di -= 1
            dj -= 1
            if di != 0 and dat.hc_up[i + 1] < di:
                eval = 0
            if dj != 0 and dat.hc_up[l + 1] < dj:
                eval = 0

    elif d == VRNA_DECOMP_ML_COAXIAL:
        if dat.mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC:
            eval = 1

    elif d == VRNA_DECOMP_ML_COAXIAL_ENC:
        if (dat.mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC and
                dat.mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC):
            eval = 1

    else:
        print(f"Warning: Unrecognized decomposition {d}")

    return eval




def hc_mb_cb_def_user_window(i, j, k, l, d, data:hc_mb_def_dat):
    dat = data  # Assume data is an instance of HC_MB_Def_Dat
    
    # Evaluate using hc_mb_cb_def_window
    eval = hc_mb_cb_def_window(i, j, k, l, d, data)

    # Apply additional user-defined function logic
    if dat.hc_f(i, j, k, l, d, dat.hc_dat):
        return eval
    else:
        return 0


def hc_mb_cb_def(i, j, k, l, d, data:hc_mb_def_dat):
    eval = 0
    di = k - i
    dj = j - l
    n = data.n

    # Switch case translation
    if d == VRNA_DECOMP_ML_ML_ML:
        u = l - k - 1
        eval = 1
        if u != 0 and data.hc_up[k + 1] < u:
            eval = 0

    elif d == VRNA_DECOMP_ML_ML:
        eval = 1
        if di != 0 and data.hc_up[i] < di:
            eval = 0
        if dj != 0 and data.hc_up[l + 1] < dj:
            eval = 0

    elif d == VRNA_DECOMP_ML_STEM:
        if data.mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC:
            eval = 1
            if di != 0 and data.hc_up[i] < di:
                eval = 0
            if dj != 0 and data.hc_up[l + 1] < dj:
                eval = 0

    elif d == VRNA_DECOMP_PAIR_ML:
        if data.mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP:
            eval = 1
            di -= 1
            dj -= 1
            if di != 0 and data.hc_up[i + 1] < di:
                eval = 0
            if dj != 0 and data.hc_up[l + 1] < dj:
                eval = 0

    elif d == VRNA_DECOMP_PAIR_ML_EXT:
        if data.mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP:
            eval = 1
            di += 1
            dj += 1
            if di != 0 and data.hc_up[k + 1] < di:
                eval = 0
            if dj != 0 and data.hc_up[j + 1] < dj:
                eval = 0

    elif d == VRNA_DECOMP_ML_ML_STEM:
        u = l - k - 1
        if data.mx[n * j + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC:
            eval = 1
        if u != 0 and data.hc_up[k + 1] < u:
            eval = 0

    elif d == VRNA_DECOMP_ML_COAXIAL:
        if data.mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC:
            eval = 1

    elif d == VRNA_DECOMP_ML_COAXIAL_ENC:
        if (data.mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC and
                data.mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC):
            eval = 1

    else:
        print("hc_mb_cb_def@multibranch_hc.inc: Unrecognized decomposition %d" % d)

    return eval


def hc_mb_cb_def_sn(i, j, k, l, d, data:hc_mb_def_dat):
    eval = hc_mb_cb_def(i, j, k, l, d, data) if hc_sn(i, j, k, l, d, data) else 0
    return eval


def hc_sn(i, j, k, l, d, data:hc_mb_def_dat):
    sn = data.sn
    eval = 0

    # Switch case translation
    if d == VRNA_DECOMP_ML_ML_ML or d == VRNA_DECOMP_ML_ML_STEM:
        if sn[k] == sn[l]:
            eval = 1

    elif d == VRNA_DECOMP_ML_STEM or d == VRNA_DECOMP_ML_ML:
        if (sn[i] == sn[k] and
            sn[l] == sn[j] and
            sn[i - 1] == sn[i] and
            sn[j + 1] == sn[j]):
            eval = 1

    elif d == VRNA_DECOMP_PAIR_ML_EXT or d == VRNA_DECOMP_PAIR_ML:
        if sn[i] == sn[k] and sn[l] == sn[j]:
            eval = 1

    elif d == VRNA_DECOMP_ML_COAXIAL:
        if (i == k - 1 and sn[i] == sn[k]) or (l + 1 == j and sn[l] == sn[j]):
            eval = 1

    elif d == VRNA_DECOMP_ML_COAXIAL_ENC:
        if sn[j] == sn[k]:
            eval = 1

    else:
        print("hc_sn@multibranch_hc.inc: Unrecognized decomposition %d" % d)

    return eval




def hc_mb_cb_def_user(i, j, k, l, d, data:hc_mb_def_dat):
    dat = data  # Assume data is an instance of HC_MB_Def_Dat
    
    # Evaluate using hc_mb_cb_def
    eval = hc_mb_cb_def(i, j, k, l, d, data)

    # Apply additional user-defined function logic
    if dat.hc_f(i, j, k, l, d, dat.hc_dat):
        return eval
    else:
        return 0

def hc_mb_cb_def_sn_user(i, j, k, l, d, data:hc_mb_def_dat):
    dat = data  # Assume data is an instance of HC_MB_Def_Dat
    
    # Evaluate using hc_mb_cb_def
    eval = hc_mb_cb_def(i, j, k, l, d, data)

    # Apply strand number function logic
    if hc_sn(i, j, k, l, d, data):
        return eval
    else:
        return 0


def prepare_hc_mb_def(fc:vrna_fold_compound_t, dat:hc_mb_def_dat):
    dat.mx         = fc.hc.mx
    dat.n          = fc.hc.n
    dat.mx_window  = fc.hc.matrix_local
    dat.hc_up      = fc.hc.up_ml
    dat.sn         = fc.strand_number
    
    if fc.hc.f:
        dat.hc_f = fc.hc.f
        dat.hc_dat = fc.hc.data
        if fc.hc.type == VRNA_HC_WINDOW:
            return hc_mb_cb_def_user_window
        elif fc.strands == 1:
            return hc_mb_cb_def_user
        else:
            return hc_mb_cb_def_sn_user
    
    if fc.hc.type == VRNA_HC_WINDOW:
        return hc_mb_cb_def_window
    elif fc.strands == 1:
        return hc_mb_cb_def
    else:
        return hc_mb_cb_def_sn


    


class sc_ext_exp_dat:
    def __init__(self):
        self.up = None  # FLT_OR_DBL **up
        self.red_ext = sc_ext_exp_cb  # sc_ext_exp_cb red_ext
        self.red_stem = sc_ext_exp_cb  # sc_ext_exp_cb red_stem
        self.red_up = sc_ext_exp_red_up  # sc_ext_exp_red_up red_up
        self.split = sc_ext_exp_split  # sc_ext_exp_split split

        self.user_cb = vrna_sc_exp_f  # vrna_sc_exp_f user_cb
        self.user_data = None  # void *user_data

        self.n_seq = 0  # int n_seq
        self.a2s = None  # unsigned int **a2s
        self.up_comparative = None  # FLT_OR_DBL ***up_comparative

        self.user_cb_comparative = vrna_sc_exp_f  # vrna_sc_exp_f *user_cb_comparative
        self.user_data_comparative = None  # void **user_data_comparative
      
sc_ext_exp_cb = Callable[[int, int, int, int, sc_ext_exp_dat], float]
sc_ext_exp_red_up = Callable[[int, int, sc_ext_exp_dat], float]
sc_ext_exp_split = Callable[[int, int, int, sc_ext_exp_dat], float]
vrna_sc_exp_f = Callable[[int, int, int, int, int, None], float]


def sc_ext_exp_cb_red(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:

  sc_up = data.up

  q_sc = 1.

  length_1  = k - i
  start_2   = l + 1
  length_2  = j - l

  if (length_1 != 0):
    q_sc *= sc_up[i][length_1]

  if (length_2 != 0):
    q_sc *= sc_up[start_2][length_2]

  return q_sc


def sc_ext_exp_cb_red_user_to_ext(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    return data.user_cb(i, j, k, l, VRNA_DECOMP_EXT_EXT, data.user_data)


def sc_ext_exp_cb_red_user_def_to_ext(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    return sc_ext_exp_cb_red(i, j, k, l, data) * sc_ext_exp_cb_red_user_to_ext(i, j, k, l, data)


def sc_ext_exp_cb_red_user_to_stem(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    return data.user_cb(i, j, k, l, VRNA_DECOMP_EXT_STEM, data.user_data)


def sc_ext_exp_cb_red_user_def_to_stem(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    return sc_ext_exp_cb_red(i, j, k, l, data) * sc_ext_exp_cb_red_user_to_stem(i, j, k, l, data)

def sc_ext_exp_cb_up(i:int, j:int, data:sc_ext_exp_dat) -> float:
    sc_up = data.up
    length = j - i + 1
    q_sc = 1.
    if length != 0:
        q_sc *= sc_up[i][length]
    return q_sc

def sc_ext_exp_cb_up_user(i:int, j:int, data:sc_ext_exp_dat) -> float:
    return data.user_cb(i, j, i, j, VRNA_DECOMP_EXT_UP, data.user_data)


def sc_ext_exp_cb_up_user_def(i:int, j:int, data:sc_ext_exp_dat) -> float:
    return sc_ext_exp_cb_up(i, j, data) * sc_ext_exp_cb_up_user(i, j, data)


def sc_ext_exp_cb_split_user(i:int, j:int, k:int, data:sc_ext_exp_dat) -> float:
    return data.user_cb(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, data.user_data)


def sc_ext_exp_cb_red_comparative(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    sc_up = data.up_comparative
    a2s   = data.a2s

    q_sc = 1.
    for s in range(data.n_seq): 
        if (sc_up[s]): 
            length_1  = a2s[s][k] - a2s[s][i]
            start_2   = a2s[s][l] + 1
            length_2  = a2s[s][j] - a2s[s][l]

        if (length_1 != 0):
            q_sc *= sc_up[s][a2s[s][i]][length_1]

        if (length_2 != 0):
            q_sc *= sc_up[s][start_2][length_2]   

    return q_sc


def sc_ext_exp_cb_red_user_to_ext_comparative(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    q_sc = 1.
    for s in range(data.n_seq):
        q_sc *= data.user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_EXT_EXT, data.user_data_comparative[s])
    return q_sc



def sc_ext_exp_cb_red_user_def_to_ext_comparative(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    return sc_ext_exp_cb_red_comparative(i, j, k, l, data) * sc_ext_exp_cb_red_user_to_ext_comparative(i, j, k, l, data)


def sc_ext_exp_cb_red_user_to_stem_comparative(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    q_sc = 1.
    for s in range(data.n_seq):
        q_sc *= data.user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_EXT_STEM, data.user_data_comparative[s])
    return q_sc

def sc_ext_exp_cb_red_user_def_to_stem_comparative(i:int, j:int, k:int, l:int, data:sc_ext_exp_dat) -> float:
    return sc_ext_exp_cb_red_comparative(i, j, k, l, data) * sc_ext_exp_cb_red_user_to_stem_comparative(i, j, k, l, data)

def sc_ext_exp_cb_up_comparative(i:int, j:int, data:sc_ext_exp_dat) -> float:
    a2s = data.a2s
    sc_up = data.up_comparative
    q_sc = 1.
    length = 0
    for s in range(data.n_seq):
        length = a2s[s][j - 1] - a2s[s][i]
        
        if length != 0:
            q_sc *= sc_up[s][a2s[s][i]][length]
    
    return q_sc


   
def sc_ext_exp_cb_up_user_comparative(i:int, j:int, data:sc_ext_exp_dat) -> float:
    q_sc = 1.
    for s in range(data.n_seq):
        q_sc *= data.user_cb_comparative[s](i, j, i, j, VRNA_DECOMP_EXT_UP, data.user_data_comparative[s])

    return q_sc
 

def sc_ext_exp_cb_up_user_def_comparative(i:int, j:int, data:sc_ext_exp_dat) -> float:
    return sc_ext_exp_cb_up_comparative(i, j, data) * sc_ext_exp_cb_up_user_comparative(i, j, data)


def sc_ext_exp_cb_split_user_comparative(i:int, j:int, k:int, data:sc_ext_exp_dat) -> float:
    q_sc = 1.
    for s in range(data.n_seq):
        q_sc *= data.user_cb_comparative[s](i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT,data.user_data_comparative[s])
    return q_sc



def init_sc_ext_exp(fc:vrna_fold_compound_t, sc_wrapper:sc_ext_exp_dat):
    sc_wrapper.up = None
    sc_wrapper.user_cb = None
    sc_wrapper.user_data = None
    sc_wrapper.n_seq = 1
    sc_wrapper.a2s = None
    sc_wrapper.up_comparative = None
    sc_wrapper.user_cb_comparative = None
    sc_wrapper.user_data_comparative = None
    sc_wrapper.red_ext = None
    sc_wrapper.red_stem = None
    sc_wrapper.red_up = None
    sc_wrapper.split = None

    # no soft constraints by default
    if fc.type == VRNA_FC_TYPE_SINGLE:
        sc = fc.sc
        if sc:
            sc_wrapper.up = sc.exp_energy_up
            sc_wrapper.user_cb = sc.exp_f
            sc_wrapper.user_data = sc.data
            if sc.exp_energy_up:
                if sc.exp_f:
                    sc_wrapper.red_ext = sc_ext_exp_cb_red_user_def_to_ext
                    sc_wrapper.red_stem = sc_ext_exp_cb_red_user_def_to_stem
                    sc_wrapper.red_up = sc_ext_exp_cb_up_user_def
                    sc_wrapper.split = sc_ext_exp_cb_split_user
                else:
                    sc_wrapper.red_ext = sc_ext_exp_cb_red
                    sc_wrapper.red_stem = sc_ext_exp_cb_red
                    sc_wrapper.red_up = sc_ext_exp_cb_up
            elif sc.exp_f:
                sc_wrapper.red_ext = sc_ext_exp_cb_red_user_to_ext
                sc_wrapper.red_stem = sc_ext_exp_cb_red_user_to_stem
                sc_wrapper.red_up = sc_ext_exp_cb_up_user

    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        scs = fc.scs
        sc_wrapper.n_seq = fc.n_seq
        sc_wrapper.a2s = fc.a2s
        if scs:
            sc_wrapper.up_comparative = []
            sc_wrapper.user_cb_comparative = []
            sc_wrapper.user_data_comparative = []
            provides_sc_up = False
            provides_sc_user_cb = False
            for s in range(fc.n_seq):
                if scs[s]:
                    sc_wrapper.up_comparative.append(scs[s].exp_energy_up)
                    sc_wrapper.user_cb_comparative.append(scs[s].exp_f)
                    sc_wrapper.user_data_comparative.append(scs[s].data)
                    if scs[s].exp_energy_up:
                        provides_sc_up = True
                    if scs[s].exp_f:
                        provides_sc_user_cb = True
                    if provides_sc_up:
                        if provides_sc_user_cb:
                            sc_wrapper.red_ext = sc_ext_exp_cb_red_user_def_to_ext_comparative
                            sc_wrapper.red_stem = sc_ext_exp_cb_red_user_def_to_stem_comparative
                            sc_wrapper.red_up = sc_ext_exp_cb_up_user_def_comparative
                            sc_wrapper.split = sc_ext_exp_cb_split_user_comparative
                        else:
                            sc_wrapper.red_ext = sc_ext_exp_cb_red_comparative
                            sc_wrapper.red_stem = sc_ext_exp_cb_red_comparative
                            sc_wrapper.red_up = sc_ext_exp_cb_up_comparative
                    elif provides_sc_user_cb:
                        sc_wrapper.red_ext = sc_ext_exp_cb_red_user_to_ext_comparative
                        sc_wrapper.red_stem = sc_ext_exp_cb_red_user_to_stem_comparative
                        sc_wrapper.red_up = sc_ext_exp_cb_up_user_comparative
    return

# constant define for init_sc_hp_exp ###########
VRNA_MX_WINDOW = 1
VRNA_SC_WINDOW = 1

# class define for init_sc_hp_exp ###########
class sc_hp_exp_dat(sc_ext_exp_dat):
    def __init__(self):
        super().__init__()
        self.n = self.idx = 0
        self.bp = self.bp_comparative = self.bp_local = self.bp_local_comparative = 0.
        self.pair = self.pair_ext = sc_hp_exp_cb

sc_hp_exp_cb = Callable[[int, int, sc_hp_exp_dat], float]

# function define for init_sc_hp_exp ###########
def sc_hp_exp_cb_ext_up(i:int, j:int, data:sc_hp_exp_dat) -> float:
    u1 = data.n - j
    u2 = i - 1
    sc = 1.0

    if u1 > 0:
        sc *= data.up[j + 1][u1]

    if u2 > 0:
        sc *= data.up[1][u2]

    return sc

def sc_hp_exp_cb_ext_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return data.user_cb(j, i, j, i,VRNA_DECOMP_PAIR_HP,data.user_data)

def sc_hp_exp_cb_ext_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    sc = 1.
    for s in range(data.n_seq):
        if (data.user_cb_comparative[s]):
            sc *= data.user_cb_comparative[s](j, i, j, i,
                                         VRNA_DECOMP_PAIR_HP,
                                         data.user_data_comparative[s])

    return sc

def sc_hp_exp_cb_ext_up_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_ext_up(i, j, data) * sc_hp_exp_cb_ext_user(i, j, data)


def sc_hp_exp_cb_up(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return data.up[i + 1][j - i -1]

def sc_hp_exp_cb_bp(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return data.bp[data.idx[j] + i]

def sc_hp_exp_cb_bp_local(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return data.bp_local[i][j - i]

def sc_hp_exp_cb_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return data.user_cb(i, j, i, j, VRNA_DECOMP_PAIR_HP, data.user_data)

def sc_hp_exp_cb_up_bp_local_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up(i, j, data) * sc_hp_exp_cb_bp_local(i, j, data) * sc_hp_exp_cb_user(i, j, data)

def sc_hp_exp_cb_up_bp_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up(i, j, data) * sc_hp_exp_cb_bp(i, j, data) * sc_hp_exp_cb_user(i, j, data)


def sc_hp_exp_cb_up_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up(i, j, data) * sc_hp_exp_cb_user(i, j, data)

def sc_hp_exp_cb_bp_local_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_bp_local(i, j, data) * sc_hp_exp_cb_user(i, j, data)

def sc_hp_exp_cb_bp_user(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_bp(i, j, data) * sc_hp_exp_cb_user(i, j, data)



def sc_hp_exp_cb_up_bp_local(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up(i, j, data) * sc_hp_exp_cb_bp_local(i, j, data)

def sc_hp_exp_cb_up_bp(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up(i, j, data) * sc_hp_exp_cb_bp(i, j, data)



def sc_hp_exp_cb_ext_up_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.up_comparative[s]:
            u1 = data.a2s[s][data.n] - data.a2s[s][j]
            u2 = data.a2s[s][i - 1]

            if u1 > 0:
                sc *= data.up[data.a2s[s][j + 1]][u1]

            if u2 > 0:
                sc *= data.up[1][u2]

    return sc


def sc_hp_exp_cb_ext_up_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_ext_up_comparative(i, j, data) * sc_hp_exp_cb_ext_user_comparative(i, j, data)


def sc_hp_exp_cb_up_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.up_comparative[s]:
            u = data.a2s[s][j - 1] - data.a2s[s][i]
            sc *= data.up_comparative[s][data.a2s[s][i + 1]][u]

    return sc


def sc_hp_exp_cb_bp_local_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.bp_local_comparative[s]:
            sc *= data.bp_local_comparative[s][i][j - i]

    return sc

def sc_hp_exp_cb_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, i, j,
                                              VRNA_DECOMP_PAIR_HP,
                                              data.user_data_comparative[s])

    return sc


def sc_hp_exp_cb_up_bp_local_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up_comparative(i, j, data) * sc_hp_exp_cb_bp_local_comparative(i, j, data) * sc_hp_exp_cb_user_comparative(i, j, data)


def sc_hp_exp_cb_bp_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    sc = 1.
    for s in range(data.n_seq):
        if data.bp_comparative[s]:
            sc *= data.bp_comparative[s][data.idx[j] + i]
    return sc

def sc_hp_exp_cb_up_bp_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up_comparative(i, j, data) * sc_hp_exp_cb_bp_comparative(i, j, data) * sc_hp_exp_cb_user_comparative(i, j, data)



def sc_hp_exp_cb_up_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up_comparative(i, j, data) * sc_hp_exp_cb_user_comparative(i, j, data)

def sc_hp_exp_cb_bp_local_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_bp_local_comparative(i, j, data) * sc_hp_exp_cb_user_comparative(i, j, data)

def sc_hp_exp_cb_bp_user_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_bp_comparative(i, j, data) * sc_hp_exp_cb_user_comparative(i, j, data)

def sc_hp_exp_cb_up_bp_local_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up_comparative(i, j, data) * sc_hp_exp_cb_bp_local_comparative(i, j, data)

def sc_hp_exp_cb_up_bp_comparative(i:int, j:int, data:sc_hp_exp_dat) -> float:
    return sc_hp_exp_cb_up_comparative(i, j, data) * sc_hp_exp_cb_bp_comparative(i, j, data)





# init_sc_hp_exp start #################
def init_sc_hp_exp(fc:vrna_fold_compound_t, sc_wrapper:sc_hp_exp_dat):
    sliding_window = 0

    if fc.exp_matrices:
        sliding_window = 1 if fc.exp_matrices.type == VRNA_MX_WINDOW else 0
    elif (fc.type == VRNA_FC_TYPE_SINGLE) and (fc.sc):
        sliding_window = 1 if fc.sc.type == VRNA_SC_WINDOW else 0
    elif fc.hc:
        sliding_window = 1 if fc.hc.type == VRNA_HC_WINDOW else 0
    else:
        sliding_window = 0

    sc_wrapper.n = int(fc.length)
    sc_wrapper.idx = fc.jindx
    sc_wrapper.n_seq = 1
    sc_wrapper.a2s = None

    sc_wrapper.up = None
    sc_wrapper.up_comparative = None
    sc_wrapper.bp = None
    sc_wrapper.bp_comparative = None
    sc_wrapper.bp_local = None
    sc_wrapper.bp_local_comparative = None

    sc_wrapper.user_cb = None
    sc_wrapper.user_data = None
    sc_wrapper.user_cb_comparative = None
    sc_wrapper.user_data_comparative = None

    sc_wrapper.pair = None
    sc_wrapper.pair_ext = None

    if fc.type == VRNA_FC_TYPE_SINGLE:
        sc = fc.sc

        if sc:
            provides_sc_up = 0
            provides_sc_bp = 0
            provides_sc_user = 0

            sc_wrapper.up = sc.exp_energy_up
            sc_wrapper.bp = None if sliding_window else sc.exp_energy_bp
            sc_wrapper.bp_local = sc.exp_energy_bp_local if sliding_window else None
            sc_wrapper.user_cb = sc.exp_f
            sc_wrapper.user_data = sc.data

            if sc.exp_energy_up:
                provides_sc_up = 1

            if sliding_window:
                if sc.exp_energy_bp_local:
                    provides_sc_bp = 1
            elif sc.exp_energy_bp:
                provides_sc_bp = 1

            if sc.exp_f:
                provides_sc_user = 1

            if provides_sc_user:
                sc_wrapper.pair_ext = sc_hp_exp_cb_ext_user
                if provides_sc_up:
                    sc_wrapper.pair_ext = sc_hp_exp_cb_ext_up_user

                    if provides_sc_bp:
                        sc_wrapper.pair = sc_hp_exp_cb_up_bp_local_user if sliding_window else sc_hp_exp_cb_up_bp_user
                    else:
                        sc_wrapper.pair = sc_hp_exp_cb_up_user
                elif provides_sc_bp:
                    sc_wrapper.pair = sc_hp_exp_cb_bp_local_user if sliding_window else sc_hp_exp_cb_bp_user
                else:
                    sc_wrapper.pair = sc_hp_exp_cb_user
            else:
                if provides_sc_up:
                    sc_wrapper.pair_ext = sc_hp_exp_cb_ext_up
                    if provides_sc_bp:
                        sc_wrapper.pair = sc_hp_exp_cb_up_bp_local if sliding_window else sc_hp_exp_cb_up_bp
                    else:
                        sc_wrapper.pair = sc_hp_exp_cb_up
                elif provides_sc_bp:
                    sc_wrapper.pair = sc_hp_exp_cb_bp_local if sliding_window else sc_hp_exp_cb_bp

    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        sc_wrapper.n_seq = fc.n_seq
        sc_wrapper.a2s = fc.a2s

        scs = fc.scs

        if scs:
            provides_sc_up = 0
            provides_sc_bp = 0
            provides_sc_user = 0

            sc_wrapper.up_comparative = [None] * fc.n_seq
            sc_wrapper.bp_comparative = [None] * fc.n_seq
            sc_wrapper.bp_local_comparative = [None] * fc.n_seq
            sc_wrapper.user_cb_comparative = [None] * fc.n_seq
            sc_wrapper.user_data_comparative = [None] * fc.n_seq

            for s in range(fc.n_seq):
                if scs[s]:
                    sliding_window = 1 if scs[s].type == VRNA_SC_WINDOW else 0
                    sc_wrapper.up_comparative[s] = scs[s].exp_energy_up
                    sc_wrapper.bp_comparative[s] = None if sliding_window else scs[s].exp_energy_bp
                    sc_wrapper.bp_local_comparative[s] = scs[s].exp_energy_bp_local if sliding_window else None
                    sc_wrapper.user_cb_comparative[s] = scs[s].exp_f
                    sc_wrapper.user_data_comparative[s] = scs[s].data

                    if scs[s].exp_energy_up:
                        provides_sc_up = 1

                    if sliding_window:
                        if scs[s].exp_energy_bp_local:
                            provides_sc_bp = 1
                    elif scs[s].exp_energy_bp:
                        provides_sc_bp = 1

                    if scs[s].exp_f:
                        provides_sc_user = 1

            if provides_sc_user:
                sc_wrapper.pair_ext = sc_hp_exp_cb_ext_user_comparative
                if provides_sc_up:
                    sc_wrapper.pair_ext = sc_hp_exp_cb_ext_up_user_comparative

                    if provides_sc_bp:
                        sc_wrapper.pair = sc_hp_exp_cb_up_bp_local_user_comparative if sliding_window else sc_hp_exp_cb_up_bp_user_comparative
                    else:
                        sc_wrapper.pair = sc_hp_exp_cb_up_user_comparative
                elif provides_sc_bp:
                    sc_wrapper.pair = sc_hp_exp_cb_bp_local_user_comparative if sliding_window else sc_hp_exp_cb_bp_user_comparative
                else:
                    sc_wrapper.pair = sc_hp_exp_cb_user_comparative
            else:
                if provides_sc_up:
                    sc_wrapper.pair_ext = sc_hp_exp_cb_ext_up_comparative
                    if provides_sc_bp:
                        sc_wrapper.pair = sc_hp_exp_cb_up_bp_local_comparative if sliding_window else sc_hp_exp_cb_up_bp_comparative
                    else:
                        sc_wrapper.pair = sc_hp_exp_cb_up_comparative
                elif provides_sc_bp:
                    sc_wrapper.pair = sc_hp_exp_cb_bp_local_comparative if sliding_window else sc_hp_exp_cb_bp_comparative

# init_sc_hp_exp end #################

# class for init_sc_int_exp start ###########
class sc_int_exp_dat(sc_hp_exp_dat):
    def __init__(self):
        super().__init__()
        self.stack = self.stack_comparative = 0.
        self.pair = self.pair_ext = sc_int_exp_cb

sc_int_exp_cb = Callable[[int, int, int, int, sc_int_exp_dat], float]

# function for init_sc_int_exp start ###########
def sc_int_exp_cb_up(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    u1 = k - i - 1
    u2 = j - l - 1
    sc = 1.0

    if u1 > 0:
        sc *= data.up[i + 1][u1]

    if u2 > 0:
        sc *= data.up[l + 1][u2]

    return sc

def sc_int_exp_cb_bp_local(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return data.bp_local[i][j - i]

def sc_int_exp_cb_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    sc = 1.0

    if i + 1 == k and l + 1 == j:
        sc *= data.stack[i] * data.stack[k] * data.stack[l] * data.stack[j]

    return sc

def sc_int_exp_cb_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return data.user_cb(i, j, k, l, VRNA_DECOMP_PAIR_IL, data.user_data)



def sc_int_exp_cb_up_bp_local_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp_local(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)


def sc_int_exp_cb_bp(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return data.bp[data.idx[j] + i]

def sc_int_exp_cb_up_bp_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)


def sc_int_exp_cb_ext_up(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    sc = 1.0

    u1 = i - 1
    u2 = k - j - 1
    u3 = data.n - l

    if u1 > 0:
        sc *= data.up[1][u1]

    if u2 > 0:
        sc *= data.up[j + 1][u2]

    if u3 > 0:
        sc *= data.up[l + 1][u3]

    return sc

def sc_int_exp_cb_ext_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    sc = 1.0

    if i == 1 and j + 1 == k and l == data.n:
        sc *= data.stack[i] * data.stack[k] * data.stack[l] * data.stack[j]

    return sc


def sc_int_exp_cb_ext_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return data.user_cb(i, j, k, l, VRNA_DECOMP_PAIR_IL, data.user_data)


def sc_int_exp_cb_up_bp_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_ext_up(i, j, k, l, data) * sc_int_exp_cb_ext_stack(i, j, k, l, data) * sc_int_exp_cb_ext_user(i, j, k, l, data)


def sc_int_exp_cb_ext_up_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_ext_up(i, j, k, l, data) * sc_int_exp_cb_ext_stack(i, j, k, l, data) * sc_int_exp_cb_ext_user(i, j, k, l, data)

def sc_int_exp_cb_up_bp_local_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp_local(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_up_bp_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_ext_up_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_ext_up(i, j, k, l, data) * sc_int_exp_cb_ext_user(i, j, k, l, data)

def sc_int_exp_cb_up_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_up_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_bp_local_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_bp_local(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_bp_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_bp(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_ext_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_ext_stack(i, j, k, l, data) * sc_int_exp_cb_ext_user(i, j, k, l, data)

def sc_int_exp_cb_bp_local_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_bp_local(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_bp_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_bp(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_stack_user(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_stack(i, j, k, l, data) * sc_int_exp_cb_user(i, j, k, l, data)

def sc_int_exp_cb_up_bp_local_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp_local(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data)

def sc_int_exp_cb_up_bp_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data)

def sc_int_exp_cb_ext_up_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_ext_up(i, j, k, l, data) * sc_int_exp_cb_ext_stack(i, j, k, l, data)

def sc_int_exp_cb_up_bp_local(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp_local(i, j, k, l, data)

def sc_int_exp_cb_up_bp(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_bp(i, j, k, l, data)

def sc_int_exp_cb_bp_local_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_bp_local(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data)

def sc_int_exp_cb_bp_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_bp(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data)

def sc_int_exp_cb_up_stack(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up(i, j, k, l, data) * sc_int_exp_cb_stack(i, j, k, l, data)


def sc_int_exp_cb_up_comparative(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.up_comparative[s]:
            u1 = data.a2s[s][k - 1] - data.a2s[s][i]
            u2 = data.a2s[s][j - 1] - data.a2s[s][l]

            if u1 > 0:
                sc *= data.up_comparative[s][data.a2s[s][i + 1]][u1]

            if u2 > 0:
                sc *= data.up_comparative[s][data.a2s[s][l + 1]][u2]

    return sc

def sc_int_exp_cb_bp_local_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    sc = 1.0
    
    for s in range(data.n_seq):
        if data.bp_local_comparative[s]:
            sc *= data.bp_local_comparative[s][i][j - i]
    
    return sc

def sc_int_exp_cb_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    sc = 1.0
    
    for s in range(data.n_seq):
        if data.stack_comparative[s]:
            if (data.a2s[s][k - 1] == data.a2s[s][i]) and (data.a2s[s][j - 1] == data.a2s[s][l]):
                sc *= (data.stack_comparative[s][data.a2s[s][i]] *
                       data.stack_comparative[s][data.a2s[s][k]] *
                       data.stack_comparative[s][data.a2s[s][l]] *
                       data.stack_comparative[s][data.a2s[s][j]])
    
    return sc

def sc_int_exp_cb_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    sc = 1.0
    
    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_PAIR_IL, data.user_data_comparative[s])
    
    return sc



def sc_int_exp_cb_up_bp_local_stack_user_comparative(i:int, j:int, k:int, l:int, data:sc_int_exp_dat) -> float:
    return sc_int_exp_cb_up_comparative(i, j, k, l, data) * sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) * sc_int_exp_cb_stack_comparative(i, j, k, l, data) * sc_int_exp_cb_user_comparative(i, j, k, l, data)

def sc_int_exp_cb_bp_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.bp_comparative[s]:
            sc *= data.bp_comparative[s][data.idx[j] + i]

    return sc



def sc_int_exp_cb_up_bp_stack_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_ext_up_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.up_comparative[s]:
            u1 = data.a2s[s][i - 1]
            u2 = data.a2s[s][k - 1] - data.a2s[s][j]
            u3 = data.a2s[s][data.n] - data.a2s[s][l]

            if u1 > 0:
                sc *= data.up_comparative[s][1][u1]

            if u2 > 0:
                sc *= data.up_comparative[s][data.a2s[s][j + 1]][u2]

            if u3 > 0:
                sc *= data.up_comparative[s][data.a2s[s][l + 1]][u3]

    return sc

def sc_int_exp_cb_ext_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.stack_comparative[s]:
            if (data.a2s[s][i] == 1 and
                data.a2s[s][j] == data.a2s[s][k - 1] and
                data.a2s[s][l] == data.a2s[s][data.n]):
                sc *= (data.stack_comparative[s][data.a2s[s][i]] *
                       data.stack_comparative[s][data.a2s[s][k]] *
                       data.stack_comparative[s][data.a2s[s][l]] *
                       data.stack_comparative[s][data.a2s[s][j]])

    return sc


def sc_int_exp_cb_ext_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    sc = 1.0

    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_PAIR_IL, data.user_data_comparative[s])

    return sc



def sc_int_exp_cb_ext_up_stack_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_ext_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_ext_stack_comparative(i, j, k, l, data) *
        sc_int_exp_cb_ext_user_comparative(i, j, k, l, data)
    )




def sc_int_exp_cb_up_bp_local_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_up_bp_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_ext_up_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_ext_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_ext_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_up_stack_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_up_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_bp_local_stack_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_bp_stack_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_ext_stack_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_ext_stack_comparative(i, j, k, l, data) *
        sc_int_exp_cb_ext_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_bp_local_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_bp_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_stack_user_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
        sc_int_exp_cb_user_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_up_bp_local_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_up_bp_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_ext_up_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_ext_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_ext_stack_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_up_bp_local_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_bp_local_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_up_bp_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_bp_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_bp_local_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data)
    )

def sc_int_exp_cb_bp_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data)
    )


def sc_int_exp_cb_up_stack_comparative(i: int, j: int, k: int, l: int, data: sc_int_exp_dat) -> float:
    return (
        sc_int_exp_cb_up_comparative(i, j, k, l, data) *
        sc_int_exp_cb_stack_comparative(i, j, k, l, data)
    )





# init_sc_int_exp start ###########
def init_sc_int_exp(fc:vrna_fold_compound_t, sc_wrapper:sc_int_exp_dat):
    sliding_window = 0
    provides_sc_up = 0
    provides_sc_bp = 0
    provides_sc_stack = 0
    provides_sc_user = 0

    if fc.exp_matrices:
        sliding_window = 1 if fc.exp_matrices.type == VRNA_MX_WINDOW else 0
    elif (fc.type == VRNA_FC_TYPE_SINGLE) and (fc.sc):
        sliding_window = 1 if fc.sc.type == VRNA_SC_WINDOW else 0
    elif fc.hc:
        sliding_window = 1 if fc.hc.type == VRNA_HC_WINDOW else 0

    sc_wrapper.n = fc.length
    sc_wrapper.n_seq = 1
    sc_wrapper.a2s = None
    sc_wrapper.idx = fc.jindx
    sc_wrapper.up = None
    sc_wrapper.up_comparative = None
    sc_wrapper.bp = None
    sc_wrapper.bp_comparative = None
    sc_wrapper.bp_local = None
    sc_wrapper.bp_local_comparative = None
    sc_wrapper.stack = None
    sc_wrapper.stack_comparative = None
    sc_wrapper.user_cb = None
    sc_wrapper.user_cb_comparative = None
    sc_wrapper.user_data = None
    sc_wrapper.user_data_comparative = None
    sc_wrapper.pair = None
    sc_wrapper.pair_ext = None

    if fc.type == VRNA_FC_TYPE_SINGLE:
        sc = fc.sc

        if sc:
            sc_wrapper.up = sc.exp_energy_up
            sc_wrapper.bp = None if sliding_window else sc.exp_energy_bp
            sc_wrapper.bp_local = sc.exp_energy_bp_local if sliding_window else None
            sc_wrapper.stack = sc.exp_energy_stack
            sc_wrapper.user_cb = sc.exp_f
            sc_wrapper.user_data = sc.data

            provides_sc_up = 1 if sc.exp_energy_up else 0
            if sliding_window:
                provides_sc_bp = 1 if sc.exp_energy_bp_local else 0
            else:
                provides_sc_bp = 1 if sc.exp_energy_bp else 0
            provides_sc_stack = 1 if sc.exp_energy_stack else 0
            provides_sc_user = 1 if sc.exp_f else 0

            if provides_sc_user:
                if provides_sc_up:
                    if provides_sc_bp:
                        if provides_sc_stack:
                            sc_wrapper.pair = sc_int_exp_cb_up_bp_local_stack_user if sliding_window else sc_int_exp_cb_up_bp_stack_user
                            sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack_user
                        else:
                            sc_wrapper.pair = sc_int_exp_cb_up_bp_local_user if sliding_window else sc_int_exp_cb_up_bp_user
                            sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_user
                    elif provides_sc_stack:
                        sc_wrapper.pair = sc_int_exp_cb_up_stack_user
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack_user
                    else:
                        sc_wrapper.pair = sc_int_exp_cb_up_user
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_user
                elif provides_sc_bp:
                    if provides_sc_stack:
                        sc_wrapper.pair = sc_int_exp_cb_bp_local_stack_user if sliding_window else sc_int_exp_cb_bp_stack_user
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack_user
                    else:
                        sc_wrapper.pair = sc_int_exp_cb_bp_local_user if sliding_window else sc_int_exp_cb_bp_user
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_user
                elif provides_sc_stack:
                    sc_wrapper.pair = sc_int_exp_cb_stack_user
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack_user
                else:
                    sc_wrapper.pair = sc_int_exp_cb_user
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_user
            elif provides_sc_bp:
                if provides_sc_up:
                    if provides_sc_stack:
                        sc_wrapper.pair = sc_int_exp_cb_up_bp_local_stack if sliding_window else sc_int_exp_cb_up_bp_stack
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack
                    else:
                        sc_wrapper.pair = sc_int_exp_cb_up_bp_local if sliding_window else sc_int_exp_cb_up_bp
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up
                elif provides_sc_stack:
                    sc_wrapper.pair = sc_int_exp_cb_bp_local_stack if sliding_window else sc_int_exp_cb_bp_stack
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack
                else:
                    sc_wrapper.pair = sc_int_exp_cb_bp_local if sliding_window else sc_int_exp_cb_bp
            elif provides_sc_up:
                if provides_sc_stack:
                    sc_wrapper.pair = sc_int_exp_cb_up_stack
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack
                else:
                    sc_wrapper.pair = sc_int_exp_cb_up
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_up
            elif provides_sc_stack:
                sc_wrapper.pair = sc_int_exp_cb_stack
                sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack

    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        sc_wrapper.n_seq = fc.n_seq
        sc_wrapper.a2s = fc.a2s
        scs = fc.scs

        if scs:
            sc_wrapper.up_comparative = [None] * fc.n_seq
            sc_wrapper.bp_comparative = [None] * fc.n_seq
            sc_wrapper.bp_local_comparative = [None] * fc.n_seq
            sc_wrapper.stack_comparative = [None] * fc.n_seq
            sc_wrapper.user_cb_comparative = [None] * fc.n_seq
            sc_wrapper.user_data_comparative = [None] * fc.n_seq

            for s in range(fc.n_seq):
                if scs[s]:
                    sliding_window = 1 if scs[s].type == VRNA_SC_WINDOW else 0
                    sc_wrapper.up_comparative[s] = scs[s].exp_energy_up
                    sc_wrapper.bp_comparative[s] = None if sliding_window else scs[s].exp_energy_bp
                    sc_wrapper.bp_local_comparative[s] = scs[s].exp_energy_bp_local if sliding_window else None
                    sc_wrapper.stack_comparative[s] = scs[s].exp_energy_stack
                    sc_wrapper.user_cb_comparative[s] = scs[s].exp_f
                    sc_wrapper.user_data_comparative[s] = scs[s].data

                    provides_sc_up = 1 if scs[s].exp_energy_up else 0
                    if sliding_window:
                        provides_sc_bp = 1 if scs[s].exp_energy_bp_local else 0
                    else:
                        provides_sc_bp = 1 if scs[s].exp_energy_bp else 0
                    provides_sc_stack = 1 if scs[s].exp_energy_stack else 0
                    provides_sc_user = 1 if scs[s].exp_f else 0

            if provides_sc_user:
                if provides_sc_up:
                    if provides_sc_bp:
                        if provides_sc_stack:
                            sc_wrapper.pair = sc_int_exp_cb_up_bp_local_stack_user_comparative if sliding_window else sc_int_exp_cb_up_bp_stack_user_comparative
                            sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack_user_comparative
                        else:
                            sc_wrapper.pair = sc_int_exp_cb_up_bp_local_user_comparative if sliding_window else sc_int_exp_cb_up_bp_user_comparative
                            sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_user_comparative
                    elif provides_sc_stack:
                        sc_wrapper.pair = sc_int_exp_cb_up_stack_user_comparative
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack_user_comparative
                    else:
                        sc_wrapper.pair = sc_int_exp_cb_up_user_comparative
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_user_comparative
                elif provides_sc_bp:
                    if provides_sc_stack:
                        sc_wrapper.pair = sc_int_exp_cb_bp_local_stack_user_comparative if sliding_window else sc_int_exp_cb_bp_stack_user_comparative
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack_user_comparative
                    else:
                        sc_wrapper.pair = sc_int_exp_cb_bp_local_user_comparative if sliding_window else sc_int_exp_cb_bp_user_comparative
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_user_comparative
                elif provides_sc_stack:
                    sc_wrapper.pair = sc_int_exp_cb_stack_user_comparative
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack_user_comparative
                else:
                    sc_wrapper.pair = sc_int_exp_cb_user_comparative
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_user_comparative
            elif provides_sc_bp:
                if provides_sc_up:
                    if provides_sc_stack:
                        sc_wrapper.pair = sc_int_exp_cb_up_bp_local_stack_comparative if sliding_window else sc_int_exp_cb_up_bp_stack_comparative
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack_comparative
                    else:
                        sc_wrapper.pair = sc_int_exp_cb_up_bp_local_comparative if sliding_window else sc_int_exp_cb_up_bp_comparative
                        sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_comparative
                elif provides_sc_stack:
                    sc_wrapper.pair = sc_int_exp_cb_bp_local_stack_comparative if sliding_window else sc_int_exp_cb_bp_stack_comparative
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack_comparative
                else:
                    sc_wrapper.pair = sc_int_exp_cb_bp_local_comparative if sliding_window else sc_int_exp_cb_bp_comparative
            elif provides_sc_up:
                if provides_sc_stack:
                    sc_wrapper.pair = sc_int_exp_cb_up_stack_comparative
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_stack_comparative
                else:
                    sc_wrapper.pair = sc_int_exp_cb_up_comparative
                    sc_wrapper.pair_ext = sc_int_exp_cb_ext_up_comparative
            elif provides_sc_stack:
                sc_wrapper.pair = sc_int_exp_cb_stack_comparative
                sc_wrapper.pair_ext = sc_int_exp_cb_ext_stack_comparative
# init_sc_int_exp end ################

# init_sc_mb_exp class start ################

    

class sc_mb_exp_dat(sc_int_exp_dat):
    def __init__(self):
        super().__init__()
        self.red_stem = self.red_ml = self.decomp_ml = sc_mb_exp_red_cb
        self.pair = self.pair_ext = sc_mb_exp_pair_cb
        
sc_mb_exp_red_cb = Callable[[int, int, int, int, sc_mb_exp_dat], float]   
sc_mb_exp_pair_cb = Callable[[int, int, sc_mb_exp_dat], float]   


# init_sc_mb_exp method start ################
def sc_mb_exp_pair_cb_bp(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return data.bp[data.idx[j] + i]


def sc_mb_exp_pair_cb_bp_comparative(i: int, j: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.bp_comparative[s]:
            sc *= data.bp_comparative[s][data.idx[j] + i]
    return sc


def sc_mb_exp_pair_cb_bp_local(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return data.bp_local[i][j - i]


def sc_mb_exp_pair_cb_bp_local_comparative(i: int, j: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.bp_local_comparative[s]:
            sc *= data.bp_local_comparative[s][i][j - i]
    return sc


def sc_mb_exp_pair_cb_user(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return data.user_cb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, data.user_data)


def sc_mb_exp_pair_ext_cb_user(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return data.user_cb(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, data.user_data)


def sc_mb_exp_pair_cb_user_comparative(i: int, j: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, data.user_data_comparative[s])
    return sc


def sc_mb_exp_pair_ext_cb_user_comparative(i: int, j: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML, data.user_data_comparative[s])
    return sc


def sc_mb_exp_pair_cb_bp_user(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_pair_cb_bp(i, j, data) * sc_mb_exp_pair_cb_user(i, j, data)


def sc_mb_exp_pair_cb_bp_user_comparative(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_pair_cb_bp_comparative(i, j, data) * sc_mb_exp_pair_cb_user_comparative(i, j, data)


def sc_mb_exp_pair_cb_bp_local_user(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_pair_cb_bp_local(i, j, data) * sc_mb_exp_pair_cb_user(i, j, data)


def sc_mb_exp_pair_cb_bp_local_user_comparative(i: int, j: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_pair_cb_bp_local_comparative(i, j, data) * sc_mb_exp_pair_cb_user_comparative(i, j, data)


def sc_mb_exp_red_cb_up(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    l1 = k - i
    l2 = j - l
    sc = 1.0
    if l1 > 0:
        sc *= data.up[i][l1]
    if l2 > 0:
        sc *= data.up[l + 1][l2]
    return sc


def sc_mb_exp_red_cb_up_comparative(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.up_comparative[s]:
            l1 = data.a2s[s][k] - data.a2s[s][i]
            l2 = data.a2s[s][j] - data.a2s[s][l]
            start2 = data.a2s[s][l] + 1
            if l1 > 0:
                sc *= data.up_comparative[s][data.a2s[s][i]][l1]
            if l2 > 0:
                sc *= data.up_comparative[s][start2][l2]
    return sc


def sc_mb_exp_red_cb_user(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    return data.user_cb(i, j, k, l, VRNA_DECOMP_ML_ML, data.user_data)


def sc_mb_exp_red_cb_user_comparative(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_ML_ML, data.user_data_comparative[s])
    return sc


def sc_mb_exp_red_cb_up_user(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_red_cb_up(i, j, k, l, data) * sc_mb_exp_red_cb_user(i, j, k, l, data)


def sc_mb_exp_red_cb_up_user_comparative(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_red_cb_up_comparative(i, j, k, l, data) * sc_mb_exp_red_cb_user_comparative(i, j, k, l, data)


def sc_mb_exp_red_cb_stem_user(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    return data.user_cb(i, j, k, l, VRNA_DECOMP_ML_STEM, data.user_data)


def sc_mb_exp_red_cb_stem_user_comparative(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_ML_STEM, data.user_data)
    return sc


def sc_mb_exp_red_cb_stem_up_user(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_red_cb_up(i, j, k, l, data) * sc_mb_exp_red_cb_stem_user(i, j, k, l, data)


def sc_mb_exp_red_cb_stem_up_user_comparative(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    return sc_mb_exp_red_cb_up_comparative(i, j, k, l, data) * sc_mb_exp_red_cb_stem_user_comparative(i, j, k, l, data)


def sc_mb_exp_split_cb_user(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    return data.user_cb(i, j, k, l, VRNA_DECOMP_ML_ML_ML, data.user_data)


def sc_mb_exp_split_cb_user_comparative(i: int, j: int, k: int, l: int, data: sc_mb_exp_dat) -> float:
    sc = 1.0
    for s in range(data.n_seq):
        if data.user_cb_comparative[s]:
            sc *= data.user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_ML_ML_ML, data.user_data_comparative[s])
    return sc












# init_sc_mb_exp start ################
def init_sc_mb_exp(fc:vrna_fold_compound_t, sc_wrapper: sc_mb_exp_dat) -> None:
    sliding_window: int
    sc: vrna_sc_t
    scs: Optional[List[vrna_sc_t]]

    sc_wrapper.n = fc.length
    sc_wrapper.n_seq = 1
    sc_wrapper.idx = fc.jindx
    sc_wrapper.a2s = None

    sc_wrapper.up = None
    sc_wrapper.up_comparative = None
    sc_wrapper.bp = None
    sc_wrapper.bp_comparative = None
    sc_wrapper.bp_local = None
    sc_wrapper.bp_local_comparative = None

    sc_wrapper.user_cb = None
    sc_wrapper.user_data = None
    sc_wrapper.user_cb_comparative = None
    sc_wrapper.user_data_comparative = None

    sc_wrapper.pair = None
    sc_wrapper.pair_ext = None
    sc_wrapper.red_stem = None
    sc_wrapper.red_ml = None
    sc_wrapper.decomp_ml = None

    sliding_window = 1 if fc.hc.type == VRNA_HC_WINDOW else 0

    if fc.type == VRNA_FC_TYPE_SINGLE:
        sc = fc.sc

        if sc:
            provides_sc_up = 0
            provides_sc_bp = 0
            provides_sc_user = 0

            sc_wrapper.up = sc.exp_energy_up
            sc_wrapper.user_cb = sc.exp_f
            sc_wrapper.user_data = sc.data

            if sliding_window:
                sc_wrapper.bp_local = sc.exp_energy_bp_local
            else:
                sc_wrapper.bp = sc.exp_energy_bp

            if sc.exp_energy_up:
                provides_sc_up = 1

            if sliding_window:
                if sc.exp_energy_bp_local:
                    provides_sc_bp = 1
            elif sc.exp_energy_bp:
                provides_sc_bp = 1

            if sc.exp_f:
                provides_sc_user = 1

            if provides_sc_user:
                sc_wrapper.decomp_ml = sc_mb_exp_split_cb_user
                sc_wrapper.red_stem = sc_mb_exp_red_cb_stem_user
                sc_wrapper.red_ml = sc_mb_exp_red_cb_user
                sc_wrapper.pair = sc_mb_exp_pair_cb_user
                if not sliding_window:
                    sc_wrapper.pair_ext = sc_mb_exp_pair_ext_cb_user

                if provides_sc_bp:
                    if sliding_window:
                        sc_wrapper.pair = sc_mb_exp_pair_cb_bp_local_user
                    else:
                        sc_wrapper.pair = sc_mb_exp_pair_cb_bp_user
                        sc_wrapper.pair_ext = sc_mb_exp_pair_ext_cb_user

                if provides_sc_up:
                    sc_wrapper.red_stem = sc_mb_exp_red_cb_stem_up_user
                    sc_wrapper.red_ml = sc_mb_exp_red_cb_up_user

            elif provides_sc_bp:
                if sliding_window:
                    sc_wrapper.pair = sc_mb_exp_pair_cb_bp_local
                else:
                    sc_wrapper.pair = sc_mb_exp_pair_cb_bp

                if provides_sc_up:
                    sc_wrapper.red_stem = sc_mb_exp_red_cb_up
                    sc_wrapper.red_ml = sc_mb_exp_red_cb_up

            elif provides_sc_up:
                sc_wrapper.red_stem = sc_mb_exp_red_cb_up
                sc_wrapper.red_ml = sc_mb_exp_red_cb_up

    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        sc_wrapper.a2s = fc.a2s
        sc_wrapper.n_seq = fc.n_seq
        scs = fc.scs

        if scs:
            provides_sc_up = 0
            provides_sc_bp = 0
            provides_sc_user = 0

            sc_wrapper.up_comparative = [[] for _ in range(fc.n_seq)]
            sc_wrapper.bp_comparative = [None] * fc.n_seq
            sc_wrapper.bp_local_comparative = [[] for _ in range(fc.n_seq)]
            sc_wrapper.user_cb_comparative = [None] * fc.n_seq
            sc_wrapper.user_data_comparative = [None] * fc.n_seq

            for s in range(fc.n_seq):
                if scs[s]:
                    sc_wrapper.up_comparative[s] = scs[s].exp_energy_up
                    sc_wrapper.bp_comparative[s] = None if sliding_window else scs[s].exp_energy_bp
                    sc_wrapper.bp_local_comparative[s] = scs[s].exp_energy_bp_local if sliding_window else None
                    sc_wrapper.user_cb_comparative[s] = scs[s].exp_f
                    sc_wrapper.user_data_comparative[s] = scs[s].data

                    if scs[s].exp_energy_up:
                        provides_sc_up = 1

                    if sliding_window:
                        if scs[s].exp_energy_bp_local:
                            provides_sc_bp = 1
                    elif scs[s].exp_energy_bp:
                        provides_sc_bp = 1

                    if scs[s].exp_f:
                        provides_sc_user = 1

            if provides_sc_user:
                sc_wrapper.decomp_ml = sc_mb_exp_split_cb_user_comparative
                sc_wrapper.red_stem = sc_mb_exp_red_cb_stem_user_comparative
                sc_wrapper.red_ml = sc_mb_exp_red_cb_user_comparative
                sc_wrapper.pair = sc_mb_exp_pair_cb_user_comparative
                if not sliding_window:
                    sc_wrapper.pair_ext = sc_mb_exp_pair_ext_cb_user_comparative

                if provides_sc_bp:
                    if sliding_window:
                        sc_wrapper.pair = sc_mb_exp_pair_cb_bp_local_user_comparative
                    else:
                        sc_wrapper.pair = sc_mb_exp_pair_cb_bp_user_comparative
                        sc_wrapper.pair_ext = sc_mb_exp_pair_ext_cb_user_comparative

                if provides_sc_up:
                    sc_wrapper.red_stem = sc_mb_exp_red_cb_stem_up_user_comparative
                    sc_wrapper.red_ml = sc_mb_exp_red_cb_up_user_comparative

            elif provides_sc_bp:
                if sliding_window:
                    sc_wrapper.pair = sc_mb_exp_pair_cb_bp_local_comparative
                else:
                    sc_wrapper.pair = sc_mb_exp_pair_cb_bp_comparative

                if provides_sc_up:
                    sc_wrapper.red_stem = sc_mb_exp_red_cb_up_comparative
                    sc_wrapper.red_ml = sc_mb_exp_red_cb_up_comparative

            elif provides_sc_up:
                sc_wrapper.red_stem = sc_mb_exp_red_cb_up_comparative
                sc_wrapper.red_ml = sc_mb_exp_red_cb_up_comparative








def get_constraints_helper(fc):
    helpers = constraints_helper()

    helpers.hc_eval_ext = prepare_hc_ext_def(fc, helpers.hc_dat_ext)
    helpers.hc_eval_hp = prepare_hc_hp_def(fc, helpers.hc_dat_hp)
    helpers.hc_eval_int = prepare_hc_int_def(fc, helpers.hc_dat_int)
    helpers.hc_eval_mb = prepare_hc_mb_def(fc, helpers.hc_dat_mb)

    init_sc_ext_exp(fc, helpers.sc_wrapper_ext)
    init_sc_hp_exp(fc, helpers.sc_wrapper_hp)
    init_sc_int_exp(fc, helpers.sc_wrapper_int)
    init_sc_mb_exp(fc, helpers.sc_wrapper_mb)

    return helpers


# compute_bpp_internal constant and class start ###############
MAXALPHA = 20
VRNA_MODEL_DEFAULT_SALT = 1.021
K0 = 273.15


class vrna_ep_t:
    def __init__(self) -> None:
        self.i = self.j = self.type = 0
        self.p = 0.

class vrna_md_t:
    def __init__(self):
        self.temperature = 0.0
        self.betaScale = 0.0
        self.pf_smooth = 0
        self.dangles = 0
        self.special_hp = 0
        self.noLP = 0
        self.noGU = 0
        self.noGUclosure = 0
        self.logML = 0
        self.circ = 0
        self.gquad = 0
        self.uniq_ML = 0
        self.energy_set = 0
        self.backtrack = 0
        self.backtrack_type = ''
        self.compute_bpp = 0
        self.nonstandards = ''
        self.max_bp_span = 0
        self.min_loop_size = 0
        self.window_size = 0
        self.oldAliEn = 0
        self.ribo = 0
        self.cv_fact = 0.0
        self.nc_fact = 0.0
        self.sfact = 0.0
        self.rtype = [0] * 8
        self.alias = [0] * (MAXALPHA + 1)
        self.pair = [[0] * (MAXALPHA + 1) for _ in range(MAXALPHA + 1)]
        self.pair_dist = [[0.0] * 7 for _ in range(7)]
        self.salt = 0.0
        self.saltMLLower = 0
        self.saltMLUpper = 0
        self.saltDPXInit = 0
        self.saltDPXInitFact = 0.0
        self.helical_rise = 0.0
        self.backbone_length = 0.0

# compute_bpp_internal method start ###############
def vrna_get_ptype_md(i: int, j: int, md: vrna_md_t) -> int:
    tt = int(md.pair[i][j])
    return 7 if tt == 0 else tt


def vrna_get_ptype(ij: int, ptype: str) -> int:
    tt = int(ptype[ij])
    return 7 if tt == 0 else tt


################################################ build in remove ########################
GASCONST = 1.98717
Eular_const = 0.58
PI = 3.141592653589793

# math utils
EUL = 0.57721566490153286060
MAXLOG = 709.782712893384
MACHEP = 2.2204460492503131e-16
MAXNUM = 1.7976931348623157e+308
MAXFAC = 31

def tau_ss(T, backbonelen):
    bjerrum_length_inv = 1 / bjerrum_length(T)
    return min(1 / backbonelen, bjerrum_length_inv)

def approx_hyper(y):
    a = 1 / (pow(y, 6.)/pow(2 * PI, 6.) + 1)
    b = pow(y, 4.)/(36 * pow(PI, 4.)) - pow(y, 3.) / (24 * PI * PI) + y * y / (2 * PI * PI) - y / 2
    c = math.log(2 * PI / y) - 1.96351
    return a * b + (1 - a) * c


def ionic_strength(rho):
    return rho


def epsilonr(T):
    return 5321 / T + 233.76 - 0.9297 * T + 1.417 * T * T / 1000 - 0.8292 * T * T * T / 1000000
    
    
def bjerrum_length(T):
    return 167100.052/(T * epsilonr(T))

def kappa(rho, T):
        return math.sqrt(bjerrum_length(T) * ionic_strength(rho)) / 8.1284

def loop_salt_aux(kmlss, L, T, backbonelen):
        a = (GASCONST / 1000.) * T * bjerrum_length(T) * L * backbonelen * tau_ss(T, backbonelen) * tau_ss(T, backbonelen)
        b = math.log(kmlss) - math.log(PI / 2) + Eular_const + approx_hyper(kmlss) + 1 / kmlss * (1 - math.exp(-kmlss) + kmlss * expn(1, kmlss))
        return a * b * 100

def vrna_salt_loop(L, rho, T, backbonelen):
        if L == 0:return 0
        kmlss_ref = kappa(VRNA_MODEL_DEFAULT_SALT, T) * L * backbonelen
        kmlss = kappa(rho, T) * L * backbonelen
        correction = loop_salt_aux(kmlss, L, T, backbonelen) - loop_salt_aux(kmlss_ref, L, T, backbonelen)
        return correction

def vrna_salt_loop_int(L, rho, T, backbonelen):
        correction = vrna_salt_loop(L, rho, T, backbonelen)
        return int(correction + 0.5 - (correction < 0))

################################################ build in remove ########################


def exp_E_IntLoop(u1: int,
                  u2: int,
                  type: int,
                  type2: int,
                  si1: int,
                  sj1: int,
                  sp1: int,
                  sq1: int,
                  P: vrna_exp_param_t) -> float:
    ul: int
    us: int
    no_close: int = 0
    z: float = 0.0
    noGUclosure: int = P.model_details.noGUclosure
    backbones: int
    salt_stack_correction: float = P.expSaltStack
    salt_loop_correction: float = 1.0

    if (noGUclosure) and ((type2 == 3) or (type2 == 4) or (type == 3) or (type == 4)):
        no_close = 1

    if u1 > u2:
        ul = u1
        us = u2
    else:
        ul = u2
        us = u1

    # salt correction for loop
    backbones = ul + us + 2

    if P.model_details.salt != VRNA_MODEL_DEFAULT_SALT:
        if backbones <= MAXLOOP + 1:
            salt_loop_correction = P.expSaltLoop[backbones]
        else:
            salt_loop_correction = math.exp(-vrna_salt_loop_int(backbones, P.model_details.salt, P.temperature + K0, P.model_details.backbone_length) * 10.0 / P.kT)

    if ul == 0:
        # stack
        z = P.expstack[type][type2] * salt_stack_correction
    elif not no_close:
        if us == 0:
            # bulge
            z = P.expbulge[ul]
            if ul == 1:
                z *= P.expstack[type][type2]
            else:
                if type > 2:
                    z *= P.expTermAU
                if type2 > 2:
                    z *= P.expTermAU
            return z * salt_loop_correction
        elif us == 1:
            if ul == 1:
                # 1x1 loop
                return P.expint11[type][type2][si1][sj1] * salt_loop_correction
            elif ul == 2:
                # 2x1 loop
                if u1 == 1:
                    return P.expint21[type][type2][si1][sq1][sj1] * salt_loop_correction
                else:
                    return P.expint21[type2][type][sq1][si1][sp1] * salt_loop_correction
            else:
                # 1xn loop
                z = P.expinternal[ul + us] * P.expmismatch1nI[type][si1][sj1] * P.expmismatch1nI[type2][sq1][sp1]
                return z * P.expninio[2][ul - us] * salt_loop_correction
        elif us == 2:
            if ul == 2:
                # 2x2 loop
                return P.expint22[type][type2][si1][sp1][sq1][sj1] * salt_loop_correction
            elif ul == 3:
                # 2x3 loop
                z = P.expinternal[5] * P.expmismatch23I[type][si1][sj1] * P.expmismatch23I[type2][sq1][sp1]
                return z * P.expninio[2][1] * salt_loop_correction

        # generic interior loop (no else here!)
        z = P.expinternal[ul + us] * P.expmismatchI[type][si1][sj1] * P.expmismatchI[type2][sq1][sp1]
        return z * P.expninio[2][ul - us] * salt_loop_correction

    return z

    
def compute_gquad_prob_internal(fc: vrna_fold_compound_t, l: int) -> None:
    n: int = fc.length
    S1: List[int] = fc.sequence_encoding
    ptype: str = fc.ptype
    my_iindx: List[int] = fc.iindx
    jindx: List[int] = fc.jindx
    pf_params: vrna_exp_param_t = fc.exp_params
    G: List[float] = fc.exp_matrices.G
    probs: List[float] = fc.exp_matrices.probs
    scale: List[float] = fc.exp_matrices.scale

    expintern: List[float] = pf_params.expinternal

    if l < n - 3:
        for k in range(2, l - VRNA_GQUAD_MIN_BOX_SIZE + 2):
            kl: int = my_iindx[k] - l
            if G[kl] == 0.0:
                continue

            tmp2: float = 0.0
            i: int = k - 1
            for j in range(min(l + MAXLOOP + 1, n), l + 3, -1):
                ij: int = my_iindx[i] - j
                type: int = ord(ptype[jindx[j] + i])
                if not type:
                    continue

                u1: int = j - l - 1
                qe: float = pf_params.expTermAU if type > 2 else 1.0
                tmp2 += probs[ij] * qe * expintern[u1] * pf_params.expmismatchI[type][S1[i + 1]][S1[j - 1]] * scale[u1 + 2]

            probs[kl] += tmp2 * G[kl]

    if l < n - 1:
        for k in range(3, l - VRNA_GQUAD_MIN_BOX_SIZE + 2):
            kl: int = my_iindx[k] - l
            if G[kl] == 0.0:
                continue

            tmp2: float = 0.0
            for i in range(max(1, k - MAXLOOP - 1), k - 1):
                u1: int = k - i - 1
                for j in range(l + 2, min(l + MAXLOOP - u1 + 1, n) + 1):
                    ij: int = my_iindx[i] - j
                    type: int = ord(ptype[jindx[j] + i])
                    if not type:
                        continue

                    u2: int = j - l - 1
                    qe: float = pf_params.expTermAU if type > 2 else 1.0
                    tmp2 += probs[ij] * qe * (expintern[u1 + u2]) * pf_params.expmismatchI[type][S1[i + 1]][S1[j - 1]] * scale[u1 + u2 + 2]

            probs[kl] += tmp2 * G[kl]

    if l < n:
        for k in range(4, l - VRNA_GQUAD_MIN_BOX_SIZE + 2):
            kl: int = my_iindx[k] - l
            if G[kl] == 0.0:
                continue

            tmp2: float = 0.0
            j: int = l + 1
            for i in range(max(1, k - MAXLOOP - 1), k - 3):
                ij: int = my_iindx[i] - j
                type: int = ord(ptype[jindx[j] + i])
                if not type:
                    continue

                u2: int = k - i - 1
                qe: float = pf_params.expTermAU if type > 2 else 1.0
                tmp2 += probs[ij] * qe * (expintern[u2]) * pf_params.expmismatchI[type][S1[i + 1]][S1[j - 1]] * scale[u2 + 2]

            probs[kl] += tmp2 * G[kl]

def compute_gquad_prob_internal_comparative(fc: vrna_fold_compound_t, l: int) -> None:
    n: int = fc.length
    n_seq: int = fc.n_seq
    S: List[List[int]] = fc.S
    S5: List[List[int]] = fc.S5
    S3: List[List[int]] = fc.S3
    a2s: List[List[int]] = fc.a2s
    my_iindx: List[int] = fc.iindx
    pf_params: vrna_exp_param_t = fc.exp_params
    G: List[float] = fc.exp_matrices.G
    qb: List[float] = fc.exp_matrices.qb
    probs: List[float] = fc.exp_matrices.probs
    scale: List[float] = fc.exp_matrices.scale
    md: vrna_md_t = pf_params.model_details

    expintern: List[float] = pf_params.expinternal

    if l < n - 3:
        for k in range(2, l - VRNA_GQUAD_MIN_BOX_SIZE + 2):
            kl: int = my_iindx[k] - l
            if G[kl] == 0.0:
                continue

            tmp2: float = 0.0
            i: int = k - 1
            for j in range(min(l + MAXLOOP + 1, n), l + 3, -1):
                ij: int = my_iindx[i] - j
                if qb[ij] == 0.0:
                    continue

                qe: float = 1.0
                u1: int = j - l - 1

                for s in range(n_seq):
                    type: int = vrna_get_ptype_md(S[s][i], S[s][j], md)
                    u1_local: int = a2s[s][j - 1] - a2s[s][l]
                    qe *= (expintern[u1_local])

                    if md.dangles == 2:
                        qe *= (pf_params.expmismatchI[type][S3[s][i]][S5[s][j]])

                    if type > 2:
                        qe *= (pf_params.expTermAU)

                tmp2 += probs[ij] * qe * scale[u1 + 2]

            probs[kl] += tmp2 * G[kl]

    if l < n - 1:
        for k in range(3, l - VRNA_GQUAD_MIN_BOX_SIZE + 2):
            kl: int = my_iindx[k] - l
            if G[kl] == 0.0:
                continue

            tmp2: float = 0.0
            for i in range(max(1, k - MAXLOOP - 1), k):
                u1: int = k - i - 1
                for j in range(l + 2, min(l + MAXLOOP - u1 + 1, n) + 1):
                    ij: int = my_iindx[i] - j
                    if qb[ij] == 0.0:
                        continue

                    qe: float = 1.0
                    u2: int = j - l - 1

                    for s in range(n_seq):
                        type: int = vrna_get_ptype_md(S[s][i], S[s][j], md)
                        u1_local: int = a2s[s][k - 1] - a2s[s][i]
                        u2_local: int = a2s[s][j - 1] - a2s[s][l]
                        qe *= (expintern[u1_local + u2_local])

                        if md.dangles == 2:
                            qe *= (pf_params.expmismatchI[type][S3[s][i]][S5[s][j]])

                        if type > 2:
                            qe *= (pf_params.expTermAU)

                    tmp2 += probs[ij] * qe * scale[u1 + u2 + 2]

            probs[kl] += tmp2 * G[kl]

    if l < n:
        for k in range(4, l - VRNA_GQUAD_MIN_BOX_SIZE + 2):
            kl: int = my_iindx[k] - l
            if G[kl] == 0.0:
                continue

            tmp2: float = 0.0
            j: int = l + 1
            for i in range(max(1, k - MAXLOOP - 1), k - 3):
                ij: int = my_iindx[i] - j
                if qb[ij] == 0.0:
                    continue

                qe: float = 1.0
                u2: int = k - i - 1

                for s in range(n_seq):
                    type: int = vrna_get_ptype_md(S[s][i], S[s][j], md)
                    u2_local: int = a2s[s][k - 1] - a2s[s][i]
                    qe *= (expintern[u2_local])

                    if md.dangles == 2:
                        qe *= (pf_params.expmismatchI[type][S3[s][i]][S5[s][j]])

                    if type > 2:
                        qe *= (pf_params.expTermAU)

                tmp2 += probs[ij] * qe * scale[u2 + 2]

            probs[kl] += tmp2 * G[kl]




def compute_bpp_internal(fc:vrna_fold_compound_t,
                         l,
                         bp_correction:list[vrna_ep_t],
                         corr_cnt,
                         corr_size,
                         Qmax,
                         ov,
                         constraints:constraints_helper):
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

# compute_bpp_multibranch constant and class start ####################
class helper_arrays:
    def __init__(self) -> None:
        self. prm_l = self.prm_l1 = self.prml = self.pmlu = self.prm_MLbu = 0.
        self.ud_max_size = 0



# compute_bpp_multibranch method start ####################
def exp_E_MLstem(type: int, si1: int, sj1: int, P: vrna_exp_param_t) -> float:
    energy: float = 1.0

    if si1 >= 0 and sj1 >= 0:
        energy = P.expmismatchM[type][si1][sj1]
    elif si1 >= 0:
        energy = P.expdangle5[type][si1]
    elif sj1 >= 0:
        energy = P.expdangle3[type][sj1]

    if type > 2:
        energy *= P.expTermAU

    energy *= P.expMLintern[type]
    return energy  # Casting to FLT_OR_DBL not necessary in Python

def rotate_ml_helper_arrays_outer(ml_helpers:helper_arrays):
    # rotate prm_l and prm_l1 arrays
    tmp = ml_helpers.prm_l1
    ml_helpers.prm_l1 = ml_helpers.prm_l
    ml_helpers.prm_l = tmp
    
    # rotate pmlu entries required for unstructured domain feature
    if ml_helpers.pmlu:
        tmp = ml_helpers.pmlu[ml_helpers.ud_max_size]
        
        for u in range(ml_helpers.ud_max_size, 0, -1):
            ml_helpers.pmlu[u] = ml_helpers.pmlu[u - 1]
        
        ml_helpers.pmlu[0] = tmp
        
        for u in range(ml_helpers.ud_max_size + 1):
            ml_helpers.prm_MLbu[u] = 0.0

def rotate_ml_helper_arrays_inner(ml_helpers:helper_arrays):
    # rotate prm_MLbu entries required for unstructured domain feature
    if ml_helpers.prm_MLbu:
        for u in range(ml_helpers.ud_max_size, 0, -1):
            ml_helpers.prm_MLbu[u] = ml_helpers.prm_MLbu[u - 1]


# compute_bpp_multibranch  start ####################
def compute_bpp_multibranch(fc:vrna_fold_compound_t,
                            l,
                            ml_helpers:helper_arrays,
                            Qmax,
                            ov,
                            constraints:constraints_helper):
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

# bppm_circ constant and class start ######################




# bppm_circ method start ######################
def numerator_single(vc: vrna_fold_compound_t, i: int, j: int) -> float:
    return 1.0

def numerator_comparative(vc: vrna_fold_compound_t, i: int, j: int) -> float:
    pscore: list[int] = vc.pscore  # precomputed array of pair types
    kTn: float = vc.exp_params.kT / 10.0  # kT in cal/mol
    jindx: list[int] = vc.jindx

    return math.exp(pscore[jindx[j] + i] / kTn)



import numpy as np

# bppm_circ  start ######################
def bppm_circ(fc, constraints):
    n = fc.length
    n_seq = 1 if fc.type == VRNA_FC_TYPE_SINGLE else fc.n_seq
    SS = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.S
    S5 = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.S5
    S3 = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.S3
    a2s = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.a2s
    pf_params = fc.exp_params
    md = pf_params.model_details
    S = fc.sequence_encoding2 if fc.type == VRNA_FC_TYPE_SINGLE else None
    S1 = fc.sequence_encoding if fc.type == VRNA_FC_TYPE_SINGLE else None
    my_iindx = fc.iindx
    jindx = fc.jindx
    ptype = fc.ptype if fc.type == VRNA_FC_TYPE_SINGLE else None
    hc = fc.hc
    matrices = fc.exp_matrices
    qb = matrices.qb
    qm = matrices.qm
    qm1 = matrices.qm1
    probs = matrices.probs
    scale = matrices.scale
    expMLbase = matrices.expMLbase
    qo = matrices.qo
    hard_constraints = hc.mx
    hc_dat_mb = constraints.hc_dat_mb
    hc_eval_mb = constraints.hc_eval_mb
    sc_dat_int = constraints.sc_wrapper_int
    sc_dat_mb = constraints.sc_wrapper_mb

    expMLclosing = pf_params.expMLclosing
    rtype = pf_params.model_details.rtype

    if fc.type == VRNA_FC_TYPE_SINGLE:
        numerator_f = numerator_single
        tt = None
    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        numerator_f = numerator_comparative
        tt = np.zeros(n_seq, dtype=int)
    else:
        numerator_f = None

    for i in range(1, n + 1):
        probs[my_iindx[i] - i] = 0
        for j in range(i + 1, n + 1):
            ij = my_iindx[i] - j
            if qb[ij] > 0:
                probs[ij] = numerator_f(fc, i, j) / qo

                if fc.type == VRNA_FC_TYPE_SINGLE:
                    type = vrna_get_ptype_md(S[i], S[j], md)
                    rt = vrna_get_ptype_md(S[j], S[i], md)
                else:
                    for s in range(n_seq):
                        tt[s] = vrna_get_ptype_md(SS[s][j], SS[s][i], md)

                tmp2 = vrna_exp_E_hp_loop(fc, j, i)

                if hard_constraints[i * n + j] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP:
                    for k in range(1, i - 1):
                        ln1 = k - 1
                        ln3 = n - j

                        if hc.up_int[j + 1] < (ln1 + ln3):
                            break

                        if (ln1 + ln3) > MAXLOOP:
                            break

                        lstart = (ln1 + ln3) + i - 1 - MAXLOOP
                        if lstart < k + 1:
                            lstart = k + 1

                        for l in range(lstart, i):
                            ln2 = i - l - 1

                            if (ln1 + ln2 + ln3) > MAXLOOP:
                                continue

                            if hc.up_int[l + 1] < ln2:
                                continue

                            if qb[my_iindx[k] - l] == 0:
                                continue

                            eval = 1 if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) else 0
                            if hc.f:
                                eval = hc.f(k, l, i, j, VRNA_DECOMP_PAIR_IL, hc.data)

                            if eval:
                                tmp = qb[my_iindx[k] - l]

                                if fc.type == VRNA_FC_TYPE_SINGLE:
                                    type_2 = vrna_get_ptype(jindx[l] + k, ptype)

                                    tmp *= exp_E_IntLoop(ln1 + ln3,
                                                         ln2,
                                                         rt,
                                                         rtype[type_2],
                                                         S1[j + 1],
                                                         S1[i - 1],
                                                         S1[k - 1],
                                                         S1[l + 1],
                                                         pf_params)
                                else:
                                    for s in range(n_seq):
                                        ln2a = a2s[s][i - 1] - a2s[s][l]
                                        ln1a = a2s[s][n] - a2s[s][j] + a2s[s][k - 1]
                                        type_2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md)
                                        tmp *= exp_E_IntLoop(ln1a, ln2a, tt[s], type_2,
                                                             S3[s][j],
                                                             S5[s][i],
                                                             S5[s][k],
                                                             S3[s][l],
                                                             pf_params)

                                if sc_dat_int.pair_ext:
                                    tmp *= sc_dat_int.pair_ext(k, l, i, j, sc_dat_int)

                                tmp2 += tmp * scale[ln1 + ln2 + ln3]

                    for k in range(j + 1, n):
                        ln1 = k - j - 1

                        if hc.up_int[j + 1] < ln1:
                            break

                        if (ln1 + i - 1) > MAXLOOP:
                            break

                        lstart = ln1 + i - 1 + n - MAXLOOP
                        if lstart < k + 1:
                            lstart = k + 1

                        for l in range(lstart, n + 1):
                            ln2 = i - 1
                            ln3 = n - l

                            if (ln1 + ln2 + ln3) > MAXLOOP:
                                continue

                            if hc.up_int[l + 1] < (ln2 + ln3):
                                continue

                            if qb[my_iindx[k] - l] == 0:
                                continue

                            eval = 1 if (hard_constraints[k * n + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) else 0
                            if hc.f:
                                eval = eval if hc.f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc.data) else 0

                            if eval:
                                tmp = qb[my_iindx[k] - l]

                                if fc.type == VRNA_FC_TYPE_SINGLE:
                                    type_2 = vrna_get_ptype(jindx[l] + k, ptype)
                                    tmp *= exp_E_IntLoop(ln2 + ln3,
                                                         ln1,
                                                         rtype[type_2],
                                                         rt,
                                                         S1[l + 1],
                                                         S1[k - 1],
                                                         S1[i - 1],
                                                         S1[j + 1],
                                                         pf_params)
                                else:
                                    for s in range(n_seq):
                                        ln1a = a2s[s][k] - a2s[s][j + 1]
                                        ln2a = a2s[s][i - 1] + a2s[s][n] - a2s[s][l]
                                        type_2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md)
                                        tmp *= exp_E_IntLoop(ln2a, ln1a, type_2, tt[s],
                                                             S3[s][l],
                                                             S5[s][k],
                                                             S5[s][i],
                                                             S3[s][j],
                                                             pf_params)

                                if sc_dat_int.pair_ext:
                                    tmp *= sc_dat_int.pair_ext(i, j, k, l, sc_dat_int)

                                tmp2 += tmp * scale[ln1 + ln2 + ln3]

                if hc_eval_mb(i, j, i - 1, j + 1, VRNA_DECOMP_PAIR_ML_EXT, hc_dat_mb):
                    sc_contrib = 1

                    if sc_dat_mb.pair_ext:
                        sc_contrib = sc_dat_mb.pair_ext(i, j, sc_dat_mb)

                    if (i > 2) and (j < n - 1):
                        tmp = qm[my_iindx[1] - i + 1] * qm[my_iindx[j + 1] - n]

                        if fc.type == VRNA_FC_TYPE_SINGLE:
                            tmp *= exp_E_MLstem(type,
                                                S1[i - 1],
                                                S1[j + 1],
                                                pf_params) * expMLclosing
                        else:
                            for s in range(n_seq):
                                tmp *= exp_E_MLstem(rtype[tt[s]],
                                                    S5[s][i],
                                                    S3[s][j],
                                                    pf_params)

                            tmp *= np.power(expMLclosing, n_seq)

                        tmp2 += tmp * sc_contrib
                    # 1.3.2 Left part
                    if hc.up_ml[j + 1] >= (n - j):
                        for k in range(2, i - 2):
                            if (hc_eval_mb(i, n, i, k, VRNA_DECOMP_PAIR_ML, hc_dat_mb)) and (hc_eval_mb(1, i - 1, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc_dat_mb)):
                                tmp = qm[my_iindx[1] - k] * qm1[jindx[k] + i - 1] * expMLbase[n - j]

                                if fc.type == VRNA_FC_TYPE_SINGLE:
                                    tmp *= exp_E_MLstem(type,
                                                        S1[i - 1],
                                                        S1[j + 1],
                                                        pf_params) * expMLclosing
                                else:
                                    for s in range(n_seq):
                                        tmp *= exp_E_MLstem(rtype[tt[s]],
                                                            S5[s][i],
                                                            S3[s][j],
                                                            pf_params)

                                    tmp *= np.power(expMLclosing, n_seq)
                                if (sc_dat_mb.red_ml):
                                    tmp *= sc_dat_mb.red_ml(i, n, i, j, sc_dat_mb)
                                if (sc_dat_mb.decomp_ml):
                                    tmp *= sc_dat_mb.decomp_ml(1, i - 1, k, k + 1, sc_dat_mb)
                                tmp2 += tmp * sc_contrib
                    # 1.3.3 Right part
                    if hc.up_ml[1] >= (i - 1):
                        for k in range(j + 2, n):
                            if (hc_eval_mb(1, j, k, n, VRNA_DECOMP_PAIR_ML, hc_dat_mb)) and (hc_eval_mb(j + 1, n, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc_dat_mb)):
                                tmp = qm[my_iindx[k] - n] * qm1[jindx[1] + j] * expMLbase[i - 1]

                                if fc.type == VRNA_FC_TYPE_SINGLE:
                                    tmp *= exp_E_MLstem(type,
                                                        S1[i - 1],
                                                        S1[j + 1],
                                                        pf_params) * expMLclosing
                                else:
                                    for s in range(n_seq):
                                        tmp *= exp_E_MLstem(rtype[tt[s]],
                                                            S5[s][i],
                                                            S3[s][j],
                                                            pf_params)

                                    tmp *= np.power(expMLclosing, n_seq)

                                if (sc_dat_mb.red_ml):
                                    tmp *= sc_dat_mb.red_ml(1, j, i, j, sc_dat_mb)
                                if (sc_dat_mb.decomp_ml):
                                    tmp *= sc_dat_mb.decomp_ml(j + 1, n, k, k + 1, sc_dat_mb)
                                tmp2 += tmp * sc_contrib
                probs[ij] *= tmp2
            else:
                probs[ij] = 0
                    

# compute_bpp_external constant and class start 




# compute_bpp_external method start 
def vrna_exp_E_ext_stem(type:int, n5d:int, n3d:int, p:vrna_exp_param_t) -> float:
    energy = 1.0
    if n5d >= 0 and n3d >= 0:
        energy = p.expmismatchExt[type][n5d][n3d]
    elif (n5d >= 0):
        energy = p.expdangle5[type][n5d]
    elif (n3d >= 0):
        energy = p.expdangle3[type][n3d]

    if (type > 2):
        energy *= p.expTermAU

    return energy

def contrib_ext_pair(fc: vrna_fold_compound_t,
                     i: int,
                     j: int,
                     constraints: constraints_helper) -> float:
    type: int
    ptype: str
    S1: List[int]
    s5: int
    s3: int
    sn: List[int]
    n: int
    jindx: List[int]
    contribution: float
    pf_params: vrna_exp_param_t
    sc: vrna_sc_t

    n = fc.length
    pf_params = fc.exp_params
    S1 = fc.sequence_encoding
    sn = fc.strand_number
    ptype = fc.ptype
    jindx = fc.jindx
    sc = fc.sc

    type = vrna_get_ptype(jindx[j] + i, ptype)
    s5 = S1[i - 1] if i > 1 and sn[i] == sn[i - 1] else -1
    s3 = S1[j + 1] if j < n and sn[j + 1] == sn[j] else -1

    contribution = vrna_exp_E_ext_stem(type, s5, s3, pf_params)

    if sc and sc.exp_f:
        contribution *= sc.exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, sc.data)

    return contribution


def contrib_ext_pair_comparative(fc: vrna_fold_compound_t, i: int, j: int, constraints: constraints_helper) -> float:
    type: int
    S: List[List[int]]
    S5: List[List[int]]
    S3: List[List[int]]
    s5: int
    s3: int
    a2s: List[List[int]]
    n: int
    s: int
    n_seq: int
    jindx: List[int]
    pscore: List[int]
    contribution: float
    kTn: float
    pf_params: vrna_exp_param_t
    md: vrna_md_t
    scs: List[vrna_sc_t]

    n = fc.length
    n_seq = fc.n_seq
    jindx = fc.jindx
    pf_params = fc.exp_params
    md = pf_params.model_details
    S = fc.S
    S5 = fc.S5
    S3 = fc.S3
    a2s = fc.a2s
    pscore = fc.pscore  # precomputed array of pair types
    scs = fc.scs
    kTn = pf_params.kT / 10.0  # kT in cal/mol

    contribution = math.exp(pscore[jindx[j] + i] / kTn)

    for s in range(n_seq):
        type = vrna_get_ptype_md(S[s][i], S[s][j], md)
        s5 = S5[s][i] if a2s[s][i] > 1 else -1
        s3 = S3[s][j] if a2s[s][j] < a2s[s][n] else -1

        contribution *= vrna_exp_E_ext_stem(type, s5, s3, pf_params)

    if scs:
        for s in range(n_seq):
            if scs[s].exp_f:
                contribution *= scs[s].exp_f(1, n, i, j, VRNA_DECOMP_EXT_STEM_OUTSIDE, scs[s].data)

    return contribution




# compute_bpp_external start 
def compute_bpp_external(fc:vrna_fold_compound_t, constraints:constraints_helper):
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


def multistrand_update_Y5(fc:vrna_fold_compound_t, l, Y5, Y5p, constraints:constraints_helper):
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


def multistrand_update_Y3(fc:vrna_fold_compound_t, l, Y3, Y3p, constraints:constraints_helper):
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


def multistrand_contrib(fc:vrna_fold_compound_t, l, Y5, Y3, constraints:constraints_helper, Qmax, ov):
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


# ud_outside_ext_loops method start
from typing import Optional, List

def vrna_nucleotide_IUPAC_identity(nt: str, mask: str) -> int:
    n1: str
    n2: str
    p: Optional[str]

    p = None
    n1 = nt.upper()
    n2 = mask.upper()

    if n1 == 'A':
        p = n2 in "ARMWDHVN"
    elif n1 == 'C':
        p = n2 in "CYMSBHVN"
    elif n1 == 'G':
        p = n2 in "GRKSBDVN"
    elif n1 == 'T':
        p = n2 in "TYKWBDHN"
    elif n1 == 'U':
        p = n2 in "UYKWBDHN"
    elif n1 == 'I':
        p = n2 in "IN"
    elif n1 == 'R':
        p = n2 in "AGR"
    elif n1 == 'Y':
        p = n2 in "CTUY"
    elif n1 == 'K':
        p = n2 in "GTUK"
    elif n1 == 'M':
        p = n2 in "ACM"
    elif n1 == 'S':
        p = n2 in "GCS"
    elif n1 == 'W':
        p = n2 in "ATUW"
    elif n1 == 'B':
        p = n2 in "GCTBU"
    elif n1 == 'D':
        p = n2 in "AGTUD"
    elif n1 == 'H':
        p = n2 in "ACTUH"
    elif n1 == 'V':
        p = n2 in "ACGV"
    elif n1 == 'N':
        p = n2 in "ACGTUN"

    return 1 if p else 0


def get_motifs(vc: vrna_fold_compound_t, i: int, loop_type: int) -> Optional[List[int]]:
    k: int
    j: int
    u: int
    n: int
    cnt: int
    guess: int
    sequence: str
    domains_up: vrna_ud_t
    motif_list: Optional[List[int]]

    sequence = vc.sequence
    n = int(vc.length)
    domains_up = vc.domains_up

    cnt = 0
    guess = domains_up.motif_count
    motif_list = [0] * (guess + 1)

    # collect list of motif numbers we find that start at position i
    for k in range(domains_up.motif_count):
        if not (domains_up.motif_type[k] & loop_type):
            continue
        
        j = i + domains_up.motif_size[k] - 1
        if j <= n:
            # only consider motif that does not exceed sequence length (does not work for circular RNAs!)
            for u in range(i, j + 1):
                if not vrna_nucleotide_IUPAC_identity(sequence[u - 1], domains_up.motif[k][u - i]):
                    break
            
            if u > j:  # got a complete motif match
                motif_list[cnt] = k
                cnt += 1

    if cnt == 0:
        return None

    motif_list = motif_list[:cnt]
    motif_list.append(-1)  # end of list marker

    return motif_list


def vrna_ud_get_motif_size_at(vc: vrna_fold_compound_t, i: int, loop_type: int):
    if vc and vc.domains_up:
        k: int
        l: int
        cnt: int
        ret: Optional[List[int]]
        ptr: Optional[List[int]]

        ret = None
        if i > 0 and i <= vc.length:
            ptr = get_motifs(vc, i, loop_type)
            if ptr:
                k = 0
                while ptr[k] != -1:
                    ptr[k] = vc.domains_up.motif_size[ptr[k]]
                    k += 1

                # make the list unique
                ret = [-1] * (k + 1)
                cnt = 0
                k = 0
                while ptr[k] != -1:
                    l = 0
                    while l < cnt:
                        if ptr[k] == ret[l]:
                            break
                        l += 1

                    if l == cnt:
                        ret[cnt] = ptr[k]
                        ret[cnt + 1] = -1
                        cnt += 1
                    k += 1

                # resize ret array
                ret = ret[:cnt + 1]


        return ret

    return None



# ud_outside_ext_loops start
def ud_outside_ext_loops(vc:vrna_fold_compound_t):
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

          
          
# ud_outside_hp_loops start


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


# ud_outside_int_loops start


def vrna_get_ptype_window(i: int, j: int, ptype: List[str]) -> int:
    tt = ord(ptype[i][j - i])

    return 7 if tt == 0 else tt


eval_hc = Callable[[int, int, int, int, None], int]

from typing import Optional, List

def exp_E_interior_loop(fc: vrna_fold_compound_t,
                        i: int,
                        j: int,
                        k: int,
                        l: int) -> float:
    sliding_window: int
    type_: int
    type2: int
    ptype: Optional[str]
    ptype_local: Optional[List[str]]
    hc_mx: Optional[List[int]]
    hc_mx_local: Optional[List[List[int]]]
    eval_loop: int
    hc_decompose_ij: int
    hc_decompose_kl: int
    qbt1: float
    q_temp: float
    scale: List[float]
    pf_params: vrna_exp_param_t
    md: vrna_md_t
    domains_up: Optional[vrna_ud_t]
    evaluate: eval_hc
    hc_dat_local: hc_int_def_dat
    sc_wrapper: sc_int_exp_dat

    sliding_window = 1 if fc.hc.type == VRNA_HC_WINDOW else 0
    n = fc.length
    n_seq = 1 if fc.type == VRNA_FC_TYPE_SINGLE else fc.n_seq
    ptype = None if fc.type == VRNA_FC_TYPE_SINGLE else (fc.ptype if not sliding_window else None)
    ptype_local = None if fc.type == VRNA_FC_TYPE_SINGLE else (fc.ptype_local if sliding_window else None)
    S1 = fc.sequence_encoding if fc.type == VRNA_FC_TYPE_SINGLE else None
    SS = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.S
    S5 = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.S5
    S3 = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.S3
    a2s = None if fc.type == VRNA_FC_TYPE_SINGLE else fc.a2s
    jindx = fc.jindx
    hc_mx = None if sliding_window else fc.hc.mx
    hc_mx_local = fc.hc.matrix_local if sliding_window else None
    hc_up = fc.hc.up_int
    pf_params = fc.exp_params
    sn = fc.strand_number
    md = pf_params.model_details
    scale = fc.exp_matrices.scale
    domains_up = fc.domains_up
    rtype = md.rtype
    qbt1 = 0.0
    u1 = k - i - 1
    u2 = j - l - 1

    if sn[k] != sn[i] or sn[j] != sn[l]:
        return qbt1

    if hc_up[l + 1] < u2:
        return qbt1

    if hc_up[i + 1] < u1:
        return qbt1

    evaluate = prepare_hc_int_def(fc, hc_dat_local)
    init_sc_int_exp(fc, sc_wrapper)

    hc_decompose_ij = hc_mx_local[i][j - i] if sliding_window else hc_mx[n * i + j]
    hc_decompose_kl = hc_mx_local[k][l - k] if sliding_window else hc_mx[n * k + l]
    eval_loop = 1 if (hc_decompose_ij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) and (hc_decompose_kl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) else 0

    if eval_loop and evaluate(i, j, k, l, hc_dat_local):
        q_temp = 0.0

        if fc.type == VRNA_FC_TYPE_SINGLE:
            type_ = vrna_get_ptype_window(i, j, ptype_local) if sliding_window else vrna_get_ptype(jindx[j] + i, ptype)
            type2 = rtype[vrna_get_ptype_window(k, l, ptype_local)] if sliding_window else rtype[vrna_get_ptype(jindx[l] + k, ptype)]

            q_temp = exp_E_IntLoop(u1,
                                   u2,
                                   type_,
                                   type2,
                                   S1[i + 1],
                                   S1[j - 1],
                                   S1[k - 1],
                                   S1[l + 1],
                                   pf_params)
        elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
            q_temp = 1.0

            for s in range(n_seq):
                u1_local = a2s[s][k - 1] - a2s[s][i]
                u2_local = a2s[s][j - 1] - a2s[s][l]
                type_ = vrna_get_ptype_md(SS[s][i], SS[s][j], md)
                type2 = vrna_get_ptype_md(SS[s][l], SS[s][k], md)

                q_temp *= exp_E_IntLoop(u1_local,
                                        u2_local,
                                        type_,
                                        type2,
                                        S3[s][i],
                                        S5[s][j],
                                        S5[s][k],
                                        S3[s][l],
                                        pf_params)

        if sc_wrapper.pair:
            q_temp *= sc_wrapper.pair(i, j, k, l, sc_wrapper)

        qbt1 += q_temp * scale[u1 + u2 + 2]

        if domains_up and domains_up.exp_energy_cb:
            qq5 = qq3 = 0.0

            if u1 > 0:
                qq5 = domains_up.exp_energy_cb(fc,
                                               i + 1, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                               domains_up.data)

            if u2 > 0:
                qq3 = domains_up.exp_energy_cb(fc,
                                               l + 1, j - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                               domains_up.data)

            qbt1 += q_temp * qq5 * scale[u1 + u2 + 2]
            qbt1 += q_temp * qq3 * scale[u1 + u2 + 2]
            qbt1 += q_temp * qq5 * qq3 * scale[u1 + u2 + 2]


    return qbt1

def vrna_exp_E_interior_loop(fc: vrna_fold_compound_t,
                             i: int,
                             j: int,
                             k: int,
                             l: int) -> float:
    if fc:
        return exp_E_interior_loop(fc, i, j, k, l)
    return 0.0


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


# vrna_exp_E_hp_loop start ############
# Python 中模拟 C 函数 hc_hp_cb_def_window 的功能
def hc_hp_cb_def_window(i: int, j: int, k: int, l: int, d: int, data: hc_hp_def_dat) -> int:
    dat = data  # 假设 data 是 HcHpDefDat 类的实例

    eval_result = 0
    u = j - i - 1

    # 检查 mx_window 和 hc_up 的条件
    if dat.mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP:
        eval_result = 1
        if dat.hc_up[i + 1] < u:
            eval_result = 0

    return eval_result


def hc_hp_cb_def_user_window(i: int, j: int, k: int, l: int, d: int, data: hc_hp_def_dat) -> int:
    # 数据类型的转换
    dat = data  # 假设 data 是 struct hc_hp_def_dat *

    # 调用 hc_hp_cb_def_window
    eval_result = hc_hp_cb_def_window(i, j, k, l, d, data)
    # 检查 hc_f 的返回值
    eval_result = eval_result if dat.hc_f(i, j, k, l, d, dat.hc_dat) else 0

    return eval_result
# 模拟私有函数 prepare_hc_hp_def_window
def prepare_hc_hp_def_window(fc: vrna_fold_compound_t, dat: hc_hp_def_dat):
    # 设置 dat 的属性
    dat.mx_window = fc.hc.matrix_local
    dat.hc_up = fc.hc.up_hp
    dat.n = fc.length
    dat.sn = fc.strand_number

    # 检查是否有 fc.hc.f 函数
    if fc.hc.f:
        dat.hc_f = fc.hc.f
        dat.hc_dat = fc.hc.data
        return hc_hp_cb_def_user_window  # 返回模拟函数指针

    return hc_hp_cb_def_window  # 返回另一个模拟函数指针


def exp_E_Hairpin(u: int, type: int, si1: int, sj1: int, string: str, P:vrna_exp_param_t) -> float:
    q = 0.0
    kT = P.kT  # kT in cal/mol
    salt_correction = 1.0

    if P.model_details.salt != VRNA_MODEL_DEFAULT_SALT:
        if u <= MAXLOOP:
            salt_correction = P.expSaltLoop[u + 1]
        else:
            salt_correction = math.exp(-vrna_salt_loop_int(u + 1, P.model_details.salt, P.temperature + K0, P.model_details.backbone_length) * 10.0 / kT)

    if u <= 30:
        q = P.exphairpin[u]
    else:
        q = P.exphairpin[30] * math.exp(-(P.lxc * math.log(u / 30.0)) * 10.0 / kT)

    q *= salt_correction

    if u < 3:
        return q  # should only be the case when folding alignments

    if string and P.model_details.special_hp:
        if u == 4:
            tl = string[:6]
            ts = P.Tetraloops.find(tl)
            if ts != -1:
                if type != 7:
                    return P.exptetra[ts // 7] * salt_correction
                else:
                    q *= P.exptetra[ts // 7]
        elif u == 6:
            tl = string[:8]
            ts = P.Hexaloops.find(tl)
            if ts != -1:
                return P.exphex[ts // 9] * salt_correction
        elif u == 3:
            tl = string[:5]
            ts = P.Triloops.find(tl)
            if ts != -1:
                return P.exptri[ts // 6] * salt_correction

            if type > 2:
                return q * P.expTermAU
            else:
                return q

    q *= P.expmismatchH[type][si1][sj1]

    return q

# Python 中模拟 C 函数 exp_eval_hp_loop 的功能
def exp_eval_hp_loop(fc: vrna_fold_compound_t, i: int, j: int) -> float:
    P = fc.exp_params
    md = P.model_details
    sn = fc.strand_number
    scale = fc.exp_matrices.scale
    domains_up = fc.domains_up

    sc_wrapper = sc_hp_exp_dat()
    init_sc_hp_exp(fc, sc_wrapper)

    q = 0.0

    if sn[j] != sn[i]:
        return q

    if fc.type == VRNA_FC_TYPE_SINGLE:
        S = fc.sequence_encoding
        S2 = fc.sequence_encoding2
        u = j - i - 1
        type = vrna_get_ptype_md(S2[i], S2[j], md)

        if sn[j] == sn[i]:
            q = exp_E_Hairpin(u, type, S[i + 1], S[j - 1], fc.sequence[i - 1:], P)
        else:
            # 处理外部 hairpin loop 的情况
            pass

    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        SS = fc.S
        S5 = fc.S5
        S3 = fc.S3
        Ss = fc.Ss
        a2s = fc.a2s
        n_seq = fc.n_seq
        qbt1 = 1.0

        for s in range(n_seq):
            u = a2s[s][j - 1] - a2s[s][i]
            if a2s[s][i] < 1:
                continue

            type = vrna_get_ptype_md(SS[s][i], SS[s][j], md)
            qbt1 *= exp_E_Hairpin(u, type, S3[s][i], S5[s][j], Ss[s][a2s[s][i] - 1:], P)

        q = qbt1

    # 添加软约束
    if sc_wrapper.pair:
        q *= sc_wrapper.pair(i, j, sc_wrapper)

    if domains_up and domains_up.exp_energy_cb:
        q += q * domains_up.exp_energy_cb(fc, i + 1, j - 1, VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP, domains_up.data)

    q *= scale[j - i + 1]


    return q


# Python 中模拟 C 函数 exp_eval_ext_hp_loop 的功能
def exp_eval_ext_hp_loop(fc: vrna_fold_compound_t, i: int, j: int) -> float:
    n = fc.length
    P = fc.exp_params
    md = P.model_details
    noGUclosure = md.noGUclosure
    scale = fc.exp_matrices.scale
    domains_up = fc.domains_up

    sc_wrapper = sc_hp_exp_dat()
    init_sc_hp_exp(fc, sc_wrapper)

    q = 0.0
    u1 = n - j
    u2 = i - 1

    if (u1 + u2) < 3:
        return q

    if fc.type == VRNA_FC_TYPE_SINGLE:
        sequence = fc.sequence
        S = fc.sequence_encoding
        S2 = fc.sequence_encoding2
        type = vrna_get_ptype_md(S2[j], S2[i], md)

        if (type == 3 or type == 4) and noGUclosure:
            return q

        loopseq = sequence[j - 1:] + sequence[:i]
        loopseq = loopseq[:u1 + u2 + 2]

        q = exp_E_Hairpin(u1 + u2, type, S[j + 1], S[i - 1], loopseq, P)

    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        SS = fc.S
        S5 = fc.S5
        S3 = fc.S3
        Ss = fc.Ss
        a2s = fc.a2s
        n_seq = fc.n_seq
        qbt1 = 1.0

        for s in range(n_seq):
            u1_local = a2s[s][n] - a2s[s][j]
            u2_local = a2s[s][i - 1]
            loopseq = Ss[s][a2s[s][j] - 1:] + Ss[s][:a2s[s][i]]
            loopseq = loopseq[:u1_local + u2_local + 2]

            type = vrna_get_ptype_md(SS[s][j], SS[s][i], md)
            qbt1 *= exp_E_Hairpin(u1_local + u2_local, type, S3[s][j], S5[s][i], loopseq, P)

        q = qbt1

    # 添加软约束
    if sc_wrapper.pair_ext:
        q *= sc_wrapper.pair_ext(i, j, sc_wrapper)

    if domains_up and domains_up.exp_energy_cb:
        q += q * domains_up.exp_energy_cb(fc, j + 1, i - 1, VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP, domains_up.data)

    q *= scale[u1 + u2]


    return q





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



            






def pf_create_bppm(vc:vrna_fold_compound_t, structure):
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

    print(f"qb is {qb}")
    print(f"probs is {probs}")
    print(f"circular is {circular}")
    print(f"q1k is {q1k}")
    print(f"qln is {qln}")
    
    # not into 
    if qb and probs and (circular or (q1k is not None and qln is not None)):
        with_gquad = pf_params.model_details.gquad
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

        # init diagonal entries unable to pair in pr matrix
        probs[my_iindx[1] - 1:my_iindx[n] - n + 1] = 0

        if circular:
            bppm_circ(vc, constraints)
        else:
            compute_bpp_external(vc, constraints)

        l = n
        compute_bpp_int(vc, l, bp_correction, corr_cnt, corr_size, Qmax, ov, constraints)

        for l in range(n - 1, 1, -1):
            compute_bpp_int(vc, l, bp_correction, corr_cnt, corr_size, Qmax, ov, constraints)
            compute_bpp_mul(vc, l, ml_helpers, Qmax, ov, constraints)

            if vc.strands > 1:
                multistrand_update_Y5(vc, l, Y5, Y5p, constraints)
                multistrand_update_Y3(vc, l, Y3, Y3p, constraints)
                multistrand_contrib(vc, l, Y5, Y3, constraints, Qmax, ov)

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
                                    if ptr and (ptr.i != 0):
                                        bp_correction[corr_cnt].i   = ptr.i
                                        bp_correction[corr_cnt].j   = ptr.j
                                        bp_correction[corr_cnt].p = probs[ij] * qhp
                                        if corr_cnt == corr_size:
                                            corr_size += 5
                                            bp_correction = np.resize(bp_correction, corr_size)
                                        # bp_correction[corr_cnt] = (ptr.i, ptr.j, probs[ij] * qhp)
                                        corr_cnt += 1

                for i in range(corr_cnt):
                    ij = my_iindx[bp_correction[i].i] - bp_correction[i].j
                    probs[ij] += bp_correction[i].p / qb[ij]

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

# vrna_fold_compound ####################
def nullify(fc: vrna_fold_compound_t):
    if fc:
        fc.length            = 0
        fc.strands           = 0
        fc.cutpoint          = -1
        fc.strand_number     = None
        fc.strand_order      = None
        fc.strand_order_uniq = None
        fc.strand_start      = None
        fc.strand_end        = None
        fc.nucleotides       = None
        fc.alignment         = None

        fc.hc            = None
        fc.matrices      = None
        fc.exp_matrices  = None
        fc.params        = None
        fc.exp_params    = None
        fc.iindx         = None
        fc.jindx         = None

        fc.stat_cb       = None
        fc.auxdata       = None
        fc.free_auxdata  = None

        fc.domains_struc = None
        fc.domains_up    = None
        fc.aux_grammar   = None

        if fc.type == VRNA_FC_TYPE_SINGLE:
            fc.sequence            = None
            fc.sequence_encoding   = None
            fc.encoding5           = None
            fc.encoding3           = None
            fc.sequence_encoding2  = None
            fc.ptype               = None
            fc.ptype_pf_compat     = None
            fc.sc                  = None

        elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
            fc.sequences         = None
            fc.n_seq             = 0
            fc.cons_seq          = None
            fc.S_cons            = None
            fc.S                 = None
            fc.S5                = None
            fc.S3                = None
            fc.Ss                = None
            fc.a2s               = None
            fc.pscore            = None
            fc.pscore_local      = None
            fc.pscore_pf_compat  = None
            fc.scs               = None
            fc.oldAliEn          = 0

        fc.maxD1         = 0
        fc.maxD2         = 0
        fc.reference_pt1 = None
        fc.reference_pt2 = None
        fc.referenceBPs1 = None
        fc.referenceBPs2 = None
        fc.bpdist        = None
        fc.mm1           = None
        fc.mm2           = None

        fc.window_size = -1
        fc.ptype_local = None
        # Uncomment the following line if VRNA_WITH_SVM is defined in your context
        # fc.zscore_data = None

def init_fc_single():
    init = vrna_fold_compound_t()
    init.type = VRNA_FC_TYPE_SINGLE
    fc = vrna_fold_compound_t()
    if fc:
        fc.type = init.type
        nullify(fc)
    return fc

def vrna_md_copy(md_to: vrna_md_t, md_from: vrna_md_t) -> vrna_md_t:
    if md_from is None:
        return None

    if md_to is None:
        md_to = vrna_md_t()

    if md_to != md_from:
        md_to.rtype = md_from.rtype.copy()
        md_to.alias = md_from.alias.copy()
        md_to.nonstandards = md_from.nonstandards.copy()
        md_to.pair = [row.copy() for row in md_from.pair]
        md_to.pair_dist = [row.copy() for row in md_from.pair_dist]
        # Copy other fields as necessary

    return md_to

def vrna_md_set_default(md:vrna_md_t):
            if md:
                defaults = vrna_md_t()
                vrna_md_copy(md, defaults)
    
    
from .remove import Duplex, Tmeasure, expn



 
def vrna_params(md=None):
    if md:
        return Duplex.get_scaled_params(md)
    else:
        md = vrna_md_t()
        vrna_md_set_default(md)
        return Duplex.get_scaled_params(md)  


def get_scaled_exp_params(md, pfs):
    pf = vrna_exp_param_t()

    if last_parameter_file() is not None:
        pf.param_file = last_parameter_file()

    pf.model_details = md
    pf.temperature = md.temperature
    pf.alpha = md.betaScale
    pf.kT = kT = md.betaScale * (md.temperature + K0) * GASCONST  # kT in cal/mol
    pf.pf_scale = pfs
    pf_smooth = md.pf_smooth
    TT = (md.temperature + K0) / Tmeasure
    salt = md.salt
    saltStandard = VRNA_MODEL_DEFAULT_SALT

    pf.lxc = lxc37 * TT
    pf.expDuplexInit = RESCALE_BF(DuplexInit37, DuplexInitdH, TT, kT)
    pf.expTermAU = RESCALE_BF(TerminalAU37, TerminalAUdH, TT, kT)
    pf.expMLbase = RESCALE_BF(ML_BASE37, ML_BASEdH, TT, kT)
    pf.expMLclosing = RESCALE_BF(ML_closing37, ML_closingdH, TT, kT)
    pf.expgquadLayerMismatch = RESCALE_BF(GQuadLayerMismatch37, GQuadLayerMismatchH, TT, kT)
    pf.gquadLayerMismatchMax = GQuadLayerMismatchMax

    for i in range(VRNA_GQUAD_MIN_STACK_SIZE, VRNA_GQUAD_MAX_STACK_SIZE + 1):
        for j in range(3 * VRNA_GQUAD_MIN_LINKER_LENGTH, 3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1):
            GQuadAlpha_T = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, TT)
            GQuadBeta_T = RESCALE_dG(GQuadBeta37, GQuadBetadH, TT)
            GT = GQuadAlpha_T * (i - 1) + GQuadBeta_T * math.log(j - 2)
            pf.expgquad[i][j] = math.exp(-TRUNC_MAYBE(GT) * 10.0 / kT)

    for i in range(31):
        pf.exphairpin[i] = RESCALE_BF(hairpin37[i], hairpindH[i], TT, kT)

    for i in range(MIN2(30, MAXLOOP) + 1):
        pf.expbulge[i] = RESCALE_BF(bulge37[i], bulgedH[i], TT, kT)
        pf.expinternal[i] = RESCALE_BF(internal_loop37[i], internal_loopdH[i], TT, kT)

    if salt == saltStandard:
        for i in range(MIN2(30, MAXLOOP) + 1):
            pf.SaltLoopDbl[i] = 0.0
            pf.expSaltLoop[i] = 1.0
        for i in range(31, MAXLOOP + 1):
            pf.SaltLoopDbl[i] = 0.0
            pf.expSaltLoop[i] = 1.0
    else:
        for i in range(MIN2(30, MAXLOOP) + 1):
            pf.SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md.backbone_length)
            saltLoop = int(pf.SaltLoopDbl[i] + 0.5 - (pf.SaltLoopDbl[i] < 0))
            pf.expSaltLoop[i] = math.exp(-saltLoop * 10.0 / kT)
        for i in range(31, MAXLOOP + 1):
            pf.SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md.backbone_length)
            saltLoop = int(pf.SaltLoopDbl[i] + 0.5 - (pf.SaltLoopDbl[i] < 0))
            pf.expSaltLoop[i] = math.exp(-saltLoop * 10.0 / kT)

    if james_rule:
        pf.expinternal[2] = math.exp(-80 * 10.0 / kT)

    GT = RESCALE_dG(bulge37[30], bulgedH[30], TT)
    for i in range(31, MAXLOOP + 1):
        pf.expbulge[i] = math.exp(-TRUNC_MAYBE(GT + (pf.lxc * math.log(i / 30.0))) * 10.0 / kT)

    GT = RESCALE_dG(internal_loop37[30], internal_loopdH[30], TT)
    for i in range(31, MAXLOOP + 1):
        pf.expinternal[i] = math.exp(-TRUNC_MAYBE(GT + (pf.lxc * math.log(i / 30.0))) * 10.0 / kT)

    GT = RESCALE_dG(ninio37, niniodH, TT)
    for j in range(MAXLOOP + 1):
        pf.expninio[2][j] = math.exp(-MIN2(MAX_NINIO, j * TRUNC_MAYBE(GT)) * 10.0 / kT)

    for i in range(0, len(Tetraloops) // 7):
        pf.exptetra[i] = RESCALE_BF(Tetraloop37[i], TetraloopdH[i], TT, kT)

    for i in range(0, len(Triloops) // 5):
        pf.exptri[i] = RESCALE_BF(Triloop37[i], TriloopdH[i], TT, kT)

    for i in range(0, len(Hexaloops) // 9):
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
                    pf.expmismatchM[i][j][k] = pf.expmismatchExt[i][j][k] = 1.0

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

    pf.Tetraloops = Tetraloops
    pf.Triloops = Triloops
    pf.Hexaloops = Hexaloops

    pf.SaltMLbase = pf.SaltMLclosing = pf.SaltDPXInit = 0.0

    if salt == saltStandard:
        pf.expSaltStack = 1.0
    else:
        pf.expSaltStack = math.exp(-vrna_salt_stack(salt, saltT, md.helical_rise) * 10.0 / kT)
        vrna_salt_ml(pf.SaltLoopDbl, md.saltMLLower, md.saltMLUpper, pf.SaltMLbase, pf.SaltMLclosing)
        if md.saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT:
            pf.SaltDPXInit = md.saltDPXInit
        elif md.saltDPXInit:
            pf.SaltDPXInit = vrna_salt_duplex_init(md)
        pf.expMLclosing *= math.exp(-pf.SaltMLbase * 10.0 / kT)
        pf.expMLclosing *= math.exp(-pf.SaltMLclosing * 10.0 / kT)
        pf.expMLbase *= math.exp(-pf.SaltMLbase * 10.0 / kT)
        for i in range(NBPAIRS + 1):
            pf.expMLintern[i] *= math.exp(-pf.SaltMLbase * 10.0 / kT)
        pf.expDuplexInit *= math.exp(-pf.SaltDPXInit * 10.0 / kT)

    return pf

   

def vrna_exp_params(md):
    if md:
        return get_scaled_exp_params(md, -1.0)
    else:
        md = vrna_md_t()
        vrna_md_set_default(md)
        return get_scaled_exp_params(md, -1.0)



def get_exp_params_ali(md: vrna_md_t, n_seq: int, pfs: float) -> vrna_exp_param_t:
    pf = vrna_exp_param_t()
    pf.model_details = md
    pf.alpha = md.betaScale
    pf.temperature = md.temperature
    pf.pf_scale = pfs
    pf.kT = kTn = float(n_seq) * md.betaScale * (md.temperature + K0) * GASCONST  # kT in cal/mol
    pf_smooth = md.pf_smooth
    TT = (md.temperature + K0) / Tmeasure
    salt = md.salt
    saltStandard = VRNA_MODEL_DEFAULT_SALT

    pf.lxc = lxc37 * TT
    pf.expDuplexInit = RESCALE_BF(DuplexInit37, DuplexInitdH, TT, kTn)
    pf.expTermAU = RESCALE_BF(TerminalAU37, TerminalAUdH, TT, kTn)
    pf.expMLbase = RESCALE_BF(ML_BASE37, ML_BASEdH, TT, kTn / n_seq)
    pf.expMLclosing = RESCALE_BF(ML_closing37, ML_closingdH, TT, kTn)
    pf.expgquadLayerMismatch = RESCALE_BF(GQuadLayerMismatch37, GQuadLayerMismatchH, TT, kTn)
    pf.gquadLayerMismatchMax = GQuadLayerMismatchMax

    for i in range(VRNA_GQUAD_MIN_STACK_SIZE, VRNA_GQUAD_MAX_STACK_SIZE + 1):
        for j in range(3 * VRNA_GQUAD_MIN_LINKER_LENGTH, 3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1):
            GQuadAlpha_T = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, TT)
            GQuadBeta_T = RESCALE_dG(GQuadBeta37, GQuadBetadH, TT)
            GT = GQuadAlpha_T * (i - 1) + GQuadBeta_T * math.log(j - 2)
            pf.expgquad[i][j] = math.exp(-TRUNC_MAYBE(GT) * 10.0 / kTn)

    for i in range(31):
        pf.exphairpin[i] = RESCALE_BF(hairpin37[i], hairpindH[i], TT, kTn)
    for i in range(3):
        GT = 600 * TT  # Penalty
        pf.exphairpin[i] = math.exp(-GT * 10.0 / kTn)

    for i in range(min(31, MAXLOOP) + 1):
        pf.expbulge[i] = RESCALE_BF(bulge37[i], bulgedH[i], TT, kTn)
        pf.expinternal[i] = RESCALE_BF(internal_loop37[i], internal_loopdH[i], TT, kTn)
        pf.SaltLoopDbl[i] = 0.0 if salt == saltStandard else vrna_salt_loop(i, salt, saltT, md.backbone_length)
        saltLoop = int(pf.SaltLoopDbl[i] + 0.5 - (pf.SaltLoopDbl[i] < 0))
        pf.expSaltLoop[i] = math.exp(-saltLoop * 10.0 / kTn)

    if salt == saltStandard:
        for i in range(min(31, MAXLOOP) + 1):
            pf.SaltLoopDbl[i] = 0.0
            pf.expSaltLoop[i] = 1.0
        for i in range(31, MAXLOOP + 1):
            pf.SaltLoopDbl[i] = 0.0
            pf.expSaltLoop[i] = 1.0
    else:
        for i in range(min(31, MAXLOOP) + 1):
            pf.SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md.backbone_length)
            saltLoop = int(pf.SaltLoopDbl[i] + 0.5 - (pf.SaltLoopDbl[i] < 0))
            pf.expSaltLoop[i] = math.exp(-saltLoop * 10.0 / kTn)
        for i in range(31, MAXLOOP + 1):
            pf.SaltLoopDbl[i] = vrna_salt_loop(i, salt, saltT, md.backbone_length)
            saltLoop = int(pf.SaltLoopDbl[i] + 0.5 - (pf.SaltLoopDbl[i] < 0))
            pf.expSaltLoop[i] = math.exp(-saltLoop * 10.0 / kTn)

    if james_rule:
        pf.expinternal[2] = math.exp(-80 * 10.0 / kTn)

    GT = RESCALE_dG(bulge37[30], bulgedH[30], TT)
    for i in range(31, MAXLOOP + 1):
        pf.expbulge[i] = math.exp(-(GT + (pf.lxc * math.log(i / 30.0))) * 10.0 / kTn)

    GT = RESCALE_dG(internal_loop37[30], internal_loopdH[30], TT)
    for i in range(31, MAXLOOP + 1):
        pf.expinternal[i] = math.exp(-(GT + (pf.lxc * math.log(i / 30.0))) * 10.0 / kTn)

    GT = RESCALE_dG(ninio37, niniodH, TT)
    for j in range(MAXLOOP + 1):
        pf.expninio[2][j] = math.exp(-min(MAX_NINIO, j * GT) * 10.0 / kTn)

    for i in range(0, len(Tetraloops) // 7):
        pf.exptetra[i] = RESCALE_BF(Tetraloop37[i], TetraloopdH[i], TT, kTn)

    for i in range(0, len(Triloops) // 5):
        pf.exptri[i] = RESCALE_BF(Triloop37[i], TriloopdH[i], TT, kTn)

    for i in range(0, len(Hexaloops) // 9):
        pf.exphex[i] = RESCALE_BF(Hexaloop37[i], HexaloopdH[i], TT, kTn)

    for i in range(NBPAIRS + 1):
        pf.expMLintern[i] = RESCALE_BF(ML_intern37, ML_interndH, TT, kTn)

    for i in range(NBPAIRS + 1):
        for j in range(5):
            if md.dangles:
                pf.expdangle5[i][j] = RESCALE_BF_SMOOTH(dangle5_37[i][j], dangle5_dH[i][j], TT, kTn)
                pf.expdangle3[i][j] = RESCALE_BF_SMOOTH(dangle3_37[i][j], dangle3_dH[i][j], TT, kTn)
            else:
                pf.expdangle3[i][j] = pf.expdangle5[i][j] = 1.0

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            pf.expstack[i][j] = RESCALE_BF(stack37[i][j], stackdH[i][j], TT, kTn)

    for i in range(NBPAIRS + 1):
        for j in range(5):
            for k in range(5):
                pf.expmismatchI[i][j][k] = RESCALE_BF(mismatchI37[i][j][k], mismatchIdH[i][j][k], TT, kTn)
                pf.expmismatch1nI[i][j][k] = RESCALE_BF(mismatch1nI37[i][j][k], mismatch1nIdH[i][j][k], TT, kTn)
                pf.expmismatchH[i][j][k] = RESCALE_BF(mismatchH37[i][j][k], mismatchHdH[i][j][k], TT, kTn)
                pf.expmismatch23I[i][j][k] = RESCALE_BF(mismatch23I37[i][j][k], mismatch23IdH[i][j][k], TT, kTn)
                if md.dangles:
                    pf.expmismatchM[i][j][k] = RESCALE_BF_SMOOTH(mismatchM37[i][j][k], mismatchMdH[i][j][k], TT, kTn)
                    pf.expmismatchExt[i][j][k] = RESCALE_BF_SMOOTH(mismatchExt37[i][j][k], mismatchExtdH[i][j][k], TT, kTn)
                else:
                    pf.expmismatchM[i][j][k] = pf.expmismatchExt[i][j][k] = 1.0

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            for k in range(5):
                for l in range(5):
                    pf.expint11[i][j][k][l] = RESCALE_BF(int11_37[i][j][k][l], int11_dH[i][j][k][l], TT, kTn)

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            for k in range(5):
                for l in range(5):
                    for m in range(5):
                        pf.expint21[i][j][k][l][m] = RESCALE_BF(int21_37[i][j][k][l][m], int21_dH[i][j][k][l][m], TT, kTn)

    for i in range(NBPAIRS + 1):
        for j in range(NBPAIRS + 1):
            for k in range(5):
                for l in range(5):
                    for m in range(5):
                        for n in range(5):
                            pf.expint22[i][j][k][l][m][n] = RESCALE_BF(int22_37[i][j][k][l][m][n], int22_dH[i][j][k][l][m][n], TT, kTn)

    pf.Tetraloops = Tetraloops
    pf.Triloops = Triloops
    pf.Hexaloops = Hexaloops

    pf.SaltMLbase = pf.SaltMLclosing = pf.SaltDPXInit = 0.0

    if salt == saltStandard:
        pf.expSaltStack = 1.0
    else:
        pf.expSaltStack = math.exp(-Duplex.vrna_salt_stack(salt, saltT, md.helical_rise) * 10.0 / kTn)
        Duplex.vrna_salt_ml(pf.SaltLoopDbl, md.saltMLLower, md.saltMLUpper, pf.SaltMLbase, pf.SaltMLclosing)

        if md.saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT:
            pf.SaltDPXInit = md.saltDPXInit
        elif md.saltDPXInit:
            pf.SaltDPXInit = Duplex.vrna_salt_duplex_init(md)

        pf.expMLclosing *= math.exp(-pf.SaltMLbase * 10.0 / kTn)
        pf.expMLclosing *= math.exp(-pf.SaltMLclosing * 10.0 / kTn)
        pf.expMLbase *= math.exp(-pf.SaltMLbase * 10.0 / kTn)
        for i in range(NBPAIRS + 1):
            pf.expMLintern[i] *= math.exp(-pf.SaltMLbase * 10.0 / kTn)

        pf.expDuplexInit *= math.exp(-pf.SaltDPXInit * 10.0 / kTn)

    return pf





def vrna_exp_params_comparative(n_seq: int, md: vrna_md_t) -> vrna_exp_param_t:
    if md:
        return get_exp_params_ali(md, n_seq, -1.0)
    else:
        md = vrna_md_t()
        vrna_md_set_default(md)
        return get_exp_params_ali(md, n_seq, -1.0)



VRNA_OPTION_PF = 1 << 1
    
def vrna_params_prepare(fc:vrna_fold_compound_t, options:int):
    if fc:
        md_p = fc.params.model_details

        if options & VRNA_OPTION_PF:
            if fc.exp_params:
                if fc.exp_params.model_details != md_p:
                    fc.exp_params = None

            if not fc.exp_params:
                if fc.type == VRNA_FC_TYPE_SINGLE:
                    fc.exp_params = vrna_exp_params(md_p)
                else:
                    fc.exp_params = vrna_exp_params_comparative(fc.n_seq, md_p)
# too long to read , chaos ##########################  
                
def add_params(fc: vrna_fold_compound_t, md_p: vrna_md_t, options: int):
    if fc.params:
        if fc.params.model_details != md_p:
            fc.params = None
    
    if not fc.params:
        fc.params = vrna_params(md_p)
    
    vrna_params_prepare(fc, options)
    
    
def sanitize_bp_span(fc:vrna_fold_compound_t, options:int):
    md = fc.params.model_details

    # 确保 min_loop_size, max_bp_span 和 window_size 合理
    if options & VRNA_OPTION_WINDOW:
        if md.window_size <= 0:
            md.window_size = fc.length
        elif md.window_size > fc.length:
            md.window_size = fc.length

        fc.window_size = md.window_size
    else:
        # 非局部折叠模式
        md.window_size = fc.length

    if md.max_bp_span <= 0 or md.max_bp_span > md.window_size:
        md.max_bp_span = md.window_size

def vrna_strsplit(string: Optional[str], delimiter: Optional[str] = '&') -> Optional[List[str]]:
    if not string:
        return None
    
    if delimiter:
        delim = delimiter[0]
    else:
        delim = '&'
    
    # Split the string using the given delimiter
    split = string.split(delim)
    
    # Remove empty strings if any
    split = [s for s in split if s]

    # Return the split list
    return split
import ctypes
def vrna_realloc(ptr, size):
    # Placeholder for actual realloc function
    return ctypes.cast(ptr, ctypes.POINTER(ctypes.c_ubyte * size)).contents

def vrna_sequence_add(vc:vrna_fold_compound_t, string:str, options:int) -> int:
    if vc and vc.type == VRNA_FC_TYPE_SINGLE and string:
        add_length = len(string)

        # Reallocate nucleotides
        vc.nucleotides = vrna_realloc(vc.nucleotides,
                                      ctypes.sizeof(ctypes.POINTER(ctypes.c_void_p)) * (vc.strands + 1))

        # Dummy placeholder for actual sequence setting
        # set_sequence(vc.nucleotides[vc.strands], string, None, vc.params.contents.model_details, options)

        # Increase strands counter
        vc.strands += 1

        # Reallocate sequence
        new_sequence_length = vc.length + add_length + 1
        vc.sequence = vrna_realloc(vc.sequence,
                                   ctypes.sizeof(ctypes.c_char) * new_sequence_length)
        vc.sequence[vc.length:vc.length + add_length] = string.encode('utf-8')

        # Reallocate encoding
        new_encoding_length = vc.length + add_length + 2
        vc.sequence_encoding = vrna_realloc(vc.sequence_encoding,
                                            ctypes.sizeof(ctypes.c_short) * new_encoding_length)

        # Dummy placeholder for actual encoding adjustment
        # vc.sequence_encoding[vc.length + 1:vc.length + add_length + 1] = vc.nucleotides[vc.strands - 1].encoding[1:]

        # Restore circular encoding
        vc.sequence_encoding[vc.length + add_length + 1] = vc.sequence_encoding[1]
        vc.sequence_encoding[0] = vc.sequence_encoding[vc.length + add_length]

        # Reallocate and adjust encoding2
        vc.sequence_encoding2 = vrna_realloc(vc.sequence_encoding2,
                                             ctypes.sizeof(ctypes.c_short) * new_encoding_length)
        
        # Dummy placeholder for actual encoding2 adjustment
        # enc = vrna_seq_encode_simple(vc.nucleotides[vc.strands - 1].string, vc.params.contents.model_details)
        # vc.sequence_encoding2[vc.length + 1:vc.length + add_length + 1] = enc[1:]
        
        vc.sequence_encoding2[vc.length + add_length + 1] = vc.sequence_encoding2[1]
        vc.sequence_encoding2[0] = vc.length + add_length

        # Increase length property of the fold compound
        vc.length += add_length

        return 1

    return 0


def vrna_idx_col_wise(length: int) -> List[int]:
    idx = [0] * (length + 1)

    for i in range(1, length + 1):
        idx[i] = (i * (i - 1)) // 2

    return idx



def vrna_ptypes(S: List[int], md: vrna_md_t) -> str:
    import numpy as np

    n = S[0]
    min_loop_size = md.min_loop_size

    if n > 1000000:
        print("vrna_ptypes@alphabet.c: sequence length of %d exceeds addressable range" % n)
        return None

    ptype = np.zeros(((n * (n + 1)) // 2 + 2,), dtype=np.uint8)
    idx = vrna_idx_col_wise(n)

    for k in range(1, n - min_loop_size):
        for l in range(1, 3):
            i = k
            j = i + min_loop_size + l
            if j > n:
                continue

            type_ = md.pair[S[i]][S[j]]
            ntype = 0
            otype = 0
            while i >= 1 and j <= n:
                if i > 1 and j < n:
                    ntype = md.pair[S[i - 1]][S[j + 1]]

                if md.noLP and not otype and not ntype:
                    type_ = 0  # i.j can only form isolated pairs

                ptype[idx[j] + i] = type_
                otype = type_
                type_ = ntype
                i -= 1
                j += 1

    return ptype.tobytes()

def wrap_get_ptypes(S: List[int], md: vrna_md_t) -> Optional[str]:
    if S:
        n = S[0]
        ptype = bytearray((n * (n + 1)) // 2 + 2)
        idx = vrna_idx_row_wise(n)
        min_loop_size = md.min_loop_size

        for k in range(1, n - min_loop_size):
            for l in range(1, 3):  # l ranges from 1 to 2 inclusive
                type_ = 0
                ntype = 0
                otype = 0
                i = k
                j = i + min_loop_size + l
                if j > n:
                    continue

                type_ = md.pair[S[i]][S[j]]
                while (i >= 1) and (j <= n):
                    if (i > 1) and (j < n):
                        ntype = md.pair[S[i - 1]][S[j + 1]]

                    if md.noLP and (otype == 0) and (ntype == 0):
                        type_ = 0  # i.j can only form isolated pairs

                    ptype[idx[i] - j] = chr(type_)
                    otype = type_
                    type_ = ntype
                    i -= 1
                    j += 1

        return ptype.decode('ascii')
    else:
        return None



def get_ptypes(S: List[int], md: vrna_md_t, idx_type: int) -> Optional[str]:
    if S:
        if len(S) > 1000000:
            print(
                f"get_ptypes@alphabet.c: sequence length of {S[0]} exceeds addressable range"
            )
            return None

        if idx_type:
            return wrap_get_ptypes(S, md)
        else:
            return vrna_ptypes(S, md)
    else:
        return None

Law_and_Order = "_ACGUTXKI"

def vrna_nucleotide_decode(enc: int, md: Optional[vrna_md_t]) -> str:
    if md:
        if md.energy_set > 0:
            return chr(enc + ord('A') - 1)
        else:
            return Law_and_Order[enc]
    else:
        return chr(0)



def vrna_aln_consensus_sequence(alignment: List[str], md_p: Optional[vrna_md_t]) -> Optional[str]:
    if alignment:
        n = len(alignment[0])
        if n > 0:
            # Check alignment for consistency
            for s in range(1, len(alignment)):
                if len(alignment[s]) != n:
                    print(
                        f"vrna_aln_consensus_sequence: Length of aligned sequence #{s + 1} does not match length of first sequence\n"
                        f"{alignment[s]}\n\n"
                    )
                    return None

            n_seq = len(alignment)

            md = md_p if md_p else vrna_md_set_default()

            consensus = [''] * n

            for i in range(n):
                freq = [0] * 8

                for s in range(n_seq):
                    nucleotide = alignment[s][i]
                    encoded = Duplex.vrna_nucleotide_encode(nucleotide, md)
                    freq[encoded] += 1

                max_freq = max(freq)
                c = freq.index(max_freq)

                if c > 4:
                    c += 1  # Skip T

                consensus[i] = vrna_nucleotide_decode(c, md)

            return ''.join(consensus)
    return None


VRNA_SEQ_RNA = 1





def vrna_sequence_prepare(fc: vrna_fold_compound_t) -> None:
    if fc:
        if fc.strand_number is not None:
            del fc.strand_number
        if fc.strand_order is not None:
            del fc.strand_order
        if fc.strand_order_uniq is not None:
            del fc.strand_order_uniq
        if fc.strand_start is not None:
            del fc.strand_start
        if fc.strand_end is not None:
            del fc.strand_end

        fc.strand_order = None
        fc.strand_order_uniq = None
        fc.strand_start = None
        fc.strand_end = None

        fc.strand_number = [0] * (fc.length + 2)

        if fc.type == VRNA_FC_TYPE_SINGLE:
            fc.strand_order_uniq = list(range(fc.strands + 1))
            fc.strand_order = list(range(fc.strands + 1))

            fc.strand_start = [0] * (fc.strands + 1)
            fc.strand_end = [0] * (fc.strands + 1)

            fc.strand_start[0] = 1
            fc.strand_end[0] = fc.strand_start[0] + fc.nucleotides[0].length - 1

            for cnt in range(1, fc.strands):
                fc.strand_start[cnt] = fc.strand_end[cnt - 1] + 1
                fc.strand_end[cnt] = fc.strand_start[cnt] + fc.nucleotides[cnt].length - 1
                for i in range(fc.strand_start[cnt], fc.strand_end[cnt] + 1):
                    fc.strand_number[i] = cnt

            fc.strand_number[0] = fc.strand_number[1]
            fc.strand_number[fc.length + 1] = fc.strand_number[fc.length]

        elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
            fc.nucleotides = [None] * (fc.strands + 1)
            fc.nucleotides[0] = vrna_seq_t(string=None, type=VRNA_SEQ_RNA, length=fc.length)

            fc.strand_order_uniq = [0, 1]
            fc.strand_order = [0, 1]

            fc.strand_start = [0, 1]
            fc.strand_end = [0, 1]
            fc.strand_start[0] = 1
            fc.strand_end[0] = fc.strand_start[0] + fc.length - 1

VRNA_SEQUENCE_RNA = 1
VRNA_OPTION_EVAL_ONLY = 1 << 3
WITH_PTYPE = 1
WITH_PTYPE_COMPAT = 2
VRNA_OPTION_WINDOW = 1 << 4


def set_fold_compound(fc:vrna_fold_compound_t, options:int, aux:int):
    sequence = None
    sequences = None
    md_p = fc.params.model_details

    if fc.type == VRNA_FC_TYPE_SINGLE:
        sequence = fc.sequence
        fc.sequence = None
        fc.length = 0

        # 将输入序列分割为默认分隔符 '&'
        sequences = vrna_strsplit(sequence, None)

        # 将单个序列添加到折叠复合物
        for seq in sequences:
            vrna_sequence_add(fc, seq, VRNA_SEQUENCE_RNA)

        if fc.strands > 1:
            fc.cutpoint = fc.nucleotides[0].length + 1

            # if md_p.min_loop_size == TURN:
            #     md_p.min_loop_size = 0

        if not options & VRNA_OPTION_EVAL_ONLY:
            if fc.strands > 1:
                min_loop_size = md_p.min_loop_size
                md_p.min_loop_size = 0
                fc.ptype = vrna_ptypes(fc.sequence_encoding2, md_p) if aux & WITH_PTYPE else None
                md_p.min_loop_size = min_loop_size
            else:
                fc.ptype = vrna_ptypes(fc.sequence_encoding2, md_p) if aux & WITH_PTYPE else None

            fc.ptype_pf_compat = get_ptypes(fc.sequence_encoding2, md_p, 1) if aux & WITH_PTYPE_COMPAT else None

    elif fc.type == VRNA_FC_TYPE_COMPARATIVE:
        sequences = fc.sequences
        fc.length = length = fc.length

        fc.cons_seq = vrna_aln_consensus_sequence(sequences, md_p)
        fc.S_cons = vrna_seq_encode_simple(fc.cons_seq, md_p)

        fc.pscore = vrna_alloc(int, (length * (length + 1)) // 2 + 2)
        fc.pscore_pf_compat = vrna_alloc(int, (length * (length + 1)) // 2 + 2) if aux & WITH_PTYPE_COMPAT else None

        oldAliEn = fc.oldAliEn = md_p.oldAliEn

        fc.S = vrna_alloc(short, fc.n_seq + 1)
        fc.S5 = vrna_alloc(short, fc.n_seq + 1)
        fc.S3 = vrna_alloc(short, fc.n_seq + 1)
        fc.a2s = vrna_alloc(int, fc.n_seq + 1)
        fc.Ss = vrna_alloc(char, fc.n_seq + 1)

        for s in range(fc.n_seq):
            vrna_aln_encode(fc.sequences[s], fc.S[s], fc.S5[s], fc.S3[s], fc.Ss[s], fc.a2s[s], md_p)

        fc.S5[fc.n_seq] = None
        fc.S3[fc.n_seq] = None
        fc.a2s[fc.n_seq] = None
        fc.Ss[fc.n_seq] = None
        fc.S[fc.n_seq] = None

    vrna_sequence_prepare(fc)

    if not options & VRNA_OPTION_WINDOW and fc.length <= 1000000:
        fc.iindx = vrna_idx_row_wise(fc.length)
        fc.jindx = vrna_idx_col_wise(fc.length)


    




def vrna_fold_compound(sequence: str, md_p: Optional[vrna_md_t], options: int) -> Optional[vrna_fold_compound_t]:
    # Check if sequence is None
    if sequence is None:
        return None

    # Sanity check for sequence length
    length = len(sequence)
    if length == 0:
        print("vrna_fold_compound@data_structures.c: sequence length must be greater 0")
        return None

    # Check sequence length against addressable range
    # if length > vrna_sequence_length_max(options):
    if length > 1000000:
        print(f"vrna_fold_compound@data_structures.c: sequence length of {length} exceeds addressable range")
        return None

    # Initialize fold compound
    fc = init_fc_single()

    # Set length and sequence in fold compound
    fc.length = length
    fc.sequence = sequence

    aux_options = 0

    # Get a copy of the model details
    md = md_p if md_p else vrna_md_set_default(md)

    # Add energy parameters
    add_params(fc, md, options)

    sanitize_bp_span(fc, options)

    if options & VRNA_OPTION_WINDOW:
        set_fold_compound(fc, options, aux_options)

        if not (options & VRNA_OPTION_EVAL_ONLY):
            # Add minimal hard constraint data structure
            vrna_hc_init_window(fc)

            # Add DP matrices
            vrna_mx_add(fc, VRNA_MX_WINDOW, options)
    else:
        # Regular global structure prediction
        aux_options |= WITH_PTYPE

        if options & VRNA_OPTION_PF:
            aux_options |= WITH_PTYPE_COMPAT

        set_fold_compound(fc, options, aux_options)

        if not (options & VRNA_OPTION_EVAL_ONLY):
            # Add default hard constraints
            vrna_hc_init(fc)

            # Add DP matrices (if required)
            vrna_mx_add(fc, VRNA_MX_DEFAULT, options)

    return fc
VRNA_OPTION_DEFAULT = 0
VRNA_OPTION_HYBRID = 1 << 2
from .remove import vrna_md_t
def generate_struc():
    struc = "(((((((((((((((((((((()))))))))))...)))))))))))"
    seq = "aactgcccactcagtacatcaa&TTGATGTACTGCCAAGTGGGCAGTT"
    # The length of the sequence (or sequence alignment) vrna_fold_compound_t
    # n = vc.length
    n = len(struc)
    # vrna_mx_pf_t *matrices
    
    # probs = matrices.probs
    # probs = [0] * n**2
    vc = vrna_fold_compound(seq,vrna_md_t(),VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID)

    return pf_create_bppm(vc, None)

if __name__ == "__main__":
    print(generate_struc())

