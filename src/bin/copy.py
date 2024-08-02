from sys import path
path.append("/home/wayne/Repository/ViennaRNA-2.6.4/src/bin")
from dimer_cofold import vrna_fold_compound
from opt import options, init_default_options

# global constant define start ##############
VRNA_OPTION_DEFAULT = 0
VRNA_OPTION_HYBRID = 1 << 2





def main():
    sequences = "CUUCCUCGGGUUCAAAGCUGGAUU&GUCCAGUUUUCCCAGGAAU"
    opts = options()
    num_input = 0
    init_default_options(opts)
    vc = vrna_fold_compound(sequences, opts.md, VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID)
    print(vc.iindx)
    
    


if __name__ == "__main__":
    main()