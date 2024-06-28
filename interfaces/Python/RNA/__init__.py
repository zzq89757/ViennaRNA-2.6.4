__license__ = '''
    Disclaimer and Copyright

The programs, library and source code of the Vienna RNA Package are free
software. They are distributed in the hope that they will be useful
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Permission is granted for research, educational, and commercial use
and modification so long as 1) the package and any derived works are not
redistributed for any fee, other than media costs, 2) proper credit is
given to the authors and the Institute for Theoretical Chemistry of the
University of Vienna.

If you want to include this software in a commercial product, please contact
the authors.

Note that the file ./src/ViennaRNA/plotting/naview.c has its own copyright attached.
'''

__doc__ = '''
A library for the prediction and comparison of RNA secondary structures.

Amongst other things, our implementations allow you to:

- predict minimum free energy secondary structures
- calculate the partition function for the ensemble of structures
- compute various equilibrium probabilities
- calculate suboptimal structures in a given energy range
- compute local structures in long sequences
- predict consensus secondary structures from a multiple sequence alignment
- predict melting curves
- search for sequences folding into a given structure
- compare two secondary structures 
- predict interactions between multiple RNA molecules
'''

__version__ = '2.6.4'

__author__ = 'Ronny Lorenz <rna@tbi.univie.ac.at>'

# Try loading RNA module by default and export all symbols to
# package namespace
try:
    from RNA.RNA import *
except Exception as E:
    raise ImportError('Failed to load ViennaRNA RNAlib Python wrapper')

