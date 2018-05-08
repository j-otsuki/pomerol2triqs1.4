from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
from pytriqs.utility.comparison_tests import *
import numpy as np
from itertools import product

# Single orbital Anderson model

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
U = 4.3                 # Coulomb repulsion
mu = 1.1                # Chemical potential

# Levels of the bath sites
epsilon = [-1.5, 1.3]
# Hopping amplitudes
V = [0.52, 0.55]

spin_names = ("up", "dn")

# Number of bosonic Matsubara frequencies for G^2 calculations
g2_n_wb = 3
# Number of fermionic Matsubara frequencies for G^2 calculations
g2_n_wf = 3

# GF structure
gf_struct = {'up' : [0], 'dn' : [0]}

# Conversion from TRIQS to Pomerol notation for operator indices
index_converter = {}
index_converter.update({(sn, 0) : ("loc", 0, "down" if sn == "dn" else "up") for sn in spin_names})
index_converter.update({("B%i_%s" % (k, sn), 0) : ("bath" + str(k), 0, "down" if sn == "dn" else "up")
                        for k, sn in product(range(len(epsilon)), spin_names)})

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Number of particles on the impurity
H_loc = -mu*(n('up', 0) + n('dn', 0)) + U * n('up', 0) * n('dn', 0)

# Bath Hamiltonian
H_bath = sum(eps*n("B%i_%s" % (k, sn), 0)
             for sn, (k, eps) in product(spin_names, enumerate(epsilon)))

# Hybridization Hamiltonian
H_hyb = Operator()
for k, v in enumerate(V):
    H_hyb += sum(        v   * c_dag("B%i_%s" % (k, sn), 0) * c(sn, 0) +
                 np.conj(v)  * c_dag(sn, 0) * c("B%i_%s" % (k, sn), 0)
                 for sn in spin_names)

# Complete Hamiltonian
H = H_loc + H_hyb + H_bath

# Diagonalize H
ed.diagonalize(H)

###########
# G^{(2)} #
###########

common_g2_params = {'channel' : "PH",
                    'gf_struct' : gf_struct,
                    'beta' : beta,}

four_indices = []
four_indices.append( (('up',0), ('dn',0), ('dn',0), ('up',0)) )  # uddu

# compute in a low-freq box
G2_iw_freq_box = ed.G2_iw_freq_box( four_indices=four_indices, n_f=g2_n_wf, n_b=g2_n_wb, **common_g2_params )[0]
print G2_iw_freq_box.shape

# compute for fixed freqs
three_freqs = [ (iwb, iwf1-g2_n_wf, iwf2-g2_n_wf) for iwb, iwf1, iwf2 in product( range(g2_n_wb), range(2*g2_n_wf), range(2*g2_n_wf) ) ]
G2_iw_freq_fix = ed.G2_iw_freq_fix( four_indices=four_indices, three_freqs=three_freqs, **common_g2_params )[0]
print G2_iw_freq_fix.shape
# reshape from 1D array to 3D array
G2_iw_freq_fix = np.reshape( G2_iw_freq_fix, (g2_n_wb, 2*g2_n_wf, 2*g2_n_wf) )
print G2_iw_freq_fix.shape

assert( np.allclose(G2_iw_freq_box, G2_iw_freq_fix) )
