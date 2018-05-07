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
U = 2.0                 # Coulomb repulsion
mu = 1.0                # Chemical potential

# Levels of the bath sites
epsilon = [-1.0, 1.0]
# Hopping amplitudes
V = [0.5, 0.5]

spin_names = ("up", "dn")

# Number of bosonic Matsubara frequencies for G^2 calculations
g2_n_wb = 5
# Number of fermionic Matsubara frequencies for G^2 calculations
g2_n_wf = 5

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
                    'beta' : beta,
                    'n_f' : g2_n_wf,
                    'n_b' : g2_n_wb, }


# # Compute G^{(2),ph}(i\omega;i\nu,i\nu'), AABB block order
# G2_iw_ph_uuuu = ed.G2_iw_legacy( index1=('up',0), index2=('up',0), index3=('up',0), index4=('up',0), **common_g2_params )
# G2_iw_ph_dddd = ed.G2_iw_legacy( index1=('dn',0), index2=('dn',0), index3=('dn',0), index4=('dn',0), **common_g2_params )
# G2_iw_ph_uudd = ed.G2_iw_legacy( index1=('up',0), index2=('up',0), index3=('dn',0), index4=('dn',0), **common_g2_params )
# G2_iw_ph_dduu = ed.G2_iw_legacy( index1=('dn',0), index2=('dn',0), index3=('up',0), index4=('up',0), **common_g2_params )
#
# # Compute G^{(2),ph}(i\omega;i\nu,i\nu'), ABBA block order
# G2_iw_ph_uddu = ed.G2_iw_legacy( index1=('up',0), index2=('dn',0), index3=('dn',0), index4=('up',0), **common_g2_params )
# G2_iw_ph_duud = ed.G2_iw_legacy( index1=('dn',0), index2=('up',0), index3=('up',0), index4=('dn',0), **common_g2_params )


four_indices = []
four_indices.append( (('up',0), ('up',0), ('up',0), ('up',0)) )  # uuuu
four_indices.append( (('dn',0), ('dn',0), ('dn',0), ('dn',0)) )  # dddd
four_indices.append( (('up',0), ('up',0), ('dn',0), ('dn',0)) )  # uudd
four_indices.append( (('dn',0), ('dn',0), ('up',0), ('up',0)) )  # dduu
four_indices.append( (('up',0), ('dn',0), ('dn',0), ('up',0)) )  # uddu
four_indices.append( (('dn',0), ('up',0), ('up',0), ('dn',0)) )  # duud

G2_iw_ph = ed.G2_iw_freq_box( vec_four_indices=four_indices, **common_g2_params )
G2_iw_ph_uuuu = G2_iw_ph[0]
G2_iw_ph_dddd = G2_iw_ph[1]
G2_iw_ph_uudd = G2_iw_ph[2]
G2_iw_ph_dduu = G2_iw_ph[3]
G2_iw_ph_uddu = G2_iw_ph[4]
G2_iw_ph_duud = G2_iw_ph[5]


# Compute G^{(2),pp}(i\omega;i\nu,i\nu'), AABB block order
# G2_iw_inu_inup_pp_AABB = ed.G2_iw_inu_inup(channel = "PP",
#                                            block_order = "AABB",
#                                            n_inu = g2_n_inu,
#                                            **common_g2_params)

# Compute G^{(2),pp}(i\omega;i\nu,i\nu'), ABBA block order
# G2_iw_inu_inup_pp_ABBA = ed.G2_iw_inu_inup(channel = "PP",
#                                            block_order = "ABBA",
#                                            n_inu = g2_n_inu,
#                                            **common_g2_params)

diff = lambda g1, g2: np.linalg.norm(g1-g2)/g1.size  # mean squared error
# diff = lambda g1, g2: np.linalg.norm(g1-g2, ord=np.inf)  # max(|x|)
gfs_are_close = lambda g1, g2: diff(g1,g2) < 1e-8

if mpi.is_master_node():
    with HDFArchive('anderson_g2_matsubara_np%i.out.h5' % mpi.size, 'w') as ar:
        ar['H'] = H
        ar['G2_iw_ph_uuuu'] = G2_iw_ph_uuuu
        ar['G2_iw_ph_dddd'] = G2_iw_ph_dddd
        ar['G2_iw_ph_uudd'] = G2_iw_ph_uudd
        ar['G2_iw_ph_dduu'] = G2_iw_ph_dduu
        ar['G2_iw_ph_uddu'] = G2_iw_ph_uddu
        ar['G2_iw_ph_duud'] = G2_iw_ph_duud

    with HDFArchive("anderson_g2_matsubara.ref.h5", 'r') as ar:
        assert (ar['H'] - H).is_zero()
        assert ( gfs_are_close(ar['G2_iw_ph_uuuu'], G2_iw_ph_uuuu) )
        assert ( gfs_are_close(ar['G2_iw_ph_dddd'], G2_iw_ph_dddd) )
        assert ( gfs_are_close(ar['G2_iw_ph_uudd'], G2_iw_ph_uudd) )
        assert ( gfs_are_close(ar['G2_iw_ph_dduu'], G2_iw_ph_dduu) )
        assert ( gfs_are_close(ar['G2_iw_ph_uddu'], G2_iw_ph_uddu) )
        assert ( gfs_are_close(ar['G2_iw_ph_duud'], G2_iw_ph_duud) )
