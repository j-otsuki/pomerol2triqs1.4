from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
from pytriqs.utility.comparison_tests import *
import numpy as np
from itertools import product

# Non-interacting impurity in a bath

####################
# Input parameters #
####################

beta = 10.0             # Inverse temperature
e_d = 0.1               # Local level
h = 0.15                # Magnetic field

# Levels of the bath sites
# epsilon = [-1.0, 1.0]
epsilon = [-1.0, 1.2]

# Hopping amplitudes
# V = [0.5, 0.5]
V = [0.5, 0.6]

spin_names = ("up", "dn")

# Number of Matsubara frequencies for GF calculation
n_iw = 200

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

# Local Hamiltonian
H_loc = e_d*(n('up', 0) + n('dn', 0)) + h*(n('up', 0) - n('dn', 0))

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

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

###########
# G^{(2)} #
###########

common_g2_params = {'channel' : "PH",
                    'gf_struct' : gf_struct,
                    'beta' : beta,
                    'n_f' : g2_n_wf,
                    'n_b' : g2_n_wb, }

G2_ph_uuuu = ed.G2_iw_freq_box( index1=('up',0), index2=('up',0), index3=('up',0), index4=('up',0), **common_g2_params )
G2_ph_dddd = ed.G2_iw_freq_box( index1=('dn',0), index2=('dn',0), index3=('dn',0), index4=('dn',0), **common_g2_params )
G2_ph_uudd = ed.G2_iw_freq_box( index1=('up',0), index2=('up',0), index3=('dn',0), index4=('dn',0), **common_g2_params )
G2_ph_dduu = ed.G2_iw_freq_box( index1=('dn',0), index2=('dn',0), index3=('up',0), index4=('up',0), **common_g2_params )
G2_ph_uddu = ed.G2_iw_freq_box( index1=('up',0), index2=('dn',0), index3=('dn',0), index4=('up',0), **common_g2_params )
G2_ph_duud = ed.G2_iw_freq_box( index1=('dn',0), index2=('up',0), index3=('up',0), index4=('dn',0), **common_g2_params )

G2_ph_uuuu_wick = G2_ph_uuuu.copy()
G2_ph_dddd_wick = G2_ph_dddd.copy()
G2_ph_uudd_wick = G2_ph_uudd.copy()
G2_ph_dduu_wick = G2_ph_dduu.copy()
G2_ph_uddu_wick = G2_ph_uddu.copy()
G2_ph_duud_wick = G2_ph_duud.copy()

G = lambda s, i: G_iw[s].data[i + n_iw, 0, 0]

for (i, j1, j2) in product(range(g2_n_wb),range(2*g2_n_wf),range(2*g2_n_wf)):
    wb = i
    wf1 = j1 - g2_n_wf
    wf2 = j2 - g2_n_wf
    d_w_0 = int(wb == 0)
    # d_w_nu_nup = int(m == i + ip + 1)
    d_nu_nup = int(j1 == j2)

    G2_ph_uudd_wick[i,j1,j2] = beta * d_w_0 * G("up", wf1) * G("dn", wf2)
    G2_ph_dduu_wick[i,j1,j2] = beta * d_w_0 * G("dn", wf1) * G("up", wf2)

    G2_ph_uuuu_wick[i,j1,j2] = beta * d_w_0 * G("up", wf1) * G("up", wf2) - beta * d_nu_nup * G("up", wf1+wb) * G("up", wf1)
    G2_ph_dddd_wick[i,j1,j2] = beta * d_w_0 * G("dn", wf1) * G("dn", wf2) - beta * d_nu_nup * G("dn", wf1+wb) * G("dn", wf1)

    G2_ph_uddu_wick[i,j1,j2] =                                            - beta * d_nu_nup * G("dn", wf1+wb) * G("up", wf1)
    G2_ph_duud_wick[i,j1,j2] =                                            - beta * d_nu_nup * G("up", wf1+wb) * G("dn", wf1)


diff = lambda g1, g2: np.linalg.norm(g1-g2)/g1.size  # mean squared error
# diff = lambda g1, g2: np.linalg.norm(g1-g2, ord=np.inf)  # max(|x|)
gfs_are_close = lambda g1, g2: diff(g1,g2) < 1e-8

print diff(G2_ph_uuuu_wick,G2_ph_uuuu)
print diff(G2_ph_dddd_wick,G2_ph_dddd)
print diff(G2_ph_uudd_wick,G2_ph_uudd)
print diff(G2_ph_dduu_wick,G2_ph_dduu)
print diff(G2_ph_uddu_wick,G2_ph_uddu)
print diff(G2_ph_duud_wick,G2_ph_duud)

assert( gfs_are_close(G2_ph_uuuu_wick, G2_ph_uuuu) )
assert( gfs_are_close(G2_ph_dddd_wick, G2_ph_dddd) )
assert( gfs_are_close(G2_ph_uudd_wick, G2_ph_uudd) )
assert( gfs_are_close(G2_ph_dduu_wick, G2_ph_dduu) )
assert( gfs_are_close(G2_ph_uddu_wick, G2_ph_uddu) )
assert( gfs_are_close(G2_ph_duud_wick, G2_ph_duud) )
