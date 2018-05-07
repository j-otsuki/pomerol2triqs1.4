from pytriqs.archive import HDFArchive
from pytriqs.gf import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.operators.util.hamiltonians import h_int_kanamori
from pytriqs.utility import mpi
from pytriqs.applications.impurity_solvers.pomerol2triqs import PomerolED
import numpy as np
from itertools import product

# 2-orbital Hubbard-Kanamori atom

####################
# Input parameters #
####################

beta = 10.0         # Inverse temperature
num_orb = 2         # Number of orbitals
U = 2.0             # Coulomb repulsion
mu = 1.5            # Chemical potential
J = 0.2             # Hund coupling

spin_names = ("up", "dn")
orb_names = range(num_orb)

# Number of Matsubara frequencies for GF calculation
n_iw = 1024
# Number of imaginary time slices for GF calculation
n_tau = 2001

# Energy window for real frequency GF calculation
energy_window = (-5, 5)
# Number of frequency points for real frequency GF calculation
n_w = 1000

# Number of bosonic Matsubara frequencies for G^2 calculations
g2_n_wb = 5
# Number of fermionic Matsubara frequencies for G^2 calculations
g2_n_wf = 10
# Number of Legendre coefficients for G^2 calculations
# g2_n_l = 10

gf_struct = {"up" : orb_names, "dn" : orb_names}
print "Block structure of single-particle Green's functions:", gf_struct

# Conversion from TRIQS to Pomerol notation for operator indices
# TRIQS: block_name, inner_index
# Pomerol: site_label, orbital_index, spin_name
index_converter = {(sn, o) : ("loc", o, "down" if sn == "dn" else "up")
                   for sn, o in product(spin_names, orb_names)}
print "index_converter:", index_converter

# Make PomerolED solver object
ed = PomerolED(index_converter, verbose = True)

# Number of particles on the impurity
N = sum(n(sn, o) for sn, o in product(spin_names, orb_names))

# Hamiltonian
H = h_int_kanamori(spin_names, orb_names,
                   np.array([[0, U-3*J], [U-3*J, 0]]),
                   np.array([[U, U-2*J], [U-2*J, U]]),
                   J, True)
H -= mu*N

# Diagonalize H
ed.diagonalize(H)

# save data
ed.save_quantum_numbers("quantum_numbers.dat")
ed.save_eigenvalues("eigenvalues.dat")

# set density-matrix cutoff
ed.set_density_matrix_cutoff(1e-10)

# Compute G(i\omega)
G_iw = ed.G_iw(gf_struct, beta, n_iw)

# Compute G(\tau)
G_tau = ed.G_tau(gf_struct, beta, n_tau)

# Compute G(\omega)
G_w = ed.G_w(gf_struct, beta, energy_window, n_w, 0.01)

###########
# G^{(2)} #
###########

common_g2_params = {'channel' : "PH",
                    'gf_struct' : gf_struct,
                    'beta' : beta,
                    'n_f' : g2_n_wf,
                    'n_b' : g2_n_wb, }

###############################
# G^{(2)}(i\omega;i\nu,i\nu') #
###############################

# G2_iw = ed.G2_iw_legacy( index1=('up',0), index2=('dn',0), index3=('dn',1), index4=('up',1), **common_g2_params )
G2_iw = ed.G2_iw_freq_box( vec_four_indices=[(('up',0), ('dn',0), ('dn',1), ('up',1)),], **common_g2_params )[0]
print type(G2_iw)
print G2_iw.shape


# # Compute G^{(2),ph}(i\omega;i\nu,i\nu'), AABB block order
# G2_iw_inu_inup_ph_AABB = ed.G2_iw_inu_inup(channel = "PH",
#                                            block_order = "AABB",
#                                            n_inu = g2_n_inu,
#                                            **common_g2_params)

#########################
# G^{(2)}(i\omega;l,l') #
#########################

# # Compute G^{(2),ph}(i\omega;l,l'), AABB block order
# G2_iw_l_lp_ph_AABB = ed.G2_iw_l_lp(channel = "PH",
#                                    block_order = "AABB",
#                                    n_l = g2_n_l,
#                                    **common_g2_params)

################
# Save results #
################

if mpi.is_master_node():
    with HDFArchive('2band.atom.h5', 'w') as ar:
        ar['G_iw'] = G_iw
        ar['G_tau'] = G_tau
        ar['G_w'] = G_w
        ar['G2_ph'] = G2_iw
