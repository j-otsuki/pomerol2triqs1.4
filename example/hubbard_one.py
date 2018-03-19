from pytriqs.applications.impurity_solvers.pomerol2triqs.hubbard_one_solver import Solver

from pytriqs.archive import HDFArchive
from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
from pytriqs.operators.util.U_matrix import U_matrix
from pytriqs.operators.util.hamiltonians import h_int_slater
from pytriqs.operators.util.observables import N_op, S_op, L_op, LS_op
from pytriqs.utility import mpi
from itertools import product
import numpy as np

# Example of using hubbard_one_solver
# 5-orbital atom with Slater interaction term

####################
# Input parameters #
####################

# Angular momentumd
L = 2

# Inverse temperature
beta = 10.0

# Slater parameters
U = 5.0
J = 0.1
# F0 = U
# F2 = J*(14.0/(1.0 + 0.625))
# F4 = F2*0.625

mu = U*0.5  # n=1
# mu = U*1.5  # n=2

# spin-orbit coupling
ls_coupling = 0.0

# Number of Matsubara frequencies for GF calculation
n_iw = 200

# Number of imaginary time slices for GF calculation
n_tau = 1001

# Number of bosonic Matsubara frequencies for G^2 calculations
g2_n_wb = 5
# Number of fermionic Matsubara frequencies for G^2 calculations
g2_n_wf = 10

#####################
# Input for Pomerol #
#####################

spin_names = ("up", "dn")
orb_names = map(str, range(-L, L+1))
# orb_names = range(-L, L+1)  # This also works

flag_so = True if ls_coupling != 0 else False

# GF and operator structure
if not flag_so:
    # ('up', '-2'), ('dn', '-2'), etc
    gf_struct = set_operator_structure(spin_names, orb_names, off_diag=True)
    map_operator_structure = { (s,o) : (s,o) for s, o in product(spin_names, orb_names) }
else:
    # ('ud', 'up_-2'), ('ud', 'dn_-2'), etc
    gf_struct = {'ud': [str(s)+'_'+str(o) for s, o in product(spin_names, orb_names)] }
    map_operator_structure = { (s,o) : ('ud', str(s)+'_'+str(o)) for s, o in product(spin_names, orb_names) }
# print gf_struct
# print map_operator_structure

# Operators
N = N_op(spin_names, orb_names, map_operator_structure=map_operator_structure)
Sz = S_op('z', spin_names, orb_names, map_operator_structure=map_operator_structure)
Lz = L_op('z', spin_names, orb_names, map_operator_structure=map_operator_structure, basis = 'spherical')
Jz = Sz + Lz

# LS coupling term
# LS = ls_op(spin_names, orb_names, off_diag=off_diag, basis = 'spherical')

# Hamiltonian
# U_mat = U_matrix(L, radial_integrals = [F0,F2,F4], basis='spherical')
U_mat = U_matrix(L, U_int=U, J_hund=J, basis='spherical')
H_int = h_int_slater(spin_names, orb_names, U_mat, map_operator_structure=map_operator_structure)

# H -= mu*N
# H += ls_coupling * LS

# Double check that we are actually using integrals of motion
# h_comm = lambda op: H*op - op*H
# str_if_commute = lambda op: "= 0" if h_comm(op).is_zero() else "!= 0"

# print "Check integrals of motion:"
# print "[H, N]", str_if_commute(N)
# print "[H, Sz]", str_if_commute(Sz)
# print "[H, Lz]", str_if_commute(Lz)
# print "[H, Jz]", str_if_commute(Jz)

E_levels = {}
if not flag_so:
    for sp in spin_names:
        E_levels[sp] = np.diag( [-mu]*len(orb_names) )
else:
    # TODO: add matrix elements of spin-orbit coupling
    E_levels['ud'] = np.diag( [-mu]*len(orb_names)*2 )

# print type(E_levels)
# print type(E_levels['up'])

# const_of_motion = [N, Sz, Lz]
const_of_motion = [N, Sz, Lz] if not flag_so else [N, Jz]

#####################
# Pomerol ED solver #
#####################

S = Solver(beta, gf_struct, n_iw, spin_orbit=flag_so, verbose=True)

S.solve(H_int, E_levels, const_of_motion=const_of_motion, file_quantum_numbers="quantum_numbers.dat", file_eigenvalues="eigenvalues.dat")

if mpi.is_master_node():
    with HDFArchive('hubbard_one.out.h5', 'w') as ar:
        ar['G_iw'] = S.G_iw
        ar['G0_iw'] = S.G0_iw
        ar['Sigma_iw'] = S.Sigma_iw
