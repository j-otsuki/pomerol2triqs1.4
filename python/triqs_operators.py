"""Additional operators to pytriqs.operators.util.observables"""

import numpy as np
from pytriqs.operators.operators import Operator, n, c_dag, c
from pytriqs.operators.util.op_struct import get_mkind
from pytriqs.operators.util.U_matrix import spherical_to_cubic
from itertools import product

def ls_op(spin_names, orb_names, off_diag = None, map_operator_structure = None, basis='spherical', T=None):
    """one-body spin-orbit coupling operator"""

    spin_range = range(len(spin_names))

    l = (len(orb_names)-1)/2.0
    orb_range = range(int(2*l+1))

    pauli_matrix = {'x' : np.array([[0,1],[1,0]]),
                    'y' : np.array([[0,-1j],[1j,0]]),
                    'z' : np.array([[1,0],[0,-1]]),
                    '+' : np.array([[0,2],[0,0]]),
                    '-' : np.array([[0,0],[2,0]])}

    L_melem_dict = {'z' : lambda m,mp: m if np.isclose(m,mp) else 0,
                    '+' : lambda m,mp: np.sqrt(l*(l+1)-mp*(mp+1)) if np.isclose(m,mp+1) else 0,
                    '-' : lambda m,mp: np.sqrt(l*(l+1)-mp*(mp-1)) if np.isclose(m,mp-1) else 0,
                    'x' : lambda m,mp: 0.5*(L_melem_dict['+'](m,mp) + L_melem_dict['-'](m,mp)),
                    'y' : lambda m,mp: -0.5j*(L_melem_dict['+'](m,mp) - L_melem_dict['-'](m,mp))}

    # define S matrix
    S_matrix = {}
    for component in ['z', '+', '-']:
        pm = pauli_matrix[component]
        S_matrix[component] = np.array([[ 0.5*pm[n1,n2] for n2 in spin_range ] for n1 in spin_range ])

    # define L matrix
    L_matrix = {}
    for component in ['z', '+', '-']:
        L_melem = L_melem_dict[component]
        L_matrix[component] = np.array([[L_melem(o1-l,o2-l) for o2 in orb_range] for o1 in orb_range])

        # Transform from spherical basis if needed
        if basis == "cubic":
            if not np.isclose(np.mod(l,1),0):
                raise ValueError("L_op: cubic basis is only defined for the integer orbital momenta.")
            T = spherical_to_cubic(l)
        if basis == "other" and T is None: raise ValueError("L_op: provide T for other bases.")
        if T is not None: L_matrix = np.einsum("ij,jk,kl",np.conj(T),L_matrix,np.transpose(T))

    # LS_matrix[n1,n2,o1,o2] = sum_x S_matrix[x][n1,n2] * L_matrix[x][o1,o2]
    LS_matrix = np.einsum("ij,kl", S_matrix['z'], L_matrix['z'])\
              + np.einsum("ij,kl", S_matrix['+'], L_matrix['-']) * 0.5\
              + np.einsum("ij,kl", S_matrix['-'], L_matrix['+']) * 0.5

    mkind  = get_mkind(off_diag,map_operator_structure)
    ls = Operator()
    for n1, n2 in product(spin_range,spin_range):
        for o1, o2 in product(orb_range,orb_range):
            ls += c_dag(*mkind(spin_names[n1],orb_names[o1])) * LS_matrix[n1,n2,o1,o2] * c(*mkind(spin_names[n2],orb_names[o2]))
    return ls
