
from pomerol2triqs import PomerolED
from pytriqs.gf.local import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
import pytriqs.utility.mpi as mpi
import numpy as np
from itertools import product

class Solver:
    def __init__(self, beta, gf_struct, spin_orbit=False, verbose=True):

        # TODO: spin_orbit
        if spin_orbit:
            print "*** hubbard_one_solver.Solver: spin_orbit not implemented yet"

        self.beta = beta
        self.gf_struct = gf_struct
        self.verbose = verbose

        if self.verbose and mpi.is_master_node():
            print "\n*** Hubbard I solver using Pomerol library"
            print "*** gf_struct =", gf_struct

        # get spin_name and orb_names from gf_struct
        self.__analyze_gf_struct(gf_struct)

        if self.verbose and mpi.is_master_node():
            print "*** spin_names =", self.spin_names
            print "*** orb_names  =", self.orb_names

        off_diag = True
        mkind = get_mkind(off_diag, None)
        # Conversion from TRIQS to Pomerol notation for operator indices
        index_converter = {mkind(sn, bn) : ("atom", bi, "down" if sn == "dn" else sn)
                           for sn, (bi, bn) in product(self.spin_names, enumerate(self.orb_names))}

        self.__ed = PomerolED(index_converter, verbose)

    def solve(self, H_int, E_levels, n_iw, const_of_motion=None, density_matrix_cutoff=1e-10, file_quantum_numbers="", file_eigenvalues=""):

        self.n_iw = n_iw
        self.__copy_E_levels(E_levels)

        H_0 = sum( self.E_levels[sn][o1,o2].real*c_dag(sn,o1)*c(sn,o2) for sn, o1, o2 in product(self.spin_names, self.orb_names, self.orb_names))
        if self.verbose and mpi.is_master_node():
            print "\n*** compute G_iw and Sigma_iw"
            print "*** n_iw =", n_iw
            print "*** E_levels =", E_levels
            print "*** H_0 =", H_0
            print "*** const_of_motion =", const_of_motion

        if const_of_motion is None:
            self.__ed.diagonalize(H_0+H_int)
        else:
            self.__ed.diagonalize(H_0+H_int, const_of_motion)

        # save data
        if file_quantum_numbers:
            self.__ed.save_quantum_numbers(file_quantum_numbers)
        if file_eigenvalues:
            self.__ed.save_eigenvalues(file_eigenvalues)

        # set density-matrix cutoff
        self.__ed.set_density_matrix_cutoff(density_matrix_cutoff)

        # Compute G(i\omega)
        self.G_iw = self.__ed.G_iw(self.gf_struct, self.beta, self.n_iw)
        # print type(self.G_iw)

        # Compute G0 and Sigma
        self.G0_iw = self.G_iw.copy()
        self.Sigma_iw = self.G_iw.copy()

        E_list = [ self.E_levels[sn] for sn in self.spin_names ]  # from dict to list
        self.G0_iw <<= iOmega_n
        self.G0_iw -= E_list
        self.Sigma_iw <<= self.G0_iw - inverse(self.G_iw)
        self.G0_iw.invert()

    def __copy_E_levels(self, E_levels):
        """Copy E_levels after checking the data structure"""

        # check data structure
        assert( isinstance(E_levels, dict) )
        for key, value in E_levels.items():
            assert( isinstance(value, np.ndarray) )
            assert( value.shape == (len(self.orb_names),len(self.orb_names)) )
        # copy
        self.E_levels = E_levels.copy()

    def __analyze_gf_struct(self, gf_struct):
        assert( isinstance(gf_struct, dict))
        self.spin_names = gf_struct.keys()
        self.orb_names = gf_struct.values()[0]
        for value in gf_struct.values():
            assert( self.orb_names == value )  # check if all spin blocks have the same orbitals
