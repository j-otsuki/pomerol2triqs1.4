
from pomerol2triqs import PomerolED
from pytriqs.gf.local import *
from pytriqs.operators import Operator, c, c_dag, n
from pytriqs.operators.util.op_struct import set_operator_structure, get_mkind
import pytriqs.utility.mpi as mpi
import numpy as np
from itertools import product

class Solver:
    def __init__(self, beta, gf_struct, n_iw, spin_orbit=False, verbose=True):

        self.beta = beta
        self.gf_struct = gf_struct
        self.n_iw = n_iw
        self.spin_orbit = spin_orbit
        self.verbose = verbose

        if self.verbose and mpi.is_master_node():
            print "\n*** Hubbard I solver using Pomerol library"
            print "*** gf_struct =", gf_struct
            print "*** n_iw =", n_iw

        # get spin_name and orb_names from gf_struct
        self.__analyze_gf_struct(gf_struct)

        if self.verbose and mpi.is_master_node():
            print "*** spin_names =", self.spin_names
            print "*** orb_names  =", self.orb_names

        # check spin_names
        for sn in self.spin_names:
            if not spin_orbit:  assert sn in ('down', 'dn', 'up')  # become either 'down' or 'up'
            if     spin_orbit:  assert sn in ('down', 'dn', 'ud')  # all become 'down'

        # Conversion from TRIQS to Pomerol notation for operator indices
        # TRIQS: ('dn', '-2') --> Pomerol: ('atom', 0, 'down')
        # NOTE: When spin_orbit is true, only spin 'down' is used.
        mkind = get_mkind(True, None)
        index_converter = {mkind(sn, bn) : ("atom", bi, "down" if sn == "dn" or sn == "ud" else sn)
                           for sn, (bi, bn) in product(self.spin_names, enumerate(self.orb_names))}

        self.__ed = PomerolED(index_converter, verbose, spin_orbit)

        # init G_iw
        glist = lambda : [ GfImFreq(indices=self.orb_names, beta=beta, n_points=n_iw) for block, inner in gf_struct.items() ]
        self.G_iw = BlockGf(name_list=self.spin_names, block_list=glist(), make_copies=False)
        self.G_iw.zero()
        self.G0_iw = self.G_iw.copy()
        self.Sigma_iw = self.G_iw.copy()

    def solve(self, H_int, E_levels, integrals_of_motion=None, density_matrix_cutoff=1e-10, file_quantum_numbers="", file_eigenvalues=""):

        self.__copy_E_levels(E_levels)

        H_0 = sum( self.E_levels[sn][oi1,oi2].real*c_dag(sn,on1)*c(sn,on2) for sn, (oi1,on1), (oi2,on2) in product(self.spin_names, enumerate(self.orb_names), enumerate(self.orb_names)) )
        if self.verbose and mpi.is_master_node():
            print "\n*** compute G_iw and Sigma_iw"
            print "*** E_levels =", E_levels
            print "*** H_0 =", H_0
            # print "*** integrals_of_motion =", integrals_of_motion

        N_op = sum( c_dag(sn,on)*c(sn,on) for sn, on in product(self.spin_names, self.orb_names) )

        if integrals_of_motion is None:
            if self.spin_orbit:
                # use only N op as integrals_of_motion when spin-orbit coupling is included
                self.__ed.diagonalize(H_0+H_int, integrals_of_motion=[N_op])
            else:
                # use N and S_z as integrals of motion by default
                self.__ed.diagonalize(H_0+H_int)
        else:
            assert isinstance(integrals_of_motion, list)

            # check commutation relation in advance (just for test)
            if self.verbose and mpi.is_master_node():
                print "*** Integrals of motion:"
                for i, op in enumerate(integrals_of_motion):
                    print " op%d =" %i, op

                print "*** commutation relations:"
                is_commute = lambda h, op: h*op - op*h
                str_if_commute = lambda h, op: "== 0" if is_commute(h, op).is_zero() else "!= 0"
                for i, op in enumerate(integrals_of_motion):
                    print " [H_0, op%d]" %i, str_if_commute(H_0, op), " [H_int, op%d]" %i, str_if_commute(H_int, op)

            self.__ed.diagonalize(H_0+H_int, integrals_of_motion)

        # save data
        if file_quantum_numbers:
            self.__ed.save_quantum_numbers(file_quantum_numbers)
        if file_eigenvalues:
            self.__ed.save_eigenvalues(file_eigenvalues)

        # set density-matrix cutoff
        self.__ed.set_density_matrix_cutoff(density_matrix_cutoff)

        # Compute G(i\omega)
        self.G_iw << self.__ed.G_iw(self.gf_struct, self.beta, self.n_iw)

        # Compute G0 and Sigma
        E_list = [ self.E_levels[sn] for sn in self.spin_names ]  # from dict to list
        self.G0_iw << iOmega_n
        self.G0_iw -= E_list
        self.Sigma_iw << self.G0_iw - inverse(self.G_iw)
        self.G0_iw.invert()

        # *********************************************************************
        # TODO: tail
        # set tail of Sigma at all zero in the meantime, because tail of G is
        #  not computed in pomerol solver. This part should be improved.
        # *********************************************************************
        for s, sig in self.Sigma_iw:
            sig.tail.zero()

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
