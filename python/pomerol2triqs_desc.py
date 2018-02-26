# Generated automatically using the command :
# c++2py.py ../c++/pomerol_ed.hpp -I../../../local/pomerol/include -I/usr/include/eigen3 -I../c++ -p -mpytriqs.applications.impurity_solvers.pomerol2triqs -o pomerol2triqs --moduledoc "TRIQS wrapper around Pomerol ED library"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.impurity_solvers.pomerol2triqs", doc = "TRIQS wrapper around Pomerol ED library", app_name = "pytriqs.applications.impurity_solvers.pomerol2triqs")

# All the triqs C++/Python modules
module.use_module('gf', 'triqs')
module.use_module('operators', 'triqs')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("pomerol_ed.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/set.hpp>
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/variant.hpp>
#include <triqs/python_tools/converters/variant_int_string.hpp>
#include <triqs/python_tools/converters/tuple.hpp>
#include <triqs/python_tools/converters/operators_real_complex.hpp>
#include <triqs/python_tools/converters/gf.hpp>
""")
module.add_using("namespace pomerol2triqs")
module.add_using("namespace Pomerol")

module.add_enum("spin",          ["down", "up"], "Pomerol", "Spin projection")
# module.add_enum("block_order_t", ["AABB", "ABBA"], "pomerol2triqs", "G^{(2)} block order")
# module.add_enum("channel_t",     ["PP", "PH", "AllFermionic"], "pomerol2triqs", "G^{(2)} channel")

# The class pomerol_ed
c = class_(
        py_type = "PomerolED",  # name of the python class
        c_type = "pomerol_ed",   # name of the C++ class
        doc = r"",   # doc of the C++ class
)

c.add_constructor("""(index_converter_t index_converter, bool verbose = false)""",
                  doc = """Create a new solver object """)

c.add_method("""void diagonalize (many_body_op_t hamiltonian, bool ignore_symmetries = false)""",
             doc = """Diagonalize Hamiltonian optionally employing conservation of N and S_z """)

c.add_method("""void diagonalize (many_body_op_t hamiltonian, std::vector<many_body_op_t> integrals_of_motion)""",
             doc = """Diagonalize Hamiltonian using provided integrals of motion """)

c.add_method("""void saveQuantumNumbers (std::string filename)""",
             doc = """Save quantum numbers and block size """)

c.add_method("""void saveEigenValues (std::string filename)""",
             doc = """Save all eigenvalues and corresponding quantum numbers """)

c.add_method("""block_gf<imfreq> G_iw (gf_struct_t gf_struct, double beta, int n_iw)""",
             doc = """Green\'s function in Matsubara frequencies """)

c.add_method("""block_gf<imtime> G_tau (gf_struct_t gf_struct, double beta, int n_tau)""",
             doc = """Green\'s function in imaginary time """)

c.add_method("""block_gf<refreq> G_w (gf_struct_t gf_struct, double beta, std::pair<double,double> energy_window, int n_w, double im_shift = 0)""",
             doc = """Retarded Green\'s function on real energy axis """)

module.add_class(c)

module.generate_code()
