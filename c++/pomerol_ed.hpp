/**
 * pomerol2triqs1.4
 *
 * Copyright (C) 2017 Igor Krivenko
 * Copyright (C) 2018 Junya Otsuki
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#include <utility>
#include <vector>
#include <memory>
#include <tuple>
#include <functional>

#include <pomerol.h>
#include <triqs/gfs.hpp>
#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
// #include <triqs/utility/optional_compat.hpp>
#include <triqs/mpi/boost.hpp>
#include <triqs/arrays.hpp>

#include "g2_parameters.hpp"

namespace pomerol2triqs {

#ifdef POMEROL_COMPLEX_MATRIX_ELEMENTS
  using h_scalar_t = std::complex<double>;
#else
  using h_scalar_t = double;
#endif

  // Operator with real or complex value
  using many_body_op_t = triqs::operators::many_body_operator_generic<h_scalar_t>;

  using namespace triqs::gfs;
  using triqs::hilbert_space::gf_struct_t;

  using indices_t         = triqs::hilbert_space::fundamental_operator_set::indices_t;
  using pomerol_indices_t = std::tuple<std::string, unsigned int, Pomerol::spin>;
  using index_converter_t = std::map<indices_t, pomerol_indices_t>;

  class pomerol_ed {

    triqs::mpi::communicator comm;

    const bool verbose;
    index_converter_t index_converter;
    Pomerol::Lattice bare_lattice;
    Pomerol::IndexClassification index_info;

    std::unique_ptr<Pomerol::Lattice> lattice;
    std::unique_ptr<Pomerol::IndexHamiltonian> storage;
    std::unique_ptr<Pomerol::Symmetrizer> symm;
    std::unique_ptr<Pomerol::StatesClassification> states_class;
    std::unique_ptr<Pomerol::Hamiltonian> matrix_h;
    std::unique_ptr<Pomerol::DensityMatrix> rho;
    std::set<Pomerol::ParticleIndex> computed_ops;
    std::unique_ptr<Pomerol::FieldOperatorContainer> ops_container;

    Pomerol::Lattice init(bool spin_orbit);
    Pomerol::ParticleIndex lookup_pomerol_index(indices_t const &i) const;
    std::set<Pomerol::ParticleIndex> gf_struct_to_pomerol_indices(gf_struct_t const &gf_struct) const;
    double diagonalize_prepare(many_body_op_t const &hamiltonian);
    void diagonalize_main(double gs_shift);
    void compute_rho(double beta);
    void compute_field_operators(gf_struct_t const &gf_struct);
    template <typename Mesh, typename Filler> block_gf<Mesh> compute_gf(gf_struct_t const &gf_struct, gf_mesh<Mesh> const &mesh, Filler filler) const;

    // using w_nu_nup_t = cartesian_product<imfreq, imfreq, imfreq>;
    // using w_l_lp_t   = cartesian_product<imfreq, legendre, legendre>;
    using g2_t = triqs::arrays::array<std::complex<double>, 3>;
    using g2_three_freqs_t = std::vector<std::complex<double> >;
    // template <typename Mesh, typename Filler>
    // block2_gf<Mesh, tensor_valued<4>> compute_g2(gf_struct_t const &gf_struct, gf_mesh<Mesh> const &mesh, block_order_t block_order,
    //                                              g2_blocks_t const &g2_blocks, Filler filler) const;
    std::vector<std::complex<double> > compute_g2(gf_struct_t const &gf_struct, double beta, channel_t channel, indices_t index1, indices_t index2, indices_t index3, indices_t index4, three_freqs_t const &three_freqs);

    double density_matrix_cutoff = 1e-15;

    public:
    /// Create a new solver object
    pomerol_ed(index_converter_t const &index_converter, bool verbose = false, bool spin_orbit = false);

    /// Diagonalize Hamiltonian optionally employing conservation of N and S_z
    void diagonalize(many_body_op_t const &hamiltonian, bool ignore_symmetries = false);

    /// Diagonalize Hamiltonian using provided integrals of motion
    void diagonalize(many_body_op_t const &hamiltonian, std::vector<many_body_op_t> const& integrals_of_motion);

    /// Save quantum numbers and block size
    void save_quantum_numbers(const std::string &filename);

    /// Save all eigenvalues and corresponding quantum numbers
    void save_eigenvalues(const std::string &filename);

    /// Set density-matrix cutoff (default: 1e-15, 0 = No cutoff). Density matrix will be recomputed.
    void set_density_matrix_cutoff(const double cutoff);

    /// Green's function in Matsubara frequencies
    block_gf<imfreq> G_iw(gf_struct_t const &gf_struct, double beta, int n_iw);

    /// Green's function in imaginary time
    block_gf<imtime> G_tau(gf_struct_t const &gf_struct, double beta, int n_tau);

    /// Retarded Green's function on real energy axis
    block_gf<refreq> G_w(gf_struct_t const &gf_struct, double beta, std::pair<double, double> const &energy_window, int n_w, double im_shift = 0);

    /// Two-particle Green's function, Matsubara frequencies
    TRIQS_WRAP_ARG_AS_DICT
    g2_t G2_iw(g2_iw_inu_inup_params_t const &p);

    TRIQS_WRAP_ARG_AS_DICT
    g2_three_freqs_t G2_iw_three_freqs(g2_three_freqs_params_t const &p);

    /// Two-particle Green's function, Matsubara frequencies
    // TRIQS_WRAP_ARG_AS_DICT
    // block2_gf<w_nu_nup_t, tensor_valued<4>> G2_iw_inu_inup(g2_iw_inu_inup_params_t const &p);

    /// Two-particle Green's function, bosonic Matsubara frequency + Legendre coefficients
    // TRIQS_WRAP_ARG_AS_DICT
    // block2_gf<w_l_lp_t, tensor_valued<4>> G2_iw_l_lp(g2_iw_l_lp_params_t const &p);
  };

  inline void pomerol_ed::set_density_matrix_cutoff(const double cutoff){
    this->density_matrix_cutoff = cutoff;
    rho.release();
  }
}
