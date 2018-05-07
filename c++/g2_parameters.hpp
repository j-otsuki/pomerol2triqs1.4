#pragma once

#include <triqs/hilbert_space/fundamental_operator_set.hpp>
using triqs::hilbert_space::gf_struct_t;
using indices_t = triqs::hilbert_space::fundamental_operator_set::indices_t;
using four_indices_t = std::tuple<indices_t, indices_t, indices_t, indices_t>;
using three_freqs_t = std::vector< std::tuple<int, int, int> >;

namespace pomerol2triqs {

  // enum block_order_t { AABB, ABBA };
  enum channel_t { PP, PH, AllFermionic };

  // using g2_blocks_t = std::set<std::pair<std::string, std::string>>;

  struct g2_iw_freq_box_params_t {

    /// Block structure of GF
    gf_struct_t gf_struct;

    /// Inverse temperature
    double beta;

    /// Channel in which Matsubara frequency representation is defined.
    channel_t channel = PH;

    /// Number of bosonic and fermionic Matsubara frequencies.
    int n_b, n_f;

    /// indices of operators in TRIQS convention: (block_name, inner_index)
    indices_t index1, index2, index3, index4;

    // g2_iw_freq_box_params_t() {}
    // g2_iw_freq_box_params_t(gf_struct_t const &gf_struct, double beta) : gf_struct(gf_struct), beta(beta) {}
  };

  struct g2_iw_freq_vec_params_t {

    /// Block structure of GF
    gf_struct_t gf_struct;

    /// Inverse temperature
    double beta;

    /// Channel in which Matsubara frequency representation is defined.
    channel_t channel = PH;

    /// three frequencies (wb, wf1, wf2).
    three_freqs_t three_freqs;

    /// set of indices of four operators in TRIQS convention: (block_name, inner_index)*4
    // indices_t index1, index2, index3, index4;
    std::vector<four_indices_t> vec_four_indices;

    // g2_iw_freq_vec_params_t() {}
    // g2_iw_freq_vec_params_t(gf_struct_t const &gf_struct, double beta) : gf_struct(gf_struct), beta(beta) {}
  };

/*
  struct g2_iw_l_lp_params_t {

    /// Structure of G^2 blocks.
    gf_struct_t gf_struct;

    /// Inverse temperature
    double beta;

    /// Channel in which Matsubara frequency representation is defined.
    channel_t channel = PH;

    /// Order of block indices in the definition of G^2.
    block_order_t block_order = AABB;

    /// List of block index pairs of G^2 to measure.
    /// default: measure all blocks
    g2_blocks_t blocks = g2_blocks_t{};

    /// Number of bosonic Matsubara frequencies.
    int n_iw = 30;

    /// Number of Legendre coefficients.
    int n_l = 20;

    /// Maximum number of positive Matsubara frequencies in summation.
    int n_inu_sum = 500;

    /// Tolerance for Matsubara frequency summation.
    double inu_sum_tol = 1e-6;

    g2_iw_l_lp_params_t() {}
    g2_iw_l_lp_params_t(gf_struct_t const &gf_struct, double beta) : gf_struct(gf_struct), beta(beta) {}
  };
*/
}
