
#pragma once

#include "error.hpp"
#include "field.hpp"
#include "gmp.hpp"
#include "large_prime_field.hpp"
#include "matrix.hpp"
#include "medium_prime_field.hpp"
#include "prime_poly.hpp"
#include "small_prime_field.hpp"
#include "vector.hpp"
#include "zech_field.hpp"
#include "zech_poly_field.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

struct Lattice;
extern std::map<integer_t, Lattice> lattices;

struct Lattice {
  const integer_t p;
  std::vector<std::shared_ptr<AbstractField>> fields;
  Lattice(const integer_t p) : p(p) { add_prime_field(); }

  std::shared_ptr<AbstractField> add_prime_field();
  std::shared_ptr<AbstractField> add_zech_field(const std::string &variable,
                                                const uint64_t k);
  std::shared_ptr<AbstractField>
  add_prime_poly_field(const std::string &variable, const uint64_t k);
  std::shared_ptr<AbstractField>
  add_zech_poly_field(AbstractField *base, const std::string &variable,
                      const uint64_t k);

  std::shared_ptr<AbstractField> add_field(const uint64_t k,
                                           const std::string &_variable = "");

  friend std::ostream &operator<<(std::ostream &os, const Lattice &lattice) {
    for (const auto &field : lattice.fields) {
      os << " - " << field << ": " << std::setw(20) << std::setfill(' ')
         << field_type_to_string(field->type()) << ": ("
         << field->characteristic() << "^" << field->degree() << ") "
         << field->to_string() << std::endl;
    }
    return os;
  }
};

inline void print_lattices(std::ostream &os) {
  os << std::endl;
  os << "------- LATTICES -------" << std::endl;
  for (const auto &[p, lattice] : lattices) {
    os << "Lattice for prime " << p << ": " << std::endl;
    os << lattice << std::endl << std::endl;
  }
  os << "------------------------" << std::endl;
}

inline std::shared_ptr<AbstractField> FiniteField(const integer_t p,
                                                  const uint64_t k = 1) {
  if (p < 0 || !gmp::is_prime(p))
    throw math_error() << "FiniteField expected positive prime, got " << p;
  if (k == 0)
    throw math_error() << "FiniteField expected positive degree, got " << k;
  if (lattices.count(p) == 0)
    lattices.emplace(p, Lattice(p));
  return lattices.at(p).add_field(k);
}

// There are four elements stored in an explicit field embedding:
// 1. The map phi_{E -> F}, as a matrix relating the generators of E and F
//    over P
// 2. The relative generator for F over E
// 3. The vector-space isomorphism between E^{(d)} and F
// 4. The minimal polynomial of the relative generator
template <class PField, class EField, class FField> struct FieldEmbedding {
  using p_field_t = PField;
  using e_field_t = EField;
  using f_field_t = FField;
  using p_field_element_t = typename PField::element_t;
  using e_field_element_t = typename EField::element_t;
  using f_field_element_t = typename FField::element_t;

  const PField &P;
  const EField &E;
  const FField &F;

  const size_t e, f, d;

  // 1. f x e matrix in P representing the map phi_{E -> F}
  Matrix<PField> phi;
  // 2. An element alpha_{F / E} of F such that E[alpha] = F
  f_field_element_t alpha_FE;
  // 3. A vector space isomorphism N encoding the transformation E^{(d)} -> F
  Matrix<PField> psi, psi_inverse;
  // 4. The minimal polynomial of alpha_{F/E} over E
  Polynomial<EField> f_FE;

  FieldEmbedding(const PField &P, const EField &E, const FField &F)
      : P(P), E(E), F(F), e(E.degree()), f(F.degree()), d(f / e), phi(P),
        alpha_FE(F.element(F.zero())), psi(P), psi_inverse(P),
        f_FE(E, "w", {E.element(E.zero())}, {}) {
    if (f % e != 0)
      throw math_error() << "Could not create embedding for incompatible "
                            "fields E and F with absolute degrees "
                         << e << " and " << f << ", respectively";
  }

  f_field_element_t apply_embedding(const e_field_element_t &elem) const {
    return F.from_vector(phi * E.to_vector(elem));
  }

  Polynomial<f_field_t> lift(const Polynomial<e_field_t> &f) const {
    std::vector<f_field_element_t> result_coeffs;
    for (const auto &coeff : f.coeffs)
      result_coeffs.push_back(apply_embedding(coeff));
    return Polynomial<f_field_t>(F, f.variable, result_coeffs, f.support);
  }

  f_field_element_t from_E_vector(const Vector<e_field_t> &vec) const {
    // 1. Apply (\phi_E^{-1})^(d) to vec
    std::vector<p_field_element_t> components;
    components.reserve(f);
    for (const e_field_element_t &component : vec.data) {
      const auto p_vector = E.to_vector(component);
      components.insert(components.end(), p_vector.data.begin(),
                        p_vector.data.end());
    }
    const auto in_vector = Vector(P, f, components);

    // 2. Apply the linear transformation in psi
    const auto out_vector = psi * in_vector;

    // 3. Apply phi_F
    return F.from_vector(out_vector);
  }

  Vector<e_field_t> to_E_vector(const f_field_element_t &elem) const {
    // 1. Apply phi_F^{-1}
    const auto out_vector = F.to_vector(elem);

    // 2. Apply the inverse of psi
    const auto in_vector = psi_inverse * out_vector;

    // 3. Apply (\psi_E)^{(d)}
    std::vector<e_field_element_t> result_coeffs;
    result_coeffs.reserve(d);
    for (size_t i = 0; i < d; ++i) {
      std::vector<p_field_element_t> component;
      component.reserve(e);
      component.insert(component.end(), in_vector.data.begin() + i * e,
                       in_vector.data.begin() + (i + 1) * e);
      const auto E_component = E.from_vector(Vector(P, e, component));
      result_coeffs.push_back(E_component);
    }
    return Vector(E, d, result_coeffs);
  }

  Polynomial<EField> compute_minimal_polynomial(const f_field_element_t &beta) {
    auto B = Matrix(E, 0, d);
    auto pow = F.element(F.one());
    for (size_t i = 0; i <= d; ++i, pow *= beta) {
      // Invariant: pow = beta^i

      // 1. Apply the inverse vector space isomorphism to get a vector in E^d
      const auto this_row = to_E_vector(pow);

      // 2. Append the vector to the matrix B and check its rank
      B.add_row(this_row.data);

      // 3. If the matrix is no longer full-rank, break, and compute the
      // resulting linear dependence and thus the minimal polynomial
      if (B.rank() != B.rows)
        break;
    }

    const auto degree_plus_one = B.rows;
    const auto neg_beta_k = -B.pop_row();
    const auto coeffs = B.transpose().solve(neg_beta_k);
    std::vector<e_field_element_t> result_coeffs = coeffs.data;
    result_coeffs.push_back(E.element(E.one()));
    assert(result_coeffs.size() == degree_plus_one);
    const auto result = Polynomial(E, "w", result_coeffs);
    assert(lift(result)(beta) == 0);
    return result;
  }
};

template <class PField, class EField, class FField>
FieldEmbedding<PField, EField, FField> Embed(const PField &P, const EField &E,
                                             const FField &F) {
  FieldEmbedding<PField, EField, FField> result(P, E, F);

  const auto d = result.d, e = result.e, f = result.f;

  // 1. Computing prime-field generators
  const auto alpha_E = E.generating_element();
  const auto alpha_F = F.generating_element();
  const auto tau = find_root(E.f.lift_from_prime_field(F));
  log() << "alpha_E = " << alpha_E << std::endl;
  log() << "alpha_F = " << alpha_F << std::endl;
  log() << "tau = " << tau << std::endl;

  if (F.g != F.g.var_poly())
    throw math_error("Currently only supports fields generated by w");

  // 2. Computing the matrix phi_{E -> F}
  auto M = Matrix(P, e, f);
  for (size_t i = 0; i < e; ++i) {
    const auto tau_i = F.to_polynomial((tau ^ i).value);
    for (size_t j = 0; j < f; ++j)
      M[i][j] = tau_i[j];
  }
  result.phi = M.transpose();
  // log() << M << std::endl;

  // 3. Finding a generator for F/E
  // If G is contained in E and we already have a generator for F over
  // G, then use that
  // TODO: Implement the above more generally once we have embeddings
  // For now, we just restrict ourselves to fields F with G = P
  const auto alpha_FE = F.generating_element();
  result.alpha_FE = alpha_FE;

  // 4. Compute the vector-space isomorphism E^{(d)} -> F
  Matrix N(P, f, f);
  for (size_t i = 0; i < d; ++i) {
    for (size_t j = 0; j < e; ++j) {
      // Compute \phi_F^{-1}(alpha_FE^i * tau^j)
      const auto alpha_FE_to_i = alpha_FE ^ i;
      const auto tau_to_j = tau ^ j;
      const auto product = alpha_FE_to_i * tau_to_j;
      const auto row = F.to_vector(product);
      N[i * e + j] = row.data;
    }
  }
  result.psi = N.transpose();
  result.psi_inverse = result.psi.inverse();

  // log() << N << std::endl;
  assert(N.rank() == f);

  // 5. Compute the minimal polynomial of alpha_FE over E
  result.f_FE = result.compute_minimal_polynomial(alpha_FE);
  return result;
}

// inline void Embed(const std::shared_ptr<AbstractField> &E,
//                   const std::shared_ptr<AbstractField> &F) {}
