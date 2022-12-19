
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

template <class PrimeField> struct Lattice;
extern std::map<integer_t, Lattice<SmallPrimeField>> small_prime_lattices;
extern std::map<integer_t, Lattice<MediumPrimeField>> medium_prime_lattices;
extern std::map<integer_t, Lattice<LargePrimeField>> large_prime_lattices;

template <class PrimeField> struct Lattice {
  using field_t = std::shared_ptr<LatticeField<PrimeField>>;

  const integer_t p;
  std::vector<field_t> fields;
  Lattice(const integer_t p) : p(p) { add_prime_field(); }

  field_t add_field(const uint64_t k, const std::string &_variable = "");

private:
  field_t add_prime_field();
  field_t add_prime_poly_field(const std::string &variable, const uint64_t k);

  field_t add_zech_field(const std::string &variable, const uint64_t k)
    requires(!std::is_same<PrimeField, LargePrimeField>::value);

public:
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

template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_prime_field() {
  const auto result = std::make_shared<PrimeField>(p);
  fields.push_back(result);
  return result;
}

template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_zech_field(const std::string &variable,
                                    const uint64_t k)
  requires(!std::is_same<PrimeField, LargePrimeField>::value)
{
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  const auto P = dynamic_cast<PrimeField &>(*prime_field);
  const auto result = std::make_shared<ZechField<PrimeField>>(P, variable, k);
  fields.push_back(result);
  return result;
}

template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_prime_poly_field(const std::string &variable,
                                          const uint64_t k) {
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  const auto P = dynamic_cast<PrimeField &>(*prime_field);
  const auto result =
      std::make_shared<PrimePolyField<PrimeField>>(P, variable, k);
  fields.push_back(result);
  return result;
}

// When creating a new field, there are a few options:
// 1. If the field is a prime field, then
//    a. If p < 2^16,                      -> SmallPrimeField
//    b. If p < 2^32,                      -> MediumPrimeField
//    c. Otherwise,                        -> LargePrimeField
// 2. If the cardinality is at most 2^20,  -> ZechField
// 3. If there exists a field of size p^l
//      with l a non-trivial factor of k,
//      and p^l <= 2^20, create a
//      ZechField S with degree l, and     - - > ZechField (recursively)
//      create a ZechPolyField of degree
//      (k/l) over it                      -> ZechPolyField
// 4. If p < 2^16                          -> PrimePolyField
// 5. Finally                              -> GeneralPolyField
template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_field(const uint64_t k, const std::string &_variable) {
  const std::string variable =
      _variable != "" ? _variable : "w" + std::to_string(k);
  assert(fields.size() > 0);

  log() << "add_field called on lattice with characteristic " << p
        << " and degree " << k << std::endl;

  // 1. Prime field
  if (k == 1)
    return fields[0];

  const integer_t q = gmp::pow(p, k);

  if constexpr (std::is_same<PrimeField, LargePrimeField>::value) {
    return add_prime_poly_field(variable, k);
  } else {
    // 2. If the cardinality is at most 2^20, ZechField
    if (q <= (1_mpz << 20))
      return add_zech_field(variable, k);

    // 3. Two-step optimized representation
    uint32_t best_ell = 1;
    for (uint32_t ell = 2; ell < k; ++ell) {
      if (gmp::pow(p, ell) > (1_mpz << 20))
        break;
      if (k % ell == 0)
        best_ell = ell;
    }

    if (best_ell == 1) {
      // None of the factors of n had p^ell <= 2^20, so return PrimePolyField
      return add_prime_poly_field(variable, k);
    } else {
      // Otherwise, create a two-step optimized representation
      assert(fields.size() > 0);
      const auto prime_field = fields[0];
      const auto P = dynamic_cast<PrimeField &>(*prime_field);
      const auto intermediate_zech = std::make_shared<ZechField<PrimeField>>(
          P, "_w" + std::to_string(best_ell), best_ell);
      fields.push_back(intermediate_zech);

      const auto result =
          std::make_shared<ZechPolyField<ZechField<PrimeField>>>(
              *intermediate_zech, variable, k / best_ell);
      fields.push_back(result);
      return result;
    }
  }
}

inline void print_lattices(std::ostream &os) {
  os << std::endl;
  os << "------- LATTICES -------" << std::endl;
  for (const auto &[p, lattice] : small_prime_lattices) {
    os << "Lattice for prime " << p << ": " << std::endl;
    os << lattice << std::endl << std::endl;
  }
  for (const auto &[p, lattice] : medium_prime_lattices) {
    os << "Lattice for prime " << p << ": " << std::endl;
    os << lattice << std::endl << std::endl;
  }
  for (const auto &[p, lattice] : large_prime_lattices) {
    os << "Lattice for prime " << p << ": " << std::endl;
    os << lattice << std::endl << std::endl;
  }
  os << "------------------------" << std::endl;
}

inline void FiniteField(const integer_t p, const uint64_t k = 1) {
  if (p < 0 || !gmp::is_prime(p))
    throw math_error() << "FiniteField expected positive prime, got " << p;
  if (k == 0)
    throw math_error() << "FiniteField expected positive degree, got " << k;
  if (p.fits_ushort_p()) {
    if (small_prime_lattices.count(p) == 0)
      small_prime_lattices.emplace(p, Lattice<SmallPrimeField>(p));
    small_prime_lattices.at(p).add_field(k);
  } else if (p.fits_uint_p()) {
    if (medium_prime_lattices.count(p) == 0)
      medium_prime_lattices.emplace(p, Lattice<MediumPrimeField>(p));
    medium_prime_lattices.at(p).add_field(k);
  } else {
    if (large_prime_lattices.count(p) == 0)
      large_prime_lattices.emplace(p, Lattice<LargePrimeField>(p));
    large_prime_lattices.at(p).add_field(k);
  }
}

// There are four elements stored in an explicit field embedding:
// 1. The map phi_{E -> F}, as a matrix relating the generators of E and F
//    over P
// 2. The relative generator for F over E
// 3. The vector-space isomorphism between E^{(d)} and F
// 4. The minimal polynomial of the relative generator
template <class EField, class FField> struct FieldEmbedding {
  using p_field_t = typename EField::prime_field_t;
  using e_field_t = EField;
  using f_field_t = FField;
  using p_field_element_t = typename p_field_t::element_t;
  using e_field_element_t = typename e_field_t::element_t;
  using f_field_element_t = typename f_field_t::element_t;

  const p_field_t &P;
  const e_field_t &E;
  const f_field_t &F;

  const size_t e, f, d;

  // 1. f x e matrix in P representing the map phi_{E -> F}
  Matrix<p_field_t> phi;
  // 2. An element alpha_{F / E} of F such that E[alpha] = F
  f_field_element_t alpha_FE;
  // 3. A vector space isomorphism N encoding the transformation E^{(d)} -> F
  Matrix<p_field_t> psi, psi_inverse;
  // 4. The minimal polynomial of alpha_{F/E} over E
  Polynomial<EField> f_FE;

  FieldEmbedding(const p_field_t &P, const EField &E, const FField &F)
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

template <class EField, class FField,
          class PField = typename EField::prime_field_t>
FieldEmbedding<EField, FField> Embed(const EField &E, const FField &F) {
  const PField &P = E.prime_field();
  FieldEmbedding<EField, FField> result(P, E, F);

  const auto d = result.d, e = result.e, f = result.f;

  // 1. Computing prime-field generators
  const auto alpha_E = E.generating_element();
  const auto alpha_F = F.generating_element();
  const auto lifted_f = E.f.lift_from_prime_field(F);
  log() << "Finding root [tau] of " << lifted_f << " in " << F << std::endl;
  const auto tau = find_root(lifted_f);
  log() << "alpha_E = " << alpha_E << std::endl;
  log() << "alpha_F = " << alpha_F << std::endl;
  log() << "tau = " << tau << std::endl;

  // 2. Compute the matrix phi_{E -> F}
  // a. First, compute the matrix B of powers of alpha_F
  auto B = Matrix(P, 0, f);
  auto alpha_F_pow = F.element(F.one());
  for (size_t i = 0; i < f; ++i, alpha_F_pow *= alpha_F) {
    B.add_row(F.to_vector(alpha_F_pow).data);
  }
  B = B.transpose().inverse();

  auto M = Matrix(P, e, f);
  auto pow = F.element(F.one());
  for (size_t i = 0; i < e; ++i, pow *= tau) {
    const auto coeffs = B * F.to_vector(pow);
    for (size_t j = 0; j < f; ++j)
      M[i][j] = coeffs[j];
  }
  result.phi = M.transpose();
  log() << M << std::endl;

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
  log() << N << std::endl;
  assert(N.rank() == f);

  result.psi = N.transpose();
  result.psi_inverse = result.psi.inverse();

  // 5. Compute the minimal polynomial of alpha_FE over E
  result.f_FE = result.compute_minimal_polynomial(alpha_FE);
  return result;
}

// inline void Embed(const std::shared_ptr<AbstractField> &E,
//                   const std::shared_ptr<AbstractField> &F) {}
