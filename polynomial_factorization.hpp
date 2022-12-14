
#pragma once

#include "error.hpp"
#include "field.hpp"
#include "polynomial.hpp"
#include <vector>

template <class Field> struct PolynomialFactor {
  Polynomial<Field> base;
  uint64_t exp;
  PolynomialFactor(const Polynomial<Field> &base, const uint64_t exp)
      : base(base), exp(exp) {}
};

template <class Field>
inline void verify_polynomial_factorization(
    const Polynomial<Field> &f,
    const std::vector<PolynomialFactor<Field>> &factorization) {
  Polynomial<Field> product = f.one_poly();
  for (const auto &factor : factorization) {
    assert(factor.base.is_irreducible_rabin());
    product *= factor.base ^ factor.exp;
  }
  assert(product == f);
}

template <class Field>
std::vector<PolynomialFactor<Field>>
equal_degree_factorization(const Polynomial<Field> &f, const uint64_t d) {
  if (f.degree() % d != 0)
    throw math_error() << "Expected f (" << f
                       << ") to consist of irreducible factors with degree d = "
                       << d << ", but d did not divide the degree of f ("
                       << f.degree() << ")";

  // TODO: Check that f is square-free
  if (f.degree() == d)
    return {PolynomialFactor(f, 1)};

  // Computes (h^1 + h^2 + h^4 + ... + h^{2^{d-1}}) % f, useful when p = 2
  const auto compute_squared_sum = [](const Polynomial<Field> &h,
                                      const uint64_t d,
                                      const Polynomial<Field> &f) {
    Polynomial<Field> result = h.zero_poly(), base = h;
    for (size_t k = 0; k < d; ++k) {
      result = (result + base) % f;
      base = (base * base) % f;
    }
    return result;
  };

  // Cantor-Zassenhaus
  // 1. Maintain a set of factors, initialized with f
  // 2. While we have less than r factors,
  //    a. Sample an element in Field[x] with degree < n = deg(f) at random
  //    b. Compute g = pow_mod(h, (q^d - 1) / 2, f) - 1
  //        i. For each factor u in factors with degree > n,
  //        ii. Compute gcd(g, u)
  //        iii. If this is non-trivial, break the factor
  const bool odd_characteristic = f.field.characteristic() % 2 == 1;
  std::vector<Polynomial<Field>> factors = {f};
  const uint64_t n = f.degree(), r = n / d;
  const integer_t q = f.field.cardinality(), exp = (gmp::pow(q, d) - 1) / 2;
  while (factors.size() < r) {
    std::vector<Polynomial<Field>> new_factors;
    const Polynomial<Field> h =
        random_polynomial<false, false>(f.field, f.variable, n - 1);
    const Polynomial<Field> g = odd_characteristic
                                    ? (pow_mod(h, exp, f) - 1) % f
                                    : compute_squared_sum(h, d, f);
    for (const Polynomial<Field> &u : factors) {
      if (u.degree() <= d) {
        new_factors.push_back(u);
        continue;
      }
      const Polynomial<Field> gcd = polynomial_gcd(g, u);
      if (gcd == u || gcd == u.one_poly()) {
        new_factors.push_back(u);
      } else {
        new_factors.push_back(gcd);
        new_factors.push_back(u / gcd);
      }
    }
    factors = new_factors;
  }

  std::vector<PolynomialFactor<Field>> result;
  for (const auto &factor : factors)
    result.emplace_back(factor, 1);
  verify_polynomial_factorization(f, result);
  return result;
}

// Given a polynomial f with degree >= 2 which splits into linear factors over a
// Field F, return an arbitrary root alpha of f in F
template <class Field>
inline typename Field::element_t find_root(const Polynomial<Field> &f) {
  // Idea: Basically a version of Cantor-Zassenhaus above, specialized to d = 1,
  // which exits after finding any linear factor

  const bool odd_characteristic = f.field.characteristic() % 2 == 1;
  std::vector<Polynomial<Field>> factors = {f};
  const uint64_t n = f.degree();
  const integer_t q = f.field.cardinality(), exp = (q - 1) / 2;
  while (!factors.empty()) {
    std::vector<Polynomial<Field>> new_factors;
    const Polynomial<Field> h =
        random_polynomial<false, false>(f.field, f.variable, n - 1);
    const Polynomial<Field> g =
        odd_characteristic ? (pow_mod(h, exp, f) - 1) % f : h;
    for (const Polynomial<Field> &u : factors) {
      if (u.degree() == 1) {
        const auto root = -u[0];
        assert(f.at(root) == 0);
        return root;
      }
      const Polynomial<Field> gcd = polynomial_gcd(g, u);
      if (gcd == u || gcd == u.one_poly()) {
        new_factors.push_back(u);
      } else {
        new_factors.push_back(gcd);
        new_factors.push_back(u / gcd);
      }
    }
    factors = new_factors;
  }
  __builtin_unreachable();
}

template <class Field>
inline void print_polynomial_factorization(
    const std::vector<PolynomialFactor<Field>> &factorization) {
  for (const auto &factor : factorization) {
    std::cout << "(" << factor.base << ")";
    if (factor.exp != 1)
      std::cout << "^" << factor.exp << " ";
  }
  std::cout << std::endl;
}
