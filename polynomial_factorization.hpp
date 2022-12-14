
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
std::vector<PolynomialFactor<Field>>
equal_degree_factorization_odd(const Polynomial<Field> &f, const uint64_t d) {
  if (f.field.characteristic() % 2 == 0)
    throw math_error() << "Expected coefficient field (" << f.field
                       << ") to have odd characteristic, but got "
                       << f.field.characteristic();
  if (f.degree() % d != 0)
    throw math_error() << "Expected f (" << f
                       << ") to consist of irreducible factors with degree d = "
                       << d << ", but d did not divide the degree of f ("
                       << f.degree() << ")";

  // TODO: Check that f is square-free
  if (f.degree() == d)
    return {PolynomialFactor(f, 1)};

  // Cantor-Zassenhaus
  // 1. Maintain a set of factors, initialized with f
  // 2. While we have less than r factors,
  //    a. Sample an element in Field[x] with degree < n = deg(f) at random
  //    b. Compute g = pow_mod(h, (q^d - 1) / 2, f) - 1
  //        i. For each factor u in factors with degree > n,
  //        ii. Compute gcd(g, u)
  //        iii. If this is non-trivial, break the factor
  std::vector<Polynomial<Field>> factors = {f};
  const uint64_t n = f.degree(), r = n / d;
  const integer_t q = f.field.cardinality(), exp = (gmp::pow(q, d) - 1) / 2;
  while (factors.size() < r) {
    std::vector<Polynomial<Field>> new_factors;
    const Polynomial<Field> h =
        random_polynomial<false, false>(f.field, f.variable, n - 1);
    const Polynomial<Field> g = (pow_mod(h, exp, f) + f - 1) % f;
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
  for (const auto &factor : factors) {
    result.emplace_back(factor, 1);
  }
  return result;
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
