
#pragma once

#include "error.hpp"
#include "field.hpp"
#include "gmp.hpp"
#include "polynomial.hpp"
#include "small_prime_field.hpp"
#include <sstream>
#include <vector>

constexpr std::string_view irreducible_polynomial_file =
    "tables/irreducible_polynomials.txt";

inline std::map<std::pair<uint64_t, uint64_t>, std::vector<uint64_t>>
    irreducible_polynomials;

inline void load_cached_irreducible_polynomials() {
  std::ifstream ifs(irreducible_polynomial_file);
  std::string line;
  while (std::getline(ifs, line)) {
    std::istringstream iss(line);
    uint64_t p, degree;
    iss >> p >> degree;
    std::vector<uint64_t> coeffs(degree + 1);
    for (size_t i = 0; i <= degree; ++i)
      iss >> coeffs[i];
    const auto key = std::make_pair(p, degree);
    irreducible_polynomials[key] = coeffs;
  }
  log() << "Finished loading " << irreducible_polynomials.size()
        << " irreducible polynomials from cache" << std::endl;
}

template <class Field>
inline void write_irreducible_polynomial(const Polynomial<Field> &f) {
  static std::ofstream ofs(irreducible_polynomial_file, std::ios::app);

  assert(f.field.degree() == 1);
  const uint64_t p = gmp::to_uint(f.field.characteristic()),
                 degree = f.degree();
  const auto key = std::make_pair(p, degree);
  if (irreducible_polynomials.count(key) > 0)
    return;
  std::vector<uint64_t> coeffs;
  for (size_t i = 0; i <= degree; ++i)
    coeffs.push_back(f.coeffs[i].value);
  irreducible_polynomials[key] = coeffs;

  ofs << p << " " << degree;
  for (size_t i = 0; i <= degree; ++i)
    ofs << " " << coeffs[i];
  ofs << std::endl;
}

template <class Field>
Polynomial<Field> cache_to_polynomial(const Field &field,
                                      const std::string &variable,
                                      const std::vector<uint64_t> coeffs) {
  std::vector<typename Field::element_t> result_coeffs;
  for (size_t i = 0; i < coeffs.size(); ++i) {
    result_coeffs.push_back(field(gmp::from_uint(coeffs[i])));
  }
  return Polynomial(field, variable, result_coeffs);
}

template <class Field>
static inline Polynomial<Field>
get_irreducible_polynomial(const Field &field, const std::string &variable,
                           const uint64_t degree) {
  if (field.degree() == 1) {
    const auto key =
        std::make_pair(gmp::to_uint(field.characteristic()), degree);
    if (irreducible_polynomials.count(key) > 0)
      return cache_to_polynomial(field, variable, irreducible_polynomials[key]);
  }

  while (true) {
    const Polynomial<Field> &f =
        random_polynomial<true, true>(field, variable, degree);
    if (f.is_irreducible_rabin())
      return f;
  }
}

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
