
#pragma once

#include "error.hpp"
#include "field.hpp"
#include "gmp.hpp"
#include "polynomial.hpp"
#include "small_prime_field.hpp"
#include <cstddef>
#include <sstream>
#include <vector>

const inline std::string irreducible_polynomial_file =
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

template <class Field> struct PolynomialFactorization {
  std::vector<PolynomialFactor<Field>> factors;

  PolynomialFactorization() = default;
  explicit PolynomialFactorization(const Polynomial<Field> &f)
      : factors(1, PolynomialFactor<Field>(f, 1)) {}

  void combine(const PolynomialFactorization &other) {
    // TODO: Combine like terms
    factors.insert(factors.end(), other.factors.begin(), other.factors.end());
  }

  PolynomialFactorization pow(const size_t exp) const {
    PolynomialFactorization result;
    for (const auto &factor : factors)
      result.factors.emplace_back(factor.base, factor.exp * exp);
    return result;
  }

  void emplace_back(const Polynomial<Field> &base, const uint64_t exp) {
    factors.emplace_back(base, exp);
  }

  Polynomial<Field> product() const {
    assert(!factors.empty());
    Polynomial<Field> product = factors[0].base.one_poly();
    for (const auto &factor : factors)
      product *= factor.base ^ factor.exp;
    return product;
  }

  bool all_irreducible() const {
    for (const auto &factor : factors)
      if (!factor.base.is_irreducible_rabin())
        return false;
    return true;
  }

  void verify(const Polynomial<Field> &f) { assert(product() == f); }

  auto begin() const { return factors.begin(); }
  auto end() const { return factors.end(); }

  friend std::ostream &
  operator<<(std::ostream &os,
             const PolynomialFactorization<Field> &factorization) {
    for (const auto &factor : factorization) {
      os << "(" << factor.base << ")";
      if (factor.exp != 1)
        os << "^" << factor.exp << " ";
    }
    return os;
  }
};

template <class Field>
PolynomialFactorization<Field>
square_free_factorization(const Polynomial<Field> &f) {
  const integer_t p = f.field.characteristic();
  const uint64_t p_uint = gmp::to_uint(p);
  PolynomialFactorization<Field> result;

  Polynomial<Field> c = polynomial_gcd(f, f.derivative()), w = f / c;

  for (size_t i = 1; w != 1; ++i) {
    const Polynomial<Field> y = polynomial_gcd(w, c);
    const Polynomial<Field> factor = w / y;
    if (factor != 1)
      result.emplace_back(factor, i);
    w = y;
    c /= y;
  }

  if (c != 1) {
    assert(c.degree() % p_uint == 0);
    const size_t new_degree = c.degree() / p_uint;
    std::vector<typename Field::element_t> new_coeffs(
        new_degree + 1, f.field.element(f.field.zero()));
    for (const size_t d : c.support) {
      assert(d % p_uint == 0);
      new_coeffs[d / p_uint] = c.coeffs[d];
    }
    const Polynomial new_c = Polynomial(c.field, c.variable, new_coeffs);
    const auto factorization = square_free_factorization(new_c);
    result.combine(factorization.pow(p_uint));
  }
  result.verify(f);
  return result;
}

template <class Field>
std::vector<std::pair<Polynomial<Field>, size_t>>
distinct_degree_factorization(const Polynomial<Field> &f) {
  // TODO: Check that f is square-free
  std::vector<std::pair<Polynomial<Field>, size_t>> result;
  const integer_t q = f.field.cardinality();
  const Polynomial<Field> x = f.var_poly();
  size_t i = 1;
  Polynomial<Field> f1 = f;
  while (f1.degree() >= 2 * i) {
    const Polynomial<Field> g =
        polynomial_gcd(f1, pow_mod(x, gmp::pow(q, i), f) - x);
    if (g != 1) {
      result.emplace_back(g, i);
      f1 /= g;
    }
    ++i;
  }
  if (f1 != 1)
    result.emplace_back(f1, f1.degree());
  if (result.empty())
    return {std::make_pair(f, 1)};
  return result;
}

template <class Field>
PolynomialFactorization<Field>
equal_degree_factorization(const Polynomial<Field> &f, const uint64_t d) {
  if (f.degree() % d != 0)
    throw math_error() << "Expected f (" << f
                       << ") to consist of irreducible factors with degree d = "
                       << d << ", but d did not divide the degree of f ("
                       << f.degree() << ")";

  // TODO: Check that f is square-free
  if (f.degree() == d) {
    return PolynomialFactorization(f);
  }

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

  PolynomialFactorization<Field> result;
  for (const auto &factor : factors)
    result.emplace_back(factor, 1);
  assert(result.all_irreducible());
  result.verify(f);
  return result;
}

template <class Field>
PolynomialFactorization<Field> factor_polynomial(const Polynomial<Field> &f) {
  assert(f != 0);
  PolynomialFactorization<Field> result;
  const auto leading_coefficient = f.coeffs.back();
  if (leading_coefficient != 1)
    result.emplace_back(leading_coefficient * f.one_poly(), 1);
  const auto square_free_factors =
      square_free_factorization(f / leading_coefficient);
  log() << "Square-free decomposition: " << square_free_factors << std::endl;
  for (const auto &factor : square_free_factors) {
    const auto distinct_degree_factors =
        distinct_degree_factorization(factor.base);
    for (const auto &[g, d] : distinct_degree_factors) {
      const auto factorization = equal_degree_factorization(g, d);
      result.combine(factorization.pow(factor.exp));
    }
  }
  assert(result.all_irreducible());
  result.verify(f);
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
