
#pragma once

#include "gmp.hpp"
#include "prime_factorization.hpp"

#include <string>
#include <tuple>
#include <vector>

template <class Field> struct Polynomial {
  using element_t = typename Field::element_t;

  const Field &field;
  char variable;
  const element_t zero, one;
  std::vector<element_t> coeffs;
  std::vector<size_t> support; // k such that coeffs[k] != zero

  // Support provided
  Polynomial(const Field &field, const char variable,
             const std::vector<element_t> &coeffs,
             const std::vector<size_t> &support)
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs), support(support) {
    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();
  }

  // No support provided
  Polynomial(const Field &field, const char variable,
             const std::vector<element_t> &coeffs)
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs) {
    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();

    // Compute support
    for (size_t i = 0; i < this->coeffs.size(); ++i) {
      if (this->coeffs[i] != zero)
        this->support.push_back(i);
    }
  }

public:
  Polynomial(const Field &field, const char variable,
             const std::vector<integer_t> &coeffs = {})
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs.size(), zero),
        support() {
    if (coeffs.empty()) {
      this->coeffs = {zero, one};
      support = {0, 1};
      return;
    }

    for (size_t i = 0; i < coeffs.size(); ++i)
      this->coeffs[i] = field(coeffs[i]);

    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();

    // Compute support
    for (size_t i = 0; i < this->coeffs.size(); ++i) {
      if (this->coeffs[i] != zero)
        support.push_back(i);
    }
  }

  Polynomial(const Polynomial &other)
      : field(other.field), variable(other.variable), zero(other.zero),
        one(other.one), coeffs(other.coeffs), support(other.support) {}

  Polynomial(Polynomial &&other)
      : field(other.field), variable(other.variable), zero(other.zero),
        one(other.one), coeffs(std::move(other.coeffs)),
        support(std::move(other.support)) {}

  Polynomial &operator=(const Polynomial &other) {
    if (field != other.field)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different field");
    if (variable != other.variable)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different variable");
    coeffs = other.coeffs;
    support = other.support;
    return *this;
  }

  Polynomial &operator=(Polynomial &&other) {
    if (field != other.field)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different field");
    if (variable != other.variable)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different variable");
    std::swap(coeffs, other.coeffs);
    std::swap(support, other.support);
    return *this;
  }

  inline size_t degree() const { return coeffs.size() - 1; }
  inline element_t operator[](const size_t idx) const {
    return idx < coeffs.size() ? coeffs[idx] : zero;
  }
  inline Polynomial zero_poly() const {
    return Polynomial(field, variable, {zero}, {});
  }
  inline Polynomial one_poly() const {
    return Polynomial(field, variable, {one}, {0});
  }
  static inline Polynomial sample(const Field &field, const char variable,
                                  const uint64_t degree) {
    std::vector<element_t> coeffs(degree + 1, field.element(field.zero()));
    for (size_t i = 0; i <= degree; ++i)
      coeffs[i] = field.random_element();
    return Polynomial(field, variable, coeffs);
  }

  inline Polynomial monic() const {
    return *this == zero_poly() ? *this : (*this / coeffs.back());
  }

  Polynomial operator-() const {
    const size_t degree_plus_one = coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, zero);
    // TOOD: Parallelize this
    for (size_t i : support)
      result_coeffs[i] = -coeffs[i];
    return Polynomial(field, variable, result_coeffs, support);
  }

  friend Polynomial operator+(const Polynomial &a, const Polynomial &b) {
    if (a.field != b.field)
      throw math_error()
          << "Cannot add polynomials with elements from different fields";
    if (a.variable != b.variable)
      throw math_error()
          << "Cannot add polynomials with two different variables '"
          << a.variable << "' and '" << b.variable << "'";
    const size_t degree_plus_one = std::max(a.coeffs.size(), b.coeffs.size());
    std::vector<element_t> result_coeffs(degree_plus_one, a.zero);
    // TODO: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = a[i] + b[i];
    }
    return Polynomial(a.field, a.variable, result_coeffs);
  }

  friend Polynomial operator+(const Polynomial &p, const element_t k) {
    std::vector<element_t> result_coeffs = p.coeffs;
    result_coeffs[0] += k;
    return Polynomial(p.field, p.variable, result_coeffs);
  }
  friend Polynomial operator+(const element_t k, const Polynomial &p) {
    return p + k;
  }
  friend Polynomial operator+(const Polynomial &p, const integer_t k) {
    return p + p.field(k);
  }
  friend Polynomial operator+(const integer_t k, const Polynomial &p) {
    return p + k;
  }

  friend Polynomial operator-(const Polynomial &a, const Polynomial &b) {
    return a + (-b);
  }
  friend Polynomial operator-(const Polynomial &p, const element_t k) {
    std::vector<element_t> result_coeffs = p.coeffs;
    result_coeffs[0] -= k;
    return Polynomial(p.field, p.variable, result_coeffs);
  }
  friend Polynomial operator-(const element_t k, const Polynomial &p) {
    return -(p - k);
  }
  friend Polynomial operator-(const Polynomial &p, const integer_t k) {
    return p - p.field(k);
  }
  friend Polynomial operator-(const integer_t k, const Polynomial &p) {
    return -(p - k);
  }

  friend Polynomial operator*(const element_t k, const Polynomial &p) {
    if (k == p.zero)
      return p.zero_poly();
    // k is non-zero, so it doesn't change the support
    const size_t degree_plus_one = p.coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, p.zero);
    // TODO: Parallelize this
    for (const size_t i : p.support)
      result_coeffs[i] = k * p.coeffs[i];
    return Polynomial(p.field, p.variable, result_coeffs, p.support);
  }
  friend Polynomial operator*(const Polynomial &p, const element_t k) {
    return k * p;
  }
  friend Polynomial operator*(const integer_t k, const Polynomial &p) {
    return p * p.field(k);
  }
  friend Polynomial operator*(const Polynomial &p, const integer_t k) {
    return k * p;
  }
  Polynomial &operator*=(const element_t k) { return *this = *this * k; }
  Polynomial &operator*=(const integer_t k) { return *this = *this * k; }

  friend Polynomial operator*(const Polynomial &a, const Polynomial &b) {
    if (a.field != b.field)
      throw math_error()
          << "Cannot multiply polynomials with elements from different fields";
    if (a.variable != b.variable)
      throw math_error()
          << "Cannot multiply polynomials with two different variables '"
          << a.variable << "' and '" << b.variable << "'";
    // TODO: Write a sub-quadratic multiplication algorithm (maybe FFT?)
    const size_t a_degree = a.coeffs.size() - 1;
    const size_t b_degree = b.coeffs.size() - 1;
    const size_t result_degree = a_degree + b_degree;
    std::vector<element_t> result_coeffs(result_degree + 1, a.zero);
    for (const size_t i : a.support)
      for (const size_t j : b.support)
        result_coeffs[i + j] += a[i] * b[j];
    return Polynomial(a.field, a.variable, result_coeffs);
  }
  Polynomial &operator*=(const Polynomial &other) {
    return *this = *this * other;
  }

  friend Polynomial operator/(const Polynomial &p, const element_t k) {
    return p * k.inv();
  }

  Polynomial operator^(uint64_t exp) const {
    const Polynomial x = Polynomial(field, variable);
    if (*this == x) {
      std::vector<element_t> result_coeffs(exp + 1, zero);
      result_coeffs[exp] = one;
      return Polynomial(field, variable, result_coeffs, {exp});
    }

    Polynomial result(field, variable, {one}, {0}), pow = *this;
    while (exp > 0) {
      if (exp % 2 == 1)
        result *= pow;
      pow *= pow;
      exp /= 2;
    }
    return result;
  }

  friend Polynomial pow_mod(const Polynomial &f, integer_t exp,
                            const Polynomial &mod) {
    Polynomial result = f.one_poly(), base = f;
    // TODO: Mod out by the order q^n - 1
    while (exp != 0) {
      if (exp % 2 == 1)
        result = (result * base) % mod;
      base = (base * base) % mod;
      exp /= 2;
    }
    return result;
  }

  friend Polynomial inv_mod(const Polynomial &f, const Polynomial &mod) {
    if (f.field != mod.field)
      throw math_error() << "Cannot compute modular inverses for polynomials "
                            "over different fields";
    if (f.variable != mod.variable)
      throw math_error() << "Cannot compute modular inverses with polynomials "
                            "with two different variables '"
                         << f.variable << "' and '" << mod.variable << "'";

    const Polynomial zero = f.zero_poly(), one = f.one_poly();
    Polynomial r0 = mod, r1 = f, s0 = one, s1 = zero, t0 = zero, t1 = one;
    while (r1 != zero) {
      const auto &[q, r2] = divmod(r0, r1);
      std::tie(r0, s0, t0, r1, s1, t1) =
          std::make_tuple(r1, s1, t1, r2, s0 - s1 * q, t0 - t1 * q);
    }
    return (t0 + mod) % mod;
  }

  friend std::pair<Polynomial, Polynomial> divmod(const Polynomial &p,
                                                  const Polynomial &q) {
    if (p.field != q.field)
      throw math_error()
          << "Cannot divide polynomials with elements from different fields";
    if (p.variable != q.variable)
      throw math_error()
          << "Cannot divide polynomials with two different variables '"
          << p.variable << "' and '" << q.variable << "'";

    const size_t p_degree = p.degree(), q_degree = q.degree();
    if (q.coeffs.back() == 0)
      throw math_error() << "Division by zero: " << p << " / " << q;
    if (p_degree < q_degree)
      return std::make_pair(Polynomial(p.field, p.variable, {p.zero}), p);
    const size_t result_degree = p_degree - q_degree;
    const element_t q_leading_coeff = q.coeffs.back(),
                    inv_q_leading_coeff = q_leading_coeff.inv();
    std::vector<element_t> result_coeffs(result_degree + 1, p.zero);
    std::vector<element_t> p_coeffs = p.coeffs;
    for (int d = result_degree; d >= 0; --d) {
      if (p_coeffs[q_degree + d] == p.zero)
        continue;
      const element_t pivot = p_coeffs[q_degree + d] * inv_q_leading_coeff;
      result_coeffs[d] = pivot;
      for (const size_t i : q.support)
        p_coeffs[i + d] -= pivot * q.coeffs[i];
    }
    return std::make_pair(Polynomial(p.field, p.variable, result_coeffs),
                          Polynomial(p.field, p.variable, p_coeffs));
  }

  friend Polynomial operator/(const Polynomial &p, const Polynomial &q) {
    return divmod(p, q).first;
  }
  friend Polynomial operator%(const Polynomial &p, const Polynomial &q) {
    return divmod(p, q).second;
  }

  friend constexpr bool operator==(const Polynomial &a, const Polynomial &b) {
    return a.variable == b.variable && a.support == b.support &&
           a.coeffs == b.coeffs;
  }

  std::string to_string() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }

  friend std::ostream &operator<<(std::ostream &os, const Polynomial &p) {
    if (p.coeffs == std::vector<element_t>({p.zero}))
      return os << "0";
    bool first = true;
    for (auto d_it = p.support.rbegin(); d_it != p.support.rend(); ++d_it) {
      const size_t d = *d_it;
      if (!first)
        os << " + ";
      first = false;

      const element_t coeff = p.coeffs[d];
      if (d != 0) {
        if (coeff != p.one)
          os << coeff;
        os << p.variable;
      } else {
        os << coeff;
      }
      if (d > 1)
        os << "^" << d;
    }
    return os;
  }

  friend Polynomial polynomial_gcd(const Polynomial &f, const Polynomial &g) {
    if (f.field != g.field)
      throw math_error()
          << "Cannot compute gcd of polynomials over different fields";
    if (f.variable != g.variable)
      throw math_error()
          << "Cannot compute gcd of polynomials with two different variables '"
          << f.variable << "' and '" << g.variable << "'";

    Polynomial a = f, b = g;
    const Polynomial zero_poly(f.field, f.variable,
                               {f.field.element(f.field.zero())});
    while (true) {
      if (b == zero_poly)
        return a.monic();
      a = a % b;
      if (a == zero_poly)
        return b.monic();
      b = b % a;
    }
  }

  bool is_irreducible_rabin() const {
    if (*this == zero_poly())
      return false;
    const Polynomial f = monic();
    const integer_t q = f.field.cardinality();
    const uint64_t n = f.degree();

    // 1. Find the prime factors of n
    // 2. For each prime factor pi
    //    a. Compute ni = n / pi
    //    b. Compute h = (x^(q^ni) - x) mod f
    //    c. If gcd(f, h) != 1, return false
    // 3. If (x^(q^n) - x) mod f = 0, return true
    // 4. Otherwise, return false

    const Polynomial x = Polynomial(f.field, f.variable);
    for (const Factor &factor : prime_factor(gmp::from_uint(n))) {
      const uint64_t pi = gmp::to_uint(factor.base), ni = n / pi;
      const Polynomial h = (pow_mod(x, gmp::pow(q, ni), f) + f - x) % f;
      const Polynomial g = polynomial_gcd(f, h);
      if (g != f.one_poly())
        return false;
    }
    const Polynomial h = (pow_mod(x, gmp::pow(q, n), f) + f - x) % f;
    return h == f.zero_poly();
  }

  // More efficient than computing the order and comparing to q - 1
  bool is_primitive(const Polynomial &f) const {
    if (*this == zero_poly())
      return false;
    const integer_t p = field.characteristic(), q = field.cardinality();
    const uint64_t k = field.degree(), n = f.degree();
    const integer_t full_order = gmp::pow(p, k * n) - 1;
    const Factorization factorization = factor_pk_minus_one(p, k * n);
    for (const Factor &factor : factorization) {
      if (pow_mod(*this, full_order / factor.base, f) == one_poly())
        return false;
    }
    return true;
  }
};
