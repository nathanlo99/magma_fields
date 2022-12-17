
#pragma once

#include "gmp.hpp"
#include "logger.hpp"
#include "prime_factorization.hpp"

#include <string>
#include <tuple>
#include <vector>

template <class Field> struct Polynomial {
  using element_t = typename Field::element_t;

  const Field &field;
  const std::string variable;
  const element_t zero, one;
  std::vector<element_t> coeffs;
  std::vector<size_t> support; // k such that coeffs[k] != zero

  // Support provided
  Polynomial(const Field &field, const std::string &variable,
             const std::vector<element_t> &coeffs,
             const std::vector<size_t> &support)
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs), support(support) {
    remove_leading_zeroes();
  }

  // No support provided
  Polynomial(const Field &field, const std::string &variable,
             const std::vector<element_t> &coeffs)
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs) {
    remove_leading_zeroes();

    // Compute support
    for (size_t i = 0; i < this->coeffs.size(); ++i) {
      if (this->coeffs[i] != zero)
        this->support.push_back(i);
    }
  }

  void remove_leading_zeroes() {
    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();
  }

public:
  Polynomial(const Field &field, const std::string &variable,
             const std::vector<integer_t> &coeffs = {})
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs.size(), zero),
        support() {
    if (coeffs.empty()) {
      this->coeffs = {zero, one};
      support = {1};
      return;
    }

    for (size_t i = 0; i < coeffs.size(); ++i)
      this->coeffs[i] = field(coeffs[i]);

    remove_leading_zeroes();

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
    coeffs = std::move(other.coeffs);
    support = std::move(other.support);
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
  inline Polynomial var_poly() const {
    return Polynomial(field, variable, {zero, one}, {1});
  }

  bool is_monic() const { return coeffs.back() == 1; }
  inline Polynomial to_monic() const {
    if (*this == 0)
      return *this;
    return *this / coeffs.back();
  }
  inline Polynomial derivative() const {
    if (degree() == 0)
      return zero_poly();
    std::vector<element_t> new_coeffs(coeffs.size() - 1, zero);
    for (size_t i = 0; i + 1 < coeffs.size(); ++i) {
      new_coeffs[i] = (i + 1) * coeffs[i + 1];
    }
    return Polynomial(field, variable, new_coeffs);
  }

  template <class TargetField>
  Polynomial<TargetField> lift(const TargetField &target_field) const {
    if (field.degree() != 1)
      throw math_error()
          << "Can currently only lift polynomials from the base field";
    std::vector<typename TargetField::element_t> result_coeffs;
    for (size_t i = 0; i < coeffs.size(); ++i) {
      const integer_t value = coeffs[i].as_integer();
      result_coeffs.push_back(target_field(value));
    }
    return Polynomial<TargetField>(target_field, variable, result_coeffs,
                                   support);
  }

  Polynomial operator-() const {
    const size_t degree_plus_one = coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, zero);
    // TOOD: Parallelize this
    for (const size_t i : support)
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
    for (size_t i = 0; i < degree_plus_one; ++i)
      result_coeffs[i] = a[i] + b[i];
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
    if (a.field != b.field)
      throw math_error()
          << "Cannot subtract polynomials with elements from different fields";
    if (a.variable != b.variable)
      throw math_error()
          << "Cannot subtract polynomials with two different variables '"
          << a.variable << "' and '" << b.variable << "'";
    const size_t degree_plus_one = std::max(a.coeffs.size(), b.coeffs.size());
    std::vector<element_t> result_coeffs(degree_plus_one, a.zero);
    // TODO: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i)
      result_coeffs[i] = a[i] - b[i];
    return Polynomial(a.field, a.variable, result_coeffs);
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
    const size_t degree_plus_one = p.coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, p.zero);
    // TODO: Parallelize this
    for (const size_t i : p.support)
      result_coeffs[i] = k * p.coeffs[i];
    // k is non-zero, so it doesn't change the support
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
  Polynomial &operator*=(const element_t k) {
    if (k == zero)
      return *this = zero_poly();
    for (const size_t i : support)
      coeffs[i] *= k;
    // If k is non-zero, it doesn't change the support, no need to re-calculate
    // or remove leading zeroes
    return *this;
  }
  Polynomial &operator*=(const integer_t k) { return *this *= field(k); }

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
    if (k == p.zero)
      throw math_error("Division by zero");
    return p * k.inv();
  }

  Polynomial operator^(uint64_t exp) const {
    const Polynomial x = Polynomial(field, variable);
    if (*this == x) {
      std::vector<element_t> result_coeffs(exp + 1, zero);
      result_coeffs[exp] = one;
      return Polynomial(field, variable, result_coeffs, {exp});
    }

    Polynomial result = one_poly(), pow = *this;
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
    while (exp != 0) {
      if (exp % 2 == 1)
        result = (result * base) % mod;
      base = (base * base) % mod;
      exp /= 2;
    }
    return result;
  }

  // Source:
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Simple_algebraic_field_extensions
  friend Polynomial inv_mod(const Polynomial &f, const Polynomial &mod) {
    if (f.field != mod.field)
      throw math_error() << "Cannot compute modular inverses for polynomials "
                            "over different fields";
    if (f.variable != mod.variable)
      throw math_error() << "Cannot compute modular inverses with polynomials "
                            "with two different variables '"
                         << f.variable << "' and '" << mod.variable << "'";

    const Polynomial zero = f.zero_poly(), one = f.one_poly();
    Polynomial t0 = zero, t1 = one, r0 = mod, r1 = f;

    while (r1 != zero) {
      const auto &[q, r2] = divmod(r0, r1);
      std::tie(r0, r1) = std::make_pair(r1, r2);
      std::tie(t0, t1) = std::make_pair(t1, t0 - q * t1);
    }

    if (r0.degree() > 0)
      throw math_error()
          << "Could not compute inverse: f and mod share a common factor "
          << r0;

    const Polynomial result = (r0.coeffs[0].inv() * t0) % mod;
    if ((result * f) % mod != 1)
      throw math_error() << "Assertion error: inv_mod did not compute inverse: "
                            "the inverse of f = "
                         << f << " mod " << mod << " was computed as " << t0
                         << " but (t0 * f) % mod = " << (result * f) % mod;
    return result;
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
      return std::make_pair(p.zero_poly(), p);
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
  Polynomial &operator/=(const Polynomial &other) {
    return *this = *this / other;
  }
  Polynomial &operator%=(const Polynomial &other) {
    return *this = *this % other;
  }

  friend constexpr bool operator==(const Polynomial &a, const Polynomial &b) {
    return a.variable == b.variable && a.support == b.support &&
           a.coeffs == b.coeffs;
  }
  friend constexpr bool operator!=(const Polynomial &a, const Polynomial &b) {
    return !(a == b);
  }
  friend bool operator==(const Polynomial &a, const integer_t elem) {
    return a.coeffs.size() == 1 && a.coeffs[0] == elem;
  }
  friend bool operator!=(const Polynomial &a, const integer_t elem) {
    return !(a == elem);
  }

  element_t at(const element_t &a) const {
    element_t result = zero;
    for (int d = coeffs.size() - 1; d >= 0; --d)
      result = result * a + coeffs[d];
    return result;
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
      const std::string coeff_str = p.field.value_to_string(coeff.value);
      const std::string display_str = coeff_str.find('+') == std::string::npos
                                          ? coeff_str
                                          : "(" + coeff_str + ")";

      if (d != 0) {
        if (coeff != p.one)
          os << display_str;
        os << p.variable;
      } else {
        os << display_str;
      }
      if (d > 1)
        os << "^" << d;
    }
    return os;
  }

  std::string to_debug_string() const {
    std::stringstream ss;
    ss << "{" << std::endl;
    ss << "  field: " << field << std::endl;
    ss << "  variable: '" << variable << "'" << std::endl;
    ss << "  coeffs: [";
    for (size_t i = 0; i < coeffs.size(); ++i) {
      ss << (i > 0 ? ", " : "");
      ss << coeffs[i];
    }
    ss << "  ]" << std::endl;
    ss << "  support: [";
    for (size_t i = 0; i < support.size(); ++i) {
      ss << (i > 0 ? ", " : "");
      ss << support[i];
    }
    ss << "  ]" << std::endl;
    ss << "}" << std::endl;
    return ss.str();
  }

  friend Polynomial polynomial_gcd(const Polynomial &f, const Polynomial &g) {
    if (f.field != g.field)
      throw math_error()
          << "Cannot compute gcd of polynomials over different fields";
    if (f.variable != g.variable)
      throw math_error()
          << "Cannot compute gcd of polynomials with two different variables '"
          << f.variable << "' and '" << g.variable << "'";

    Polynomial a = f.to_monic(), b = g.to_monic();
    while (true) {
      if (b == 0)
        return a;
      a = (a % b).to_monic();
      if (a == 0)
        return b;
      b = (b % a).to_monic();
    }
  }

  bool is_irreducible_rabin() const {
    if (*this == 0)
      return false;
    const Polynomial f = to_monic();
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
      const Polynomial h = (pow_mod(x, gmp::pow(q, ni), f) - x) % f;
      const Polynomial g = polynomial_gcd(f, h);
      if (g != 1)
        return false;
    }
    const Polynomial h = (pow_mod(x, gmp::pow(q, n), f) - x) % f;
    return h == 0;
  }

  // More efficient than computing the order and comparing to q - 1
  bool is_primitive_mod(const Polynomial &f) const {
    if (f % *this == 0)
      return false;
    const integer_t p = field.characteristic(), q = field.cardinality();
    const uint64_t k = field.degree(), n = f.degree();
    const integer_t full_order = gmp::pow(p, k * n) - 1;
    const Factorization factorization = factor_pk_minus_one(p, k * n);
    for (const Factor &factor : factorization) {
      if (pow_mod(*this, full_order / factor.base, f) == 1)
        return false;
    }
    return true;
  }
};

template <bool exact_degree, bool monic, class Field>
static inline Polynomial<Field> random_polynomial(const Field &field,
                                                  const std::string &variable,
                                                  const uint64_t degree) {
  std::vector<typename Field::element_t> coeffs(degree + 1,
                                                field.element(field.zero()));
  for (size_t i = 0; i < degree; ++i)
    coeffs[i] = field.random_element();

  // monic is stronger than exact_degree, so handle that first
  if constexpr (monic) {
    coeffs[degree] = field.element(field.one());
  } else if constexpr (exact_degree) {
    while (coeffs[degree] == field.element(field.zero()))
      coeffs[degree] = field.random_element();
  } else {
    coeffs[degree] = field.random_element();
  }
  return Polynomial(field, variable, coeffs);
}
