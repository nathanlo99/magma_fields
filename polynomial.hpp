
#pragma once

#include "gmp.hpp"

#include <string>
#include <tuple>
#include <vector>

template <class Field = integer_t> struct Polynomial {
  using element_t = typename Field::element_t;

  const Field &field;
  char variable;
  const element_t zero, one;
  std::vector<element_t> coeffs;

  Polynomial(const Field &field, const char variable,
             const std::vector<element_t> &coeffs = {})
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs) {
    if (coeffs.empty()) {
      this->coeffs = {zero, one};
      return;
    }

    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();
  }

public:
  Polynomial(const Field &field, const char variable,
             const std::vector<integer_t> &coeffs)
      : field(field), variable(variable), zero(field.element(field.zero())),
        one(field.element(field.one())), coeffs(coeffs.size(), zero) {
    if (coeffs.empty()) {
      this->coeffs = {zero, one};
      return;
    }

    for (size_t i = 0; i < coeffs.size(); ++i)
      this->coeffs[i] = field(coeffs[i]);

    // Remove leading zeroes
    while (this->coeffs.size() > 1 && this->coeffs.back() == zero)
      this->coeffs.pop_back();
  }

  Polynomial(const Polynomial &other)
      : field(other.field), variable(other.variable),
        zero(field.element(field.zero())), one(field.element(field.one())),
        coeffs(other.coeffs) {}

  Polynomial &operator=(const Polynomial &other) {
    if (field != other.field)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different field");
    if (variable != other.variable)
      throw math_error(
          "Cannot assign Polynomial to Polynomial with different variable");
    coeffs = other.coeffs;
    return *this;
  }

  inline size_t degree() const { return coeffs.size() - 1; }
  inline element_t operator[](const size_t idx) const {
    return idx < coeffs.size() ? coeffs[idx] : zero;
  }

  Polynomial operator-() const {
    const size_t degree_plus_one = coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, zero);
    // TOOD: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i) {
      result_coeffs[i] = -coeffs[i];
    }
    return Polynomial(field, variable, result_coeffs);
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

  friend Polynomial operator*(const element_t k, const Polynomial &p) {
    const size_t degree_plus_one = p.coeffs.size();
    std::vector<element_t> result_coeffs(degree_plus_one, p.zero);
    // TODO: Parallelize this
    for (size_t i = 0; i < degree_plus_one; ++i)
      result_coeffs[i] = k * p.coeffs[i];
    return Polynomial(p.field, p.variable, result_coeffs);
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
    for (size_t i = 0; i <= a_degree; ++i) {
      for (size_t j = 0; j <= b_degree; ++j) {
        result_coeffs[i + j] += a[i] * b[j];
      }
    }
    return Polynomial(a.field, a.variable, result_coeffs);
  }
  Polynomial &operator*=(const Polynomial &other) {
    return *this = *this * other;
  }

  Polynomial operator^(uint64_t exp) const {
    Polynomial result(field, variable, {one}), pow = *this;
    while (exp > 0) {
      if (exp % 2 == 1)
        result *= pow;
      pow *= pow;
      exp /= 2;
    }
    return result;
  }

  friend Polynomial pow_mod(const Polynomial &f, uint64_t exp,
                            const Polynomial &mod) {
    Polynomial result = Polynomial(f.field, 'x', {f.one}), base = f;
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
    const Polynomial zero = Polynomial(f.field, f.variable, {f.zero}),
                     one = Polynomial(f.field, f.variable, {f.one});
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
      const element_t pivot = p_coeffs[q_degree + d] * inv_q_leading_coeff;
      result_coeffs[d] = pivot;
      p_coeffs[q_degree + d] = p.zero;
      for (size_t i = 0; i < q_degree; ++i)
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
    return a.variable == b.variable && a.coeffs == b.coeffs;
  }

  std::string to_string() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }

  friend std::ostream &operator<<(std::ostream &os, const Polynomial &p) {
    if (p.coeffs.size() == 1 && p[0] == p.zero)
      return os << "0";
    bool first = true;
    for (int d = static_cast<int>(p.coeffs.size()) - 1; d >= 0; --d) {
      if (p[d] == p.zero)
        continue;
      if (first) {
        first = false;
      } else {
        os << " + ";
      }

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
};
