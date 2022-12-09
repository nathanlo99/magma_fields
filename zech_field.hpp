
#pragma once

#include "field.hpp"
#include "gmp.h"
#include "gmp.hpp"
#include "polynomial.hpp"

#include <iostream>

template <class BaseField>
uint32_t get_polynomial_index(const Polynomial<BaseField> &f) {
  // Assumes the BaseField has degree 1
  const auto base_field = f.field;
  const uint32_t p = to_uint(base_field.characteristic());
  uint32_t result = 0;
  for (int i = f.degree(); i >= 0; --i) {
    result = result * p + base_field.as_integer(f.coeffs[i].value);
  }
  return result;
}

template <class BaseField>
std::vector<uint32_t> compute_zech_table(const uint32_t p, const uint32_t k,
                                         const uint32_t q,
                                         const Polynomial<BaseField> &f) {
  // Compute the Zech logarithm table, find s[r] such that x^s[r] = x^r + 1
  std::vector<uint32_t> result(q, -1);
  const auto base_field = f.field;
  const auto one = Polynomial(base_field, 'x', {1}),
             x = Polynomial(base_field, 'x', {0, 1});

  Polynomial<BaseField> pow = one;
  std::vector<uint32_t> log_x(q, q);
  std::vector<uint32_t> next_index(q, q);
  for (uint32_t exp = 0; exp + 1 < q; ++exp) {
    if (exp > 0 && pow == one)
      throw math_error() << "Given polynomial was not primitive, x^" << exp
                         << " = 1";
    const uint32_t idx = get_polynomial_index(pow);
    next_index[exp] = get_polynomial_index(pow + one);
    log_x[idx] = exp;
    pow = (pow * x) % f;
  }

  // Special cases
  result[q - 1] = 0;
  const uint32_t log_negative_one = p == 2 ? 0 : (q - 1) / 2;
  result[log_negative_one] = q - 1;

  for (uint32_t i = 0; i < q; ++i) {
    if (result[i] >= q)
      result[i] = log_x[next_index[i]];
  }

  return result;
}

// For prime extension fields with cardinality at most 2^20
//
// Elements are represented as powers of a primitive element, and a table of
// Zech logarithms is stored for use in addition.
template <class BaseField> struct ZechField : Field<uint32_t> {
  using value_t = uint32_t;
  using element_t = FieldElement<ZechField>;

  uint32_t p, k, q;
  const BaseField &base_field;
  std::vector<value_t> zech_table;
  const Polynomial<BaseField> &f;

  ZechField(const BaseField &base_field, const Polynomial<BaseField> &f)
      : p(to_uint(base_field.characteristic())), k(f.degree()),
        base_field(base_field), f(f) {
    if (base_field != f.field)
      throw math_error()
          << "Polynomial base field did not match provided base_field";
    if (base_field.degree() != 1)
      throw math_error()
          << "Zech field expected base field with prime order, got " << f.field;

    integer_t p_tmp = p, q_tmp, max_size = 1 << 20;
    mpz_pow_ui(q_tmp.get_mpz_t(), p_tmp.get_mpz_t(), k); // q = p^k
    if (q_tmp > max_size)
      throw math_error() << "Zech field expected cardinality at most 2^20, got "
                         << p << "^" << k << " = " << q;
    q = to_uint(q_tmp);

    zech_table = compute_zech_table(this->p, k, q, f);
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return k; }
  integer_t cardinality() const override { return q; }

  value_t zero() const override { return q - 1; }
  value_t one() const override { return 0; }
  value_t integer(const integer_t number) const override {
    value_t result = zero(), base = one();
    // This conversion will fit since p fits in a uint32_t
    uint32_t num = to_uint(unsigned_mod(number, p));
    while (num != 0) {
      if (num % 2 == 1)
        result = add(result, base);
      base = add(base, base);
      num /= 2;
    }
    return result;
  }
  element_t operator()(const integer_t num) const {
    return element_t(*this, integer(num));
  }
  element_t element(const value_t value) const {
    return element_t(*this, value);
  }
  element_t generator() const { return element(1); } // g^1

  value_t neg(const value_t a) const override {
    if (p == 2 || a == zero())
      return a;
    const value_t negative_one = (q - 1) / 2;
    return (a + negative_one) % (q - 1);
  }
  value_t add(const value_t a, const value_t b) const override {
    if (a == zero())
      return b;
    if (b == zero())
      return a;
    if (a < b) {
      // g^a + g^b = g^a (1 + g^{b-a}) = g^a ( g^s[b-a] )
      return mul(a, zech_table[b - a]);
    } else {
      // Similarly
      return mul(b, zech_table[a - b]);
    }
  }
  value_t sub(const value_t a, const value_t b) const override {
    return add(a, neg(b));
  }

  value_t inv(const value_t a) const override {
    if (a == zero())
      throw math_error("Division by zero");
    return a == 0 ? 0 : q - 1 - a;
  }
  value_t mul(const value_t a, const value_t b) const override {
    if (a == zero() || b == zero())
      return zero();
    return (a + b) % (q - 1);
  }
  value_t div(const value_t a, const value_t b) const override {
    if (b == zero())
      throw math_error("Division by zero");
    if (a == zero())
      return zero();
    return (a + (q - 1) - b) % (q - 1);
  }
  value_t pow(const value_t a, const integer_t exp) const override {
    if (a == zero())
      return zero();
    const integer_t modded_exp = unsigned_mod(exp, q - 1);
    return (a * to_uint(modded_exp)) % (q - 1);
  }

  virtual bool eq(const value_t a, const value_t b) const override {
    return a == b;
  }

  std::string to_string(const value_t a) const override {
    if (a == q - 1)
      return "0";
    const Polynomial<BaseField> x = Polynomial(base_field, 'x', {0, 1});
    const Polynomial<BaseField> actual = pow_mod(x, a, f);
    return actual.to_string();
  }

  friend std::ostream &operator<<(std::ostream &os, const ZechField &field) {
    return os << "ZechField over [" << field.base_field
              << "] with generating polynomial " << field.f;
  }
};
