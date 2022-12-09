
#pragma once

#include "field.hpp"
#include "util.hpp"
#include <cstdint>
#include <vector>

// LargePrimeField: FiniteField(p) for primes 2^32 <= p
struct LargePrimeField : Field<integer_t> {
  using value_t = integer_t;
  using element_t = FieldElement<LargePrimeField>;

  const value_t p;

  LargePrimeField(const integer_t p) : p(p) {
    if (p <= 0 || !is_prime(p))
      throw math_error() << "LargePrimeField expects a positive prime, got "
                         << p;
    if (p.fits_uint_p())
      throw math_error() << "Cannot instantiate LargePrimeField with modulus "
                         << p << ", expects a number larger than 2^32";
  }

  integer_t characteristic() const override { return p; }
  integer_t degree() const override { return 1; }
  integer_t cardinality() const override { return p; }

  value_t zero() const override { return 0; }
  value_t one() const override { return 1; }
  value_t integer(const integer_t number) const override {
    return ((number % p) + p) % p;
  }
  uint32_t as_integer(const value_t number) const { return to_uint(number); }
  element_t operator()(const integer_t num) const {
    return element_t(*this, integer(num));
  }
  element_t element(const value_t value) const {
    return element_t(*this, value);
  }

  value_t neg(const value_t a) const override { return a == 0 ? a : p - a; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % p;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a + p - b) % p;
  }

  // TODO: Use Montgomery modular representation to make this faster
  value_t inv(const value_t a) const override {
    value_t r0 = p, r1 = a, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
    while (r1 != 0) {
      const value_t q = r0 / r1, r2 = r0 % r1;
      std::tie(r0, s0, t0, r1, s1, t1) =
          std::make_tuple(r1, s1, t1, r2, s0 - s1 * q, t0 - t1 * q);
    }
    return t0 >= 0 ? t0 : t0 + p;
  }
  value_t mul(const value_t a, const value_t b) const override {
    return (a * b) % p;
  }

  // NOTE: We take the default implementations of div and pow

  bool eq(const value_t a, const value_t b) const override { return a == b; }

  std::string to_string(const value_t a) const override { return a.get_str(); }

  friend std::ostream &operator<<(std::ostream &os,
                                  const LargePrimeField &field) {
    return os << "Finite field of order " << field.p;
  }
};
