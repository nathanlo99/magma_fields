
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include <cstdint>
#include <tuple>
#include <vector>

// SmallPrimeField: FiniteField(p) for primes 256 <= p < 65536
struct SmallPrimeField : Field<uint32_t> {
  using value_t = uint32_t;
  using element_t = FieldElement<SmallPrimeField>;

  const value_t p;

  SmallPrimeField(const integer_t p) : p(gmp::to_int(p)) {
    if (p <= 0 || !gmp::is_prime(p))
      throw math_error() << "SmallPrimeField expects a positive prime, got "
                         << p;
    if (!p.fits_ushort_p())
      throw math_error()
          << "Cannot instantiate SmallPrimeField with modulus " << p
          << ", expects a number which fits in an unsigned short";
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return 1; }
  integer_t cardinality() const override { return p; }
  FieldType type() const override { return SmallPrimeFieldType; }

  value_t zero() const override { return 0; }
  value_t one() const override { return 1; }
  value_t integer(const integer_t number) const override {
    return gmp::to_uint(gmp::unsigned_mod(number, p));
  }
  uint32_t as_integer(const value_t number) const { return number; }
  element_t operator()(const integer_t num) const {
    return element_t(*this, integer(num));
  }
  element_t element(const value_t value) const {
    return element_t(*this, value);
  }
  element_t primitive_element() const {
    for (value_t c = 1; c < p; ++c) {
      if (is_primitive(c))
        return element_t(*this, c);
    }
    __builtin_unreachable();
  }

  value_t neg(const value_t a) const override { return a == 0 ? a : p - a; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % p;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a + p - b) % p;
  }

  value_t inv(const value_t a) const override {
    int32_t r0 = p, r1 = a, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
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

  std::string to_string(const value_t a) const override {
    return std::to_string(a);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const SmallPrimeField &field) {
    return os << "Finite field of order " << field.p;
  }
};
