

#pragma once

#include "field.hpp"
#include "util.hpp"
#include <cstdint>
#include <vector>

// MediumPrimeField: FiniteField(p) for primes 2^16 <= p < 2^32
struct MediumPrimeField : FiniteField<uint32_t> {
  using value_t = uint32_t;
  using element_t = FieldElement<MediumPrimeField>;

  const value_t p;

  MediumPrimeField(const integer_t p) : p(to_int(p)) {
    if (p <= 0 || !is_prime(p))
      throw math_error() << "MediumPrimeField expects a positive prime, got "
                         << p;
    if (!p.fits_uint_p())
      throw math_error() << "Cannot instantiate MediumPrimeField with modulus "
                         << p
                         << ", expects a number which fits in an unsigned int";
  }

  integer_t characteristic() const override { return p; }
  integer_t degree() const override { return 1; }
  integer_t cardinality() const override { return p; }

  value_t zero() const override { return 0; }
  value_t one() const override { return 1; }
  value_t integer(const integer_t number) const override {
    const int64_t num = to_int(number);
    return ((num % p) + p) % p;
  }
  element_t operator()(const integer_t num) const {
    return element_t(*this, integer(num));
  }

  value_t neg(const value_t a) const override { return a == 0 ? a : p - a; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % p;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a + p - b) % p;
  }

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
    return (static_cast<uint64_t>(a) * b) % p;
  }

  // NOTE: We take the default implementations of div and pow

  bool eq(const value_t a, const value_t b) const override { return a == b; }
};
