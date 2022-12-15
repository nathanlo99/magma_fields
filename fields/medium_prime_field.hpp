
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "random.hpp"

#include <cstdint>
#include <tuple>
#include <vector>

// MediumPrimeField: FiniteField(p) for primes 2^16 <= p < 2^32
struct MediumPrimeField : Field<uint64_t> {
  using value_t = uint64_t;
  using element_t = FieldElement<MediumPrimeField>;

  const value_t p;

  MediumPrimeField(const integer_t p) : p(gmp::to_uint(p)) {
    if (p <= 0 || !gmp::is_prime(p))
      throw math_error() << "MediumPrimeField expects a positive prime, got "
                         << p;
    if (!p.fits_uint_p())
      throw math_error() << "Cannot instantiate MediumPrimeField with modulus "
                         << p
                         << ", expects a number which fits in an unsigned int";
  }

  integer_t characteristic() const override { return gmp::from_uint(p); }
  uint32_t degree() const override { return 1; }
  integer_t cardinality() const override { return gmp::from_uint(p); }
  FieldType type() const override { return FieldType::MediumPrime; }

  value_t zero() const override { return 0; }
  value_t one() const override { return 1; }
  value_t integer(const integer_t number) const override {
    return gmp::to_uint(gmp::unsigned_mod(number, gmp::from_uint(p)));
  }
  integer_t as_integer(const value_t number) const {
    return gmp::from_uint(number);
  }
  element_t operator()(const integer_t num) const {
    return element_t(*this, integer(num));
  }
  element_t element(const value_t value) const {
    return element_t(*this, value);
  }
  element_t primitive_element() const {
    static value_t result = 0;
    if (result != 0)
      return element_t(*this, result);
    for (value_t c = 2; c < p; ++c) {
      if (is_primitive(c)) {
        result = c;
        return element_t(*this, result);
      }
    }
    __builtin_unreachable();
  }
  element_t generating_element() const { return element_t(*this, 1); }
  // NOTE: This is not cryptographically secure or even uniform.
  element_t random_element() const {
    return element_t(*this, random_uint64() % p);
  }

  value_t neg(const value_t a) const override { return a == 0 ? a : p - a; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % p;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a + p - b) % p;
  }

  value_t inv(const value_t a) const override {
    int64_t r0 = p, r1 = a, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
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

  std::string value_to_string(const value_t a) const override {
    return std::to_string(a);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const MediumPrimeField &field) {
    return os << "Prime field with order " << field.p;
  }

  std::string to_string() const override {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }
};
