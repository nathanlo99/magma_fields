
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "logger.hpp"
#include "random.hpp"
#include "vector.hpp"

#include <cstdint>
#include <tuple>
#include <vector>

// SmallPrimeField: FiniteField(p) for primes p < 65536
struct SmallPrimeField : Field<uint32_t, SmallPrimeField> {
  using value_t = uint32_t;
  using element_t = FieldElement<SmallPrimeField>;
  using prime_field_t = SmallPrimeField;
  using ground_field_t = SmallPrimeField;
  using vector_t = Vector<prime_field_t>;

  const value_t p;
  std::vector<value_t> inverses;

  SmallPrimeField(const integer_t p)
      : p(gmp::to_uint(p)), inverses(this->p, 0) {
    if (p <= 0 || !gmp::is_prime(p))
      throw math_error() << "SmallPrimeField expects a positive prime, got "
                         << p;
    if (!p.fits_ushort_p())
      throw math_error()
          << "Cannot instantiate SmallPrimeField with modulus " << p
          << ", expects a number which fits in an unsigned short";

    const auto compute_inverse_mod = [](const value_t a, const value_t p) {
      int32_t r0 = p, r1 = a, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
      while (r1 != 0) {
        const value_t q = r0 / r1, r2 = r0 % r1;
        std::tie(r0, s0, t0, r1, s1, t1) =
            std::make_tuple(r1, s1, t1, r2, s0 - s1 * q, t0 - t1 * q);
      }
      return t0 >= 0 ? t0 : t0 + p;
    };
    for (value_t val = 1; val < p; ++val) {
      inverses[val] = compute_inverse_mod(val, this->p);
    }
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return 1; }
  integer_t cardinality() const override { return p; }
  FieldType type() const override { return FieldType::SmallPrime; }
  const prime_field_t &prime_field() const override { return *this; }
  const ground_field_t &ground_field() const { return *this; }

  value_t zero() const override { return 0; }
  value_t one() const override { return 1; }
  value_t integer(const integer_t number) const override {
    return gmp::to_uint(gmp::unsigned_mod(number, p));
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
    for (value_t c = 1; c < p; ++c) {
      if (is_primitive(c)) {
        result = c;
        return element_t(*this, result);
      }
    }
    __builtin_unreachable();
  }
  // NOTE: This is not cryptographically secure or even uniform.
  element_t random_element() const {
    return element_t(*this, random_uint64() % p);
  }

  // Vector space structure
  element_t generating_element() const { return element_t(*this, 1); }
  vector_t to_vector(const element_t elem) const {
    return vector_t(*this, 1, {elem});
  }
  element_t from_vector(const vector_t vec) const { return vec[0]; }

  // Field operations
  value_t neg(const value_t a) const override { return a == 0 ? a : p - a; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % p;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a + p - b) % p;
  }

  value_t inv(const value_t a) const override {
    if (a == 0)
      throw math_error("Division by zero");
    return inverses[a];
  }
  value_t mul(const value_t a, const value_t b) const override {
    return (a * b) % p;
  }

  // NOTE: We take the default implementations of div and pow

  bool eq(const value_t a, const value_t b) const override { return a == b; }

  std::string value_to_string(const value_t a) const override {
    return std::to_string(a);
  }

  std::string to_string() const override {
    std::stringstream ss;
    ss << "Prime field with order " << p;
    return ss.str();
  }
};
