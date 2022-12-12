
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "logger.hpp"
#include "random.hpp"

#include <cstdint>
#include <tuple>
#include <vector>

// LargePrimeField: FiniteField(p) for primes 2^32 <= p
struct LargePrimeField : Field<integer_t> {
  using value_t = integer_t;
  using element_t = FieldElement<LargePrimeField>;

  const value_t p;

  LargePrimeField(const integer_t p) : p(p) {
    if (p <= 0 || !gmp::is_prime(p))
      throw math_error() << "LargePrimeField expects a positive prime, got "
                         << p;
    if (p.fits_uint_p())
      throw math_error() << "Cannot instantiate LargePrimeField with modulus "
                         << p << ", expects a number larger than 2^32";
    log() << "Done constructing LargePrimeField with cardinality " << p
          << std::endl;
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return 1; }
  integer_t cardinality() const override { return p; }
  FieldType type() const override { return FieldType::LargePrime; }

  value_t zero() const override { return 0; }
  value_t one() const override { return 1; }
  value_t integer(const integer_t number) const override {
    return gmp::unsigned_mod(number, p);
  }
  uint32_t as_integer(const value_t number) const {
    return gmp::to_uint(number);
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
  // NOTE: This is not cryptographically secure
  element_t random_element() const {
    integer_t random_number;
    mpz_urandomm(random_number.get_mpz_t(), random_state, p.get_mpz_t());
    return element_t(*this, random_number);
  }

  value_t neg(const value_t a) const override { return a == 0 ? a : p - a; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % p;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a + p - b) % p;
  }

  // TODO: Use Montgomery modular representation to make this even faster
  value_t inv(const value_t a) const override {
    value_t result;
    mpz_invert(result.get_mpz_t(), a.get_mpz_t(), p.get_mpz_t());
    return result;
  }
  value_t mul(const value_t a, const value_t b) const override {
    return (a * b) % p;
  }

  // NOTE: We take the default implementations of div and pow

  bool eq(const value_t a, const value_t b) const override { return a == b; }

  friend std::ostream &operator<<(std::ostream &os,
                                  const LargePrimeField &field) {
    return os << "Prime field with order " << field.p;
  }

  std::string value_to_string(const value_t value) const override {
    return value.get_str();
  }

  std::string to_string() const override {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }
};
