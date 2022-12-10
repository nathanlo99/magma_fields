
#pragma once

#include <gmpxx.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "error.hpp"
#include "gmp.h"

using integer_t = mpz_class;
using rational_t = mpq_class;

inline bool is_prime(const integer_t n) {
  return mpz_probab_prime_p(n.get_mpz_t(), 50);
}
inline void next_prime(integer_t &n) {
  mpz_nextprime(n.get_mpz_t(), n.get_mpz_t());
}

inline int64_t to_int(const integer_t n) {
  if (!n.fits_sint_p())
    throw math_error() << "Could not convert " << n << " to int: too large";
  return mpz_get_si(n.get_mpz_t());
}

inline uint64_t to_uint(const integer_t n) {
  if (!n.fits_uint_p())
    throw math_error() << "Could not convert " << n << " to uint: too large";
  return mpz_get_ui(n.get_mpz_t());
}

inline integer_t from_uint(const uint64_t n) {
  integer_t result;
  mpz_set_ui(result.get_mpz_t(), n);
  return result;
}

inline integer_t unsigned_mod(const integer_t n, const integer_t p) {
  return ((n % p) + p) % p;
}