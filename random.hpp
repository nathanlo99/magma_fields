
#pragma once

#include "gmp.hpp"

#include <cstdint>
#include <random>

extern gmp_randstate_t random_state;
inline void init_gmp_random_state() { gmp_randinit_default(random_state); }

inline uint64_t random_uint64() {
  static std::mt19937 gen(127);
  static std::uniform_int_distribution<uint64_t> dist;
  return dist(gen);
}

// Return a random number between 0 and n-1 inclusive
inline integer_t random_integer_t(const integer_t n) {
  integer_t random_number;
  mpz_urandomm(random_number.get_mpz_t(), random_state, n.get_mpz_t());
  return random_number;
}
