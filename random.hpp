
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
