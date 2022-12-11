
#pragma once

#include <cstdint>
#include <random>

inline uint64_t random_uint64() {
  static std::mt19937 gen(127);
  static std::uniform_int_distribution<uint64_t> dist;
  return dist(gen);
}
