
#pragma once

#include <chrono>
#include <iostream>

inline long long get_time_ns() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}

inline long long get_time_ms() {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}

inline long long get_elapsed_ms() {
  static long long start_time = get_time_ms();
  return get_time_ms() - start_time;
}

template <typename Func>
inline void timeit(const std::string_view &msg, const Func &f) {
  std::cout << "[" << msg << "] | Starting" << std::endl;
  const auto start_ms = get_time_ns();
  f();
  const auto end_ms = get_time_ns();
  std::cout << "[" << msg << "] | Done in " << (end_ms - start_ms) / 1e9
            << " seconds" << std::endl
            << std::endl;
}

template <typename Func>
inline void notimeit(const std::string_view &msg, const Func &f) {}
