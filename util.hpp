
#pragma once

#include <gmpxx.h>
#include <sstream>
#include <stdexcept>

using integer_t = mpz_class;
using rational_t = mpq_class;

inline bool is_prime(const integer_t n) {
  return mpz_probab_prime_p(n.get_mpz_t(), 50);
}

inline int64_t to_int(const integer_t n) { return mpz_get_si(n.get_mpz_t()); }
inline uint64_t to_uint(const integer_t n) { return mpz_get_ui(n.get_mpz_t()); }

class math_error : public std::runtime_error {
  std::stringstream ss;
  mutable std::string str;

public:
  math_error(const std::string &val = "") : std::runtime_error(""), ss(val) {}
  math_error(const math_error &other)
      : std::runtime_error(""), ss(other.ss.str()) {}
  template <typename T> math_error &operator<<(const T &val) {
    ss << val;
    return *this;
  }
  const char *what() const throw() override {
    str = ss.str();
    return str.data();
  }
};

template <typename T> struct Indexed {
  static inline size_t next_id = 0;
  size_t m_id;
  Indexed() { m_id = next_id++; }

  bool operator==(const Indexed<T> &other) const { return m_id == other.m_id; }
};
