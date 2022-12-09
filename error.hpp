
#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

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
