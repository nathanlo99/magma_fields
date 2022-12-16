
#pragma once

#include "field.hpp"

#include <iostream>
#include <vector>

template <class Field> struct Vector {
  using element_t = typename Field::element_t;
  const Field &field;
  size_t size;
  std::vector<element_t> data;

  Vector(const Field &field, const size_t size)
      : field(field), size(size), data(size, field.element(field.zero())) {}

  Vector(const Field &field, const size_t size,
         const std::vector<element_t> &data)
      : field(field), size(size), data(data) {
    assert(data.size() == size);
  }

  Vector(const Field &field, const size_t size,
         const std::vector<integer_t> &input_data)
      : field(field), size(size) {
    data.reserve(size);
    for (size_t i = 0; i < size; ++i)
      data.push_back(field(input_data[i]));
  }

  element_t &operator[](const size_t i) { return data[i]; }
  const element_t &operator[](const size_t i) const { return data[i]; }

  bool operator==(const Vector &other) const { return data == other.data; }

  friend std::ostream &operator<<(std::ostream &os, const Vector &vec) {
    os << "[";
    for (size_t i = 0; i < vec.size; ++i) {
      os << (i == 0 ? " " : ", ") << vec[i];
    }
    os << " ]";
    return os;
  }
};
