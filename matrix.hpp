
#pragma once

#include "error.hpp"
#include "field.hpp"

template <typename Field> struct Matrix {
  using element_t = typename Field::element_t;

  const Field &field;
  const size_t rows, cols;
  std::vector<std::vector<element_t>> data;

  Matrix(const Field &field, const size_t rows, const size_t cols)
      : field(field), rows(rows), cols(cols) {
    data = std::vector<std::vector<element_t>>(
        rows, std::vector<element_t>(cols, field.element(field.zero())));
    check_invariants();
  }

  Matrix(const Field &field, const size_t rows, const size_t cols,
         const std::vector<std::vector<element_t>> &data)
      : field(field), rows(rows), cols(cols), data(data) {
    check_invariants();
  }

  bool is_square() const { return rows == cols; }
  void check_invariants() const {
    assert(data.size() == rows);
    for (const auto &row : data)
      assert(row.size() == cols);
  }

  std::vector<element_t> &operator[](const size_t i) { return data[i]; }
  const std::vector<element_t> &operator[](const size_t i) const {
    return data[i];
  }

  // Addition
  Matrix &operator+=(const Matrix &other) {
    if (field != other.field)
      throw math_error("Cannot add matrices with different base fields");
    if (rows != other.rows)
      throw math_error("Cannot add matrices with different numbers of rows");
    if (cols != other.cols)
      throw math_error("Cannot add matrices with different numbers of cols");

    for (size_t row = 0; row < rows; ++row)
      for (size_t col = 0; col < cols; ++col)
        data[row][col] += other[row][col];
    return *this;
  }
  friend Matrix operator+(const Matrix &a, const Matrix &b) {
    Matrix result = a;
    result += b;
    return result;
  }

  // Subtraction
  Matrix &operator-=(const Matrix &other) {
    if (field != other.field)
      throw math_error("Cannot subtract matrices with different base fields");
    if (rows != other.rows)
      throw math_error(
          "Cannot subtract matrices with different numbers of rows");
    if (cols != other.cols)
      throw math_error(
          "Cannot subtract matrices with different numbers of cols");

    for (size_t row = 0; row < rows; ++row)
      for (size_t col = 0; col < cols; ++col)
        data[row][col] -= other[row][col];
    return *this;
  }
  friend Matrix operator-(const Matrix &a, const Matrix &b) {
    Matrix result = a;
    result -= b;
    return result;
  }

  // Matrix multiplication
  friend Matrix operator*(const Matrix &a, const Matrix &b) {
    if (a.field != b.field)
      throw math_error("Cannot multiply matrices with different base fields");
    if (a.cols != b.rows)
      throw math_error("Cannot multiply matrices with incompatible dimensions");

    Matrix result(a.field, a.rows, b.cols);
    for (size_t i = 0; i < result.rows; ++i)
      for (size_t j = 0; j < result.cols; ++j)
        for (size_t k = 0; k < a.cols; ++k)
          result[i][j] += a[i][k] * b[k][j];
    return result;
  }

  // Element multiplication
  Matrix &operator*=(const element_t &value) {
    for (size_t row = 0; row < rows; ++row)
      for (size_t col = 0; col < cols; ++col)
        data[row][col] *= value;
    return *this;
  }
  friend Matrix operator*(const Matrix &m, const element_t &k) {
    Matrix result = m;
    result *= k;
    return result;
  }
  friend Matrix operator*(const element_t &k, const Matrix &m) { return m * k; }

  // Inversion (for square matrices only)

  // Printing
  friend std::ostream &operator<<(std::ostream &os, const Matrix &m) {
    const size_t width = m.field.characteristic().get_str().size();
    for (size_t i = 0; i < m.rows; ++i) {
      os << (i > 0 ? "  " : "[ ");
      for (size_t j = 0; j < m.cols; ++j) {
        os << (j > 0 ? ", " : "") << std::setw(width) << m.data[i][j];
      }
      os << (i + 1 == m.rows ? " ]" : "\n");
    }
    return os;
  }
};
