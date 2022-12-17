
#pragma once

#include "error.hpp"
#include "field.hpp"
#include "vector.hpp"

template <typename Field> struct Matrix {
  using element_t = typename Field::element_t;

  const Field &field;
  size_t rows, cols;
  std::vector<std::vector<element_t>> data;

  Matrix(const Field &field) : field(field), rows(0), cols(0), data() {}
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

  Matrix(const Field &field, const size_t rows, const size_t cols,
         const std::vector<std::vector<integer_t>> &input_data)
      : field(field), rows(rows), cols(cols) {
    for (size_t row = 0; row < rows; ++row) {
      data.emplace_back();
      for (size_t col = 0; col < cols; ++col) {
        data[row].push_back(field(input_data[row][col]));
      }
    }
    check_invariants();
  }

  Matrix(const Matrix &other) = default;
  Matrix(Matrix &&other) = default;

  Matrix &operator=(const Matrix &other) {
    if (field != other.field)
      throw math_error(
          "Cannot assign matrix to matrix with different underlying field");
    rows = other.rows;
    cols = other.cols;
    data = other.data;
    return *this;
  }

  Matrix &operator=(Matrix &&other) {
    if (field != other.field)
      throw math_error(
          "Cannot assign matrix to matrix with different underlying field");
    rows = other.rows;
    cols = other.cols;
    data = std::move(other.data);
    return *this;
  }

  Matrix eye(const size_t size) const {
    const element_t one = field.element(field.one());
    Matrix result(field, size, size);
    for (size_t i = 0; i < size; ++i)
      result[i][i] = one;
    return result;
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

  Matrix transpose() const {
    Matrix result(field, cols, rows);
    for (size_t row = 0; row < rows; ++row)
      for (size_t col = 0; col < cols; ++col)
        result.data[col][row] = data[row][col];
    result.check_invariants();
    return result;
  }

  void add_row(const std::vector<element_t> &row) {
    if (row.size() != cols)
      throw math_error()
          << "Cannot add row of incompatible size, width of matrix was " << cols
          << " but got row with " << row.size() << " elements";
    data.push_back(row);
    ++rows;
    check_invariants();
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

  friend Vector<Field> operator*(const Matrix &a, const Vector<Field> &v) {
    if (a.cols != v.size)
      throw math_error()
          << "Cannot multiply matrix with vector: incompatible sizes ("
          << a.rows << " x " << a.cols << ") and (" << v.size << " x " << 1
          << ")";
    Vector<Field> result(a.field, a.rows);
    for (size_t row = 0; row < a.rows; ++row)
      for (size_t col = 0; col < a.cols; ++col)
        result[row] += a[row][col] * v[col];
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

  // Row-reduces the matrix and returns the rank
  size_t row_reduce() {
    const element_t zero = field.element(field.zero());

    // Finds the pivot row and column and swaps it to the top left, and return
    // the pivot column
    const auto find_pivot = [&](const size_t row, const size_t column) {
      for (size_t col = column; col < cols; ++col) {
        size_t pivot_row = rows;
        for (size_t lead_row = row; lead_row < rows; ++lead_row) {
          if (data[lead_row][col] != zero) {
            pivot_row = lead_row;
            break;
          }
        }
        if (pivot_row != rows) {
          std::swap(data[row], data[pivot_row]);
          return col;
        }
      }
      // Otherwise, every column past [column] is all zeroes: no more pivots
      return cols;
    };

    const auto divide_row = [&](const size_t row, const element_t elem) {
      assert(elem != zero);
      const element_t inv = elem.inv();
      for (size_t col = 0; col < cols; ++col)
        data[row][col] *= inv;
    };

    const auto subtract_multiple_of_row = [&](const size_t other_row,
                                              const element_t multiple,
                                              const size_t row) {
      if (multiple == zero)
        return;
      for (size_t col = 0; col < cols; ++col)
        data[other_row][col] -= multiple * data[row][col];
    };

    size_t col = 0;
    for (size_t row = 0; row < rows; ++row) {
      const size_t pivot_column = find_pivot(row, col);
      if (pivot_column == cols)
        return row;
      assert(data[row][pivot_column] != zero);
      divide_row(row, data[row][pivot_column]);
      for (size_t other_row = 0; other_row < rows; ++other_row) {
        if (other_row != row)
          subtract_multiple_of_row(other_row, data[other_row][pivot_column],
                                   row);
      }
    }
    return rows; // Full-rank!
  }

  size_t rank() const {
    Matrix tmp = *this;
    return tmp.row_reduce();
  }

  // Given b, solves Ax = b and returns x
  Vector<Field> solve(const Vector<Field> &b) const {
    if (b.size != rows)
      throw math_error() << "Incompatible matrix sizes: matrix had " << rows
                         << " rows while supplied column vector had " << b.size
                         << " elements";
    const element_t zero = field.element(field.zero()),
                    one = field.element(field.one());
    std::vector<std::vector<element_t>> augmented_coeffs(
        rows, std::vector<element_t>(cols + 1, zero));
    for (size_t row = 0; row < rows; ++row) {
      for (size_t col = 0; col < cols; ++col)
        augmented_coeffs[row][col] = data[row][col];
      augmented_coeffs[row][cols] = b[row];
    }
    Matrix augmented(field, rows, cols + 1, augmented_coeffs);
    const size_t rank = augmented.row_reduce();
    if (rank != rows || augmented[rows - 1][rows - 1] != one)
      throw math_error("Could not solve equation, matrix was singular");

    Vector result(field, b.size);
    for (size_t row = 0; row < rows; ++row)
      result[row] = augmented[row][cols];
    assert((*this) * result == b);
    return result;
  }

  Matrix inverse() const {
    if (rows != cols)
      throw math_error("Cannot invert a non-square matrix");
    const element_t zero = field.element(field.zero()),
                    one = field.element(field.one());

    std::vector<std::vector<element_t>> augmented_coeffs(
        rows, std::vector<element_t>(2 * cols, zero));
    for (size_t row = 0; row < rows; ++row) {
      for (size_t col = 0; col < cols; ++col)
        augmented_coeffs[row][col] = data[row][col];
      augmented_coeffs[row][cols + row] = one;
    }
    Matrix augmented(field, rows, 2 * cols, augmented_coeffs);

    const size_t rank = augmented.row_reduce();
    if (rank != rows || augmented[rows - 1][rows - 1] != one)
      throw math_error(
          "Supplied matrix was not full-rank and thus not invertible");
    std::vector<std::vector<element_t>> result_coeffs(
        rows, std::vector<element_t>(cols, zero));
    for (size_t row = 0; row < rows; ++row)
      for (size_t col = 0; col < cols; ++col)
        result_coeffs[row][col] = augmented[row][cols + col];

    const auto inverse = Matrix(field, rows, cols, result_coeffs);
    assert((*this) * inverse == eye(rows));
    assert(inverse * (*this) == eye(rows));
    return inverse;
  }

  bool operator==(const Matrix &other) const { return data == other.data; }
  bool operator!=(const Matrix &other) const { return !(*this == other); }

  // Printing
  friend std::ostream &operator<<(std::ostream &os, const Matrix &m) {
    os << "\n";
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
