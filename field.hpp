
#pragma once

#include "gmp.hpp"

#include <iostream>

template <class Field> struct FieldElement;

template <typename T> struct Indexed {
  static inline size_t next_id = 0;
  size_t m_id;
  Indexed() { m_id = next_id++; }

  bool operator==(const Indexed<T> &other) const { return m_id == other.m_id; }
};

template <class Value> struct Field : Indexed<Field<Value>> {
  using value_t = Value;

  // Field properties
  virtual integer_t characteristic() const = 0;
  virtual integer_t degree() const = 0;
  virtual integer_t cardinality() const = 0;

  // Element constructors
  virtual value_t zero() const = 0;
  virtual value_t one() const = 0;
  virtual value_t integer(const integer_t number) const = 0;

  // Unary and binary operations
  // Addition
  virtual value_t neg(const value_t a) const = 0;
  virtual value_t add(const value_t a, const value_t b) const = 0;
  virtual value_t sub(const value_t a, const value_t b) const {
    return add(a, neg(b));
  }

  // Multiplication
  virtual value_t inv(const value_t a) const = 0;
  virtual value_t mul(const value_t a, const value_t b) const = 0;
  virtual value_t div(const value_t a, const value_t b) const {
    return mul(a, inv(b));
  }
  virtual value_t pow(const value_t a, const integer_t exp) const {
    value_t base = a, result = one();
    integer_t expt = exp;
    if (expt < 0) {
      base = inv(base);
      expt = -expt;
    }
    while (expt > 0) {
      if (expt % 2 == 1)
        result = mul(result, base);
      base = mul(base, base);
      expt /= 2;
    }
    return result;
  }

  // Element equality
  virtual bool eq(const value_t a, const value_t b) const = 0;

  virtual std::string to_string(const value_t a) const = 0;
};

template <typename Field> struct FieldElement {
  using element_t = FieldElement<Field>;
  using value_t = typename Field::value_t;

  const Field &field;
  value_t value;

  FieldElement(const Field &field, const value_t value)
      : field(field), value(value) {}
  FieldElement(const FieldElement &other)
      : field(other.field), value(other.value) {}
  FieldElement &operator=(const FieldElement &other) {
    if (field != other.field)
      throw math_error("Field element cannot be assigned to field "
                       "element with different base field");
    value = other.value;
    return *this;
  }

  element_t zero(const Field &field) const {
    return FieldElement(field, field.zero());
  }
  element_t one(const Field &field) const {
    return FieldElement(field, field.one());
  }

  element_t operator-() const { return element_t(field, field.neg(value)); }
  element_t operator+(const element_t &other) const {
    return element_t(field, field.add(value, other.value));
  }
  element_t operator-(const element_t &other) const {
    return element_t(field, field.sub(value, other.value));
  }
  element_t &operator+=(const element_t &other) {
    return *this = *this + other;
  }
  element_t &operator-=(const element_t &other) {
    return *this = *this - other;
  }
  friend element_t operator+(const element_t &a, const integer_t b) {
    return a + a.field(b);
  }
  friend element_t operator+(const integer_t a, const element_t &b) {
    return b + a;
  }

  element_t inv() const { return element_t(field, field.inv(value)); }
  element_t operator*(const element_t &other) const {
    return element_t(field, field.mul(value, other.value));
  }
  element_t operator/(const element_t &other) const {
    return element_t(field, field.div(value, other.value));
  }
  element_t operator^(const integer_t exp) const {
    return element_t(field, field.pow(value, exp));
  }

  bool operator==(const element_t &other) const {
    return field.eq(value, other.value);
  }

  friend element_t operator*(const integer_t a, const element_t &b) {
    return b.field(a) * b;
  }
  friend element_t operator*(const element_t &a, const integer_t b) {
    return a * a.field(b);
  }

  friend bool operator==(const integer_t a, const element_t &b) {
    return b.field.eq(b.field.integer(a), b.value);
  }
  friend bool operator==(const element_t &a, const integer_t b) {
    return a.field.eq(a.value, a.field.integer(b));
  }

  friend std::ostream &operator<<(std::ostream &os, const FieldElement &val) {
    return os << val.field.to_string(val.value);
  }
};
