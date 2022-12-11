
#pragma once

#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "polynomial.hpp"

template <class BaseField>
struct PrimePolyField : Field<Polynomial<BaseField>> {
  using value_t = Polynomial<BaseField>;
  using base_element_t = FieldElement<BaseField>;
  using element_t = FieldElement<PrimePolyField>;

  const BaseField &base_field;
  const Polynomial<BaseField> &f;
  const integer_t p;
  const uint32_t k;
  const integer_t q;

  PrimePolyField(const BaseField &base_field, const Polynomial<BaseField> &f)
      : base_field(base_field), f(f), p(base_field.characteristic()),
        k(base_field.degree() * f.degree()), q(gmp::pow(p, k)) {
    // TODO: Check that f is irreducible

    if (base_field.type() != SmallPrimeFieldType)
      throw math_error()
          << "PrimePolyField expected SmallPrimeField as base field, got "
          << field_type_to_string(base_field.type());
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return k; }
  integer_t cardinality() const override { return q; }
  FieldType type() const override { return PrimePolyFieldType; }

  value_t zero() const override {
    return value_t(base_field, f.variable, {base_field.element(0)});
  }
  value_t one() const override {
    return value_t(base_field, f.variable, {base_field.element(1)});
  }
  value_t integer(const integer_t number) const override {
    const base_element_t coeff = base_field(number);
    return value_t(base_field, f.variable, {coeff});
  }
  element_t operator()(const integer_t num) const {
    return element_t(*this, integer(num));
  }
  element_t element(const value_t value) const {
    return element_t(*this, value);
  }
  element_t primitive_element() const {
    return element_t(*this, value_t(base_field, f.variable));
  }

  value_t neg(const value_t a) const override { return (f - a) % f; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % f;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a + f - b) % f;
  }

  value_t inv(const value_t a) const override { return inv_mod(a, f); }
  value_t mul(const value_t a, const value_t b) const override {
    return (a * b) % f;
  }
  value_t div(const value_t a, const value_t b) const override {
    return mul(a, inv(b));
  }
  // value_t pow(const value_t a, const integer_t exp) const; // Default

  bool eq(const value_t a, const value_t b) const override { return a == b; }

  std::string to_string(const value_t a) const override {
    return a.to_string();
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const PrimePolyField &field) {
    return os << "PrimePolyField over [" << field.base_field << "] defined by "
              << field.f;
  }
};
