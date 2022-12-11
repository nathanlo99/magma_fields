
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "polynomial.hpp"

template <class BaseField> struct ZechPolyField : Field<Polynomial<BaseField>> {
  using value_t = Polynomial<BaseField>;
  using base_element_t = FieldElement<BaseField>;
  using element_t = FieldElement<ZechPolyField>;

  const BaseField &base_field;
  Polynomial<BaseField> f;
  const integer_t p;
  const uint32_t k;
  const integer_t q;

  ZechPolyField(const BaseField &base_field, const Polynomial<BaseField> &f)
      : base_field(base_field), f(f), p(base_field.characteristic()),
        k(base_field.degree() * f.degree()), q(gmp::pow(p, k)) {
    if (base_field.type() != ZechFieldType)
      throw math_error()
          << "ZechPolyField expected SmallPrimeField as base field, got "
          << field_type_to_string(base_field.type());
    if (!f.is_irreducible_rabin())
      throw math_error()
          << "ZechPolyField expected irreducible polynomial, got " << f;
  }

  ZechPolyField(const BaseField &base_field, const char variable,
                const uint64_t k)
      : base_field(base_field), f(base_field, variable),
        p(base_field.characteristic()), k(base_field.degree() * k),
        q(gmp::pow(p, this->k)) {
    if (base_field.type() != ZechFieldType)
      throw math_error()
          << "ZechPolyField expected SmallPrimeField as base field, got "
          << field_type_to_string(base_field.type());
    do {
      f = Polynomial<BaseField>::sample(base_field, variable, k).monic();
    } while (f.coeffs[k] == f.zero || !f.is_irreducible_rabin());
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return k; }
  integer_t cardinality() const override { return q; }
  FieldType type() const override { return ZechPolyFieldType; }

  value_t zero() const override {
    return value_t(base_field, f.variable,
                   {base_field.element(base_field.zero())});
  }
  value_t one() const override {
    return value_t(base_field, f.variable,
                   {base_field.element(base_field.one())});
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
  element_t random_element() const {
    const auto random_poly =
        Polynomial<BaseField>::sample(base_field, f.variable, k);
    return element_t(*this, random_poly);
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
                                  const ZechPolyField &field) {
    return os << "ZechPolyField over [" << field.base_field << "] defined by "
              << field.f;
  }
};
