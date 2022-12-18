
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "polynomial.hpp"
#include "polynomial_factorization.hpp"

template <class BaseField>
struct ZechPolyField
    : Field<Polynomial<BaseField>, typename BaseField::prime_field_t> {
  using value_t = Polynomial<BaseField>;
  using base_element_t = FieldElement<BaseField>;
  using element_t = FieldElement<ZechPolyField>;
  using prime_field_t = typename BaseField::prime_field_t;
  using ground_field_t = BaseField;
  using vector_t = Vector<prime_field_t>;

  const BaseField &base_field;
  Polynomial<BaseField> f;
  const integer_t p;
  const uint32_t k;
  const integer_t q;

  ZechPolyField(const BaseField &base_field, const Polynomial<BaseField> &f)
      : base_field(base_field), f(f), p(base_field.characteristic()),
        k(base_field.degree() * f.degree()), q(gmp::pow(p, k)) {
    if (base_field.type() != FieldType::Zech)
      throw math_error()
          << "ZechPolyField expected SmallPrimeField as base field, got "
          << field_type_to_string(base_field.type());
    if (!f.is_irreducible_rabin())
      throw math_error()
          << "ZechPolyField expected irreducible polynomial, got " << f;
  }

  ZechPolyField(const BaseField &base_field, const std::string &variable,
                const uint64_t k)
      : base_field(base_field), f(base_field, variable),
        p(base_field.characteristic()), k(base_field.degree() * k),
        q(gmp::pow(p, this->k)) {
    if (base_field.type() != FieldType::Zech)
      throw math_error()
          << "ZechPolyField expected SmallPrimeField as base field, got "
          << field_type_to_string(base_field.type());
    log() << "Creating degree " << k << " extension over " << base_field
          << std::endl;
    f = get_irreducible_polynomial(base_field, variable, k);
    log() << "Found the irreducible polynomial " << f << std::endl;
    log() << "Constructed PrimePolyField over " << base_field << " with degree "
          << k << std::endl;
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return k; }
  integer_t cardinality() const override { return q; }
  FieldType type() const override { return FieldType::ZechPoly; }
  const prime_field_t &prime_field() const override {
    return base_field.prime_field();
  }
  const ground_field_t &ground_field() const { return base_field; }

  // Element constructors
  value_t zero() const override { return f.zero_poly(); }
  value_t one() const override { return f.one_poly(); }
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
    static value_t cache = f.zero_poly();
    if (cache == f.zero_poly()) {
      do {
        cache = random_element().value;
      } while (this->is_primitive(cache));
    }
    return element_t(*this, cache);
  }
  element_t random_element() const {
    const auto random_poly =
        random_polynomial<false, false>(base_field, f.variable, k - 1);
    return element_t(*this, random_poly);
  }

  // Vector space structure
  element_t generating_element() const { return primitive_element(); }
  // TODO: Implement to_vector and from_vector

  // Field structure
  value_t neg(const value_t a) const override { return (-a) % f; }
  value_t add(const value_t a, const value_t b) const override {
    return (a + b) % f;
  }
  value_t sub(const value_t a, const value_t b) const override {
    return (a - b) % f;
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

  std::string value_to_string(const value_t a) const override {
    return a.to_string();
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const ZechPolyField &field) {
    return os << "Degree " << field.f.degree() << " extension with cardinality "
              << field.characteristic() << "^" << field.degree() << " = "
              << field.cardinality();
  }

  std::string to_string() const override {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }
};
