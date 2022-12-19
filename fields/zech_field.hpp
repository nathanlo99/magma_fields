
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "logger.hpp"
#include "polynomial.hpp"
#include "polynomial_factorization.hpp"
#include "random.hpp"
#include "vector.hpp"

#include <iostream>

template <class BaseField>
uint32_t get_polynomial_index(const Polynomial<BaseField> &f) {
  const auto base_field = f.field;
  const uint32_t p = gmp::to_uint(base_field.cardinality());
  uint32_t result = 0;
  for (int i = f.degree(); i >= 0; --i) {
    const uint32_t idx = gmp::to_uint(f.coeffs[i].as_integer());
    result = result * p + idx;
  }
  return result;
}

template <class BaseField>
std::vector<uint32_t>
compute_zech_table(const uint32_t p, const uint32_t k, const uint32_t q,
                   const Polynomial<BaseField> &f,
                   const Polynomial<BaseField> &g,
                   std::vector<uint32_t> &constants, uint32_t &var) {
  log() << "Computing Zech table with f = " << f << " and g = " << g
        << std::endl;
  // Compute the Zech logarithm table, find s[r] such that x^s[r] = x^r + 1
  std::vector<uint32_t> result(q, -1);
  const auto base_field = f.field;
  const auto one = f.one_poly();

  Polynomial<BaseField> pow = one;
  std::vector<uint32_t> log_x(q, q);
  std::vector<uint32_t> next_index(q, q);
  constants = std::vector<uint32_t>(p);
  constants[0] = q - 1;

  for (uint32_t exp = 0; exp + 1 < q; ++exp) {
    if (exp > 0 && pow == one)
      throw math_error() << "Given polynomial was not primitive, g^" << exp
                         << " = 1";
    const uint32_t idx = get_polynomial_index(pow);

    // Cache special cases
    if (pow == f.var_poly())
      var = exp;
    if (pow.degree() == 0)
      constants[pow[0].value] = exp;

    next_index[exp] = get_polynomial_index(pow + one);
    log_x[idx] = exp;
    pow = (pow * g) % f;
  }

  // Special cases
  result[q - 1] = 0;
  const uint32_t log_negative_one = p == 2 ? 0 : (q - 1) / 2;
  result[log_negative_one] = q - 1;

  for (uint32_t i = 0; i < q; ++i) {
    if (result[i] >= q)
      result[i] = log_x[next_index[i]];
  }

  log() << "Done computing Zech table" << std::endl;
  return result;
}

// For prime extension fields with cardinality at most 2^20
//
// Elements are represented as powers of a primitive element, and a table of
// Zech logarithms is stored for use in addition.
template <class BaseField>
struct ZechField : Field<uint32_t, typename BaseField::prime_field_t> {
  using value_t = uint32_t;
  using element_t = FieldElement<ZechField>;
  using prime_field_t = typename BaseField::prime_field_t;
  using ground_field_t = BaseField;
  using vector_t = Vector<prime_field_t>;

  uint32_t p, k, q;
  const BaseField &base_field;
  std::vector<value_t> zech_table;
  std::vector<value_t> constants;
  Polynomial<BaseField> f, g;
  value_t var; // The value of the polynomial 'x'

  ZechField(const BaseField &base_field, const Polynomial<BaseField> &f)
      : p(gmp::to_uint(base_field.characteristic())), k(f.degree()),
        q(gmp::to_uint(gmp::pow(gmp::from_uint(p), k))), base_field(base_field),
        f(f), g(f.field, f.variable) {
    if (base_field != f.field)
      throw math_error()
          << "Polynomial base field did not match provided base_field";
    if (base_field.degree() != 1)
      throw math_error()
          << "Zech field expected base field with prime order, got " << f.field;
    if (q > (1_mpz << 20))
      throw math_error() << "Zech field expected cardinality at most 2^20, got "
                         << p << "^" << k << " = " << q;
    zech_table = compute_zech_table(this->p, k, q, f, g, constants, var);
  }

  ZechField(const BaseField &base_field, const std::string &variable,
            const uint64_t k)
      : p(gmp::to_uint(base_field.characteristic())), k(k),
        q(gmp::to_uint(gmp::pow(gmp::from_uint(p), k))), base_field(base_field),
        f(base_field, variable), g(base_field, variable) {
    if (base_field != f.field)
      throw math_error()
          << "Polynomial base field did not match provided base_field";
    if (base_field.degree() != 1)
      throw math_error()
          << "Zech field expected base field with prime order, got " << f.field;
    if (q > (1_mpz << 20))
      throw math_error() << "Zech field expected cardinality at most 2^20, got "
                         << p << "^" << k << " = " << q;

    log() << "Creating a ZechField extension of degree " << k
          << " over base_field " << base_field << std::endl;
    f = get_irreducible_polynomial(base_field, variable, k);
    log() << "Found an irreducible polynomial of degree " << k << ": '" << f
          << "'" << std::endl;

    // 2. Generate a polynomial g with full order in Z_q[x]/<f>
    do {
      g = random_polynomial<false, true>(base_field, variable, k);
    } while (!g.is_primitive_mod(f));
    log() << "Found the primitive polynomial '" << g << "'" << std::endl;

    // 3. Generate the Zech Table from it
    zech_table = compute_zech_table(this->p, k, q, f, g, constants, var);
  }

  integer_t characteristic() const override { return p; }
  uint32_t degree() const override { return k; }
  integer_t cardinality() const override { return q; }
  FieldType type() const override { return FieldType::Zech; }
  const prime_field_t &prime_field() const override { return base_field; }
  const ground_field_t &ground_field() const { return base_field; }

  // Element constructors
  value_t zero() const override { return q - 1; }
  value_t one() const override { return 0; }
  value_t integer(const integer_t number) const override {
    value_t result = zero(), base = one();
    // This conversion will fit since p fits in a uint32_t
    uint32_t num = gmp::to_uint(gmp::unsigned_mod(number, p));
    while (num != 0) {
      if (num % 2 == 1)
        result = add(result, base);
      base = add(base, base);
      num /= 2;
    }
    return result;
  }
  element_t operator()(const integer_t num) const {
    return element_t(*this, integer(num));
  }
  element_t element(const value_t value) const {
    return element_t(*this, value);
  }
  element_t primitive_element() const { return element_t(*this, 1); }
  // NOTE: This is not cryptographically secure or even uniform.
  element_t random_element() const {
    return element_t(*this, random_uint64() % q);
  }

  // Vector space structure
  element_t generating_element() const { return element_t(*this, var); }
  vector_t to_vector(const element_t elem) const {
    return vector_t(f.field, degree(), to_polynomial(elem.value).coeffs);
  }
  element_t from_vector(const vector_t vec) const {
    const element_t x = element_t(*this, var);
    element_t result = element_t(*this, zero());
    for (int d = vec.size - 1; d >= 0; --d) {
      const uint32_t coeff = vec.data[d].value;
      result = result * x + element_t(*this, constants[coeff]);
    }
    return result;
  }

  // Field operations
  value_t neg(const value_t a) const override {
    if (p == 2 || a == zero())
      return a;
    const value_t negative_one = (q - 1) / 2;
    return (a + negative_one) % (q - 1);
  }
  value_t add(const value_t a, const value_t b) const override {
    if (a == zero())
      return b;
    if (b == zero())
      return a;
    if (a < b) {
      // g^a + g^b = g^a (1 + g^{b-a}) = g^a ( g^s[b-a] )
      return mul(a, zech_table[b - a]);
    } else {
      // Similarly
      return mul(b, zech_table[a - b]);
    }
  }
  value_t sub(const value_t a, const value_t b) const override {
    return add(a, neg(b));
  }

  value_t inv(const value_t a) const override {
    if (a == zero())
      throw math_error("Division by zero");
    return a == 0 ? 0 : q - 1 - a;
  }
  value_t mul(const value_t a, const value_t b) const override {
    if (a == zero() || b == zero())
      return zero();
    return (a + b) % (q - 1);
  }
  value_t div(const value_t a, const value_t b) const override {
    if (b == zero())
      throw math_error("Division by zero");
    if (a == zero())
      return zero();
    return (a + (q - 1) - b) % (q - 1);
  }
  value_t pow(const value_t a, const integer_t exp) const override {
    if (a == zero())
      return zero();
    const integer_t modded_exp = gmp::unsigned_mod(exp, q - 1);
    return (a * gmp::to_uint(modded_exp)) % (q - 1);
  }

  virtual bool eq(const value_t a, const value_t b) const override {
    return a == b;
  }

  Polynomial<BaseField> to_polynomial(const value_t a) const {
    if (a == q - 1)
      return f.zero_poly();
    return pow_mod(g, a, f);
  }

  std::string value_to_string(const value_t a) const override {
    return to_polynomial(a).to_string();
  }

  std::string to_string() const override {
    std::stringstream ss;
    ss << "Degree " << f.degree() << " extension over " << f.field.to_string()
       << " with defining polynomial '" << f << "'";
    return ss.str();
  }
};
