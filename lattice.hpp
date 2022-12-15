
#pragma once

#include "error.hpp"
#include "field.hpp"
#include "gmp.hpp"
#include "large_prime_field.hpp"
#include "medium_prime_field.hpp"
#include "prime_poly.hpp"
#include "small_prime_field.hpp"
#include "zech_field.hpp"
#include "zech_poly_field.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

struct Lattice;
extern std::map<integer_t, Lattice> lattices;

struct Lattice {
  const integer_t p;
  std::vector<std::shared_ptr<AbstractField>> fields;
  Lattice(const integer_t p) : p(p) { add_prime_field(); }

  std::shared_ptr<AbstractField> add_prime_field();
  std::shared_ptr<AbstractField> add_zech_field(const std::string &variable,
                                                const uint64_t k);
  std::shared_ptr<AbstractField>
  add_prime_poly_field(const std::string &variable, const uint64_t k);
  std::shared_ptr<AbstractField>
  add_zech_poly_field(AbstractField *base, const std::string &variable,
                      const uint64_t k);

  std::shared_ptr<AbstractField> add_field(const uint64_t k,
                                           const std::string &_variable = "");

  friend std::ostream &operator<<(std::ostream &os, const Lattice &lattice) {
    for (const auto &field : lattice.fields) {
      os << " - " << field << ": " << std::setw(20) << std::setfill(' ')
         << field_type_to_string(field->type()) << ": ("
         << field->characteristic() << "^" << field->degree() << ") "
         << field->to_string() << std::endl;
    }
    return os;
  }
};

inline void print_lattices(std::ostream &os) {
  os << std::endl;
  os << "------- LATTICES -------" << std::endl;
  for (const auto &[p, lattice] : lattices) {
    os << "Lattice for prime " << p << ": " << std::endl;
    os << lattice << std::endl << std::endl;
  }
  os << "------------------------" << std::endl;
}

inline std::shared_ptr<AbstractField> FiniteField(const integer_t p,
                                                  const uint64_t k = 1) {
  if (p < 0 || !gmp::is_prime(p))
    throw math_error() << "FiniteField expected positive prime, got " << p;
  if (k == 0)
    throw math_error() << "FiniteField expected positive degree, got " << k;
  if (lattices.count(p) == 0)
    lattices.emplace(p, Lattice(p));
  return lattices.at(p).add_field(k);
}

inline void Embed(const std::shared_ptr<AbstractField> &E,
                  const std::shared_ptr<AbstractField> &F) {
  //
}
