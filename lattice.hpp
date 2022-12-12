
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

struct Lattice {
  const integer_t p;
  std::vector<std::shared_ptr<AbstractField>> fields;
  Lattice(const integer_t p) : p(p) { add_field(1); }

  // When creating a new field, there are a few options:
  // 1. If the field is a prime field, then
  //    a. If p < 2^16,                      -> SmallPrimeField
  //    b. If p < 2^32,                      -> MediumPrimeField
  //    c. Otherwise,                        -> LargePrimeField
  // 2. If the cardinality is at most 2^20,  -> ZechField
  // 3. If there exists a field of size p^l
  //      with l a non-trivial factor of k,
  //      and p^l <= 2^20, create a
  //      ZechField S with degree l, and     - - > ZechField (recursively)
  //      create a ZechPolyField of degree
  //      (k/l) over it                      -> ZechPolyField
  // 4. If p < 2^16                          -> PrimePolyField
  // 5. Finally                              -> GeneralPolyField
  std::shared_ptr<AbstractField> add_field(const uint64_t k) {
    std::cout << "add_field called on lattice with characteristic " << p
              << " and degree " << k << std::endl;

    // 1. Prime fields
    if (k == 1) {
      // If we already have a prime field, just return that
      if (fields.size() > 0)
        return fields[0];
      if (p < (1_mpz << 16))
        fields.push_back(std::make_shared<SmallPrimeField>(p));
      else if (p < (1_mpz << 32))
        fields.push_back(std::make_shared<MediumPrimeField>(p));
      else
        fields.push_back(std::make_shared<LargePrimeField>(p));
      return fields.back();
    }

    // Otherwise, if we get here, we better already have a prime field
    assert(fields.size() >= 1);
    AbstractField *prime_field = fields[0].get();

    const integer_t q = gmp::pow(p, k);
    // 2. If the cardinality is at most 2^20, ZechField
    if (q <= (1_mpz << 20)) {
      switch (prime_field->type()) {
      case FieldType::SmallPrime: {
        const SmallPrimeField &P =
            dynamic_cast<SmallPrimeField &>(*prime_field);
        fields.push_back(
            std::make_shared<ZechField<SmallPrimeField>>(P, "z", k));
        return fields.back();
      }

      case FieldType::MediumPrime: {
        const MediumPrimeField &P =
            dynamic_cast<MediumPrimeField &>(*prime_field);
        fields.push_back(
            std::make_shared<ZechField<MediumPrimeField>>(P, "z", k));
        return fields.back();
      }

      case FieldType::LargePrime: {
        const LargePrimeField &P =
            dynamic_cast<LargePrimeField &>(*prime_field);
        fields.push_back(
            std::make_shared<ZechField<LargePrimeField>>(P, "z", k));
        return fields.back();
      }

      default:
        throw math_error()
            << "Prime field was not actually of prime field type";
      }
      return nullptr;
    }

    // 3. Two-step optimized representation
    uint32_t best_ell = 1;
    for (uint32_t ell = 2; ell < k; ++ell) {
      if (gmp::pow(p, ell) > (1_mpz << 20))
        break;
      if (k % ell == 0)
        best_ell = ell;
    }

    if (best_ell == 1) {
      // None of the factors of n had p^ell <= 2^20, so return PrimePolyField
      switch (prime_field->type()) {
      case FieldType::SmallPrime: {
        const SmallPrimeField &S =
            dynamic_cast<SmallPrimeField &>(*prime_field);
        fields.push_back(
            std::make_shared<PrimePolyField<SmallPrimeField>>(S, "z", k));
        return fields.back();
      }

      case FieldType::MediumPrime: {
        const MediumPrimeField &S =
            dynamic_cast<MediumPrimeField &>(*prime_field);
        fields.push_back(
            std::make_shared<PrimePolyField<MediumPrimeField>>(S, "z", k));
        return fields.back();
      }

      case FieldType::LargePrime: {
        const LargePrimeField &S =
            dynamic_cast<LargePrimeField &>(*prime_field);
        fields.push_back(
            std::make_shared<PrimePolyField<LargePrimeField>>(S, "z", k));
        return fields.back();
      }

      default:
        throw math_error()
            << "Prime field was not actually of prime field type";
      }
    } else {
      // Otherwise, create a two-step optimized representation
      const auto S_tmp = add_field(best_ell);
      assert(S_tmp->type() == FieldType::Zech);

      switch (prime_field->type()) {
      case FieldType::SmallPrime: {
        const ZechField<SmallPrimeField> &S =
            dynamic_cast<ZechField<SmallPrimeField> &>(*S_tmp);
        fields.push_back(
            std::make_shared<ZechPolyField<ZechField<SmallPrimeField>>>(
                S, "z", k / best_ell));
        return fields.back();
      }

      case FieldType::MediumPrime: {
        const ZechField<MediumPrimeField> &S =
            dynamic_cast<ZechField<MediumPrimeField> &>(*S_tmp);
        fields.push_back(
            std::make_shared<ZechPolyField<ZechField<MediumPrimeField>>>(
                S, "z", k / best_ell));
        return fields.back();
      }

      case FieldType::LargePrime: {
        const ZechField<LargePrimeField> &S =
            dynamic_cast<ZechField<LargePrimeField> &>(*S_tmp);
        fields.push_back(
            std::make_shared<ZechPolyField<ZechField<LargePrimeField>>>(
                S, "z", k / best_ell));
        return fields.back();
      }

      default:
        throw math_error()
            << "Prime field was not actually of prime field type";
      }
    }

    throw math_error()
        << "Unimplemented code path: could not add field with cardinality " << p
        << "^" << k;
  }

  friend std::ostream &operator<<(std::ostream &os, const Lattice &lattice) {
    for (const auto &field : lattice.fields) {
      std::cout << " - " << field << ": " << field->to_string() << std::endl;
    }
    return os;
  }
};

struct LatticeManager {
  std::map<integer_t, std::shared_ptr<Lattice>> lattices;

  std::shared_ptr<AbstractField> FiniteField(const integer_t p,
                                             const uint64_t k = 1) {
    if (p < 0 || !gmp::is_prime(p))
      throw math_error() << "FiniteField expected positive prime, got " << p;
    if (k == 0)
      throw math_error() << "FiniteField expected positive degree, got " << k;
    if (lattices.count(p) == 0)
      lattices[p] = std::make_shared<Lattice>(p);
    const auto lattice = lattices.at(p);
    return lattice->add_field(k);
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const LatticeManager &manager) {
    os << std::endl;
    os << "------- LATTICES -------" << std::endl;
    for (const auto &[p, lattice] : manager.lattices) {
      os << "Lattice for prime " << p << ": " << std::endl;
      os << *lattice << std::endl << std::endl;
    }
    os << "------------------------" << std::endl;
    return os;
  }
};
