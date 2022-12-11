
#pragma once

#include "field.hpp"
#include "gmp.hpp"
#include "large_prime_field.hpp"
#include "medium_prime_field.hpp"
#include "small_prime_field.hpp"
#include "zech_field.hpp"
#include "zech_poly_field.hpp"

#include <array>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

// Returns the largest k such that p^k <= 2^20 = (1024)^2
inline uint32_t get_largest_power_for_zech_field(const integer_t p) {
  if (p >= 1024)
    return 1;
  if (p >= 103)
    return 2;
  if (p >= 37)
    return 3;
  if (p >= 17)
    return 4;
  if (p >= 11)
    return 5;
  constexpr std::array<uint32_t, 8> lookup_table = {0, 0, 20, 12, 0, 8, 0, 7};
  return lookup_table[gmp::to_uint(p)];
}

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
    // TODO: Construct and add a field with degree k
    std::cout << "add_field called on lattice with characteristic " << p
              << " and degree " << k << std::endl;

    // 1. Prime fields
    if (k == 1) {
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
    const auto prime_field = fields[0];

    const integer_t q = gmp::pow(p, k);
    // 2.
    if (q <= (1_mpz << 20)) {
      // TODO: Fix this or write a ZechField constructor which finds an
      // irreducible polynomial automatically
      // fields.push_back(std::make_shared<ZechField>(*prime_field, k));
      // return fields.back();
      return nullptr;
    }

    uint32_t best_ell = 1;
    for (uint32_t ell = 2; ell < k; ++ell) {
      if (gmp::pow(p, ell) > (1_mpz << 20))
        break;
      if (k % ell == 0)
        best_ell = ell;
    }
    if (best_ell == 1) {
      // TODO: Return PrimePolyField
    } else {
      const auto S = add_field(best_ell);
      // assert(S->type() == ZechFieldType);
      // TODO: Return ZechPolyField of degree k / best_ell over S
    }
    return nullptr;
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
};
