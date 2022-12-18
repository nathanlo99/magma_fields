
#include "lattice.hpp"
#include "field.hpp"
#include "gmp.hpp"
#include "large_prime_field.hpp"

#include <utility>

std::map<integer_t, Lattice<SmallPrimeField>> small_prime_lattices;
std::map<integer_t, Lattice<MediumPrimeField>> medium_prime_lattices;
std::map<integer_t, Lattice<LargePrimeField>> large_prime_lattices;

template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_prime_field() {
  fields.push_back(std::make_shared<PrimeField>(p));
  return fields.back();
}

template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_zech_field(const std::string &variable,
                                    const uint64_t k)
  requires(!std::is_same<PrimeField, LargePrimeField>::value)
{
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  const auto P = dynamic_cast<PrimeField &>(*prime_field);
  fields.push_back(std::make_shared<ZechField<PrimeField>>(P, variable, k));
  return fields.back();
}

template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_prime_poly_field(const std::string &variable,
                                          const uint64_t k) {
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  const auto P = dynamic_cast<PrimeField &>(*prime_field);
  fields.push_back(
      std::make_shared<PrimePolyField<PrimeField>>(P, variable, k));
  return fields.back();
}

template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_zech_poly_field(LatticeField<PrimeField> *base,
                                         const std::string &variable,
                                         const uint64_t k)
  requires(!std::is_same<PrimeField, LargePrimeField>::value)
{
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  const auto S = dynamic_cast<ZechField<PrimeField> &>(*base);
  fields.push_back(
      std::make_shared<ZechPolyField<ZechField<PrimeField>>>(S, variable, k));
  return fields.back();
}

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
template <class PrimeField>
std::shared_ptr<LatticeField<PrimeField>>
Lattice<PrimeField>::add_field(const uint64_t k, const std::string &_variable) {
  const std::string variable =
      _variable != "" ? _variable : "w" + std::to_string(k);
  assert(fields.size() > 0);

  log() << "add_field called on lattice with characteristic " << p
        << " and degree " << k << std::endl;

  // 1. Prime field
  if (k == 1)
    return fields[0];

  const integer_t q = gmp::pow(p, k);

  if constexpr (std::is_same<PrimeField, LargePrimeField>::value) {
    return add_prime_poly_field(variable, k);
  } else {
    // 2. If the cardinality is at most 2^20, ZechField
    if (q <= (1_mpz << 20))
      return add_zech_field(variable, k);

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
      return add_prime_poly_field(variable, k);
    } else {
      // Otherwise, create a two-step optimized representation
      const auto S_tmp = add_field(best_ell);
      assert(S_tmp->type() == FieldType::Zech);
      return add_zech_poly_field(S_tmp.get(), variable, k / best_ell);
    }
  }
}
