
#include "lattice.hpp"
#include "field.hpp"
#include "gmp.hpp"

#include <utility>

std::map<integer_t, Lattice> lattices;

std::shared_ptr<AbstractField> Lattice::add_prime_field() {
  if (p.fits_ushort_p())
    fields.push_back(std::make_shared<SmallPrimeField>(p));
  else if (p.fits_uint_p())
    fields.push_back(std::make_shared<MediumPrimeField>(p));
  else
    fields.push_back(std::make_shared<LargePrimeField>(p));
  return fields.back();
}

std::shared_ptr<AbstractField>
Lattice::add_zech_field(const std::string &variable, const uint64_t k) {
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  switch (prime_field->type()) {
  case FieldType::SmallPrime: {
    const SmallPrimeField &P = dynamic_cast<SmallPrimeField &>(*prime_field);
    fields.push_back(
        std::make_shared<ZechField<SmallPrimeField>>(P, variable, k));
    return fields.back();
  }

  case FieldType::MediumPrime: {
    const MediumPrimeField &P = dynamic_cast<MediumPrimeField &>(*prime_field);
    fields.push_back(
        std::make_shared<ZechField<MediumPrimeField>>(P, variable, k));
    return fields.back();
  }

  case FieldType::LargePrime: {
    const LargePrimeField &P = dynamic_cast<LargePrimeField &>(*prime_field);
    fields.push_back(
        std::make_shared<ZechField<LargePrimeField>>(P, variable, k));
    return fields.back();
  }

  default:
    throw math_error() << "Prime field was not actually prime";
  }
}

std::shared_ptr<AbstractField>
Lattice::add_prime_poly_field(const std::string &variable, const uint64_t k) {
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  switch (prime_field->type()) {
  case FieldType::SmallPrime: {
    const SmallPrimeField &P = dynamic_cast<SmallPrimeField &>(*prime_field);
    fields.push_back(
        std::make_shared<PrimePolyField<SmallPrimeField>>(P, variable, k));
    return fields.back();
  }

  case FieldType::MediumPrime: {
    const MediumPrimeField &P = dynamic_cast<MediumPrimeField &>(*prime_field);
    fields.push_back(
        std::make_shared<PrimePolyField<MediumPrimeField>>(P, variable, k));
    return fields.back();
  }

  case FieldType::LargePrime: {
    const LargePrimeField &P = dynamic_cast<LargePrimeField &>(*prime_field);
    fields.push_back(
        std::make_shared<PrimePolyField<LargePrimeField>>(P, variable, k));
    return fields.back();
  }

  default:
    throw math_error() << "Prime field was not actually of prime field type";
  }
}

std::shared_ptr<AbstractField>
Lattice::add_zech_poly_field(AbstractField *base, const std::string &variable,
                             const uint64_t k) {
  assert(fields.size() > 0);
  const auto prime_field = fields[0];
  switch (prime_field->type()) {
  case FieldType::SmallPrime: {
    const ZechField<SmallPrimeField> &S =
        dynamic_cast<ZechField<SmallPrimeField> &>(*base);
    fields.push_back(
        std::make_shared<ZechPolyField<ZechField<SmallPrimeField>>>(S, variable,
                                                                    k));
    return fields.back();
  }

  case FieldType::MediumPrime: {
    const ZechField<MediumPrimeField> &S =
        dynamic_cast<ZechField<MediumPrimeField> &>(*base);
    fields.push_back(
        std::make_shared<ZechPolyField<ZechField<MediumPrimeField>>>(
            S, variable, k));
    return fields.back();
  }

  case FieldType::LargePrime: {
    const ZechField<LargePrimeField> &S =
        dynamic_cast<ZechField<LargePrimeField> &>(*base);
    fields.push_back(
        std::make_shared<ZechPolyField<ZechField<LargePrimeField>>>(S, variable,
                                                                    k));
    return fields.back();
  }

  default:
    throw math_error() << "Prime field was not actually of prime field type";
  }
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
std::shared_ptr<AbstractField>
Lattice::add_field(const uint64_t k, const std::string &_variable) {
  const std::string variable =
      _variable != "" ? _variable : "w" + std::to_string(k);
  assert(fields.size() > 0);

  log() << "add_field called on lattice with characteristic " << p
        << " and degree " << k << std::endl;

  // 1. Prime field
  if (k == 1)
    return fields[0];

  const integer_t q = gmp::pow(p, k);
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

  throw math_error()
      << "Unimplemented code path: could not add field with cardinality " << p
      << "^" << k;
}
