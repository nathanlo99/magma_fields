
#pragma once

#include "gmp.h"
#include "gmp.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <vector>

struct Factor {
  integer_t base;
  uint64_t exp;
  Factor(const integer_t base, const uint64_t exp = 0) : base(base), exp(exp) {}
};
using Factorization = std::vector<Factor>;
using FactorizationMemo =
    std::map<std::pair<uint64_t, uint64_t>, Factorization>;

inline integer_t pk_minus_one(const integer_t p, const uint64_t k) {
  integer_t n;
  mpz_pow_ui(n.get_mpz_t(), p.get_mpz_t(), k);
  n -= 1;
  return n;
}

inline integer_t pk_plus_one(const integer_t p, const uint64_t k) {
  integer_t n;
  mpz_pow_ui(n.get_mpz_t(), p.get_mpz_t(), k);
  n += 1;
  return n;
}

inline void verify_factorization(const integer_t num,
                                 const Factorization &factorization) {
  integer_t product = 1;
  for (const Factor &factor : factorization) {
    assert(is_prime(factor.base));
    integer_t power;
    mpz_pow_ui(power.get_mpz_t(), factor.base.get_mpz_t(), factor.exp);
    product *= power;
  }
  assert(product == num);
}

constexpr std::string_view factorization_table_filename =
    "tables/prime_factorizations.txt";
inline FactorizationMemo factorization_memo;
inline void load_cache_factorizations() {
  std::ifstream ifs(factorization_table_filename);
  std::string line;
  while (std::getline(ifs, line)) {
    Factorization factorization;
    std::stringstream iss(line);
    uint64_t p, k, exp;
    std::string base;
    iss >> p >> k;
    std::pair<uint64_t, uint64_t> key = std::make_pair(p, k);
    while (iss >> base && iss >> exp) {
      factorization.push_back(Factor(integer_t(base), exp));
    }
    factorization_memo[key] = factorization;
    verify_factorization(pk_minus_one(from_uint(p), k), factorization);
  }
}

inline void cache_factorization(const uint64_t p, const uint64_t k,
                                const Factorization &factorization) {
  static std::ofstream ofs(factorization_table_filename, std::ios::app);
  ofs << p << " " << k;
  for (const Factor &factor : factorization) {
    ofs << " " << factor.base << " " << factor.exp;
  }
  ofs << std::endl;
  factorization_memo[std::make_pair(p, k)] = factorization;
}

inline Factorization combine_factorizations(const Factorization &a,
                                            const Factorization &b) {
  Factorization result;
  size_t ai = 0, bi = 0;
  while (ai < a.size() && bi < b.size()) {
    const integer_t a_base = a[ai].base, b_base = b[bi].base;
    const uint64_t a_exp = a[ai].exp, b_exp = b[bi].exp;
    if (a_base == b_base) {
      result.push_back(Factor(a_base, a_exp + b_exp));
      ai++;
      bi++;
    } else if (a_base < b_base) {
      result.push_back(a[ai]);
      ai++;
    } else {
      result.push_back(b[bi]);
      bi++;
    }
  }
  if (ai < a.size())
    result.insert(result.end(), a.begin() + ai, a.end());
  if (bi < b.size())
    result.insert(result.end(), b.begin() + bi, b.end());

  return result;
}

inline Factorization prime_factor(const integer_t num) {
  integer_t n = num;
  std::cout << "Factoring " << n << std::endl;

  if (is_prime(n))
    return {Factor(n, 1)};
  Factorization result;
  for (integer_t p = 2; p * p <= n; next_prime(p)) {
    if (n % p != 0)
      continue;
    const uint64_t exp =
        mpz_remove(n.get_mpz_t(), n.get_mpz_t(), p.get_mpz_t());
    result.push_back(Factor(p, exp));
    std::cout << "      ... " << n << std::endl;
    if (is_prime(n))
      break;
  }
  if (n != 1)
    result.push_back(Factor(n, 1));
  return result;
}

// Return a prime factorization of the number p^k - 1
inline Factorization factor_pk_minus_one(const integer_t p, const uint64_t k) {
  std::cout << "Factoring " << p << "^" << k << " - 1" << std::endl;
  const bool use_memo = p.fits_ulong_p();
  if (use_memo) {
    const std::pair<uint64_t, uint64_t> key = std::make_pair(to_uint(p), k);
    if (factorization_memo.count(key) > 0) {
      std::cout << "Grabbed factorization of " << p << "^" << k
                << " - 1 from cache" << std::endl;
      return factorization_memo[key];
    }
  }

  Factorization result;
  if (k % 2 == 0) {
    result = factor_pk_minus_one(p, k / 2);
    result =
        combine_factorizations(result, prime_factor(pk_plus_one(p, k / 2)));
  } else if (k % 3 == 0) {
    const integer_t x = pk_minus_one(p, k / 3) + 1; // p^{k/3}
    result = factor_pk_minus_one(p, k / 3);
    result = combine_factorizations(result, prime_factor(x * x + x + 1));
  } else if (k % 5 == 0) {
    const integer_t x = pk_minus_one(p, k / 5) + 1; // p^{k/5}
    result = factor_pk_minus_one(p, k / 5);
    result = combine_factorizations(
        result, prime_factor(x * x * x * x + x * x * x + x * x + x + 1));
  } else {
    result = prime_factor(pk_minus_one(p, k));
  }

  if (use_memo) {
    cache_factorization(to_uint(p), k, result);
  }

  verify_factorization(pk_minus_one(p, k), result);
  return result;
}

inline void print_factorization(const Factorization &factorization) {
  bool first = true;
  for (const Factor &factor : factorization) {
    if (first)
      first = false;
    else
      std::cout << " * ";
    std::cout << factor.base;
    if (factor.exp != 1)
      std::cout << "^" << factor.exp;
  }
  std::cout << std::endl;
}