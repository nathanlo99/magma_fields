
#include "gmp.hpp"
#include "large_prime_field.hpp"
#include "medium_prime_field.hpp"
#include "prime_factorization.hpp"
#include "small_prime_field.hpp"
#include "timing.hpp"
#include "zech_field.hpp"
#include "zech_poly_field.hpp"

#include <gmpxx.h>
#include <iostream>

inline void init_all() { load_cache_factorizations(); }

int main(int argc, char *argv[]) {
  init_all();

  timeit("Prime field of size near 2^16", []() {
    const auto F = SmallPrimeField(65521);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  });

  timeit("Prime field of medium size", []() {
    const auto F = MediumPrimeField(100003);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  });

  timeit("Prime field of size over 2^32", []() {
    const auto F = LargePrimeField(4294967311);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  });

  timeit("Field of cardinality 3^3", []() {
    const auto F3 = SmallPrimeField(3);
    const auto f = Polynomial(F3, 'x', {1, 2, 0, 1});
    const auto Z = ZechField(F3, f);
    std::cout << "Z = " << Z << std::endl;
    const auto x = Z.generator();
    for (int i = 0; i < Z.cardinality(); ++i) {
      std::cout << "x^" << i << " = " << (x ^ i) << std::endl;
    }
  });

  // timeit("Field of cardinality 2^20", []() {
  //   const auto F2 = SmallPrimeField(2);
  //   const auto x = Polynomial(F2, 'x');
  //   const auto f = (x ^ 20) + (x ^ 3) + 1;
  //   const auto Z = ZechField(F2, f);
  //   std::cout << "Z = " << Z << std::endl;
  //   std::cout << "Z has cardinality " << Z.cardinality() << std::endl;
  // });

  // Untested
  timeit("Field of cardinality 2^32", []() {
    const auto F2 = SmallPrimeField(2);
    const auto x = Polynomial(F2, 'x');
    const auto f = (x ^ 16) + (x ^ 12) + (x ^ 3) + x + 1;
    const auto F2_16 = ZechField(F2, f);
    std::cout << F2_16 << " has cardinality " << F2_16.cardinality()
              << std::endl;
    const auto z = Polynomial(F2_16, 'z');
    const auto g = (z ^ 2) + F2_16.generator() * z + 1;
    const auto F2_32 = ZechPolyField(F2_16, g);
    std::cout << F2_32 << " has cardinality " << F2_32.cardinality()
              << std::endl;
    const auto zg = F2_32.generator();
  });

  timeit("Factoring", []() {
    for (integer_t p = 2; p < 13; next_prime(p)) {
      for (uint64_t k = 1; k < 100; ++k) {
        std::cout << "Factoring " << p << "^" << k << " - 1" << std::endl;
        print_factorization(factor_pk_minus_one(p, k));
      }
    }
  });
}
