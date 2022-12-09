
#include "gmp.hpp"
#include "large_prime_field.hpp"
#include "medium_prime_field.hpp"
#include "small_prime_field.hpp"
#include "timing.hpp"
#include "zech_field.hpp"

#include <gmpxx.h>
#include <iostream>

int main(int argc, char *argv[]) {
  {
    const auto F = SmallPrimeField(65521);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  }

  {
    const auto F = MediumPrimeField(100003);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  }

  {
    const auto F = LargePrimeField(4294967311);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  }

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

  timeit("Field of cardinality 2^20", []() {
    const auto F2 = SmallPrimeField(2);
    const auto x = Polynomial(F2, 'x', {0, 1});
    const auto f = (x ^ 20) + (x ^ 3) + 1;
    const auto Z = ZechField(F2, f);
    std::cout << "Z = " << Z << std::endl;
    std::cout << "Z has cardinality " << Z.cardinality() << std::endl;
  });
}
