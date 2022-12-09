
#include "large_prime_field.hpp"
#include "medium_prime_field.hpp"
#include "small_prime_field.hpp"
#include "util.hpp"
#include "zech_field.hpp"
#include <iostream>

#include <gmpxx.h>

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

  {
    const auto F3 = SmallPrimeField(3);
    const auto f = Polynomial(F3, 'x', {1, 2, 0, 1});
    const auto Z = ZechField(F3, 3, f);
    std::cout << "Z = " << Z << std::endl;
    const auto x = Z.generator();
    for (int i = 0; i < Z.cardinality(); ++i) {
      std::cout << "x^" << i << " = " << (x ^ i) << std::endl;
    }
  }
}
