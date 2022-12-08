
#include "large_prime_field.hpp"
#include "medium_prime_field.hpp"
#include "small_prime_field.hpp"
#include "util.hpp"
#include <iostream>

#include <gmpxx.h>

using FiniteField =
    std::variant<SmallPrimeField, MediumPrimeField, LargePrimeField>;

mpz_class fact(const mpz_class n) {
  return n <= 1 ? 1_mpz : n * factorial(n - 1);
}

int main(int argc, char *argv[]) {
  const mpz_class n = 100_mpz;
  std::cout << "factorial(" << n << ") = " << factorial(n) << std::endl;

  const auto F = SmallPrimeField(65521);
  const auto two = F(20000), three = F(30000);
  std::cout << two << " * " << three << " = " << two * three << std::endl;
}
