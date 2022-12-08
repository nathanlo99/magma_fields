
#include "small_prime_field.hpp"
#include <iostream>

#include <gmpxx.h>

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
