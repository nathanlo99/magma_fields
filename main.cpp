
#include <iostream>

#include <gmpxx.h>

mpz_class fact(const mpz_class n) {
  return n <= 1 ? mpz_class(1) : n * factorial(n - 1);
}

int main(int argc, char *argv[]) {
  const mpz_class n = 100;
  std::cout << "factorial(" << n << ") = " << factorial(n) << std::endl;
  std::cout << "Hello, world!" << std::endl;
}
