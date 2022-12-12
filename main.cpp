
#include "gmp.hpp"
#include "large_prime_field.hpp"
#include "lattice.hpp"
#include "medium_prime_field.hpp"
#include "prime_factorization.hpp"
#include "prime_poly.hpp"
#include "random.hpp"
#include "small_prime_field.hpp"
#include "timing.hpp"
#include "zech_field.hpp"
#include "zech_poly_field.hpp"

#include <gmpxx.h>
#include <iostream>

inline void init_all() {
  load_cached_factorizations();
  init_gmp_random_state();
}

int main(int argc, char *argv[]) {
  init_all();

  // timeit("Prime field of size near 2^16", []() {
  //   const auto F = SmallPrimeField(65521);
  //   const auto a = F(20000), b = F(30000);
  //   std::cout << "In " << F << ", " << std::endl;
  //   std::cout << a << " * " << b << " = " << a * b << std::endl;
  // });

  // timeit("Prime field of medium size", []() {
  //   const auto F = MediumPrimeField(100003);
  //   const auto a = F(20000), b = F(30000);
  //   std::cout << "In " << F << ", " << std::endl;
  //   std::cout << a << " * " << b << " = " << a * b << std::endl;
  // });

  // timeit("Prime field of size over 2^32", []() {
  //   const auto F = LargePrimeField(4294967311);
  //   const auto a = F(20000), b = F(30000);
  //   std::cout << "In " << F << ", " << std::endl;
  //   std::cout << a << " * " << b << " = " << a * b << std::endl;
  // });

  // timeit("Field of cardinality 3^3", []() {
  //   const auto F3 = SmallPrimeField(3);
  //   const auto f = Polynomial(F3, "x", {1, 2, 0, 1});
  //   const auto Z = ZechField(F3, f);
  //   std::cout << "Z = " << Z << std::endl;
  //   const auto x = Z.primitive_element();
  //   for (int i = 0; i < Z.cardinality(); ++i) {
  //     std::cout << "x^" << i << " = " << (x ^ i) << std::endl;
  //   }
  // });

  // timeit("Field of cardinality 2^5", []() {
  //   const auto F2 = SmallPrimeField(2);
  //   const auto x = Polynomial(F2, "x");
  //   const auto f = (x ^ 5) + (x ^ 2) + 1;
  //   const auto GF = PrimePolyField(F2, f);
  //   std::cout << "GF = " << GF << std::endl;
  //   std::cout << "GF has cardinality " << GF.cardinality() << std::endl;

  //   const auto x_GF = GF.element(x);
  //   for (int i = 0; i < GF.cardinality(); ++i) {
  //     std::cout << "x^" << i << " = " << (x_GF ^ i) << std::endl;
  //   }
  // });

  // Untested
  // timeit("Field of cardinality 2^32", []() {
  //   const auto F2 = SmallPrimeField(2);
  //   const auto x = Polynomial(F2, "x");
  //   const auto f = (x ^ 16) + (x ^ 12) + (x ^ 3) + x + 1;
  //   const auto F2_16 = ZechField(F2, f);
  //   std::cout << F2_16 << " has cardinality " << F2_16.cardinality()
  //             << std::endl;
  //   const auto z = Polynomial(F2_16, "z");
  //   const auto g = (z ^ 2) + F2_16.primitive_element() * z + 1;
  //   const auto F2_32 = ZechPolyField(F2_16, g);
  //   std::cout << F2_32 << " has cardinality " << F2_32.cardinality()
  //             << std::endl;
  //   const auto zg = F2_32.primitive_element();
  // });

  // timeit("Factoring",
  //        []() { print_factorization(factor_pk_minus_one(2, 127)); });

  // timeit("Order-finding", []() {
  //   const auto F = MediumPrimeField(1000000007);
  //   std::cout << "Primitive element in " << F << " is " <<
  //   F.primitive_element() << std::endl;
  // });

  timeit("Debug demo", []() {
    LatticeManager manager;
    manager.FiniteField(5, 4);   // Zech
    manager.FiniteField(3, 24);  // Two-step: ZechPoly over a Zech
    manager.FiniteField(2, 8);   // ZechField
    manager.FiniteField(2, 103); // PrimePolyField
    manager.FiniteField(2, 120); // ZechPoly over a ZechField of degree 20
    for (integer_t p = 2; p < 10000000000000_mpz; p *= 2, gmp::next_prime(p)) {
      manager.FiniteField(p, 10);
    }
    std::cout << manager << std::endl;
  });

  // timeit("Polynomial gcd's", []() {
  //   const auto F = SmallPrimeField(127);
  //   const auto x = Polynomial(F, "x");
  //   const auto f = (x ^ 2) + 7 * x + 6;
  //   std::cout << f << std::endl;
  //   const auto g = (x ^ 2) - 5 * x - 6;
  //   std::cout << g << std::endl;
  //   const auto gcd = polynomial_gcd(f, g);
  //   std::cout << "The gcd of " << f << " and " << g << " is " << gcd
  //             << std::endl;
  // });

  // timeit("Rabin irreducibility algorithm", []() {
  //   const integer_t p = 2;
  //   const auto F = SmallPrimeField(p);
  //   for (uint64_t degree = 1; degree <= 11; ++degree) {
  //     std::cout << "For degree " << degree << "... " << std::flush;
  //     uint64_t least_support = degree + 2;
  //     Polynomial best_polynomial = Polynomial(F, "x");
  //     uint64_t num_irreducible = 0;
  //     for (integer_t iter = 0; iter < gmp::pow(p, degree + 1); ++iter) {
  //       std::vector<integer_t> coeffs(degree + 1);
  //       integer_t coeff_idx = iter;
  //       for (size_t i = 0; i <= degree; ++i) {
  //         coeffs[i] = coeff_idx % p;
  //         coeff_idx /= p;
  //       }
  //       if (coeffs[degree] == 0 || coeffs[0] == 0)
  //         continue;
  //       const Polynomial f = Polynomial(F, "x", coeffs);
  //       if (!f.is_irreducible_rabin())
  //         continue;
  //       num_irreducible++;
  //       if (f.support.size() < least_support) {
  //         least_support = f.support.size();
  //         best_polynomial = f;
  //       }
  //     }
  //     std::cout << num_irreducible
  //               << " irreducible polynomials, with shortest = "
  //               << best_polynomial << std::endl;
  //   }
  // });

  // timeit("Generating random polynomials", []() {
  //   const auto F = SmallPrimeField(5);
  //   for (int i = 0; i < 10; ++i) {
  //     const Polynomial p = random_polynomial<false, true>(F, "x", 4);
  //     std::cout << p << std::endl;
  //   }
  // });
}
