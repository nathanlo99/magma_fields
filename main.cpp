
#include "field.hpp"
#include "gmp.hpp"
#include "large_prime_field.hpp"
#include "lattice.hpp"
#include "matrix.hpp"
#include "medium_prime_field.hpp"
#include "polynomial.hpp"
#include "polynomial_factorization.hpp"
#include "prime_factorization.hpp"
#include "prime_poly.hpp"
#include "random.hpp"
#include "small_prime_field.hpp"
#include "timing.hpp"
#include "vector.hpp"
#include "zech_field.hpp"
#include "zech_poly_field.hpp"

#include <gmpxx.h>
#include <iostream>

inline void init_all() {
  init_gmp_random_state();
  load_cached_factorizations();
  load_cached_irreducible_polynomials();
}

int main(int argc, char *argv[]) {
  init_all();

  notimeit("Prime field of size near 2^16", []() {
    const auto F = SmallPrimeField(65521);
    const auto a = F(20000), b = F(30000);
    const auto product = a * b;
    const auto sum = a + b;
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << product << std::endl;
    std::cout << a << " + " << b << " = " << sum << std::endl;
  });

  notimeit("Prime field of medium size", []() {
    const auto F = MediumPrimeField(100003);
    const auto a = F(20000), b = F(30000);
    const auto product = a * b;
    const auto sum = a + b;
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << product << std::endl;
    std::cout << a << " + " << b << " = " << sum << std::endl;
  });

  notimeit("Prime field of size over 2^32", []() {
    const auto F = LargePrimeField(4294967311);
    const auto a = F(20000), b = F(30000);
    const auto product = a * b;
    const auto sum = a + b;
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << product << std::endl;
    std::cout << a << " + " << b << " = " << sum << std::endl;
  });

  notimeit("Field of cardinality 3^3", []() {
    const auto F3 = SmallPrimeField(3);
    const auto f = Polynomial(F3, "x", {1, 2, 0, 1});
    const auto Z = ZechField(F3, f);
    std::cout << "Z = " << Z << std::endl;
    const auto x = Z.primitive_element();
    for (int i = 0; i < Z.cardinality(); ++i) {
      std::cout << "x^" << i << " = " << (x ^ i) << std::endl;
    }
  });

  notimeit("Field of cardinality 2^5", []() {
    const auto F2 = SmallPrimeField(2);
    const auto x = Polynomial(F2, "x");
    const auto f = (x ^ 5) + (x ^ 2) + 1;
    const auto GF = PrimePolyField(F2, f);
    std::cout << "GF = " << GF << std::endl;
    std::cout << "GF has cardinality " << GF.cardinality() << std::endl;

    const auto x_GF = GF.element(x);
    for (int i = 0; i < GF.cardinality(); ++i) {
      std::cout << "x^" << i << " = " << (x_GF ^ i) << std::endl;
    }
  });

  notimeit("Field of cardinality 2^32", []() {
    const auto F2 = SmallPrimeField(2);
    const auto x = Polynomial(F2, "x");
    const auto f = (x ^ 16) + (x ^ 12) + (x ^ 3) + x + 1;
    const auto F2_16 = ZechField(F2, f);
    std::cout << F2_16 << " has cardinality " << F2_16.cardinality()
              << std::endl;
    const auto z = Polynomial(F2_16, "z");
    const auto g = (z ^ 2) + F2_16.primitive_element() * z + 1;
    const auto F2_32 = ZechPolyField(F2_16, g);
    std::cout << F2_32 << " has cardinality " << F2_32.cardinality()
              << std::endl;
    const auto zg = F2_32.primitive_element();
  });

  // NOTE: Currently buggy due to template and polymorphism hell
  notimeit("Finite field construction", []() {
    FiniteField(5, 4);   // Zech
    FiniteField(3, 24);  // Two-step: ZechPoly over a Zech
    FiniteField(2, 8);   // ZechField
    FiniteField(2, 103); // PrimePolyField
    FiniteField(2, 120); // ZechPoly over a ZechField of degree 20
    FiniteField(2, 1);   // Repeating the prime field
    FiniteField(2, 53);  // PrimePolyField
    FiniteField(2, 530); // Example from source
    FiniteField(2, 128); //
    print_lattices(std::cout);
  });

  timeit("Polynomial gcd's", []() {
    const auto F = SmallPrimeField(127);
    const auto x = Polynomial(F, "x");
    const auto f = (x ^ 2) + 7 * x + 6;
    std::cout << f << std::endl;
    const auto g = (x ^ 2) - 5 * x - 6;
    std::cout << g << std::endl;
    const auto gcd = polynomial_gcd(f, g);
    std::cout << "The gcd of " << f << " and " << g << " is " << gcd
              << std::endl;
  });

  timeit("Odd characteristic polynomial factorization", []() {
    const auto P = SmallPrimeField(3);
    const auto x = Polynomial(P, "x");
    const auto F = PrimePolyField(P, (x ^ 4) + x + 2);

    const auto w = Polynomial(F, "w");
    const auto f = (w ^ 2) + w + 2;

    std::cout << equal_degree_factorization(f, 1) << std::endl;
  });

  timeit("p = 2 polynomial factorization", []() {
    const auto P = SmallPrimeField(2);
    const auto x = Polynomial(P, "x");
    const auto F = PrimePolyField(P, (x ^ 4) + (x ^ 3) + 1);

    const auto w = Polynomial(F, "w");
    const auto f = (w ^ 2) + w + 1;

    std::cout << equal_degree_factorization(f, 1) << std::endl;

    const auto r = find_root(f);
    std::cout << "Root of " << f << ": " << r << std::endl;
  });

  notimeit("Root-finding fuzz test", []() {
    const std::array<integer_t, 2> primes = {2, 3};
    for (const integer_t &prime : primes) {
      const auto P = SmallPrimeField(prime);
      for (int base_degree = 3; base_degree <= 5; ++base_degree) {
        for (int extension_degree = 2 * base_degree; extension_degree <= 20;
             extension_degree += base_degree) {
          // Create an extension field of degree [extension_degree] and try and
          // find roots polynomials of irreducible polynomials in P[x]
          const auto f = get_irreducible_polynomial(P, "x", extension_degree);
          const auto F = PrimePolyField(P, f);
          log() << "Creating the field extension over " << P
                << " defined by f = " << f << std::endl;

          for (int i = 0; i < 4; ++i) {
            const auto h = get_irreducible_polynomial(P, "w", base_degree);

            // Lift h to F and factor it in F
            std::vector<integer_t> coeffs(h.degree() + 1, 0);
            for (size_t i : h.support)
              coeffs[i] = h.coeffs[i].value;
            const auto lifted_h = Polynomial(F, "w", coeffs);
            const auto root = find_root(lifted_h);
            std::cout << "Root of " << h << " is " << root << std::endl;
          }
        }
      }
    }
  });

  timeit("General polynomial factorization", []() {
    const auto P = SmallPrimeField(3);
    const auto x = Polynomial(P, "x");
    const auto f = (((x ^ 11) + 2 * (x ^ 9) + 2 * (x ^ 8) + (x ^ 6) + (x ^ 5) +
                     2 * (x ^ 3) + 2 * (x ^ 2) + 1) ^
                    20) +
                   1;
    // f is a 20-th power, plus 1, and has degree 220
    // On my computer, it factors in around 100ms
    std::cout << f << std::endl;
    std::cout << factor_polynomial(f) << std::endl;
  });

  timeit("Polynomial factorization fuzz-test", []() {
    // 'If it works for primes under 50, it works for all primes'
    const size_t max_prime = 50, max_degree = 60, num_iterations = 10,
                 degree_multiplier = max_degree / num_iterations;
    for (integer_t p = 2; p < max_prime; gmp::next_prime(p)) {
      const std::string lead_string = "Testing prime " + p.get_str() + "... ";
      std::cout << lead_string << std::flush;
      const auto P = SmallPrimeField(p);
      for (size_t i = 0; i < num_iterations; ++i) {
        std::cout << "\r" << lead_string << i << "/" << num_iterations
                  << std::flush;
        const auto g =
            random_polynomial<true, false>(P, "x", degree_multiplier * i + 1);
        const auto result = factor_polynomial(g);
      }
      std::cout << "\r" << lead_string << num_iterations << "/"
                << num_iterations << std::endl;
    }
  });

  timeit("Matrix row-reduction", []() {
    const auto P = SmallPrimeField(5);

    {
      std::vector<std::vector<integer_t>> coeffs = {
          {1, 2, -1, -4}, {2, 3, -1, -11}, {-2, 0, -3, 22}};

      auto M = Matrix(P, 3, 4, coeffs);
      std::cout << M << std::endl;
      const size_t rank = M.row_reduce();
      std::cout << M << std::endl;
      std::cout << "Rank " << rank << std::endl;
    }

    {
      std::vector<std::vector<integer_t>> coeffs = {
          {1, 2, -1}, {2, 3, -1}, {-2, 0, -3}};
      const auto M2 = Matrix(P, 3, 3, coeffs);
      const auto M_inv = M2.inverse();
      std::cout << M2 << std::endl;
      std::cout << M_inv << std::endl;
    }

    {
      std::vector<std::vector<integer_t>> coeffs = {{1, 2}, {3, 7}};
      const auto M = Matrix(P, 2, 2, coeffs);
      const auto M_inv = M.inverse();
      std::cout << M << std::endl;
      std::cout << M_inv << std::endl;
    }

    {
      std::vector<std::vector<integer_t>> A_coeffs = {{1, 2}, {3, 7}};
      const auto A = Matrix(P, 2, 2, A_coeffs);
      std::vector<integer_t> b_coeffs = {5, 10};
      const auto b = Vector(P, 2, b_coeffs);
      const auto x = A.solve(b);
      std::cout << "A = \n" << A << std::endl;
      std::cout << "b = " << b << std::endl;
      std::cout << "x = " << x << std::endl;
    }
  });

  timeit("To and from vectors", []() {
    const auto P = SmallPrimeField(5);
    const auto F = PrimePolyField(P, "x", 8ULL);
    std::cout << F << std::endl;
    for (int i = 0; i < 10'000; ++i) {
      const auto poly1 = F.random_element();
      const auto vec1 = F.to_vector(poly1);
      const auto poly2 = F.from_vector(vec1);
      const auto vec2 = F.to_vector(poly2);
      assert(poly1 == poly2);
      assert(vec1 == vec2);
    }
  });

  timeit("Embedding workspace", []() {
    // 0. Setup
    const auto P = SmallPrimeField(2);
    const auto x = Polynomial(P, "x");
    const auto E = ZechField(P, (x ^ 2) + x + 1);
    const auto w = Polynomial(P, "w");
    const auto F = ZechField(P, (w ^ 4) + w + 1);
    std::cout << P << std::endl;
    std::cout << E << std::endl;
    std::cout << F << std::endl;

    const auto embedding = Embed(E, F);
    const auto expected_phi = Matrix(P, 4, 2, {{1, 1}, {0, 1}, {0, 1}, {0, 0}});
    const auto expected_psi = Matrix(
        P, 4, 4, {{1, 1, 0, 0}, {0, 1, 1, 1}, {0, 1, 0, 1}, {0, 0, 0, 1}});
    const auto expected_alpha_FE = F.generating_element();
    const auto expected_f_FE_string = "w^2 + w + (x + 1)";

    assert(embedding.phi == expected_phi);
    assert(embedding.psi == expected_psi);
    assert(embedding.alpha_FE == expected_alpha_FE);
    assert(embedding.f_FE.to_string() == expected_f_FE_string);

    for (int i = 0; i < 100; ++i) {
      const auto elem = F.random_element();
      const auto vec_E = embedding.to_E_vector(elem);
      const auto elem2 = embedding.from_E_vector(vec_E);
      const auto vec_E2 = embedding.to_E_vector(elem2);
      assert(vec_E == vec_E2);
      assert(elem == elem2);
    }
  });

  timeit("More embedding tests", []() {
    const auto P = SmallPrimeField(5);
    const auto E = ZechField(P, "x", 3);
    const auto F = PrimePolyField(P, "w", 12);
    std::cout << P << std::endl;
    std::cout << E << std::endl;
    std::cout << F << std::endl;

    const auto embedding = Embed(E, F);
  });

  timeit("Square roots in finite field extensions", []() {
    const auto P = SmallPrimeField(11);
    const auto E = ZechField(P, "x", 3);
    const auto F = PrimePolyField(P, "w", 12);
    std::cout << P << std::endl;
    std::cout << E << std::endl;
    std::cout << F << std::endl;

    for (int i = 0; i < 10; ++i) {
      const auto elem = P.random_element();
      const auto &[has_root, root] = nth_root(P, elem, 2);
      if (has_root)
        std::cout << "A square root of " << elem << " in P is " << root
                  << std::endl;
    }

    for (int i = 0; i < 10; ++i) {
      const auto elem = E.random_element();
      const auto &[has_root, root] = nth_root(E, elem, 2);
      if (has_root)
        std::cout << "A square root of " << elem << " in E is " << root
                  << std::endl;
    }

    for (int i = 0; i < 10; ++i) {
      const auto elem = F.random_element();
      const auto &[has_root, root] = nth_root(F, elem, 2);
      if (has_root)
        std::cout << "A square root of " << elem << " in F is " << root
                  << std::endl;

      const auto &[has_cube_root, cube_root] = nth_root(F, elem, 3);
      if (has_cube_root)
        std::cout << "A cube root of " << elem << " in F is " << cube_root
                  << std::endl;
    }
  });
}
