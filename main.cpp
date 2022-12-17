
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
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  });

  notimeit("Prime field of medium size", []() {
    const auto F = MediumPrimeField(100003);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  });

  notimeit("Prime field of size over 2^32", []() {
    const auto F = LargePrimeField(4294967311);
    const auto a = F(20000), b = F(30000);
    std::cout << "In " << F << ", " << std::endl;
    std::cout << a << " * " << b << " = " << a * b << std::endl;
  });

  timeit("Field of cardinality 3^3", []() {
    const auto F3 = SmallPrimeField(3);
    const auto f = Polynomial(F3, "x", {1, 2, 0, 1});
    const auto Z = ZechField(F3, f);
    std::cout << "Z = " << Z << std::endl;
    const auto x = Z.primitive_element();
    for (int i = 0; i < Z.cardinality(); ++i) {
      std::cout << "x^" << i << " = " << (x ^ i) << std::endl;
    }
  });

  timeit("Field of cardinality 2^5", []() {
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

  notimeit("Debug demo", []() {
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

  notimeit("Polynomial gcd's", []() {
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

    print_polynomial_factorization(equal_degree_factorization(f, 1));
  });

  timeit("p = 2 polynomial factorization", []() {
    const auto P = SmallPrimeField(2);
    const auto x = Polynomial(P, "x");
    const auto F = PrimePolyField(P, (x ^ 4) + (x ^ 3) + 1);

    const auto w = Polynomial(F, "w");
    const auto f = (w ^ 2) + w + 1;

    print_polynomial_factorization(equal_degree_factorization(f, 1));

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

  notimeit("Matrix row-reduction", []() {
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

    const auto e = E.degree(), f = F.degree();
    assert(f % e == 0);
    const auto d = f / e;

    // 1. Computing prime-field generators
    const auto alpha_E = E.generating_element();
    const auto alpha_F = F.generating_element();
    const auto tau = find_root(E.f.lift(F));
    std::cout << "alpha_E = " << alpha_E << std::endl;
    std::cout << "alpha_F = " << alpha_F << std::endl;
    std::cout << "tau = " << tau << std::endl;

    if (F.g != F.g.var_poly())
      throw math_error("Currently only supports fields generated by w");

    // 2. Computing the matrix phi_{E -> F}
    auto M = Matrix(P, e, f);
    for (size_t i = 0; i < e; ++i) {
      const auto tau_i = F.to_polynomial((tau ^ i).value);
      std::cout << "tau^" << i << " = " << tau_i << std::endl;
      for (size_t j = 0; j < f; ++j) {
        M[i][j] = tau_i[j];
      }
    }
    std::cout << M << std::endl;

    // 3. Finding a generator for F/E
    // If G is contained in E and we already have a generator for F over
    // G, then use that
    // TODO: Implement the above more generally once we have embeddings
    // For now, we just restrict ourselves to fields F with G = P
    const auto alpha_FE = F.generating_element();

    // 3. Compute the vector-space isomorphism E^{(d)} -> F
    Matrix N(P, f, f);
    for (size_t i = 0; i < d; ++i) {
      for (size_t j = 0; j < e; ++j) {
        // Compute \phi_F^{-1}(alpha_FE^i * tau^j)
        const auto alpha_FE_to_i = alpha_FE ^ i;
        const auto tau_to_j = tau ^ j;
        const auto product = alpha_FE_to_i * tau_to_j;
        const auto row = F.to_vector(product);
        N[i * e + j] = row.data;
      }
    }
    std::cout << N << std::endl;
  });
}
