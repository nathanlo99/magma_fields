
#include "lattice.hpp"
#include "field.hpp"
#include "gmp.hpp"
#include "large_prime_field.hpp"

#include <utility>

std::map<integer_t, Lattice<SmallPrimeField>> small_prime_lattices;
std::map<integer_t, Lattice<MediumPrimeField>> medium_prime_lattices;
std::map<integer_t, Lattice<LargePrimeField>> large_prime_lattices;
