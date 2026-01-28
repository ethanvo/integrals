# Quantum Chemistry Integral Engine

A C library for computing quantum chemistry integrals, designed for educational purposes and practice.

## Features

This library implements the calculation of fundamental quantum chemistry integrals using Gaussian basis functions:

- **Overlap integrals (S)**: Measure the overlap between two basis functions
- **Kinetic energy integrals (T)**: Compute the kinetic energy operator matrix elements
- **Nuclear attraction integrals (V)**: Calculate the electron-nucleus attraction
- **Two-electron repulsion integrals (ERI)**: Compute the electron-electron repulsion integrals

The implementation supports:
- Primitive Gaussian functions
- Contracted Gaussian functions (e.g., STO-3G basis sets)
- s-type orbitals (spherical Gaussians)

## Building

Build the library and tests:

```bash
make          # Build the static library (libintegrals.a)
make tests    # Build test programs
make examples # Build example programs
make all      # Build everything
```

Run tests:

```bash
make run-tests
```

Clean build artifacts:

```bash
make clean
```

## Usage

### Basic Example

```c
#include "integrals.h"

int main() {
    // Compute overlap integral between two s-type Gaussians
    double s = overlap_1s(
        1.0, 0.0, 0.0, 0.0,  // alpha1, x1, y1, z1
        1.0, 1.0, 0.0, 0.0   // alpha2, x2, y2, z2
    );
    printf("Overlap: %f\n", s);
    
    return 0;
}
```

### Using Contracted Gaussians

```c
#include "integrals.h"

int main() {
    // Create a contracted Gaussian with 3 primitives (STO-3G style)
    ContractedGaussian *cg = create_contracted_gaussian(3, 0.0, 0.0, 0.0);
    
    // Set exponents and coefficients
    cg->primitives[0].exponent = 3.42525091;
    cg->primitives[0].coeff = 0.15432897;
    // ... set other primitives ...
    
    // Compute integrals
    double s = overlap_contracted(cg, cg);
    
    // Cleanup
    free_contracted_gaussian(cg);
    
    return 0;
}
```

### Complete H2 Molecule Example

See `examples/h2_molecule.c` for a complete example that computes all integrals for an H2 molecule using the STO-3G basis set.

```bash
make examples
./bin/h2_molecule
```

## Library Structure

```
integrals/
├── include/
│   └── integrals.h         # Header file with API
├── src/
│   ├── overlap.c           # Overlap integral implementations
│   ├── kinetic.c           # Kinetic energy integral implementations
│   ├── nuclear.c           # Nuclear attraction integrals and Boys function
│   ├── eri.c               # Two-electron repulsion integrals
│   └── utils.c             # Utility functions
├── tests/
│   ├── test_overlap.c      # Overlap integral tests
│   ├── test_kinetic.c      # Kinetic energy tests
│   ├── test_nuclear.c      # Nuclear attraction tests
│   └── test_eri.c          # ERI tests
└── examples/
    └── h2_molecule.c       # H2 molecule example
```

## Theory

### Gaussian Basis Functions

The library uses Gaussian-type orbitals (GTOs) of the form:

```
φ(r) = exp(-α|r-R|²)
```

where α is the Gaussian exponent and R is the center position.

### Contracted Gaussians

Contracted Gaussian functions are linear combinations of primitive Gaussians:

```
χ(r) = Σᵢ cᵢ exp(-αᵢ|r-R|²)
```

This is used in standard basis sets like STO-3G, where each basis function is a contraction of 3 primitive Gaussians.

### Integral Types

1. **Overlap Integral**: S = ⟨φₐ|φᵦ⟩
2. **Kinetic Energy**: T = ⟨φₐ|(-1/2)∇²|φᵦ⟩
3. **Nuclear Attraction**: V = ⟨φₐ|(-Z/|r-C|)|φᵦ⟩
4. **Electron Repulsion**: (ab|cd) = ⟨φₐφᵦ|(1/r₁₂)|φ꜀φ_d⟩

### Boys Function

The Boys function F_n(x) is used in nuclear attraction and electron repulsion integrals:

```
F_n(x) = ∫₀¹ t²ⁿ exp(-xt²) dt
```

## References

This implementation is based on standard quantum chemistry textbooks:
- Szabo & Ostlund, "Modern Quantum Chemistry"
- Helgaker, Jorgensen & Olsen, "Molecular Electronic-Structure Theory"

## License

This is an educational project for practice purposes.

## Author

Built as a practice project for learning quantum chemistry computational methods.