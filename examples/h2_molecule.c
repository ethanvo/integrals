#include <stdio.h>
#include <stdlib.h>
#include "../include/integrals.h"

/*
 * Example: Computing integrals for H2 molecule
 * 
 * This example calculates the one-electron integrals (overlap, kinetic, nuclear)
 * and two-electron repulsion integrals for a simple H2 molecule using STO-3G basis.
 */

int main() {
    printf("======================================\n");
    printf("Quantum Chemistry Integral Engine\n");
    printf("Example: H2 Molecule (STO-3G basis)\n");
    printf("======================================\n\n");
    
    /* H2 molecule: Two H atoms separated by 1.4 bohr along x-axis */
    double h1_x = -0.7, h1_y = 0.0, h1_z = 0.0;  /* First hydrogen */
    double h2_x = 0.7, h2_y = 0.0, h2_z = 0.0;   /* Second hydrogen */
    
    /* STO-3G basis set for hydrogen 1s orbital */
    /* Exponents and coefficients from Pople's original paper */
    double sto3g_exponents[] = {3.42525091, 0.62391373, 0.16885540};
    double sto3g_coeffs[] = {0.15432897, 0.53532814, 0.44463454};
    
    printf("Setting up basis functions...\n");
    printf("H1 center: (%.2f, %.2f, %.2f)\n", h1_x, h1_y, h1_z);
    printf("H2 center: (%.2f, %.2f, %.2f)\n", h2_x, h2_y, h2_z);
    printf("Basis: STO-3G (3 primitives per atom)\n\n");
    
    /* Create contracted Gaussians for each hydrogen atom */
    ContractedGaussian *h1_basis = create_contracted_gaussian(3, h1_x, h1_y, h1_z);
    ContractedGaussian *h2_basis = create_contracted_gaussian(3, h2_x, h2_y, h2_z);
    
    /* Set up the STO-3G basis parameters */
    for (int i = 0; i < 3; i++) {
        h1_basis->primitives[i].exponent = sto3g_exponents[i];
        h1_basis->primitives[i].coeff = sto3g_coeffs[i];
        h2_basis->primitives[i].exponent = sto3g_exponents[i];
        h2_basis->primitives[i].coeff = sto3g_coeffs[i];
    }
    
    /* Compute overlap integrals */
    printf("Computing overlap integrals (S matrix)...\n");
    double s11 = overlap_contracted(h1_basis, h1_basis);
    double s12 = overlap_contracted(h1_basis, h2_basis);
    double s22 = overlap_contracted(h2_basis, h2_basis);
    
    printf("  S[1,1] = %.10f (should be ~1.0 for normalized basis)\n", s11);
    printf("  S[1,2] = %.10f\n", s12);
    printf("  S[2,2] = %.10f (should be ~1.0 for normalized basis)\n", s22);
    
    /* Compute kinetic energy integrals */
    printf("\nComputing kinetic energy integrals (T matrix)...\n");
    double t11 = kinetic_contracted(h1_basis, h1_basis);
    double t12 = kinetic_contracted(h1_basis, h2_basis);
    double t22 = kinetic_contracted(h2_basis, h2_basis);
    
    printf("  T[1,1] = %.10f\n", t11);
    printf("  T[1,2] = %.10f\n", t12);
    printf("  T[2,2] = %.10f\n", t22);
    
    /* Compute nuclear attraction integrals */
    printf("\nComputing nuclear attraction integrals (V matrix)...\n");
    printf("  V[1,1] with nucleus 1:\n");
    double v11_n1 = nuclear_attraction_contracted(h1_basis, h1_basis, h1_x, h1_y, h1_z);
    printf("    %.10f\n", v11_n1);
    
    printf("  V[1,1] with nucleus 2:\n");
    double v11_n2 = nuclear_attraction_contracted(h1_basis, h1_basis, h2_x, h2_y, h2_z);
    printf("    %.10f\n", v11_n2);
    
    printf("  V[1,2] with nucleus 1:\n");
    double v12_n1 = nuclear_attraction_contracted(h1_basis, h2_basis, h1_x, h1_y, h1_z);
    printf("    %.10f\n", v12_n1);
    
    printf("  V[1,2] with nucleus 2:\n");
    double v12_n2 = nuclear_attraction_contracted(h1_basis, h2_basis, h2_x, h2_y, h2_z);
    printf("    %.10f\n", v12_n2);
    
    printf("  V[2,2] with nucleus 1:\n");
    double v22_n1 = nuclear_attraction_contracted(h2_basis, h2_basis, h1_x, h1_y, h1_z);
    printf("    %.10f\n", v22_n1);
    
    printf("  V[2,2] with nucleus 2:\n");
    double v22_n2 = nuclear_attraction_contracted(h2_basis, h2_basis, h2_x, h2_y, h2_z);
    printf("    %.10f\n", v22_n2);
    
    /* Compute some two-electron repulsion integrals */
    printf("\nComputing two-electron repulsion integrals (ERIs)...\n");
    printf("  (11|11) = ");
    double eri_1111 = electron_repulsion_contracted(h1_basis, h1_basis, h1_basis, h1_basis);
    printf("%.10f\n", eri_1111);
    
    printf("  (11|12) = ");
    double eri_1112 = electron_repulsion_contracted(h1_basis, h1_basis, h1_basis, h2_basis);
    printf("%.10f\n", eri_1112);
    
    printf("  (11|22) = ");
    double eri_1122 = electron_repulsion_contracted(h1_basis, h1_basis, h2_basis, h2_basis);
    printf("%.10f\n", eri_1122);
    
    printf("  (12|12) = ");
    double eri_1212 = electron_repulsion_contracted(h1_basis, h2_basis, h1_basis, h2_basis);
    printf("%.10f\n", eri_1212);
    
    /* Core Hamiltonian matrix elements */
    printf("\nCore Hamiltonian (H = T + V):\n");
    double h11 = t11 + v11_n1 + v11_n2;
    double h12 = t12 + v12_n1 + v12_n2;
    double h22 = t22 + v22_n1 + v22_n2;
    printf("  H[1,1] = %.10f\n", h11);
    printf("  H[1,2] = %.10f\n", h12);
    printf("  H[2,2] = %.10f\n", h22);
    
    /* Cleanup */
    free_contracted_gaussian(h1_basis);
    free_contracted_gaussian(h2_basis);
    
    printf("\n======================================\n");
    printf("Example completed successfully!\n");
    printf("======================================\n");
    
    return 0;
}
