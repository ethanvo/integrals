#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/integrals.h"

#define TOLERANCE 1e-6

int test_count = 0;
int test_passed = 0;

void assert_near(const char *test_name, double actual, double expected, double tol) {
    test_count++;
    if (fabs(actual - expected) < tol) {
        test_passed++;
        printf("[PASS] %s: %.10f ≈ %.10f\n", test_name, actual, expected);
    } else {
        printf("[FAIL] %s: %.10f != %.10f (diff: %.10e)\n", 
               test_name, actual, expected, fabs(actual - expected));
    }
}

int main() {
    printf("=== Overlap Integral Tests ===\n\n");
    
    /* Test 1: Overlap of two identical s-type Gaussians at same position */
    printf("Test 1: Identical Gaussians at same position\n");
    double s1 = overlap_1s(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    /* Expected: (π/2)^(3/2) ≈ 1.244 */
    double expected1 = pow(M_PI / 2.0, 1.5);
    assert_near("S(1s,1s) same center", s1, expected1, TOLERANCE);
    
    /* Test 2: Overlap of two s-type Gaussians separated by 1 bohr */
    printf("\nTest 2: Gaussians separated by 1 bohr\n");
    double s2 = overlap_1s(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    /* Expected: (π/2)^(3/2) * exp(-0.5) ≈ 0.754 */
    double expected2 = pow(M_PI / 2.0, 1.5) * exp(-0.5);
    assert_near("S(1s,1s) 1 bohr apart", s2, expected2, TOLERANCE);
    
    /* Test 3: Overlap with different exponents */
    printf("\nTest 3: Different exponents\n");
    double s3 = overlap_1s(0.5, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0);
    /* Expected: (π/2)^(3/2) ≈ 1.244 (same center) */
    double gamma = 0.5 + 1.5;
    double expected3 = pow(M_PI / gamma, 1.5);
    assert_near("S(1s,1s) different alpha", s3, expected3, TOLERANCE);
    
    /* Test 4: Contracted Gaussian overlap */
    printf("\nTest 4: Contracted Gaussians (STO-3G style)\n");
    ContractedGaussian *cg1 = create_contracted_gaussian(3, 0.0, 0.0, 0.0);
    ContractedGaussian *cg2 = create_contracted_gaussian(3, 0.0, 0.0, 0.0);
    
    /* STO-3G exponents and coefficients for H 1s */
    double sto3g_exponents[] = {3.42525091, 0.62391373, 0.16885540};
    double sto3g_coeffs[] = {0.15432897, 0.53532814, 0.44463454};
    
    for (int i = 0; i < 3; i++) {
        cg1->primitives[i].exponent = sto3g_exponents[i];
        cg1->primitives[i].coeff = sto3g_coeffs[i];
        cg2->primitives[i].exponent = sto3g_exponents[i];
        cg2->primitives[i].coeff = sto3g_coeffs[i];
    }
    
    double s4 = overlap_contracted(cg1, cg2);
    /* For contracted Gaussians, just check that it's positive and reasonable */
    /* The raw STO-3G coefficients are not normalized to unity */
    if (s4 > 0.0 && s4 < 20.0) {
        test_passed++;
        test_count++;
        printf("[PASS] S(STO-3G,STO-3G) > 0 and reasonable: %.10f\n", s4);
    } else {
        test_count++;
        printf("[FAIL] S(STO-3G,STO-3G) unreasonable: %.10f\n", s4);
    }
    
    free_contracted_gaussian(cg1);
    free_contracted_gaussian(cg2);
    
    /* Summary */
    printf("\n=== Test Summary ===\n");
    printf("Tests passed: %d/%d\n", test_passed, test_count);
    
    return (test_passed == test_count) ? 0 : 1;
}
