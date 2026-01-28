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
    printf("=== Kinetic Energy Integral Tests ===\n\n");
    
    /* Test 1: Kinetic energy of two identical s-type Gaussians at same position */
    printf("Test 1: Identical Gaussians at same position\n");
    double t1 = kinetic_1s(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    double s1 = overlap_1s(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
    /* T = α₁α₂/(α₁+α₂) * 3 * S when R=0 */
    double expected1 = 0.5 * 3.0 * s1;
    assert_near("T(1s,1s) same center", t1, expected1, TOLERANCE);
    
    /* Test 2: Kinetic energy with different exponents */
    printf("\nTest 2: Different exponents\n");
    double alpha1 = 0.5, alpha2 = 1.5;
    double t2 = kinetic_1s(alpha1, 0.0, 0.0, 0.0, alpha2, 0.0, 0.0, 0.0);
    double s2 = overlap_1s(alpha1, 0.0, 0.0, 0.0, alpha2, 0.0, 0.0, 0.0);
    double term = alpha1 * alpha2 / (alpha1 + alpha2);
    double expected2 = term * 3.0 * s2;
    assert_near("T(1s,1s) different alpha", t2, expected2, TOLERANCE);
    
    /* Test 3: Kinetic energy of separated Gaussians */
    printf("\nTest 3: Gaussians separated by 1 bohr\n");
    double t3 = kinetic_1s(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    double s3 = overlap_1s(1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
    double rab2 = 1.0;
    term = 0.5;
    double expected3 = term * (3.0 - 2.0 * term * rab2) * s3;
    assert_near("T(1s,1s) 1 bohr apart", t3, expected3, TOLERANCE);
    
    /* Test 4: Contracted Gaussian kinetic energy */
    printf("\nTest 4: Contracted Gaussians (STO-3G style)\n");
    ContractedGaussian *cg1 = create_contracted_gaussian(2, 0.0, 0.0, 0.0);
    ContractedGaussian *cg2 = create_contracted_gaussian(2, 0.0, 0.0, 0.0);
    
    /* Simple test with 2 primitives */
    cg1->primitives[0].exponent = 1.0;
    cg1->primitives[0].coeff = 0.6;
    cg1->primitives[1].exponent = 0.3;
    cg1->primitives[1].coeff = 0.4;
    
    cg2->primitives[0].exponent = 1.0;
    cg2->primitives[0].coeff = 0.6;
    cg2->primitives[1].exponent = 0.3;
    cg2->primitives[1].coeff = 0.4;
    
    double t4 = kinetic_contracted(cg1, cg2);
    /* Just verify it's positive and reasonable */
    if (t4 > 0.0 && t4 < 10.0) {
        test_passed++;
        test_count++;
        printf("[PASS] T(contracted,contracted) > 0: %.10f\n", t4);
    } else {
        test_count++;
        printf("[FAIL] T(contracted,contracted) unreasonable: %.10f\n", t4);
    }
    
    free_contracted_gaussian(cg1);
    free_contracted_gaussian(cg2);
    
    /* Summary */
    printf("\n=== Test Summary ===\n");
    printf("Tests passed: %d/%d\n", test_passed, test_count);
    
    return (test_passed == test_count) ? 0 : 1;
}
