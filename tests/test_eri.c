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
    printf("=== Two-Electron Repulsion Integral Tests ===\n\n");
    
    /* Test 1: ERI for four identical Gaussians at same position */
    printf("Test 1: Four identical Gaussians at same position\n");
    double eri1 = electron_repulsion_1s(1.0, 0.0, 0.0, 0.0,
                                        1.0, 0.0, 0.0, 0.0,
                                        1.0, 0.0, 0.0, 0.0,
                                        1.0, 0.0, 0.0, 0.0);
    /* Should be positive */
    if (eri1 > 0.0) {
        test_passed++;
        test_count++;
        printf("[PASS] ERI(same position) > 0: %.10f\n", eri1);
    } else {
        test_count++;
        printf("[FAIL] ERI should be positive: %.10f\n", eri1);
    }
    
    /* Test 2: ERI should decrease with separation */
    printf("\nTest 2: ERI decreases with separation\n");
    double eri2 = electron_repulsion_1s(1.0, 0.0, 0.0, 0.0,
                                        1.0, 0.0, 0.0, 0.0,
                                        1.0, 2.0, 0.0, 0.0,
                                        1.0, 2.0, 0.0, 0.0);
    if (eri1 > eri2 && eri2 > 0.0) {
        test_passed++;
        test_count++;
        printf("[PASS] ERI(close) > ERI(far): %.10f > %.10f\n", eri1, eri2);
    } else {
        test_count++;
        printf("[FAIL] ERI should decrease with separation\n");
    }
    
    /* Test 3: ERI symmetry - (ab|cd) = (ba|cd) */
    printf("\nTest 3: ERI symmetry (ab|cd) = (ba|cd)\n");
    double eri3a = electron_repulsion_1s(1.0, 0.0, 0.0, 0.0,
                                         1.5, 1.0, 0.0, 0.0,
                                         0.8, 0.5, 0.5, 0.0,
                                         1.2, 1.5, 0.5, 0.0);
    double eri3b = electron_repulsion_1s(1.5, 1.0, 0.0, 0.0,
                                         1.0, 0.0, 0.0, 0.0,
                                         0.8, 0.5, 0.5, 0.0,
                                         1.2, 1.5, 0.5, 0.0);
    assert_near("ERI symmetry (ab|cd)=(ba|cd)", eri3a, eri3b, TOLERANCE);
    
    /* Test 4: ERI symmetry - (ab|cd) = (ab|dc) */
    printf("\nTest 4: ERI symmetry (ab|cd) = (ab|dc)\n");
    double eri4a = electron_repulsion_1s(1.0, 0.0, 0.0, 0.0,
                                         1.5, 1.0, 0.0, 0.0,
                                         0.8, 0.5, 0.5, 0.0,
                                         1.2, 1.5, 0.5, 0.0);
    double eri4b = electron_repulsion_1s(1.0, 0.0, 0.0, 0.0,
                                         1.5, 1.0, 0.0, 0.0,
                                         1.2, 1.5, 0.5, 0.0,
                                         0.8, 0.5, 0.5, 0.0);
    assert_near("ERI symmetry (ab|cd)=(ab|dc)", eri4a, eri4b, TOLERANCE);
    
    /* Test 5: Contracted Gaussian ERI */
    printf("\nTest 5: Contracted Gaussians\n");
    ContractedGaussian *cg1 = create_contracted_gaussian(1, 0.0, 0.0, 0.0);
    ContractedGaussian *cg2 = create_contracted_gaussian(1, 0.0, 0.0, 0.0);
    ContractedGaussian *cg3 = create_contracted_gaussian(1, 1.0, 0.0, 0.0);
    ContractedGaussian *cg4 = create_contracted_gaussian(1, 1.0, 0.0, 0.0);
    
    cg1->primitives[0].exponent = 1.0;
    cg1->primitives[0].coeff = 1.0;
    cg2->primitives[0].exponent = 1.0;
    cg2->primitives[0].coeff = 1.0;
    cg3->primitives[0].exponent = 1.0;
    cg3->primitives[0].coeff = 1.0;
    cg3->primitives[0].x = 1.0;
    cg4->primitives[0].exponent = 1.0;
    cg4->primitives[0].coeff = 1.0;
    cg4->primitives[0].x = 1.0;
    
    double eri5 = electron_repulsion_contracted(cg1, cg2, cg3, cg4);
    if (eri5 > 0.0) {
        test_passed++;
        test_count++;
        printf("[PASS] ERI(contracted) > 0: %.10f\n", eri5);
    } else {
        test_count++;
        printf("[FAIL] ERI(contracted) should be positive: %.10f\n", eri5);
    }
    
    free_contracted_gaussian(cg1);
    free_contracted_gaussian(cg2);
    free_contracted_gaussian(cg3);
    free_contracted_gaussian(cg4);
    
    /* Summary */
    printf("\n=== Test Summary ===\n");
    printf("Tests passed: %d/%d\n", test_passed, test_count);
    
    return (test_passed == test_count) ? 0 : 1;
}
