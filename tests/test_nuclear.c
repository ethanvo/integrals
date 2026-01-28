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
    printf("=== Nuclear Attraction Integral Tests ===\n\n");
    
    /* Test 1: Boys function at x=0 */
    printf("Test 1: Boys function F_0(0) = 1\n");
    double f0_0 = boys_function(0, 0.0);
    assert_near("F_0(0)", f0_0, 1.0, TOLERANCE);
    
    /* Test 2: Boys function F_1(0) = 1/3 */
    printf("\nTest 2: Boys function F_1(0) = 1/3\n");
    double f1_0 = boys_function(1, 0.0);
    assert_near("F_1(0)", f1_0, 1.0/3.0, TOLERANCE);
    
    /* Test 3: Nuclear attraction with nucleus at center of Gaussian */
    printf("\nTest 3: Nuclear attraction, nucleus at Gaussian center\n");
    double v1 = nuclear_attraction_1s(1.0, 0.0, 0.0, 0.0, 
                                      1.0, 0.0, 0.0, 0.0,
                                      0.0, 0.0, 0.0);
    /* When nucleus is at center, this should be negative */
    if (v1 < 0.0) {
        test_passed++;
        test_count++;
        printf("[PASS] V(nucleus at center) < 0: %.10f\n", v1);
    } else {
        test_count++;
        printf("[FAIL] V(nucleus at center) should be negative: %.10f\n", v1);
    }
    
    /* Test 4: Nuclear attraction should be more negative when closer */
    printf("\nTest 4: Nuclear attraction increases with proximity\n");
    double v2 = nuclear_attraction_1s(1.0, 0.0, 0.0, 0.0, 
                                      1.0, 0.0, 0.0, 0.0,
                                      2.0, 0.0, 0.0);
    if (v1 < v2 && v2 < 0.0) {
        test_passed++;
        test_count++;
        printf("[PASS] V(closer) < V(farther): %.10f < %.10f\n", v1, v2);
    } else {
        test_count++;
        printf("[FAIL] Nuclear attraction should decrease with distance\n");
    }
    
    /* Test 5: Contracted Gaussian nuclear attraction */
    printf("\nTest 5: Contracted Gaussians\n");
    ContractedGaussian *cg1 = create_contracted_gaussian(2, 0.0, 0.0, 0.0);
    ContractedGaussian *cg2 = create_contracted_gaussian(2, 1.0, 0.0, 0.0);
    
    cg1->primitives[0].exponent = 1.0;
    cg1->primitives[0].coeff = 0.6;
    cg1->primitives[1].exponent = 0.3;
    cg1->primitives[1].coeff = 0.4;
    
    cg2->primitives[0].exponent = 1.0;
    cg2->primitives[0].coeff = 0.6;
    cg2->primitives[1].exponent = 0.3;
    cg2->primitives[1].coeff = 0.4;
    
    /* Update primitive centers */
    for (int i = 0; i < 2; i++) {
        cg2->primitives[i].x = 1.0;
    }
    
    double v3 = nuclear_attraction_contracted(cg1, cg2, 0.5, 0.0, 0.0);
    /* Should be negative */
    if (v3 < 0.0) {
        test_passed++;
        test_count++;
        printf("[PASS] V(contracted) < 0: %.10f\n", v3);
    } else {
        test_count++;
        printf("[FAIL] V(contracted) should be negative: %.10f\n", v3);
    }
    
    free_contracted_gaussian(cg1);
    free_contracted_gaussian(cg2);
    
    /* Summary */
    printf("\n=== Test Summary ===\n");
    printf("Tests passed: %d/%d\n", test_passed, test_count);
    
    return (test_passed == test_count) ? 0 : 1;
}
