#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gaussian_overlap.h"

static int nearly_equal(double a, double b, double tol) {
    return fabs(a - b) < tol;
}

int main(void) {
    const double A[3] = {0.0, 0.0, 0.0};
    const double B[3] = {0.0, 0.0, 0.0};
    const double alpha = 0.5;
    const double beta = 0.5;
    const double pi = acos(-1.0);
    const double expected = pow(pi / (alpha + beta), 1.5);
    const double computed = gaussian_overlap_s(alpha, beta, A, B);

    if (!nearly_equal(computed, expected, 1e-12)) {
        fprintf(stderr, "Expected %.15f but got %.15f\n", expected, computed);
        return EXIT_FAILURE;
    }

    const double B_shifted[3] = {1.0, 0.0, 0.0};
    const double computed_shifted = gaussian_overlap_s(alpha, beta, A, B_shifted);
    if (computed_shifted >= computed) {
        fprintf(stderr, "Expected shifted overlap to be smaller.\n");
        return EXIT_FAILURE;
    }

    printf("All gaussian overlap tests passed.\n");
    return EXIT_SUCCESS;
}
