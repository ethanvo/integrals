#include <math.h>
#include "gaussian_overlap.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Obara-Saika s-type Gaussian overlap integral between centers A and B.
 * alpha and beta are Gaussian exponents and must be positive.
 */
double gaussian_overlap_s(double alpha, double beta, const double A[3], const double B[3]) {
    if (alpha <= 0.0 || beta <= 0.0) {
        return 0.0;
    }
    const double gamma = alpha + beta;
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double rab2 = dx * dx + dy * dy + dz * dz;
    return pow(M_PI / gamma, 1.5) * exp(-(alpha * beta / gamma) * rab2);
}
