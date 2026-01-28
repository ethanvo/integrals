#include <math.h>
#include "gaussian_overlap.h"

double gaussian_overlap_s(double alpha, double beta, const double A[3], const double B[3]) {
    const double pi = acos(-1.0);
    const double gamma = alpha + beta;
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    const double rab2 = dx * dx + dy * dy + dz * dz;
    return pow(pi / gamma, 1.5) * exp(-(alpha * beta / gamma) * rab2);
}
