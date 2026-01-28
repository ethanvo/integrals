#ifndef GAUSSIAN_OVERLAP_H
#define GAUSSIAN_OVERLAP_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Computes the s-type Gaussian overlap integral between centers A and B:
 * S = (pi / (alpha + beta))^(3/2) * exp(-(alpha * beta / (alpha + beta)) * |A - B|^2)
 * alpha and beta are Gaussian exponents and must be positive.
 * Returns 0.0 if alpha or beta are non-positive.
 */
double gaussian_overlap_s(double alpha, double beta, const double A[3], const double B[3]);

#ifdef __cplusplus
}
#endif

#endif
