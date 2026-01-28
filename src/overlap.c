#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/integrals.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* 
 * Overlap integral between two primitive s-type Gaussian functions
 * 
 * S = ∫ exp(-alpha1*r1²) exp(-alpha2*r2²) dr
 * 
 * For s-type Gaussians centered at A and B:
 * S = (π/(alpha1+alpha2))^(3/2) * exp(-alpha1*alpha2/(alpha1+alpha2)*|A-B|²)
 */
double overlap_1s(double alpha1, double ax, double ay, double az,
                  double alpha2, double bx, double by, double bz) {
    double rab2, gamma, prefactor;
    
    /* Distance squared between centers */
    rab2 = (ax - bx) * (ax - bx) + 
           (ay - by) * (ay - by) + 
           (az - bz) * (az - bz);
    
    /* gamma = alpha1 + alpha2 */
    gamma = alpha1 + alpha2;
    
    /* Prefactor: (π/γ)^(3/2) */
    prefactor = pow(M_PI / gamma, 1.5);
    
    /* Exponential factor: exp(-α₁α₂/γ * R²) */
    return prefactor * exp(-alpha1 * alpha2 / gamma * rab2);
}

/*
 * Overlap integral between two contracted Gaussian functions
 */
double overlap_contracted(ContractedGaussian *ga, ContractedGaussian *gb) {
    double result = 0.0;
    int i, j;
    
    for (i = 0; i < ga->nprim; i++) {
        for (j = 0; j < gb->nprim; j++) {
            double contrib = ga->primitives[i].coeff * gb->primitives[j].coeff *
                           overlap_1s(ga->primitives[i].exponent,
                                     ga->primitives[i].x,
                                     ga->primitives[i].y,
                                     ga->primitives[i].z,
                                     gb->primitives[j].exponent,
                                     gb->primitives[j].x,
                                     gb->primitives[j].y,
                                     gb->primitives[j].z);
            result += contrib;
        }
    }
    
    return result;
}
