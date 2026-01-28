#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/integrals.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Kinetic energy integral between two primitive s-type Gaussian functions
 * 
 * T = ∫ exp(-alpha1*r1²) (-1/2 ∇²) exp(-alpha2*r2²) dr
 * 
 * For s-type Gaussians:
 * T = alpha2*(3 - 2*alpha2*R²) * S
 * where S is the overlap integral and R is the distance between centers
 */
double kinetic_1s(double alpha1, double ax, double ay, double az,
                  double alpha2, double bx, double by, double bz) {
    double rab2, s_integral, term;
    
    /* Distance squared between centers */
    rab2 = (ax - bx) * (ax - bx) + 
           (ay - by) * (ay - by) + 
           (az - bz) * (az - bz);
    
    /* Get overlap integral */
    s_integral = overlap_1s(alpha1, ax, ay, az, alpha2, bx, by, bz);
    
    /* Kinetic energy expression */
    term = alpha1 * alpha2 / (alpha1 + alpha2);
    return term * (3.0 - 2.0 * term * rab2) * s_integral;
}

/*
 * Kinetic energy integral between two contracted Gaussian functions
 */
double kinetic_contracted(ContractedGaussian *ga, ContractedGaussian *gb) {
    double result = 0.0;
    int i, j;
    
    for (i = 0; i < ga->nprim; i++) {
        for (j = 0; j < gb->nprim; j++) {
            double contrib = ga->primitives[i].coeff * gb->primitives[j].coeff *
                           kinetic_1s(ga->primitives[i].exponent,
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
