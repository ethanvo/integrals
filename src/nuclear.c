#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/integrals.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Boys function F_n(x)
 * 
 * F_n(x) = ∫[0,1] t^(2n) exp(-x*t²) dt
 * 
 * For small x, use Taylor series
 * For large x, use asymptotic expansion
 */
double boys_function(int n, double x) {
    double result;
    int i;
    
    if (x < 1e-10) {
        /* For very small x, use limiting value */
        return 1.0 / (2.0 * n + 1.0);
    }
    
    if (n == 0) {
        /* F_0(x) = (√π/2) * erf(√x) / √x */
        if (x < 30.0) {
            return 0.5 * sqrt(M_PI / x) * erf(sqrt(x));
        } else {
            /* Asymptotic expansion for large x */
            return 0.5 * sqrt(M_PI / x);
        }
    }
    
    /* Recursive relation for n > 0:
     * F_n(x) = ((2n-1) * F_(n-1)(x) - exp(-x)) / (2x)
     */
    if (x < 30.0) {
        result = boys_function(0, x);
        for (i = 1; i <= n; i++) {
            result = ((2.0 * i - 1.0) * result - exp(-x)) / (2.0 * x);
        }
        return result;
    } else {
        /* For large x, use asymptotic form */
        double f = 1.0;
        for (i = 1; i <= n; i++) {
            f *= (2.0 * i - 1.0) / (2.0 * x);
        }
        return f * 0.5 * sqrt(M_PI / x);
    }
}

/*
 * Nuclear attraction integral for s-type Gaussians
 * 
 * V = ∫ exp(-alpha1*r1²) (-Z/|r-C|) exp(-alpha2*r2²) dr
 * 
 * where C is the nuclear center
 */
double nuclear_attraction_1s(double alpha1, double ax, double ay, double az,
                             double alpha2, double bx, double by, double bz,
                             double cx, double cy, double cz) {
    double gamma, px, py, pz, xpc, ypc, zpc, rpc2;
    double prefactor, t, boys_val;
    
    /* gamma = alpha1 + alpha2 */
    gamma = alpha1 + alpha2;
    
    /* Product center P = (α₁A + α₂B) / γ */
    px = (alpha1 * ax + alpha2 * bx) / gamma;
    py = (alpha1 * ay + alpha2 * by) / gamma;
    pz = (alpha1 * az + alpha2 * bz) / gamma;
    
    /* Distance squared from P to nuclear center C */
    xpc = px - cx;
    ypc = py - cy;
    zpc = pz - cz;
    rpc2 = xpc * xpc + ypc * ypc + zpc * zpc;
    
    /* Argument for Boys function */
    t = gamma * rpc2;
    
    /* Get overlap integral */
    double s_integral = overlap_1s(alpha1, ax, ay, az, alpha2, bx, by, bz);
    
    /* Prefactor: -2π/γ */
    prefactor = -2.0 * M_PI / gamma;
    
    /* Boys function F_0(t) */
    boys_val = boys_function(0, t);
    
    return prefactor * s_integral * boys_val;
}

/*
 * Nuclear attraction integral between two contracted Gaussians
 */
double nuclear_attraction_contracted(ContractedGaussian *ga, ContractedGaussian *gb,
                                     double cx, double cy, double cz) {
    double result = 0.0;
    int i, j;
    
    for (i = 0; i < ga->nprim; i++) {
        for (j = 0; j < gb->nprim; j++) {
            double contrib = ga->primitives[i].coeff * gb->primitives[j].coeff *
                           nuclear_attraction_1s(ga->primitives[i].exponent,
                                                ga->primitives[i].x,
                                                ga->primitives[i].y,
                                                ga->primitives[i].z,
                                                gb->primitives[j].exponent,
                                                gb->primitives[j].x,
                                                gb->primitives[j].y,
                                                gb->primitives[j].z,
                                                cx, cy, cz);
            result += contrib;
        }
    }
    
    return result;
}
