#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/integrals.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Two-electron repulsion integral (ERI) for s-type Gaussians
 * 
 * (ab|cd) = ∫∫ φ_a(r1) φ_b(r1) (1/r12) φ_c(r2) φ_d(r2) dr1 dr2
 * 
 * For s-type Gaussians, using the Gaussian product theorem and
 * the Coulomb operator expansion
 */
double electron_repulsion_1s(double alpha1, double ax, double ay, double az,
                             double alpha2, double bx, double by, double bz,
                             double alpha3, double cx, double cy, double cz,
                             double alpha4, double dx, double dy, double dz) {
    double gamma1, gamma2;
    double px, py, pz, qx, qy, qz;
    double xpq, ypq, zpq, rpq2;
    double rab2, rcd2;
    double prefactor, t, boys_val;
    
    /* gamma1 = alpha1 + alpha2, gamma2 = alpha3 + alpha4 */
    gamma1 = alpha1 + alpha2;
    gamma2 = alpha3 + alpha4;
    
    /* Product centers P and Q */
    px = (alpha1 * ax + alpha2 * bx) / gamma1;
    py = (alpha1 * ay + alpha2 * by) / gamma1;
    pz = (alpha1 * az + alpha2 * bz) / gamma1;
    
    qx = (alpha3 * cx + alpha4 * dx) / gamma2;
    qy = (alpha3 * cy + alpha4 * dy) / gamma2;
    qz = (alpha3 * cz + alpha4 * dz) / gamma2;
    
    /* Distance squared between P and Q */
    xpq = px - qx;
    ypq = py - qy;
    zpq = pz - qz;
    rpq2 = xpq * xpq + ypq * ypq + zpq * zpq;
    
    /* Distance squared between A and B, C and D */
    rab2 = (ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz);
    rcd2 = (cx - dx) * (cx - dx) + (cy - dy) * (cy - dy) + (cz - dz) * (cz - dz);
    
    /* Argument for Boys function */
    t = gamma1 * gamma2 / (gamma1 + gamma2) * rpq2;
    
    /* Prefactor */
    prefactor = 2.0 * pow(M_PI, 2.5) / (gamma1 * gamma2 * sqrt(gamma1 + gamma2));
    prefactor *= exp(-alpha1 * alpha2 / gamma1 * rab2);
    prefactor *= exp(-alpha3 * alpha4 / gamma2 * rcd2);
    
    /* Boys function F_0(t) */
    boys_val = boys_function(0, t);
    
    return prefactor * boys_val;
}

/*
 * Two-electron repulsion integral between four contracted Gaussians
 */
double electron_repulsion_contracted(ContractedGaussian *ga, ContractedGaussian *gb,
                                     ContractedGaussian *gc, ContractedGaussian *gd) {
    double result = 0.0;
    int i, j, k, l;
    
    for (i = 0; i < ga->nprim; i++) {
        for (j = 0; j < gb->nprim; j++) {
            for (k = 0; k < gc->nprim; k++) {
                for (l = 0; l < gd->nprim; l++) {
                    double contrib = ga->primitives[i].coeff *
                                   gb->primitives[j].coeff *
                                   gc->primitives[k].coeff *
                                   gd->primitives[l].coeff *
                                   electron_repulsion_1s(
                                       ga->primitives[i].exponent,
                                       ga->primitives[i].x,
                                       ga->primitives[i].y,
                                       ga->primitives[i].z,
                                       gb->primitives[j].exponent,
                                       gb->primitives[j].x,
                                       gb->primitives[j].y,
                                       gb->primitives[j].z,
                                       gc->primitives[k].exponent,
                                       gc->primitives[k].x,
                                       gc->primitives[k].y,
                                       gc->primitives[k].z,
                                       gd->primitives[l].exponent,
                                       gd->primitives[l].x,
                                       gd->primitives[l].y,
                                       gd->primitives[l].z);
                    result += contrib;
                }
            }
        }
    }
    
    return result;
}
