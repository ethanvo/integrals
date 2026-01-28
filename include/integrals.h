#ifndef INTEGRALS_H
#define INTEGRALS_H

#include <math.h>

/* Gaussian primitive structure */
typedef struct {
    double exponent;   /* Gaussian exponent (alpha) */
    double coeff;      /* Contraction coefficient */
    double x, y, z;    /* Center coordinates */
} GaussianPrimitive;

/* Contracted Gaussian function */
typedef struct {
    int nprim;                    /* Number of primitives */
    GaussianPrimitive *primitives; /* Array of primitives */
    double x, y, z;               /* Center coordinates */
} ContractedGaussian;

/* Function prototypes */

/* Overlap integral between two s-type Gaussians */
double overlap_1s(double alpha1, double ax, double ay, double az,
                  double alpha2, double bx, double by, double bz);

/* Overlap integral between two contracted Gaussians */
double overlap_contracted(ContractedGaussian *ga, ContractedGaussian *gb);

/* Kinetic energy integral between two s-type Gaussians */
double kinetic_1s(double alpha1, double ax, double ay, double az,
                  double alpha2, double bx, double by, double bz);

/* Kinetic energy integral between two contracted Gaussians */
double kinetic_contracted(ContractedGaussian *ga, ContractedGaussian *gb);

/* Nuclear attraction integral for s-type Gaussians */
double nuclear_attraction_1s(double alpha1, double ax, double ay, double az,
                             double alpha2, double bx, double by, double bz,
                             double cx, double cy, double cz);

/* Nuclear attraction integral between two contracted Gaussians */
double nuclear_attraction_contracted(ContractedGaussian *ga, ContractedGaussian *gb,
                                     double cx, double cy, double cz);

/* Two-electron repulsion integral for s-type Gaussians */
double electron_repulsion_1s(double alpha1, double ax, double ay, double az,
                             double alpha2, double bx, double by, double bz,
                             double alpha3, double cx, double cy, double cz,
                             double alpha4, double dx, double dy, double dz);

/* Two-electron repulsion integral between four contracted Gaussians */
double electron_repulsion_contracted(ContractedGaussian *ga, ContractedGaussian *gb,
                                     ContractedGaussian *gc, ContractedGaussian *gd);

/* Utility functions */
double boys_function(int n, double x);
ContractedGaussian* create_contracted_gaussian(int nprim, double x, double y, double z);
void free_contracted_gaussian(ContractedGaussian *cg);

#endif /* INTEGRALS_H */
