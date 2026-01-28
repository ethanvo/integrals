#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/integrals.h"

/*
 * Create a contracted Gaussian function
 */
ContractedGaussian* create_contracted_gaussian(int nprim, double x, double y, double z) {
    ContractedGaussian *cg = (ContractedGaussian*)malloc(sizeof(ContractedGaussian));
    if (cg == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for ContractedGaussian\n");
        return NULL;
    }
    
    cg->nprim = nprim;
    cg->x = x;
    cg->y = y;
    cg->z = z;
    
    cg->primitives = (GaussianPrimitive*)malloc(nprim * sizeof(GaussianPrimitive));
    if (cg->primitives == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for primitives\n");
        free(cg);
        return NULL;
    }
    
    /* Initialize primitives to the same center */
    for (int i = 0; i < nprim; i++) {
        cg->primitives[i].x = x;
        cg->primitives[i].y = y;
        cg->primitives[i].z = z;
        cg->primitives[i].exponent = 0.0;
        cg->primitives[i].coeff = 0.0;
    }
    
    return cg;
}

/*
 * Free a contracted Gaussian function
 */
void free_contracted_gaussian(ContractedGaussian *cg) {
    if (cg != NULL) {
        if (cg->primitives != NULL) {
            free(cg->primitives);
        }
        free(cg);
    }
}
