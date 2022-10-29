#ifndef VEGAFEM_GETINTEGRATORSOLVER_H
#define VEGAFEM_GETINTEGRATORSOLVER_H

// returns the string corresponding to the selected integrator solver
// "solver" must be pre-allocated
// result: PARDISO, SPOOLES or PCG
void GetIntegratorSolver(char * solver);
 
#endif

