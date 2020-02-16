#ifndef OPERATORSBULK_H
#define OPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"
#include "BulkDatum.h"

/*! @file OperatorsBulk.h
    @brief This file assembles different methods related to the bulk that can be employed in several contexts.

*/

/*! method to compute the stiffness matrix for the elliptic scalar problem

@param M pointer to the matrix that stores the computed values
@param FemSol reference to the FEM space for the test and trial functions
@param FemCoef reference to the FEM space for the diffusion coefficient
@param Diff reference to the BulkDatum object for the diffusion coefficient
@param im reference to the object describing the integration method*/
void stiffness ( sparseMatrixPtr_Type M, const FEM& FemSol, const FEM& FemCoef, BulkDatum& Diff, const getfem::mesh_im& im);


/*! method to compute the stiffness matrix for the linear elasticity problem

@param M pointer to the matrix that stores the computed values
@param FemSol reference to the FEM space for the test and trial functions
@param FemCoef reference to the FEM space for the Lamè coefficients
@param Mu reference to the BulkDatum object for the Lamè coefficient
@param Lambda reference to the BulkDatum object for the Lamè coefficient
@param im reference to the object describing the integration method*/
void linearElasticity(sparseMatrixPtr_Type M, const FEM& FemSol, const FEM& FemCoef, BulkDatum& Mu, BulkDatum& Lambda, const getfem::mesh_im& im);


/*! method to compute the volumetric source term

@param V pointer to the vector that stores the computed values
@param FemSol reference to the FEM space for the test functions
@param FemSource reference to the FEM space for the volumetric source term
@param Source reference to the BulkDatum object for the volumetric source term
@param im reference to the object describing the integration method*/
void bulkLoad (scalarVectorPtr_Type V, const FEM& FemSol, const FEM& FemSource, BulkDatum& Source, const getfem::mesh_im& im);


/*! method to evaluate the exact solution in every degree of freedom of the FEM space

@param V pointer to the vector that stores the computed values
@param FemSol reference to the FEM space for the trial functions
@param Solution reference to the BulkDatum object for the exact solution
*/
void exactSolution (scalarVectorPtr_Type V, const FEM& FemSol, BulkDatum& Solution);

/*! method to evaluate the jump q0 in every degree of freedom of the FEM space

@param V pointer to the vector that stores the computed values
@param FemSol reference to the FEM space for the trial functions
@param Jump reference to the BulkDatum object for the interface jump q0
*/
void jump(scalarVectorPtr_Type V, const FEM& FemSol, BulkDatum& Jump);

#endif
