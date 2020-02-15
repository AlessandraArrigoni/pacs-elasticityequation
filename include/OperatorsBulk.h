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

/*! method to compute the stiffness matrix for the elliptic scalar problem */
void stiffness ( sparseMatrixPtr_Type M, const FEM& femSol, const FEM& femCoef, BulkDatum& Diff, const getfem::mesh_im& im);

/*! method to compute the stiffness matrix for the linear elasticity problem */
void linearElasticity(sparseMatrixPtr_Type M, const FEM& femSol, const FEM& femCoef, BulkDatum& Mu, BulkDatum& Lambda, const getfem::mesh_im& im);

/*! method to compute the volumetric source term */
void bulkLoad (scalarVectorPtr_Type V, const FEM& FemSol, const FEM& FemSource, BulkDatum& Source, const getfem::mesh_im& im);  

/*! method to evaluate the exact solution in every degree of freedom of the FEM */
void exactSolution (scalarVectorPtr_Type V, const FEM& FemSol, BulkDatum& Solution); 

/*! method to evaluate the jump in every degree of freedom of the FEM */
void jump(scalarVectorPtr_Type V, const FEM& FemSol, BulkDatum& Jump); 

#endif
