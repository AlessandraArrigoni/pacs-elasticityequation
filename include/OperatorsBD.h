#ifndef OPERATORSBD_H
#define OPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

/*! @file OperatorsBD.h
    @brief This file includes the method for the evaluation of natural boundary conditions.
*/

/*! method to evaluate the Neumann boundary conditions using the function getfem::asm_source_term

@param V pointer to the vector that stores the computed values
@param medium reference to the Bulk object (computational domain)
@param bcRef reference to the BC object associated to medium
@param femSol reference to the FEM space for the test functions
@param femDatum reference to the FEM space for the Neumann datum
@param im reference to the object describing the integration method
*/
void stressRHS( scalarVectorPtr_Type V, const Bulk& medium, BC& bcRef, const FEM& femSol, const FEM& femDatum, const getfem::mesh_im& im);

#endif
