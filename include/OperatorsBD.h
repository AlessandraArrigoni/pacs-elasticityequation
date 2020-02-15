#ifndef OPERATORSBD_H
#define OPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

/*! @file OperatorsBD.h
    @brief This file includes the method for the evaluation of natural boundary conditions.
*/

/*! method to evaluate the Neumann boundary conditions using the function getfem::asm_source_term */
void stressRHS( scalarVectorPtr_Type V, const Bulk& medium, BC& bcRef, const FEM& femSol, const FEM& femDatum, const getfem::mesh_im& im); 

#endif
