#ifndef OPERATORSBD_H
#define OPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

//condizioni al contorno naturali

void stressRHS( scalarVectorPtr_Type V, const Bulk& medium, BC& bcRef, const FEM& femSol, const FEM& femDatum, const getfem::mesh_im& im);  //condizione di sforzo

#endif
