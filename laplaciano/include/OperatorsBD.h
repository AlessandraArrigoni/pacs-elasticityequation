#ifndef OPERATORSBD_H
#define OPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

//condizioni al contorno naturali

void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium,  BC* bcPtr,  FEM& femSol, FEM& femDatum, getfem::mesh_im& im);  //condizione di sforzo

#endif
