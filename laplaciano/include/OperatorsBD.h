#ifndef OPERATORSBD_H
#define OPERATORSBD_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

//condizioni al contorno naturali

void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time, BC* bcPtr,  FEM& femV, FEM& femP, getfem::mesh_im& im);  //condizione di sforzo

#endif
