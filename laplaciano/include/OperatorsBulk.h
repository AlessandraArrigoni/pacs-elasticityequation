#ifndef OPERATORSBULK_H
#define OPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

//operatori elasticità nel bulk

void stiffness ( sparseMatrixPtr_Type M, Bulk* medium, FEM& femV, FEM& femP, getfem::mesh_im& im);  //matrice di stiffness
void bulkLoad (scalarVectorPtr_Type V, Bulk* medium,  FEM& FemD, FEM& FemC, getfem::mesh_im& im);  //termine forzante volumetrico
void massMatrix(sparseMatrixPtr_Type M,  Bulk* medium,  FEM& FemD, getfem::mesh_im& im);  //matrice di massa
void exactSolution(scalarVectorPtr_Type V, Bulk* medium, FEM& FemD); //vettore con la soluzione esatta in ogni nodo

#endif
