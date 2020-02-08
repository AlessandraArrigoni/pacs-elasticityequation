#ifndef OPERATORSBULK_H
#define OPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"

//operatori elasticità nel bulk

void stiffness ( sparseMatrixPtr_Type M, Bulk* medium, FEM& FemSol, FEM& FemCoef, getfem::mesh_im& im);  //matrice di stiffness
void bulkLoad (scalarVectorPtr_Type V, Bulk* medium,  FEM& FemSol, FEM& FemSource, getfem::mesh_im& im);  //termine forzante volumetrico
void massMatrix(sparseMatrixPtr_Type M,  Bulk* medium,  FEM& FemD, getfem::mesh_im& im);  //matrice di massa
void exactSolution(scalarVectorPtr_Type V, Bulk* medium, FEM& FemSol); //vettore con la soluzione esatta in ogni nodo
void jump(scalarVectorPtr_Type V, Bulk* medium, FEM& FemSol);

#endif
