#ifndef OPERATORSBULK_H
#define OPERATORSBULK_H

#include "Core.h"
#include "FEM.h"
#include "Bulk.h"
#include "BC.h"
#include "BulkDatum.h"

//operatori nel bulk

void stiffness ( sparseMatrixPtr_Type M, const FEM& femSol, const FEM& femCoef, BulkDatum& Diff, const getfem::mesh_im& im);
void linearElasticity(sparseMatrixPtr_Type M, const FEM& femSol, const FEM& femCoef, BulkDatum& Mu, BulkDatum& Lambda, const getfem::mesh_im& im);
void bulkLoad (scalarVectorPtr_Type V, const FEM& FemSol, const FEM& FemSource, BulkDatum& Source, const getfem::mesh_im& im);  //termine forzante volumetrico
void exactSolution (scalarVectorPtr_Type V, const FEM& FemSol, BulkDatum& Solution); //vettore con la soluzione esatta in ogni nodo

#endif
