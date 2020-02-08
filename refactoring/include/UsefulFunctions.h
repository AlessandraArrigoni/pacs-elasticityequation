#ifndef USEFULFUNCTIONS_H
#define USEFULFUNCTIONS_H

// AGGIUNGEREI QUA LA FUNZIONE INDEX_INSIDE facendole restituire un booleano e ricevendo in input un sizeVector_Type

#include "Core.h"
#include "FEM.h"

void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution );

void exportSolutionInCell ( const std::string& fileName,
                            const std::string& solutionName,
                            const getfem::mesh_fem& meshFEM,
                            const scalarVector_Type& solution );

void exportMesh ( const std::string& fileName, const getfem::mesh& mesh );

void massLumping ( sparseMatrix_Type& matrix );

void fromBitVectorToStdVector ( dal::bit_vector& bitVector,
                                std::vector < size_type >& stdVector );

char intToChar ( const size_type& integer );

scalar_type pointDistance ( const scalar_type& x0,
                            const scalar_type& x1,
                            const scalar_type& y0,
                            const scalar_type& y1 );


#endif /* USEFULFUNCTIONS_H_ */
