#ifndef USEFULFUNCTIONS_H
#define USEFULFUNCTIONS_H

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


std::pair < std::string, size_type > comparaSegni ( const std::string& region,
                                                    const scalarVector_Type& signs);

bool isInTriangle ( const getfem::mesh& mesh,
                    const size_type& elementID,
                    const base_node& node,
                    const scalar_type& toll = 1e-7 );

bool isInTriangle ( const base_node& v1, const base_node& v2, const base_node& v3,
                    const base_node& node,
                    const scalar_type& toll = 1e-7 );

bool intersectSegments(const base_node& a1, const base_node& a2, const base_node& b1, const base_node& b2, base_node& sol);

bool intersectSegmentTriangle(const base_node& a1, const base_node& a2, const base_node& v1, const base_node& v2, const base_node& v3,  base_node& s1, base_node& s2);



#endif /* USEFULFUNCTIONS_H_ */
