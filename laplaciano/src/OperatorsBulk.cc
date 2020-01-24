#include "../include/OperatorsBulk.h"

void stiffness ( sparseMatrixPtr_Type M, Bulk* medium, FEM& FemSol, FEM& FemCoef, getfem::mesh_im& im)
{
    getfem::mesh_fem femSol(*(FemSol.getFEM()));
    getfem::mesh_fem femCoef(*(FemCoef.getFEM()));

    std::vector<scalar_type> mu(femCoef.nb_dof(),0.0);

    // Compute non costant diffusion coefficient
    for (int i=0; i<femCoef.nb_dof();++i)
    {
    	mu[i]=medium->getData()->getDiff(i);
    }

    getfem::asm_stiffness_matrix_for_laplacian(*M, im, femSol, femCoef, mu);

    std::cout << "operator a(volume)   :    [OK]" << std::endl;

}

void bulkLoad( scalarVectorPtr_Type V,
               Bulk* medium, FEM& FemSol, FEM& FemSource, getfem::mesh_im& im)
{
  getfem::mesh_fem femSol(*(FemSol.getFEM()));
  getfem::mesh_fem femSource(*(FemSource.getFEM()));

  scalarVector_Type sourceVett(femSource.nb_dof(),0.0);

  // Compute source term in each dof
  for (size_type i=0; i<femSource.nb_dof();++i)
  {
  	sourceVett[ i ] = medium->getData()->bulkLoad(femSource.point_of_basic_dof(i))[0]; // chiama il metodo bulkLoad contenuto nella classe getData che è un membro della classe Bulk
  }

  getfem::asm_source_term(*V, im, femSol, femSource, sourceVett);

  std::cout << "operator bulkLoad(volume)   :    [OK]" << std::endl;

}


void exactSolution(scalarVectorPtr_Type V, Bulk* medium, FEM& FemSol)
{
  for (size_type i=0; i<FemSol.nb_dof(); i++){
    V->at(i) = medium->getData()->exactSolution(FemSol.point_of_basic_dof(i))[0];
  }
}
