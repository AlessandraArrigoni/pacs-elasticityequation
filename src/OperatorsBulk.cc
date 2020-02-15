#include "../include/OperatorsBulk.h"

void stiffness ( sparseMatrixPtr_Type M, const FEM& FemSol, const FEM& FemCoef, BulkDatum& Diff, const getfem::mesh_im& im)
{
    getfem::mesh_fem femSol(FemSol.getFEM());
    getfem::mesh_fem femCoef(FemCoef.getFEM());

    std::vector<scalar_type> diff(femCoef.nb_dof(),0.0);

    // Compute non costant diffusion coefficient
    for (int i=0; i<femCoef.nb_dof();++i)
    {
    	diff[i] = Diff.getValue(femCoef.point_of_basic_dof(i), 0);
    }

    getfem::asm_stiffness_matrix_for_laplacian(*M, im, femSol, femCoef, diff);


    std::cout << "operator stiffness   :    [OK]" << std::endl;

}


void linearElasticity(sparseMatrixPtr_Type M, const FEM& FemSol, const FEM& FemCoef, BulkDatum& Mu, BulkDatum& Lambda, const getfem::mesh_im& im)
{
  getfem::mesh_fem femSol(FemSol.getFEM());
  getfem::mesh_fem femCoef(FemCoef.getFEM());

  scalarVector_Type lambda(femCoef.nb_dof(),0.0);
  scalarVector_Type mu(femCoef.nb_dof(),0.0);

  // Compute non costant Lamè coefficients
  for (size_type i=0; i<femCoef.nb_dof();++i)
  {
    lambda[i] = Lambda.getValue(femCoef.point_of_basic_dof(i), 0);
    mu[i] = Mu.getValue(femCoef.point_of_basic_dof(i), 0);
  }

  getfem::asm_stiffness_matrix_for_linear_elasticity(*M, im, femSol, femCoef, lambda, mu);

  std::cout << "operator linearElasticity      [OK]" << std::endl;
}



void bulkLoad(scalarVectorPtr_Type V, const FEM& FemSol, const FEM& FemSource, BulkDatum& Source, const getfem::mesh_im& im)
{
  size_type Qdim = FemSol.getFEM().get_qdim(); // I need it here since the Qdim property is defined in the Problem class, not here!

  getfem::mesh_fem femSol(FemSol.getFEM());
  getfem::mesh_fem femSource(FemSource.getFEM());

  scalarVector_Type sourceVett(femSource.nb_dof(), 0.0);

  // Compute source term in each dof
  for (size_type i = 0; i < femSource.nb_dof(); i += Qdim)
  {
    for (size_type j = 0; j < Qdim; j++){
      sourceVett[i + j] = Source.getValue(femSource.point_of_basic_dof(i), j);
    }
  }

  getfem::asm_source_term(*V, im, femSol, femSource, sourceVett);

  
  std::cout << "operator bulkLoad      [OK]" << std::endl;
}



void exactSolution(scalarVectorPtr_Type V, const FEM& FemSol, BulkDatum& Solution)
{
  size_type Qdim = FemSol.getFEM().get_qdim(); // I need it here since the Qdim property is defined in the Problem class, not here!

  for (size_type i=0; i<FemSol.nb_dof(); i+= Qdim)
  {
    for (size_type j = 0; j < Qdim; j++){
      V->at(i+j) = Solution.getValue(FemSol.point_of_basic_dof(i), j);

    }
  }
}

void jump(scalarVectorPtr_Type V, const FEM& FemSol, BulkDatum& Jump )
{
  size_type Qdim = FemSol.getFEM().get_qdim();

  for (size_type i=0; i<FemSol.nb_dof(); i+= Qdim)
  {
    for (size_type j = 0; j < Qdim; j++){
      V->at(i+j) = Jump.getValue(FemSol.point_of_basic_dof(i), j);

    }
  }
}
