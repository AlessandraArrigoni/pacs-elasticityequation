#include "../include/LinearElasticityPb.h"

LinearElasticityPb::LinearElasticityPb( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys):
  Problem(dataFile, "elasticity", bulk1, bulk2, LinearElasticityPb::Qdim, extSys),
  mu1(dataFile, "bulkData/", "elasticity","1", "mu"),
  lambda1(dataFile, "bulkData/", "elasticity","1", "lambda"),
  mu2(dataFile, "bulkData/", "elasticity","2", "mu"),
  lambda2(dataFile, "bulkData/", "elasticity","2", "lambda")
  {
    std::cout<<"Created derived class LinearElasticityPb\n"<<std::endl;
  }


void LinearElasticityPb::assembleMatrix()
{
  sparseMatrixPtr_Type A1 = std::make_shared<sparseMatrix_Type> (M_nbDOF1, M_nbDOF1);
  linearElasticity( A1, M_uFEM1, M_CoeffFEM1, mu1, lambda1, M_intMethod1);

  sparseMatrixPtr_Type A2 = std::make_shared<sparseMatrix_Type> (M_nbDOF2, M_nbDOF2);
  linearElasticity( A2, M_uFEM2, M_CoeffFEM2, mu2, lambda2, M_intMethod2);

  M_Sys.addSubMatrix(A1, 0, 0);
  M_Sys.addSubMatrix(A2, M_nbDOF1, M_nbDOF1);

  // Get all the rows associated to the dofs on the interface and move them
  for (size_type k = 0; k < dof_IFace1.size(); k++)
  {
    sparseMatrixPtr_Type curRow = std::make_shared<sparseMatrix_Type>(1, M_nbTotDOF);

    M_Sys.extractSubMatrix( curRow, dof_IFace1[k], 1, 0, M_nbTotDOF);
    M_Sys.addSubMatrix(curRow, dof_IFace2[k] + M_nbDOF1, 0);
  }

  std::cout<< "Global matrix assembled     [OK]"<<std::endl;
}
