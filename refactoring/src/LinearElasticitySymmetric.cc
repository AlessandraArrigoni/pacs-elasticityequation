#include "../include/LinearElasticitySymmetric.h"

LinearElasticitySymmetric::LinearElasticitySymmetric( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys):
  SymmetricMethod(dataFile, "elasticity", bulk1, bulk2, LinearElasticitySymmetric::Qdim, extSys),
  mu1(dataFile, "bulkData/", "elasticity","1", "mu"),
  lambda1(dataFile, "bulkData/", "elasticity","1", "lambda"),
  mu2(dataFile, "bulkData/", "elasticity","2", "mu"),
  lambda2(dataFile, "bulkData/", "elasticity","2", "lambda")
  {
    std::cout<<"Created derived class LinearElasticitySymmetric\n"<<std::endl;
  }


void LinearElasticitySymmetric::assembleMatrix()
{
	sparseMatrixPtr_Type A1 = std::make_shared<sparseMatrix_Type> (M_nbDOF1, M_nbDOF1);
	linearElasticity( A1, M_uFEM1, M_CoeffFEM1, mu1, lambda1, M_intMethod1);

	sparseMatrixPtr_Type A2 = std::make_shared<sparseMatrix_Type> (M_nbDOF2, M_nbDOF2);
	linearElasticity(A2,  M_uFEM2, M_CoeffFEM2, mu2, lambda2, M_intMethod2);


  // dal blocco [A_1Gamma A_1GammaGamma 0 0] metto a zero i termini A_1GammaGamma e li sommo a A_2GammaGamma:
	for (size_type k=0; k < M_nbDOFIFace; k++)
	{
	  size_type idx1=dof_IFace1[k];
    scalar_type value= (*A1)(idx1, idx1);
    size_type idx2=dof_IFace2[k];
    (*A2)(idx2, idx2)+= value;
	}


	M_Sys.addSubMatrix(A1, 0, 0);
	M_Sys.addSubMatrix(A2, M_nbDOF1, M_nbDOF1);

	// collego i termini legati all'interfaccia del dominio di sinistra a quelli dell'interfaccia del dominio di destra

	for (size_type k=0; k < M_nbDOFIFace; k++)
	{
		sparseMatrixPtr_Type curRow1 = std::make_shared<sparseMatrix_Type>(1,M_nbTotDOF);
		sparseMatrixPtr_Type curCol2 = std::make_shared<sparseMatrix_Type>(M_nbTotDOF,1);
	  M_Sys.extractSubMatrix( curRow1, dof_IFace1[k], 1, 0, M_nbTotDOF );
		M_Sys.extractSubMatrix( curCol2, 0, M_nbTotDOF, dof_IFace1[k], 1 );
	  M_Sys.addSubMatrix(curRow1, dof_IFace2[k] + M_nbDOF1, 0);
		M_Sys.addSubMatrix(curCol2, 0, dof_IFace2[k] + M_nbDOF1);
	}

  std::cout<< "Global matrix assembled     [OK]"<<std::endl;
}
