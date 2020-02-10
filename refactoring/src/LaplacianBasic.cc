#include "../include/LaplacianBasic.h"


LaplacianBasic::LaplacianBasic (const GetPot& dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys):
  BasicMethod(dataFile, "laplacian", bulk1, bulk2, LaplacianBasic::Qdim, extSys ),
  diff1(dataFile, "bulkData/", "laplacian", "1", "diff"),
  diff2(dataFile, "bulkData/", "laplacian", "2", "diff")
  {
    std::cout<<"Created derived class LaplacianBasic\n"<<std::endl;
  }


void LaplacianBasic::assembleMatrix()
{
  sparseMatrixPtr_Type A1 = std::make_shared<sparseMatrix_Type> (M_nbDOF1, M_nbDOF1);
  stiffness( A1, M_uFEM1, M_CoeffFEM1, diff1, M_intMethod1);

  sparseMatrixPtr_Type A2 = std::make_shared<sparseMatrix_Type> (M_nbDOF2, M_nbDOF2);
  stiffness( A2, M_uFEM2, M_CoeffFEM2, diff2, M_intMethod2);

  //std::cout<<"Dimensioni matrici in assembleMatrix: A1 = "<<A1->nrows()<<" A2 = "<<A2->ncols()<<std::endl;

  M_Sys.addSubMatrix(A1, 0, 0);
  M_Sys.addSubMatrix(A2, M_nbDOF1, M_nbDOF1);

  sparseMatrixPtr_Type curRow = std::make_shared<sparseMatrix_Type>(1, M_nbTotDOF);

  // Get all the rows associated to the dofs on the interface and move them
  for (size_type k = 0; k < M_nbDOFIFace; k++)
  {
    M_Sys.extractSubMatrix( curRow, dof_IFace1[k], 1, 0, M_nbTotDOF );
    //std::cout << "Estratta in curRow la riga "<< dof_IFace1[k] <<std::endl;
    //std::cout<<"Aggiunta la curRow alla riga "<< dof_IFace2[k] + M_nbDOF1 <<std::endl;
    M_Sys.addSubMatrix(curRow, dof_IFace2[k] + M_nbDOF1, 0);

  }

  std::cout<< "Global matrix assembled     [OK]"<<std::endl;
}
