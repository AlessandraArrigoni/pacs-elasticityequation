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
  // Assemble the matrices separately on each subdomain
  sparseMatrixPtr_Type A1 = std::make_shared<sparseMatrix_Type> (M_nbDOF1, M_nbDOF1);
  stiffness( A1, M_uFEM1, M_CoeffFEM1, diff1, M_intMethod1);

  sparseMatrixPtr_Type A2 = std::make_shared<sparseMatrix_Type> (M_nbDOF2, M_nbDOF2);
  stiffness( A2, M_uFEM2, M_CoeffFEM2, diff2, M_intMethod2);

      #ifdef DEBUG
      std::cout << "In LaplacianBasic::assembleMatrix() the matrices have dimensions: A1 = " << A1->nrows() << "x" << A1->ncols() << " and A2 = " << A2->nrows() << "x" << A2->ncols()<<std::endl;
      #endif

  // Place the matrices in the global system
  M_Sys.addSubMatrix(A1, 0, 0);
  M_Sys.addSubMatrix(A2, M_nbDOF1, M_nbDOF1);

  sparseMatrixPtr_Type curRow = std::make_shared<sparseMatrix_Type>(1, M_nbTotDOF);

  // For each interface node...
  for (size_type k = 0; k < M_nbDOFIFace; k++)
  {
    // Get the associated row in the block from Omega1
    M_Sys.extractSubMatrix( curRow, dof_IFace1[k], 1, 0, M_nbTotDOF );
    // Add the row to the terms from Omega2 associated to the dof linked to the same grid node
    M_Sys.addSubMatrix(curRow, dof_IFace2[k] + M_nbDOF1, 0);

  }

  std::cout<< "In LaplacianBasic::assembleMatrix(), global matrix assembled     [OK]"<<std::endl;
}
