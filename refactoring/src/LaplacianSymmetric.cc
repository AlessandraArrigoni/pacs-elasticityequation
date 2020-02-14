#include "../include/LaplacianSymmetric.h"


LaplacianSymmetric::LaplacianSymmetric (const GetPot& dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys):
  SymmetricMethod(dataFile, "laplacian", bulk1, bulk2, LaplacianSymmetric::Qdim, extSys ),
  diff1(dataFile, "bulkData/", "laplacian", "1", "diff"),
  diff2(dataFile, "bulkData/", "laplacian", "2", "diff")
  {
    std::cout<<"Created derived class LaplacianSymmetric\n"<<std::endl;
  }


void LaplacianSymmetric::assembleMatrix()
{
  // Assemble the matrices separately on each subdomain
	sparseMatrixPtr_Type A1 = std::make_shared<sparseMatrix_Type> (M_nbDOF1, M_nbDOF1);
	stiffness( A1, M_uFEM1, M_CoeffFEM1, diff1, M_intMethod1);

  sparseMatrixPtr_Type A2 = std::make_shared<sparseMatrix_Type> (M_nbDOF2, M_nbDOF2);
	stiffness( A2, M_uFEM2, M_CoeffFEM2, diff2, M_intMethod2);

      #ifdef DEBUG
      std::cout << "In LaplacianSymmetric::assembleMatrix() the matrices have dimensions: A1 = " << A1->nrows() << "x" << A1->ncols() << " and A2 = " << A2->nrows() << "x" << A2->ncols()<<std::endl;
      #endif

  // Construct the term A_Gamma,Gamma: for each interface node...
	for (size_type k = 0; k < M_nbDOFIFace; k++)
  {
    // Get the values in A_Gamma1,Gamma1
	  size_type idx1 = dof_IFace1[k];
    scalar_type value = (*A1)(idx1, idx1);

    // Sum them to A_Gamma2,Gamma2
    size_type idx2 = dof_IFace2[k];
    (*A2)(idx2, idx2) += value;
  }

  // Place the matrices in the global system
	M_Sys.addSubMatrix(A1, 0, 0);
	M_Sys.addSubMatrix(A2, M_nbDOF1, M_nbDOF1);

  sparseMatrixPtr_Type curRow1 = std::make_shared<sparseMatrix_Type>(1, M_nbTotDOF);
  sparseMatrixPtr_Type curCol2 = std::make_shared<sparseMatrix_Type>(M_nbTotDOF, 1);

  // Link the interface dofs from Omega1 to the interface dofs from Omega2
  // For each interface node...
  for (size_type k=0; k < M_nbDOFIFace; k++)
	{
    // Get the row from A_Gamma,1
    M_Sys.extractSubMatrix( curRow1, dof_IFace1[k], 1, 0, M_nbTotDOF );
    // Get the column from A_1,Gamma
 	  M_Sys.extractSubMatrix( curCol2, 0, M_nbTotDOF, dof_IFace1[k], 1 );
    // Add the row to the equation associated to the same node from Omega2
	  M_Sys.addSubMatrix(curRow1, dof_IFace2[k] + M_nbDOF1, 0);
    // Add the column to include the dependence on the unknown average
		M_Sys.addSubMatrix(curCol2, 0, dof_IFace2[k] + M_nbDOF1);
  }

  std::cout<< "In LaplacianSymmetric::assembleMatrix() global matrix assembled     [OK]"<<std::endl;
}
