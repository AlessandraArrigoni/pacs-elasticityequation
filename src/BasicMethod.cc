#include "../include/BasicMethod.h"


// constructor 
BasicMethod::BasicMethod(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys):
  Problem(dataFile, problem, bulk1, bulk2, dim, extSys) {}


void BasicMethod::assembleRHS()
{
  // Volumetric source term
  scalarVectorPtr_Type source1 = std::make_shared<scalarVector_Type> (M_nbDOF1);
	bulkLoad(source1, M_uFEM1, M_uFEM1, M_source1, M_intMethod1);

  scalarVectorPtr_Type source2 = std::make_shared<scalarVector_Type> (M_nbDOF2);
	bulkLoad(source2, M_uFEM2, M_uFEM2, M_source2, M_intMethod2);

  // Put the source for Omega1 at the top of the rhs
	M_Sys.addSubVector(source1, 0);
  // Put the source for Omega2 at the bottom of the rhs
	M_Sys.addSubVector(source2, M_nbDOF1);

	// For each interface node...
	for (size_type k = 0; k < M_nbDOFIFace; k++)
  {
    // Sum the values from Omega1 and Omega2
		scalar_type newVal = (M_Sys.getRHS())->at(dof_IFace1[k]) + (M_Sys.getRHS())->at(dof_IFace2[k] + M_nbDOF1);
    // Set the new value of the interface dof on Omega2
		M_Sys.setRHSValue(dof_IFace2[k] + M_nbDOF1, newVal);
	}

  // Set Neumann boundary conditions
	scalarVectorPtr_Type BCvec1 = std::make_shared<scalarVector_Type> (M_nbDOF1);
	stressRHS( BCvec1, M_Bulk1, M_BC1, M_uFEM1, M_uFEM1, M_intMethod1);

	scalarVectorPtr_Type BCvec2 = std::make_shared<scalarVector_Type> (M_nbDOF2);
	stressRHS( BCvec2, M_Bulk2, M_BC2, M_uFEM2, M_uFEM2, M_intMethod2);

	M_Sys.addSubVector(BCvec1, 0);
	M_Sys.addSubVector(BCvec2, M_nbDOF1);
}


void BasicMethod::enforceStrongBC(size_type const domainIdx)
{
		// Set local variables according to the current domain
		FEM * M_uFEMcur;
		BC * M_BCcur;
		sizeVector_Type  M_rowsStrongBCcur;
		sizeVector_Type  M_rowsStrongBCFlagscur;

		if (domainIdx == 1){
			M_uFEMcur = &M_uFEM1;
			M_BCcur = &M_BC1;
			M_rowsStrongBCcur = M_rowsStrongBC1;
			M_rowsStrongBCFlagscur = M_rowsStrongBCFlags1;
		}
		else if (domainIdx == 2){
			M_uFEMcur = &M_uFEM2;
			M_BCcur = &M_BC2;
			M_rowsStrongBCcur = M_rowsStrongBC2;
			M_rowsStrongBCFlagscur = M_rowsStrongBCFlags2;
		}

    size_type Qdim = M_uFEMcur->getFEM().get_qdim();

    // Store the dofs on all the Dirichlet sides in M_rowsStrongBCcur
    // Loop on the Dirichlet sides
		for ( size_type bndID = 0; bndID < M_BCcur->getDiriBD().size(); bndID++ )
		{
      // Get the dofs on the current side
			dal::bit_vector quali = M_uFEMcur->getFEM().dof_on_region( M_BCcur->getDiriBD()[bndID]);
			for(dal::bv_visitor i(quali); !i.finished(); ++i)
			{
				M_rowsStrongBCcur.push_back(i);
				M_rowsStrongBCFlagscur.push_back(bndID);
			}
		}

    // For each dof on the Dirichlet boundaries
		for (size_type i = 0; i < M_rowsStrongBCcur.size(); i += Qdim)
		{
			size_type ii = M_rowsStrongBCcur[i] ;

      // For each dimension of the problem
			for(size_type j = 0; j < Qdim; j++)
			{
        // Get coordinates of the node
				bgeot::base_node where = M_uFEMcur->getFEM().point_of_basic_dof(ii + j);
        // Evaluate the Dirichlet condition
        scalar_type value = M_BCcur->BCDiri(where, j, M_BCcur->getDiriBD()[M_rowsStrongBCFlagscur[i]]);

				M_Sys.setNullRow(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ));
        // Set the diagonal value to 1
        M_Sys.setMatrixValue(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), 1);

        // Set value on the rhs
				M_Sys.setRHSValue(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), value);

			}

		}

}


void BasicMethod::treatIFaceDofs(){

  size_type Qdim = M_uFEM1.getFEM().get_qdim();

  // Get index of the interface (order: bottom, right, top, left)
	size_type idxIFaceInDiriBD;
 	for (size_type j=0; j<M_BC1.getDiriBD().size(); j++){
 		if (M_BC1.getDiriBD()[j] == interfaceIdx1) {idxIFaceInDiriBD = j;}
 	}

	// For each interface node...
	for (size_type k = 0 ; k < M_nbDOFIFace ; k += Qdim){


		size_type idx1 = dof_IFace1[k];
		size_type idx2 = dof_IFace2[k] + M_nbDOF1;

    // For each dimension of the problem
		for (size_type j = 0; j < Qdim; j++)
		{
			// Get coordinates of the node
			bgeot::base_node where;
			where = M_uFEM1.getFEM().point_of_basic_dof(idx1 + j);

			// Change the matrix
			M_Sys.setNullRow(idx1 + j);
			M_Sys.setMatrixValue(idx1 + j, idx1 + j, 1);
			M_Sys.setMatrixValue(idx1 + j, idx2 + j, -1);

			// Evaluate the interface condition q0 and set the rhs value
			scalar_type value = M_BC1.BCDiri(where, j, M_BC1.getDiriBD()[idxIFaceInDiriBD]);
			M_Sys.setRHSValue(idx1 + j, value);

		}

	}
}


void BasicMethod::solve()
{
  M_Sys.solve();

  //Store the solution contained in LinearSystem in the vector M_uSol of Problem
  M_Sys.extractSubVector(M_uSol, 0, "sol");

        #ifdef DEBUG
        std::cout << "Size of the solution in Problem: "<< M_uSol.size() << std::endl;
        std::cout << "\nSolution in BasicMethod::solve() "<<std::endl;
        for (size_type k=0; k < M_uSol.size(); k++)
        {
          std::cout << M_uSol[k]<< "\t"<<std::endl;
        }
        std::cout << "\nEND Solution in BasicMethod::solve()\n"<<std::endl;
        #endif

}
