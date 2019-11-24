#include "../include/Problem.h"

Problem::Problem (const GetPot& dataFile, Bulk* bulk1, Bulk* bulk2):
			M_Bulk1(bulk1), M_Bulk2(bulk2),
			M_BC1(dataFile, "laplacian1/"), // TODO : change the dataFile accordingly and set it
			M_BC1(dataFile, "laplacian2/"), // TODO : change the dataFile accordingly and set it
			/*I_Face1(dataFile, "laplacian/"), // TODO : change the dataFile accordingly and set it
			I_Face2(dataFile, "laplacian/"), // TODO : change the dataFile accordingly and set it */
			M_uFEM1( bulk1->getMesh(), dataFile, "laplacian1/", "Sol", "bulkData1/"),// TODO : change the dataFile accordingly and set it
			M_uFEM12( bulk2->getMesh(), dataFile, "laplacian2/", "Sol", "bulkData2/"),// TODO : change the dataFile accordingly and set it
			M_CoeffFEM1( bulk1->getMesh(), dataFile, "laplacian1/", "Coeff", "bulkData1/"), // TODO : change the dataFile accordingly and set it
			M_CoeffFEM2( bulk2->getMesh(), dataFile, "laplacian2/", "Coeff", "bulkData2/"), // TODO : change the dataFile accordingly and set it
			M_intMethod1(*(bulk1->getMesh()) ), M_intMethod2(*(bulk2->getMesh()) ),
			M_Sys()
{

   M_nbDOF1 = M_uFEM1.nb_dof();
	 M_nbDOF2 = M_uFEM2.nb_dof();
	 M_nbTotDOF = M_nbDOF1 + M_nbDOF2;
	 // M_nbDOFIFace = ; TODO

	 // Per ora faccio tutto sdoppiato, anche se poi probabilmente il metodo di integrazione sarà lo stesso sui 2 domini. Per ora ho solo copiato e incollato, non so ancora come modificare il file data!
   std::string intMethod1(dataFile ( std::string("bulkData1/laplacian1/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );
	 std::string intMethod2(dataFile ( std::string("bulkData2/laplacian2/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );

   M_intMethod1.set_integration_method(bulk1->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod1) );

	 M_intMethod2.set_integration_method(bulk2->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod2) );

   M_Bulk1->getData()->setDiff(M_CoeffFEM1.getDOFpoints());  //valuta i coefficienti negli elementi della mesh
	 M_Bulk1->getData()->setDiff(M_CoeffFEM2.getDOFpoints());  //valuta i coefficienti negli elementi della mesh

	 // Con questa funzione sto assegnando le facce della frontiera e dell'interfaccia a M_bulk1, M_bulk2 e alle sue regions quindi potrei cercare di recuperare il numero di dofs associati alla region che so essere l'interfaccia. Per adesso posso dire di sapere che è quella a destra per 1 e a sinistra per 2, poi magari posso scrivere una funzione che la determini in base alla valutazione dell'insieme di livello sui nodi.
   M_BC1.setBoundaries(M_Bulk1->getMesh());
	 M_BC2.setBoundaries(M_Bulk2->getMesh());

}

void Problem::initialize()
{
         	M_uSol.reset(new scalarVector_Type (M_uFEM.nb_dof()));
         	gmm::clear(*M_uSol);
}

//restituisce un puntatore ai FEM dl displacement.

FEM* Problem::getFEM()
{
 		return &M_uFEM;
}

// dimensiona il sistema lineare per contenere il problema
void Problem::addToSys(LinearSystem* sys)
{
	M_Sys=sys;
        M_Sys->addToMatrix(M_nbTotDOF);

}

// CONTO DEI DOF
// Il valore M_nbDOFIFace è quello relativo alle funzioni che vanno raddoppiate, quindi in questo caso M_nbTotDOF è dato dalla somma M_nbDOF1 + M_nbDOF2 per le trial functions visto che suppongo che siano 2 domini separati; poi per le test il numero di dofs sarà M_nbTotDOF - M_nbDOFIFace.
size_type Problem::getNDOF(std::string variable)
{
	if (variable=="u1")
	{		return M_nbDOF1;	}

	if (variable=="u2")
	{		return M_nbDOF2;	}

	if (variable=="v") // Meaning the test functions (smaller space)
	{		return M_nbTotDOF - M_nbDOFIFace;	}

	if (variable=="u") // Meaning the global solution (trial function)
	{		return M_nbTotDOF;	}

	if (variable=="all") // Equivalente a "u"
	{		return M_nbTotDOF; }
}

// assembla la matrice
void Problem::assembleMatrix(LinearSystem* sys)
{

	sparseMatrixPtr_Type A;
	A.reset(new sparseMatrix_Type (M_uFEM.nb_dof(),M_uFEM.nb_dof()));
	gmm::clear(*A);
	stiffness( A, M_Bulk, M_uFEM, M_CoeffFEM, M_intMethod);
	M_Sys->addSubMatrix(A, 0,0);

}

//assemblo il RHS

void Problem::assembleRHS(LinearSystem* sys)
{

	scalarVectorPtr_Type  source;
	source.reset(new scalarVector_Type (M_uFEM.nb_dof()));
	gmm::clear(*source);

	bulkLoad(source, M_Bulk, M_uFEM, M_CoeffFEM, M_intMethod);

	M_Sys->addSubVector(source,0);

	scalarVectorPtr_Type  BCvec;
	BCvec.reset(new scalarVector_Type (M_uFEM.nb_dof()));

	stressRHS( BCvec, M_Bulk,  0 , &M_BC, M_uFEM, M_uFEM, M_intMethod);

	M_Sys->addSubVector(BCvec, 0);

}


// risolve il sistema e estrae la soluzione

void Problem::solve()
{
			size_type ii;
			// Indice global
       	M_Sys->solve();

       	M_uSol.reset(new scalarVector_Type (M_uFEM.nb_dof()));
	gmm::clear(*M_uSol);

       	M_Sys->extractSubVector(M_uSol, 0, "sol");

}


void Problem::extractSol(scalarVectorPtr_Type sol)
{
	M_uSol.reset(new scalarVector_Type (M_uFEM.nb_dof()));

	gmm::clear(*M_uSol);
	gmm::copy ( gmm::sub_vector (*sol, gmm::sub_interval (0,  M_uFEM.nb_dof()) ), *M_uSol);
}


void Problem::exportVtk(std::string folder, std::string what)
{
	std::cout << "export"<<std::endl;
	getfem::vtk_export exp(folder + "Solution.vtk" );
	exp.exporting( *(M_uFEM.getFEM()));
	std::vector<scalar_type> disp(M_uFEM.getFEM()->nb_dof(),0.0);
	gmm::copy(*M_uSol, disp);
	exp.write_mesh();
	exp.write_point_data( *(M_uFEM.getFEM()), disp, "u");
}



void Problem::enforceStrongBC(bool firstTime)
{

	for ( size_type bndID = 0; bndID < M_BC.getDiriBD().size(); bndID++ )
	{
		dal::bit_vector quali = M_uFEM.getFEM()->dof_on_region( M_BC.getDiriBD()[bndID]);
		size_type counter=0;
		for(dal::bv_visitor i(quali); !i.finished(); ++i)
		{
			// Sono vettori di size_type cioè di indici: contengono gli indici dei dof interessati dalle BC e gli indici delle frontiere con Dirichlet, credo
			M_rowsStrongBC.push_back(i);
			M_rowsStrongBCFlags.push_back(bndID);
			counter++;
		}

	}

		// Loop sui dofs
		for (int i=0;i<M_rowsStrongBC.size(); ++i)
		{
			size_type ii;

			// Indice globale del dof
			ii=M_rowsStrongBC[i];

			bgeot::base_node where;
			where=M_uFEM.getFEM()->point_of_basic_dof(ii);

			// Cambia la matrice direttamente! Come dovrei fare io!
			M_Sys->setNullRow(ii);
			M_Sys->setMatrixValue(ii,ii,1);

			scalar_type value= M_BC.BCDiri(where, M_BC.getDiriBD()[M_rowsStrongBCFlags[i]]);
			M_Sys->setRHSValue(ii,value);

		}



}
