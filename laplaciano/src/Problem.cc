#include "../include/Problem.h"

Problem::Problem (const GetPot& dataFile, Bulk* bulk1, Bulk* bulk2):
			M_Bulk1(bulk1), M_Bulk2(bulk2),
			M_BC1(dataFile, "laplacian1/"), // change the dataFile accordingly and set it
			M_BC2(dataFile, "laplacian2/"), // change the dataFile accordingly and set it
			/*I_Face1(dataFile, "laplacian/"), // TODO : change the dataFile accordingly and set it
			I_Face2(dataFile, "laplacian/"), // TODO : change the dataFile accordingly and set it */
			interfaceIdx(),// change the dataFile accordingly and set it: devo capire come leggere il valore direttamente dal file!
			M_uFEM1( bulk1->getMesh(), dataFile, "laplacian1/", "Sol", "bulkData1/"),// change the dataFile accordingly and set it
			M_uFEM2( bulk2->getMesh(), dataFile, "laplacian2/", "Sol", "bulkData2/"),// change the dataFile accordingly and set it
			M_CoeffFEM1( bulk1->getMesh(), dataFile, "laplacian1/", "Coeff", "bulkData1/"), // change the dataFile accordingly and set it
			M_CoeffFEM2( bulk2->getMesh(), dataFile, "laplacian2/", "Coeff", "bulkData2/"), // change the dataFile accordingly and set it
			M_intMethod1(*(bulk1->getMesh()) ), M_intMethod2(*(bulk2->getMesh()) ),
			M_Sys()
{

   M_nbDOF1 = M_uFEM1.nb_dof();
	 M_nbDOF2 = M_uFEM2.nb_dof();
	 M_nbTotDOF = M_nbDOF1 + M_nbDOF2;
	 // M_nbDOFIFace = ; TODO : UTILE /*dal::bit_vector quali = M_uFEMcur.getFEM()->dof_on_region( M_BCcur.getDiriBD()[bndID]);*/



	 // Per ora faccio tutto sdoppiato, anche se poi probabilmente il metodo di integrazione sarà lo stesso sui 2 domini. Per ora ho solo copiato e incollato.
   std::string intMethod1(dataFile ( std::string("bulkData1/laplacian1/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );
	 std::string intMethod2(dataFile ( std::string("bulkData2/laplacian2/integrationMethod" ).data (), "IM_TRIANGLE(6)" ) );

   M_intMethod1.set_integration_method(bulk1->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod1) );
	 M_intMethod2.set_integration_method(bulk2->getMesh()->convex_index(),getfem::int_method_descriptor(intMethod2) );

   M_Bulk1->getData()->setDiff(M_CoeffFEM1.getDOFpoints());  //valuta i coefficienti negli elementi della mesh
	 M_Bulk2->getData()->setDiff(M_CoeffFEM2.getDOFpoints());  //valuta i coefficienti negli elementi della mesh

	 // Con questa funzione sto assegnando le facce della frontiera e dell'interfaccia a M_bulk1, M_bulk2 e alle sue regions quindi potrei cercare di recuperare il numero di dofs associati alla region che so essere l'interfaccia. Per adesso posso dire di sapere che è quella a destra per 1 e a sinistra per 2, poi magari posso scrivere una funzione che la determini in base alla valutazione dell'insieme di livello sui nodi.
   M_BC1.setBoundaries(M_Bulk1->getMesh());
	 M_BC2.setBoundaries(M_Bulk2->getMesh());

}

void Problem::initialize()
{
         	M_uSol.reset(new scalarVector_Type (M_nbTotDOF);
         	gmm::clear(*M_uSol);
}

//restituisce un puntatore ai FEM del displacement.

FEM* Problem::getFEM(size_type const idx)
{
	if (idx == 1)
 		{return &M_uFEM1;}
	if (idx == 2)
		{return &M_uFEM2;}
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

// assembla la matrice: decido di mettere prima tutti i dof relativi a Omega1 e poi quelli relativi a Omega2
void Problem::assembleMatrix(LinearSystem* sys)
{
	sparseMatrixPtr_Type A1, A2;

	A1.reset(new sparseMatrix_Type (M_nbDOF1, M_nbDOF1);
	gmm::clear(*A1);
	stiffness( A1, M_Bulk1, M_uFEM1, M_CoeffFEM1, M_intMethod1);

	A2.reset(new sparseMatrix_Type (M_nbDOF2, M_nbDOF2);
	gmm::clear(*A2);
	stiffness( A2, M_Bulk2, M_uFEM2, M_CoeffFEM2, M_intMethod2);

	M_Sys->addSubMatrix(A1, 0, 0);
	M_Sys->addSubMatrix(A2, M_nbDOF1, M_nbDOF1)

}

//assemblo il RHS
void Problem::assembleRHS(LinearSystem* sys)
{

	scalarVectorPtr_Type  source1;
	source.reset(new scalarVector_Type (M_nbDOF1);
	gmm::clear(*source1);
	bulkLoad(source1, M_Bulk1, M_uFEM1, M_CoeffFEM1, M_intMethod1);

	scalarVectorPtr_Type  source2;
	source.reset(new scalarVector_Type (M_nbDOF2);
	gmm::clear(*source2);
	bulkLoad(source2, M_Bulk2, M_uFEM2, M_CoeffFEM2, M_intMethod2);

	M_Sys->addSubVector(source1,0);
	M_sys->addSubVector(source2,M_nbDOF2);

  // To include the Neumann boundary conditions (maybe then we will need to change this part adding the jump too - or defining it in the dataFile in a proper way)
	scalarVectorPtr_Type  BCvec1;
	BCvec1.reset(new scalarVector_Type (M_nbDOF1);
	scalarVectorPtr_Type  BCvec2;
	BCvec2.reset(new scalarVector_Type (M_nbDOF2);

	stressRHS( BCvec1, M_Bulk1,  0 , &M_BC1, M_uFEM1, M_uFEM1, M_intMethod1);
	stressRHS( BCvec2, M_Bulk2,  0 , &M_BC2, M_uFEM2, M_uFEM2, M_intMethod2);

	M_Sys->addSubVector(BCvec1, 0);
	M_Sys->addSubVector(BCvec2, M_nbDOF1);

}


// risolve il sistema e estrae la soluzione
void Problem::solve()
{
			size_type ii;
			// Indice global
       	M_Sys->solve();

       	M_uSol.reset(new scalarVector_Type (M_nbTotDOF);

				gmm::clear(*M_uSol);

       	M_Sys->extractSubVector(M_uSol, 0, "sol");
				// The content of the variable M_sol contained in the LinearSystem M_Sys (actually it is all a matter of pointers...) is copied into the variable M_uSol of our class Problem.
}

// Provo a modificarla per estrarre le due soluzioni su Omega1 e Omega2 separatamente e poi anche per tenere quella globale. Definisco le due variabili M_uSol1 e M_uSol2 in questa classe Problem ed estraggo le parti da quella globale (poi eventualmente avrò bisogno di un'altro "layer" di funzioni per restituirle in output se serve)
void Problem::extractSol(scalarVectorPtr_Type destination, std::string variable="all")
{
	if (variable == "u1"){
			destination.reset(new scalarVector_Type (M_nbDOF1));
			gmm::clear(*destination);
			gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (0,  M_nbDOF1 ), *destination); //gmm::copy(source, destination)
	}
	else if (variable == "u2"){
			destination.reset(new scalarVector_Type (M_nbDOF2));
			gmm::clear(*destination);
			gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (M_nbDOF1,  M_nbTotDOF ), *destination); //gmm::copy(source, destination)
	}
	else {
			destination.reset(new scalarVector_Type (M_nbTotDOF));
			gmm::clear(*destination);
			gmm::copy ( gmm::sub_vector (*M_uSol, gmm::sub_interval (0,  M_nbTotDOF ), *destination);
	}
}

// We keep the two solutions always separated since we don't have the global FEM space for the whole solution (I guess)
void Problem::exportVtk(std::string folder, std::string what)
{
	std::cout << "export "+ what <<std::endl;
	getfem::vtk_export exp(folder + "Solution_" + what + ".vtk" );
	if (what == "u1"){
		exp.exporting( *(M_uFEM1.getFEM()));
		/* // provo a utilizzare la funzione extract che ho definito per ottenere la soluzione u1 qui direttamente; altrimenti dovrei chiamarla prima per "riempire" la variabile M_uSol1 e poi fare la copy
		std::vector<scalar_type> disp(M_nbDOF1,0.0);
		gmm::copy(*M_uSol1, disp); */
		scalarVectorPtr_Type dispPtr = new scalarVector_Type(M_nbDOF1,0.0);
		extractSol(dispPtr, what)
		exp.write_mesh();
		exp.write_point_data( *(M_uFEM1.getFEM()), *dispPtr, what);
	}
	if (what == "u2"){
		exp.exporting( *(M_uFEM2.getFEM()));
		/* // provo a utilizzare la funzione extract che ho definito per ottenere la soluzione u1 qui direttamente; altrimenti dovrei chiamarla prima per "riempire" la variabile M_uSol2 e poi fare la copy
		std::vector<scalar_type> disp(M_nbDOF2,0.0);
		gmm::copy(*M_uSol2, disp); */
		scalarVectorPtr_Type dispPtr = new scalarVector_Type(M_nbDOF2,0.0);
		extractSol(dispPtr, what)
		exp.write_mesh();
		exp.write_point_data( *(M_uFEM2.getFEM()), *dispPtr, what);
	}

	/* ORIGINAL VERSION
	exp.exporting( *(M_uFEM.getFEM()));
	std::vector<scalar_type> disp(M_uFEM.getFEM()->nb_dof(),0.0);
	gmm::copy(*M_uSol, disp);
	exp.write_mesh();
	exp.write_point_data( *(M_uFEM.getFEM()), disp, "u");*/
}


// Sets the DIRICHLET boundary conditions for the two domains in the global matrix contained in Sys. Actually it sets the rows of 0 and 1 also for the interface in the left (1) domain, since we defined the boundary condition in the Datafile as a Dirichlet one. Then we must modify these lines by adding the -1 in the columns that correspond to the basis functions on the right (2) domain.
void Problem::enforceStrongBC(bool firstTime, size_type const domainIdx)
{
	// Set local variables according to the current domain
	if (domainIdx == 1){
		FEM & M_uFEMcur = M_uFEM1;
		BC & M_BCcur = M_BC1;
		sizeVector_Type & M_rowsStrongBCcur = M_rowsStrongBC1;
		sizeVector_Type & M_rowsStrongBCFlagscur = M_rowsStrongBCFlags1;
	}
	else if (domainIdx == 2){
		FEM & M_uFEMcur = M_uFEM2;
		BC & M_BCcur = M_BC2;
		sizeVector_Type & M_rowsStrongBCcur = M_rowsStrongBC2;
		sizeVector_Type & M_rowsStrongBCFlagscur = M_rowsStrongBCFlags2;
	}

	for ( size_type bndID = 0; bndID < M_BCcur.getDiriBD().size(); bndID++ )
	{
		dal::bit_vector quali = M_uFEMcur.getFEM()->dof_on_region( M_BCcur.getDiriBD()[bndID]); // QUESTO PUÒ ESSERE UTILE PER OTTENERE I DOFS SULL'INTERFACCIA SE RIESCO A CAPIRE COME! getDiriBD restituisce un vettore che contiene gli indici (1,2,3,4) dei lati interessati da una BC di Dirichlet.
		size_type counter=0;
		for(dal::bv_visitor i(quali); !i.finished(); ++i)
		{
			// Sono vettori di size_type cioè di indici: contengono gli indici dei dof interessati dalle BC e gli indici delle frontiere con Dirichlet, credo
			M_rowsStrongBCcur.push_back(i); // aggiungo gli indici dei dofs, che però non sono necessariamente in ordine!
			M_rowsStrongBCFlagscur.push_back(bndID);
			counter++;
		}
	}

		// Loop sui dofs
		for (int i=0;i<M_rowsStrongBCcur.size(); ++i)
		{
			size_type ii;

			// Indice "locale" del dof
			ii = M_rowsStrongBCcur[i] ;

			bgeot::base_node where;
			where=M_uFEMcur.getFEM()->point_of_basic_dof(ii);

			// Cambia la matrice direttamente! Qui abbiamo bisogno di usare l'indice "globale" del dof : se consideriamo il secondo dominio, assumendo che la matrice in Sys sia ordinata con tutti i dof del primo dominio all'inizio e tutti gli altri alla fine, per recuperare gli indici globali delle funzioni interessate dalle BC devo aggiungere il numero di righe occupate dal primo blocco relativo al dominio1.
			M_Sys->setNullRow(ii + (domainIdx=1 ? 0 : M_nbDOF1 ));
			M_Sys->setMatrixValue(ii + (domainIdx=1 ? 0 : M_nbDOF1 ), ii + (domainIdx=1 ? 0 : M_nbDOF1 ), 1);

			scalar_type value= M_BCcur.BCDiri(where, M_BCcur.getDiriBD()[M_rowsStrongBCFlagscur[i]]);
			M_Sys->setRHSValue(ii + (domainIdx=1 ? 0 : M_nbDOF1 ), value);

		}
}

// To be called AFTER the enforceStrongBC since that function sets the whole row to 0 (even the one associated to the interface!)
// devo recuperare gli indici relativi all'interfaccia di entrambi i domini e poi modificare le righe della matrice relative a omega1: in realtà potrei solo aggiungere i -1 in corrispondenza delle colonne di omega2, visto che enforceStrongBC modifica già le righe corrispondenti, ma meglio rifare tutto.
void Problem::enforceInterfaceJump(){

}
