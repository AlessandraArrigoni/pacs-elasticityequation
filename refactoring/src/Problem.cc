#include "../include/Problem.h"

// NB: NON SO SE È GIUSTO PERCHÈ NEL CONSTRUCTOR DEL PROBLEM DOVREI PASSARE LA DIMENSIONE CHE PERÒ È CONTENUTA NELLE CLASSI DERIVATE! COMUNQUE È UN MEMBRO STATICO QUINDI IN TEORIA POSSO USARLO ANCHE SENZA AVERE NESSUN OGGETTO DI QUELLA CLASSE, MA NON SONO SICURA!

// Constructor
Problem::Problem( GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, size_type const dim, LinearSystem & extSys):
  M_Bulk1(bulk1), M_Bulk2(bulk2),
  M_BC1(dataFile, problem + "1/", "bulkData/"), // problem = laplacian or elasticity
  M_BC2(dataFile, problem + "2/", "bulkData/"),
  interfaceIdx1(dataFile( ("bulkData/ " + problem + "1/interfaceIdx").data(),1) ),
  interfaceIdx2(dataFile( ("bulkData/ " + problem + "2/interfaceIdx").data(),3) ),
  M_uFEM1( bulk1.getMeshRef(), dataFile, "femspaces/", "Sol", "bulkData/", dim),// the last parameter is the dimension of the space (2 = vector)
  M_uFEM2( bulk2.getMeshRef(), dataFile, "femspaces/", "Sol", "bulkData/", dim),// the last parameter is the dimesion of the space (2 = vector)
  M_CoeffFEM1(  bulk1.getMeshRef(), dataFile, "femspaces/", "Coeff", "bulkData/"), // change the dataFile accordingly and set it
  M_CoeffFEM2(  bulk2.getMeshRef(), dataFile, "femspaces/", "Coeff", "bulkData/"), // change the dataFile accordingly and set it
  M_intMethod1( bulk1.getMeshRef() ),
  M_intMethod2( bulk2.getMeshRef() ),
  M_Sys(extSys),
  M_exact_sol1(dataFile, "bulkData/", problem, "1", "exact_sol"),
  M_exact_sol2(dataFile, "bulkData/", problem, "2", "exact_sol"),
  M_source1(dataFile, "bulkData/", problem, "1", "bulkload"),
  M_source2(dataFile, "bulkData/", problem, "2", "bulkload") //,
//  M_nbDOF1(M_uFEM1.nb_dof()),
//  M_nbDOF2(M_uFEM2.nb_dof()),
//  M_nbTotDOF(M_nbDOF1 + M_nbDOF2)
{

  // Set dofs number
  M_nbDOF1 = M_uFEM1.nb_dof();
  M_nbDOF2 = M_uFEM2.nb_dof();
  M_nbTotDOF = M_nbDOF1 + M_nbDOF2;

  std::cout<<"In the Problem constructor: total number of dofs in 1 = "<< M_nbDOF1 <<std::endl;
  std::cout<<"In the Problem constructor: total number of dofs in 2 = "<< M_nbDOF2 <<std::endl;
  std::cout<<"In the Problem constructor: total number of dofs = "<< M_nbTotDOF <<std::endl;

  // Set solution vectors
  M_uSol.resize(M_nbTotDOF, 0.0);
  M_uSol1.resize(M_nbDOF1, 0.0);
  M_uSol2.resize(M_nbDOF2, 0.0);


  // Set matrix dimensions in the linear system
  M_Sys.addToMatrix(M_nbTotDOF);

  std::cout<< "Matrix dimension in Problem constructor = "<<M_Sys.getMatrix()->nrows()<<" x " <<M_Sys.getMatrix()->ncols()<<std::endl;

  // Set integration method (in realtà non sono sicura che la stringa venga letta)
  std::string intMethod(dataFile(std::string("bulkData/femspaces/integrationMethod").data(), "IM_TRIANGLE(6)" ) );

  M_intMethod1.set_integration_method(bulk1.getMesh().convex_index(), getfem::int_method_descriptor(intMethod) );
  M_intMethod2.set_integration_method(bulk2.getMesh().convex_index(), getfem::int_method_descriptor(intMethod) );

  // Set boundaries in BC
  M_BC1.setBoundaries(M_Bulk1.getMeshRef());
  M_BC2.setBoundaries(M_Bulk2.getMeshRef());

  // Set interface related quantities: Omega1
  size_type idxIFaceInDiriBD;
 	for (size_type j=0; j<M_BC1.getDiriBD().size(); j++){
 		if (M_BC1.getDiriBD()[j] == interfaceIdx1) {idxIFaceInDiriBD = j;}
 	}

	dal::bit_vector dal_dof_IFace1 = M_uFEM1.getFEM().dof_on_region(M_BC1.getDiriBD()[idxIFaceInDiriBD]);
 	for(dal::bv_visitor i(dal_dof_IFace1); !i.finished(); ++i)
 	{		 M_rowsIFace1.push_back(i);  }

  fromBitVectorToStdVector(dal_dof_IFace1, dof_IFace1);

  // Set interface related quantities: Omega2
  size_type idxIFaceInNeumBD;
 	for (size_type j=0; j<M_BC2.getNeumBD().size(); j++){
 			if (M_BC2.getNeumBD()[j] == interfaceIdx2) {idxIFaceInNeumBD = j;}
 	}

  dal::bit_vector dal_dof_IFace2 = M_uFEM2.getFEM().dof_on_region(M_BC2.getNeumBD()[idxIFaceInNeumBD]);
 	for(dal::bv_visitor i(dal_dof_IFace2); !i.finished(); ++i)
 	{		M_rowsIFace2.push_back(i);	 }

  fromBitVectorToStdVector ( dal_dof_IFace2, dof_IFace2 );

  // Set number of interface dofs: since the mesh is conforming, we can take one of the two vectors (same dimension)
  M_nbDOFIFace = dof_IFace1.size();

  // Output info
  std::cout<<"In the Problem constructor: total number of dofs = "<< M_nbTotDOF << ", number of interface dofs from domain Omega2 = "<< dof_IFace2.size() << std::endl;
}


FEM Problem::getFEM(size_type const idx) const
{
  if (idx == 1)
 		{return M_uFEM1;}
	if (idx == 2)
		{return M_uFEM2;}
}


size_type Problem::getNDOF(std::string variable) const
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

void Problem::assembleRHS()
{
  // Set the volumic source term
  scalarVectorPtr_Type source1 = std::make_shared<scalarVector_Type> (M_nbDOF1); // Quello che metto tra parentesi viene passato al constructor del tipo tra parentesi angolate, quindi credo che sia la dimensione del vettore.
	bulkLoad(source1, M_uFEM1, M_uFEM1, M_source1, M_intMethod1);

  scalarVectorPtr_Type source2 = std::make_shared<scalarVector_Type> (M_nbDOF2);
	bulkLoad(source2, M_uFEM2, M_uFEM2, M_source2, M_intMethod2);

	M_Sys.addSubVector(source1, 0);
	M_Sys.addSubVector(source2, M_nbDOF1);

	// Come per la matrice, "sposto" (sommandoli) anche i valori del rhs relativi alle funzioni definite sull'interfaccia di Omega1, perchè poi questi valori vengono persi quando impongo la condizione q0 con enforceInterfaceJump.
	for (size_type k = 0; k < dof_IFace1.size(); k++)
  {
		scalar_type newVal = (M_Sys.getRHS())->at(dof_IFace1[k]) + (M_Sys.getRHS())->at(dof_IFace2[k] + M_nbDOF1);
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


void Problem::solve()
{
  M_Sys.solve();

  M_Sys.extractSubVector(M_uSol, 0, "sol"); // Passo il vettore per indirizzo perchè la funzione vuole un pointer: in realtà è uno shared pointer: funziona lo stesso?
	std::cout << "Size of the solution in Problem: "<< M_uSol.size() << std::endl;
	// The content of the variable M_sol contained in the LinearSystem M_Sys (actually it is all a matter of pointers...) is copied into the variable M_uSol of our class Problem.
}

void Problem::extractSol(scalarVector_Type & destination, std::string const variable)
{
	if (variable == "u1")
  {
		gmm::copy ( gmm::sub_vector (M_uSol, gmm::sub_interval (0,  M_nbDOF1 )), destination); //gmm::copy(source, destination)
	}
	else if (variable == "u2")
  {
		gmm::copy ( gmm::sub_vector (M_uSol, gmm::sub_interval (M_nbDOF1,  M_nbDOF2 )), destination); //gmm::copy(source, destination)
	}
	else
  {
		gmm::copy ( gmm::sub_vector (M_uSol, gmm::sub_interval (0,  M_nbTotDOF )), destination);
	}
}


void Problem::exportVtk( std::string const folder, std::string const what)
{
	std::cout << "export "+ what <<std::endl;
	getfem::vtk_export exp(folder + "Solution_" + what + ".vtk" );
	if (what == "u1"){
		exp.exporting( M_uFEM1.getFEM());
		extractSol(M_uSol1, what);
		std::cout<< "Solution extracted! "<< std::endl;
		exp.write_mesh();
		exp.write_point_data( M_uFEM1.getFEM(), M_uSol1, what);
	}
	if (what == "u2"){
		exp.exporting( M_uFEM2.getFEM());
		extractSol(M_uSol2, what);
		std::cout<< "Solution extracted! "<< std::endl;
		exp.write_mesh();
		exp.write_point_data( M_uFEM2.getFEM(), M_uSol2, what);
	}

}

// Computes the L2 and H1 errors with respect to the exact solution and saves them in the associated variables.
// NB: per prima cosa estraggo le due sottosoluzioni, anche se magari è già stato fatto dall'exportVtk.
void Problem::computeErrors(){

	extractSol(M_uSol1, "u1"); // Sovrascrivono quelle già estratte dall'exportVtk!
	extractSol(M_uSol2, "u2");

  scalarVectorPtr_Type DIFF1 = std::make_shared<scalarVector_Type> (M_nbDOF1);
  exactSolution(DIFF1, M_uFEM1, M_exact_sol1); // calcolo la soluzione esatta

  scalarVectorPtr_Type DIFF2 = std::make_shared<scalarVector_Type> (M_nbDOF2);
	exactSolution(DIFF2, M_uFEM2, M_exact_sol2);

	// Compute difference with the numerical solution
	for (size_type i = 0; i < M_nbDOF1; i++){
		 DIFF1->at(i) -= M_uSol1.at(i);
	}

	for (size_type i = 0; i < M_nbDOF2; i++){
		 DIFF2->at(i) -= M_uSol2.at(i);
	}

	scalar_type errL2sx = getfem::asm_L2_norm(M_intMethod1, M_uFEM1.getFEM(), *DIFF1);
	scalar_type errH1sx = getfem::asm_H1_norm(M_intMethod1, M_uFEM1.getFEM(), *DIFF1);
	scalar_type errL2dx = getfem::asm_L2_norm(M_intMethod2, M_uFEM2.getFEM(), *DIFF2);
	scalar_type errH1dx = getfem::asm_H1_norm(M_intMethod2, M_uFEM2.getFEM(), *DIFF2);

	errL2 = errL2sx + errL2dx;
	errH1 = errH1sx + errH1dx;
}


void Problem::enforceStrongBC(size_type const domainIdx)
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

    size_type Qdim = M_uFEMcur->getFEM().get_qdim(); // DEVO RECUPERARE IL DATO DAL FEM PERCHÈ LA VARIABILE STATICA È NELLE CLASSI DERIVATE, MENTRE QUI SONO NELLA BASE!

		for ( size_type bndID = 0; bndID < M_BCcur->getDiriBD().size(); bndID++ )
		{
			dal::bit_vector quali = M_uFEMcur->getFEM().dof_on_region( M_BCcur->getDiriBD()[bndID]);
			size_type counter=0;
			for(dal::bv_visitor i(quali); !i.finished(); ++i)
			{
				// Sono vettori di size_type cioè di indici: contengono gli indici dei dof interessati dalle BC e gli indici delle frontiere con Dirichlet, credo
				M_rowsStrongBCcur.push_back(i); // aggiungo gli indici dei dofs, che però non sono necessariamente in ordine!
				M_rowsStrongBCFlagscur.push_back(bndID);
				counter++;
			}
		}

		for (size_type i = 0; i < M_rowsStrongBCcur.size(); i += Qdim)
		{
			size_type ii; // é l'indice "locale" del primo dof associato a un punto fisico (spero!)
			ii = M_rowsStrongBCcur[i] ;

			for(size_type j = 0; j < Qdim; j++)
			{
				bgeot::base_node where;
				where = M_uFEMcur->getFEM().point_of_basic_dof(ii + j);
				scalar_type value = M_BCcur->BCDiri(where, j, M_BCcur->getDiriBD()[M_rowsStrongBCFlagscur[i]]);

				M_Sys.setNullRow(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ));
				M_Sys.setMatrixValue(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), 1);

				M_Sys.setRHSValue(ii + j + (domainIdx==1 ? 0 : M_nbDOF1 ), value);

			} // fine loop su j (2 dofs associati a un punto fisico)

		}

}


void Problem::enforceInterfaceJump(){

  size_type Qdim = M_uFEM1.getFEM().get_qdim(); // DEVO RECUPERARE IL DATO DAL FEM PERCHÈ LA VARIABILE STATICA È NELLE CLASSI DERIVATE, MENTRE QUI SONO NELLA BASE!

	size_type idxIFaceInDiriBD;
 	for (size_type j=0; j<M_BC1.getDiriBD().size(); j++){
 		if (M_BC1.getDiriBD()[j] == interfaceIdx1) {idxIFaceInDiriBD = j;}
 	}

	// Loop on the interface dofs and change the matrix
	for (size_type k = 0 ; k < M_nbDOFIFace ; k += Qdim){

		// Global indices in the two domains: vectors dof_IFace1 and dof_IFace2 are initialized by the constructor. I hope, once again, that the 2 dofs associated to the same physical node occupy two consecutive places in the vector!
		size_type idx1 = dof_IFace1[k];
		size_type idx2 = dof_IFace2[k] + M_nbDOF1;

		for (size_type j = 0; j < Qdim; j++)
		{
			// Get the coordinates of the points
			bgeot::base_node where;
			where = M_uFEM1.getFEM().point_of_basic_dof(idx1 + j);

			// Change the matrix in M_Sys
			M_Sys.setNullRow(idx1 + j);
			M_Sys.setMatrixValue(idx1 + j, idx1 + j, 1);
			M_Sys.setMatrixValue(idx1 + j, idx2 + j, -1);

			// Change the RHS
			scalar_type value = M_BC1.BCDiri(where, j, M_BC1.getDiriBD()[idxIFaceInDiriBD]);
			M_Sys.setRHSValue(idx1 + j, value);

		} // fine loop su j (2 dofs associati a un punto fisico)

	}
}
