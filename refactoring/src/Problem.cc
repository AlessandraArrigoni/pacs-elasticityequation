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



