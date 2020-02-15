#include "../include/Problem.h"


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
  M_source2(dataFile, "bulkData/", problem, "2", "bulkload")
{

  // Set dofs number
  M_nbDOF1 = M_uFEM1.nb_dof();
  M_nbDOF2 = M_uFEM2.nb_dof();
  M_nbTotDOF = M_nbDOF1 + M_nbDOF2;

  // Set solution vectors
  M_uSol.resize(M_nbTotDOF, 0.0);
  M_uSol1.resize(M_nbDOF1, 0.0);
  M_uSol2.resize(M_nbDOF2, 0.0);


  // Set matrix dimensions in the linear system
  M_Sys.addToMatrix(M_nbTotDOF);

        #ifdef DEBUG
        std::cout<<"In the Problem constructor: total number of dofs in 1 = "<< M_nbDOF1 <<std::endl;
        std::cout<<"In the Problem constructor: total number of dofs in 2 = "<< M_nbDOF2 <<std::endl;
        std::cout<<"In the Problem constructor: total number of dofs = "<< M_nbTotDOF <<std::endl;
        std::cout<< "Matrix dimension in Problem constructor = "<<M_Sys.getMatrix()->nrows()<<" x " <<M_Sys.getMatrix()->ncols()<<std::endl;
        #endif

  // Set integration method 
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
  std::cout<<"In the Problem constructor: total number of dofs = "<< M_nbTotDOF << ", number of interface dofs = "<< M_nbDOFIFace << std::endl;
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

	if (variable=="all") // Equivalent to "u"
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
        #ifdef DEBUG
  			scalarVectorPtr_Type DIFF1 = std::make_shared<scalarVector_Type> (M_nbDOF1);
  			scalarVectorPtr_Type DIFF2 = std::make_shared<scalarVector_Type> (M_nbDOF2);

  			exactSolution(DIFF1, M_uFEM1, M_exact_sol1);  // compute te exact solution
  			exactSolution(DIFF2, M_uFEM2, M_exact_sol2);
  			#endif

	std::cout << "export "+ what <<std::endl;
	getfem::vtk_export exp(folder + "Solution_" + what + ".vtk" );
	if (what == "u1"){
		exp.exporting( M_uFEM1.getFEM());
		extractSol(M_uSol1, what);
		std::cout<< "Solution extracted! "<< std::endl;

        // DEBUG : Compute difference with the numerical solution
        #ifdef DEBUG
        for (size_type i=0; i<M_nbDOF1; i++){
           DIFF1->at(i) -= M_uSol1.at(i);
        }
        #endif

		exp.write_mesh();
		exp.write_point_data( M_uFEM1.getFEM(), M_uSol1, what);

        #ifdef DEBUG
        exp.write_point_data(M_uFEM1.getFEM(), *DIFF1, "error1");
        #endif
	}
	if (what == "u2"){
		exp.exporting( M_uFEM2.getFEM());
		extractSol(M_uSol2, what);
		std::cout<< "Solution extracted! "<< std::endl;

        // DEBUG : Compute difference with the numerical solution
        #ifdef DEBUG
        for (size_type i=0; i<M_nbDOF2; i++){
           DIFF2->at(i) -= M_uSol2.at(i);
        }
        #endif

		exp.write_mesh();
		exp.write_point_data( M_uFEM2.getFEM(), M_uSol2, what);

        #ifdef DEBUG
        exp.write_point_data(M_uFEM2.getFEM(), *DIFF2, "error2");
        #endif
	}

}

// Computes the L2 and H1 errors with respect to the exact solution and saves them in the associated variables.

void Problem::computeErrors(){

	extractSol(M_uSol1, "u1"); // Sovrascrivono quelle già estratte dall'exportVtk!
	extractSol(M_uSol2, "u2");

  scalarVectorPtr_Type DIFF1 = std::make_shared<scalarVector_Type> (M_nbDOF1);
  exactSolution(DIFF1, M_uFEM1, M_exact_sol1); // compute te exact solution

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

void Problem::printErrors(std::string const filename1, std::string const filename2, std::string const test )
{
std::ofstream f1,f2;
f1.open(filename1, std::ios_base::app);
if (M_Bulk1.nSubY()== 20) {f1<<"Test:"<<test<<"\n";}
f1<<errL2<<"\n";
if (M_Bulk1.nSubY()== 160) {f1<<"\n";}

f2.open(filename2, std::ios_base::app);
if (M_Bulk1.nSubY()== 20) {f2<<"Test:"<<test<<"\n";}
f2<<errH1<<"\n";
if (M_Bulk1.nSubY()== 160) {f2<<"\n";}
}

#ifdef DEBUG

// Prints the coordinates of the nodes on the interface following the order of the vectors (dof_IFace1, dof_IFace2) that contain them, to check that they are the same (so we don't need to check every time if the coordinates are the same, but just take the element in the same position)
void Problem::printCoordinatesInterfaceNodes() const {

	for (size_type k=0; k < M_nbDOFIFace; k++){
			bgeot::base_node whereSX = M_uFEM1.getFEM().point_of_basic_dof(dof_IFace1[k]);
			bgeot::base_node whereDX = M_uFEM2.getFEM().point_of_basic_dof(dof_IFace2[k]);

			std::cout<<" Coordinate dofs interface LEFT: "<<whereSX[0]<<" "<<whereSX[1]<<" RIGHT: "<<whereDX[0]<< " "<<whereDX[1]<<std::endl;
	}

}

// Prints the values of the two solutions on the interface, associating the ones with the same coordinates
void Problem::printInterfaceValues() const {

	std::vector<base_node> coordsSx = M_uFEM1.getDOFpoints();
	std::vector<base_node> coordsDx = M_uFEM2.getDOFpoints();
	std::cout<< " The solution on the LEFT and RIGHT interface is:" <<std::endl;

	for (size_type k=0; k < M_nbDOFIFace; k++){
		for (size_type h=0; h < M_nbDOFIFace; h++){

			if (gmm::vect_norm2(coordsSx[dof_IFace1[k]] - coordsDx[dof_IFace2[h]])<1.0e-9){
				std::cout<<"Coords LEFT: ("<<coordsSx[dof_IFace1[k]][0]<<" , "<<coordsSx[dof_IFace1[k]][1]<<")\tCoords RIGHT: ("<<coordsDx[dof_IFace2[h]][0]<<" , "<<coordsDx[dof_IFace2[h]][1]<<")"<<std::endl;
				std::cout<<"Sol LEFT: "<< M_uSol1.at(dof_IFace1[k])<<" Sol RIGHT: "<<M_uSol2.at(dof_IFace2[h])<<std::endl;
			}
		}
	}

}
#endif
