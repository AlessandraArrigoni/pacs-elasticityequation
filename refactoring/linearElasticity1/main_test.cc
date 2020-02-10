
/**
 * Petrov-Galerkin method for interface problem
 *
*/

// #include <utility>
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_config.h"
#include "getfem/getfem_assembling.h" // import assembly methods (and comp. of
#include "getfem/getfem_import.h"
#include "getfem/getfem_export.h"     // export functions (save the solution in a file)
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"
#include "gmm/gmm_MUMPS_interface.h"
#include "gmm/gmm_superlu_interface.h"

#include "../include/UsefulFunctions.h"
#include "../include/LinearElasticityBasic.h"

/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */      // System Matrix
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif


int main(int argc, char *argv[]) {

    GetPot command_line(argc, argv);

    const std::string data_file_name = command_line.follow("inputData/linearq0", 2, "-f","--file");

    GetPot dataFile(data_file_name.data());
    std::cout<< "File name : "<< data_file_name << std::endl;

    const std::string section = "";

    const std::string vtkFolder = "output_vtk/";

    // Anche qui i nomi del problema potrei leggerli da terminale per avere lo stesso main per entrambi i problemi.
    Bulk myDomainLeft(dataFile, "bulkData/", "domain/", "elasticity","1");    //creo i domini
    Bulk myDomainRight(dataFile, "bulkData/", "domain/", "elasticity","2");

    // Initial export
    myDomainLeft.exportMesh(vtkFolder+"meshLeft.vtk");
    myDomainRight.exportMesh(vtkFolder+"meshRight.vtk");

    // Create the System
    LinearSystem mySys;

    // Create the Problem and link it to the system
    LinearElasticityBasic myProblem(dataFile, myDomainLeft, myDomainRight, mySys);

    // Assemble matrix and RHS
    myProblem.assembleMatrix();
    myProblem.assembleRHS();

    // Enforce Dirichet boundary conditions
    myProblem.enforceStrongBC(1);
    myProblem.enforceStrongBC(2);
    std::cout << "Dirichlet boundary conditions      [OK]" << std::endl;

    // Enforce interface conditions
    myProblem.treatIFaceDofs();
    std::cout << "Interface conditions      [OK]" << std::endl;

    std::string outputMatrixName = "matrixElasticity_Nref" + std::to_string(int(myDomainLeft.nSubY())) + ".mm";
    std::cout<<"Stringa nome matrice output "<<outputMatrixName<<std::endl;
    mySys.saveMatrix(outputMatrixName.c_str());   //esporto la matrice per guardarla in matlab, se serve

    // Solve the problem
    myProblem.solve();
    std::cout << "Solved system      [OK]" << std::endl;

    // Export the two solutions to Paraview
    myProblem.exportVtk(vtkFolder,"u1");  // esporto la soluzione per paraview separando i due domini
    myProblem.exportVtk(vtkFolder,"u2");

    // Compute and print errors
    myProblem.computeErrors();
    myProblem.printErrors("linearElasticityBasic_errorsL2","linearElasticityBasic_errorsH1",data_file_name);

    std::cout<<"L2 error TOTAL : " << myProblem.getL2ERR() << std::endl;
    std::cout<<"H1 error TOTAL : " << myProblem.getH1ERR() << std::endl;

}
