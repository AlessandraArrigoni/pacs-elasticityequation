
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

//#include "import_neu.h"                // the .neu mesh reader

// Level Set and Xfem stuff:
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"

#include "include/UsefulFunctions.h"

#include "include/Problem.h"


/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */      // System Matrix
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector;  /* special class for small (dim < 16) vectors */
using bgeot::base_node;          /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type;        /* = double */
using bgeot::size_type;          /* = unsigned long */
using bgeot::short_type;

/* definition of some matrix/vector types. These ones are built
 * using the predefined types in Gmm++
 */
typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;

int main(int argc, char *argv[]) {

    GetPot command_line(argc, argv);

    const std::string data_file_name1 = command_line.follow("inputData/dataSX_seno", 2, "-f",    "--file");
    const std::string data_file_name2 = command_line.follow("inputData/dataDX_seno", 2, "-f",    "--file");

    GetPot dataFile1(data_file_name1.data());
    std::cout<< "File name 1: "<< data_file_name1 << std::endl;

    GetPot dataFile2(data_file_name2.data());
    std::cout<< "File name 2: " << data_file_name2 << std::endl;

    const std::string section = "";

    const std::string vtkFolder = "output_vtk/";

    Bulk myDomainLeft(dataFile1, "bulkData/", "domain/", "laplacian/");    //creo i domini
    Bulk myDomainRight(dataFile2, "bulkData/", "domain/", "laplacian/");

    myDomainLeft.exportMesh(vtkFolder+"meshLeft.vtk");       //initial export
    myDomainRight.exportMesh(vtkFolder+"meshRight.vtk");

    Problem myProblem(dataFile1, dataFile2, & myDomainLeft, & myDomainRight);   //creo il problema
//    myProblem.printCoordinatesInterfaceNodes(); // Controllo che i vettori che contengono i dof sull'interfaccia seguano lo stesso ordine-->YES

    LinearSystem mySys;   //il sistema lineare corrispondente

    myProblem.addToSys(&mySys);   // collego problema e sistema lineare

    myProblem.initialize();   // inizializzo il problema (azzero la soluzione, forse inutile se non è tempo dip.)

    myProblem.assembleMatrix(&mySys);   //assemblo matrice e termine noto

    myProblem.assembleRHS(&mySys);

    myProblem.enforceStrongBC(1);   // impongo le condizioni essenziali lavorando sulle righe della matrice
    myProblem.enforceStrongBC(2);
    std::cout << "Dirichlet boundary conditions      [OK]" << std::endl;

    myProblem.enforceInterfaceJump(); // Impongo le condizioni di interfaccia (solo di DIRICHLET per adesso, il salto sulla derivata dovrebbe essere già assegnato da assembleRHS nel caso)
    std::cout << "Dirichlet interface conditions      [OK]" << std::endl;

    mySys.saveMatrix();   //esporto la matrice per guardarla in matlab, se serve

    myProblem.solve();   // risolvo
    std::cout << "Solved system      [OK]" << std::endl;

    myProblem.exportVtk(vtkFolder,"u1");  // esporto la soluzione per paraview
    myProblem.exportVtk(vtkFolder,"u2");

    // Compute and print errors
    myProblem.computeErrors();

    std::cout<<"L2 error LEFT : "<< myProblem.errL2sx<< "\tL2 error RIGHT : "<<myProblem.errL2dx<<"\tL2 error TOTAL : "<<myProblem.getL2Err()<<std::endl;
    std::cout<<"H1 error LEFT : "<< myProblem.errH1sx<< "\tH1 error RIGHT : "<<myProblem.errH1dx<<"\tH1 error TOTAL : "<<myProblem.getH1Err()<<std::endl;

    //myProblem.printInterfaceValues(); // per stampare i valori sull'interfaccia (associando quelli sullo stesso nodo fisico)
}


// Check number of total dofs: considera solo le frontiere come regioni! A me interessa sapere se la numerazione dei dof inizia da 0 o da 1, penso da 1
/*
std::cout<<"the total number of dofs in omega1 is : "<<std::endl;
for (int u = 0; u<40; u++){
  if (myProblem.getBulk(1)->getMesh()->has_region(u)){
    std::cout<<"\nla regione ha numero : "<<u<<std::endl;
    dal::bit_vector quali = myProblem.getFEM(1)->getFEM()->dof_on_region(u);
    for( dal::bv_visitor i(quali); !i.finished(); ++i){
      std::cout<< i << " ";
    }
  }
}*/
