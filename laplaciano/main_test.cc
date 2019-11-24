
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
    const std::string data_file_name = command_line.follow("data", 2, "-f",
            "--file");

    GetPot dataFile(data_file_name.data());

    const std::string section = "";

    const std::string vtkFolder = "output_vtk/";

    Bulk myDomain(dataFile);    //creo il dominio

    myDomain.exportMesh(vtkFolder+"mesh.vtk");       //initial export

    Problem myProblem(dataFile, &myDomain);   //creo il problema

    LinearSystem mySys;   //il sistema lineare corrispondente

    myProblem.addToSys(&mySys);   // collego problema e sistema lineare

    myProblem.initialize();   // inizializzo il problema (azzero la soluzione, forse inutile se non è tempo dip.)

    myProblem.assembleMatrix(&mySys);   //assemblo matrice e termine noto

    myProblem.assembleRHS(&mySys);

    myProblem.enforceStrongBC(true);   // impongo le condizioni essenziali lavorando sulle righe della matrice

    mySys.saveMatrix();   //esporto la matrice per guardarla in matlab, se serve

    myProblem.solve();   // risolvo

    myProblem.exportVtk(vtkFolder,"all");  // esporto la soluzione per paraview

}
