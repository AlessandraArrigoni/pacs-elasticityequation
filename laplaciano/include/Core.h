#ifndef CORE_H
#define CORE_H

//include e typedef necessari

#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_derivatives.h>
#include <getfem/getfem_config.h>
#include <getfem/getfem_assembling.h> // import assembly methods (and comp. of
#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>     // export functions (save the solution in a file)
#include <gmm/gmm.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_MUMPS_interface.h>
#include <gmm/gmm_superlu_interface.h>

#include <getfem/bgeot_mesh.h>

#include "GetPot"
#include <string>

//#include <getfem/getfem_mesh_fem_product.h>
//#include <getfem/getfem_mesh_fem_global_function.h>

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <iomanip>

/* some Getfem++ types that we will be using */
using bgeot::base_small_vector;
/* special class for small (dim < 16) vectors */
using bgeot::base_node;
/* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type;
/* = double */
using bgeot::size_type;
/* = unsigned long */
using std::string;

/* definition of some matrix/vector types. These ones are built
 using the predefined types in Gmm++ */

typedef gmm::rsvector<scalar_type> sparseVector_Type;

typedef gmm::row_matrix<sparseVector_Type> sparseMatrix_Type;
typedef boost::shared_ptr<sparseMatrix_Type> sparseMatrixPtr_Type;

typedef std::vector<scalar_type> scalarVector_Type;
typedef boost::shared_ptr<scalarVector_Type> scalarVectorPtr_Type; // pointer a un vettore di double

typedef std::vector<size_type> sizeVector_Type;
typedef boost::shared_ptr<sizeVector_Type> sizeVectorPtr_Type;

typedef std::vector < std::string > stringContainer_Type;

namespace LifeV{
typedef scalar_type Double;
typedef size_type UInt;
typedef size_type ID;
typedef int Int;
}


#endif
