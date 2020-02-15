#ifndef LAPLACIANBASIC_H
#define LAPLACIANBASIC_H

#include "BasicMethod.h"

/*! @file LaplacianBasic.h
    @brief This is the class for the management of a Laplacian problem with the basic formulation
 
   @details It is a final class of the hierarchy, where the method for the assembly of the matrix has been overriden.
*/

class LaplacianBasic final : public BasicMethod
{
public:

   /*! static member for the dimension of the solution and the necessary FEM (1 because the problem is scalar) */
  static const size_type Qdim = 1; 

  /*! constructor */
  LaplacianBasic ( GetPot const & dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys);

  /*! overriden method for the assembly of the Matrix in the case of an elliptic scalar problem */
  void assembleMatrix() override;


private:

  /*! diffusion coefficients stored as objects of the class BulkDatum */
  BulkDatum diff1, diff2;

};


#endif
