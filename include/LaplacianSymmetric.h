#ifndef LAPLACIANSYMMETRIC_H
#define LAPLACIANSYMMETRIC_H

#include "SymmetricMethod.h"

/*! @file LaplacianSymmetric.h
    @brief This is the class for the management of a Laplacian problem with the symmetric formulation

   @details It is a final class of the hierarchy, where the method for the assembly of the matrix has been overridden
*/

class LaplacianSymmetric final : public SymmetricMethod
{
public:

  /*! static member for the dimension of the solution and the necessary FEM  (1 because the problem is scalar) */
  static const size_type Qdim = 1;

  //! constructor
  /*! calls the constructor of its base class SymmetricMethod

  @param dataFile reference to the GetPot object that reads the input data
  @param bulk1 reference to the Bulk object for the left subdomain
  @param bulk2 reference to the Bulk object for the right subdomain
  @param extSys reference to the LinearSystem linked to the problem
  */
  LaplacianSymmetric ( GetPot const & dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys);

   /*! overridden method for the assembly of the Matrix in the case of an elliptic scalar problem */
  void assembleMatrix() override;


private:

  /*! diffusion coefficients stored as objects of the class BulkDatum */
  BulkDatum diff1, diff2;

};


#endif
