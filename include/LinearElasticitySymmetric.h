#ifndef LINEARELASTICITYSYMMETRIC_H
#define LINEARELASTICITYSYMMETRIC_H

#include "SymmetricMethod.h"

/*! @file LinearElasticitySymmetric.h
    @brief This is the class for the management of a linear Elasticity problem with the symmetric formulation

   @details It is a final class of the hierarchy, where the method for the assembly of the matrix has been overridden.
*/

class LinearElasticitySymmetric final : public SymmetricMethod
{
public:

  /*! static member for the dimension of the solution and the necessary FEM (2 because the problem is vectorial) */
  static const size_type Qdim = 2;

  //! constructor
  /*! calls the constructor of its base class SymmetricMethod

  @param dataFile reference to the GetPot object that reads the input data
  @param bulk1 reference to the Bulk object for the left subdomain
  @param bulk2 reference to the Bulk object for the right subdomain
  @param extSys reference to the LinearSystem linked to the problem
  */
  LinearElasticitySymmetric ( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys);

  /*! overridden method for the assembly of the Matrix in the case of the linear elasticity problem */
  void assembleMatrix() override;


private:

   /*! Lamé parameters stored as objects of the class BulkDatum */
  BulkDatum mu1, lambda1, mu2, lambda2;

};



#endif
