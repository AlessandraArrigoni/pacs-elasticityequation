#ifndef LINEARELASTICITYSYMMETRIC_H
#define LINEARELASTICITYSYMMETRIC_H

#include "SymmetricMethod.h"

/*! @file LinearElasticitySymmetric.h
    @brief This is the class for the management of a linear Elasticity problem with the symmetric formulation
 
   @details It is a final class of the hierarchy, where the method for the assembly of the matrix has been overriden.
*/

class LinearElasticitySymmetric final : public SymmetricMethod
{
public:

  /*! static member for the dimension of the solution and the necessary FEM (2 because the problem is vectorial) */
  static const size_type Qdim = 2; 

  /*! constructor */
  LinearElasticitySymmetric ( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys);

  /*! overriden method for the assembly of the Matrix in the case of a linear Elasticity problem */
  void assembleMatrix() override;


private:

   /*! Lamè parameters stored as objects of the class BulkDatum */
  BulkDatum mu1, lambda1, mu2, lambda2;

};



#endif
