#ifndef LINEARELASTICITYSYMMETRIC_H
#define LINEARELASTICITYSYMMETRIC_H

#include "SymmetricMethod.h"


class LinearElasticitySymmetric final : public SymmetricMethod
{
public:
  static const size_type Qdim = 2; // Dimensione della soluzione e degli spazi FEM necessari

  LinearElasticitySymmetric ( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum mu1, lambda1, mu2, lambda2;

};



#endif
