#ifndef LINEARELASTICITYBASIC_H
#define LINEARELASTICITYBASIC_H

#include "BasicMethod.h"


class LinearElasticityBasic final : public BasicMethod
{
public:
  static const size_type Qdim = 2; // Dimensione della soluzione e degli spazi FEM necessari

  LinearElasticityBasic ( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum mu1, lambda1, mu2, lambda2;

};



#endif
