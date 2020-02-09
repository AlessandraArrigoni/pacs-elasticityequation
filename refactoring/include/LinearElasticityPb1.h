#ifndef LINEARELASTICITYPB1_H
#define LINEARELASTICITYPB1_H

#include "Method1.h"


class LinearElasticityPb1: public Method1
{
public:
  static const size_type Qdim = 2; // Dimensione della soluzione e degli spazi FEM necessari

  LinearElasticityPb1 ( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum mu1, lambda1, mu2, lambda2;
  
};



#endif
