#ifndef LINEARELASTICITYPB2_H
#define LINEARELASTICITYPB2_H

#include "Method2.h"


class LinearElasticityPb2: public Method2
{
public:
  static const size_type Qdim = 2; // Dimensione della soluzione e degli spazi FEM necessari

  LinearElasticityPb2 ( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum mu1, lambda1, mu2, lambda2;
  
};



#endif
