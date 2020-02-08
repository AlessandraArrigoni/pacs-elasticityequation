#ifndef LINEARELASTICITYPB_H
#define LINEARELASTICITYPB_H

#include "Problem.h"


class LinearElasticityPb: public Problem
{
public:
  static const size_type Qdim = 2; // Dimensione della soluzione e degli spazi FEM necessari

  LinearElasticityPb ( GetPot const & dataFile, Bulk  & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum mu1, lambda1, mu2, lambda2;
  
};



#endif
