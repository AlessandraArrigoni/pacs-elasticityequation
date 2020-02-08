#ifndef LAPLACIANPB_H
#define LAPLACIANPB_H

#include "Problem.h"

class LaplacianPb: public Problem
{
public:
  static const size_type Qdim = 1; // Dimensione della soluzione e degli spazi FEM necessari

  LaplacianPb ( GetPot const & dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum diff1, diff2;

};


#endif
