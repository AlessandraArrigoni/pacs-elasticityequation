#ifndef LAPLACIANPB2_H
#define LAPLACIANPB2_H

#include "Method2.h"

class LaplacianPb2: public Method2
{
public:
  static const size_type Qdim = 1; // Dimensione della soluzione e degli spazi FEM necessari

  LaplacianPb2 ( GetPot const & dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum diff1, diff2;

};


#endif
