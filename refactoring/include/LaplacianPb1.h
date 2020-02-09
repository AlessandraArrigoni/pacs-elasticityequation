#ifndef LAPLACIANPB1_H
#define LAPLACIANPB1_H

#include "Method1.h"

class LaplacianPb1: public Method1
{
public:
  static const size_type Qdim = 1; // Dimensione della soluzione e degli spazi FEM necessari

  LaplacianPb1 ( GetPot const & dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum diff1, diff2;

};


#endif
