#ifndef LAPLACIANBASIC_H
#define LAPLACIANBASIC_H

#include "BasicMethod.h"

class LaplacianBasic final : public BasicMethod
{
public:
  static const size_type Qdim = 1; // Dimensione della soluzione e degli spazi FEM necessari

  LaplacianBasic ( GetPot const & dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum diff1, diff2;

};


#endif
