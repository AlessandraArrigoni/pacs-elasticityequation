#ifndef LAPLACIANSYMMETRIC_H
#define LAPLACIANSYMMETRIC_H

#include "SymmetricMethod.h"

class LaplacianSymmetric final : public SymmetricMethod
{
public:
  static const size_type Qdim = 1; // Dimensione della soluzione e degli spazi FEM necessari

  LaplacianSymmetric ( GetPot const & dataFile, Bulk & bulk1, Bulk & bulk2, LinearSystem & extSys);

  void assembleMatrix() override;


private:
  BulkDatum diff1, diff2;

};


#endif
