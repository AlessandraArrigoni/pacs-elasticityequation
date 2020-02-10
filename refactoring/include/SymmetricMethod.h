#ifndef SYMMETRICMETHOD_H
#define SYMMETRICMETHOD_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "Problem.h"


class SymmetricMethod: public Problem
{

public:
	SymmetricMethod(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys);

  virtual void assembleMatrix() = 0;

  void assembleRHS() override;

  void enforceStrongBC(size_type const domainIdx) override;

  void treatIFaceDofs() override;

  void solve() override;

  //size_type index_inside(size_type k, std::vector<long unsigned int>& vec);

  virtual ~SymmetricMethod(){};


protected:
BulkDatum M_jump1, M_jump2;
scalarVector_Type M_q01, M_q02;

};

#endif
