#ifndef BASICMETHOD_H
#define BASICMETHOD_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "Problem.h"


class BasicMethod: public Problem
{

public:
	BasicMethod(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys);

  virtual void assembleMatrix() = 0;

  void assembleRHS() override;

  void enforceStrongBC(size_type const domainIdx) override;

  void treatIFaceDofs() override;

  void solve() override;

  virtual ~BasicMethod(){};
};


#endif
