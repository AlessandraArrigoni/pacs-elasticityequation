#ifndef METHOD1_H
#define METHOD1_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "Problem.h" 


class Method1: public Problem
{

public:
	Method1(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys);		
 
  virtual void assembleMatrix() = 0;

  void assembleRHS() override;

  void enforceStrongBC(size_type const domainIdx) override;

  void treatIFaceDofs() override;

  void solve() override;

  virtual ~Method1(){};
};


#endif
