#ifndef METHOD2_H
#define METHOD2_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "Problem.h" 


class Method2: public Problem
{

public:
	Method2(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys);		
 
  virtual void assembleMatrix() = 0;

  void assembleRHS() override;

  void enforceStrongBC(size_type const domainIdx) override;

  void treatIFaceDofs() override;

  void solve() override;

  size_type index_inside(size_type k, std::vector<long unsigned int>& vec);
  
  virtual ~Method2(){};


protected:
BulkDatum M_jump1, M_jump2;

};

#endif
