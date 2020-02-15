#ifndef BASICMETHOD_H
#define BASICMETHOD_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"
#include "Problem.h"

/*! @file BasicMethod.h
    @brief An abstract class to group the common features of the basic formulation.

    @details This class inherits from the "Problem" class and specializes the methods for the assembly of the matrix, the assembly of the right hand side, the imposition of the boundary and interface conditions and the function for the resolution of the linear system.
*/

class BasicMethod: public Problem
{

public:

	/*! constructor */
	BasicMethod(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys); 

  /*! a pure virtual method for the assembly of the matrix*/
  virtual void assembleMatrix() = 0;

   /*! overriden method for the assembly of the right hand side*/
  void assembleRHS() override;

  //! overriden method for the strong Dirichlet boundary conditions 
  /*!
      This methods takes one argument
      @param domainIdx index to identify if we are in the first or in the second   domain
  */
  void enforceStrongBC(size_type const domainIdx) override;

  /*! overriden method for the management of the interface degrees of freedom */
  void treatIFaceDofs() override;

  /*! overriden method for the resolution of the linear system */
  void solve() override;

  /*! destructor*/
  virtual ~BasicMethod(){}; 

};
#endif
