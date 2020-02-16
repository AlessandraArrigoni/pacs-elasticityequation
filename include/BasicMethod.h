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

	//! constructor
	/*! calls the constructor of its base class Problem

		@param dataFile reference to the GetPot object that reads the input data
		@param problem name of the problem to solve, as contained in the input file
		@param bulk1 reference to the Bulk object for the left subdomain
		@param bulk2 reference to the Bulk object for the right subdomain
		@param dim dimension of the problem (1=scalar, 2=vectorial)
		@param extSys reference to the LinearSystem linked to the problem
	 */
	BasicMethod(GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys);

  /*! a pure virtual method for the assembly of the matrix*/
  virtual void assembleMatrix() = 0;

   /*! overridden method for the assembly of the right hand side*/
  void assembleRHS() override;

  //! overridden method for the strong imposition of Dirichlet boundary conditions
  /*!
      This methods takes one argument
      @param domainIdx index to identify the domain (1=left, 2=right)
  */
  void enforceStrongBC(size_type const domainIdx) override;

  /*! overridden method for the management of the interface degrees of freedom */
  void treatIFaceDofs() override;

  /*! overridden method for the resolution of the linear system */
  void solve() override;

  /*! destructor*/
  virtual ~BasicMethod(){};

};
#endif
