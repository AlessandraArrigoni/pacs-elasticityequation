#ifndef PROBLEM_H
#define PROBLEM_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"

/*! @file Problem.h
    @brief This is the base abstract class. It contains all the methods andattributes that both the “symmetric” and the “basic” approach need without specialization.

    @details This class contains the functions that extract the solution and compute or print errors in the L^2 and H^1 norm. The methods which deeply change according to the PDEs system and the algebraic formulation have been defined virtual. These are the ones for the construction of the matrix, of the right-hand-sideand the treatment of the boundary and interface conditions and the resolution of the linear system.

*/

class Problem
{
public:

  /*! constructor

    @param dataFile reference to the GetPot object that reads the input data
    @param problem name of the problem to solve, as contained in the input file
    @param bulk1 reference to the Bulk object for the left subdomain
    @param bulk2 reference to the Bulk object for the right subdomain
    @param dim dimension of the problem (1=scalar, 2=vectorial)
    @param extSys reference to the LinearSystem linked to the problem
   */
  Problem(  GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys);


	/*! returns the FEM object for the solution

		@param idx index of the subdomain (1=left, 2=right)
	 */
  FEM getFEM(size_type const idx) const;

  inline LinearSystem& getSYS() const
  {
    return M_Sys;
  }

  /*! returns the number of dofs of the discrete spaces in the Petrov-Galerkin method

		@param variable string to select the discrete space ("u1" = dofs on Omega1, "u2" = dofs on Omega2, "u" = dofs in the broken space , "v" = dofs in the continuous space, "all" = total number of dofs for the problem, as "u")
	*/
  size_type getNDOF(std::string const variable = "all") const ;

  inline scalar_type getL2ERR() const
  {
    return errL2;

  }

  inline scalar_type getH1ERR() const
  {
    return errH1;
  }

  /*! pure virtual method for the assembly of the matrix */
  virtual void assembleMatrix() = 0;

  /*! pure virtual method for the assembly of the RHS */
  virtual void assembleRHS() = 0;

  /*! pure virtual method for the imposition of the Dirichlet boundary conditions */
  virtual void enforceStrongBC(size_type const domainIdx) = 0;

  /*! pure virtual method for the treatment of the interface degrees of freedom */
  virtual void treatIFaceDofs() = 0;

  /*! pure virtual method for the resolution of the linear system */
  virtual void solve() = 0;

  /*! method to extract in "destSol" the values of the variable specified by "variable", taking them from the global solution M_uSol; the default value of the std::string implies that the whole solutions M_uSol is extracted.

  @param destSol vector to store the desired solution by copy
  @param variable string to select the solution ("u1" = solution on Omega1, "u2" = solution on Omega2, any other string = solution on the whole Omega)
  */
  void extractSol(scalarVector_Type & destSol, std::string const variable = "all");

  /*! method to export the solution in ./vtk extension

  @param folder name of the folder where to store the .vtk output files
  @param what string to select the solution ("u1" = solution on Omega1, "u2" = solution on Omega2)
  */
  void exportVtk(std::string const folder = "./vtk", std::string const what = "all");

  /*! method to compute errors in L^2 and H^1 norm */
  void computeErrors();

  /*! method to print errors of the "test" in L^2 norm in "filename1" and the errors in H^1 norm in "filename2"

  @param filename1 name of the file with L^2 norm errors
  @param filename2 name of the file with H^1 norm errors
  @param test name of the test to be run
  */
  void printErrors(std::string const filename1, std::string const filename2, std::string const test );

  /*! destructor */
  virtual ~Problem(){};

  #ifdef DEBUG
  void printCoordinatesInterfaceNodes() const ;

  void printInterfaceValues() const ;
  #endif


protected:
  Bulk & M_Bulk1; // reference to the left bulk
  Bulk & M_Bulk2; // reference to the right bulk

  BC M_BC1, M_BC2; // objects of the BC class
  size_type interfaceIdx1, interfaceIdx2;

  FEM M_uFEM1, M_uFEM2;
  FEM M_CoeffFEM1, M_CoeffFEM2;
  getfem::mesh_im M_intMethod1, M_intMethod2; // integration methods

  LinearSystem& M_Sys;

  BulkDatum M_exact_sol1, M_exact_sol2, M_source1, M_source2;
  scalarVector_Type M_uSol, M_uSol1, M_uSol2;


  size_type  M_nbDOF1, M_nbDOF2, M_nbTotDOF, M_nbDOFIFace;
  sizeVector_Type dof_IFace1, dof_IFace2;

  sizeVector_Type M_rowsStrongBC1, M_rowsStrongBC2, M_rowsIFace1, M_rowsIFace2;
  sizeVector_Type M_rowsStrongBCFlags1, M_rowsStrongBCFlags2, M_rowsIFace;

  scalar_type errL2, errH1;
};





#endif
