#ifndef BC_H
#define BC_H

#include "Core.h"
#include "Parser.h"

/*!  @file BC.h
  @brief This is a class for the management of the boundary conditions

  @details This class contains the functions to evaluate the boundary conditions and is generic with respect to the PDEs system we want to solve.
   */


class BC
{
public:
 	/*! constructor */
	BC ( const GetPot& dataFile,
	      const std::string& problem,
        const std::string& section ); 

	void setBoundaries( getfem::mesh& meshRef);

		
	inline std::vector<size_type> getNeumBD() const
	{
		return M_NeumRG;
	};

	inline std::vector<size_type> getDiriBD() const
	{
		return M_DiriRG;
	};

//! method to to evaluate the Neumann boundary conditions.
  /*!
    This method takes three input parameters and returns a scalar value.

    @param x coordinates of the point where we want to evaluate the function.
    @param flag index indicating the side of the domain
    @param what index 0 if the datum is a scalar or if we want the first component of the vector; index 1 if we want to evaluate the second component.
  */
	scalar_type BCNeum(const base_node& x, const size_type what, const size_type& flag) ; 

//! method to to evaluate the Dirichlet boundary conditions.
  /*!
    This method takes three input parameters and returns a scalar value.

    @param x coordinates of the point where we want to evaluate the function.
    @param flag index indicating the side of the domain
    @param what index 0 if the datum is a scalar or if we want the first component of the vector; index 1 if we want to evaluate the second component.
  */

	scalar_type BCDiri(const base_node& x, const size_type what, const size_type& flag) ; 

private:

    std::string M_section;
    std::string M_BCstring;

     /*! string for functions defining the Neumann boundary conditions */
    std::string M_BCNeum; 
     /*! string for functions defining the Dirichlet boundary conditions */
    std::string M_BCDiri; 

    /*! vector containing the indices of the degrees of freedom on Neumann boundaries */
    std::vector<size_type>  M_NeumRG; 
    /*! vector containing the indices of the degrees of freedom on Dirichlet boundaries */
    std::vector<size_type>  M_DiriRG;

    /*! vector containing 0-Dirichlet or 1-Neumann for each edge of the domain in this order: down-right-up-left */
    std::vector<size_type> M_BC; 

    LifeV::Parser M_parser;

    int M_nBoundaries;

};

#endif
