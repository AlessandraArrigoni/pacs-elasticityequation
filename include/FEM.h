#ifndef FEM_H
#define FEM_H

#include "Core.h"

/*! @file FEM.h
    @brief This class contains all the necessary features for a generic finite element method
	
    @details This class is actually a wrapper of the already existing getfem::fem class in the library Getfem++. It incorporates only the instruments that were really needed, such as the problem's dimension, the definition of the degrees of freedom and the selection f the functional space for the trial and test functions. 

*/


class FEM
{
public:

	/*! constructor based on the input file */
	FEM (const getfem::mesh& mesh,
	      const GetPot& dataFile,
	      const std::string& problem,
	      const std::string& variable,
        const std::string& section = "bulkData/",
        const size_type qdim = 1);

	/*! constructor that is not based on the input file */
 	FEM (const getfem::mesh& mesh,
	     const std::string femType,
   		 const size_type spaceDim );

  inline size_type nb_dof() const //number of degrees of freedom
	{
		return M_FEM.nb_dof();
	};


	inline std::string type() const
	{
	 return M_femType;
	}

  inline getfem::mesh_fem getFEM() const
  {
  	return M_FEM;
  }

  inline std::vector<base_node> getDOFpoints() const //physical point corresponding to the degrees of freedom
  {
   	return M_DOFpoints;
  }


  inline base_node point_of_basic_dof(size_type i) const
  {
   	return M_FEM.point_of_basic_dof(i);
  }

private:

    std::string M_section;
    std::string M_femType;
    size_type M_SpaceDim;
    getfem::mesh_fem M_FEM; 
    std::vector<base_node> M_DOFpoints;

};

#endifo
