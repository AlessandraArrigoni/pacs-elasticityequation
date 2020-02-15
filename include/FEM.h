#ifndef FEM_H
#define FEM_H

#include "Core.h"


// Classe wrapper per la classe fem di getFEM

class FEM
{
public:

	FEM (const getfem::mesh& mesh,
	      const GetPot& dataFile,
	      const std::string& problem,
	      const std::string& variable,
        const std::string& section = "bulkData/",
        const size_type qdim = 1);

 	FEM (const getfem::mesh& mesh,
	     const std::string femType,
   		 const size_type spaceDim );

  inline size_type nb_dof() const //numero gradi di libertà
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

  inline std::vector<base_node> getDOFpoints() const //punti corrispondenti ai dof
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

#endif
