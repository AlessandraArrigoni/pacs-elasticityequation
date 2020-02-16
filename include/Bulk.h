#ifndef BULK_H
#define BULK_H

#include "Core.h"
#include "Parser.h"

/*! @file Bulk.h
    @brief This class is for the management of a 2-dimensional domain.

    @details This class contains all the data related to the domain and its triangulation. It does not include data regarding the specific problem to solve.
*/


class Bulk
{
public:
 	/*! constructor*/
	 Bulk ( const GetPot& dataFile,
                   const std::string& section = "bulkData/",
                   const std::string& sectionDomain = "domain/",
									 const std::string& sectionProblem = "laplacian",
								   const std::string& domainNumber = "1");

  //! method to export the mesh
  /*!
      This method takes one argument as input.
	@param nomefile a string containing the name of the file where we want to save the exported mesh
   */
  void exportMesh(std::string const nomefile) const;

  inline getfem::mesh& getMeshRef() // In order to modify it in the BC
	{
		return M_mesh;
  }

	inline getfem::mesh getMesh() const // In order not to modify it in Bulk
	{
		return M_mesh;
	}

	/*! length of the domain (x axis)*/
  inline scalar_type Lx() const
  {
  	return M_Lx;
  }

	/*! length of the domain (y axis)*/
  inline scalar_type Ly() const
  {
  	return M_Ly;
  }

	/*! number of subdivisions (x axis)*/
	inline scalar_type nSubX() const
	{
		return M_Nx;
	}

	/*! number of subdivisions (y axis)*/
	inline scalar_type nSubY() const
	{
		return M_Ny;
	}

private:

    // Attributes
    std::string M_section;
    std::string M_sectionDomain;
    std::string M_sectionProblem;
    std::string M_meshFile;
    std::string M_meshFolder;
    GetPot M_datafile;

    size_type M_Nx;
    size_type M_Ny;
    scalar_type M_Lx;
    scalar_type M_Ly;

    std::string M_meshType;
	/*! variable for the mesh */
    getfem::mesh M_mesh;

};

#endif
