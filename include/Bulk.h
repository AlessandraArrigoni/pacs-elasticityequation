#ifndef BULK_H
#define BULK_H

#include "Core.h"
#include "Parser.h"

//gestione del dominio 2D "bulk": non contiene i dati (coefficienti, soluzione e fornzante)

class Bulk
{
public:
	 Bulk ( const GetPot& dataFile,
                   const std::string& section = "bulkData/",
                   const std::string& sectionDomain = "domain/",
									 const std::string& sectionProblem = "laplacian",
								   const std::string& domainNumber = "1");

  void exportMesh(std::string const nomefile) const;

  inline getfem::mesh& getMeshRef() // In order to modify it in the BC
	{
		return M_mesh;
  }

	inline getfem::mesh getMesh() const // In order not to modify it in Bulk
	{
		return M_mesh;
	}

  inline scalar_type Lx() const
  {
  	return M_Lx;
  }

  inline scalar_type Ly() const
  {
  	return M_Ly;
  }

	inline scalar_type nSubX() const
	{
		return M_Nx;
	}

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
    getfem::mesh M_mesh;          // the mesh

};

#endif
