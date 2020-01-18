#ifndef BULK_H
#define BULK_H

#include "Core.h"
#include "Parser.h"
#include "BulkData.h"
#include "BC.h"

//gestione del dominio 2D "bulk"

class Bulk
{
public:
	 Bulk ( const GetPot& dataFile,
                   const std::string& section = "bulkData/",
                   const std::string& sectionDomain = "domain/",

									 const std::string& sectionProblem = "elasticity/");
  void exportMesh(std::string nomefile);

  inline getfem::mesh* getMesh()
  {
	return &M_mesh;
  }

  inline std::vector<getfem::mesh_region*>* getNeumBoundaries()
  {
	return &M_NeumBoundaries;
  }

  inline BulkData* getData()  //restituisce puntatore alla classe che riunisce i dati (letti da input file)
  {
	return M_DataPtr;
  }

  inline scalar_type Lx()
  {
  	return M_Lx;
  }

   inline scalar_type Ly()
  {
  	return M_Ly;
  }



private:

    // Attributes
    std::string M_section;
    std::string M_sectionDomain;
    std::string M_sectionProblem;
    std::string M_meshFile;
    std::string M_meshFolder;
    GetPot M_datafile;

    BulkData M_Data;
    BulkData* M_DataPtr;

    size_type M_Nx;
    size_type M_Ny;
    scalar_type M_Lx;
    scalar_type M_Ly;
		scalar_type M_origx;
		scalar_type M_origy;

    std::string M_meshType;
    getfem::mesh M_mesh;          // the mesh

    std::vector<getfem::mesh_region*> M_NeumBoundaries;

};

#endif
