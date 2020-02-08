#include "../include/Bulk.h"

Bulk::Bulk ( const GetPot& dataFile,
                             const std::string& section,
                             const std::string& sectionDomain,
                             const std::string& sectionProblem,
			     const std::string& domainNumber
                             ) :
            M_datafile(dataFile),
            M_section ( section ),
             M_sectionDomain ( section + sectionDomain ),
            M_sectionProblem ( section + sectionProblem + domainNumber + "/"),
	          M_Data( dataFile, section, sectionProblem, domainNumber),
            // domain
	    M_meshType( dataFile ( ( M_sectionDomain + "meshType" ).data (), "GT_PK(2,1)" ) ),
	    M_meshFile( dataFile ( ( M_sectionDomain + "meshExternal" ).data (), "none" )   ),
	    M_meshFolder( dataFile ( ( M_sectionDomain + "meshFolder" ).data (), "" )   ),
      M_Nx ( dataFile ( ( M_sectionDomain + "spatialDiscretizationX" ).data (), 10 ) ),
      M_Ny ( dataFile ( ( M_sectionDomain + "spatialDiscretizationY" ).data (), 10 ) ),
	    M_Lx ( dataFile ( ( M_sectionDomain + "lengthAbscissa" ).data (), 1. ) ),
      M_Ly ( dataFile ( ( M_sectionDomain + "lengthOrdinate" ).data (), 1. ) )
      //M_origx( dataFile ( ( M_sectionDomain + "startingAbscissa" ).data(), 0. ) ),
     // M_origy( dataFile ( ( M_sectionDomain + "startingOrdinate" ).data(), 0. ) )
{
   std::cout <<  M_sectionDomain <<std::endl;
    bgeot::pgeometric_trans pgt;

    pgt =  bgeot::geometric_trans_descriptor(M_meshType);
    size_type N = pgt->dim();

    std::vector<size_type> nsubdiv(N);

    if (M_meshFile.compare("none")==0)
    {
    nsubdiv[0]= M_Nx;
    nsubdiv[1]= M_Ny;
    getfem::regular_unit_mesh(M_mesh, nsubdiv, pgt);  //creates a semi-structured mesh

    // Translate and rescale the mesh
    bgeot::base_matrix M(N,N);  // transformation matrix (scaling, shearing, etc)
    M(0,0)=M_Lx;  	        // scale the unit mesh to [LX,LY]
    M(1,1)=M_Ly;

    bgeot::base_small_vector trasl(N);
    trasl[0]= dataFile(( M_sectionDomain + "startingAbscissa" + domainNumber ).data(), 0.);
    trasl[1]= dataFile(( M_sectionDomain + "startingOrdinate" + domainNumber ).data(), 0.);

    M_mesh.transformation(M);
    M_mesh.translation(trasl);

    }
    else
    {
       std::cout << M_meshFile<<std::endl;
       getfem::import_mesh(M_meshFolder + M_meshFile, M_mesh);
    }
    M_DataPtr=&M_Data;

}

void Bulk::exportMesh(std::string nomefile)
{
   getfem::vtk_export vtkmesh(nomefile);
   vtkmesh.exporting(M_mesh);
   vtkmesh.write_mesh();

}
