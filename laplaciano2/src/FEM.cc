#include "../include/FEM.h"

FEM::FEM ( const getfem::mesh* mesh,
	   const GetPot& dataFile,
	   const std::string& femspaces,
	   const std::string& variable,
	   const std::string& section) :
           M_section ( section + femspaces ),
           M_femType ( ),
           M_SpaceDim( ),
           M_FEM(*mesh), // the 2 means that we want a vectorial FEM space on this mesh with 2 components
           M_meshPtr(mesh)
{
	 M_femType = dataFile((M_section+ "FEMType"+variable ).data(),"FEM_PK(2,1)");
   M_SpaceDim= dataFile((M_section+ "spaceDimension" ).data (), 2 ); // qui credo che prenda il valore di default, ma va bene comunque

   getfem::pfem pf_v;
   pf_v = getfem::fem_descriptor(M_femType);

   M_FEM.set_finite_element(mesh->convex_index(), pf_v);
	 // convex_index returns the list of all the valid convex elements on the mesh

   for (size_type i=0; i<M_FEM.nb_dof();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(i));
   }

}
// questo secondo costruttore non si basa sul file di input

FEM::FEM (const getfem::mesh* mesh,
	         std::string femType,
   		 size_type spaceDim

              ):
		M_femType ( ),
          	 M_SpaceDim( ),
         	  M_FEM(*mesh, 2),
         	  M_meshPtr(mesh)
{
   getfem::pfem pf_v;
   pf_v = getfem::fem_descriptor(femType);
   M_FEM.set_qdim(spaceDim);
   M_SpaceDim=spaceDim;
   M_femType=femType;


   M_FEM.set_finite_element(mesh->convex_index(), pf_v);

   for (size_type i=0; i<M_FEM.nb_dof();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(i));
   }

}



size_type FEM::nb_dof(std::string which)
{
	return M_FEM.nb_dof();
}
