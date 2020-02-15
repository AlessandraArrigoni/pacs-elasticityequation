#include "../include/FEM.h"

FEM::FEM ( const getfem::mesh& mesh,
	   				const GetPot& dataFile,
	   				const std::string& problem,
	   				const std::string& variable,
	   				const std::string& section,
	 	 				const size_type qdim) :
           			M_section ( section + problem ),
           			M_femType ( ),
           			M_SpaceDim( qdim ),
           			M_FEM(mesh)
{

   M_femType = dataFile ( ( M_section+ "FEMType"+variable ).data (), "FEM_PK(2,1)" );

   getfem::pfem pf_v;
   pf_v = getfem::fem_descriptor(M_femType);

	 M_FEM.set_qdim(qdim);
   M_FEM.set_finite_element(mesh.convex_index(), pf_v);
	 // convex_index returns the list of all the valid convex elements on the mesh

   for (size_type i=0; i<M_FEM.nb_dof();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(i));
   }

			 #ifdef DEBUG
			 std::cout<<"In constructor FEM the dimension of the space is "<< M_SpaceDim <<", the type is "<<M_femType<<" and the number of dofs is "<< M_FEM.nb_dof() << std::endl;
			 #endif
}



// questo secondo costruttore non si basa sul file di input

FEM::FEM (const getfem::mesh& mesh, const std::string femType, const size_type spaceDim ):
		M_femType ( femType ),
    M_SpaceDim( spaceDim ),
    M_FEM(mesh)
{
   getfem::pfem pf_v;
   pf_v = getfem::fem_descriptor(femType);
   M_FEM.set_qdim(spaceDim);

   M_FEM.set_finite_element(mesh.convex_index(), pf_v);

   for (size_type i=0; i<M_FEM.nb_dof();++i)
   {
   	M_DOFpoints.push_back(M_FEM.point_of_basic_dof(i));
   }

			 #ifdef DEBUG
			 std::cout<<"In constructor FEM the dimension of the space is "<< M_SpaceDim <<", the type is "<<M_femType<<" and the number of dofs is "<< M_FEM.nb_dof() << std::endl;
			 #endif

}
