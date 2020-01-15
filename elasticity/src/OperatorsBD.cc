#include "../include/OperatorsBD.h"

void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, scalar_type time, BC* bcPtr,  FEM& FemD, FEM& FemC, getfem::mesh_im& im)

{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));

    scalarVector_Type V_(femD.nb_dof());

    getfem::generic_assembly assem_surf;

    assem_surf.set("datax=data$1(#2);"
        "t=comp(Base(#1).Base(#2));"
        "V(#1)+=t(:,k).datax(k);");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(im);

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(femD);

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(femC);

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < bcPtr->getNeumBD().size(); bndID++ )
    {

	    scalarVector_Type datax(femC.nb_dof());

      for (size_type i=0; i<femC.nb_dof();++i)
	    {
	    	datax [ i ] = bcPtr->BCNeum(femC.point_of_basic_dof(i), bcPtr->getNeumBD()[bndID]);
	    }


    	    // Assign the coefficients
   	    assem_surf.push_data(datax);

  	    // Set the matrices to save the evaluations
   	    assem_surf.push_vec(V_);

        assem_surf.assembly(medium->getMesh()->region(bcPtr->getNeumBD()[bndID])); // Sets the region where we want to assembly the term (only the boundary)

        // Add the Neumann (or q1) terms we computed to the rhs of the problem (that already contains the source term)
        for ( size_type i = 0; i < femD.nb_dof(); ++i )
 	      {
       		(*V)[i]+=V_[i];
 	      }
	    gmm::clear(V_);
    }
    std::cout << "stress boundary    [OK]" << std::endl;

}
