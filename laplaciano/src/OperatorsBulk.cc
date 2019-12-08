#include "../include/OperatorsBulk.h"

void stiffness( sparseMatrixPtr_Type M,
               Bulk* medium, FEM& FemD, FEM& FemC, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));

    getfem::generic_assembly assem;

    assem.set("mu=data$1(#2);"
			   "t=comp(Grad(#1).Grad(#1).Base(#2));"
			   "M(#1,#1)+=  t(:,i,:,i,k).mu(k)");
     // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femD);

    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);

    std::vector<scalar_type> mu(femC.nb_dof(),0.0);

    for (int i=0; i<femC.nb_dof();++i)
    {
    	mu[i]=medium->getData()->getDiff(i);
    }

    assem.push_data(mu);

    // Set the matrices to save the evaluations
    assem.push_mat(*M);

    assem.assembly(-1);

    std::cout << "operator a(volume)      [OK]" << std::endl;

}

void bulkLoad( scalarVectorPtr_Type V,
               Bulk* medium, FEM& FemD, FEM& FemC, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));

    getfem::generic_assembly assem;

    assem.set("datax=data$1(#2);"
        "t=comp(Base(#1).Base(#2));"
        "V(#1)+=t(:,k).datax(k);");

      // Assign the mesh integration method
    assem.push_mi(im);
    // Assign the mesh finite element space
    assem.push_mf(femD);

    // Assign the mesh finite element space for the coefficients
    assem.push_mf(femC);

    scalarVector_Type datax(femC.nb_dof());
    for (size_type i=0; i<femC.nb_dof();++i)
    {
    	datax [ i ] = medium->getData()->bulkLoad(femC.point_of_basic_dof(i))[0]; // chiama il metodo bulkLoad contenuto nella classe getData che è un membro della classe Bulk
    }

    assem.push_data(datax);

    // Set the matrices to save the evaluations
    assem.push_vec(*V);
    assem.assembly(-1);

}
