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

    // Compute non costant diffusion coefficient
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


// Uso la funzione già fornita da getFEM per l'assemblaggio, e qui valuto solo i coefficienti nel punto giusto. Alla fine è un wrapper per chiamare l'altra, quindi è abbastanza inutile e potrei toglierla per semplificare il codice : valuto i coefficienti in Problem::assembleMatrix e bon, però poi devo ripensare alla struttura di tutto il codice completo, altrimenti elasticità e laplaciano risultano completamente diversi.
void linearElasticity(sparseMatrixPtr_Type M, Bulk* medium, FEM& FemSol, FEM& FemCoef, getfem::mesh_im& im)
{
  getfem::mesh_fem femSol(*(FemSol.getFEM()));
  getfem::mesh_fem femCoef(*(FemCoef.getFEM()));

  std::vector<scalar_type> lambda(femCoef.nb_dof(),0.0);
  std::vector<scalar_type> mu(femCoef.nb_dof(),0.0);

  // Compute non costant Lamè coefficients
  for (int i=0; i<femCoef.nb_dof();++i)
  {
    lambda[i] = medium->getData()->getLambda(i); // To be defined in BULKDATA!
    mu[i] = medium->getData()->getMu(i);
  }

  getfem::asm_stiffness_matrix_for_linear_elasticity(M, im, femSol, femCoef, lambda, mu);

  std::cout << "operator linearElasticity(volume)      [OK]" << std::endl;
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


void exactSolution(scalarVectorPtr_Type V, Bulk* medium, FEM& FemD)
{
  for (size_type i=0; i<FemD.nb_dof(); i++){
    V->at(i) = medium->getData()->exactSolution(FemD.point_of_basic_dof(i))[0];
  }
}
