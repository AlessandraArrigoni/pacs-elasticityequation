#include "../include/OperatorsBulk.h"


// Uso la funzione già fornita da getFEM per l'assemblaggio, e qui valuto solo i coefficienti nel punto giusto. Alla fine è un wrapper per chiamare l'altra, quindi è abbastanza inutile e potrei toglierla per semplificare il codice : valuto i coefficienti in Problem::assembleMatrix e bon, però poi devo ripensare alla struttura di tutto il codice completo, altrimenti elasticità e laplaciano risultano completamente diversi.
void linearElasticity(sparseMatrixPtr_Type M, Bulk* medium, FEM& FemSol, FEM& FemCoef, getfem::mesh_im& im)
{
  getfem::mesh_fem femSol(*(FemSol.getFEM()));
  getfem::mesh_fem femCoef(*(FemCoef.getFEM()));

  scalarVector_Type lambda(femCoef.nb_dof(),0.0);
  scalarVector_Type mu(femCoef.nb_dof(),0.0);

  // Compute non costant Lamè coefficients
  for (size_type i=0; i<femCoef.nb_dof();++i)
  {
    lambda[i] = medium->getData()->getLambda(i); // Defined in BULKDATA!
    mu[i] = medium->getData()->getMu(i);
  }

            // DEBUG:
            std::cout<<"The dimension of the mesh is "<<medium->getMesh()->dim()<<std::endl;
            std::cout<<"The dimension of the space is "<<FemSol.getFEM()->get_qdim()<<std::endl;

  getfem::asm_stiffness_matrix_for_linear_elasticity(*M, im, femSol, femCoef, lambda, mu);

  std::cout << "operator linearElasticity(volume)      [OK]" << std::endl;
}

// Metodo con la funzione della libreria a cui passo il vettore di dati con Fx e Fy associate allo stesso nodo fisico una in seguito all'altra
// getfem::asm_source_term(B, mim, mfu, mfd, V);
// NB: Here there isn't any femCoef, but we specify the FEM space where the datum (source) is defined; in general it is the same as the FEM space where the solution is defined.
void bulkLoad(scalarVectorPtr_Type V,
               Bulk* medium, FEM& FemSol, FEM& FemSource, getfem::mesh_im& im)
{
  size_type Qdim = FemSol.getFEM()->get_qdim(); // I need it here since the Qdim property is defined in the Problem class, not here!

  getfem::mesh_fem femSol(*(FemSol.getFEM()));
  getfem::mesh_fem femSource(*(FemSource.getFEM()));

  scalarVector_Type sourceVett(femSource.nb_dof(), 0.0);
      //DEBUG
      //std::cout << "in PROBLEM::ASSEMBLERHS::BULKLOAD, dimensione vettore source: "<<femSource.nb_dof()<<" dimensione spazio soluzione: "<<femSol.nb_dof()<<" QDIM = "<<Qdim<<std::endl;

  // Compute source term in each dof
  for (size_type i = 0; i < femSource.nb_dof(); i += Qdim)
  {
    for (size_type j = 0; j < Qdim; j++){
      sourceVett[i + j] = medium->getData()->bulkLoad(femSource.point_of_basic_dof(i))[j]; // chiama il metodo bulkLoad contenuto nella classe getData che è un membro della classe Bulk: restituisce un vettore con i due valori, se il parser fa quello che dico io.
    }
  }

  getfem::asm_source_term(*V, im, femSol, femSource, sourceVett);

  std::cout << "operator bulkLoad(volume)      [OK]" << std::endl;
}



void exactSolution(scalarVectorPtr_Type V, Bulk* medium, FEM& FemSol)
{
  size_type Qdim = FemSol.getFEM()->get_qdim(); // I need it here since the Qdim property is defined in the Problem class, not here!

  for (size_type i=0; i<FemSol.nb_dof(); i+= Qdim)
  {
    for (size_type j = 0; j < Qdim; j++){
      V->at(i+j) = medium->getData()->exactSolution(FemSol.point_of_basic_dof(i))[j];

          //DEBUG
          //std::cout<<"punto: ("<<FemSol.point_of_basic_dof(i)[0]<<" , "<<FemSol.point_of_basic_dof(i)[1]<<") e valore della componente "<<j<<": "<< V->at(i+j)<<std::endl;
    }
  }
}

void jump(scalarVectorPtr_Type V, Bulk* medium, FEM& FemSol)
{
  size_type Qdim = FemSol.getFEM()->get_qdim(); 

  for (size_type i=0; i<FemSol.nb_dof(); i+= Qdim)
  {
    for (size_type j = 0; j < Qdim; j++){
      V->at(i+j) = medium->getData()->jump(FemSol.point_of_basic_dof(i))[j];

    }
  }
}




/////////////////////////////////////////////////////////////////////
// VERSIONI CON LE STRINGHE PASSATE A MANO, NON VALIDE PER VETTORI //
/////////////////////////////////////////////////////////////////////


/*
void bulkLoad( scalarVectorPtr_Type V,
               Bulk* medium, FEM& FemD, FEM& FemC, getfem::mesh_im& im)
{
    getfem::mesh_fem femD(*(FemD.getFEM()));
    getfem::mesh_fem femC((*FemC.getFEM()));

    getfem::generic_assembly assem;

    // Versione scalare per il Laplaciano
    //assem.set("datax=data$1(#2);"
    //    "t=comp(Base(#1).Base(#2));"
    //    "V(#1)+=t(:,k).datax(k);");

    assem.set("datax=data$1(#2);"
        "t=comp(vBase(#1).vBase(#2));" // uso vBase per dire che fanno parte di uno spazio FEM vettoriale!
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

} */

/*
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

}*/
