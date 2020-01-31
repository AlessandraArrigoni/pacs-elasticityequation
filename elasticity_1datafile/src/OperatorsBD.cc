#include "../include/OperatorsBD.h"

// Versione nuova con la funzione di getFem come per il termine sorgente, passando l'ID della frontiera da trattare
void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, BC* bcPtr,  FEM& FemSol, FEM& FemDatum, getfem::mesh_im& im)

{
    size_type Qdim = FemSol.getFEM()->get_qdim(); // I need it here since the Qdim property is defined in the Problem class, not here!


    getfem::mesh_fem femSol(*(FemSol.getFEM()));
    getfem::mesh_fem femDatum((*FemDatum.getFEM()));


    // Assemble in each sub region:
    // In modo brutale, sperando non sia un problema, definisco un nuovo vettore dei dati per ogni frontiera (contiene il valore del dato però su tutti i dof del dominio credo) e poi sommo i risultati (automaticamente con add_source_term) nel vettore di output, sperando che quando specifico la regione su cui voglio inegrare capisca da solo che deve considerare solo i valori del dato associati! Nel nostro caso non dovrebbe essere un problema visto che abbiamo solo una frontiera di Neumann in generale.
    for ( size_type bndID = 0; bndID < bcPtr->getNeumBD().size(); bndID++ )
    {
      scalarVector_Type neumVett(femDatum.nb_dof(), 0.0);

      // Compute datum vector in each dof: spero che non sia un problema, altrimenti provo poi a distinguere i nodi interni (li lascio a 0) da quelli sulla frontiera che davvero modifico. In questo caso però le BC vanno scritte separatamente nel file di dati perchè ho già definito il metodo che distingue X e Y e il parser non mi restituisce un vettore ma già degli scalari.
      for (size_type i = 0; i < femDatum.nb_dof(); i += Qdim)
      {
        for (size_type j = 0; j < Qdim; j++){
          neumVett[i + j] = bcPtr->BCNeum(femDatum.point_of_basic_dof(i), j, bcPtr->getNeumBD()[bndID]);
        }
      }

      getfem::asm_source_term(*V, im, femSol, femDatum, neumVett, medium->getMesh()->region(bcPtr->getNeumBD()[bndID]) );

      std::cout << "operator NeumannBC on boundary "<<bcPtr->getNeumBD()[bndID]<<"      [OK]" << std::endl;

    }


}
















/*

// Versione vecchia con le stringhe per costruire a mano l'integrazione, ma non funziona per le funzioni vettoriali (dovrei almeno mettere vBase)
void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, BC* bcPtr,  FEM& FemD, FEM& FemC, getfem::mesh_im& im)

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

}*/
