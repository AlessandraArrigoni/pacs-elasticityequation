#include "../include/OperatorsBD.h"


void stressRHS( scalarVectorPtr_Type V,
               Bulk* medium, BC* bcPtr,  FEM& FemSol, FEM& FemDatum, getfem::mesh_im& im)

{

    getfem::mesh_fem femSol(*(FemSol.getFEM()));
    getfem::mesh_fem femDatum((*FemDatum.getFEM()));


    // Assemble in each sub region:
    for ( size_type bndID = 0; bndID < bcPtr->getNeumBD().size(); bndID++ )
    {
      scalarVector_Type neumVett(femDatum.nb_dof(), 0.0);

      // Compute datum vector in each dof: spero che non sia un problema, altrimenti provo poi a distinguere i nodi interni (li lascio a 0) da quelli sulla frontiera che davvero modifico.
      for (size_type i = 0; i < femDatum.nb_dof(); ++i)
      {
        neumVett[i] = bcPtr->BCNeum(femDatum.point_of_basic_dof(i), bcPtr->getNeumBD()[bndID]);
      }

      getfem::asm_source_term(*V, im, femSol, femDatum, neumVett, medium->getMesh()->region(bcPtr->getNeumBD()[bndID]) );

      std::cout << "operator NeumannBC on boundary "<< bcPtr->getNeumBD()[bndID] <<"      [OK]" << std::endl;

    }


}
