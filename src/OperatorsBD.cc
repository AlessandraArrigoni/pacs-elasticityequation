#include "../include/OperatorsBD.h"


void stressRHS( scalarVectorPtr_Type V, const Bulk& medium, BC & bcRef, const FEM& FemSol, const FEM& FemDatum, const getfem::mesh_im& im)

{
    size_type Qdim = FemSol.getFEM().get_qdim(); // I need it here since the Qdim property is defined in the Problem class, not here!

    getfem::mesh_fem femSol(FemSol.getFEM());
    getfem::mesh_fem femDatum(FemDatum.getFEM());


    // Assemble in each sub region:
    for ( size_type bndID = 0; bndID < bcRef.getNeumBD().size(); bndID++ )
    {
      scalarVector_Type neumVett(femDatum.nb_dof(), 0.0);

      // Compute datum vector in each dof
      for (size_type i = 0; i < femDatum.nb_dof(); i += Qdim)
      {
        for (size_type j = 0; j < Qdim; j++){
          neumVett[i + j] = bcRef.BCNeum(femDatum.point_of_basic_dof(i), j, bcRef.getNeumBD()[bndID]);
        }
      }

      getfem::asm_source_term(*V, im, femSol, femDatum, neumVett, medium.getMesh().region(bcRef.getNeumBD()[bndID]) );

      std::cout << "operator NeumannBC on boundary "<<bcRef.getNeumBD()[bndID]<<"      [OK]" << std::endl;

    }


}
