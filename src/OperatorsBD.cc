#include "../include/OperatorsBD.h"

// Versione nuova con la funzione di getFem come per il termine sorgente, passando l'ID della frontiera da trattare
void stressRHS( scalarVectorPtr_Type V, const Bulk& medium, BC & bcRef, const FEM& FemSol, const FEM& FemDatum, const getfem::mesh_im& im)

{
    size_type Qdim = FemSol.getFEM().get_qdim(); // I need it here since the Qdim property is defined in the Problem class, not here!

    getfem::mesh_fem femSol(FemSol.getFEM());
    getfem::mesh_fem femDatum(FemDatum.getFEM());


    // Assemble in each sub region:
    // In modo brutale, sperando non sia un problema, definisco un nuovo vettore dei dati per ogni frontiera (contiene il valore del dato però su tutti i dof del dominio credo) e poi sommo i risultati (automaticamente con add_source_term) nel vettore di output, sperando che quando specifico la regione su cui voglio inegrare capisca da solo che deve considerare solo i valori del dato associati! Nel nostro caso non dovrebbe essere un problema visto che abbiamo solo una frontiera di Neumann in generale.
    for ( size_type bndID = 0; bndID < bcRef.getNeumBD().size(); bndID++ )
    {
      scalarVector_Type neumVett(femDatum.nb_dof(), 0.0);

      // Compute datum vector in each dof: spero che non sia un problema, altrimenti provo poi a distinguere i nodi interni (li lascio a 0) da quelli sulla frontiera che davvero modifico. In questo caso però le BC vanno scritte separatamente nel file di dati perchè ho già definito il metodo che distingue X e Y e il parser non mi restituisce un vettore ma già degli scalari.
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
