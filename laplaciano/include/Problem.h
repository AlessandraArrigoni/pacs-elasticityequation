//classe per la gestione del problema dell'elasticità

#ifndef PROBLEM_H
#define PROBLEM_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"

class Problem
{
public:
	 Problem ( const GetPot& dataFile, Bulk* bulk=NULL );

         void initialize();

	 FEM* getFEM();

         void addToSys(LinearSystem* sys);

         void assembleMatrix(LinearSystem* sys);

         void assembleRHS(LinearSystem* sys);

	 inline LinearSystem* getSys()
	 {
		return M_Sys;
	 }

         void solve();

         void updateSol();

         void extractSol(scalarVectorPtr_Type sol);

         void exportVtk(std::string folder="./vtk", std::string what="all");

         size_type getNDOF(std::string variable="all");

         void enforceStrongBC(bool firstTime);


 private:

    Bulk* M_Bulk1, M_Bulk2;
    BC M_BC1, M_BC2, IFace1, IFace2 ;
    FEM M_uFEM1, M_uFEM2;
    FEM M_CoeffFEM1, M_CoeffFEM2;

		// Posso provare a lasciare un solo oggetto LinearSystem e una sola soluzione e vedere se riesco a cavarmela modificando la funzione assembleMatrix in modo che utilizzi il addSubMatrix o simili.
    LinearSystem* M_Sys;

    scalarVectorPtr_Type M_uSol;

		// Il valore M_nbDOFIFace è quello relativo alle funzioni che vanno raddoppiate, quindi in questo caso M_nbTotDOF è dato dalla somma M_nbDOF1 + M_nbDOF2 per le trial functions visto che suppongo che siano 2 domini separati; poi per le test il numero di dofs sarà M_nbTotDOF - M_nbDOFIFace.
    size_type M_nbTotDOF, M_nbDOF1, M_nbDOF2, M_nbDOFIFace;

    getfem::mesh_im M_intMethod1, M_intMethod;

		// Anche qui consideriamo di trattare l'interfaccia come una BC di Dirichlet per uno dei 2 domini mentre per l'altro viene lasciata "libera", nel senso che il valore di u si trova risolvendo il sistema; ciò equivale ad assegnare BC di Neumann.
    std::vector<size_type> M_rowsStrongBC1, M_rowsStrongBC2, M_rowsIFace;
    std::vector<size_type> M_rowsStrongBCFlags1, M_rowsStrongBCFlags2, M_rowsIFace;

   // mutable LifeV::Parser M_parser;
};


#endif
