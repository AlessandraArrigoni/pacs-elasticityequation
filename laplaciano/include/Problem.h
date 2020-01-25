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
	 Problem ( const GetPot& dataFile, Bulk* bulk1=NULL , Bulk* bulk2=NULL);

         void initialize();

	 		 	 FEM* getFEM(size_type const idx);

				 Bulk* getBulk(size_type const idx);

         void addToSys(LinearSystem* sys);

         void assembleMatrix(LinearSystem* sys);

         void assembleRHS(LinearSystem* sys);

	 		 	 inline LinearSystem* getSys()
	 			 {
					return M_Sys;
	 		   }

         void solve();

         void updateSol();

				 // This function extracts in destSol the values of the variable specified by variable, taking them from the global solution M_uSol; the default value of the string implies that the whole solutions M_uSol is extracted.
         void extractSol(scalarVectorPtr_Type destSol, std::string variable="all");

         void exportVtk(std::string folder="./vtk", std::string what="all");

         size_type getNDOF(std::string variable="all");

         void enforceStrongBC(size_type const domainIdx);
				 void enforceInterfaceJump();

				 // ERRORS
				 // solo per adesso li metto pubblici, poi eventualmente definisco i getter
				 scalar_type errL2dx, errL2sx, errH1dx, errH1sx;

				 void computeErrors();

				 inline scalar_type const getL2Err(){
					 return errL2;
				 }

				 inline scalar_type const getH1Err(){
					 return errH1;
				 }



				 // DEBUG
				 void printInterfaceValues();
				 void printCoordinatesInterfaceNodes();

 private:

    Bulk* M_Bulk1;
		Bulk* M_Bulk2;
    BC M_BC1, M_BC2;
		size_type interfaceIdx1, interfaceIdx2;
    FEM M_uFEM1, M_uFEM2;
    FEM M_CoeffFEM1, M_CoeffFEM2;

		// Posso provare a lasciare un solo oggetto LinearSystem e una sola soluzione e vedere se riesco a cavarmela modificando la funzione assembleMatrix in modo che utilizzi il addSubMatrix o simili.
    LinearSystem* M_Sys;

    scalarVectorPtr_Type M_uSol, M_uSol1, M_uSol2;

		// Il valore M_nbDOFIFace è quello relativo alle funzioni che vanno raddoppiate, quindi in questo caso M_nbTotDOF è dato dalla somma M_nbDOF1 + M_nbDOF2 per le trial functions visto che suppongo che siano 2 domini separati; poi per le test il numero di dofs sarà M_nbTotDOF - M_nbDOFIFace.
    size_type M_nbTotDOF, M_nbDOF1, M_nbDOF2, M_nbDOFIFace;
		std::vector<size_type> dof_IFace1, dof_IFace2;

		size_type Qdim;

    getfem::mesh_im M_intMethod1, M_intMethod2;

		// Anche qui consideriamo di trattare l'interfaccia come una BC di Dirichlet per uno dei 2 domini mentre per l'altro viene lasciata "libera", nel senso che il valore di u si trova risolvendo il sistema; ciò equivale ad assegnare BC di Neumann.
    std::vector<size_type> M_rowsStrongBC1, M_rowsStrongBC2, M_rowsIFace1, M_rowsIFace2;
    std::vector<size_type> M_rowsStrongBCFlags1, M_rowsStrongBCFlags2, M_rowsIFace;

		scalar_type errL2, errH1; //, errL2sx, errL2dx, errH1sx, errH1dx;

   // mutable LifeV::Parser M_parser;
};


#endif
