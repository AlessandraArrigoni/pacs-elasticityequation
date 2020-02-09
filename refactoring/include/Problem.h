#ifndef PROBLEM_H
#define PROBLEM_H

#include "LinearSystem.h"
#include "Operators.h"
#include "UsefulFunctions.h"
#include "StringUtility.h"

// Copio da Problem della Scotti ma invece dei puntatori metto le references!

class Problem
{
public:
  Problem(  GetPot const & dataFile, std::string const problem, Bulk & bulk1, Bulk & bulk2, const size_type dim, LinearSystem & extSys);

  //void initialize(); // Ho messo questa parte nel constructor

  FEM getFEM(size_type const idx) const;

  //const Bulk& getBULK(size_type const idx) const; //solo metodi const possono trattare variabili const;

  inline LinearSystem& getSYS() const // non sono molto sicura di questa funzione: posso restituire una reference così?
  {
    return M_Sys;
  }

  size_type getNDOF(std::string const variable = "all") const ;

  inline scalar_type getL2ERR() const
  {
    return errL2;
  }

  inline scalar_type getH1ERR() const
  {
    return errH1;
  }

  //void addToSys( LinearSystem & extSys ); // Già fatto nel constructor

  virtual void assembleMatrix() = 0;

  virtual void assembleRHS() = 0;

  virtual void enforceStrongBC(size_type const domainIdx) = 0;

  virtual void treatIFaceDofs() = 0;

  virtual void solve() = 0;

  // This function extracts in destSol the values of the variable specified by variable, taking them from the global solution M_uSol; the default value of the std::string implies that the whole solutions M_uSol is extracted.
  void extractSol(scalarVector_Type & destSol, std::string const variable = "all");

  void exportVtk(std::string const folder = "./vtk", std::string const what = "all"); //non la metto come const perchè chiama extractSol che modifica la classe!

  void computeErrors();

  virtual ~Problem(){};



protected:
  Bulk & M_Bulk1;
  Bulk & M_Bulk2; // Qui ho delle REFERENCES perchè i due domini sono creati nel main prima del problem, quindi gli oggetti esistono già, invece le BC e gli spazi FEM sono creati allo stesso momento del Problem, direttamente come membri di questa classe. Ho definito questi Bulk come const perchè il problem non li deve modificare.
  BC M_BC1, M_BC2;
  size_type interfaceIdx1, interfaceIdx2;
  FEM M_uFEM1, M_uFEM2;
  FEM M_CoeffFEM1, M_CoeffFEM2;
  getfem::mesh_im M_intMethod1, M_intMethod2;

  LinearSystem& M_Sys;

  BulkDatum M_exact_sol1, M_exact_sol2, M_source1, M_source2;
  scalarVector_Type M_uSol, M_uSol1, M_uSol2;

  // Il valore M_nbDOFIFace è quello relativo alle funzioni che vanno raddoppiate, quindi in questo caso M_nbTotDOF è dato dalla somma M_nbDOF1 + M_nbDOF2 per le trial functions visto che suppongo che siano 2 domini separati; poi per le test il numero di dofs sarà M_nbTotDOF - M_nbDOFIFace.
  size_type  M_nbDOF1, M_nbDOF2, M_nbTotDOF, M_nbDOFIFace;
  sizeVector_Type dof_IFace1, dof_IFace2; // sono inizializzati dal constructor!

  sizeVector_Type M_rowsStrongBC1, M_rowsStrongBC2, M_rowsIFace1, M_rowsIFace2;
  sizeVector_Type M_rowsStrongBCFlags1, M_rowsStrongBCFlags2, M_rowsIFace;

  scalar_type errL2, errH1;
};





#endif
