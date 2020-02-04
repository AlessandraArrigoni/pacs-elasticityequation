#ifndef BULKDATUM_H
#define BULKDATUM_H

#include "Core.h"
#include "Parser.h"
//#include "UsefulFunctions.h"

//Uno dei possibili dati (scalare o vettoriale associati al dominio): per adesso non metto il vettore con tutti i valori in ogni punto della griglia perchè non saprei bene come trattare il caso scalare e vettoriale: vuol dire che poi dovrò semplicemente valutare punto per punto (come già viene fatto nel codice).
// Metto entrambe le funzioni scalare e vettoriale come protected, così poi posso accedervi dalle due classi derivate che definiranno il metodo getValue(x) usando l'una o l'altra e scegliendo il valore di ritorno opportuno perchè altrimenti è un casino con il parser passare da scalare a vettoriale.


class BulkDatum
{
public:
  BulkDatum(const GetPot& dataFile,
             const std::string& section , // "bulkData/",
             const std::string& sectionProblem, // "laplacian",
             const std::string& domainNumber, //  "1",
             const std::string& datum); //  "exact_sol"

/*protected:
  scalar_type getScalarValue(const base_node& x) const;
  base_node getVectorValue(const base_node& x) const;
*/

scalar_type getValue(const base_node & x, const size_type what); // Non posso mettere const altrimenti fa casino con il parser
// Se il bulkdatum è SCALARE metterò sempre 0 come what, se è vettoriale metterò 0 o 1 a seconda della componente di cui voglio ottenere il valore.

private:
  std::string M_section;
  std::string M_sectionProblem;

  std::string M_datum;

  LifeV::Parser M_parser;


};






#endif
