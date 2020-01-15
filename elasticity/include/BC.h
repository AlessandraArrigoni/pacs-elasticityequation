#ifndef BC_H
#define BC_H

#include "Core.h"
#include "Parser.h"

//gestione delle condizioni al contorno. Generico rispetto al problema

class BC
{
public:
	 BC ( const GetPot& dataFile,
	      const std::string& problem,
                          const std::string& section1 );

	scalar_type BCNeum(const base_node& x, const size_type what, const size_type& flag);  //valuta la funzione che restituisce il valore delle BC
	scalar_type BCDiri(const base_node& x, const size_type what, const size_type& flag);

	void setBoundaries(getfem::mesh* meshPrt);

	std::vector<size_type> getNeumBD(std::string where="bulk");

	std::vector<size_type> getDiriBD(std::string where="bulk");

private:

    std::string M_section;

    // Attributes
    std::string M_BCstring;

		// Contengono le stringhe che descrivono le funzioni che definiscono le BC
    std::string M_BCNeumX, M_BCNeumY;
    std::string M_BCDiriX, M_BCDiriY;

		// Contengono gli indici delle frontiere su cui ho definito le due diverse BC
    std::vector<size_type>  M_NeumRG;
    std::vector<size_type>  M_DiriRG;

		// Contiene 0-Dirichlet o 1-Neumann nell'ordine sotto-dx-sopra-sx
    std::vector<size_type> M_BC;

    LifeV::Parser M_parser;

    int M_nBoundaries;
    std::vector<scalar_type> M_bdNodesX;
    std::vector<scalar_type> M_bdNodesY;

};

#endif
