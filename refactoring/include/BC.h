#ifndef BC_H
#define BC_H

#include "Core.h"
#include "Parser.h"

//gestione delle condizioni al contorno. Generico rispetto al problema: contiene le stringhe (riferite agli scalari o ai vettori) e i dati relativi a quali frontiere sono interessate (lettura da file); poi la valutazione nei punti verrà fatta dalle 2 classi derivate come per il bulkdata.

class BC
{
public:
	BC ( const GetPot& dataFile,
	      const std::string& problem,
        const std::string& section );

	//scalar_type BCNeum(const base_node& x, const size_type& flag);  //valuta la funzione che restituisce il valore delle BC
	//scalar_type BCDiri(const base_node& x, const size_type& flag);

	void setBoundaries( getfem::mesh& meshRef);

	inline std::vector<size_type> getNeumBD() const
	{
		return M_NeumRG;
	};

	inline std::vector<size_type> getDiriBD() const
	{
		return M_DiriRG;
	};

	scalar_type BCNeum(const base_node& x, const size_type what, const size_type& flag) ; // non posso metterla come const perchè altrimenti fa casino con il parser!

	scalar_type BCDiri(const base_node& x, const size_type what, const size_type& flag) ; // non posso metterla come const perchè altrimenti fa casino con il parser!

/*
protected:
	scalar_type getScalarValue(const std::string function, const base_node& x, const size_type flag) const ;
  base_node getVectorValue(const std::string function, const base_node& x, const size_type flag) const;
*/

private:

    std::string M_section;
    std::string M_BCstring;

		// Contengono le stringhe che descrivono le funzioni che definiscono le BC
    std::string M_BCNeum;
    std::string M_BCDiri;

		// Contengono gli indici delle frontiere su cui ho definito le due diverse BC
    std::vector<size_type>  M_NeumRG;
    std::vector<size_type>  M_DiriRG;

		// Contiene 0-Dirichlet o 1-Neumann nell'ordine sotto-dx-sopra-sx
    std::vector<size_type> M_BC;

    LifeV::Parser M_parser;

    int M_nBoundaries;

};

#endif
