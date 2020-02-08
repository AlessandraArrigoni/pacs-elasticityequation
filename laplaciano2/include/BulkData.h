#ifndef BULKDATA_H
#define BULKDATA_H

#include "Core.h"
#include "Parser.h"

//dati meccanica relativi al bulk

class BulkData
{
public:
	 BulkData ( const GetPot& dataFile,
              const std::string& section = "bulkData/",
              const std::string& sectionProblem = "nothing",
							const std::string& domainNumber = "nothing");

	scalar_type Diff(const base_node& x );

    void setDiff(std::vector<base_node> nodes);

 	inline scalar_type getDiff(size_type i)

	{
		return (*M_DiffVector)[i];  //restituisce la diffusività in uno specifico DOF
	}

  inline std::vector<scalar_type> getDiff()
	{
		return (*M_DiffVector); //restituisce la diffusività in tutti i DOF
	}
   


	bgeot::base_node bulkLoad(bgeot::base_node x);   // la forzante
	bgeot::base_node exactSolution(const base_node & x); // la soluzione esatta
	scalar_type jump(const base_node & x);

private:

    // Attributes
    std::string M_section;
    std::string M_sectionProblem;

    std::string M_diff;
    std::string M_load;
    std::string M_jump;
		std::string M_solExact;

    scalarVectorPtr_Type M_DiffVector; 

    LifeV::Parser M_parser;
};

#endif
