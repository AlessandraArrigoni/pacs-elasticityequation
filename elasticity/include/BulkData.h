#ifndef BULKDATA_H
#define BULKDATA_H

#include "Core.h"
#include "Parser.h"

//dati meccanica relativi al bulk per il caso elastico

class BulkData
{
public:
	 BulkData ( const GetPot& dataFile,
                          const std::string& section = "bulkData/",
                          const std::string& sectionProblem = "elasticity/");

	scalar_type Lambda(const base_node& x );
	scalar_type Mu(const base_node& x );

  void setLambda(std::vector<base_node> nodes);
	void setMu(std::vector<base_node> nodes);

	inline scalar_type getLambda(size_type i)
	{
		return (*M_LambdaVector)[i];  //restituisce il coefficiente lambda in uno specifico DOF
	}

  inline std::vector<scalar_type> getLambda()
	{
		return (*M_LambdaVector); //restituisce il coefficiente lambda in tutti i DOF
	}

	inline scalar_type getMu(size_type i)
	{
		return (*M_MuVector)[i];  //restituisce il coefficiente mu in uno specifico DOF
	}

  inline std::vector<scalar_type> getMu()
	{
		return (*M_MuVector); //restituisce il coefficiente mu in tutti i DOF
	}

	bgeot::base_node bulkLoad(bgeot::base_node x);   // la forzante
	bgeot::base_node exactSolution(const base_node & x); // la soluzione esatta


private:

    // Attributes
    std::string M_section;
    std::string M_sectionProblem;

    std::string M_lambda;
		std::string M_mu;
    std::string M_load;
		std::string M_solExact;

    scalarVectorPtr_Type M_LambdaVector;
		scalarVectorPtr_Type M_MuVector

    LifeV::Parser M_parser;
};

#endif
