#include "../include/BulkData.h"

BulkData::BulkData ( const GetPot& dataFile,
                                 const std::string& section,
                                 const std::string& sectionPb) :
            M_section ( section ),
            M_sectionProblem ( M_section + sectionPb ),
	    M_lambda( dataFile ( ( M_sectionProblem + "lambda" ).data (), "1." ) ),
      M_mu( dataFile ( ( M_sectionProblem + "mu" ).data (), "1." ) ),
	    M_load( dataFile ( ( M_sectionProblem + "bulkload" ).data (), "1." ) ),
      M_solExact(dataFile( ( M_sectionProblem + "exact_sol" ).data (), "1." ) )
  {
    std::cout<<"la forzante è "<<M_load<<std::endl;
    std::cout<<"La soluzione esatta è "<<M_solExact<<std::endl;

  }


scalar_type BulkData::Lambda(const base_node& x )
	 {
		M_parser.setString ( M_lambda );
    M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
 	  return M_parser.evaluate ();
	 }

scalar_type BulkData::Mu(const base_node& x )
 	 {
 		M_parser.setString ( M_mu );
    M_parser.setVariable ( "x", x [ 0 ] );
		M_parser.setVariable ( "y", x [ 1 ] );
	  return M_parser.evaluate ();
 	 }

void BulkData::setLambda(std::vector<base_node> nodes)
	{

		M_LambdaVector.reset(new scalarVector_Type (nodes.size()));
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_LambdaVector)[i]=Lambda(nodes[i]);
			}
	}

  void BulkData::setMu(std::vector<base_node> nodes)
  	{

  		M_MuVector.reset(new scalarVector_Type (nodes.size()));
  		for (size_type i=0;i<nodes.size();++i)
  			{
  				(*M_MuVector)[i]=Mu(nodes[i]);
  			}
  	}

bgeot::base_node BulkData::bulkLoad(bgeot::base_node x)
	{
    //std::cout<<"Calcolo il termine sorgente nel punto ("<<x[0]<<" , "<<x[1]<<")"<<std::endl;

		bgeot::base_node sol(1,0); // crea un vettore 2d con valori 1 e 0

	  for ( size_type i = 0; i < 2; ++i )
    {
 	 		M_parser.setString ( M_load );
   	  M_parser.setVariable ( "x", x [ 0 ] );
    	M_parser.setVariable ( "y", x [ 1 ] );
    	sol[i] = M_parser.evaluate (i) ;
   	}

	  return sol;
	}

  // Faccio un base_node perchè mi aspetto che poi per l'elasticità avrò due valori in x e y, anche se devo capire come definire la stringa della soluzione esatta in quel caso (potrei fare un vettore di stringhe e accedere alla prima per la soluzione in x e alla seconda per la soluzione in y ma vediamo). Nel caso del laplaciano invece devo poi prendere solo il primo elemento perchè è l'unico che sto modificando!
  bgeot::base_node BulkData::exactSolution(const base_node & x)
  {
    bgeot::base_node sol(0,0);

    for (size_type i = 0; i < 1; ++i ){
      M_parser.setString(M_solExact);
      M_parser.setVariable ( "x", x [ 0 ] );
      M_parser.setVariable ( "y", x [ 1 ] );
      sol[i]= M_parser.evaluate ();
    }
    return sol;
  }
