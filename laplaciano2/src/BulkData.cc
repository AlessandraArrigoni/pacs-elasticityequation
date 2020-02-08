#include "../include/BulkData.h"

BulkData::BulkData ( const GetPot& dataFile,
                     const std::string& section,
                     const std::string& sectionPb,
                     const std::string& domainNumber) :
      M_section ( section ),
      M_sectionProblem ( section + sectionPb + domainNumber + "/"),
	    M_diff( dataFile ( ( M_sectionProblem + "diff" ).data (), "1." ) ),
	    M_load( dataFile ( ( M_sectionProblem + "bulkload" ).data (), "1." ) ),
      M_jump( dataFile ( ( M_sectionProblem + "qzero" ).data (), "1." ) ),
      M_solExact(dataFile( ( M_sectionProblem + "exact_sol" ).data (), "1." ) )
      
  {
    std::cout<<"la forzante è "<<M_load<<std::endl;
    std::cout<<"La soluzione esatta è "<<M_solExact<<std::endl;
    std::cout<<"il salto è "<<M_jump<<std::endl;
  }


scalar_type BulkData::Diff(const base_node& x )
	 {
		M_parser.setString ( M_diff);
    		M_parser.setVariable ( "x", x [ 0 ] );
 		M_parser.setVariable ( "y", x [ 1 ] );
 	        return M_parser.evaluate ();
	 }


void BulkData::setDiff(std::vector<base_node> nodes)
	{

		M_DiffVector.reset(new scalarVector_Type (nodes.size()));
		for (size_type i=0;i<nodes.size();++i)
			{
				(*M_DiffVector)[i]=Diff(nodes[i]);
			}
	}


bgeot::base_node BulkData::bulkLoad(bgeot::base_node x)
	{

		bgeot::base_node sol(1,0);

	        for ( size_type i = 0; i < 1; ++i )
    		{
 	 		M_parser.setString ( M_load);
   		        M_parser.setVariable ( "x", x [ 0 ] );
    	 		M_parser.setVariable ( "y", x [ 1 ] );
     	 		sol[i]=M_parser.evaluate (i) ;
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


  scalar_type BulkData::jump(const base_node & x)
  {
    scalar_type jump;   
      M_parser.setString(M_jump);
      M_parser.setVariable ( "x", x [ 0 ] );
      M_parser.setVariable ( "y", x [ 1 ] );
      jump= M_parser.evaluate ();
    
    return jump;
  }
