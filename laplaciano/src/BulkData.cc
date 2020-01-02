#include "../include/BulkData.h"

BulkData::BulkData ( const GetPot& dataFile,
                                 const std::string& section,
                                 const std::string& sectionPb) :
            M_section ( section ),
            M_sectionProblem ( M_section + sectionPb ),
	    M_diff( dataFile ( ( M_sectionProblem + "diff" ).data (), "1." ) ),
	    M_load( dataFile ( ( M_sectionProblem + "bulkload" ).data (), "1." ) )
{
    std::cout<<"la forzante è "<<M_load<<std::endl;

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
