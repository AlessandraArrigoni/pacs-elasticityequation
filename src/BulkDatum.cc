#include  "../include/BulkDatum.h"

BulkDatum::BulkDatum(const GetPot& dataFile, const std::string& section, const std::string& sectionProblem,
 const std::string& domainNumber, const std::string& datum):
  M_section ( section ),
  M_sectionProblem ( section + sectionProblem + domainNumber + "/"),
  M_datum( dataFile ( ( M_sectionProblem + datum ).data (), "1." ) )
  {
    #ifdef DEBUG
    std::cout<< "In constructor BulkDatum the datum is "<<M_datum<<std::endl;
    #endif
  }


scalar_type BulkDatum::getValue(const base_node & x, const size_type what)
{
  M_parser.setString( M_datum );
  M_parser.setVariable ( "x", x [ 0 ] );
  M_parser.setVariable ( "y", x [ 1 ] );

  return M_parser.evaluate(what);
}
