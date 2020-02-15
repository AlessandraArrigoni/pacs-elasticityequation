#ifndef BULKDATUM_H
#define BULKDATUM_H

#include "Core.h"
#include "Parser.h"

/*!  @file BulkDatum.h
  @brief This is a class for any kind of data related to the problem.

  @details This class is used to store the string describing the parameters (diffusion, Lamè...) and the functions (forcing, exact solution) related to the problem we are interested in. They can be both scalars and vectors.
   */


class BulkDatum
{
public:
 /*!constructor */
  BulkDatum(const GetPot& dataFile,
             const std::string& section , // "bulkData/",
             const std::string& sectionProblem, // "laplacian",
             const std::string& domainNumber, //  "1",
             const std::string& datum); //  "exact_sol"


//! method to evaluate the string.
  /*!
    This method takes two input parameters and returns a scalar value.

    @param x coordinates of the point where we want to evaluate the function.
    @param what index 0 if the datum is a scalar or if we want the first component of the vector; index 1 if we want to evaluate the second component.
  */
scalar_type getValue(const base_node & x, const size_type what); 

private:
  std::string M_section;
  std::string M_sectionProblem;

  std::string M_datum;

  LifeV::Parser M_parser;


};






#endif
