#include "../include/UsefulFunctions.h"

void massLumping ( sparseMatrix_Type& matrix )
{

    const scalar_type matrixSize = gmm::mat_nrows(matrix);

    for ( size_type i = 0; i < matrixSize; ++i )
    {
        scalar_type lumped = 0;

        for ( size_type j = 0; j < matrixSize; ++j )
        {
            lumped += matrix(i, j);
            matrix(i, j) = 0;
        }
        matrix(i, i) = lumped;
    }
}

//chiamate all'exporter di getfem
void exportSolution ( const std::string& fileName,
                      const std::string& solutionName,
                      const getfem::mesh_fem& meshFEM,
                      const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(meshFEM);

    exporter.write_mesh();

    exporter.write_point_data(meshFEM, solution, solutionName.data());

}

void exportSolutionInCell ( const std::string& fileName,
                            const std::string& solutionName,
                            const getfem::mesh_fem& meshFEM,
                            const scalarVector_Type& solution )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(meshFEM);

    exporter.write_mesh();

    exporter.write_cell_data(solution, solutionName.data());

}

void exportMesh ( const std::string& fileName, const getfem::mesh& mesh )
{
    getfem::vtk_export exporter(fileName.data());

    exporter.exporting(mesh);

    exporter.write_mesh();
}

// calcola la distanza fra 2 punti
scalar_type pointDistance ( const scalar_type& x0,
                            const scalar_type& x1,
                            const scalar_type& y0,
                            const scalar_type& y1 )
{
    return std::sqrt(std::pow(x0 - x1, 2) + std::pow(y0 - y1, 2));
}

void fromBitVectorToStdVector ( dal::bit_vector& bitVector, std::vector < size_type >& stdVector )
{
        size_type i_cv = 0;

        for ( i_cv << bitVector; i_cv != size_type(-1); i_cv << bitVector )
        {
                stdVector.push_back ( i_cv );
        }

} // fromBitVectorToStdVector

char intToChar ( const size_type& integer )
{
        return static_cast<char>( integer + 97 );
} // intToChar


/*scalar_type getScalarValue(const std::string function, const base_node& x) const
{
  M_parser.setString ( function );
  M_parser.setVariable ( "x", x [ 0 ] );
  M_parser.setVariable ( "y", x [ 1 ] );

  return M_parser.evaluate ();
}


base_node getVectorValue(const std::string function, const base_node& x) const
{
  base_node value(0,0);

  for (size_type i = 0; i < 2; ++i ){
    M_parser.setString( function );
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    value[i]= M_parser.evaluate (i);
  }

  return value;
};*/
