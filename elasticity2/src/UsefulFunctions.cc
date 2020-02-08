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



//determina se un punto è o meno in un triangolo (2 versioni)

bool isInTriangle ( const getfem::mesh& mesh, const size_type& elementID, const base_node& node, const scalar_type& toll )
{
    const base_node x0 = node - mesh.points_of_convex(elementID) [ 0 ];
    const base_node x1 = node - mesh.points_of_convex(elementID) [ 1 ];
    const base_node x2 = node - mesh.points_of_convex(elementID) [ 2 ];

    const scalar_type A = 0.5 * ( gmm::abs(x0 [ 0 ] * x1 [ 1 ] - x0 [ 1 ] * x1 [ 0 ]) +
                                  gmm::abs(x0 [ 0 ] * x2 [ 1 ] - x0 [ 1 ] * x2 [ 0 ]) +
                                  gmm::abs(x2 [ 0 ] * x1 [ 1 ] - x2 [ 1 ] * x1 [ 0 ]) );

    if ( gmm::abs( A - mesh.convex_area_estimate(elementID) ) <= toll )
    {
        return true;
    }
    else
    {
        return false;
    }

} // isInTriangle


bool isInTriangle ( const base_node& v1,const base_node& v2,const base_node& v3,  const base_node& node, const scalar_type& toll )
{
    const base_node x0 = node - v1;
    const base_node x1 = node - v2;
    const base_node x2 = node - v3;

    const scalar_type A = 0.5 * ( gmm::abs(x0 [ 0 ] * x1 [ 1 ] - x0 [ 1 ] * x1 [ 0 ]) +
                                  gmm::abs(x0 [ 0 ] * x2 [ 1 ] - x0 [ 1 ] * x2 [ 0 ]) +
                                  gmm::abs(x2 [ 0 ] * x1 [ 1 ] - x2 [ 1 ] * x1 [ 0 ]) );

    const scalar_type area = 0.5*gmm::abs((v1[0]-v2[0])*(v1[1]-v3[1]) - (v1[1]-v2[1])*(v1[0]-v3[0]));

    if ( gmm::abs( A -area )/area <= toll )
    {
        return true;
    }
    else
    {
        return false;
    }

} // isInTriangle

bool intersectSegments(const base_node& a1, const base_node& a2, const base_node& b1, const base_node& b2, base_node& sol){

scalar_type a11(a2[0]-a1[0]);
scalar_type a12(b1[0]-b2[0]);
scalar_type a21(a2[1]-a1[1]);
scalar_type a22(b1[1]-b2[1]);
scalar_type r1(b1[0]-a1[0]);
scalar_type r2(b1[1]-a1[1]);

scalar_type det=a11*a22-a12*a21;

if (gmm::abs(det)>1.0e-12){
scalar_type t=(a22*r1-a12*r2)/det;
scalar_type s=(-a21*r1+ a11*r2)/det;

if (t>0 && t<1 && s>0 && s<1){
sol[0]=a1[0]+(a2[0]-a1[0])*t;
sol[1]=a1[1]+(a2[1]-a1[1])*t;
sol[2]=0;

return true;
}
return false;
}
else
{
return false;
}
}

bool intersectSegmentTriangle(const base_node& a1, const base_node& a2, const base_node& v1, const base_node& v2, const base_node& v3,  base_node& s1, base_node& s2){
		
scalar_type toll(1.0e-7);
if (isInTriangle (v1,v2, v3, a1, toll ) && isInTriangle (v1,v2,v3,a2, toll )){
	s1=a1;
	s2=a2;
	return true;
}

if (isInTriangle (v1,v2,v3, a1, toll ) && !isInTriangle (v1,v2,v3,a2, toll )){
       base_node altro(0,0,0);
       if (intersectSegments(v1,v2, a1, a2, altro)){
            s1=a1;
	    s2=altro;
	
	    return true;
       }
       if (intersectSegments(v2,v3, a1, a2, altro)){
            s1=a1;
	    s2=altro;
	
		return true;
       }
       if (intersectSegments(v1,v3, a1, a2, altro)){
            s1=a1;
	    s2=altro;
	
	return true;
       }

}

if (isInTriangle (v1,v2,v3, a2, toll) && !isInTriangle (v1,v2,v3,a1, toll )){
       base_node altro(0,0,0);
       if (intersectSegments(v1,v2, a1, a2, altro)){
            s1=a2;
	    s2=altro;
		return true;
       }
       if (intersectSegments(v2,v3, a1, a2, altro)){
            s1=a2;
	    s2=altro;
	return true;
       }
       if (intersectSegments(v1,v3, a1, a2, altro)){
            s1=a2;
	    s2=altro;
	return true;
       }
       
}

if (!isInTriangle (v1,v2,v3,a1, toll ) && !isInTriangle (v1,v2,v3,a2, toll )){	
	base_node altro1(0,0,0),altro2(0,0,0);
	if (intersectSegments(v1,v2, a1, a2, altro1)&&intersectSegments(v2,v3, a1, a2, altro2)){
		s1=altro1;
		s2=altro2;
		return true;
 	}

	if (intersectSegments(v1,v2, a1, a2, altro1)&&intersectSegments(v1,v3, a1, a2, altro2)){
		s1=altro1;
		s2=altro2;
		return true;
 	}

	if (intersectSegments(v1,v3, a1, a2, altro1)&&intersectSegments(v2,v3, a1, a2, altro2)){
		s1=altro1;
		s2=altro2;
		return true;
 	}
return false;
}
return false;
}



std::pair < std::string, size_type > comparaSegni ( const std::string& region, const scalarVector_Type& signs)
{
    std::pair < std::string, size_type > compara;
    compara.first = "Base";

    scalar_type flag = 0;

    for ( size_type i = 0; i < signs.size(); ++i )
    {
        if ( !( ( region[i] == '+' && signs[i] > 0 ) || ( region[i] == '-' && signs[i] < 0 ) ) )
        {
            ++flag;
            compara.first = "Extended";
            compara.second = i;
        }
    }

    if ( flag == signs.size() )
    {
        compara.first = "Extra";
    }

    return compara;

}


