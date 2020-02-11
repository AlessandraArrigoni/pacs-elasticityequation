#include "../include/BC.h"

BC::BC ( const GetPot& dataFile, const std::string& problem, const std::string& section) :
      M_section ( section + problem ),
	    M_nBoundaries( dataFile ( ( M_section + "nBoundaries" ).data (), 4) ),
	    M_BCstring( dataFile ( ( M_section + "bcflag" ).data (), "1") ),
 	    M_BCNeum( dataFile ( ( M_section + "du_BC" ).data (), "1") ),
 	    M_BCDiri( dataFile ( ( M_section + "u_BC" ).data (), "1") )
{
  M_BC.resize(M_nBoundaries,0);

  M_parser.setString ( M_BCstring );

  for ( size_type i = 0; i < M_nBoundaries; ++i )
  {
    M_BC [ i ] = M_parser.evaluate ( i ); // vettore di indici: suppongo che contenga 0 per Dirichlet e 1 per Neumann nell'ordine dato dal Datafile (sotto, dx, sopra, sx)

    if (M_BC[i]==0)
    {
    	M_DiriRG.push_back(i);
    }
    if (M_BC[i]==1)
    {
    	M_NeumRG.push_back(i);
    }
  }

        // Print values of the boundary conditions
        #ifdef DEBUG

        std::cout << "In constructor BC, string section: " << M_section << std::endl;
      	std::cout << "string BC: " << M_BCstring << std::endl;
      	std::cout << "string Dirichlet condition: "<< M_BCDiri << std::endl;

        std::cout<<"\nDirichlet boundaries"<<std::endl;
        for (size_type k = 0; k<M_DiriRG.size(); k++)
        {
          std::cout << M_DiriRG[k]<<"\t";
        }
        std::cout<<"\nNeumann boundaries"<<std::endl;
        for (size_type k = 0; k<M_NeumRG.size(); k++)
        {
          std::cout << M_NeumRG[k]<<"\t";
        }
        #endif
}

scalar_type BC::BCNeum(const base_node & x, const size_type what, const size_type & flag)
{
	  M_parser.setString( M_BCNeum );
    M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "n", flag );

    return M_parser.evaluate (what);
}

scalar_type BC::BCDiri(const base_node& x, const size_type what, const size_type& flag)
{
		M_parser.setString ( M_BCDiri );
	  M_parser.setVariable ( "x", x [ 0 ] );
    M_parser.setVariable ( "y", x [ 1 ] );
    M_parser.setVariable ( "n", flag );

    return M_parser.evaluate (what);
}


void BC::setBoundaries( getfem::mesh &  meshRef)
{
    getfem::mesh_region border_faces;
    getfem::outer_faces_of_mesh(meshRef, border_faces);

    for (getfem::mr_visitor i(border_faces); !i.finished(); ++i)
    {
    	assert(i.is_face());
      base_node un = meshRef.normal_of_face_of_convex(i.cv(), i.f());
      un /= gmm::vect_norm2(un);

			if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1]+1)<1.0e-7)
      {
      	meshRef.region(0).add(i.cv(),i.f());
      }
      if (gmm::abs(un[0]-1)<1.0e-7 && gmm::abs(un[1])<1.0e-7)
      {
      	meshRef.region(1).add(i.cv(),i.f());
      }
      if (gmm::abs(un[0])<1.0e-7 && gmm::abs(un[1]-1)<1.0e-7)
     	{
 	     meshRef.region(2).add(i.cv(),i.f());
      }
      if (gmm::abs(un[0]+1)<1.0e-7 && gmm::abs(un[1])<1.0e-7)
      {
		 	 meshRef.region(3).add(i.cv(),i.f());
      }
    }
}
