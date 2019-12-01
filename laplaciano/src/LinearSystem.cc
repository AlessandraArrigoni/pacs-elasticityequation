
#include "../include/LinearSystem.h"

LinearSystem::LinearSystem( ):
		M_Matrix(),
		M_RHS()
{
   M_gotInverse=false;
}

void LinearSystem::addToMatrix(int ndof)
{
	M_ndof=ndof;
	if (M_RHS!=NULL)
	{
		int provv_size(M_RHS->size());
		M_RHS->resize(provv_size+ndof);
		M_Sol->resize(provv_size+ndof);

		M_Matrix->resize(provv_size+ndof,provv_size+ndof);
	}
	else
	{
		M_RHS.reset(new scalarVector_Type (ndof));
		M_Sol.reset(new scalarVector_Type (ndof));

		M_Matrix.reset(new sparseMatrix_Type (ndof,ndof));
	}
}

void LinearSystem::copySubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column, scalar_type scale, bool transpose)
{
	if (transpose)
	{
	gmm::copy ( gmm::transposed(gmm::scaled(*M,scale)), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_ncols(*M)),
                        gmm::sub_interval (first_column, gmm::mat_nrows(*M)) ) );
	}
	else
	{
	gmm::copy ( gmm::scaled(*M,scale), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_nrows(*M)),
                        gmm::sub_interval (first_column, gmm::mat_ncols(*M)) ) );
	}
}

void LinearSystem::addSubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column, scalar_type scale, bool transpose)
{
	if (transpose)
	{
	gmm::add ( gmm::transposed(gmm::scaled(*M,scale)), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_ncols(*M)),
                        gmm::sub_interval (first_column, gmm::mat_nrows(*M)) ) );
	}
	else
	{
	gmm::add ( gmm::scaled(*M,scale), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (first_row,gmm::mat_nrows(*M)),
                        gmm::sub_interval (first_column, gmm::mat_ncols(*M)) ) );
	}
}

void LinearSystem::copySubVector(scalarVectorPtr_Type M, int first_row, scalar_type scale)
{

	gmm::copy ( gmm::scaled(*M,scale), gmm::sub_vector (*M_RHS,
                        gmm::sub_interval (first_row,M->size()) ));
}

void LinearSystem::addSubVector(scalarVectorPtr_Type M, int first_row, scalar_type scale)
{
	gmm::add ( gmm::scaled(*M,scale), gmm::sub_vector (*M_RHS,
                        gmm::sub_interval (first_row, M->size()) ));
}

void LinearSystem::extractSubVector(scalarVectorPtr_Type M, int first_row, std::string where)
{
	if (where=="sol")
	{
		gmm::copy ( gmm::sub_vector (*M_Sol, gmm::sub_interval (first_row, M->size()) ), *M);
	}
	else
	{
		gmm::copy ( gmm::sub_vector (*M_RHS, gmm::sub_interval (first_row,M->size()) ), *M);
	}
}

void LinearSystem::addSubSystem(LinearSystem* small, size_type shiftRows, size_type shiftColumns)
{
	gmm::add ( *(small->getMatrix()), gmm::sub_matrix (*M_Matrix,
                        gmm::sub_interval (shiftRows,gmm::mat_nrows(*(small->getMatrix()))),
                        gmm::sub_interval (shiftColumns, gmm::mat_ncols(*(small->getMatrix()))) ) );

        gmm::add ( *(small->getRHS()), gmm::sub_vector (*M_RHS,
                        gmm::sub_interval (shiftRows,gmm::mat_nrows(*(small->getMatrix()))) ) );

}


void LinearSystem::multAddToRHS(scalarVectorPtr_Type V, int first_row, int first_column, int nrows, int ncols)
{
	int length((*V).size());
	if (ncols==(*V).size()){
	gmm::mult_add(gmm::sub_matrix(*M_Matrix, gmm::sub_interval (first_row,nrows),gmm::sub_interval (first_column,length)),*V, gmm::sub_vector(*M_RHS,gmm::sub_interval(first_row,nrows)));
	}
	else
	{
		std::cout << "dimension mismatch"<<std::endl;
	}
}

void LinearSystem::multAddToRHS(sparseMatrixPtr_Type M, scalarVectorPtr_Type V,  int first_rowVector, int first_rowRHS, scalar_type scale,bool transposed)
{
	if (transposed)
	{
        gmm::mult_add(gmm::scaled(gmm::transposed(*M),scale), gmm::sub_vector(*V,gmm::sub_interval(first_rowVector,gmm::mat_nrows(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval		(first_rowRHS,gmm::mat_ncols(*M))));
	}
	else
	{
	gmm::mult_add(gmm::scaled(*M,scale), gmm::sub_vector(*V,gmm::sub_interval(first_rowVector,gmm::mat_ncols(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval(first_rowRHS,gmm::mat_nrows(*M))));
	}
}

void LinearSystem::multAddToRHS(sparseMatrixPtr_Type M, scalarVector_Type& V,  int first_rowVector, int first_rowRHS, scalar_type scale,bool transposed)
{
	if (transposed)
	{
	gmm::mult_add(gmm::scaled(gmm::transposed(*M),scale),gmm::sub_vector(V,gmm::sub_interval(first_rowVector,gmm::mat_nrows(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval(first_rowRHS,gmm::mat_ncols(*M))));
	}
	else
	{
	gmm::mult_add(gmm::scaled(*M,scale),gmm::sub_vector(V,gmm::sub_interval(first_rowVector,gmm::mat_ncols(*M))), gmm::sub_vector(*M_RHS,gmm::sub_interval(first_rowRHS,gmm::mat_nrows(*M))));
	}
}

void LinearSystem::solve()
{
   if (M_gotInverse)
   {
        gmm::clear(*M_Sol);
	gmm::mult(*M_InverseMatrix, *M_RHS, *M_Sol);
   }
   else
   {
   	scalar_type rcond;
  	gmm::clear(*M_Sol);
	SuperLU_solve(*M_Matrix, *M_Sol, *M_RHS, rcond);
   }
}

void LinearSystem::computeInverse()  //lento
{
    std::vector<scalar_type> RHS(M_ndof,0.0);
    std::vector<scalar_type> provv(M_ndof,0.0);
    M_InverseMatrix.reset(new sparseMatrix_Type (M_ndof,M_ndof));

    for (int i=0; i<M_ndof;++i)
    {
	gmm::clear(RHS);
	gmm::clear(provv);

	RHS[i]=1;
        scalar_type rcond;
        SuperLU_solve(*M_Matrix, provv, RHS, rcond);
	gmm::copy(gmm::col_vector(provv), gmm::sub_matrix(*M_InverseMatrix, gmm::sub_interval(0, M_ndof),  gmm::sub_interval(0,1)));

    }

    M_gotInverse=true;

}

void LinearSystem::saveMatrix(const char* nomefile)
{
	gmm::MatrixMarket_IO::write(nomefile , *M_Matrix);
}
