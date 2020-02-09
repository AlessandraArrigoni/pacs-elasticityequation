
#include "../include/LinearSystem.h"

LinearSystem::LinearSystem( ): M_Matrix(), M_RHS(), M_Sol() // Con il constructor ho creato dei puntatori vuoti, poi userò initialize e addToMatrix per riempirli
{
   M_gotInverse=false;
}

void LinearSystem::addToMatrix(int ndof)
{
	M_ndof = ndof;
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

void LinearSystem::copySubMatrix(sparseMatrixPtr_Type source, int first_row, int first_column, scalar_type scale, bool transpose)
{
	if (transpose)
	{
	gmm::copy ( gmm::transposed(gmm::scaled(*source, scale)),
							gmm::sub_matrix (*M_Matrix, gmm::sub_interval (first_row,gmm::mat_ncols(*source)), gmm::sub_interval (first_column, gmm::mat_nrows(*source)) ) );
	}
	else
	{
	gmm::copy ( gmm::scaled(*source, scale),
							gmm::sub_matrix (*M_Matrix, gmm::sub_interval (first_row,gmm::mat_nrows(*source)), gmm::sub_interval (first_column, gmm::mat_ncols(*source)) ) );
	}
}


void LinearSystem::addSubMatrix(sparseMatrixPtr_Type source, int first_row, int first_column, scalar_type scale, bool transpose)
{
	if (transpose)
	{
	gmm::add ( gmm::transposed(gmm::scaled(*source, scale)),
						 gmm::sub_matrix (*M_Matrix, gmm::sub_interval (first_row,gmm::mat_ncols(*source)), gmm::sub_interval (first_column, gmm::mat_nrows(*source)) ) );
	}
	else
	{
	gmm::add ( gmm::scaled(*source, scale),
	 					 gmm::sub_matrix (*M_Matrix, gmm::sub_interval (first_row,gmm::mat_nrows(*source)), gmm::sub_interval (first_column, gmm::mat_ncols(*source)) ) );
	}
}


void LinearSystem::extractSubMatrix(sparseMatrixPtr_Type destination, int first_row, int number_rows, int first_column, int number_cols ) const
{

	gmm::copy(gmm::sub_matrix(*M_Matrix, gmm::sub_interval(first_row, number_rows), gmm::sub_interval(first_column, number_cols)), *destination);

}


void LinearSystem::copySubVector(scalarVectorPtr_Type source, int first_row, scalar_type scale)
{
	gmm::copy ( gmm::scaled(*source,scale),
							gmm::sub_vector (*M_RHS, gmm::sub_interval (first_row, source->size()) ));
}


void LinearSystem::addSubVector(scalarVectorPtr_Type source, int first_row, scalar_type scale)
{
	gmm::add ( gmm::scaled(*source, scale),
						 gmm::sub_vector (*M_RHS, gmm::sub_interval (first_row, source->size()) ));
}

// Con il PUNTATORE
void LinearSystem::extractSubVector(scalarVectorPtr_Type destination, int first_row, std::string where) const
{
	if (where=="sol")
	{
		gmm::copy ( gmm::sub_vector (*M_Sol, gmm::sub_interval (first_row, destination->size()) ), *destination);
	}
	else
	{
		gmm::copy ( gmm::sub_vector (*M_RHS, gmm::sub_interval (first_row, destination->size()) ), *destination);
	}
}

// Con la REFERENCE
void LinearSystem::extractSubVector(scalarVector_Type & destination, int first_row, std::string where) const
{
  if (where=="sol")
	{
		gmm::copy ( gmm::sub_vector (*M_Sol, gmm::sub_interval (first_row, destination.size()) ), destination);
	}
	else
	{
		gmm::copy ( gmm::sub_vector (*M_RHS, gmm::sub_interval (first_row, destination.size()) ), destination);
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

   /*
   std::cout << "\nSoluzione nel solve del sistema lineare "<<std::endl;
   for (size_type k=0; k<M_Sol->size(); k++)
   {
     std::cout << M_Sol->at(k)<< "\t"<<std::endl;
   }
   std::cout << "\nFINE Soluzione nel solve del sistema lineare "<<std::endl;
   */
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

void LinearSystem::saveMatrix(const char* nomefile) const
{
	gmm::MatrixMarket_IO::write(nomefile , *M_Matrix);
}



void LinearSystem::eliminateRowsColumns(std::vector<size_type> indexes)
{
// ordina vettore in oridine crescente
sort(indexes.begin(),indexes.end());


// eliminare le righe e le colonne una alla volta, a partire da quella "più bassa" e quella "più a destra" nella matrice 
// contemporaneamente elimina le stesse righe dal rhs 

for (size_type k=0; k<indexes.size(); k++)
{
  size_type idx=indexes[indexes.size()-1-k]; // indice riga e colonna che sono da eliminare in ordine decrescente
 // std::cout<< "idx è:"<< idx<<std::endl;

sparseMatrixPtr_Type A1 = std::make_shared<sparseMatrix_Type>(idx,idx);
sparseMatrixPtr_Type A2 = std::make_shared<sparseMatrix_Type>(idx,M_ndof-idx-1);
sparseMatrixPtr_Type A3 = std::make_shared<sparseMatrix_Type> (M_ndof-idx-1,idx);
sparseMatrixPtr_Type A4 = std::make_shared<sparseMatrix_Type>(M_ndof-idx-1,M_ndof-idx-1);

extractSubMatrix(A1,0,idx,0,idx);
extractSubMatrix(A2,0,idx,idx+1,M_ndof-idx-1);
extractSubMatrix(A3,idx+1,M_ndof-idx-1,0,idx);
extractSubMatrix(A4,idx+1,M_ndof-idx-1,idx+1,M_ndof-idx-1);

cleanMAT();
M_Matrix->resize(M_ndof-1,M_ndof-1);
addSubMatrix(A1, 0, 0);
addSubMatrix(A2, 0, idx);
addSubMatrix(A3, idx, 0);
addSubMatrix(A4, idx, idx);

scalarVectorPtr_Type RHS1= std::make_shared<scalarVector_Type> (idx);
scalarVectorPtr_Type RHS2 = std::make_shared<scalarVector_Type> (M_ndof-1-idx);


extractSubVector(RHS1, 0, "M_RHS");
extractSubVector(RHS2, idx+1, "M_RHS");

cleanRHS();
M_RHS->resize(M_ndof-1);
addSubVector(RHS1, 0);
addSubVector(RHS2, idx);

M_Sol->resize(M_ndof-1);

M_ndof=M_ndof-1;

}

}
