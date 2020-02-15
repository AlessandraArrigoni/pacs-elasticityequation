#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Core.h"

/*! @file LinearSystem.h
    @brief This is the class for the management of a linear system.  

   @details This class is endowed with different methods that let us copy, add or extract parts of the matrix, of the right hand side or of the solution vector. It also provides functions to get and set values, to save the matrix and algorithms that invert the matrix and solve the linear system.

*/


class LinearSystem
{
public:
	 LinearSystem ();

	 void addToMatrix(int ndof);

	 /*! method to copy the given matrix into the existing system starting from first_row, first_column */
	 void copySubMatrix(sparseMatrixPtr_Type source, int first_row, int first_column, scalar_type scale=1.0, bool transpose=false) ; 

	/*! method to add the given matrix to the existing system starting from first_row, first_column */
	 void addSubMatrix(sparseMatrixPtr_Type source, int first_row, int first_column, scalar_type scale=1.0, bool transpose=false) ; 

         /*! method that takes the submatrix of the existing system associated to the interval of indices specified by the ints and copies it into the destination */
	 void extractSubMatrix(sparseMatrixPtr_Type destination, int first_row, int number_rows, int first_column, int number_cols ) const ; 

	/*! method to copy the given vector into the existing system starting from first_row */
	 void copySubVector(scalarVectorPtr_Type source, int first_row, scalar_type scale=1.0) ; 

	/*! method to add the given vector to existing system starting from first_row */
	 void addSubVector(scalarVectorPtr_Type source, int first_row, scalar_type scale=1.0); 

	/*! method that takes the subvector of the existing system (either the solution or the rhs, according to "where") from first_row to the dimension of destination and copies it into the destination */
	 void extractSubVector(scalarVectorPtr_Type destination, int first_row, std::string where="sol") const ; 

	 void extractSubVector(scalarVector_Type & destination, int first_row, std::string where="sol") const ;

	 sparseMatrixPtr_Type inline getMatrix() const
	 {
	 	return M_Matrix;
	 }

	 scalarVectorPtr_Type inline getRHS() const
	 {
	 	return M_RHS;
	 }

   scalarVectorPtr_Type inline getSol() const
	 {
	 	return M_Sol;
	 }


	 void addSubSystem(LinearSystem* small, size_type shiftRows, size_type shiftColumns);

 	/*! resolution of the linear system */
	 void solve();

	 void computeInverse();

	 void saveMatrix(const char* nomefile="Matrix.mm") const;

	 void multAddToRHS(scalarVectorPtr_Type V, int first_row, int first_column, int nrows, int ncols);  //method that multiplies the matrix of the system for V and adds the result to RHS 
	 void multAddToRHS(sparseMatrixPtr_Type M, scalarVectorPtr_Type V, int first_rowVector, int first_rowRHS, scalar_type scale=1.0, bool transposed=false);  //method that multiplies the matrix M for V and adds the result to RHS 
	 void multAddToRHS(sparseMatrixPtr_Type M, scalarVector_Type& V, int first_rowVector, int first_rowRHS, scalar_type scale=1.0, bool transposed=false);

	/*! method that eliminates the rows and columns of the matrix of the system. It takes one argument as input.
	@param indexes indexes of rows and columns we want to discard from the original matrix of the system
        */
	 void eliminateRowsColumns(std::vector<size_type> indexes);

	 inline void cleanRHS()
	 {
	 	gmm::clear(*M_RHS);
	 }

	 inline void cleanMAT()
	 {
	 	gmm::clear(*M_Matrix);
	 }

	 inline void setNullRow(size_type which)
	 {
		for (size_type j=0; j<M_ndof;++j)  (*M_Matrix)(which,j)=0;
	}

	 inline void setNullColumn(size_type which)
	 {
		 for (size_type j = 0; j<M_ndof; ++j) (*M_Matrix)(j, which) = 0;
	 }

	 inline void setMatrixValue(size_type i, size_type j, scalar_type value)
	 {
		 (*M_Matrix)(i,j)=value;
	 }

	 inline void setRHSValue(size_type i, scalar_type value)
	 {
		 (*M_RHS)[i]=value;
	 }


private:

	sparseMatrixPtr_Type M_Matrix; // Pointer to a sparse matrix
	sparseMatrixPtr_Type M_InverseMatrix; // Pointer to a sparse matrix
	scalarVectorPtr_Type M_RHS; // Pointer to a vector of scalars
	scalarVectorPtr_Type M_Sol; // Pointer to a vector of scalars
	bool M_gotInverse;
	int M_ndof;

};



#endif
