#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Core.h"


// classe per i sistemi lineari collegati ai vari problemi
class LinearSystem
{
public:
	 LinearSystem ();
	 void addToMatrix(int ndof);
	 void copySubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column, scalar_type scale=1.0, bool transpose=false);
	 void addSubMatrix(sparseMatrixPtr_Type M, int first_row, int first_column, scalar_type scale=1.0, bool transpose=false);

	 void copySubVector(scalarVectorPtr_Type M, int first_row,scalar_type scale=1.0);
	 void addSubVector(scalarVectorPtr_Type M, int first_row,scalar_type scale=1.0);
	 void extractSubVector(scalarVectorPtr_Type M, int first_row, std::string where="sol");

	 sparseMatrixPtr_Type inline getMatrix()
	 {
	 	return M_Matrix;
	 }

	 scalarVectorPtr_Type inline getRHS()
	 {
	 	return M_RHS;
	 }

  	 scalarVectorPtr_Type inline getSol()
	 {
	 	return M_Sol;
	 }

         //questo per problemi a blocchi
	 void addSubSystem(LinearSystem* small, size_type shiftRows, size_type shiftColumns);

	 void solve();

	 void computeInverse();

	 void saveMatrix(const char* nomefile="Matrix.mm");

	 void multAddToRHS(scalarVectorPtr_Type V, int first_row, int first_column, int nrows, int ncols);  //moltiplica la matrice del sistema per V e aggiunge il risultato a RHSe
	 void multAddToRHS(sparseMatrixPtr_Type M, scalarVectorPtr_Type V, int first_rowVector, int first_rowRHS, scalar_type scale=1.0, bool transposed=false);  //moltiplica M per V e aggiunge il risultato a RHS
	 void multAddToRHS(sparseMatrixPtr_Type M, scalarVector_Type& V, int first_rowVector, int first_rowRHS, scalar_type scale=1.0, bool transposed=false);

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
	 inline void setMatrixValue(size_type i, size_type j, scalar_type value)
	 {
		 (*M_Matrix)(i,j)=value;
	 }

	 inline void setRHSValue(size_type i,scalar_type value)
	 {
		 (*M_RHS)[i]=value;
	 }
private:

	sparseMatrixPtr_Type M_Matrix;
	sparseMatrixPtr_Type M_InverseMatrix;
	scalarVectorPtr_Type M_RHS; // Puntatore a un vettore di scalari
	scalarVectorPtr_Type M_Sol; // Puntatore a un vettore di scalari
	bool M_gotInverse;
	int M_ndof;

};



#endif
