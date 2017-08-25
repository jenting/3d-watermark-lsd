#ifndef __MY_MATRIX_H_
#define __MY_MATRIX_H_


#include "stdafx.h"


class MyMatrix {

public:
	MyMatrix();
	virtual ~MyMatrix();

	/** Calculate the cofactor of element (row,col) */
	int getMinor(double** src, double** dest, int row, int col, int order);

	/** Calculate the determinant recursively. */
	double calcDeterminant(double** mat, int order);

	/**Matrix inversion, the result is put in Y */
	void matrixInversion(double** A, int order, double** Y);

	/** Matrix multiplication */ 
	double** matrixMultiplication(double** a, double** b, int order);

	/** Matrix subtraction */
	double** matrixSubtraction(double** a, double** b, int order);

	/** Calculate Forbenius norm */
	double frobeniusNorm(double** m, int order);

	/** Calculate inverse matrix using Gaussian-Jordan elimination */
	bool getInverseMatrixGaussianJordanElimination(double** m, int order, double** mInv);
	bool getInverseMatrixGaussianJordanElimination(TNT::Array2D<double>& m, int order, TNT::Array2D<double>& mInv);

};

#endif