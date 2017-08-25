#include "MyMatrix.h"


MyMatrix::MyMatrix()
{
	//ctor
}

MyMatrix::~MyMatrix()
{
	//dtor
}


/**
  *calculate the cofactor of element (row,col)
  */
int MyMatrix::getMinor(double **src, double **dest, int row, int col, int order)
{
	// indicate which col and row is being copied to dest
	int colCount=0, rowCount=0;

	for(int i = 0; i < order; i++ )
	{
		if( i != row )
		{
			colCount = 0;
			for(int j = 0; j < order; j++ )
			{
				// when j is not the element
				if( j != col )
				{
					dest[rowCount][colCount] = src[i][j];
					colCount++;
				}
			}
			rowCount++;
		}
	}

	return 1;
}


/** 
  * Calculate the determinant recursively.
  */
double MyMatrix::calcDeterminant(double **mat, int order)
{
	//cout << "order= " << order << endl;

	// order must be >= 0
	/// stop the recursion when Matr is a single element ///
	if( order == 1 )
		return mat[0][0];

	/// the determinant value ///
	double det = 0;

	/// allocate the cofactor Matr ///
	double** minor = new double* [order-1];
	for (int i = 0; i < (order-1); i++)
		minor[i] = new double[order-1];

	for(int i = 0; i < order; i++ )
	{
		// get minor of element (0,i)
		getMinor(mat, minor, 0, i, order);

		// the recusion is here!
		det += pow( -1.0, i ) * mat[0][i] * calcDeterminant( minor,order-1 );
	}

	/// release memory ///
	for(int i=0;i<order-1;i++)
		delete [] minor[i];
	delete [] minor;

	return det;
}

/*
double MyMatrix::calcDeterminant(double **mat, int order)
{
	double det = 0.0;
	for (int i = 0; i < order; i++)
	{
		double a = 1.0, b = 1.0;
		for (int row = 0; row < order; row++)
		{
			a *= mat[row][(i+row)%order];
			b *= mat[row][(order-1) - (i+row)%order];
		}

		det += a - b;
	}

	return det;
}
*/


/**
  * Matr inversioon, the result is put in Y
  */
void MyMatrix::matrixInversion(double** A, int order, double** Y)
{
	/// get the determinant of a ///
	double det = 1.0 / calcDeterminant(A,order);

	//cout << "calculate det done\n";

	/// memory allocation ///
	//cout << "start to allocate memory...\n";

	double* temp = new double[(order-1) * (order-1)];
	
	//cout << "temp memory allocation done\n";

	double** minor = new double*[order-1];

	//cout << "temp minor allocation done\n";

	for(int i = 0; i < (order - 1); i++)
		minor[i] = temp+(i * (order-1));

	//cout << "memory allocation done\n";

	for(int j = 0; j < order; j++)
	{
		for(int i = 0; i < order; i++)
		{
			// get the co-factor (Matr) of A(j,i)
			getMinor(A, minor, j, i, order);
			Y[i][j] = det * calcDeterminant(minor, order-1);

			if( (i+j) % 2 == 1)
				Y[i][j] = -Y[i][j];
		}
	}

	/// release memory ///
	delete [] minor[0];
	delete [] minor;
}


/**
  * Matrix multiplication
  */ 
double** MyMatrix::matrixMultiplication(double** a, double** b, int order)
{
	double** result = new double*[order];
	for(int i = 0; i < order; i++)
		result[i] = new double[order];

	for(int i = 0; i < order; i++)
		for(int j = 0; j < order; j++)
		{
			double sum = 0;
			for(int k = 0; k < order; k++)
				sum += a[i][k] * b[k][j];
			result[i][j] = sum;
		}
				
	return result;
}


/**
  * Matr subtraction
  */
double** MyMatrix::matrixSubtraction(double** a, double** b, int order)
{
	double** result = new double*[order];
	for(int i = 0; i < order; i++)
		result[i] = new double[order];

	for(int i = 0; i < order; i++)
		for(int j = 0; j < order; j++)
			result[i][j] = a[i][j] - b[i][j];

	return result;
}


/**
  * calculate Forbenius norm
  */
double MyMatrix::frobeniusNorm(double** m, int order)
{
	double sum = 0;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			sum += pow(m[i][j], 2);

	return sqrt(sum);
}


/**
  * calculate inverse matrix using Gaussian-Jordan elimination
  */
bool MyMatrix::getInverseMatrixGaussianJordanElimination(double** m, int order, double** mInv)
{
	int j, i, k, L;
	double s, t;

	for (j = 0; j < order; j++)
	{

		for (i = j; i < order; i++)
		{

			if ( m[i][j] != 0 )
			{

				for (k = 0; k < order; k++)
				{
					s = m[j][k]; m[j][k] = m[i][k]; m[i][k] = s;
					s = mInv[j][k]; mInv[j][k] = mInv[i][k]; mInv[i][k] = s;
				}

				t = 1 / m[j][j];

				for (k = 0; k < order; k++)
				{
					m[j][k] = t * m[j][k];
					mInv[j][k] = t * mInv[j][k];
				}

				for (L = 0; L < order; L++)
				{
	                if ( L != j )
					{
						t = -m[L][j];

						for (k = 0; k < order; k++)
						{
	                        m[L][k] = m[L][k] + t * m[j][k];
							mInv[L][k] = mInv[L][k] + t * mInv[j][k];
						}

					}
				}

			}
			
			break;
		}

		/// Display warning if a row full of zeros is found ///
		if ( m[i][j] == 0 )
		{
			cout << "Warning: Singular Matrix\n";
			return true;
		}

	}

	return false;
}


/**
  * calculate inverse matrix using Gaussian-Jordan elimination
  */
bool MyMatrix::getInverseMatrixGaussianJordanElimination(TNT::Array2D<double>& m, int order, TNT::Array2D<double>& mInv)
{
	int j, i, k, L;
	double s, t;

	for (j = 0; j < order; j++)
	{

		for (i = j; i < order; i++)
		{

			if ( m[i][j] != 0 )
			{

				for (k = 0; k < order; k++)
				{
					s = m[j][k]; m[j][k] = m[i][k]; m[i][k] = s;
					s = mInv[j][k]; mInv[j][k] = mInv[i][k]; mInv[i][k] = s;
				}

				t = 1.0 / m[j][j];

				for (k = 0; k < order; k++)
				{
					m[j][k] = m[j][k] * t;
					mInv[j][k] = mInv[j][k] * t;
				}

				for (L = 0; L < order; L++)
				{
	                if ( L != j )
					{
						t = -m[L][j];

						for (k = 0; k < order; k++)
						{
	                        m[L][k] = m[L][k] + t * m[j][k];
							mInv[L][k] = mInv[L][k] + t * mInv[j][k];
						}

					}
				}

			}
			
			break;
		}

		/// Display warning if a row full of zeros is found ///
		if ( m[i][j] == 0 )
		{
			cout << "Warning: Singular Matrix\n";
			return true;
		}

	}

	return false;
}