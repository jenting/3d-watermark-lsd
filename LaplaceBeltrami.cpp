#include "LaplaceBeltrami.h"

/*
double getTheta(Point_3 p1, Point_3 p2, Point_3 p3)
{
	Vector_3 p21(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
	Vector_3 p31(p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z());

	int p21_length = sqrt(p21.squared_length());
	int p31_length = sqrt(p31.squared_length());

	double theta = acos(p21 * p31 / (p21_length* p31_length));

	return theta;
}

double calculateTriangleArea(Point_3 p1, Point_3 p2, Point_3 p3)
{
	Vector_3 p21(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());
	Vector_3 p31(p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z());

	int p21_length = sqrt(p21.squared_length());
	int p31_length = sqrt(p31.squared_length());

	double theta = acos(p21 * p31 / (p21_length* p31_length));
	double area = p21_length * p31_length * sin(theta) / 2;

	return area;
}

void LaplaceBeltrami::go(Mesh* mesh)
{
	int no_pts = mesh->size_of_vertices();

    // initialize
	vector<vector<double> > S(no_pts, vector<double>(no_pts, 0));
	vector<vector<double> > M(no_pts, vector<double>(no_pts, 0));
	for (int i = 0; i < no_pts; i++) {
		for (int j = 0; j < no_pts; j++) {
			S[i][j] = 0;
			M[i][j] = 0;
		}
	}

//	printf("initialize done!\n");

	// construct Laplace-Beltrami matrix
	Mesh::Vertex_iterator vit = mesh->vertices_begin();
	Mesh::Vertex_iterator vit_end = mesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		int v_index = vit->index();
		Point_3 p = vit->point();

//		printf("v_index = %d\n", v_index);

		Mesh::Vertex_handle pVertex = vit;

		Mesh::Halfedge_around_vertex_circulator pHalfEdge = pVertex->vertex_begin();
		Mesh::Halfedge_around_vertex_circulator end = pHalfEdge;

		int size = 0;
		CGAL_For_all(pHalfEdge, end)
		{
			size++;
		}

//		printf("size = %d\n", size);

		// get 1-ring neighbors
		vector<Point_3> v_n_points(size);
		vector<int> v_n_indexs(size);
		int ind = 0;
		CGAL_For_all(pHalfEdge, end)
		{
			int v_n_index = pHalfEdge->opposite()->vertex()->index();	
			Point_3 p_n_point = pHalfEdge->opposite()->vertex()->point();

			v_n_indexs[ind] = v_n_index;
			v_n_points[ind] = p_n_point;

			ind++;
		}
		
//		printf("get 1-ring neighbors done!\n");

		// calculate facet area
		// calculate cot(alpha) and cot(beta)
		double sum_area = 0;
		vector<double> cot_alpha(size);
		vector<double> cot_beta(size);
		ind = 0;
		CGAL_For_all(pHalfEdge, end)
		{
			sum_area += calculateTriangleArea(p, v_n_points[ind], v_n_points[(ind+1) % size]); // area

			double theta_alpha = getTheta(v_n_points[ind], v_n_points[(ind+1) % size], p);
			double theta_beta = getTheta(v_n_points[(ind+1) % size], v_n_points[ind], p);
			cot_alpha[ind] = 1 / tan(theta_alpha); // cot(alpha)
			cot_beta[ind] = 1 / tan(theta_beta); // cot(beta)

			ind++;
		}
		
		S[v_index][v_index] = sum_area; // matrix S

		double sum_val = 0;
		for (int i = 0; i < size; i++)
		{
			int v_n_index = v_n_indexs[i];
			//double weight = (1.5) * (cot_alpha[(i+size-1)%size] + cot_beta[i]); // matrix M, i and j adjacent
			double weight = (cot_alpha[(i+size-1)%size] + cot_beta[i]) / 2; // matrix M, i and j adjacent
			
			M[v_index][v_n_index] = (-1) * weight;
			sum_val += weight; // i == j
		}
		M[v_index][v_index] = sum_val;

	}

	printf("All Done!!!\n");
}
*/
// #include <math.h>
// #include "arlsmat.h"
// #include "arlgsym.h"

/*
  * Prints eigenvalues and eigenvectors of symmetric eigen-problems
  * on standard "cout" stream.
  */
/*
template<class FLOAT, class INT>
void Solution(INT nconv, INT n, INT nnz, FLOAT A[], INT irow[], INT pcol[],
              char uplo, FLOAT EigVal[], FLOAT* EigVec = 0)
{
  INT                  i;
  FLOAT*               Ax;
  FLOAT*               ResNorm;
  ARluSymMatrix<FLOAT> matrix(n, nnz, A, irow, pcol, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric eigenvalue problem: A*x - lambda*x \n \n";

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Finding the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrix.MultMv(&EigVec[i*n], Ax);
      axpy(n, -EigVal[i], &EigVec[i*n], 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution.
*/

/*
  * Prints eigenvalues and eigenvectors of symmetric generalized
  * eigen-problem on standard "cout" stream.
  */
/*
template<class FLOAT, class INT>
void Solution(INT nconv, INT n, INT nnzA, FLOAT A[], INT irowA[],
              INT pcolA[], INT nnzB, FLOAT B[], INT irowB[], INT pcolB[],
              char uplo, FLOAT EigVal[], FLOAT* EigVec = 0)
{

  INT                  i;
  FLOAT                *Ax, *Bx;
  FLOAT                *ResNorm;
  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  cout << endl << endl << "Testing ARPACK++ function AREig" << endl;
  cout << "Real symmetric generalized eigenvalue problem: A*x - lambda*B*x";
  cout << endl << endl;

  cout << "Dimension of the system            : " << n     << endl;
  cout << "Number of 'converged' eigenvalues  : " << nconv << endl << endl;

  // Printing eigenvalues.

  cout << "Eigenvalues:" << endl;

  for (i=0; i<nconv; i++) {
    cout << "  lambda[" << (i+1) << "]: " << EigVal[i] << endl;
  }
  cout << endl;

  // Printing eigenvectors.

  if (EigVec != 0) {

    // Printing the residual norm || A*x - lambda*B*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new FLOAT[n];
    Bx      = new FLOAT[n];
    ResNorm = new FLOAT[nconv+1];

    for (i=0; i<nconv; i++) {
      matrixA.MultMv(&EigVec[i*n], Ax);
      matrixB.MultMv(&EigVec[i*n], Bx);
      axpy(n, -EigVal[i], Bx, 1, Ax, 1);
      ResNorm[i] = nrm2(n, Ax, 1)/fabs(EigVal[i]);
    }

    for (i=0; i<nconv; i++) {
      cout << "||A*x(" << i << ") - lambda(" << i;
      cout << ")*B*x(" << i << ")||: " << ResNorm[i] << endl;
    }
    cout << endl;

    delete[] Ax;
    delete[] Bx;
    delete[] ResNorm;

  }

} // Solution.
*/

/**
  * @param matrix symmetric matrix
  * @param n dimension of symmetric matrix
  * @param nnz number of non-zero elements in matrix
  * @param val pointer to an array that sotres the non-zero elements fo matrix
  * @param irow pointer to an array that stores the row indices of the non-zero in matrix
  * @param pcol pointer to an array of pointers to the beginning of each column of matrix in val
  * @param uplo variable that indicates whether the upper (uplo='U') or the lower (uplo='L') part of val will be supplied to AREig
  */
/*
void symmMatrixToCSCFormat(vector<vector<double> > matrix, int n, int& nnz, double* &val, int* &irow, int* &pcol, char uplo = 'U')
{
	// calculate nnz
	nnz = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (matrix[i][j] != 0)
				nnz++;

	val = new double[nnz];
	irow = new int[nnz];
	pcol = new int[n+1];

	// Defining matrix val
	int ind_val = 0, ind_irow = 0, ind_pcol = 0;
	int sum_pcol = 0;
	bool start = false;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (matrix[j][i] !=0)
			{
				val[ind_val] = matrix[j][i];	irow[ind_irow] = j;	
				ind_val++;	ind_irow++;

				if (!start)
				{
					start = true;
					pcol[ind_pcol] = sum_pcol;
					ind_pcol++;
				}
				sum_pcol++;
			}
		}

		start = false;
	}
	pcol[ind_pcol] = sum_pcol; // the last

}

template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char uplo, int nev, char* which = "LM",
          int ncv = 0, FLOAT tol = 0.0, int maxit = 0,
          FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(nev, matrixA, matrixB, which,
                            ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors, regular mode.

template <class FLOAT>
int AREig(FLOAT EigVal[], int n, int nnzA, FLOAT A[], int irowA[],
          int pcolA[], int nnzB, FLOAT B[], int irowB[], int pcolB[],
          char uplo, char InvertMode, FLOAT sigma, int nev,
          char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                            which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues.

  return prob.Eigenvalues(EigVal);

} // real symmetric generalized problem, only eigenvalues,
  // shift-and-invert, buckling and Cayley modes.

template <class FLOAT>
int AREig(FLOAT EigVal[], FLOAT EigVec[], int n, int nnzA, FLOAT A[],
          int irowA[], int pcolA[], int nnzB, FLOAT B[], int irowB[],
          int pcolB[], char uplo, char InvertMode, FLOAT sigma,
          int nev, char* which = "LM", int ncv = 0, FLOAT tol = 0.0,
          int maxit = 0, FLOAT* resid = 0, bool AutoShift = true)
{

  // Creating two matrices in ARPACK++ format.

  ARluSymMatrix<FLOAT> matrixA(n, nnzA, A, irowA, pcolA, uplo);
  ARluSymMatrix<FLOAT> matrixB(n, nnzB, B, irowB, pcolB, uplo);

  // Defining the eigenvalue problem.

  ARluSymGenEig<FLOAT> prob(InvertMode, nev, matrixA, matrixB, sigma,
                            which, ncv, tol, maxit, resid, AutoShift);

  // Finding eigenvalues and eigenvectors.

  return prob.EigenValVectors(EigVec, EigVal);

} // real symmetric generalized problem, values and vectors,
  // shift-and-invert, buckling and Cayley modes.
*/

/*
void arpack()
{
	// solve the eigenvalues and eigenvectors of Laplace-Beltrami operator
	// ARluSymGenEig, Classes that require matrices, real symmetric, generalized
	// ARSymGenEig, Classes that require user-defined matrix-vector, real symmetric, generalized
	// ARluSymMatrix, SuperLU(CSC format), symmetric
	// ARumSymMatrix, UMFPACK(CSC format), symmetric
	// ARbdSymMatrix, LAPACK(band format), symmetric
	// ARdsSymMatrix, LAPACK(dense format), symmetric

	// declaring variables needed to store M_CSC, S_CSC in compressed sparse column(CSC) format

	// Defining variables;


	int n; // Dimension of the problem.
	int nnzA,   nnzB; // Number of nonzero elements in A and B.
	int *irowA, *irowB; // pointer to an array that stores the row indices of the nonzeros in A and B.
	int *pcolA, *pcolB; // pointer to an array of pointers to the beginning of each column of A (B) in valA (valB).
    double *valA,  *valB; // pointer to an array that stores the nonzero elements of A and B.
    double EigVal[101]; // Eigenvalues.
    double EigVec[1001]; // Eigenvectors stored sequentially.
    char uplo; // Variable that indicates whether the upper (uplo='U') ot the lower (uplo='L') part of A and B will be supplied to AREig.

    // Creating matrices A and B.
    n = 100;
    uplo = 'U';

    //SymmetricMatrixC(n, nnzA, valA, irowA, pcolA, uplo);
    //SymmetricMatrixD(n, nnzB, valB, irowB, pcolB, uplo);

    // Finding the four eigenvalues of A nearest to 1.0/150.0/0.0 and the related eigenvectors.

	int nconv; // Number of "converged" eigenvalues.
    //nconv = AREig(EigVal, EigVec, n, nnzA, valA, irowA, pcolA, nnzB, valB, irowB, pcolB, uplo, 'B', 1.0, 4); /// buckling mode
	
	nconv = AREig(EigVal, EigVec, n, nnzA, valA, irowA, pcolA, nnzB, valB, irowB, pcolB, uplo, 'C', 150.0, 4); /// Cayley mode

	//nconv = AREig(EigVal, EigVec, n, nnzA, valA, irowA, pcolA, nnzB, valB, irowB, pcolB, uplo, 'S', 0.0, 4); /// shift and invert mode

    // Printing solution.

    Solution(nconv, n, nnzA, valA, irowA, pcolA, nnzB,
             valB, irowB, pcolB, uplo, EigVal, EigVec);

}
*/

/*
void calculateFacetArea(Mesh* mesh)
{
	//area
	
	Mesh::Facet_const_iterator fit = mesh->facets_begin();
	Mesh::Facet_const_iterator fit_end = mesh->facets_end();
}


void getFeature(Mesh* mesh, double* featureArray)
{
	int ind = 0;
	Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
	Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		featureArray[ind] = vit->diff_length();
		ind++;
	}
}


void laplaceBeltramiValue(Mesh* mesh, int v_index, double* featureArray)
{

}
*/