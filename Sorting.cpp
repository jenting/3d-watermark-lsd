#include "Sorting.h"


/// Swap:  Swap two item(call-by-reference). ///
void QuickSort::swapVertexHandle(Mesh::Vertex_handle& v1, Mesh::Vertex_handle& v2)
{
	Mesh::Vertex_handle tmpVal = v1;
    v1 = v2;
    v2 = tmpVal;
}


///  Quicksort:  Sort an array a, using the quicksort algorithm. ///
void QuickSort::sortVolumeSDFAscending(vector<Mesh::Vertex_handle>& a, int first, int last )
{
    int p;

    if( first < last )
	{
        p = pivotVolumeSDFAscending( a, first, last );
        sortVolumeSDFAscending( a, first, p - 1 );
        sortVolumeSDFAscending( a, p + 1, last );
    }
}


/// Pivot:  Find and return the index of pivot element. ///
int QuickSort::pivotVolumeSDFAscending(vector<Mesh::Vertex_handle>& a, int first, int last ) 
{
    int  p = first;
	double pivot = a[first]->volumeSDF();

    for( int i = first+1 ; i <= last ; i++ )
	{
        if( a[i]->volumeSDF() <= pivot )
		{
            p++;
			swapVertexHandle( a[i], a[p] );
        }
    }

    swapVertexHandle( a[p], a[first] );

    return p;
}


///  Quicksort:  Sort an array a, using the quicksort algorithm. ///
void QuickSort::sortVolumeSDFDescending(vector<Mesh::Vertex_handle>& a, int first, int last )
{
    int p;

    if( first < last )
	{
        p = pivotVolumeSDFDescending( a, first, last );
        sortVolumeSDFDescending( a, first, p - 1 );
        sortVolumeSDFDescending( a, p + 1, last );
    }
}


/// Pivot:  Find and return the index of pivot element. ///
int QuickSort::pivotVolumeSDFDescending(vector<Mesh::Vertex_handle>& a, int first, int last ) 
{
    int  p = first;
	double pivot = a[first]->volumeSDF();

    for( int i = first+1 ; i <= last ; i++ )
	{
        if( a[i]->volumeSDF() >= pivot )
		{
            p++;
			swapVertexHandle( a[i], a[p] );
        }
    }

    swapVertexHandle( a[p], a[first] );

    return p;
}


/// PrintArray:  Print contents of an array. ///
void  QuickSort::printVolumeSDF(const vector<Mesh::Vertex_handle>& A, int nElements)
{
    cout << "[ ";

    for( int i = 0 ; i < nElements ; i++ )
    {
        cout << A[i]->volumeSDF();
        if( i < nElements-1 )
           cout << ", ";
    }

    cout << " ] " << endl;
}


///////////////////////////////////////////////////////////////////


///  Quicksort:  Sort an array a, using the quicksort algorithm. ///
void QuickSort::sortSquareErrorAscending(vector<Mesh::Vertex_handle>& a, int first, int last )
{
    int p;

    if( first < last )
	{
        p = pivotSquareErrorAscending( a, first, last );
        sortSquareErrorAscending( a, first, p - 1 );
        sortSquareErrorAscending( a, p + 1, last );
    }
}


/// Pivot:  Find and return the index of pivot element. ///
int QuickSort::pivotSquareErrorAscending(vector<Mesh::Vertex_handle>& a, int first, int last ) 
{
    int  p = first;
	double pivot = a[first]->squareErrors();

    for( int i = first+1 ; i <= last ; i++ )
	{
        if( a[i]->squareErrors() <= pivot )
		{
            p++;
			swapVertexHandle( a[i], a[p] );
        }
    }

    swapVertexHandle( a[p], a[first] );

    return p;
}



/// PrintArray:  Print contents of an array. ///
void  QuickSort::printSquareError(const vector<Mesh::Vertex_handle>& A, int nElements)
{
    cout << "[ ";

    for( int i = 0 ; i < nElements ; i++ )
    {
        cout << A[i]->squareErrors();
        if( i < nElements-1 )
           cout << ", ";
    }

    cout << " ] " << endl;
}
