#ifndef __SORTING_H
#define __SORTING_H

#include "stdafx.h"

class QuickSort {

private:

protected:

	void swapVertexHandle(Mesh::Vertex_handle& v1, Mesh::Vertex_handle& v2);

	int pivotSquareErrorAscending(vector<Mesh::Vertex_handle>& a, int first, int last);

	int pivotVolumeSDFDescending(vector<Mesh::Vertex_handle>& a, int first, int last);
	int pivotVolumeSDFAscending(vector<Mesh::Vertex_handle>& a, int first, int last );

public:

	void printSquareError(const vector<Mesh::Vertex_handle>& a, int nElements);
	void printVolumeSDF(const vector<Mesh::Vertex_handle>& a, int nElements);
	
	void sortSquareErrorAscending(vector<Mesh::Vertex_handle>& a, int first, int last);

	void sortVolumeSDFDescending(vector<Mesh::Vertex_handle>& a, int first, int last);
	void sortVolumeSDFAscending(vector<Mesh::Vertex_handle>& a, int first, int last);
	
};

#endif