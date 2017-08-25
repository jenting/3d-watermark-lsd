#ifndef __MEAN_VALUE_COORDINATES_H
#define __MEAN_VALUE_COORDINATES_H

#include "stdafx.h"

#include <math.h> 

#include <q3textstream.h>
#include <qfile.h>

#include "AppObject.h"

//#include "Mesh/Polyhedron_Copy.h"

//template <class Polyhedron, class kernel>
class MeanValueCoordinate/* : public QObject*/
{
	//Q_OBJECT
private:
	//Mesh* m_MeshOriginal;
	//Mesh* m_MeshDeformed;

	//CCopyPoly<Polyhedron, kernel> copier;

private:

	/** 
	 * calculate the distance between two points
	 * 
	 * @param pt1 point1
	 * @param pt2 point2
	 */
	float distance(float pt1[3], float pt2[3]);

	/** 
	 * calculate the length of a vector
	 *
	 * @param pt vector
	 */
	float length(float pt[3]);

	/** 
	 * return sign(det[v1 v2 v3])
	 *
	 * @param v0 1th column vector
	 * @param v1 2th column vector
	 * @param v2 3th column vector
	 */
	float sign_det(float v0[3], float v1[3], float v2[3]);

	Point_3 Point_3_plus(Point_3 p1, Point_3 p2);

	Point_3 Point_3_scalar_multi(float scalar, Point_3 p1);

	Point_3 Point_3_scalar_div(float scalar, Point_3 p1);

protected:
	/**
	 *
	 */
	bool saveMVC(const QString& exportFileName, float* mvc, int no_inter_pts, int no_cage_pts);

	/**
	 * convert mesh vertices into 2d origins array
	 * @param mesh mesh
	 * @param origins the 2d array that save all the mesh vertices coordinate
	 */
	void preprocess_vertex(Mesh* mesh, float (*origins)[3]);

	vector<Point_3> preprocess_vertex(Mesh* mesh);

	/**
	 * convert mesh facets into 3d origins array
	 * @param mesh mesh
	 * @param origins the 2d array that save the vertex index of all facets
	 */
	void preprocess_facet(Mesh *mesh, int (*origins)[3]);

	vector<vector<int> > preprocess_facet(Mesh *mesh);

	/**
	 * calculate interior points coordinate according to local_coordinate and deform_cage_point
	 * v` = summation(wi * vi) / summation(wi)
	 * @param no_inter_pts no of interior points
	 * @param no_cage_pts  no of cage points
	 * @param local_coordinate the 1d array that save the coordinate of each vertices, [i, i+1, i+2] means (x,y,z) coordinate
	 * @param deform_cage_point the 2d array that save the coordinate of cage points
	 */
	//float* postprocess(int no_inter_pts, int no_cage_pts, float *local_coordinate, float (*deform_cage_point)[3]);
	
	vector<Point_3> postprocess(int no_inter_pts, int no_cage_pts, float *local_coordinate, vector<Point_3>deform_cage_point);
	
	/**
	 * this does the job for one point p, returns the value of the interpolant at this point
     * the output is an array of size 1 x nofuncs
     * mesh_funcs[i*nopts + j] contains the value of i-th function at vertex j
	 *
	 * @param x interior point
	 * @param mesh_func function defined on vertices
	 * @param nofunc no of mesh_func
	 * @param mesh_coord
	 * @param nopts no of mesh_coord
	 * @param mesh_triang
	 * @param notrg no of mesh_triang
	 */
	float* mvc(float x[3], int point_index, /*float mesh_funcs[],  int nofuncs,*/ float mesh_coord[][3], int nopts, int mesh_triang[][3], int notrg);


	/**
	 * this does the job for all interior points, returns an array contains the values of the interpolant of all the points
	 * @param interior_point interior points(original mesh vertices)
	 * @param no_inter_pts no of interior points(original mesh vertices)
	 * @param mesh_funcs functions defined on all the original mesh vertices
	 * @param nofuncs no of mesh_funcs
	 * @param mesh_coord cage mesh vertices coordinates
	 * @param nopts no of mesh_coord
	 * @param mesh_triang cage mesh facets
	 * @param notrg no of mesh_triang
	 */
	float* mvc_all(float (*interior_points)[3], int no_inter_pts, /*float *mesh_funcs, int nofuncs,*/ float (*mesh_coord)[3], int nopts, int (*mesh_triang)[3], int notrg);

	//float* mvc_all(vector<Point_3> interior_points, int no_inter_pts, vector<float> mesh_funcs, int nofuncs, vector<Point_3> mesh_coord, int nopts, vector<vector<int> > mesh_triang, int notrg);

	/*void PreProcessingDeformMesh()
	{
		m_MeshDeformed->clear();
		//copier.copy(m_MeshOriginal, m_MeshDeformed);
		m_MeshDeformed->computeNormals();
		m_MeshDeformed->compute_type();
	}*/

public:

	MeanValueCoordinate();

    virtual ~MeanValueCoordinate();

	bool go(AppObject* appObject, Mesh* origin_mesh, Mesh* cage_mesh, Mesh* deform_cage_mesh);

	/*MeanValueCoordinate(Mesh* _mesh)
	{
		MeshOriginal(_mesh);

		m_MeshDeformed = new Mesh;

		m_MeshOriginal->computeNormals();
		m_MeshOriginal->compute_type();
	}

	~MeanValueCoordinate()
	{
		delete m_MeshDeformed;
		m_MeshDeformed = 0;
	}*/

	/*const Mesh* MeshOriginal() const { return m_MeshOriginal; }
	Mesh* MeshOriginal() { return m_MeshOriginal; }
	void MeshOriginal(Mesh* m) { m_MeshOriginal = m; }

	const Mesh* MeshDeformed() const { return m_MeshDeformed; }
	Mesh* MeshDeformed() { return m_MeshDeformed; }
	void MeshDeformed(Mesh* m) { m_MeshDeformed = m; }*/

};

#endif