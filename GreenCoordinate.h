#ifndef __GREEN_COORDINATES_H
#define __GREEN_COORDINATES_H

#include "stdafx.h"

#include <math.h>

#include <q3textstream.h>
#include <qfile.h>

#include "AppObject.h"

class GreenCoordinate {

private:
	Point_3 Point_3_plus(Point_3 p1, Point_3 p2);

	Point_3 Point_3_minus(Point_3 p1, Point_3 p2);

	Point_3 Point_3_scalar_multi(float scalar, Point_3 p1);

	Point_3 Point_3_scalar_div(float scalar, Vector_3 p1);

	int sign_func(float value)
	{
		if(value >= 0)
			return 1;
		else
			return -1;
	}

protected:

	bool saveGC(
		const QString& vertexGCTxtFileName, 
		const QString& facetGCTxtFileName, 
		float* gc_coord, 
		int no_inter_pts, int no_cage_pts, int no_cage_trg);

	bool saveVertexGC(const QString& exportFileName, float* gc_v, int no_inter_pts, int no_cage_pts);

	bool saveFacetGC(const QString& exportFileName, float* gc_f, int no_inter_pts, int no_cage_trg);

	/**
	 * convert mesh vertices into 2d origins array
	 * @param mesh mesh
	 */
	vector<Point_3> preprocess_vertex(Mesh* mesh);

	/**
	 * convert mesh facets into 2d origins array
	 * @param mesh mesh
	 */
	vector<vector<int> > preprocess_facet(Mesh *mesh);

	/**
	 * convert the normal vector of all facets into 1d origins array
	 * @param mesh mesh
	 */
	vector<Vector_3> preprocess_facet_norm(Mesh* mesh);

	/**
	 * calculate the scaling factor = ||n`|| / ||n||
	 */
	vector<float> get_scaling_factor(vector<Vector_3> cage_trg_norm, vector<Vector_3> deform_cage_trg_norm);

	/**
	 * calculate interior points coordinate according to local_coordinate and deform_cage_point
	 * v` = summation(wi * vi) / summation(wi)
	 * @param no_inter_pts no of interior points
	 * @param no_cage_pts  no of cage points
	 * @param local_coordinate the 1d array that save the coordinate of each vertices, [i, i+1, i+2] means (x,y,z) coordinate
	 * @param deform_cage_point the 2d array that save the coordinate of cage points
	 */
	vector<Point_3> postprocess(
		float* gc_coord, int no_inter_pts, 
		vector<Point_3> deform_cage_vertices, 
		vector<Vector_3> deform_cage_trg_norm, 
		vector<float> deform_scaling_factor);


	float TriInt(Point_3 p, Point_3 v1, Point_3 v2, Point_3 eta);

	float TriIntTheta(float lambda, float c, float theta);

	float* gc_all(
		vector<Point_3> interior_points, 
		vector<Point_3> mesh_coord, 
		vector<vector<int> > mesh_trg,
		vector<Vector_3> mesh_trg_norm);

	float* gc(
		Point_3 x, int point_index, 
		vector<Point_3> mesh_coord, int no_cage_pts, 
		vector<vector<int> >mesh_trg, 
		vector<Vector_3> mesh_trg_norm, int no_cage_trgs);


public:
	GreenCoordinate();
    virtual ~GreenCoordinate();

	bool go(AppObject* appObject, Mesh* origin_mesh, Mesh* cage_mesh, Mesh* deform_cage_mesh);

};

#endif