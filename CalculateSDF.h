#ifndef __CALCULATE_SDF_H_
#define __CALCULATE_SDF_H_

#include "stdafx.h"

#include <qobject.h>
#include <qthread.h>

#include <QApplication.h>

class QProgressDialog;

enum NormalizeTypes {None, MinMax, Log};


class CalculateSDF : public QObject 
{
	Q_OBJECT

private:
	RayIntersect m_rayIntersect;

public:

	double logit(const double& value, const int& range);

	void preprocess(
		Mesh* mesh, 
		const bool onVertices, 
		double* origins);

	/**
	 *	write results to mesh
	 *
	 * @param mesh pointer to mesh
	 * @param onVertices if true calculates on vertices, otherwise on facets
	 * @param results array in which results are stored (not normalized)
	 * @param normalize if true normalize from 0 to 1
	 */
	bool postprocess(
		Mesh* mesh, const bool onVertices, double* results, 
		const NormalizeTypes normalize, const int logAlpha,
		const bool smoothing, const int smoothingIterations);

	bool postprocess(
		Mesh* mesh, Mesh::Vertex_handle& v,
		const NormalizeTypes normalize, const int logAlpha,
		const bool smoothing, const int smoothingIterations);

	/**
	 *	calculate shape diameter function
	 *
	 * @param mesh pointer to mesh
	 * @param onVertices if true calculate on vertices, otherwise on facets
	 * @param multiThreaded if true run multi-threaded
	 */
	bool go(
		Mesh* mesh, 
		const bool onVertices, 
		const bool debugMode,
		const bool multiThreaded,
		const int threadsNum,
		const int gridSize,
		const int numCones,
		const float coneSeparationDegree,
		const int raysInCone,
		const bool removeSameDirectionIntersection,
		const bool gaussianWeights,
		const bool removeOutlier,
		const NormalizeTypes normalize,
		const int logAlpha,
		const bool smoothing,
		const int smoothingIterations,
		const bool showDialog = true);



	Aff_transformation_3 rotationMatrixAroundVertex(const Vector_3& v, const double degrees);

	void computeRayIntersect(Mesh mesh, const int gridSize = 40, const bool debugMode = false);

	/**
	 * Check if a given segment intersects with this mesh, relies on the triangle searching data structure
	 * @param segment check if this segment intersects the mesh
	 * @param distance return the distance to the nearest intersection
	 * @return true if intersects, false otherwise
	 */
	bool isIntersectsWithRay(const Ray_3 &ray, double* distance, FILE* debugFile = NULL, const bool debugMode = false);
	
	/** Given a vertex, calculate a segment starting from the vertex and going in a direction opposite to the normal */
	Ray_3 getNormalOpposite(const Mesh::Vertex_handle& v);

	/** Given a point and a normal, calculate a ray from that point, going opposite to the normal */
	Ray_3 getNormalOpposite(const Point_3& p, const Vector_3& n);

	/** given a vertex, calculate a cone of rays to shoot out */
	std::list<Ray_3> getConeRays(const Mesh::Vertex_handle& v, const float degree, const int spread);

	/** given a point and a normal, calculate a cone of rays to shoot out */
	std::list<Ray_3> getConeRays(const Point_3& p, const Vector_3& n, const float degree, const int spread);

	Ray_3 getParallelRay(const Ray_3& targetRay, const Point_3& p, const Vector_3& n, const int numCones, const float coneSeperation, const int raysInCone);


	/**
	 * Compute approximate volume for all vertices in mesh                                                                     
	 */
	/*
	void computeVolumeSDF(
		int gridSize, 
		int numCones, 
		float coneSeparationDegree,
		double coneSeparationRadian, 
		int raysInCone, 
		bool gaussianWeights,
		bool removeOutlier,
		const bool debugMode = false);
	*/

	/**
	 * Compute approximate volume for specified vertex                                                                     
	 */
	/*
	void computeVolumeSDF(
		Mesh::Vertex_handle& v,
		int gridSize, 
		int numCones,
		float coneSeparationDegree,
		double coneSeparationRadian,
		int raysInCone,
		bool gaussianWeights,
		bool removeOutlier,
		const bool debugMode, 
		std::list<std::pair<Ray_3, double> >* resultRays = NULL,
		double* debugMedian = NULL,
		double* debugRange = NULL);
	*/

	/**
	 *	Compute approximate volume for a point and its normal
	 *	@return the calculated volume
	 */
	/*
	double computeVolumeSDF(
		const Point_3& p,
		const Vector_3& n,
		int gridSize, 
		int numCones,
		float coneSeparationDegree,
		double coneSeparationRadian,
		int raysInCone,
		bool gaussianWeights,
		bool removeOutlier,
		const bool debugMode,
		std::list<std::pair<Ray_3, double> >* resultRays = NULL,
		double* debugMedian = NULL,
		double* debugRange = NULL,
		FILE* debugFile = NULL);
	*/

	double computeVolumeSDF(
		const Point_3& p, const Vector_3& n,
		const int& numCones, const float& coneSeparationDegree,	const double& coneSeparationRadian, const int& raysInCone, 
		const bool& removeSameDirectionIntersection, const bool& gaussianWeights, const bool& removeOutlier,
		const bool& debugMode = false, std::list<std::pair<Ray_3, double> >* resultRays = NULL,
		double* debugMedian = NULL, double* debugRange = NULL, FILE* debugFile = NULL);

	double computeVolumeSDF(
		const Point_3& p, const Vector_3& n, const int& gridSize, 
		const int& numCones, const float& coneSeparationDegree, const double& coneSeparationRadian, const int& raysInCone, 
		const bool& removeSameDirectionIntersection, const bool& gaussianWeights, const bool& removeOutlier,
		const TNT::Array3D<RayIntersectCell*>& grid, const triangleVector_t& triangleVector,
		const double& xmin, const double& xmax, const double& ymin, const double& ymax, const double& zmin, const double& zmax,
		const double& xspread, const double& yspread, const double& zspread,
		const bool& debugMode = false, std::list<std::pair<Ray_3, double> >* resultRays = NULL,
		double* debugMedian = NULL, double* debugRange = NULL, FILE* debugFile = NULL);

};


#endif