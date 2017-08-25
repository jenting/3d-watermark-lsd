#ifndef __WATERMARK_H_
#define __WATERMARK_H_

#include "stdafx.h"

#include <cmath>

#include <qobject.h>

enum WatermarkType {SDF, Laplacian};
enum SDFWatermarkType {Mean, Relation};
enum SearchSpace {MIN_DIST, AVG_DIST, MAX_DIST};

#include "AppManager.h"
#include "AppObject.h"
#include "Exporter.h"
#include "Importer.h"
#include "MyMatrix.h"
#include "Mesh/Polyhedron_Copy.h"
#include "Metric.h"
#include "Sorting.h"
#include "Utility.h"


class PseudoRandomSequenceGenerator
{

public:

	/** psuedo-random generator to generator a watermark sequence with length = totalWatermarkBits and seed */
	int* generate(long seed, int totalWatermarkBits);
	void generate(long seed, TNT::Array1D<int>& watermarkSequence);
};


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


class SDFWatermarkingAlgorithm
{

private:

	/** Calculate the Euclidean distance between two point*/
	double distance(Point_3 p1, Point_3 p2);

	/** Generate a pseudo-random number with range [lower, upper] */
	double randomNumber(double lower, double upper);

	/** Generate a pesudo-random point in a ball with radius and center */
	Point_3 randomPointInSphere(Point_3 center, double radius);

	/**  */
	void getMinMaxSDF(
		const Mesh* mesh, const bool onVertices, 
		double& minSDF, double& maxSDF);

public:

	SDFWatermarkingAlgorithm();
	virtual ~SDFWatermarkingAlgorithm();

	bool embed(
		AppObject* appObject, Mesh* originalMesh, Mesh* watermarkedMesh, 
		const bool onVertices, const bool multiThreaded, const int threadsNum,
		const long seed, const int totalWatermarkBits, const int robustnessFactor, const SearchSpace searchSpace,
		const int maxCycles, const int particleNumbers,
		const int gridSize, const int numCones, const float coneSeparationDegree, const int raysInCone, 
		const bool removeSameDirectionIntersection, const bool gaussianWeights, const bool removeOutlier,
		const NormalizeTypes normalizeTypes, const int logAlpha, 
		const bool smoothing, const int smoothingIterations,
		const bool debugMode, const bool showDialog);

	bool meshEmbeddingProcedure(
		AppObject* appObject, Mesh* watermarkedMesh, bool onVertices,
		const bool multiThreaded, const int threadsNum,
		const int* originalWatermarkSequence, const long seed, const int totalWatermarkBits, const int robustnessFactor, const SearchSpace searchSpace,
		const int maxCycles, const int particleNumbers, const int numBins,
		const int gridSize, const int numCones, const float coneSeparationDegree, const double coneSeparationRadian, const int raysInCone, 
		const bool removeSameDirectionIntersection, const bool gaussianWeights, const bool removeOutlier, 
		const NormalizeTypes normalizeTypes, const int logAlpha,
		const bool smoothing, const int smoothingIterations, 
		const bool debugMode, const bool showDialog);

	bool detect(
		AppObject* appObject, Mesh* watermarkedMesh, const bool onVertices, 
		const long seed, const int totalWatermarkBits,
		const int gridSize, const int numCones, const float coneSeparationDegree, const int raysInCone, 
		const bool gaussianWeights, const bool removeOutlier, const NormalizeTypes normalizeTypes, const int logAlpha,
		const bool smoothing, const int smoothingIterations,
		const bool debugMode, const bool showDialog = true);

	void meshDetectigProcedure(
		Mesh* watermarkedMesh, const int interval,
		vector<vector<int> > vertexIndexBins, TNT::Array1D<int>& detectedWatermarkSequence);

};


#endif