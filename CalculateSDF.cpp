#include "stdafx.h"

#include "CalculateSDF.h"
#include "RayIntersect.h"
#include "WorldManager.h"

#include <qdatetime.h>
#include <QApplication>
#include <QProgressDialog>

const double reducedWeight = 0.2;
const double windowSize = 0.1;

double CalculateSDF::logit(const double& value, const int& range)
{
	return logf(value * range + 1) / logf(range + 1);
}

void CalculateSDF::preprocess(
	Mesh* mesh, 
	const bool onVertices, 
	double* origins)
{
	if (onVertices)
	{
		int ind = 0;
		Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			origins[ind*6  ] = vit->point().x();
			origins[ind*6+1] = vit->point().y();
			origins[ind*6+2] = vit->point().z();
			origins[ind*6+3] = vit->normal().x();
			origins[ind*6+4] = vit->normal().y();
			origins[ind*6+5] = vit->normal().z();
			ind++;
		}
	}
	else
	{
		int ind = 0;
		Mesh::Facet_const_iterator fit = mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			Point_3 p = mesh->getFacetCenter(fit);
			origins[ind*6  ] = p.x();
			origins[ind*6+1] = p.y();
			origins[ind*6+2] = p.z();
			origins[ind*6+3] = fit->normal().x();
			origins[ind*6+4] = fit->normal().y();
			origins[ind*6+5] = fit->normal().z();
			ind++;
		}
	}
}


bool CalculateSDF::postprocess(
	Mesh* mesh, Mesh::Vertex_handle& v,
	const NormalizeTypes normalizeTypes, const int logAlpha,
	const bool smoothing, const int smoothingIterations)
{
	const int numPoints = mesh->size_of_vertices();

	//////////////////////////////////////////////////////////////////////////
	// Find the missing values, get average of neighbors
	//////////////////////////////////////////////////////////////////////////
	if (v->volumeSDF() == 0.0)
	{
		double smoothedValue = 0.0;
		double smoothedWeights = 0.0;
	
		Mesh::Halfedge_around_vertex_circulator c = v->vertex_begin();
		do {
			smoothedValue += c->opposite()->vertex()->volumeSDF();
			smoothedWeights += 1.0;
		} while (++c != v->vertex_begin());

		v->volumeSDF() = safeDiv(smoothedValue, smoothedWeights);
	}


	//////////////////////////////////////////////////////////////////////////
	// Find max and min value
	//////////////////////////////////////////////////////////////////////////
	const double minValue = mesh->getMinVolume();
	const double maxValue = mesh->getMaxVolume();
	const double avgValue = mesh->getAvgVolume();

	if (minValue == maxValue)
		return false; // something is wrong. abort all.	

	if (normalizeTypes != None)
	{		
		//const double theDivider = 1 / (maxValue - minValue);
		const double theDivider = 1 / (avgValue - minValue);
	
		v->volumeNSDF() = (v->volumeSDF() - minValue) * theDivider;
		
		if (normalizeTypes == Log)
			v->volumeNSDF() = logit(v->volumeNSDF(), logAlpha);
	}

/*
	//////////////////////////////////////////////////////////////////////////
	// Place the values in the mesh
	//////////////////////////////////////////////////////////////////////////
	if (smoothing)
	{
		double smoothedValue = v->volumeSDF();
		double smoothedWeights = (v->volumeSDF() == 0.0 ? 0.0 : 1.0);				
		Mesh::Halfedge_around_vertex_circulator c = v->vertex_begin();
		do {
			smoothedValue += (c->opposite()->vertex()->volumeSDF()) * reducedWeight;
			smoothedWeights += reducedWeight;
		} while (++c != v->vertex_begin());

		v->volumeSDF() = safeDiv(smoothedValue, smoothedWeights);
	}
*/

	return true;
}


bool CalculateSDF::postprocess(
	Mesh* mesh, const bool onVertices, double* results, 
	const NormalizeTypes normalizeTypes, const int logAlpha,
	const bool smoothing, const int smoothingIterations)
{
	const int numPoints = (onVertices ? mesh->size_of_vertices() : mesh->size_of_facets());

	//////////////////////////////////////////////////////////////////////////
	// Find the missing values, get average of neighbors
	//////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < numPoints; i++)
	{
		if (results[i] == 0.0) // missing values
		{
			double smoothedValue = 0.0;
			double smoothedWeights = 0.0;

			if (onVertices) 
			{
				Mesh::Vertex_handle v = mesh->findVertex(i);
				Mesh::Halfedge_around_vertex_circulator c = v->vertex_begin();
				do {
					if (c->opposite()->vertex() != NULL)
					{
						smoothedValue += results[c->opposite()->vertex()->index()];
						smoothedWeights += 1.0;
					}
				} while (++c != v->vertex_begin());
			} 
			else 
			{
				Mesh::Facet_handle f = mesh->findFacet(i);
				Mesh::Halfedge_around_facet_circulator c = f->facet_begin();
				do {
					if (c->opposite()->facet() != NULL)
					{
						smoothedValue += results[c->opposite()->facet()->index()];
						smoothedWeights += 1.0;
					}
				} while (++c != f->facet_begin());
			}

			results[i] = safeDiv(smoothedValue, smoothedWeights);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Find max and min value
	//////////////////////////////////////////////////////////////////////////
	double minValue = FLT_MAX;
	double maxValue = -FLT_MAX;
	double avgValue = 0.0;

	int ind = 0;
	for (Mesh::Vertex_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++)
	{
		if (results[ind] < minValue) minValue = results[ind];
		if (results[ind] > maxValue) maxValue = results[ind];
		avgValue += results[ind];

		vit->volumeSDF() = results[ind];

		ind++;
	}

	avgValue = avgValue / numPoints;
	
	mesh->setMinVolume(minValue);
	mesh->setMaxVolume(maxValue);
	mesh->setAvgVolume(avgValue);

	if (minValue == maxValue)
		return false; // something is wrong. abort all.

	/// calculate standard deviation ///
	double stdDevValue = 0.0;
	for (Mesh::Vertex_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++)
	{
		stdDevValue += pow(vit->volumeSDF() - avgValue, 2);
	}
	stdDevValue = sqrt(stdDevValue / numPoints);

	double minStdDevValue = FLT_MAX;
	double maxStdDevValue = -FLT_MAX;
	double avgStdDevValue = 0.0;
	int stdDevCounter = 0;
	for (Mesh::Vertex_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++)
	{
		double sdf = vit->volumeSDF();

		if (fabs(sdf - avgValue) <= 2 * stdDevValue)
		{
			if (sdf < minStdDevValue) minStdDevValue = sdf;
			if (sdf > maxStdDevValue) maxStdDevValue = sdf;
			avgStdDevValue += sdf;
			stdDevCounter++;
		}
	}

	avgStdDevValue = avgStdDevValue / stdDevCounter;

	mesh->setMinStdDevVolume(minStdDevValue);
	mesh->setMaxStdDevVolume(maxStdDevValue);
	mesh->setAvgStdDevVolume(avgStdDevValue);

	if (normalizeTypes != None)
	{		
		//double theDivider = 1 / (maxValue - minValue);
		//const double theDivider = 1 / (avgValue - minValue);
		const double theDivider = 1 / (avgStdDevValue - minStdDevValue);
		
		for (Mesh::Vertex_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++)
		{
			//vit->volumeNSDF() = (vit->volumeSDF() - minValue) * theDivider;
			vit->volumeNSDF() = (vit->volumeSDF() - minStdDevValue) * theDivider;
			
			if (normalizeTypes == Log)
				vit->volumeNSDF() = logit(vit->volumeNSDF(), logAlpha);
		}
	}

	double minNValue = FLT_MAX;
	double maxNValue = -FLT_MAX;
	double avgNValue = 0.0;
	double minStdDevNValue = FLT_MAX;
	double maxStdDevNValue = -FLT_MAX;
	double avgStdDevNValue = 0.0;
	stdDevCounter = 0;
	for (Mesh::Vertex_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++)
	{
		if (vit->volumeNSDF() < minNValue) minNValue = vit->volumeNSDF();
		if (vit->volumeNSDF() > maxNValue) maxNValue = vit->volumeNSDF();
		avgNValue += vit->volumeNSDF();

		if (fabs(vit->volumeSDF() - avgValue) <= 2 * stdDevValue)
		{
			if (vit->volumeNSDF() < minStdDevNValue) minStdDevNValue = vit->volumeNSDF();
			if (vit->volumeNSDF() > maxStdDevNValue) maxStdDevNValue = vit->volumeNSDF();
			avgStdDevNValue += vit->volumeNSDF();
			stdDevCounter++;
		}
	}

	avgNValue = avgNValue / numPoints;
	avgStdDevNValue = avgStdDevNValue / stdDevCounter;

	mesh->setMinNVolume(minNValue);
	mesh->setMaxNVolume(maxNValue);
	mesh->setAvgNVolume(avgNValue);

	mesh->setMinStdDevNVolume(minStdDevNValue);
	mesh->setMaxStdDevNVolume(maxStdDevNValue);
	mesh->setAvgStdDevNVolume(avgStdDevNValue);

	cout.precision(10);
	cout << "minSDF = " << minValue << endl;
	cout << "maxSDF = " << maxValue << endl;
	cout << "avgSDF = " << avgValue << endl;
	cout << "minNSDF = " << minNValue << endl;
	cout << "maxNSDF = " << maxNValue << endl;
	cout << "avgNSDF = " << avgNValue << endl;

	cout << "minStdDevValue = " << mesh->getMinStdDevVolume() << endl;
	cout << "maxStdDevValue = " << mesh->getMaxStdDevVolume() << endl;
	cout << "avgStdDevValue = " << mesh->getAvgStdDevVolume() << endl;
	cout << "minStdDevNValue = " << mesh->getMinStdDevNVolume() << endl;
	cout << "maxStdDevNValue = " << mesh->getMaxStdDevNVolume() << endl;
	cout << "avgStdDevNValue = " << mesh->getAvgStdDevNVolume() << endl;

/*
	//////////////////////////////////////////////////////////////////////////
	// Place the values in the mesh
	//////////////////////////////////////////////////////////////////////////
	if (onVertices)
	{
		int ind = 0;
		Mesh::Vertex_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			if (!smoothing)
			{
				vit->volumeSDF(results[ind]);
			}
			else
			{
				double smoothedValue = results[ind];
				double smoothedWeights = (results[ind] == 0.0 ? 0.0 : 1.0);				
				Mesh::Halfedge_around_vertex_circulator c = vit->vertex_begin();
				do {
					smoothedValue += results[c->opposite()->vertex()->index()] * reducedWeight;
					smoothedWeights += reducedWeight;
				} while (++c != vit->vertex_begin());

				vit->volumeSDF() = safeDiv(smoothedValue, smoothedWeights);
			}

			ind++;			
		}
	}
	else
	{
		int ind = 0;
		Mesh::Facet_iterator fit = mesh->facets_begin();
		Mesh::Facet_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			if (!smoothing)
			{
				fit->volumeSDF(results[ind]);
			}
			else
			{
				double initialValue = results[ind];
				for (int i = 0; i < smoothingIterations; i++)
				{
					double smoothedValue = initialValue;
					double smoothedWeights = (initialValue == 0.0 ? 0.0 : 1.0);
					Mesh::Halfedge_around_facet_circulator c = fit->facet_begin();
					do {
						if (c->opposite()->facet() != NULL)
						{
							double sneighbor = results[c->opposite()->facet()->index()];
							//if (fabs(sneighbor-results[ind])<window_size) {
								smoothedValue += results[c->opposite()->facet()->index()] * reducedWeight;
								smoothedWeights += reducedWeight;
							//}						
						}					
					} while (++c != fit->facet_begin());

					fit->volumeSDF() = safeDiv(smoothedValue, smoothedWeights);
					initialValue = fit->volumeSDF();
				}				
			}

			ind++;
		}
	}
*/
	return true;
}


bool CalculateSDF::go(
	Mesh* mesh, const bool onVertices, const bool debugMode,
	const bool multiThreaded, const int threadsNum,
	const int gridSize, const int numCones, const float coneSeparationDegree, const int raysInCone,
	const bool removeSameDirectionIntersection, const bool gaussianWeights, const bool removeOutlier, 
	const NormalizeTypes normalizeTypes, const int logAlpha,
	const bool smoothing, const int smoothingIterations, const bool showDialog)
{
	QTime timer;

	timer.start();

	const int originsSize = (onVertices ? mesh->size_of_vertices() : mesh->size_of_facets());

	double* origins = new double[originsSize * 6];

	preprocess(mesh, onVertices, origins);
	
	const int numOfThreads = (multiThreaded? threadsNum : 1);

	if (multiThreaded)
		printf("Work on multi-thread with %d threads\n", threadsNum);
	else
		printf("Work on single thread\n");

	int timerPreprocess = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Prepare search structure on mesh (Initialize RayIntersect)
	//////////////////////////////////////////////////////////////////////////
	if (debugMode) cout << "Going to compute ray intersect...";	

	TNT::Array3D<RayIntersectCell*> grid;
	triangleVector_t triangleVector;
	double xmin, xmax, ymin, ymax, zmin, zmax, xspread, yspread, zspread, diagonal;

	RayIntersect rayIntersect;
	rayIntersect.Init(*mesh, grid, triangleVector,
						xmin, xmax, ymin, ymax, zmin, zmax,
						xspread, yspread, zspread, 
						gridSize, debugMode);

	if (debugMode) cout << "Done!" << endl;

	int timerComputeSearchStructure = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Initialize the threads
	//////////////////////////////////////////////////////////////////////////
	const double coneSeparationRadian = coneSeparationDegree * DEG2RAD;

	cout << "Calculate SDF...";

	double* results = new double[originsSize];

	#pragma omp parallel num_threads(threadsNum)  /// multi-thread using OpenMP ///
	{
		#pragma omp for
		
		for (int i = 0; i < originsSize; i++)
		{
			Point_3 p(origins[i*6], origins[i*6+1], origins[i*6+2]);
			Vector_3 n(origins[i*6+3], origins[i*6+4], origins[i*6+5]);

			results[i] = computeVolumeSDF(
							p, n,  gridSize, 
							numCones, coneSeparationDegree, coneSeparationRadian, raysInCone, 
							removeSameDirectionIntersection, gaussianWeights, removeOutlier, 
							grid, triangleVector,
							xmin, xmax, ymin, ymax, zmin, zmax, 
							xspread, yspread, zspread,
							debugMode);
		}

	}

	printf("Done!\n");

	int timerAllThreadsFinished = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Put the results back in the mesh (only apply smoothing if on vertices)
	//////////////////////////////////////////////////////////////////////////
	printf("normalizeTypes = ");
	if (normalizeTypes == None)
		printf("None\n");
	else if (normalizeTypes == MinMax)
		printf("MinMax\n");
	else if (normalizeTypes == Log)
		printf("Log, alpha = %d\n", logAlpha);

	if (!postprocess(mesh, onVertices, results, normalizeTypes, logAlpha, smoothing && onVertices, smoothingIterations))
		return false;

	//////////////////////////////////////////////////////////////////////////
	// perform smoothing on the mesh itself, and store the nvolume results back in the mesh
	//////////////////////////////////////////////////////////////////////////
	/*
	if (onVertices)
		mesh->fillNormalizedVertexVolume();
	else
		mesh->fillNormalizedFacetVolume();
	*/

	/*
	if (onVertices)
		mesh->fillNormalizedVertexVolume();
	else
		mesh->makeFacesNVolume(smoothing, smoothingAnisotropic, smoothingIterations);
	*/

	int timerPostprocess = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Clean up
	//////////////////////////////////////////////////////////////////////////
	delete[] origins;
	delete[] results;

	//delete m_grid
	for (int x = 0; x < grid.dim1(); x++)
		for (int y = 0; y < grid.dim2(); y++)
			for (int z = 0; z < grid.dim3(); z++)
				delete grid[x][y][z];	

	for (int i = 0; i < triangleVector.size(); i++)
		delete triangleVector[i];

	triangleVector.clear();
	
	//////////////////////////////////////////////////////////////////////////

	int timerTotal = timer.elapsed();

	/*
	double minSDF = FLT_MAX;
	double maxSDF = -FLT_MAX;
	double avgSDF = 0.0;
	double minNSDF = FLT_MAX;
	double maxNSDF = -FLT_MAX;
	double avgNSDF = 0.0;
	std::vector<double> sdfValues(originsSize);
	std::vector<double> nsdfValues(originsSize);
	for (Mesh::Vertex_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++)
	{
		double sdf = vit->volumeSDF();
		double nsdf = vit->volumeNSDF();
		if (sdf < minSDF) minSDF = sdf;
		if (sdf > maxSDF) maxSDF = sdf;
		if (nsdf < minNSDF) minNSDF = nsdf;
		if (nsdf > maxNSDF) maxNSDF = nsdf;

		avgSDF += sdf;
		avgNSDF += nsdf;

		sdfValues[vit->index()] = sdf;
		nsdfValues[vit->index()] = nsdf;
	}
	avgSDF = avgSDF / originsSize;
	avgNSDF = avgNSDF / originsSize;

	std::nth_element(sdfValues.begin(), sdfValues.begin() + (sdfValues.size() / 2), sdfValues.end());
	double medianSDF = sdfValues[originsSize / 2];
	double sdfStandardDeviation = 0.0;
	double nsdfStandardDeviation = 0.0;
	for (int i = 0; i < originsSize; i++)
	{
		sdfStandardDeviation += pow(sdfValues[i] - avgSDF, 2);
		nsdfStandardDeviation += pow (nsdfValues[i] -avgNSDF, 2);
	}

	sdfStandardDeviation = sqrt(sdfStandardDeviation / originsSize);
	nsdfStandardDeviation = sqrt(nsdfStandardDeviation / originsSize);

	double sdfMin = FLT_MAX;
	double sdfMax = -FLT_MAX;
		
	double sdfMedianDistanceCounter = 0.0;
	int sdfMedianIntersectionCounter = 0;
	double sdfMeanDistanceCounter = 0.0;
	int sdfManIntersectionCounter = 0;
	for (int i = 0; i < originsSize; i++)
	{
		if (fabs(sdfValues[i] - medianSDF) <= 2 * sdfStandardDeviation)
		{
			sdfMedianDistanceCounter += sdfValues[i];
			sdfMedianIntersectionCounter++;
		}
		if (fabs(sdfValues[i] - avgSDF) <= 2 * sdfStandardDeviation)
		{
			if (sdfValues[i] < sdfMin) sdfMin = sdfValues[i];
			if (sdfValues[i] > sdfMax) sdfMax = sdfValues[i];
				
			sdfMeanDistanceCounter += sdfValues[i];
			sdfManIntersectionCounter++;
		}
	}
	sdfMedianDistanceCounter = sdfMedianDistanceCounter / sdfMedianIntersectionCounter;
	sdfMeanDistanceCounter = sdfMeanDistanceCounter / sdfManIntersectionCounter;

	cout.precision(10);
	cout << "minSDF = " << minSDF << endl;
	cout << "maxSDF = " << maxSDF << endl;
	cout << "avgSDF = " << avgSDF << endl;
	cout << "minNSDF = " << minNSDF << endl;
	cout << "maxNSDF = " << maxNSDF << endl;
	cout << "avgNSDF = " << avgNSDF << endl;

	cout << "sdfMedianIntersectionCounter = " << sdfMedianIntersectionCounter << endl;
	cout << "sdfMedianDistanceCounter = " << sdfMedianDistanceCounter << endl;
	cout << "sdfManIntersectionCounter = " << sdfManIntersectionCounter << endl;
	cout << "sdfMeanDistanceCounter = " << sdfMeanDistanceCounter << endl;
	cout << "sdfMin = " << sdfMin << endl;
	cout << "sdfMax = " << sdfMax << endl;
	
	printf("Total time %f seconds\n\n", (float)timerTotal / 1000);
	*/
	return true;
}

/*
void CalculateSDF::computeVolumeSDF(
	Mesh* mesh,
	int gridSize, 
	int numCones,
	float coneSeparationDegree,
	double coneSeparationRadian,
	int raysInCone,
	bool gaussianWeights,
	bool removeOutlier,
	const bool debugMode)
{
	TNT::Array3D<RayIntersectCell*> grid;
	triangleVector_t triangleVector;
	double xmin, xmax, ymin, ymax, zmin, zmax, xspread, yspread, zspread, diagonal;

	if (!m_initialized)
	{
		RayIntersect rayIntersect;

		rayIntersect.Init(*mesh, initialized, grid, triangleVector,
							xmin, xmax, ymin, ymax, zmin, zmax,
							xspread, yspread, zspread, diagonal, 
							gridSize, debugMode);

		m_initialized = true;
	}		


	double minVolume = FLT_MAX;
	double maxVolume = -FLT_MAX;
	for (Mesh::Vertex_iterator it = mesh->vertices_begin(); it != mesh->vertices_end(); it++)
	{
		computeVolumeSDF(
			it, gridSize, numCones, coneSeparationDegree, coneSeparationRadian, raysInCone, 
			gaussianWeights, removeOutlier, debugMode);

		if (it->volumeSDF() < m_minVolume) minVolume = it->volumeSDF();
		if (it->volumeSDF() > m_maxVolume) maxVolume = it->volumeSDF();
	}

	mesh->setMinVolume(minVolume);
	mesh->setMaxVolume(maxVolume);

	/// normalize
	double nsdf;
	for (Mesh::Vertex_iterator it = mesh->vertices_begin(); it != mesh->vertices_end(); it++)
	{
		nsdf = (it->volumeSDF() - minVolume) / (maxVolume - minVolume);
		it->volumeNSDF(nsdf);
	}

	minNVolume = 0;
	maxNVolume = 1;

	mesh->setMinNVolume(minNVolume);
	mesh->setMaxNVolume(maxNVolume);
}


void CalculateSDF::computeVolumeSDF(
	Vertex_handle& v,
	int gridSize,
	int numCones,
	float coneSeparationDegree,
	double coneSeparationRadian,
	int raysInCone,
	bool gaussianWeights,
	bool removeOutlier,
	const bool debugMode, 
	std::list<std::pair<Ray_3, double> >* resultRays,
	double* debugMedian,
	double* debugRange)
{
	assert(m_initialized);

	double vsdf = computeVolumeSDF(
							v->point(), v->normal(), gridSize, numCones, 
							coneSeparationDegree, coneSeparationRadian, raysInCone, gaussianWeights, removeOutlier,
							debugMode, resultRays, debugMedian, debugRange);

	v->volumeSDF(vsdf);
}*/


double CalculateSDF::computeVolumeSDF(
	const Point_3& p, const Vector_3& n, const int& gridSize, 
	const int& numCones, const float& coneSeparationDegree, const double& coneSeparationRadian, const int& raysInCone, 
	const bool& removeSameDirectionIntersection, const bool& gaussianWeights, const bool& removeOutlier,
	const TNT::Array3D<RayIntersectCell*>& grid, const triangleVector_t& triangleVector,
	const double& xmin, const double& xmax, const double& ymin, const double& ymax, const double& zmin, const double& zmax,
	const double& xspread, const double& yspread, const double& zspread,
	const bool& debugMode, std::list<std::pair<Ray_3, double> >* resultRays,
	double* debugMedian, double* debugRange, FILE* debugFile)
{
	if (debugFile) { fprintf(debugFile, "Beginning compute volume\n"); fflush(debugFile); }
	if (debugMode) { cout << "Beginning compute volume" << endl; }

	//basically shoot rays opposite the normal and do something with the distances we receive

	double result = 0.0;

	Ray_3 ray = getNormalOpposite(p, n);

	double distanceCounter = 0.0;
	float intersectionCounter = 0.0;

	double distance = FLT_MAX;

	std::vector<double> values;
	double maxValue = 0.0;
	std::vector<float> weights;
	values.reserve(numCones * raysInCone + 1); // Request a change in capacity

	if (debugFile) { fprintf(debugFile, "Found normal opposite\n"); fflush(debugFile);}
	if (debugMode) { cout << "Found normal opposite" << endl; }

	// first check straight opposite to normal //
	RayIntersect rayIntersect;

	distance = rayIntersect.Query(
		ray, /*m_initialized,*/ grid, triangleVector, 
		xmin, ymin, zmin, xspread, yspread, zspread, removeSameDirectionIntersection,
		NULL, NULL, debugFile, debugMode);

	if (distance != FLT_MAX)
	{
		//intersectionCounter+= 1.0;
		//distanceCounter += distance;

		if (distance > maxValue) 
			maxValue = distance;

		values.push_back(distance);
		weights.push_back(1.0);

		if (resultRays)
			resultRays->push_back(std::pair<Ray_3, double>(ray, distance));
	}
	else if (resultRays)
	{
		resultRays->push_back(std::pair<Ray_3, double>(ray, FLT_MAX/*diagonalLength()*/));
	}

	if (debugFile) { fprintf(debugFile, "Going to compute cones\n"); fflush(debugFile); }
	if (debugMode) { cout << "Going to compute cones" << endl; }

	//const double gaussianVar = 120.0 * DEG2RAD;
	const double gaussianVar = numCones * coneSeparationDegree * DEG2RAD;
	const double gaussianDiv2 = 1 / (2 * pow(gaussianVar,2));

	for (int i = 0; i < numCones; i++)
	{
		//now check cone
		double radianDegree = coneSeparationRadian * (i+1);
		std::list<Ray_3> list = getConeRays(p, n, radianDegree, raysInCone);

		for (std::list<Ray_3>::iterator it = list.begin(); it != list.end(); it++)
		{
			distance = rayIntersect.Query(
				*it, /*m_initialized,*/ grid, triangleVector, 
				xmin, ymin, zmin, xspread, yspread, zspread, removeSameDirectionIntersection,
				NULL, NULL, debugFile, debugMode);

			if (distance != FLT_MAX)
			{
				double weight;
				if (gaussianWeights)
				{
					//weight = (double) 1 / cos(coneSeparationRadian);
					weight = expf(-pow(radianDegree,2) * gaussianDiv2);
					//weight = gaussianDiv1 * expf(-pow(radianDegree,2) * gaussianDiv2);

					weight = 1 / (1 - weight);
				}
				else
				{
					weight = 1.0;
				}

				//intersectionCounter += weight;
				//distanceCounter += weight * distance;

				if (distance > maxValue) 
					maxValue = distance;

				values.push_back(distance);
				weights.push_back(weight);

				if (resultRays)
					resultRays->push_back(std::pair<Ray_3, double>(*it, distance));
			}
		}
	}

	//calculate the average
	for (int i = 0; i < values.size(); i++)
	{
		distanceCounter += values[i] * weights[i];
		//intersectionCounter += weights[i];
		intersectionCounter += 1.0;
	}

	if (intersectionCounter)
	{
		if (!removeOutlier)
			result = distanceCounter / intersectionCounter;
		else
		{
			result = distanceCounter / intersectionCounter;

			//fix really insane values
			std::nth_element(values.begin(), values.begin()+(values.size()/2), values.end());
			double median = values[values.size()/2];

			float sdfStandardDeviation = 0.0;
			for (int i=0; i<values.size();i++)
				sdfStandardDeviation += pow(values[i] - result, 2);

			sdfStandardDeviation = sqrt(sdfStandardDeviation / values.size());
			distanceCounter = intersectionCounter = 0.0;

			if (debugMedian) *debugMedian = median;
			if (debugRange) *debugRange = 0.75 * sdfStandardDeviation;

			for (int i=0; i<values.size(); i++)
			{
				//if (fabs(values[i] - result) <= 0.5 * sdfStandardDeviation) {
				if (fabs(values[i] - median) <= 0.5 * sdfStandardDeviation)
				{
					distanceCounter += values[i] * weights[i];
					intersectionCounter += weights[i];
				}
			}

			if (intersectionCounter)
				result = distanceCounter / intersectionCounter;
			else
				result = 0.0;
		}
	}
	else
	{
		result = 0.0;
	}


	return result;
}


/*
void Mesh::computeVolume() 
{
	if (!m_rayIntersect.isInitialized())
		computeRayIntersect();

	m_minVolume = FLT_MAX;
	m_maxVolume = -FLT_MAX;
	int counter = 0;
	for (Mesh::Vertex_iterator it = vertices_begin(); it != vertices_end(); it++)
	{
		compute_volume(it, 3, 20, true);
		//compute_volume(it, m_renderingParams->m_volumeConeNumber,
		//	m_renderingParams->m_volumeConeSeperation,
		//	m_renderingParams->m_volumeRaysInCone,
		//	m_renderingParams->m_volumeUseGaussianWeights);

		if (it->volume() < m_minVolume) m_minVolume = it->volume();
		if (it->volume() > m_maxVolume) m_maxVolume = it->volume();

		//f << counter << ": " << it->volume() << std::endl;

		counter++;
	}

	//normalize
	float mod = 100;
	float logmod = logf(mod);
	for (Mesh::Vertex_iterator it = vertices_begin(); it != vertices_end(); it++)
	{
		it->volume() = (it->volume() - m_minVolume) / (m_maxVolume - m_minVolume);
	}

	m_minVolume = 0;
	m_maxVolume = 1;
}

void CalculateSDF::computeRayIntersect(Mesh* mesh, const int gridSize, const bool debugMode)
{
	if (!m_rayIntersect.isInitialized())
		m_rayIntersect.Init(*mesh, gridSize, debugMode);
}

bool CalculateSDF::isIntersectsWithRay(const Ray_3 &ray, double* distance, FILE* debugFile, const bool debugMode)
{
	*distance = m_rayIntersect.Query(ray, NULL, NULL, debugFile, debugMode);

	return *distance != FLT_MAX;
}
*/

Ray_3 CalculateSDF::getNormalOpposite(const Point_3& p, const Vector_3& n)
{
	//direction opposite the normal
	Vector_3 direction = n * -1;

	//construct segment from vertex pointing towards direction
	Ray_3 ray(p + direction * 1e-5, direction);

	return ray;
}

Ray_3 CalculateSDF::getNormalOpposite(const Mesh::Vertex_handle& v)
{
	//direction opposite the normal
	Vector_3 direction = v->normal() * -1;

	//construct segment from vertex pointing towards direction
	Ray_3 ray(v->point() + direction * 1e-5, direction);

	return ray;
}

Aff_transformation_3 CalculateSDF::rotationMatrixAroundVertex(const Vector_3& v, const double degrees)
{
	double c = cos(degrees);
	double s = sin(degrees);
	double t = 1-c;
	Enriched_kernel::Aff_transformation_3 aff2(
		t*v.x()*v.x() + c, t*v.x()*v.y() + s*v.z(), t*v.x()*v.z() - s*v.y(), 0,
		t*v.x()*v.y() - s*v.z(), t*v.y()*v.y() + c, t*v.y()*v.z() + s*v.x(), 0,
		t*v.x()*v.z() + s*v.y(), t*v.y()*v.z() - s * v.x(), t*v.z()*v.z() + c, 0);

	return aff2;
}

std::list<Ray_3> CalculateSDF::getConeRays(
	const Point_3& p, const Vector_3& n,
	const float degree, const int spread)
{
	std::list<Ray_3> list;

	Ray_3 ray = getNormalOpposite(p, n);
	Vector_3 rayv = ray.to_vector();
	rayv = rayv / sqrt(rayv.squared_length());

	//now check a cone
	Plane_3 plane(ray.source(), ray.to_vector());
	Vector_3 base = plane.base1() / sqrt(plane.base1().squared_length());

	Enriched_kernel::Aff_transformation_3 aff = rotationMatrixAroundVertex(base, degree);
	Enriched_kernel::Aff_transformation_3 aroundRay = rotationMatrixAroundVertex(rayv, 2 * M_PI / spread);
	Enriched_kernel::Aff_transformation_3 around = rotationMatrixAroundVertex(rayv, degree);

	Vector_3 spinMe = rayv;
	spinMe = spinMe / sqrt(spinMe.squared_length());
	spinMe = spinMe.transform(aff);
	spinMe = spinMe.transform(around); ///Gaussian sphere
	for (int i = 0; i < spread; i++)
	{
		Ray_3 coneRay(p + spinMe * 1e-5, spinMe);

		list.push_back(coneRay);
		spinMe = spinMe.transform(aroundRay);
	}

	return list;
}

std::list<Ray_3> CalculateSDF::getConeRays(const Mesh::Vertex_handle& v, const float degree, const int spread)
{
	return getConeRays(v->point(), v->normal(), degree, spread);
}

Ray_3 CalculateSDF::getParallelRay(
	const Ray_3& targetRay, 
	const Point_3& p, 
	const Vector_3& n, 
	const int numCones, 
	const float coneSeperation, 
	const int raysInCone)
{
	Vector_3 direction = targetRay.to_vector();
	double targetRay_x = direction.x();
	double targetRay_y = direction.y();
	double targetRay_z = direction.z();

	Ray_3 ray = getNormalOpposite(p, n);
	Vector_3 ray_direction = ray.to_vector();
	double ray_x = ray_direction.x();
	double ray_y = ray_direction.y();
	double ray_z = ray_direction.z();

	if ((targetRay_x * ray_y - targetRay_y * ray_x) == 0)
		if ((targetRay_y * ray_z - targetRay_z * ray_y) == 0)
			if ((targetRay_x * ray_z - targetRay_z * ray_x) == 0)
				return ray;
	
	for (int i=0; i<numCones; i++)
	{
		//now check cone
		float degree = coneSeperation * (i+1);
		std::list<Ray_3> list = getConeRays(p, n, degree, raysInCone);

		for (std::list<Ray_3>::iterator it = list.begin(); it != list.end(); it++)
		{
			ray = *it;
			ray_direction = ray.to_vector();
			ray_x = ray_direction.x();
			ray_y = ray_direction.y();
			ray_z = ray_direction.z();
			
			if ((targetRay_x * ray_y - targetRay_y * ray_x) == 0)
				if ((targetRay_y * ray_z - targetRay_z * ray_y) == 0)
					if ((targetRay_x * ray_z - targetRay_z * ray_x) == 0)
						return ray;
		}
	}
}