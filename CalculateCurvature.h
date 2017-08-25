#ifndef __CALCULATE_CURVATURE_H_
#define __CALCULATE_CURVATURE_H_

#include "stdafx.h"

#include <qobject.h>
#include <qthread.h>

class CalculateCurvatureThread : public QThread
{
	Q_OBJECT;
private:
	Mesh* m_mesh;

	double* m_results;
	const double* m_origins;
	int m_size;

public:
	/**
	 *	constructor which accepts the target vertices/facets
	 */
	CalculateCurvatureThread(
		Mesh* mesh/*,
		const double* origins,
		const int size,
		double* results,
		int numCones,
		float coneSeperation,
		int raysInCone,
		bool gaussianWeights,
		const QString& name = QString::null*/);

	/**
	 *	runs the thread
	 */
	virtual void run();
};

class CalculateCurvature : public QObject 
{
	Q_OBJECT
private:
	
public:
	/**
	 *
	 */
	bool go(
		Mesh* mesh, 
		/*const bool onVertices, */
		const bool multiThreaded/*,
		const int numCones = 3,
		const float coneSeperation = 20,
		const int raysInCone = 4,
		const bool gaussianWeights = true,
		const NormalizeTypes normalize = None,
		const bool smoothing = false,
		const bool smoothingAnisotropic = false,
		const int smoothingIterations = 1,
		const bool showDialog = true*/);
};

#endif