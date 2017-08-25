#ifndef __DEFORMATION_ANALYSIS_H_
#define __DEFORMATION_ANALYSIS_H_

#include "stdafx.h"

#include <qobject.h>
#include <qthread.h>

class CalculateDeformationGradientThread : public QThread 
{
	Q_OBJECT;
private:
	QString m_name; // debug file

	QString m_dirPath; // working directory path

	Mesh* m_mesh;

	float m_deformationGradient[3][3]; // deformation gradient matrix

	//double* m_results;
	//const double* m_origins;
	//int m_size;
	//int m_numCones;
	//float m_coneSeperation;
	//int m_raysInCone;
	//bool m_gaussianWeights;

signals:
	/// number of origins completed since the last singal
	void advancedIn(int numOfOrigins);

public:
	/**
	 *	constructor which accepts the target vertices/facets
	 */

	CalculateDeformationGradientThread(
		QString dirPath,
		Mesh* mesh,
		float deformationGradient[3][3],
		const QString& name = QString::null);

	//float[][] getOrientationMatrix(const Mesh::Facet_const_handle& f);

	/*
	CalculateSDFThread(
		Mesh* mesh,
		const double* origins,
		const int size,
		double* results,
		int numCones,
		float coneSeperation,
		int raysInCone,
		bool gaussianWeights,
		const QString& name = QString::null);
	*/

	/**
	 *	runs the thread
	 */
	virtual void run();

};

#endif