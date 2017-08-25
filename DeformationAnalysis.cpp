#include "stdafx.h"

#include "DeformationAnalysis.h"

const int DEFORMATION_THREADS_NUM = 4;

/*
CalculateDeformationGradientThread::CalculateDeformationGradientThread(
		const QString& dirPath, Mesh* mesh, float deformationGradient[3][3], const bool multiThreaded, const QString& name) :
		m_dirPath(dirPath), m_mesh(mesh), m_deformationGradient(deformationGradient), m_name(name) {}
*/
/*
bool CalculateDeformationGradientThread::go(const QString& dirPath, Mesh* mesh, float deformationGradient[3][3], const bool multiThreaded)
{
	int originsSize = mesh->size_of_facets();
	int numOfThreads = (multiThreaded? DEFORMATION_THREADS_NUM : 1);
	int dataSize = originsSize / numOfThreads;
	std::vector<CalculateSDFThread*> threads;

	for (int i=0; i<numOfThreads; i++)
	{
		if (i == numOfThreads-1 && originsSize % dataSize != 0)
			dataSize += originsSize % dataSize;

		CalculateDeformationGradientThread* calculateDeformationGradientThread =
			new CalculateDeformationGradientThread(dirPath, mesh, deformationGradient);

		//connect(cst, SIGNAL(advancedIn(int)), this, SLOT(addToProgress(int)));

		//dataPtr += dataSize * 6;
		//resultsPtr += dataSize;
		threads.push_back(calculateDeformationGradientThread);
	}

	// Run the threads and wait for them to finish
	if (multiThreaded) {
		for (int i=0; i<threads.size(); i++) {
			threads[i]->start();
		}

		for (int i=0; i<threads.size(); i++) {
			threads[i]->wait();
		}
	} else {
		threads[0]->run();
	}

	return true;
}
*/

void CalculateDeformationGradientThread::run()
{
	/*
	Mesh::Facet_const_iterator fit = mesh->facets_begin();
	Mesh::Facet_const_iterator fit_end = mesh->facets_end();
	for (; fit != fit_end; fit++) {
		
	}
	*/
}
/*
float[][] CalculateDeformationGradientThread::getOrientationMatrix(const Mesh::Facet_const_handle& f) {
	Mesh::Halfedge_around_facet_const_circulator c = f->facet_begin();

	Point_3 p1 = c->vertex()->point();
	Point_3 p2 = ++c->vertex()->point();
	Point_3 p3 = ++c->vertex()->point();

	Vector_3 v1 = p2-p1;
	Vector_3 v2 = p3-p1;

	Point_3 p4 = p1 + v1*v2 / sqrt(fabs(v1*v2), 0.5);

	Vector_3 v3 = p4-p1;

	float orientationMatrix[3][3];

	orientationMatrix[0][0] = v1.x();
	orientationMatrix[1][0] = v1.y();
	orientationMatrix[2][0] = v1.z();

	orientationMatrix[0][1] = v2.x();
	orientationMatrix[1][1] = v2.y();
	orientationMatrix[2][1] = v2.z();

	orientationMatrix[0][2] = v3.x();
	orientationMatrix[1][2] = v3.y();
	orientationMatrix[2][2] = v3.z();

	return orientationMatrix;
}
*/