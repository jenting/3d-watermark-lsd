#include "stdafx.h"

#include "CalculateCurvature.h"

const int CURVATURE_THREADS_NUM = 4;

bool CalculateCurvature::go(
	Mesh* mesh,
	/*const bool onVertices,*/
	const bool multiThreaded/*,
	const int numCones,
	const float coneSeperation,
	const int raysInCone,
	const bool gaussianWeights,
	const NormalizeTypes normalize,
	const bool smoothing,
	const bool smoothingAnisotropic,
	const int smoothingIterations,
	const bool showDialog*/)
{
	//QTime timer;

	//timer.start();

	//int originsSize = (onVertices?mesh->size_of_vertices() : mesh->size_of_facets());
	//double* origins = new double[originsSize * 6];

	int originsSize = mesh->size_of_vertices();
	double* origins = new double[originsSize];

	//preprocess(mesh, /*onVertices,*/ origins);

	int numOfThreads = (multiThreaded ? CURVATURE_THREADS_NUM : 1);

	//int timerPreprocess = timer.elapsed();


	/*if (showDialog)
		m_progDlg = new QProgressDialog("Calculating SDF", "cancel", 0, originsSize - 1, g_main);
	else
		m_progDlg = NULL;*/

	//////////////////////////////////////////////////////////////////////////
	// Prepare search structure on mesh
	//////////////////////////////////////////////////////////////////////////
	//mesh->computerayIntersectError();

	//int timerComputeSearchStructure = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Initialize the threads
	//////////////////////////////////////////////////////////////////////////
	
	double* results = new double[originsSize];
	double* resultsPtr = results;

	double* dataPtr = origins;
	int dataSize = originsSize / numOfThreads;
	std::vector<CalculateCurvatureThread*> threads;
	for (int i=0; i<numOfThreads; i++)
	{
		/*if (i == numOfThreads-1 && originsSize % dataSize != 0)
			dataSize += originsSize % dataSize;*/

		/*CalculateCurvatureThread* cst =
			new CalculateCurvatureThread(mesh/*, dataPtr, dataSize, resultsPtr, numCones, coneSeperation, raysInCone, gaussianWeights, QString::number(i));*/

		//connect(cst, SIGNAL(advancedIn(int)), this, SLOT(addToProgress(int)));

		/*dataPtr += dataSize * 6;*/
		//resultsPtr += dataSize;

		//threads.push_back(cst);
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

	//int timerAllThreadsFinished = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Put the results back in the mesh
	//////////////////////////////////////////////////////////////////////////
	//only apply smoothing if on vertices
	/*if (!postprocess(mesh, onVertices, results, normalize, smoothing && onVertices))
		return false;*/

	//perform smoothing on the mesh itself, and store the nvolume results back in the mesh
	/*if (onVertices)
		mesh->fillNormalizedVertexVolume();
	else
		mesh->makeFacesNVolume(smoothing, smoothingAnisotropic, smoothingIterations);*/


	//int timerPostprocess = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Clean up
	//////////////////////////////////////////////////////////////////////////
	/*for (int i=0; i<threads.size(); i++) {
		delete threads[i];
	}
	threads.clear();
	delete[] origins;
	delete[] results;
	m_progDlg = NULL;*/

	//int timerTotal = timer.elapsed();

	/*if (showDialog) {
		qDebug("Calculated SDF total time %dms. preprocess %dms, compute %dms, threads %dms, postprocess %dms",
			timerTotal, timerPreprocess, timerComputeSearchStructure, timerAllThreadsFinished, timerPostprocess);

		QString workingon = (onVertices?"vertices" : "facets");
		QString mt = (multiThreaded?"multithreaded" : "not multithreaded");
		printf("Mesh size v=%d f=%d\n", mesh->size_of_vertices(), mesh->size_of_facets());
		cout << "Working on ";
		if(onVertices)
            cout << "vertices";
		else
            cout <<"facets";
		if(multiThreaded)
            cout << " with multithreaded" << endl;
        else
            cout << " with not multithreaded" << endl;
		printf("Total time %f seconds\n", (float)timerTotal / 1000);
	}*/

	return true;
}

void CalculateCurvatureThread::run()
{
	//QTime timer;

	//timer.start();

	int lastReported = 0;
	int i;
	for (i = 0; i < m_size; i++)
	{
		//m_results[i] = m_mesh->average_curvature_around();

		/*if ((i % 20) == 0)
		{
			emit advancedIn(i - lastReported);
			qApp->processEvents();
			sleep(0);
			lastReported = i;
		}*/
	}

	//emit advancedIn(i - lastReported);

	//int totalTime = timer.elapsed();
}