#include "stdafx.h"

#include "DirRunner.h"

#include <QVector>

//#include "ui_DirRunner.h"
/*
DirRunner::DirRunner(QWidget *parent) : QWidget(parent)
{
    //ui->setupUi(this);
	ui.setupUi(this);
}

DirRunner::~DirRunner()
{
    //delete ui;
}

void DirRunner::dirRunnerAvgDDegree(AppObject* appObject, const QString& fileDir)
{
	// 1. calculate orientation matrix of the first frame
	// 2. for 2~n frames
	//     a. calculate orientation matrix of each facet
	//     b. calculate deformation distance of each facet
	//     c. use a facet array to record current max deformation distance of each facet
	// 3. paint this result

	const int order = 3;

	Mesh* referenceMesh = appObject->mesh;

	QString ReferenceFramePath = appObject->inFileName;

	cout << "Reference frame path: " << qPrintable(ReferenceFramePath) << endl;
	
	// reference frame orientation matrix
	int numOfFacets = referenceMesh->size_of_facets();

	QVector<float**> referenceOrientationMatrixs(numOfFacets);

	for(Mesh::Facet_iterator fit = referenceMesh->facets_begin(); fit != referenceMesh->facets_end(); fit++)
	{
		Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();

		Point_3 p1 = c->vertex()->point();
		++c;
		Point_3 p2 = c->vertex()->point();
		++c;
		Point_3 p3 = c->vertex()->point();

		Vector_3 v1 = p2-p1;
		Vector_3 v2 = p3-p1;
		Vector_3 v3 = cross_product(v1, v2);
		v3 = v3 / sqrt(v3.squared_length());
		
		float** orientationMatrix = new float*[order];
		for(int i = 0; i < order; i++)
			orientationMatrix[i] = new float[order];

		float** inverseOrientatoinMatrix = new float*[order];
		for(int i = 0; i < order; i++)
			inverseOrientatoinMatrix[i] = new float[order];

		orientationMatrix[0][0] = v1.x();
		orientationMatrix[1][0] = v1.y();
		orientationMatrix[2][0] = v1.z();

		orientationMatrix[0][1] = v2.x();
		orientationMatrix[1][1] = v2.y();
		orientationMatrix[2][1] = v2.z();

		orientationMatrix[0][2] = v3.x();
		orientationMatrix[1][2] = v3.y();
		orientationMatrix[2][2] = v3.z();

		MatrixInversion(orientationMatrix, order, inverseOrientatoinMatrix);

		int fit_index = fit->index();
		referenceOrientationMatrixs[fit_index] = inverseOrientatoinMatrix;

	    // release memory
		for(int i = 0; i < order; i++)
			delete [] orientationMatrix[i];
		delete [] orientationMatrix;
	}

	printf("Reference frame done!\n");

	Mesh* mesh = new Mesh;

	QVector<float> maxDDegree(numOfFacets);
	for(int i = 0; i < numOfFacets; i++)
		maxDDegree[i] = FLT_MIN;

	// list dir files
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for(QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		QFileInfo fi(filePath);
		if(fi.suffix() != "obj") continue;

		cout << qPrintable(filePath) << ":  ";

		bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);

		if(!success) continue;

		// for each facet, calculate its deformation gradient matrix
		// deformation gradient matrix Dt,i = Ot,i * (Or,i)^-1
		// Or,i : reference frame orientation matrix
		// Ot,i : current frame orientation
		// * : matrix product
		QVector<float**> deformationGradientMatrixs(numOfFacets);

		for(int i = 0; i < numOfFacets; i++)
		{
			float** m = new float*[order];
			for(int j = 0; j < order; j++)
				m[j] = new float[order];
			deformationGradientMatrixs[i] = m;
		}

		for(Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++)
		{
			Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();

			Point_3 p1 = c->vertex()->point();
			++c;
			Point_3 p2 = c->vertex()->point();
			++c;
			Point_3 p3 = c->vertex()->point();

			Vector_3 v1 = p2-p1;
			Vector_3 v2 = p3-p1;
			Vector_3 v3 = cross_product(v1, v2);
			v3 = v3 / sqrt(v3.squared_length());

			float** orientationMatrix = new float*[order];
			for(int i = 0; i < order; i++)
				orientationMatrix[i] = new float[order];

			orientationMatrix[0][0] = v1.x();
			orientationMatrix[1][0] = v1.y();
			orientationMatrix[2][0] = v1.z();

			orientationMatrix[0][1] = v2.x();
			orientationMatrix[1][1] = v2.y();
			orientationMatrix[2][1] = v2.z();

			orientationMatrix[0][2] = v3.x();
			orientationMatrix[1][2] = v3.y();
			orientationMatrix[2][2] = v3.z();

			int fit_index = fit->index();
			float** referenceOrientatoinMatrix = referenceOrientationMatrixs.at(fit_index);
			float** deformationGradient = matrixProduct(orientationMatrix, referenceOrientatoinMatrix, order);

			deformationGradientMatrixs[fit_index] = deformationGradient;
		}
		
		QVector<float> deforDist(numOfFacets);
		for(int i = 0; i < numOfFacets; i++)
			deforDist[i] = 0;

		// sum all the deformation degree of current facet
		for (Mesh::Edge_const_iterator eit = mesh->edges_begin(); eit != mesh->edges_end(); eit++) {
			if (eit->facet() == NULL || eit->opposite()->facet() == NULL)
				continue;

			Mesh::Facet_const_handle f1 = eit->facet();
			Mesh::Facet_const_handle f2 = eit->opposite()->facet();

			int f_index1 = f1->index();
			int f_index2 = f2->index();

			float** m1 = deformationGradientMatrixs.at(f_index1);
			float** m2 = deformationGradientMatrixs.at(f_index2);

			float dd = deformationDistance(m1, m2, order);

			deforDist[f_index1] += dd;
			deforDist[f_index2] += dd;
		}

		 // avg deformation degree
		for(Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++)
		{
			int totalEdge = 0;
			int f_index = fit->index();

			Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();
			do { 
				if (c->facet() == NULL || c->opposite()->facet() == NULL)
					continue;
				totalEdge++; 
			} while (++c != fit->facet_begin());

			deforDist[f_index] /= totalEdge;

			if(deforDist.at(f_index) > maxDDegree.at(f_index))
				maxDDegree[f_index] = deforDist[f_index];
		}

		// release memory
		for(int i = 0; i < numOfFacets; i++)
		{
			float** defor = deformationGradientMatrixs.at(i);
			for(int i = 0; i < order; i++)
				delete [] defor[i];
			delete [] defor;
		}

		mesh->clear();

		cout << "ok!" << endl;
	}
	
	// release memory //
	delete mesh;

	// save the dd to reference frame mesh //
	float maxDD = FLT_MIN;
	float minDD = FLT_MAX;
	for(Mesh::Facet_iterator fit = referenceMesh->facets_begin(); fit != referenceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		float mdd = maxDDegree[f_index];

		fit->deformationDegree(mdd);

		if(mdd > maxDD) maxDD = mdd;
		if(mdd < minDD) minDD = mdd;
	}

	// normalize
	for(Mesh::Facet_iterator fit = referenceMesh->facets_begin(); fit != referenceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		float mdd = fit->deformationDegree();
		float n_mdd = (mdd - minDD) / (maxDD - minDD);
		
		fit->normalizeDeformationDegree(n_mdd);
	}

	// release memory //
	for(int i = 0; i < numOfFacets; i++)
	{
		float** ref = referenceOrientationMatrixs[i];
		for(int i = 0; i < order; i++)
			delete [] ref[i];
		delete [] ref;
	}
}

void AppManager::dirRunnerMaxDDegree(AppObject* appObject, const QString& fileDir)
{
	// 1. calculate orientation matrix of the first frame
	// 2. for 2~n frames
	//     a. calculate orientation matrix of each facet
	//     b. calculate deformation distance of each facet
	//     c. use a facet array to record current max deformation distance of each facet
	// 3. paint this result

	const int order = 3;

	Mesh* referenceMesh = appObject->mesh;

	QString ReferenceFramePath = appObject->inFileName;

	cout << "Reference frame path: " << qPrintable(ReferenceFramePath) << endl;

	// reference frame orientation matrix
	int numOfFacets = referenceMesh->size_of_facets();

	QVector<float**> referenceOrientationMatrixs(numOfFacets);

	for(Mesh::Facet_iterator fit = referenceMesh->facets_begin(); fit != referenceMesh->facets_end(); fit++)
	{
		Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();

		Point_3 p1 = c->vertex()->point();
		++c;
		Point_3 p2 = c->vertex()->point();
		++c;
		Point_3 p3 = c->vertex()->point();

		Vector_3 v1 = p2-p1;
		Vector_3 v2 = p3-p1;
		Vector_3 v3 = cross_product(v1, v2);
		v3 = v3 / sqrt(v3.squared_length());
		
		float** orientationMatrix = new float*[order];
		for(int i = 0; i < order; i++)
			orientationMatrix[i] = new float[order];

		float** inverseOrientatoinMatrix = new float*[order];
		for(int i = 0; i < order; i++)
			inverseOrientatoinMatrix[i] = new float[order];

		orientationMatrix[0][0] = v1.x();
		orientationMatrix[1][0] = v1.y();
		orientationMatrix[2][0] = v1.z();

		orientationMatrix[0][1] = v2.x();
		orientationMatrix[1][1] = v2.y();
		orientationMatrix[2][1] = v2.z();

		orientationMatrix[0][2] = v3.x();
		orientationMatrix[1][2] = v3.y();
		orientationMatrix[2][2] = v3.z();

		MatrixInversion(orientationMatrix, order, inverseOrientatoinMatrix);

		int fit_index = fit->index();
		referenceOrientationMatrixs[fit_index] = inverseOrientatoinMatrix;

	    // release memory
		for(int i = 0; i < order; i++)
			delete [] orientationMatrix[i];
		delete [] orientationMatrix;
	}

	printf("Reference frame done!\n");

	Mesh* mesh = new Mesh;

	QVector<float> maxDDegree(numOfFacets);
	for(int i = 0; i < numOfFacets; i++)
		maxDDegree[i] = FLT_MIN;

	// list dir files
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for(QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		QFileInfo fi(filePath);
		if(fi.suffix() != "obj") continue;

		cout << qPrintable(filePath) << ":  ";

		bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);

		if(!success) continue;

		// for each facet, calculate its deformation gradient matrix
		// deformation gradient matrix Dt,i = Ot,i * (Or,i)^-1
		// Or,i : reference frame orientation matrix
		// Ot,i : current frame orientation
		// * : matrix product
		QVector<float**> deformationGradientMatrixs(numOfFacets);

		for(int i = 0; i < numOfFacets; i++)
		{
			float** m = new float*[order];
			for(int j = 0; j < order; j++)
				m[j] = new float[order];
			deformationGradientMatrixs[i] = m;
		}

		for(Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++)
		{
			Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();

			Point_3 p1 = c->vertex()->point();
			++c;
			Point_3 p2 = c->vertex()->point();
			++c;
			Point_3 p3 = c->vertex()->point();

			Vector_3 v1 = p2-p1;
			Vector_3 v2 = p3-p1;
			Vector_3 v3 = cross_product(v1, v2);
			v3 = v3 / sqrt(v3.squared_length());

			float** orientationMatrix = new float*[order];
			for(int i = 0; i < order; i++)
				orientationMatrix[i] = new float[order];

			orientationMatrix[0][0] = v1.x();
			orientationMatrix[1][0] = v1.y();
			orientationMatrix[2][0] = v1.z();

			orientationMatrix[0][1] = v2.x();
			orientationMatrix[1][1] = v2.y();
			orientationMatrix[2][1] = v2.z();

			orientationMatrix[0][2] = v3.x();
			orientationMatrix[1][2] = v3.y();
			orientationMatrix[2][2] = v3.z();

			int fit_index = fit->index();
			float** referenceOrientatoinMatrix = referenceOrientationMatrixs.at(fit_index);
			float** deformationGradient = matrixProduct(orientationMatrix, referenceOrientatoinMatrix, order);

			deformationGradientMatrixs[fit_index] = deformationGradient;
		}

		// find the largest deformation distance between two adjacent faces
		for(Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++)
		{
			float maxDeforDist = FLT_MIN;

			int f_index1 = fit->index();

			Mesh::Halfedge_around_facet_circulator pHalfedge = fit->facet_begin();
			do { 
				if(pHalfedge->opposite()->facet() == NULL)
					continue;
						
				int f_index2 = pHalfedge->opposite()->facet()->index();

				float** m1 = deformationGradientMatrixs.at(f_index1);
				float** m2 = deformationGradientMatrixs.at(f_index2);

				float dd = deformationDistance(m1, m2, order);

				if( dd > maxDeforDist )
					maxDeforDist = dd;
			} while(++pHalfedge != fit->facet_begin());
			
			if(maxDeforDist > maxDDegree[f_index1])
				maxDDegree[f_index1] = maxDeforDist;
		}

		// release memory
		for(int i = 0; i < numOfFacets; i++)
		{
			float** defor = deformationGradientMatrixs.at(i);
			for(int i = 0; i < order; i++)
				delete [] defor[i];
			delete [] defor;
		}

		mesh->clear();

		cout << "ok!" << endl;
	}
	
	// release memory //
	delete mesh;

	// save the dd to reference frame mesh //
	float maxDD = FLT_MIN;
	float minDD = FLT_MAX;
	for(Mesh::Facet_iterator fit = referenceMesh->facets_begin(); fit != referenceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		float mdd = maxDDegree[f_index];

		fit->deformationDegree(mdd);

		if(mdd > maxDD) maxDD = mdd;
		if(mdd < minDD) minDD = mdd;
	}

	// normalize
	for(Mesh::Facet_iterator fit = referenceMesh->facets_begin(); fit != referenceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		float mdd = fit->deformationDegree();
		float n_mdd = (mdd - minDD) / (maxDD - minDD);
		
		fit->normalizeDeformationDegree(n_mdd);
	}

	// release memory //
	for(int i = 0; i < numOfFacets; i++)
	{
		float** ref = referenceOrientationMatrixs[i];
		for(int i = 0; i < order; i++)
			delete [] ref[i];
		delete [] ref;
	}
}

void DirRunner::dirRunnerNoiseAttack(const QString& fileDir)
{
	Mesh *mesh = new Mesh;

	float noiseIntensity[2] = {0.0005, 0.001};
	bool preserveBoundaries = true;

	for(int noiseIndex = 0; noiseIndex < 2; noiseIndex++) {
		// check whether target directory is exist or not. if not, mkdir it
		float noise = noiseIntensity[noiseIndex];

		QString stringIntensity;
		stringIntensity.setNum(noise, 'g', 6);

		QString targetDir = fileDir + "\\" + "noise_intensity" + stringIntensity;

		QDir td(targetDir);
		if(!td.exists())
			td.mkpath(targetDir);

		// list dir files
		QDir fd(fileDir);
		fd.setFilter(QDir::Files);

		const QFileInfoList list = fd.entryInfoList();

		QFileInfoList::const_iterator iterator = list.begin();

		while ( iterator != list.end() ) {
			QString inFileName = (*iterator).fileName();
			QString filePath = fileDir + "\\" + inFileName;

			QFileInfo fi(inFileName);

			QString fileName = fi.completeBaseName();
			QString fileSuffix = "." + fi.completeSuffix();
		
			if(fi.completeSuffix() != "obj") continue;

			cout << qPrintable(filePath) << ":  ";
			++iterator;

			QString boundary = (preserveBoundaries==true) ? "keepboundary" : "notkeepboundary";
			QString stringOutFilename = targetDir + "\\" + fileName + "_noise_intensity" + stringIntensity + "_" + boundary + fileSuffix;
			
			QFileInfo outFile(stringOutFilename);
			if(outFile.exists()) { cout << "objExists!" << endl; continue; }

			bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);
			
			if(!success) continue;

			DirRunnerAttack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(mesh, stringOutFilename);

			attack.NoiseAdditionUniform(noise, preserveBoundaries);

			// release memory //
			mesh->clear();

			cout << "ok!" << endl;
			
		}
	}

	delete mesh; // release memory
}

void DirRunner::dirRunnerQuantizationAttack(const QString& fileDir)
{
	Mesh *mesh = new Mesh;

	int quantizationBit[2] = {10, 11};

	for(int quantizationIndex = 0; quantizationIndex < 2; quantizationIndex++) {
		// check whether target directory is exist or not. if not, mkdir it
		int bitDepth = quantizationBit[quantizationIndex];

		QString stringBitDepth;
		stringBitDepth.setNum(bitDepth, 10);

		QString targetDir = fileDir + "\\" + "quantization_bitdepth" + stringBitDepth;

		QDir td(targetDir);
		if(!td.exists())
			td.mkpath(targetDir);

		// list dir files
		QDir fd(fileDir);
		fd.setFilter(QDir::Files);

		const QFileInfoList list = fd.entryInfoList();
		QFileInfoList::const_iterator iterator = list.begin();

		while ( iterator != list.end() ) {
			QString inFileName = (*iterator).fileName();
			QString filePath = fileDir + "\\" + inFileName;

			QFileInfo fi(inFileName);
			QString fileName = fi.completeBaseName();
			QString fileSuffix = "." + fi.completeSuffix();	
			
			if(fi.completeSuffix() != "obj") continue;

			cout << qPrintable(filePath) << ":  ";
			++iterator;

			QString stringOutFilename = targetDir + "\\" + fileName + "_quantization_bitdepth" + stringBitDepth + fileSuffix;
			
			QFileInfo outFile(stringOutFilename);
			if(outFile.exists()) { cout << "objExists!" << endl; continue; }

			bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);
			
			if(!success) continue;

			DirRunnerAttack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(mesh, stringOutFilename);

			attack.CoordinateQuantization(bitDepth);

			// release memory //
			mesh->clear();

			cout << "ok!" << endl;
		}
	}

	delete mesh; // release memory
}

void DirRunner::dirRunnerReorderingAttack(const QString& fileDir)
{
	Mesh *mesh = new Mesh;

	// check whether target directory is exist or not. if not, mkdir it
	QString targetDir = fileDir + "\\" + "elementreordering";

	QDir td(targetDir);
	if(!td.exists())
		td.mkpath(targetDir);

	// list dir files
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	QFileInfoList::const_iterator iterator = list.begin();

	while ( iterator != list.end() ) {
		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		QFileInfo fi(inFileName);
		QString fileName = fi.completeBaseName();
		QString fileSuffix = "." + fi.completeSuffix();

		if(fi.completeSuffix() != "obj") continue;
		
		cout << qPrintable(filePath) << ":  ";
		++iterator;

		QString stringOutFilename = targetDir + "\\" + fileName + "_elementreordering" + fileSuffix;	
		
		QFileInfo outFile(stringOutFilename);
		if(outFile.exists()) { cout << "objExists!" << endl; continue;}

		bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);
		
		if(!success) continue;

		DirRunnerAttack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(mesh, stringOutFilename);

		attack.ElementReordering();

		// release memory //
		mesh->clear();

		cout << "ok!" << endl;
	}

	delete mesh; // release memory
}

void DirRunner::dirRunnerSimilarityTransformAttack(const QString& fileDir)
{
	Mesh *mesh = new Mesh;

	// check whether target directory is exist or not. if not, mkdir it
	QString targetDir = fileDir + "\\" + "similaritytransform";

	QDir td(targetDir);
	if(!td.exists())
		td.mkpath(targetDir);

	// list dir files
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	QFileInfoList::const_iterator iterator = list.begin();

	while ( iterator != list.end() ) {
		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		QFileInfo fi(inFileName);
		QString fileName = fi.completeBaseName();
		QString fileSuffix = "." + fi.completeSuffix();	
		
		if(fi.completeSuffix() != "obj") continue;

		cout << qPrintable(filePath) << ":  ";
		++iterator;

		QString stringOutFilename = targetDir + "\\" + fileName + "_similaritytransform" + fileSuffix;
		
		QFileInfo outFile(stringOutFilename);
		if(outFile.exists()) { cout << "objExists!" << endl; continue; }

		bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);
		
		if(!success) continue;

		DirRunnerAttack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(mesh, stringOutFilename);

		attack.SimilarityTransformation();

		// release memory //
		mesh->clear();

		cout << "ok!" << endl;
	}

	delete mesh; // release memory
}

void DirRunner::dirRunnerSimplificationAttack(const QString& fileDir)
{
	Mesh *mesh = new Mesh;

	float simplificationRatio[2] = {10, 30};
	SimplificationType simpType = LINDSTROMTURK;
	//SimplificationType simplificationType = EDGELENGTHMIDPOINT;

	for(int simplificationIndex = 0; simplificationIndex < 2; simplificationIndex++) {
		// check whether target directory is exist or not. if not, mkdir it
		float ratio = simplificationRatio[simplificationIndex];

		QString stringType;
		if(simpType==EDGELENGTHMIDPOINT)
			stringType = "EDGELENGTHMIDPOINT";
		else if(simpType==LINDSTROMTURK)
			stringType = "LINDSTROMTURK";
		stringType = stringType.toLower();

		QString stringRatio;
		stringRatio.setNum(ratio, 'g', 6);

		QString targetDir = fileDir + "\\" + "simplification_policy" + stringType +"_ratio" + stringRatio;

		QDir td(targetDir);
		if(!td.exists())
			td.mkpath(targetDir);

		// list dir files
		QDir fd(fileDir);
		fd.setFilter(QDir::Files);

		const QFileInfoList list = fd.entryInfoList();
		QFileInfoList::const_iterator iterator = list.begin();

		while ( iterator != list.end() ) {
			QString inFileName = (*iterator).fileName();
			QString filePath = fileDir + "\\" + inFileName;

			QFileInfo fi(inFileName);
			QString fileName = fi.completeBaseName();
			QString fileSuffix = "." + fi.completeSuffix();	
			
			if(fi.completeSuffix() != "obj") continue;

			cout << qPrintable(filePath) << ":  ";
			++iterator;

			QString stringOutFilename = targetDir + "\\" + fileName + "_simplification_policy" + stringType +"_ratio" + stringRatio + fileSuffix;
			
			QFileInfo outFile(stringOutFilename);
			if(outFile.exists()) { cout << "objExists!" << endl; continue; }

			bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);
			
			if(!success) continue;

			DirRunnerAttack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(mesh, stringOutFilename);

			attack.Simplification(simpType, ratio);

			// release memory //
			mesh->clear();

			cout << "ok!" << endl;
		}
	}

	delete mesh; // release memory
}

void DirRunner::dirRunnerSmoothingAttack(const QString& fileDir)
{
	Mesh* mesh = new Mesh;

	float deformFactor = 0.1;
	int iterationList[2] = {5, 10};
	bool preserveBoundaries = true;
	
	for(int iterationIndex = 0; iterationIndex < 2; iterationIndex++) {
		// check whether target directory is exist or not. if not, mkdir it
		int iteraNum = iterationList[iterationIndex];

		QString stringLambda, stringIteration;
		stringLambda.setNum(deformFactor, 'g', 6 );
		stringIteration.setNum(iteraNum, 10);

		QString targetDir = fileDir + "\\" + "smoothing_lambda" + stringLambda + "_iteration" + stringIteration;

		QDir td(targetDir);
		if(!td.exists())
			td.mkpath(targetDir);

		// list dir files
		QDir fd(fileDir);
		fd.setFilter(QDir::Files);

		const QFileInfoList list = fd.entryInfoList();
		QFileInfoList::const_iterator iterator = list.begin();

		while ( iterator != list.end() ) {
			QString inFileName = (*iterator).fileName();
			QString filePath = fileDir + "\\" + inFileName;

			QFileInfo fi(inFileName);
			QString fileName = fi.completeBaseName();
			QString fileSuffix = "." + fi.completeSuffix();	
			
			if(fi.completeSuffix() != "obj") continue;

			cout << qPrintable(filePath) << ":  ";
			++iterator;

			QString boundary = (preserveBoundaries==true) ? "keepboundary" : "notkeepboundary";
			QString stringOutFilename = targetDir + "\\" + fileName + "_smoothing_lambda" + stringLambda + "_iteration" + stringIteration + "_" + boundary + fileSuffix;

			QFileInfo outFile(stringOutFilename);
			if(outFile.exists()) { cout << "objExists!" << endl; continue; }

			bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);
			
			if(!success) continue;

			DirRunnerAttack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(mesh, stringOutFilename);

			attack.LaplacianSmoothing(deformFactor, iteraNum, preserveBoundaries);

			// release memory //
			mesh->clear();

			cout << "ok!" << endl;
		}
	}

	delete mesh; // release memory
}

void DirRunner::dirRunnerSubdivisionAttack(const QString& fileDir)
{
	Mesh *mesh = new Mesh;

	int subdivisionIteration = 1;
	SubdivisionType subdivisionType = MIDPOINT;

	// check whether target directory is exist or not. if not, mkdir it
	QString stringType;
	if(subdivisionType==CATMULLCLARK)
		stringType = "CATMULLCLARK";
	else if(subdivisionType==LOOP)
		stringType = "LOOP";
	else if(subdivisionType==DOOSABIN)
		stringType = "DOOSABIN";
	else if(subdivisionType==SQRT3)
		stringType = "SQRT3";
	else if(subdivisionType==MIDPOINT)
		stringType = "MIDPOINT";
	stringType = 	stringType.toLower();

	QString targetDir = fileDir + "\\" + "subdivision_" + stringType;

	QDir td(targetDir);
	if(!td.exists())
		td.mkpath(targetDir);

	// list dir files
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	QFileInfoList::const_iterator iterator = list.begin();

	while ( iterator != list.end() ) {
		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		QFileInfo fi(inFileName);
		QString fileName = fi.completeBaseName();
		QString fileSuffix = "." + fi.completeSuffix();	
		
		if(fi.completeSuffix() != "obj") continue;

		cout << qPrintable(filePath) << ":  ";
		++iterator;

		QString stringDepth;
		stringDepth.setNum(subdivisionIteration, 10);

		QString stringOutFilename = targetDir + "\\" + fileName + "_subdivision_" + stringType + "_iteration"+ stringDepth + fileSuffix;

		QFileInfo outFile(stringOutFilename);
		if(outFile.exists()) { cout << "objExists!" << endl; continue; }

		bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);
		
		if(!success) continue;

		DirRunnerAttack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(mesh, stringOutFilename);

		attack.Subdivision(subdivisionType, subdivisionIteration);

		// release memory //
		mesh->clear();

		cout << "ok!" << endl;
	}

	delete mesh; // release memory
}
*/