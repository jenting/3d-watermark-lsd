#include "Importer.h"


bool FileParser::read_OffPly2(const char* filename, Mesh* mesh, int err)
{
	//SimpleMeshBuilder b(mesh);

	cout << "read off: " << filename << endl;

	RepeatingMeshBuilder b(mesh);
	Ply2Reader r(&b);
	return r.read(filename);
}

bool FileParser::read_Obj(const char* filename, Mesh* mesh, int err)
{
	//SimpleMeshBuilder b(mesh);

	cout << "read obj: " << filename << endl;

	RepeatingMeshBuilder b(mesh);
	ObjReader r(&b);
	return r.read(filename);
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


bool ImporterMesh::loadMesh(const QString& fileName, Mesh* mesh, double err)
{
	QFileInfo fi(fileName);

	QString ext = fi.extension(FALSE);

	bool success = false;

	//cout << "fileName : " << qPrintable(fileName) << endl;
	//cout << "ext : " << qPrintable(ext) << endl;

	try
	{
		if (ext == "simple" || ext == "ply2" || ext == "off")
		{
			success = FileParser::read_OffPly2(fileName.toAscii(), mesh, err);
		}
		else if (ext == "obj")
		{
			success = FileParser::read_Obj(fileName.toAscii(), mesh, err);
		}

		if (success)
		{
			// deal with quads properly
			for (Mesh::Facet_iterator it = mesh->facets_begin(); it != mesh->facets_end(); it++)
				mesh->triangulateFacet(it);

			int i=0;
			for (Mesh::Vertex_iterator vit = mesh->vertices_begin(); vit != mesh->vertices_end(); vit++)
			{
				vit->index(i);
				i++;
			}
			i=0;
			for (Mesh::Facet_iterator it = mesh->facets_begin(); it != mesh->facets_end(); it++)
			{
				it->index(i);
				i++;
			}
			success = (mesh->size_of_facets()!=0);

			try
			{
				i = 0;
				for (Mesh::Facet_const_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++)
				{
					Mesh::Halfedge_const_handle c = fit->halfedge();
					Mesh::Vertex_const_handle v1 = c->vertex();
					++i;
				}
			}
			catch(...)
			{
				printf("cought exception in pass.\n");
				return false;
			}
			
		}
		else
		{
			printf("FAILED!");
		}

	}
	catch (...)
	{
		printf("Cought exception in loadMesh.\n");
		return false;
	}

	if (success)
	{
		//AppManager::doMeshPrecalc(mesh);

		mesh->computeNormals();
		mesh->computeBoundingBox();
		mesh->computeTriangleSurfaces();
		mesh->centerOfMass();
		mesh->compute_type();
	}

	return success;
}


bool ImporterMesh::loadMesh(const QString& fileName, AppObject* appObject, double err)
{
	QFileInfo fi(fileName);

	QString ext = fi.extension(FALSE);

	bool success = false;

	try
	{
		/*	if (ext == "off")
		{
			std::ifstream stream(fileName.toAscii());
			if(!stream)
			{
				printf("failed to open file in loadMesh.\n");
				return false;
			}
			//stream >> *mesh;
			CGAL::scan_OFF(stream, *appObject->mesh, true);

			stream.close();
		}
		else*/ if (ext == "simple" || ext == "ply2" || ext == "off")
		{
			success = FileParser::read_OffPly2(fileName.toAscii(), appObject->mesh, err);
		}
		else if (ext == "obj")
		{
			//success = loadObjFile(fileName, appObject->mesh);
			//Wavefront_Parser_obj parser;
			//success = parser.read(fileName.toAscii(), appObject->mesh);

			success = FileParser::read_Obj(fileName.toAscii(), appObject->mesh, err);
		}

		if (success)
		{
			// deal with quads properly
			for (Mesh::Facet_iterator it = appObject->mesh->facets_begin(); it != appObject->mesh->facets_end(); it++)
				appObject->mesh->triangulateFacet(it);
			int i=0;
			for (Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin(); vit != appObject->mesh->vertices_end(); vit++) {
				vit->index(i);
				i++;
			}
			i=0;
			for (Mesh::Facet_iterator it = appObject->mesh->facets_begin(); it != appObject->mesh->facets_end(); it++)
			{
				it->index(i);
				i++;
			}
			success = (appObject->mesh->size_of_facets()!=0);

			try
			{
				i = 0;
				for (Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin(); fit != appObject->mesh->facets_end(); fit++)
				{
					Mesh::Halfedge_const_handle c = fit->halfedge();
					Mesh::Vertex_const_handle v1 = c->vertex();
					++i;
				}
			}
			catch(...)
			{
				printf("cought exception in pass.\n");
				return false;
			}
			
		}
		else
		{
			printf("FAILED!");
		}

	}
	catch (...)
	{
		printf("Cought exception in loadMesh.\n");
		return false;
	}


	if (success)
		AppManager::doMeshPrecalc(appObject->mesh);

	return success;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


bool ImporterMeshFeatures::loadVertexNormals(const QString& featuresFileName, double** vertexNormals)
{
	QFile file(featuresFileName);
	if (!file.open(QIODevice::ReadOnly))
		return false;
	QTextStream in(&file);

	int pointIndex = 0;
	QString line;
	do {
		line = in.readLine();
		if (line.length() == 0)
			continue;

		QStringList flds = line.split(' ', QString::SkipEmptyParts);

		vertexNormals[pointIndex][0] = flds[0].toDouble();
		vertexNormals[pointIndex][1] = flds[1].toDouble();
		vertexNormals[pointIndex][2] = flds[2].toDouble();

		pointIndex++;

	} while (!line.isNull());

	return true;
}


bool ImporterMeshFeatures::loadSDFVertices(const QString& featuresFileName, const bool normalize, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_vertices();
	double* sdf = new double[sdf_size];

	QString description;
	bool success = FileUtils::load(featuresFileName, FEATURES_FILE_SDF_VERTICES, description, (char*)sdf, sdf_size);

	if (!success)
	{
		delete[] sdf;
		return false;
	}

	Mesh::Vertex_iterator it = appObject->mesh->vertices_begin();
	Mesh::Vertex_iterator it_end = appObject->mesh->vertices_end();
	int i=0;
	for (;it != it_end; it++)
	{
		if (sdf[i]==sdf[i])
			it->volumeSDF(sdf[i]);
		else
			it->volumeSDF(0.0);
		i++;
	}

	delete[] sdf;

	return true;
}


bool ImporterMeshFeatures::loadSDFFacets(const QString& featuresFileName, const bool normalize, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_facets();
	double* sdf = new double[sdf_size];

	QString description;
	bool success = FileUtils::load(featuresFileName, FEATURES_FILE_SDF_FACETS, description, (char*)sdf, sdf_size);

	if (!success)
	{
		delete[] sdf;
		return false;
	}

	int i=0;
	for (Mesh::Facet_iterator it = appObject->mesh->facets_begin(); it != appObject->mesh->facets_end(); ++it)
	{
		if (sdf[i]==sdf[i])
			it->volumeSDF(sdf[i]);
		else
			it->volumeSDF(0.0);
		i++;
	}

//	appObject->mesh->fillNormalizedFacetVolume();

	/*appObject->mesh->makeFacesNVolume(m_appParameters.lsd_smoothing,
	m_appParameters.lsd_smoothing_anisotropic, m_appParameters.lsd_smoothing_iterations);*/

	delete[] sdf;
	return true;
}


bool ImporterMeshFeatures::diffSDFFacets(const QString& featuresFileName, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_facets();
	double* sdf = new double[sdf_size];

	QString description;
	bool success = FileUtils::load(featuresFileName, FEATURES_FILE_SDF_FACETS, description, (char*)sdf, sdf_size);

	if (!success)
	{
		delete[] sdf;
		return false;
	}

	Mesh::Facet_iterator it = appObject->mesh->facets_begin();
	Mesh::Facet_iterator it_end = appObject->mesh->facets_end();
	int i=0;
	for (;it != it_end; it++) {
		it->volumeSDF() = fabs(sdf[i] - it->volumeSDF());
		i++;
	}

	appObject->mesh->fillNormalizedFacetVolume();

	delete[] sdf;
	return true;
}


bool ImporterMeshFeatures::loadSDFText(const QString& importFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(importFileName);
	if (!file.open(QIODevice::ReadOnly))
		return false;
	QTextStream in(&file);

	unsigned int sdf_size = appObject->mesh->size_of_vertices();
	TNT::Array1D<double> sdf = TNT::Array1D<double>(sdf_size);

	int i = 0;
	QString line;	
	do {
		line = in.readLine();
		sdf[i] = line.toDouble();
		i++;
	} while (!line.isNull() && i < sdf_size);

	i=0;
	double minSDF = FLT_MAX;
	double maxSDF = -FLT_MAX;
	double avgSDF = 0.0;
	for (Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin(); vit != appObject->mesh->vertices_end(); vit++)
	{
		vit->volumeSDF(sdf[i]);

		if (sdf[i] < minSDF) minSDF = sdf[i];
		if (sdf[i] > maxSDF) maxSDF = sdf[i];
		avgSDF += sdf[i];

		i++;
	}
	avgSDF = avgSDF / sdf_size;

	appObject->mesh->setMinVolume(minSDF);
	appObject->mesh->setMaxVolume(maxSDF);
	appObject->mesh->setAvgVolume(avgSDF);

	double stdDevSDF = 0.0;
	for (Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin(); vit != appObject->mesh->vertices_end(); vit++)
	{
		stdDevSDF += pow(vit->volumeSDF() - avgSDF, 2);
	}
	stdDevSDF = sqrt(stdDevSDF / sdf_size);

	double minStdDevSDF = FLT_MAX;
	double maxStdDevSDF = -FLT_MAX;
	double avgStdDevSDF = 0.0;
	int stdDevCounter = 0;
	for (Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin(); vit != appObject->mesh->vertices_end(); vit++)
	{
		double sdf = vit->volumeSDF();

		if (fabs(sdf - avgSDF) <= 2 * stdDevSDF)
		{
			if (sdf < minStdDevSDF) minStdDevSDF = sdf;
			if (sdf > maxStdDevSDF) maxStdDevSDF = sdf;
			avgStdDevSDF += sdf;

			stdDevCounter++;
		}
	}
	avgStdDevSDF = avgStdDevSDF / stdDevCounter;

	appObject->mesh->setMinStdDevVolume(minStdDevSDF);
	appObject->mesh->setMaxStdDevVolume(maxStdDevSDF);
	appObject->mesh->setAvgStdDevVolume(avgStdDevSDF);

	cout << "load sdf: " << qPrintable(importFileName) << " ...done!" << endl;

	/*
	cout << "minSDF = " << minSDF << " , maxSDF = " << maxSDF << " , avgSDF = " << avgSDF << endl;
	cout << "minStdDevSDF = " << minStdDevSDF << " , maxStdDevSDF = " << maxStdDevSDF << " , avgStdDevSDF = " << avgStdDevSDF << endl;
	*/

	return true;
}

bool ImporterMeshFeatures::loadNSDFText(const QString& importFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(importFileName);
	if (!file.open(QIODevice::ReadOnly))
		return false;
	QTextStream in(&file);

	unsigned int sdf_size = appObject->mesh->size_of_vertices();
	TNT::Array1D<double> nsdf = TNT::Array1D<double>(sdf_size);

	int i = 0;
	QString line;	
	do {
		line = in.readLine();
		nsdf[i] = line.toDouble();
		i++;
	} while (!line.isNull() && i < sdf_size);

	i=0;
	double minNSDF = FLT_MAX;
	double maxNSDF = -FLT_MAX;
	double avgNSDF = 0.0;
	for (Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin(); vit != appObject->mesh->vertices_end(); vit++)
	{
		vit->volumeNSDF(nsdf[i]);

		if (nsdf[i] < minNSDF) minNSDF = nsdf[i];
		if (nsdf[i] > maxNSDF) maxNSDF = nsdf[i];
		avgNSDF += nsdf[i];

		i++;
	}

	cout << "sum NSDF = " << avgNSDF << endl;
	avgNSDF = avgNSDF / sdf_size;
	cout << "avg NSDF = " << avgNSDF << endl;

	appObject->mesh->setMinNVolume(minNSDF);
	appObject->mesh->setMaxNVolume(maxNSDF);
	appObject->mesh->setAvgNVolume(avgNSDF);
	
	double stdDevNSDF = 0.0;
	for (Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin(); vit != appObject->mesh->vertices_end(); vit++)
	{
		stdDevNSDF += pow(vit->volumeNSDF() - avgNSDF, 2);
	}
	stdDevNSDF = sqrt(stdDevNSDF / sdf_size);

	double minStdDevNSDF = FLT_MAX;
	double maxStdDevNSDF = -FLT_MAX;
	double avgStdDevNSDF = 0.0;
	int stdDevCounter = 0;
	for (Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin(); vit != appObject->mesh->vertices_end(); vit++)
	{
		double nsdf = vit->volumeNSDF();

		if (fabs(nsdf - avgNSDF) <= 2 * stdDevNSDF)
		{
			if (nsdf < minStdDevNSDF) minStdDevNSDF = nsdf;
			if (nsdf > maxStdDevNSDF) maxStdDevNSDF = nsdf;
			avgStdDevNSDF += nsdf;

			stdDevCounter++;
		}
	}
	avgStdDevNSDF = avgStdDevNSDF / stdDevCounter;

	appObject->mesh->setMinStdDevNVolume(minStdDevNSDF);
	appObject->mesh->setMaxStdDevNVolume(maxStdDevNSDF);
	appObject->mesh->setAvgStdDevNVolume(avgStdDevNSDF);

	cout << "load nsdf: " << qPrintable(importFileName) << " ...done!" << endl;

	/*
	cout.precision(10);
	cout << "minNSDF = " << minNSDF << " , maxNSDF = " << maxNSDF << " , avgNSDF = " << avgNSDF << endl;
	cout << "minStdDevNSDF = " << minStdDevNSDF << " , maxStdDevNSDF = " << maxStdDevNSDF << " , avgStdDevNSDF = " << avgStdDevNSDF << endl;
	*/

	return true;
}