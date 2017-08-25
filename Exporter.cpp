#include "Exporter.h"

#include "Watermark.h"


bool FileUtils::load(const QString& fileName, const int fileType, QString& description, char* data, unsigned int dataSize)
{
	QFile file(fileName);
	file.open(QIODevice::ReadOnly);
	QDataStream stream(&file);

	Q_UINT32 fileTypeInFile;
	stream >> fileTypeInFile >> description;

	if (fileType != fileTypeInFile) return false;

	Q_UINT32 dataSizeInFile;
	stream >> dataSizeInFile;

	dataSize = dataSizeInFile;

	stream.readRawBytes(data, dataSize);

	file.close();

	return true;
}

bool FileUtils::save(const QString& fileName, const int fileType, const QString& description, const char* data, const unsigned int dataSize)
{
	QFile file(fileName);
	file.open(QIODevice::WriteOnly);
	QDataStream stream(&file);

	stream << (Q_UINT32)fileType;
	stream << description;

	stream << (Q_UINT32)dataSize;

	stream.writeRawBytes(data, dataSize);

	file.close();

	return true;
}


bool ExporterMeshFile::write_obj(const QString& filename, AppObject* appObject)
{
	return appObject->mesh->write_obj(filename, 1);
}

bool ExporterMeshFile::write_off(const QString& filename, AppObject* appObject)
{
	return appObject->mesh->write_off(filename);
}

/*
bool ExporterMeshFile::saveObjWithFaceSDF(const QString& filename, AppObject* appObject)
{
	return saveObjWithFaceColor(filename, appObject, VolumeFaceColorizer(m_worldManager), false);
}

bool ExporterMeshFile::saveObjWithFaceColor(const QString& filename, AppObject* appObject, FaceColorizer& cizr, bool groups)
{
	QFile data(filename);
	if (!data.open(QFile::WriteOnly | QFile::Truncate))
		return false;
	data.setTextModeEnabled(true);

	QTextStream out(&data);
	out << "# SDF for facets.\n\n";
	out << "mtllib sdfDefault.mtl\n";
	Mesh *mesh = appObject->mesh;
	for (Mesh::Vertex_iterator ivertex = mesh->vertices_begin(); ivertex != mesh->vertices_end(); ivertex++)
	{
		Mesh::Vertex_handle v = ivertex;
		Mesh::Point p = v->point();
		out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
	}
	out << "\n";

	typedef QList<Mesh::Facet_handle> TFacetsList;
	typedef QHash<QColor, TFacetsList> TFacetsHash;
	TFacetsHash bins;

	for (Mesh::Facet_iterator fcit = mesh->facets_begin(); fcit != mesh->facets_end(); ++fcit)
	{
		Mesh::Facet_handle f = fcit;
		QColor c = cizr.getColor(f);
		bins[c].push_back(f);
	}

	for (TFacetsHash::const_iterator it = bins.constBegin(); it != bins.constEnd(); ++it)
	{
		QColor c = it.key();
		if (groups)
		{
			out << "g " << "group_" << c.red() << "_" << c.green() << "_" << c.blue() << "\n";
		}
		out << "\nusemtl " << "sdfcol_" << c.red() << "_" << c.green() << "_" << c.blue() << "\n";
		const TFacetsList &flst = it.value();
		for (TFacetsList::const_iterator fit = flst.constBegin(); fit != flst.constEnd(); ++fit)
		{
			Mesh::Facet_handle f = *fit;
			out << "f ";
			Mesh::Halfedge_around_facet_circulator	fvait = f->halfedge()->facet_begin();
			Mesh::Halfedge_around_facet_circulator	begin = fvait;
			CGAL_For_all(fvait, begin)
			{
				//Mesh::Vertex_handle v = fvait;
				out << fvait->vertex()->index() + 1<< " ";
			}
			out << "\n";
		}

	}

	return true;
}
*/

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


bool ExporterMeshFeature::saveSDFText(const QString& exportFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(exportFileName);
	if (!file.open(QIODevice::WriteOnly))
		return false;

	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = appObject->mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = appObject->mesh->vertices_end();
		for (; vit != vit_end; vit++)
			ts << vit->volumeSDF() << '\n';
	}
	else
	{
		Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = appObject->mesh->facets_end();
		for (; fit != fit_end; fit++)
			ts << fit->volumeSDF() << '\n';
	}

	file.close();

	cout << "export sdf text => " << qPrintable(exportFileName) << " done!" << endl;

	return true;
}

bool ExporterMeshFeature::saveSDFDiffText(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_facets();

	QFile file(featuresFileName);
	if (!file.open(QIODevice::WriteOnly))
		return false;

	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = appObject->mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = appObject->mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			double feature = vit->volumeSDF();
			double maxNeighborDiff = 0.0;
			Mesh::Halfedge_around_vertex_const_circulator c = vit->vertex_begin();
			do {
				Mesh::Vertex_const_handle neighbor = c->opposite()->vertex();
				if (neighbor == NULL) continue;

				double nvalue = neighbor->volumeSDF();

				if (fabs(nvalue-feature)>maxNeighborDiff) maxNeighborDiff = fabs(nvalue-feature);
			} while (++c != vit->vertex_begin());

			//double diff = MIN(maxNeighborDiff / (feature + 1e-4),1);

			ts << maxNeighborDiff << '\n';
		}
	}
	else
	{
		Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = appObject->mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			double feature = fit->volumeSDF();
			double maxNeighborDiff = 0.0;
			Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();
			do {
				Mesh::Facet_const_handle neighbor = c->opposite()->facet();
				if (neighbor == NULL) continue;

				double nvalue = neighbor->volumeSDF();

				if (fabs(nvalue-feature)>maxNeighborDiff) maxNeighborDiff = fabs(nvalue-feature);
			} while (++c != fit->facet_begin());

			//double diff = MIN(maxNeighborDiff / (feature + 1e-4),1);

			ts << maxNeighborDiff << '\n';
		}
	}

	file.close();

	return true;
}

bool ExporterMeshFeature::saveNSDFText(const QString& exportFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(exportFileName);
	if (!file.open(QIODevice::WriteOnly))
		return false;

	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = appObject->mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = appObject->mesh->vertices_end();
		for (; vit != vit_end; vit++)
			ts << vit->volumeNSDF() << '\n';
	}
	else
	{
		Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = appObject->mesh->facets_end();
		for (; fit != fit_end; fit++)
			ts << fit->volumeNSDF() << '\n';
	}

	file.close();

	cout << "export nsdf text => " << qPrintable(exportFileName) << " done!" << endl;

	return true;
}

bool ExporterMeshFeature::saveNSDFDiffText(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	//unsigned int sdf_size = appObject->mesh->size_of_facets();

	QFile file(featuresFileName);
	if (!file.open(QIODevice::WriteOnly))
		return false;
	
	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = appObject->mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = appObject->mesh->vertices_end();

		for (; vit != vit_end; vit++)
		{
			double feature = vit->volumeNSDF();
			double maxNeighborDiff = 0.0;
			Mesh::Halfedge_around_vertex_const_circulator c = vit->vertex_begin();
			do {
				Mesh::Vertex_const_handle neighbor = c->opposite()->vertex();
				if (neighbor == NULL) continue;

				double nvalue = neighbor->volumeNSDF();

				if (fabs(nvalue-feature)>maxNeighborDiff) maxNeighborDiff = fabs(nvalue-feature);

			} while (++c != vit->vertex_begin());

			//double diff = MIN(maxNeighborDiff / (feature + 1e-4),1);

			ts << maxNeighborDiff << '\n';
		}
	}
	else
	{
		Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = appObject->mesh->facets_end();

		for (; fit != fit_end; fit++)
		{
			double feature = fit->volumeNSDF();
			double maxNeighborDiff = 0.0;
			Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();
			do {
				Mesh::Facet_const_handle neighbor = c->opposite()->facet();
				if (neighbor == NULL) continue;

				double nvalue = neighbor->volumeNSDF();

				if (fabs(nvalue-feature)>maxNeighborDiff) maxNeighborDiff = fabs(nvalue-feature);

			} while (++c != fit->facet_begin());

			//double diff = MIN(maxNeighborDiff / (feature + 1e-4),1);

			ts << maxNeighborDiff << '\n';
		}
	}

	file.close();

	return true;
}


bool ExporterMeshFeature::saveVSIRefPoints(const QString& featuresFileName, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		for (Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin(); fit != appObject->mesh->facets_end(); fit++)
			ts << fit->referencePointX() << ' ' << fit->referencePointY() << ' ' << fit->referencePointZ() << '\n';

		file.close();
		return true;
	}
	
	return false;
}


bool ExporterMeshFeature::saveVSIText(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		if (onVertices)
		{
			Mesh::Vertex_const_iterator vit = appObject->mesh->vertices_begin();
			Mesh::Vertex_const_iterator vit_end = appObject->mesh->vertices_end();
			for (; vit != vit_end; vit++)
				ts << vit->volumeVSI() << '\n';
		}
		else
		{
			Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin();
			Mesh::Facet_const_iterator fit_end = appObject->mesh->facets_end();
			for(; fit != fit_end; fit++)
				ts << fit->volumeVSI() << '\n';
		}

		file.close();
		return true;
	}

	return false;
}

bool ExporterMeshFeature::saveMaximalDiameter(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		if (onVertices)
		{
			Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin();
			Mesh::Vertex_iterator vit_end = appObject->mesh->vertices_end();
			for (; vit != vit_end; vit++)
				ts << vit->maxDist() << '\n';
		}
		else
		{
			Mesh::Facet_iterator fit = appObject->mesh->facets_begin();
			Mesh::Facet_iterator fit_end = appObject->mesh->facets_end();
			for (; fit != fit_end; fit++)
				ts << fit->maxDist() << '\n';
		}

		file.close();
		return true;
	}

	return false;
}

bool ExporterMeshFeature::saveMinimalDiameter(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		if (onVertices)
		{
			Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin();
			Mesh::Vertex_iterator vit_end = appObject->mesh->vertices_end();
			for (; vit != vit_end; vit++)
				ts << vit->minDist() << '\n';
		}
		else
		{
			Mesh::Facet_iterator fit = appObject->mesh->facets_begin();
			Mesh::Facet_iterator fit_end = appObject->mesh->facets_end();
			for (; fit != fit_end; fit++)
				ts << fit->minDist() << '\n';
		}

		file.close();
		return true;
	}

	return false;
}

bool ExporterMeshFeature::saveMaximalDiameterVolume(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		if (onVertices)
		{
			Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin();
			Mesh::Vertex_iterator vit_end = appObject->mesh->vertices_end();
			for (; vit != vit_end; vit++)
				ts << vit->volumeMaxDist() << '\n';
		}
		else
		{
			Mesh::Facet_iterator fit = appObject->mesh->facets_begin();
			Mesh::Facet_iterator fit_end = appObject->mesh->facets_end();
			for (; fit != fit_end; fit++)
				ts << fit->volumeMaxDist() << '\n';
		}

		file.close();
		return true;
	}

	return false;
}

bool ExporterMeshFeature::saveMinimalDiameterVolume(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		if (onVertices)
		{
			Mesh::Vertex_iterator vit = appObject->mesh->vertices_begin();
			Mesh::Vertex_iterator vit_end = appObject->mesh->vertices_end();
			for (; vit != vit_end; vit++)
				ts << vit->volumeMinDist() << '\n';
		}
		else
		{
			Mesh::Facet_iterator fit = appObject->mesh->facets_begin();
			Mesh::Facet_iterator fit_end = appObject->mesh->facets_end();
			for (; fit != fit_end; fit++)
				ts << fit->volumeMinDist() << '\n';
		}

		file.close();
		return true;
	}

	return false;
}

bool ExporterMeshFeature::saveDDegreeText(const QString& featuresFileName, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		for (Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin(); fit != appObject->mesh->facets_end(); fit++)
			ts << fit->deformationDegree() << '\n';

		file.close();
		return true;
	}

	return false;
}


bool ExporterMeshFeature::saveNDDegreeText(const QString& featuresFileName, AppObject* appObject)
{
	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		for (Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin(); fit != appObject->mesh->facets_end(); fit++)
			ts << fit->normalizeDeformationDegree() << '\n';

		file.close();
		return true;
	}

	return false;
}

bool ExporterMeshFeature::saveNormals(const QString& featuresFileName, const bool onVertices, AppObject* appObject)
{
	if (onVertices)
	{
		QFile file(featuresFileName);
		if (file.open(QIODevice::WriteOnly))
		{
			QTextStream ts(&file);
			ts.setRealNumberPrecision(10);

			Mesh::Vertex_const_iterator vit = appObject->mesh->vertices_begin();
			Mesh::Vertex_const_iterator vit_end = appObject->mesh->vertices_end();
			for (; vit != vit_end; vit++)
				ts << vit->normal().x() << " " << vit->normal().y() << " " << vit->normal().z() << '\n';

			file.close();
			return true;
		}

		return false;
	}
	else
	{
		QFile file(featuresFileName);
		if (file.open(QIODevice::WriteOnly))
		{
			QTextStream ts(&file);
			ts.setRealNumberPrecision(10);

			Mesh::Facet_const_iterator fit = appObject->mesh->facets_begin();
			Mesh::Facet_const_iterator fit_end = appObject->mesh->facets_end();
			for (; fit != fit_end; fit++)
				ts << fit->normal().x() << " " << fit->normal().y() << " " << fit->normal().z() << '\n';

			file.close();
			return true;
		}

		return false;
	}
}


bool ExporterMeshFeature::saveNormals(const QString& featuresFileName, vector<Vector_3> normals)
{
	cout << "size = " << normals.size() << endl;

	QFile file(featuresFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);

		for (int index = 0; index < normals.size(); index++)
		{
			Vector_3 n = normals[index];
			ts << n.x() << " " << n.y() << " " << n.z() << '\n';
		}

		file.close();
		return true;
	}

	return false;
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

bool ExporterArray::exportArray1D(const TNT::Array1D<int>& m, const QString& exportFileName)
{
	QFile file(exportFileName);

	if (!file.open(QIODevice::WriteOnly))
		return false;

	QTextStream ts(&file);

	int rowDim = m.dim1();
	for (int rowIndex = 0; rowIndex < rowDim; rowIndex++)
	{
		ts << m[rowIndex];
		if (rowIndex != (rowDim - 1))
			ts << '\n';
	}

	file.close();
	return true;
}

bool ExporterArray::exportArray1D(const TNT::Array1D<double>& m, const QString& exportFileName)
{
	QFile file(exportFileName);

	if (!file.open(QIODevice::WriteOnly))
		return false;

	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	int rowDim = m.dim1();

	for (int rowIndex = 0; rowIndex < rowDim; rowIndex++)
	{
		ts << m[rowIndex];
		if (rowIndex != (rowDim - 1))
			ts << '\n';
	}

	file.close();
	return true;
}

bool ExporterArray::exportArray2D(const TNT::Array2D<double>& m, const QString& exportFileName)
{
	QFile file(exportFileName);

	if (!file.open(QIODevice::WriteOnly))
		return false;

	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	int rowDim = m.dim1();
	int colDim = m.dim2();

	for (int rowIndex = 0; rowIndex < rowDim; rowIndex++)
	{
		for (int columnIndex = 0; columnIndex < colDim; columnIndex++)
		{
			ts << m[rowIndex][columnIndex];
			if (columnIndex != (colDim - 1))
				ts << ' ';
		}
		ts << '\n';
	}

	file.close();
	return true;
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


bool ExporterLog::exportRelationSDFWatermarkEmbeddingLog(
	const QString& exportLogFileName, 
	const int gridSize, const int numCones, const float& coneSeparationDegree, const double coneSeparationRadian, const int raysInCone, 
	const bool gaussianWeights, const bool removeOutlier, const NormalizeTypes& normalizeTypes, const int logAlpha,
	const double minSDF, const double maxSDF, const double avgSDF, 
	const double minNSDF, const double maxNSDF, const double avgNSDF,
	const double minStdDevSDF, const double maxStdDevSDF, const double avgStdDevSDF, 
	const double minStdDevNSDF, const double maxStdDevNSDF, const double avgStdDevNSDF, 
	const double range, const double interval,
	const double minDeformedSDF, const double maxDeformedSDF, const double avgDeformedSDF, 
	const double minDeformedNSDF, const double maxDeformedNSDF, const double avgDeformedNSDF,
	const double minDeformedStdDevSDF, const double maxDeformedStdDevSDF, const double avgDeformedStdDevSDF, 
	const double minDeformedStdDevNSDF, const double maxDeformedStdDevNSDF, const double avgDeformedStdDevNSDF,
	const double deformedRange, const double deformedBinInterval,
	const long& seed, const int totalWatermarkBits, const int numBins, const int robustnessFactor, const SearchSpace& searchSpace,
	const int maxCycles, const int particleNumbers,
	const int totalVertices,  const int totalFacets, const int totalEdges,
	const vector<vector<int> > originalVertexIndexBins, const vector<vector<int> > deformedVertexIndexBins,
	const TNT::Array1D<int>& precheckWatermarkSequence, 
	const TNT::Array1D<int>& originalWatermarkSequence, 
	const TNT::Array1D<int>& detectedWatermarkSequence, 
	const double preBitErrorRate, const double preCorrelationCoefficient,
	const double bitErrorRate, const double correlationCoefficient,
	const double SNR, const double PSNR, 
	const double MRMS, const double MRMSwrtBB, 
	const double Hausdorff, const double HausdorffwrtBB,
	const int timerEmbed)
{
	QFile file(exportLogFileName);
	if (!file.open(QIODevice::WriteOnly))
		return false;

	cout << "\nStart to export log file" << endl;
	cout << qPrintable(exportLogFileName) << endl;

	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	ts << "#Calculate shape diameter function parameters: " << endl;
	ts << "\tGrid size: " << gridSize << endl;
	ts << "\tNum cones: " << numCones << endl;
	ts << "\tCones seperation degree: " << coneSeparationDegree << endl;
	ts << "\tCones seperation radian: " << coneSeparationRadian << endl;
	ts << "\tRays in cone: " << raysInCone << endl;

	ts << "\tUse gaussian weights: ";
	if (gaussianWeights)
		ts << "true" << endl;
	else
		ts << "false" << endl;

	ts << "\tRemove outlier: ";
	if (removeOutlier)
		ts << "true" << endl;
	else
		ts << "false" << endl;

	ts << "\tNormalize type: ";
	if (normalizeTypes == None)
		ts << "None" << endl;
	else if (normalizeTypes == MinMax)
		ts << "MinMax" << endl;
	else if (normalizeTypes == Log)
		ts << "Log" << endl;

	ts << "\tLog Alpha: " << logAlpha << endl;

	ts << "\tMin SDF: " << minSDF << endl;
	ts << "\tMax SDF: " << maxSDF << endl;
	ts << "\tAvg SDF: " << avgSDF << endl;
	ts << "\tMin NSDF: " << minNSDF << endl;
	ts << "\tMax NSDF: " << maxNSDF << endl;
	ts << "\tAvg NSDF: " << avgNSDF << endl;
	ts << "\tMin Std Dev SDF: " << minStdDevSDF << endl;
	ts << "\tMax Std Dev SDF: " << maxStdDevSDF << endl;
	ts << "\tAvg Std Dev SDF: " << avgStdDevSDF << endl;
	ts << "\tMin Std Dev SDF: " << minStdDevNSDF << endl;
	ts << "\tMax Std Dev SDF: " << maxStdDevNSDF << endl;
	ts << "\tAvg Std Dev SDF: " << avgStdDevNSDF << endl;
	ts << "\tRange: " << range << endl;
	ts << "\tBin Interval: " << interval << endl;

	ts << "\tMin Deformed SDF: " << minDeformedSDF << endl;
	ts << "\tMax Deformed SDF: " << maxDeformedSDF << endl;
	ts << "\tAvg Deformed SDF: " << avgDeformedSDF << endl;
	ts << "\tMin Deformed NSDF: " << minDeformedNSDF << endl;
	ts << "\tMax Deformed NSDF: " << maxDeformedNSDF << endl;
	ts << "\tAvg Deformed NSDF: " << avgDeformedNSDF << endl;
	ts << "\tMin Deformed Std Dev SDF: " << minDeformedStdDevSDF << endl;
	ts << "\tMax Deformed Std Dev SDF: " << maxDeformedStdDevSDF << endl;
	ts << "\tAvg Deformed Std Dev SDF: " << avgDeformedStdDevSDF << endl;
	ts << "\tMin Deformed Std Dev SDF: " << minDeformedStdDevNSDF << endl;
	ts << "\tMax Deformed Std Dev SDF: " << maxDeformedStdDevNSDF << endl;
	ts << "\tAvg Deformed Std Dev SDF: " << avgDeformedStdDevNSDF << endl;
	ts << "\tDeformed Range: " << deformedRange << endl;
	ts << "\tDeformed Bin Interval: " << deformedBinInterval << endl;

	ts << "#Embedding parameters: " << endl;
	ts << "\tSeed: " << seed << endl;
	ts << "\tTotal bits: " << totalWatermarkBits << endl;		
	ts << "\tNum bins: " << numBins << endl;
	ts << "\tRobustness factor: " << robustnessFactor << endl;
	ts << "\tSearch space : ";

	if (searchSpace == MIN_DIST)
		ts << "MIN_DIST" << endl;
	else if (searchSpace == AVG_DIST)
		ts << "AVG_DIST" << endl;
	else if (searchSpace == MAX_DIST)
		ts << "MAX_DIST" << endl;

	ts << "\tMax cycles: " << maxCycles << endl;
	ts << "\tParticle numbers: " << particleNumbers << endl;

	ts << "#Model information: " << endl;
	ts << "\tTotal vertices: " << totalVertices << endl;
	ts << "\tTotal facets: " << totalFacets << endl;
	ts << "\tTotal edges: " << totalEdges << endl;
	
	ts << "#Results of watermark embedding procedure: " << endl;
	ts << "\tPrecheck watermark sequence: ";
	ts << "[";
	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
	{
		ts << precheckWatermarkSequence[bitIndex];
		if (bitIndex != (totalWatermarkBits - 1))
			ts << ",";
	}
	ts << "]" << endl;

	ts << "\tEmbedded watermark sequence: ";
	ts << "[";
	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
	{
		ts << originalWatermarkSequence[bitIndex];
		if (bitIndex != (totalWatermarkBits - 1))
			ts << ",";
	}
	ts << "]" << endl;

	ts << "\tDetected watermark sequence: ";
	ts << "[";
	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
	{
		ts << detectedWatermarkSequence[bitIndex];
		if (bitIndex != (totalWatermarkBits - 1))
			ts << ",";
	}
	ts << "]" << endl;

	ts << "\tpre-BER = " << preBitErrorRate * 100 << "%" << " (total bits = " << totalWatermarkBits << ")" << endl;
	ts << "\tpre-CC = " << preCorrelationCoefficient << endl;
	ts << "\tBER = " << bitErrorRate * 100 << "%" << " (total bits = " << totalWatermarkBits << ")" << endl;
	ts << "\tCC = " << correlationCoefficient << endl;

	ts << "#Quality information: " << endl;
	ts << "\tSNR: " << SNR << " (db)" << endl;
	ts << "\tPSNR: " << PSNR << " (db)" << endl;
	ts << "\tMRMS: " << MRMS << endl;
	ts << "\tMRMSwrtBB: " << MRMSwrtBB << endl;
	ts << "\tHausdorff: " << Hausdorff << endl;
	ts << "\tHausdorffwrtBB: " << HausdorffwrtBB << endl;

	ts << "#Other related information: " << endl;
	ts << "\tDuration time (sec): " << timerEmbed / 1000.0 << endl;

	file.close();

	cout << "Export log file done!\n\n";
	
	return true;
}


bool ExporterLog::exportRelationSDFWatermarkDetectingLog(
	const QString& exportLogFileName, 
	const int gridSize, const int numCones, const float& coneSeparationDegree, const double coneSeparationRadian, const int raysInCone, 
	const bool gaussianWeights, const bool removeOutlier, const NormalizeTypes& normalizeTypes, const int logAlpha,
	const double minSDF, const double maxSDF, const double avgSDF, 
	const double minNSDF, const double maxNSDF, const double avgNSDF, 
	const double minStdDevSDF, const double maxStdDevSDF, const double avgStdDevSDF,
	const double minStdDevNSDF, const double maxStdDevNSDF, const double avgStdDevNSDF,
	const int numBins, const double interval,
	const long& seed, const int totalWatermarkBits,
	const int totalVertices,  const int totalFacets, const int totalEdges,
	const vector<vector<int> > vertexIndexBins,
	const int* originalWatermarkSequence, const int* detectedWatermarkSequence, 
	const double bitErrorRate, const double correlationCoefficient,
	const int timerDetect)
{
	QFile file(exportLogFileName);
	if (!file.open(QIODevice::WriteOnly))
		return false;

	cout << "\nStart to export log file" << endl;
	cout << qPrintable(exportLogFileName) << endl;

	QTextStream ts(&file);
	ts.setRealNumberPrecision(10);

	ts << "#Calculate shape diameter function parameters: " << endl;
	ts << "\tGrid size: " << gridSize << endl;
	ts << "\tNum cones: " << numCones << endl;
	ts << "\tCones seperation degree: " << coneSeparationDegree << endl;
	ts << "\tCones seperation radian: " << coneSeparationRadian << endl;
	ts << "\tRays in cone: " << raysInCone << endl;

	ts << "\tUse gaussian weights: ";
	if (gaussianWeights)
		ts << "true" << endl;
	else
		ts << "false" << endl;

	ts << "\tRemove outlier: ";
	if (removeOutlier)
		ts << "true" << endl;
	else
		ts << "false" << endl;

	ts << "\tNormalize type: ";
	if (normalizeTypes == None)
		ts << "None" << endl;
	else if (normalizeTypes == MinMax)
		ts << "MinMax" << endl;
	else if (normalizeTypes == Log)
		ts << "Log" << endl;

	ts << "\tLog Alpha: " << logAlpha << endl;

	ts << "\tMin SDF: " << minSDF << endl;
	ts << "\tMax SDF: " << maxSDF << endl;
	ts << "\tAvg SDF: " << avgSDF << endl;
	ts << "\tMin NSDF: " << minNSDF << endl;
	ts << "\tMax NSDF: " << maxNSDF << endl;
	ts << "\tAvg NSDF: " << avgNSDF << endl;
	ts << "\tMin Std Dev SDF: " << minStdDevSDF << endl;
	ts << "\tMax Std Dev SDF: " << maxStdDevSDF << endl;
	ts << "\tAvg Std Dev SDF: " << avgStdDevSDF << endl;
	ts << "\tMin Std Dev NSDF: " << minStdDevNSDF << endl;
	ts << "\tMax Std Dev NSDF: " << maxStdDevNSDF << endl;
	ts << "\tAvg Std Dev NSDF: " << avgStdDevNSDF << endl;
	
	ts << "#Embedding parameters: " << endl;
	ts << "\tNum bins: " << numBins << endl;
	ts << "\tInterval: " << interval << endl;
	ts << "\tSeed: " << seed << endl;
	ts << "\tTotal bits: " << totalWatermarkBits << endl;		

	ts << "#Model information: " << endl;
	ts << "\tTotal vertices: " << totalVertices << endl;
	ts << "\tTotal facets: " << totalFacets << endl;
	ts << "\tTotal edges: " << totalEdges << endl;
	
	ts << "#Results of watermark detecting procedure: " << endl;
	ts << "\tAll Bins size: " << endl;
	ts << "\t[" << vertexIndexBins[0].size();
	for (int binIndex = 1; binIndex < (numBins-1); binIndex++)
	{
		if ((binIndex + 1) == (numBins - 1))
			break;
		
		ts << " , " << vertexIndexBins[binIndex].size();
	}
	ts << "]" << endl;

	ts << "\tEmbedded watermark sequence: ";
	ts << "[";
	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
	{
		ts << originalWatermarkSequence[bitIndex];
		if (bitIndex != (totalWatermarkBits - 1))
			ts << ",";
	}
	ts << "]" << endl;

	ts << "\tDetected watermark sequence: ";
	ts << "[";
	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
	{
		ts << detectedWatermarkSequence[bitIndex];
		if (bitIndex != (totalWatermarkBits - 1))
			ts << ",";
	}
	ts << "]" << endl;

	ts << "\tBER = " << bitErrorRate * 100 << "%" << " (total bits = " << totalWatermarkBits << ")" << endl;
	ts << "\tCC = " << correlationCoefficient << endl;

	ts << "#Other related information: " << endl;
	ts << "\tDuration time (sec): " << timerDetect / 1000.0 << endl;

	file.close();

	cout << "Export log file done!\n\n";
	
	return true;
}