#ifndef __EXPORTER_H_
#define __EXPORTER_H_

#include "stdafx.h"

#include <q3textstream.h>

#include "AppManager.h"
#include "AppObject.h"
#include "Watermark.h"

class FileUtils{

public:

	static bool save(const QString& fileName, const int fileType, const QString& description, const char* data, const unsigned int dataSize);
	static bool load(const QString& fileName, const int fileType, QString& description, char* data, unsigned int dataSize);

};

class ExporterMeshFile {

public:

	static bool write_obj(const QString& filename, AppObject* appObject);
	static bool write_off(const QString& filename, AppObject* appObject);

	//bool saveObjWithFaceSDF(const QString& filename, AppObject* appObject);
	//bool saveObjWithFaceColor(const QString& filename, AppObject* appObject, FaceColorizer& cizr, bool group);

};

class ExporterMeshFeature {

public:

	bool saveSDFVertices(const QString& featuresFileName, AppObject* appObject);
	bool saveSDFFacets(const QString& featuresFileName, AppObject* appObject);

	bool saveNSDFVertices(const QString& featuresFileName, AppObject* appObject);
	bool saveNSDFFacets(const QString& featuresFileName, AppObject* appObject);

	bool saveSDFText(const QString& exportFileName, const bool onVertices, AppObject* appObject);
	bool saveSDFDiffText(const QString& featuresFileName, const bool onVertices, AppObject* appObject);

	bool saveNSDFText(const QString& exportFileName, const bool onVertices, AppObject* appObject);
	bool saveNSDFDiffText(const QString& featuresFileName, const bool onVertices, AppObject* appObject);

	bool saveVSIRefPoints(const QString& featuresFileName, AppObject* appObject);
	bool saveVSIText(const QString& featuresFileName, const bool onVertices, AppObject* appObject);

	bool saveMaximalDiameter(const QString& featuresFileName, const bool onVertices, AppObject* appObject);
	bool saveMinimalDiameter(const QString& featuresFileName, const bool onVertices, AppObject* appObject);
	bool saveMaximalDiameterVolume(const QString& featuresFileName, const bool onVertices, AppObject* appObject);
	bool saveMinimalDiameterVolume(const QString& featuresFileName, const bool onVertices, AppObject* appObject);

	bool saveDDegreeText(const QString& featuresFileName, AppObject* appObject);
	bool saveNDDegreeText(const QString& featuresFileName, AppObject* appObject);	

	bool saveNormals(const QString& featuresFileName, const bool onVertices, AppObject* appObject);
	bool saveNormals(const QString& featuresFileName, vector<Vector_3> normals);

};

class ExporterArray {

public:

	bool exportArray1D(const TNT::Array1D<int>& m, const QString& exportFileName);
	bool exportArray1D(const TNT::Array1D<double>& m, const QString& exportFileName);
	bool exportArray2D(const TNT::Array2D<double>& m, const QString& exportFileName);

};

class ExporterLog {

public:

	bool exportRelationSDFWatermarkEmbeddingLog(
		const QString& exportLogFileName, 
		const int gridSize, const int numCones, const float& coneSeparationDegree, const double coneSeparationRadian, const int raysInCone, 
		const bool gaussianWeights, const bool removeOutlier, const NormalizeTypes& normalizeTypes, const int logAlpha,
		const double minSDF, const double maxSDF, const double meanSDF,
		const double minNSDF, const double maxNSDF, const double meanNSDF,
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
		const int timerEmbed);

	bool exportRelationSDFWatermarkDetectingLog(
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
		const int timerDetect);

};

#endif