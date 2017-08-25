#include "Watermark.h"

#include <qdatetime.h>

#include <QProgressDialog>

#include <omp.h>

#include <iomanip>

#include "RayIntersect.h"

///// Particle Swarm Optimization /////
const int dim = 3;
const double c1 = 2.05;
const double c2 = 2.05;
/// POS-IW ///
const double weightMax = 0.9;
const double weightMin = 0.4;
/// POS-CF ///
const double factorSum = c1+c2;
const double constrictionFactor = 2 / abs(2 - factorSum - sqrt(factorSum * factorSum - 4 * factorSum));


/**
  * Psuedo-random generator to generator a watermark sequence with length = totalBits and seed
  */
int* PseudoRandomSequenceGenerator::generate(long seed, int totalWatermarkBits)
{
	srand(seed);

	int* watermarkSequence = new int[totalWatermarkBits];

	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
		watermarkSequence[bitIndex] = (rand() % 2 == 0 ? 1 : 0);
	
	return watermarkSequence;
}

void PseudoRandomSequenceGenerator::generate(long seed, TNT::Array1D<int>& watermarkSequence)
{
	srand(seed);

	const int totalWatermarkBits = watermarkSequence.dim();
	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
		watermarkSequence[bitIndex] = (rand() % 2 == 0 ? 1 : 0);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


SDFWatermarkingAlgorithm::SDFWatermarkingAlgorithm()
{
	//ctor
}

SDFWatermarkingAlgorithm::~SDFWatermarkingAlgorithm()
{
	//dtor
}


/** 
  * Calculate the Euclidean distance between two point
  */
double SDFWatermarkingAlgorithm::distance(Point_3 p1, Point_3 p2)
{
	double x = p1.x() - p2.x();
	double y = p1.y() - p2.y();
	double z = p1.z() - p2.z();

	return sqrt(x * x + y * y + z * z);
}

/**
  * Generate a pseudo-random number with range [lower, upper]
  */
double SDFWatermarkingAlgorithm::randomNumber(double lower, double upper)
{
	/// to [0,1] binInterval
	double value = (double) rand() / RAND_MAX;

	/// to [lower, upper] binInterval
	double binInterval = upper - lower;
	double new_value = binInterval * value + lower;

	return new_value;
}


/** 
  * Generate a pesudo-random point in a ball with radius and center
  */

Point_3 SDFWatermarkingAlgorithm::randomPointInSphere(Point_3 center, double radius)
{
	double u = randomNumber(0, 1);
	double v = randomNumber(0, 1);
	double r = randomNumber(0, radius);

	double theta = 2 * M_PI * u;
	double phi = acos(2 * v - 1);

	double x = r * cos(theta) * sin(phi);
	double y = r * sin(theta) * sin(phi);
	double z = r * cos(phi);

	return Point_3(center.x() + x, center.y() + y, center.z() + z);
}


void SDFWatermarkingAlgorithm::getMinMaxSDF(
	const Mesh* mesh,
	const bool onVertices,
	double& minSDF, double& maxSDF)
{
	minSDF = FLT_MAX;
	maxSDF = 0.0;

	if (onVertices)
	{
		Mesh::Vertex_const_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_const_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
		{
			if (vit->volumeSDF() < minSDF) minSDF = vit->volumeSDF();
			if (vit->volumeSDF() > maxSDF) maxSDF = vit->volumeSDF();
		}
	}
	else
	{
		Mesh::Facet_const_iterator fit = mesh->facets_begin();
		Mesh::Facet_const_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
		{
			if (fit->volumeSDF() < minSDF) minSDF = fit->volumeSDF();
			if (fit->volumeSDF() > maxSDF) maxSDF = fit->volumeSDF();
		}
	}
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


bool SDFWatermarkingAlgorithm::embed(
	AppObject* appObject, Mesh* originalMesh, Mesh* watermarkedMesh, 
	const bool onVertices, const bool multiThreaded, const int threadsNum,
	const long seed, const int totalWatermarkBits, const int robustnessFactor, const SearchSpace searchSpace,
	const int maxCycles, const int particleNumbers,
	const int gridSize, const int numCones, const float coneSeparationDegree, const int raysInCone, 
	const bool removeSameDirectionIntersection, const bool gaussianWeights, const bool removeOutlier,
	const NormalizeTypes normalizeTypes, const int logAlpha, 
	const bool smoothing, const int smoothingIterations,
	const bool debugMode, const bool showDialog)
{
	cout << endl << endl;

	QTime timer;

	timer.start();

	const double coneSeparationRadian = coneSeparationDegree * DEG2RAD;
	const int numBins = 2 * totalWatermarkBits + 2;

	//////////////////////////////////////////////////////////////////////////
	// Basic model information
	//////////////////////////////////////////////////////////////////////////
	const int totalVertices = watermarkedMesh->size_of_vertices();
	const int totalFacets = watermarkedMesh->size_of_facets();
	const int totalEdges = watermarkedMesh->size_of_halfedges() / 2; 

	//////////////////////////////////////////////////////////////////////////
	// Find minimum SDF and maximum SDF value
	//////////////////////////////////////////////////////////////////////////
	const double minOriginalSDF = watermarkedMesh->getMinVolume();
	const double maxOriginalSDF = watermarkedMesh->getMaxVolume();
	const double avgOriginalSDF = watermarkedMesh->getAvgVolume();
	const double minOriginalNSDF = watermarkedMesh->getMinNVolume();
	const double maxOriginalNSDF = watermarkedMesh->getMaxNVolume();
	const double avgOriginalNSDF = watermarkedMesh->getAvgNVolume();

	const double minOriginalStdDevSDF = watermarkedMesh->getMinStdDevVolume();
	const double maxOriginalStdDevSDF = watermarkedMesh->getMaxStdDevVolume();
	const double avgOriginalStdDevSDF = watermarkedMesh->getAvgStdDevVolume();
	const double minOriginalStdDevNSDF = watermarkedMesh->getMinStdDevNVolume();
	const double maxOriginalStdDevNSDF = watermarkedMesh->getMaxStdDevNVolume();
	const double avgOriginalStdDevNSDF = watermarkedMesh->getAvgStdDevNVolume();

	const double originalRange = maxOriginalStdDevNSDF - minOriginalStdDevNSDF;
	const double originalBinInterval = originalRange / numBins;

	cout.precision(10);
	cout << "minOriginalSDF = " << minOriginalSDF << " , maxOriginalSDF = " << maxOriginalSDF << " , avgOriginalSDF = " << avgOriginalSDF << endl;
	cout << "minOriginalNSDF = " << minOriginalNSDF << " , maxOriginalNSDF = " << maxOriginalNSDF << " , avgOriginalNSDF = " << avgOriginalNSDF << endl;
	cout << "minOriginalStdDevSDF = " << minOriginalStdDevSDF << " , maxOriginalStdDevSDF = " << maxOriginalStdDevSDF << " , avgStdDevSDF = " << avgOriginalStdDevSDF << endl;
	cout << "minStdDevNSDF = " << minOriginalStdDevNSDF << " , maxStdDevNSDF = " << maxOriginalStdDevNSDF << " , avgStdDevNSDF = " << avgOriginalStdDevNSDF << endl;
	cout << "binInterval = " << originalBinInterval << endl;

	//////////////////////////////////////////////////////////////////////////
	// Allocate all the vertices into bins
	//////////////////////////////////////////////////////////////////////////
	vector<vector<int> > originalVertexIndexBins(numBins);
	for (Mesh::Vertex_iterator vit = watermarkedMesh->vertices_begin(); vit != watermarkedMesh->vertices_end(); vit++)
	{
		double originalNSDF = vit->volumeNSDF();
		if (originalNSDF < minOriginalStdDevNSDF || originalNSDF > maxOriginalStdDevNSDF)
			continue;

		int originalBinIntervalIndex = (originalNSDF - minOriginalStdDevNSDF) / originalBinInterval;
		if (originalBinIntervalIndex < 0 || originalBinIntervalIndex > numBins)
			continue;

		if (originalBinIntervalIndex == numBins)
			originalVertexIndexBins[originalBinIntervalIndex-1].push_back(vit->index());
		else
			originalVertexIndexBins[originalBinIntervalIndex].push_back(vit->index());
	}

	//////////////////////////////////////////////////////////////////////////
	// Generate pseudo random watermark sequence
	//////////////////////////////////////////////////////////////////////////
	TNT::Array1D<int> originalWatermarkSequence = TNT::Array1D<int>(totalWatermarkBits);

	PseudoRandomSequenceGenerator generator;
	generator.generate(seed, originalWatermarkSequence);

	//////////////////////////////////////////////////////////////////////////
	// Precheck
	//////////////////////////////////////////////////////////////////////////
	TNT::Array1D<int> precheckWatermarkSequence = TNT::Array1D<int>(totalWatermarkBits);
	int watermarkIndex = 0;
	for (int binIndex = 1; binIndex < (numBins-1); binIndex+=2)
	{
		if ((binIndex + 1) == (numBins - 1))
			continue;

		int firstBinIndex = binIndex;
		int secondBinIndex = firstBinIndex + 1;

		int firstBinSize = originalVertexIndexBins[firstBinIndex].size();
		int secondBinSize = originalVertexIndexBins[secondBinIndex].size();

		if (firstBinSize <= secondBinSize)
		{
			precheckWatermarkSequence[watermarkIndex] = 1;
			//originalWatermarkSequence[watermarkIndex] = 1;
		}
		else
		{
			precheckWatermarkSequence[watermarkIndex] = 0;
			//originalWatermarkSequence[watermarkIndex] = 0;
		}

		watermarkIndex++;
	}

	cout << "original watermark sequence\n";
	cout << "[";
	for (int bitIndex = 0; bitIndex < totalWatermarkBits; bitIndex++)
	{
		cout << originalWatermarkSequence[bitIndex];
		if (bitIndex != (totalWatermarkBits - 1))
			cout << ",";
	}
	cout << "]\n";

	Statistic statistic;

	double preBitErrorRate = statistic.calculateBitErrorRate(originalWatermarkSequence, precheckWatermarkSequence);
	double preCorrelationCoefficient = statistic.calcuateCorrelationCoefficient(originalWatermarkSequence, precheckWatermarkSequence);

	printf("pre-BER = %g, pre-CC = %g\n", preBitErrorRate, preCorrelationCoefficient);

	//////////////////////////////////////////////////////////////////////////
	// Watermark embedding procedure
	//////////////////////////////////////////////////////////////////////////
	cout << "\n\nStart to embed watermark with numBins = " << numBins << endl;
	cout << "maxCycles = " << maxCycles << " , particleNumbers = " << particleNumbers << endl;

	meshEmbeddingProcedure(
			appObject, watermarkedMesh, onVertices,
			multiThreaded, threadsNum,
			originalWatermarkSequence, seed, totalWatermarkBits, robustnessFactor, searchSpace,
			maxCycles, particleNumbers, numBins,
			gridSize, numCones, coneSeparationDegree, coneSeparationRadian, raysInCone, 
			removeSameDirectionIntersection, gaussianWeights, removeOutlier, 
			normalizeTypes, logAlpha,
			smoothing, smoothingIterations,
			debugMode, showDialog);

	printf("Mesh Embedding Procedure Done!!!\n");

	int timerEmbed = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Export mesh
	//////////////////////////////////////////////////////////////////////////
	cout << "\tStart to export watermarked mesh model...\n";

	QString qBits;
	qBits.setNum(totalWatermarkBits);

	QString watermarkedMeshExportFileName = appObject->fileDir + appObject->fileSeparator + "Watermark_SDF_Relation_" + qBits + "bits_" + appObject->fileName + appObject->fileSuffix;
	watermarkedMesh->write_obj(watermarkedMeshExportFileName);

	cout << "\tExport watermarked mesh model done!\n\n";

	//////////////////////////////////////////////////////////////////////////
	// Compute SDF
	//////////////////////////////////////////////////////////////////////////
	watermarkedMesh->computeNormals();

	CalculateSDF csdf;
	csdf.go(
		watermarkedMesh, onVertices, debugMode, multiThreaded, threadsNum, 
		gridSize, numCones, coneSeparationDegree, raysInCone, 
		removeSameDirectionIntersection, gaussianWeights, removeOutlier, 
		normalizeTypes, logAlpha,
		smoothing, smoothingIterations);

	//////////////////////////////////////////////////////////////////////////
	// Allocate all the vertices into bins
	//////////////////////////////////////////////////////////////////////////
	const double minDeformedSDF = watermarkedMesh->getMinVolume();
	const double maxDeformedSDF = watermarkedMesh->getMaxVolume();
	const double avgDeformedSDF = watermarkedMesh->getAvgVolume();
	const double minDeformedNSDF = watermarkedMesh->getMinNVolume();
	const double maxDeformedNSDF = watermarkedMesh->getMaxNVolume();
	const double avgDeformedNSDF = watermarkedMesh->getAvgNVolume();

	const double minDeformedStdDevSDF = watermarkedMesh->getMinStdDevVolume();
	const double maxDeformedStdDevSDF = watermarkedMesh->getMaxStdDevVolume();
	const double avgDeformedStdDevSDF = watermarkedMesh->getAvgStdDevVolume();
	const double minDeformedStdDevNSDF = watermarkedMesh->getMinStdDevNVolume();
	const double maxDeformedStdDevNSDF = watermarkedMesh->getMaxStdDevNVolume();
	const double avgDeformedStdDevNSDF = watermarkedMesh->getAvgStdDevNVolume();

	const double deformedRange = maxDeformedStdDevNSDF - minDeformedStdDevNSDF;
	const double deformedBinInterval = deformedRange / numBins;

	vector<vector<int> > deformedVertexIndexBins(numBins);
	for (Mesh::Vertex_iterator vit = watermarkedMesh->vertices_begin(); vit != watermarkedMesh->vertices_end(); vit++)
	{
		double deformedNSDF = vit->volumeNSDF();
		if (deformedNSDF < minDeformedStdDevNSDF || deformedNSDF > maxDeformedStdDevNSDF)
			continue;
		
		int deformedBinIntervalIndex = (deformedNSDF - minDeformedStdDevNSDF) / deformedBinInterval;
		if (deformedBinIntervalIndex < 0 || deformedBinIntervalIndex > numBins)
			continue;

		if (deformedBinIntervalIndex == numBins)	
			deformedVertexIndexBins[deformedBinIntervalIndex-1].push_back(vit->index());
		else
			deformedVertexIndexBins[deformedBinIntervalIndex].push_back(vit->index());
	}

	//////////////////////////////////////////////////////////////////////////
	// Posterior-detecting procedure
	//////////////////////////////////////////////////////////////////////////
	TNT::Array1D<int> detectedWatermarkSequence = TNT::Array1D<int>(totalWatermarkBits);

	meshDetectigProcedure(watermarkedMesh, originalBinInterval, deformedVertexIndexBins, detectedWatermarkSequence);

	//////////////////////////////////////////////////////////////////////////
	// Statistical analysis
	//////////////////////////////////////////////////////////////////////////
	cout << "\nOriginal Watermark Sequence\n";
	cout << "[";
	for (int watermarkBitIndex = 0; watermarkBitIndex < totalWatermarkBits; watermarkBitIndex++)
	{
		cout << originalWatermarkSequence[watermarkBitIndex];
		if (watermarkBitIndex != (totalWatermarkBits - 1))
			cout << ",";
	}
	cout << "]" << endl;

	cout << "\nDetected Watermark Sequence\n";
	cout << "[";
	for (int watermarkBitIndex = 0; watermarkBitIndex < totalWatermarkBits; watermarkBitIndex++)
	{
		cout << detectedWatermarkSequence[watermarkBitIndex];
		if (watermarkBitIndex != (totalWatermarkBits - 1))
			cout << ",";
	}
	cout << "]" << endl;

	double bitErrorRate = statistic.calculateBitErrorRate(originalWatermarkSequence, detectedWatermarkSequence);
	double correlationCoefficient = statistic.calcuateCorrelationCoefficient(originalWatermarkSequence, detectedWatermarkSequence);

	cout << "BER = " << bitErrorRate << " (total totalWatermarkBits = " << totalWatermarkBits << ")\n";
	cout << "CC = " << correlationCoefficient << endl;

	//////////////////////////////////////////////////////////////////////////
	// Quality estimate
	//////////////////////////////////////////////////////////////////////////
	cout << "\nStart to estimate quality between original model and embedded model...";

	QString originalMeshFileName = appObject->fileDir + appObject->fileSeparator + appObject->fileName + appObject->fileSuffix;

	Quality quality;

	double SNR = quality.measureSNR(originalMesh, watermarkedMesh);
	double PSNR = quality.measurePSNR(originalMesh, watermarkedMesh);
	double MRMS;
	double MRMSwrtBB;
	double Hausdorff;
	double HausdorffwrtBB;

	quality.calculateGeoDistances(originalMeshFileName, watermarkedMeshExportFileName, MRMS, MRMSwrtBB, Hausdorff, HausdorffwrtBB);

	cout << "Done!" << endl;

	cout << "SNR = " << SNR << " , PSNR = " << PSNR << endl;
	cout << "MRMS = " << MRMS << " , MRMSwrtBB = " << MRMSwrtBB << endl;
	cout << "Hausdorff = " << Hausdorff << " , HausdorffwrtBB = " << HausdorffwrtBB << endl;

	//////////////////////////////////////////////////////////////////////////
	// Export log information
	//////////////////////////////////////////////////////////////////////////
	QString exportLogFileName = appObject->fileDir + appObject->fileSeparator + "WE_Watermark_SDF_Relation_" + qBits + "Bits_" + appObject->fileName + "_log.txt";
	
	ExporterLog exporterLog;
	exporterLog.exportRelationSDFWatermarkEmbeddingLog(
		exportLogFileName, 
		gridSize, numCones, coneSeparationDegree, coneSeparationRadian, raysInCone, 
		gaussianWeights, removeOutlier, normalizeTypes, logAlpha,
		minOriginalSDF, maxOriginalSDF, avgOriginalSDF, minOriginalNSDF, maxOriginalNSDF, avgOriginalNSDF, 
		minOriginalStdDevSDF, maxOriginalStdDevSDF, avgOriginalStdDevSDF, minOriginalStdDevNSDF, maxOriginalStdDevNSDF, avgOriginalStdDevNSDF, 
		originalRange, originalBinInterval,
		minDeformedSDF, maxDeformedSDF, avgDeformedSDF, minDeformedNSDF, maxDeformedNSDF, avgDeformedNSDF,
		minDeformedStdDevSDF, maxDeformedStdDevSDF, avgDeformedStdDevSDF, minDeformedStdDevNSDF, maxDeformedStdDevNSDF, avgDeformedStdDevNSDF,
		deformedRange, deformedBinInterval,
		seed, totalWatermarkBits, numBins, robustnessFactor, searchSpace,
		maxCycles, particleNumbers,
		totalVertices, totalFacets, totalEdges,
		originalVertexIndexBins, deformedVertexIndexBins,
		precheckWatermarkSequence, originalWatermarkSequence, detectedWatermarkSequence, 
		preBitErrorRate, preCorrelationCoefficient, 
		bitErrorRate, correlationCoefficient,
		SNR, PSNR, MRMS, MRMSwrtBB, Hausdorff, HausdorffwrtBB,
		timerEmbed);

	return true;
}


bool SDFWatermarkingAlgorithm::meshEmbeddingProcedure(
	AppObject* appObject, Mesh* watermarkedMesh, bool onVertices,
	const bool multiThreaded, const int threadsNum,
	const int* originalWatermarkSequence, const long seed, const int totalWatermarkBits, const int robustnessFactor, const SearchSpace searchSpace,
	const int maxCycles, const int particleNumbers, const int numBins,
	const int gridSize, const int numCones, const float coneSeparationDegree, const double coneSeparationRadian, const int raysInCone, 
	const bool removeSameDirectionIntersection, const bool gaussianWeights, const bool removeOutlier, 
	const NormalizeTypes normalizeTypes, const int logAlpha,
	const bool smoothing, const int smoothingIterations, 
	const bool debugMode, const bool showDialog)
{
	srand(seed);

	const int totalVertices = watermarkedMesh->size_of_vertices();
	const double delta = 1e-5;

	int threadsNumber;
	if (multiThreaded)
		threadsNumber = threadsNum;
	else
		threadsNumber = 1;

	cout << "\nStart to embed watermark with " << totalWatermarkBits << " bits...\n";

	vector<vector<int> > vertexIndexBins(numBins);
	vector<int> vertexIndexBin;
	int vertexIndexBinSize;
	for (int watermarkIndex = 0; watermarkIndex < totalWatermarkBits; watermarkIndex++)
	{
		const int watermarkBit = originalWatermarkSequence[watermarkIndex];

		const int firstBinIndex = 2 * watermarkIndex + 1;
		const int secondBinIndex = firstBinIndex + 1;

		///////////////////////////// Calculate SDF //////////////////////////////		
		watermarkedMesh->computeNormals();

		CalculateSDF csdf;
		csdf.go(
			watermarkedMesh, onVertices, debugMode, multiThreaded, threadsNum, 
			gridSize, numCones, coneSeparationDegree, raysInCone, 
			removeSameDirectionIntersection, gaussianWeights, removeOutlier, 
			normalizeTypes, logAlpha,
			smoothing, smoothingIterations);
		//////////////////////////////////////////////////////////////////////////

		const double minOriginalStdDevSDF = watermarkedMesh->getMinStdDevVolume();
		const double maxOriginalStdDevSDF = watermarkedMesh->getMaxStdDevVolume();
		const double avgOriginalStdDevSDF = watermarkedMesh->getAvgStdDevVolume();
		const double minOriginalStdDevNSDF = watermarkedMesh->getMinStdDevNVolume();
		const double maxOriginalStdDevNSDF = watermarkedMesh->getMaxStdDevNVolume();
		const double avgOriginalStdDevNSDF = watermarkedMesh->getAvgStdDevNVolume();

		const double originalRange = maxOriginalStdDevNSDF - minOriginalStdDevNSDF;
		const double originalBinInterval = originalRange / numBins;

		///////////////////////// Vertex set partition ///////////////////////////
		for (int index = 0; index < numBins; index++)
			vertexIndexBins[index].clear();

		for (Mesh::Vertex_iterator vit = watermarkedMesh->vertices_begin(); vit != watermarkedMesh->vertices_end(); vit++)
		{
			double originalNSDF = vit->volumeNSDF();
			if (originalNSDF < minOriginalStdDevNSDF || originalNSDF > maxOriginalStdDevNSDF)
				continue;

			int originalBinIntervalIndex = (originalNSDF - minOriginalStdDevNSDF) / originalBinInterval;
			if (originalBinIntervalIndex > numBins || originalBinIntervalIndex < 0)
				continue;

			if (originalBinIntervalIndex == numBins)
				vertexIndexBins[originalBinIntervalIndex-1].push_back(vit->index());
			else
				vertexIndexBins[originalBinIntervalIndex].push_back(vit->index());
		}
		//////////////////////////////////////////////////////////////////////////


		/////////////////////////// Terminal condition ///////////////////////////
		const int firstVertexIndexBinSize = vertexIndexBins[firstBinIndex].size();
		const int secondVertexIndexBinSize = vertexIndexBins[secondBinIndex].size();
		const int totalVertexIndexBinsSize = firstVertexIndexBinSize + secondVertexIndexBinSize;
		const int differenceSize = (int) ceil((firstVertexIndexBinSize + secondVertexIndexBinSize) / (double) robustnessFactor);

		if ((watermarkBit == 1) && ((secondVertexIndexBinSize - firstVertexIndexBinSize) >= differenceSize))
			continue;
		else if ((watermarkBit == 0) && ((firstVertexIndexBinSize - secondVertexIndexBinSize) >= differenceSize))
			continue;	
		//////////////////////////////////////////////////////////////////////////


		/////////////////////// Deform the vertices in bin ///////////////////////
		int shiftPointsNumber;

		int toShiftVertexIndexBinSize;
		int shiftedVertexIndexBinSize;

		vector<int> toShiftVertexIndexBin;
		vector<int> shiftedVertexIndexBin;

		double targetBinbinIntervalLeftValue;
		double targetBinbinIntervalRightValue;
		double targetLeftValue;
		double targetRightValue;

		if (watermarkBit == 1)
		{
			targetBinbinIntervalLeftValue = originalBinInterval * secondBinIndex;
			targetBinbinIntervalRightValue = originalBinInterval * (secondBinIndex + 1);

			//targetLeftValue = targetBinbinIntervalLeftValue;
			//targetRightValue = targetBinbinIntervalRightValue;
			targetLeftValue = targetBinbinIntervalLeftValue + (originalBinInterval / 4.0);
			targetRightValue = targetBinbinIntervalRightValue - (originalBinInterval / 4.0);

			if (firstVertexIndexBinSize <= secondVertexIndexBinSize)
			{
				printf("  current bit index = %d, current bit = 1, watermark bit = %d...Deformed!\n", watermarkIndex, watermarkBit);
				
				int difference = secondVertexIndexBinSize - firstVertexIndexBinSize;
				shiftPointsNumber = (int) ceil((differenceSize - difference) / 2.0);				
			}
			else
			{
				printf("  current bit index = %d, current bit = 0, watermark bit = %d...Deformed!\n", watermarkIndex, watermarkBit);

				shiftPointsNumber = firstVertexIndexBinSize - (int) ceil(totalVertexIndexBinsSize / 2.0) + differenceSize;
			}

			toShiftVertexIndexBin = vertexIndexBins[firstBinIndex];
			shiftedVertexIndexBin = vertexIndexBins[secondBinIndex];

			toShiftVertexIndexBinSize = toShiftVertexIndexBin.size();
			shiftedVertexIndexBinSize = shiftedVertexIndexBin.size();
		}
		else
		{
			targetBinbinIntervalLeftValue = originalBinInterval * firstBinIndex;
			targetBinbinIntervalRightValue = originalBinInterval * (firstBinIndex + 1);

			//targetLeftValue = targetBinbinIntervalLeftValue;
			//targetRightValue = targetBinbinIntervalRightValue;
			targetLeftValue = targetBinbinIntervalLeftValue + (originalBinInterval / 4.0);
			targetRightValue = targetBinbinIntervalRightValue - (originalBinInterval / 4.0);

			if (firstVertexIndexBinSize > secondVertexIndexBinSize)
			{
				printf("  current bit index = %d, current bit = 0, watermark bit = %d...Deformed!\n", watermarkIndex, watermarkBit);

				int difference = firstVertexIndexBinSize - secondVertexIndexBinSize;
				shiftPointsNumber = (int) ceil((differenceSize - difference) / 2.0);
			}
			else
			{
				printf("  current bit index = %d, current bit = 1, watermark bit = %d...Deformed!\n", watermarkIndex, watermarkBit);

				shiftPointsNumber = secondVertexIndexBinSize - (int) ceil(totalVertexIndexBinsSize / 2.0) + differenceSize;
			}
			
			toShiftVertexIndexBin = vertexIndexBins[secondBinIndex];
			shiftedVertexIndexBin = vertexIndexBins[firstBinIndex];

			toShiftVertexIndexBinSize = toShiftVertexIndexBin.size();
			shiftedVertexIndexBinSize = shiftedVertexIndexBin.size();
		}

		printf("  firstVertexIndexBinSize = %d, secondVertexIndexBinSize = %d\n\n", firstVertexIndexBinSize, secondVertexIndexBinSize);

		//////////////////////////////////////////////////////////////////////////
		// Find the shiftPointsNumber largest number
		//////////////////////////////////////////////////////////////////////////
		int sizeOfSortedSquareErrors = toShiftVertexIndexBinSize;
		vector<Mesh::Vertex_handle> sortedVertex;

		int sortedIndex = 0;
		for (int pointIndex = 0; pointIndex < toShiftVertexIndexBinSize; pointIndex++)
		{
			int vertexIndex = toShiftVertexIndexBin[pointIndex];
			Mesh::Vertex_handle vertex = watermarkedMesh->findVertex(vertexIndex);
			sortedVertex.push_back(vertex);
		}

		QuickSort qs;
		if (watermarkBit == 1)
			qs.sortVolumeSDFAscending(sortedVertex, 0, sizeOfSortedSquareErrors - 1);
		else
			qs.sortVolumeSDFDescending(sortedVertex, 0, sizeOfSortedSquareErrors - 1);

		cout << "toShiftVertexIndexBinSize = " << toShiftVertexIndexBinSize << " , shiftedVertexIndexBinSize = " << shiftedVertexIndexBinSize << endl;
		cout << "2.  shiftPointsNumber = " << shiftPointsNumber << endl;
		cout << "2.a toShiftVertexIndexBinSize = " << toShiftVertexIndexBinSize << endl;

		//////////////////////////////////////////////////////////////////////////
		// Deforming vertices
		//////////////////////////////////////////////////////////////////////////
		TNT::Array1D<int> sortedVertexIndex = TNT::Array1D<int>(sizeOfSortedSquareErrors);
		for (int index = 0; index < sizeOfSortedSquareErrors; index++)
			sortedVertexIndex[index] = sortedVertex[index]->index();

		//////////////////////////////////////////////////////////////////////////
		// Variables
		//////////////////////////////////////////////////////////////////////////
		TNT::Array2D<double> velocity = TNT::Array2D<double>(dim,particleNumbers); // current velocity of particle
		TNT::Array2D<double> position = TNT::Array2D<double>(dim,particleNumbers); // current position of particle

		TNT::Array2D<double> individualOptimalPosition = TNT::Array2D<double>(dim,particleNumbers);
		TNT::Array1D<double> individualOptimalFitness = TNT::Array1D<double>(particleNumbers);

		TNT::Array1D<double> globalOptimalPosition = TNT::Array1D<double>(dim);
		double globalOptimalFitness;

		TNT::Array1D<double> individualOptimalNSDF = TNT::Array1D<double>(particleNumbers);
		double globalOptimalNSDF;

		TNT::Array1D<double> sdfArray = TNT::Array1D<double>(particleNumbers);
		TNT::Array1D<double> nsdfArray = TNT::Array1D<double>(particleNumbers);

		TNT::Array1D<Point_3> deformedPointArray = TNT::Array1D<Point_3>(particleNumbers);
		TNT::Array1D<Vector_3> deformedNormalArray = TNT::Array1D<Vector_3>(particleNumbers);
		//////////////////////////////////////////////////////////////////////////
			
		int hasShiftPointsNumber = 0;
		for (int pointIndex = 0; pointIndex < sizeOfSortedSquareErrors; pointIndex++)
		{
			Mesh::Vertex_handle v = sortedVertex[pointIndex];

			const Point_3 originalPoint = v->point();		
			const double originalNSDF = v->volumeNSDF();

			//////////////////////////////////////////////////////////////////////////
			// Get min and max distance of it`s 1-ring neighbor
			//////////////////////////////////////////////////////////////////////////
			double dist;
			double minDist = FLT_MAX;
			double maxDist = -FLT_MAX;
			double avgDist = 0;
			Mesh::Halfedge_around_vertex_circulator c = v->vertex_begin();
			do {
				Mesh::Vertex_handle vNeighbor = c->opposite()->vertex();
				Point_3 pNeighborPoint = vNeighbor->point();

				dist = distance(originalPoint, pNeighborPoint);

				avgDist += dist;

				if (dist < minDist) minDist = dist;
				if (dist > maxDist) maxDist = dist;
			} while (++c != v->vertex_begin());

			avgDist /= v->degree();
			
			double searchSpaceDist;
			switch (searchSpace)
			{
				case MIN_DIST: searchSpaceDist = minDist; break;
				case AVG_DIST: searchSpaceDist = avgDist; break;
				case MAX_DIST: searchSpaceDist = maxDist; break;
			}

			//printf("minDist = %g, avgDist = %g, maxDist = %g\n", minDist, avgDist, maxDist);
		
			//////////////////////////////////////////////////////////////////////////
			// Initialization
			//////////////////////////////////////////////////////////////////////////
			#pragma omp parallel num_threads(threadsNumber)
			{
				#pragma omp for
				for (int j = 0; j < particleNumbers; j++) /// Initialize Points ///
				{
					Point_3 randomPoint = randomPointInSphere(originalPoint, searchSpaceDist);
					
					deformedPointArray[j] = randomPoint;

					position[0][j] = individualOptimalPosition[0][j] = randomPoint.x();
					position[1][j] = individualOptimalPosition[1][j] = randomPoint.y();
					position[2][j] = individualOptimalPosition[2][j] = randomPoint.z();

					velocity[0][j] = randomNumber(-searchSpaceDist / 20.0, searchSpaceDist / 20.0);
					velocity[1][j] = randomNumber(-searchSpaceDist / 20.0, searchSpaceDist / 20.0);
					velocity[2][j] = randomNumber(-searchSpaceDist / 20.0, searchSpaceDist / 20.0);

					double newRadius = distance(originalPoint, randomPoint);

					if (newRadius > searchSpaceDist)
						cout << "\tError on initialize points: " << newRadius << "," << searchSpaceDist << endl;
				}
			} /// end OpenMP

			TNT::Array1D<TNT::Array3D<RayIntersectCell*>> gridArray = TNT::Array1D<TNT::Array3D<RayIntersectCell*>>(particleNumbers);
			TNT::Array1D<triangleVector_t> triangleVectorArray = TNT::Array1D<triangleVector_t>(particleNumbers);

			TNT::Array1D<double> xminArray = TNT::Array1D<double>(particleNumbers);
			TNT::Array1D<double> xmaxArray = TNT::Array1D<double>(particleNumbers);
			TNT::Array1D<double> yminArray = TNT::Array1D<double>(particleNumbers);
			TNT::Array1D<double> ymaxArray = TNT::Array1D<double>(particleNumbers);
			TNT::Array1D<double> zminArray = TNT::Array1D<double>(particleNumbers);
			TNT::Array1D<double> zmaxArray = TNT::Array1D<double>(particleNumbers);

			TNT::Array1D<double> xspreadArray = TNT::Array1D<double>(particleNumbers);
			TNT::Array1D<double> yspreadArray = TNT::Array1D<double>(particleNumbers);
			TNT::Array1D<double> zspreadArray = TNT::Array1D<double>(particleNumbers);

			//cout << "define parameters done!" << endl;		

			/// Initialize RayIntersect ///
			for (int j = 0; j < particleNumbers; j++)
			{
				//cout << "j = " << j << " , particleNumbers = " << particleNumbers << endl;

				v->point() = deformedPointArray[j]; // set vertex position
				watermarkedMesh->computeNormals();
				deformedNormalArray[j] = v->normal();

				RayIntersect rayIntersect;
				rayIntersect.Init(*watermarkedMesh, gridArray[j], triangleVectorArray[j],
									xminArray[j], xmaxArray[j], yminArray[j], ymaxArray[j], zminArray[j], zmaxArray[j],
									xspreadArray[j], yspreadArray[j], zspreadArray[j], 
									gridSize, debugMode);

				//cout << "j = " << j << " , particleNumbers = " << particleNumbers << endl;
			}

			//cout << "rayIntersect done!" << endl;

			#pragma omp parallel num_threads(threadsNumber)
			{
				#pragma omp for
				for (int j = 0; j < particleNumbers; j++) /// Calculate SDF ///
				{	
					CalculateSDF csdf;

					sdfArray[j] = csdf.computeVolumeSDF(
						deformedPointArray[j], deformedNormalArray[j], gridSize, 
						numCones, coneSeparationDegree, coneSeparationRadian, raysInCone,
						removeSameDirectionIntersection, gaussianWeights, removeOutlier,
						gridArray[j], triangleVectorArray[j],
						xminArray[j], xmaxArray[j], yminArray[j], ymaxArray[j], zminArray[j], zmaxArray[j],
						xspreadArray[j], yspreadArray[j], zspreadArray[j], 
						debugMode);

					if (sdfArray[j] == 0)
					{
						double smoothedValue = 0.0;
						double smoothedWeights = 0.0;
					
						Mesh::Halfedge_around_vertex_circulator c = v->vertex_begin();
						do {
							smoothedValue += c->opposite()->vertex()->volumeSDF();
							smoothedWeights += 1.0;
						} while (++c != v->vertex_begin());

						sdfArray[j] = safeDiv(smoothedValue, smoothedWeights);
					}

					if (normalizeTypes != None)
					{
						//const double theDivider = 1 / (maxOriginalSDF - minOriginalSDF);
						const double theDivider = 1 / (avgOriginalStdDevSDF - minOriginalStdDevSDF);

						nsdfArray[j] = (sdfArray[j] - minOriginalStdDevSDF) * theDivider;

						if (normalizeTypes == Log)
							nsdfArray[j] = csdf.logit(nsdfArray[j], logAlpha);
					}
					
					double fitness;
					if (nsdfArray[j] >= targetLeftValue && nsdfArray[j] < targetRightValue)
						fitness = pow(distance(deformedPointArray[j], originalPoint), 2);
					else
						fitness = FLT_MAX;

					individualOptimalFitness[j] = fitness;
					individualOptimalNSDF[j] = nsdfArray[j];

					individualOptimalPosition[0][j] = deformedPointArray[j].x();
					individualOptimalPosition[1][j] = deformedPointArray[j].y();
					individualOptimalPosition[2][j] = deformedPointArray[j].z();
				}

			} /// OpenMP

			//cout << "calculate SDF done!" << endl;

			globalOptimalFitness = FLT_MAX;
			for (int j = 0; j < particleNumbers; j++)
			{
				if (individualOptimalFitness[j] >= globalOptimalFitness)
					continue;

				globalOptimalFitness = individualOptimalFitness[j];
				globalOptimalNSDF = individualOptimalNSDF[j];

				globalOptimalPosition[0] = individualOptimalPosition[0][j];
				globalOptimalPosition[1] = individualOptimalPosition[1][j];
				globalOptimalPosition[2] = individualOptimalPosition[2][j];
			}

			//cout << "initialization done!" << endl;			

			//////////////////////////////////////////////////////////////////////////
			// Release memory
			//////////////////////////////////////////////////////////////////////////
			#pragma omp parallel num_threads(threadsNumber)
			{
				#pragma omp for
				for (int j = 0; j < particleNumbers; j++)
				{
					TNT::Array3D<RayIntersectCell*> grid = gridArray[j];

					//delete m_grid
					for (int x = 0; x < grid.dim1(); x++)
						for (int y = 0; y < grid.dim2(); y++)
							for (int z = 0; z < grid.dim3(); z++)
								delete grid[x][y][z];

					triangleVector_t triangleVector = triangleVectorArray[j];

					for (int i = 0; i < triangleVector.size(); i++)
						delete triangleVector[i];

					triangleVector.clear();
				}
			}

			//cout << "release memory done!\n";

			//////////////////////////////////////////////////////////////////////////
			// Optimization
			//////////////////////////////////////////////////////////////////////////
			for (int cycle = 0; cycle < maxCycles; cycle++)
 			{
				if (cycle % 10 == 9)
					cout << "cycle = " << cycle+1 << " , maxCycles = " << maxCycles << endl;

				TNT::Array1D<TNT::Array3D<RayIntersectCell*>> gridArray = TNT::Array1D<TNT::Array3D<RayIntersectCell*>>(particleNumbers);
				TNT::Array1D<triangleVector_t> triangleVectorArray = TNT::Array1D<triangleVector_t>(particleNumbers);

				TNT::Array1D<double> xminArray = TNT::Array1D<double>(particleNumbers);
				TNT::Array1D<double> xmaxArray = TNT::Array1D<double>(particleNumbers);
				TNT::Array1D<double> yminArray = TNT::Array1D<double>(particleNumbers);
				TNT::Array1D<double> ymaxArray = TNT::Array1D<double>(particleNumbers);
				TNT::Array1D<double> zminArray = TNT::Array1D<double>(particleNumbers);
				TNT::Array1D<double> zmaxArray = TNT::Array1D<double>(particleNumbers);

				TNT::Array1D<double> xspreadArray = TNT::Array1D<double>(particleNumbers);
				TNT::Array1D<double> yspreadArray = TNT::Array1D<double>(particleNumbers);
				TNT::Array1D<double> zspreadArray = TNT::Array1D<double>(particleNumbers);

				for (int j = 0; j < particleNumbers; j++)
				{
					Point_3 deformedPoint(position[0][j], position[1][j], position[2][j]);

					v->point() = deformedPointArray[j] = deformedPoint; // set vertex position
					watermarkedMesh->computeNormals();
					deformedNormalArray[j] = v->normal();

					RayIntersect rayIntersect;
					rayIntersect.Init(*watermarkedMesh, gridArray[j], triangleVectorArray[j],
										xminArray[j], xmaxArray[j], yminArray[j], ymaxArray[j], zminArray[j], zmaxArray[j],
										xspreadArray[j], yspreadArray[j], zspreadArray[j], 
										gridSize, debugMode);
				}

				#pragma omp parallel num_threads(threadsNumber)
				{
					#pragma omp for
					for (int j = 0; j < particleNumbers; j++)
					{
						CalculateSDF csdf;

						sdfArray[j] = csdf.computeVolumeSDF(
							deformedPointArray[j], deformedNormalArray[j], gridSize, 
							numCones, coneSeparationDegree, coneSeparationRadian, raysInCone,
							removeSameDirectionIntersection, gaussianWeights, removeOutlier,
							gridArray[j], triangleVectorArray[j],
							xminArray[j], xmaxArray[j], yminArray[j], ymaxArray[j], zminArray[j], zmaxArray[j],
							xspreadArray[j], yspreadArray[j], zspreadArray[j],
							debugMode);

						if (sdfArray[j] == 0)
						{
							double smoothedValue = 0.0;
							double smoothedWeights = 0.0;
						
							Mesh::Halfedge_around_vertex_circulator c = v->vertex_begin();
							do {
								smoothedValue += c->opposite()->vertex()->volumeSDF();
								smoothedWeights += 1.0;
							} while (++c != v->vertex_begin());

							sdfArray[j] = safeDiv(smoothedValue, smoothedWeights);
						}

						if (normalizeTypes != None)
						{
							//const double theDivider = 1 / (maxOriginalSDF - minOriginalSDF);
							const double theDivider = 1 / (avgOriginalStdDevSDF - minOriginalStdDevSDF);

							nsdfArray[j] = (sdfArray[j] - minOriginalStdDevSDF) * theDivider;

							if (normalizeTypes == Log)
								nsdfArray[j] = csdf.logit(nsdfArray[j], logAlpha);
						}
						
						double fitness;
						if (nsdfArray[j] >= targetLeftValue && nsdfArray[j] < targetRightValue)							
							fitness = pow(distance(deformedPointArray[j], originalPoint), 2);
						else
							fitness = FLT_MAX;
						
						if (fitness >= individualOptimalFitness[j])
							continue;
						
						individualOptimalFitness[j] = fitness;
						individualOptimalNSDF[j] = nsdfArray[j];

						individualOptimalPosition[0][j] = deformedPointArray[j].x();
						individualOptimalPosition[1][j] = deformedPointArray[j].y();
						individualOptimalPosition[2][j] = deformedPointArray[j].z();
					}
					
				} /// end OpenMP


				//////////////////////////////////////////////////////////////////////////
				// Individual optimal done
				//////////////////////////////////////////////////////////////////////////
				for (int j = 0; j < particleNumbers; j++)
				{
					if (individualOptimalFitness[j] >= globalOptimalFitness)
						continue;
					
					globalOptimalFitness = individualOptimalFitness[j];
					globalOptimalNSDF = individualOptimalNSDF[j];

					globalOptimalPosition[0] = individualOptimalPosition[0][j];
					globalOptimalPosition[1] = individualOptimalPosition[1][j];
					globalOptimalPosition[2] = individualOptimalPosition[2][j];
				}

				//////////////////////////////////////////////////////////////////////////
				// Modify the velocity and position of particle
				//////////////////////////////////////////////////////////////////////////
				#pragma omp parallel num_threads(threadsNumber)
				{
					#pragma omp for
					for (int j = 0; j < particleNumbers; j++)
					{
						velocity[0][j] = constrictionFactor * (velocity[0][j] + c1 * randomNumber(0,1) * (individualOptimalPosition[0][j] - position[0][j]) + c2 * randomNumber(0,1) * (globalOptimalPosition[0] - position[0][j]));
						velocity[1][j] = constrictionFactor * (velocity[1][j] + c1 * randomNumber(0,1) * (individualOptimalPosition[1][j] - position[1][j]) + c2 * randomNumber(0,1) * (globalOptimalPosition[1] - position[1][j]));
						velocity[2][j] = constrictionFactor * (velocity[2][j] + c1 * randomNumber(0,1) * (individualOptimalPosition[2][j] - position[2][j]) + c2 * randomNumber(0,1) * (globalOptimalPosition[2] - position[2][j]));

						position[0][j] = position[0][j] + velocity[0][j];
						position[1][j] = position[1][j] + velocity[1][j];
						position[2][j] = position[2][j] + velocity[2][j];
					}
				} /// end OpenMP

				//////////////////////////////////////////////////////////////////////////
				// Check outliers
				//////////////////////////////////////////////////////////////////////////
				#pragma omp parallel num_threads(threadsNumber)
				{
					#pragma omp for
					for (int j = 0; j < particleNumbers; j++)
					{						
						Point_3 deformedPoint(position[0][j], position[1][j], position[2][j]);

						double newRadius = distance(originalPoint, deformedPoint);

						if (newRadius <= searchSpaceDist)
							continue;

						position[0][j] = position[0][j] + ((newRadius - searchSpaceDist + 2 * delta) * (originalPoint.x() - position[0][j]) / newRadius);
						position[1][j] = position[1][j] + ((newRadius - searchSpaceDist + 2 * delta) * (originalPoint.y() - position[1][j]) / newRadius);
						position[2][j] = position[2][j] + ((newRadius - searchSpaceDist + 2 * delta) * (originalPoint.z() - position[2][j]) / newRadius);

						Point_3 newdeformedPoint(position[0][j], position[1][j], position[2][j]);

						double updatedNewRadius = distance(originalPoint, newdeformedPoint);
						
						//cout << "updatedNewRadius= " << updatedNewRadius << ", avgDist= " << searchSpaceDist << endl;

						if (updatedNewRadius > searchSpaceDist)
						{
							cout << "\t#Error on re-assign points: updated new radius: " << updatedNewRadius  << " , radius: " << searchSpaceDist << endl;
							position[0][j] = originalPoint.x();
							position[1][j] = originalPoint.y();
							position[2][j] = originalPoint.z();
						}
					}

				} /// end OpenMP


				//////////////////////////////////////////////////////////////////////////
				// Release memory
				//////////////////////////////////////////////////////////////////////////
				#pragma omp parallel num_threads(threadsNumber)
				{
					#pragma omp for
					for (int j = 0; j < particleNumbers; j++)
					{
						TNT::Array3D<RayIntersectCell*> grid = gridArray[j];

						//delete m_grid
						for (int x = 0; x < grid.dim1(); x++)
							for (int y = 0; y < grid.dim2(); y++)
								for (int z = 0; z < grid.dim3(); z++)
									delete grid[x][y][z];

						triangleVector_t triangleVector = triangleVectorArray[j];

						for (int i = 0; i < triangleVector.size(); i++)
							delete triangleVector[i];

						triangleVector.clear();
					}
				}

			} /// end-for-loop(maxCycles)


			cout << "Optimization process done!\n";

			//////////////////////////////////////////////////////////////////////////
			// The vertex cannot deformed to the target SDF value, give it up
			//////////////////////////////////////////////////////////////////////////
			if (globalOptimalFitness == FLT_MAX)
			{
				cout << "    global optimal fitness is FLT_MAX" << endl;
				v->point() = originalPoint;
				continue;
			}

			//////////////////////////////////////////////////////////////////////////
			// Check whether the minOriginalStdDevSDF/minOriginalStdDevNSDF/minOriginalStdDevNSDF had changed
			// If not, deforming the vertex
			// Else, do not deforming the vertex
			//////////////////////////////////////////////////////////////////////////
			Point_3 updatedVertexPosition(globalOptimalPosition[0], globalOptimalPosition[1], globalOptimalPosition[2]);
			v->point() = updatedVertexPosition;

			watermarkedMesh->computeNormals();

			CalculateSDF csdf;
			csdf.go(
				watermarkedMesh, onVertices, debugMode, multiThreaded, threadsNum, 
				gridSize, numCones, coneSeparationDegree, raysInCone, 
				removeSameDirectionIntersection, gaussianWeights, removeOutlier, 
				normalizeTypes, logAlpha,
				smoothing, smoothingIterations);
			
			/// keep minStdDevSDF ///
			const double currentMinStdDevSDF = watermarkedMesh->getMinStdDevVolume();
			/*const double currentMinStdDevNSDF = watermarkedMesh->getMinStdDevNVolume();
			const double currentMaxStdDevNSDF = watermarkedMesh->getMaxStdDevNVolume();*/
			if (currentMinStdDevSDF != minOriginalStdDevSDF)
			{
				cout << "    currentMinStdDevSDF: " << currentMinStdDevSDF << ", minOriginalStdDevSDF: " << minOriginalStdDevSDF << endl;
				v->point() = originalPoint;
				continue;
			}
			/*
			if (currentMinStdDevNSDF != minOriginalStdDevNSDF)
			{
				cout << "    currentMinStdDevNSDF: " << currentMinStdDevNSDF << ", minOriginalStdDevNSDF: " << minOriginalStdDevNSDF << endl;
				v->point() = originalPoint;
				continue;
			}
			if (currentMaxStdDevNSDF != maxOriginalStdDevNSDF)
			{
				cout << "    currentMaxStdDevNSDF: " << currentMaxStdDevNSDF << ", maxOriginalStdDevNSDF: " << maxOriginalStdDevNSDF << endl;
				v->point() = originalPoint;
				continue;
			}
			if (currentMaxSDF != maxOriginalSDF)
			{
				cout << "    currentMaxSDF: " << currentMaxSDF << ", maxOriginalSDF: " << maxOriginalSDF << endl;
				v->point() = originalPoint;
				continue;
			}
			*/
			//////////////////////////////////////////////////////////////////////////

			/// deforming vertex position ///
			if (globalOptimalNSDF >= targetLeftValue && globalOptimalNSDF < targetRightValue)
			{
				hasShiftPointsNumber++;

				cout << "    current watermark index = " << watermarkIndex << " , hasShiftPointsNumber = " << hasShiftPointsNumber << " , shiftPointsNumber = " << shiftPointsNumber << endl;

				//////////////////////////////////////////////////////////////////////////
				// Update the deformed vertex position
				//////////////////////////////////////////////////////////////////////////
				Point_3 deformedPoint(globalOptimalPosition[0], globalOptimalPosition[1], globalOptimalPosition[2]);

				v->point() = deformedPoint;
				v->volumeNSDF() = globalOptimalNSDF;

				cout.precision(10);
				cout << "    [" << targetLeftValue << " <= " << v->volumeNSDF() << " < " << targetRightValue << ")\n";
				cout << "    original point = " << "(" << originalPoint.x() << "," << originalPoint.y() << "," << originalPoint.z() << ")\n";
				cout << "    deformed point = " << "(" << deformedPoint.x() << "," << deformedPoint.y() << "," << deformedPoint.z() << ")\n";
				cout << "    originalNSDF = " << originalNSDF << " , globalOptimalNSDF = " << globalOptimalNSDF << endl;
			}
			else
			{
				cout << "failed" << endl;
				v->point() = originalPoint;
			}
			
			cout << endl;

			//////////////////////////////////////////////////////////////////////////
			
			/// Terminate condition ///
			if (hasShiftPointsNumber == shiftPointsNumber)
				break;
			//////////////////////////////////////////////////////////////////////////

		} // end-for-loop(sizeOfSortedSquareErrors)		

		cout << "2.c hasShiftPointsNumber = " << hasShiftPointsNumber << endl;
		
		//////////////////////////////////////////////////////////////////////////
		// Recalculate SDF values of the whole mesh
		//////////////////////////////////////////////////////////////////////////
		watermarkedMesh->computeNormals();

		//CalculateSDF csdf;
		csdf.go(
			watermarkedMesh, onVertices, debugMode, multiThreaded, threadsNum, 
			gridSize, numCones, coneSeparationDegree, raysInCone, 
			removeSameDirectionIntersection, gaussianWeights, removeOutlier, 
			normalizeTypes, logAlpha,
			smoothing, smoothingIterations);

		//////////////////////////////////////////////////////////////////////////
		// Reallocate vertex index bins
		//////////////////////////////////////////////////////////////////////////
		const double minDeformedStdDevSDF = watermarkedMesh->getMinStdDevVolume();
		const double maxDeformedStdDevSDF = watermarkedMesh->getMaxStdDevVolume();
		const double avgDeformedStdDevSDF = watermarkedMesh->getAvgStdDevVolume();
		const double minDeformedStdDevNSDF = watermarkedMesh->getMinStdDevNVolume();
		const double maxDeformedStdDevNSDF = watermarkedMesh->getMaxStdDevNVolume();
		const double avgDeformedStdDevNSDF = watermarkedMesh->getAvgStdDevNVolume();

		const double deformedRange = maxDeformedStdDevNSDF - minDeformedStdDevNSDF;
		const double deformedBinInterval = deformedRange / numBins;

		for (int index = 0; index < numBins; index++)
			vertexIndexBins[index].clear();

		for (Mesh::Vertex_iterator vit = watermarkedMesh->vertices_begin(); vit != watermarkedMesh->vertices_end(); vit++)
		{
			double deformedNSDF = vit->volumeNSDF();
			if (deformedNSDF < minDeformedStdDevNSDF || deformedNSDF > maxDeformedStdDevNSDF)
				continue;

			int deformedBinIntervalIndex = (deformedNSDF - minDeformedStdDevNSDF) / deformedBinInterval;
			if (deformedBinIntervalIndex < 0 || deformedBinIntervalIndex > numBins)
				continue;

			if (deformedBinIntervalIndex == numBins)
				vertexIndexBins[deformedBinIntervalIndex-1].push_back(vit->index());
			else
				vertexIndexBins[deformedBinIntervalIndex].push_back(vit->index());
		}

		//////////////////////////////////////////////////////////////////////////
		// Check updated bin
		//////////////////////////////////////////////////////////////////////////
		cout << "\tCheck updated bin...";

		const int deformedFirstBinSize = vertexIndexBins[firstBinIndex].size();
		const int deformedSecondBinSize = vertexIndexBins[secondBinIndex].size();

		if (watermarkBit == 1)
		{
			cout << "Embed watermark bit '1' (totalBits = " << totalWatermarkBits << ")...";
			if (deformedFirstBinSize <= deformedSecondBinSize)
				cout << "success!\n";
			else
				cout << "failed!\n";
		}
		else
		{
			cout << "Embed watermark bit '0' (totalBits = " << totalWatermarkBits << ")...";
			if (deformedFirstBinSize > deformedSecondBinSize)
				cout << "success!\n";
			else
				cout << "failed!\n";
		}

		//////////////////////////////////////////////////////////////////////////
		// Export mesh
		//////////////////////////////////////////////////////////////////////////
		cout << "\tStart to export watermarked mesh model...\n";

		QString bitIndex;
		QString qBit;
		bitIndex.setNum(watermarkIndex);
		qBit.setNum(totalWatermarkBits);

		QString watermarkedMeshExportFileName = appObject->fileDir + appObject->fileSeparator + "Watermark_SDF_Relation_" + qBit + "Bits_bitIndex_" + bitIndex + "_" + appObject->fileName + appObject->fileSuffix;
		watermarkedMesh->write_obj(watermarkedMeshExportFileName);

		cout << "\tExport watermarked mesh model done!\n\n";

	} // end-for-loop(watermarkIndex)

	cout << "\tWatermark Embedding Procedure done!\n\n";

	return true;	
}


bool SDFWatermarkingAlgorithm::detect(
	AppObject* appObject, Mesh* watermarkedMesh, const bool onVertices, 
	const long seed, const int totalWatermarkBits,
	const int gridSize, const int numCones, const float coneSeparationDegree, const int raysInCone, 
	const bool gaussianWeights, const bool removeOutlier, const NormalizeTypes normalizeTypes, const int logAlpha,
	const bool smoothing, const int smoothingIterations,
	const bool debugMode, const bool showDialog)
{
	QTime timer;

	timer.start();

	const double coneSeparationRadian = coneSeparationDegree * DEG2RAD;
	const int numBins = 2 * totalWatermarkBits + 2;

	//////////////////////////////////////////////////////////////////////////
	// Basic model information
	//////////////////////////////////////////////////////////////////////////
	const int totalVertices = watermarkedMesh->size_of_vertices();
	const int totalFacets = watermarkedMesh->size_of_facets();
	const int totalEdges = watermarkedMesh->size_of_halfedges() / 2; 

	//////////////////////////////////////////////////////////////////////////
	// Find minimum SDF and maximum SDF value
	//////////////////////////////////////////////////////////////////////////
	const double minSDF = watermarkedMesh->getMinVolume();
	const double maxSDF = watermarkedMesh->getMaxVolume();
	const double avgSDF = watermarkedMesh->getAvgVolume();
	const double minNSDF = watermarkedMesh->getMinNVolume();
	const double maxNSDF = watermarkedMesh->getMaxNVolume();
	const double avgNSDF = watermarkedMesh->getAvgNVolume();

	const double minStdDevSDF = watermarkedMesh->getMinStdDevVolume();
	const double maxStdDevSDF = watermarkedMesh->getMaxStdDevVolume();
	const double avgStdDevSDF = watermarkedMesh->getAvgStdDevVolume();
	const double minStdDevNSDF = watermarkedMesh->getMinStdDevNVolume();
	const double maxStdDevNSDF = watermarkedMesh->getMaxStdDevNVolume();
	const double avgStdDevNSDF = watermarkedMesh->getAvgStdDevNVolume();

	//////////////////////////////////////////////////////////////////////////
	// Allocate all the vertices into bins
	//////////////////////////////////////////////////////////////////////////
	const double range = maxStdDevNSDF - minStdDevNSDF;
	const double binInterval = range / numBins;
	
	cout.precision(10);
	cout << "minSDF = " << minSDF << " , maxSDF = " << maxSDF << " , avgSDF = " << avgSDF << endl;
	cout << "minNSDF = " << minNSDF << " , maxNSDF = " << maxNSDF << " , avgNSDF = " << avgNSDF << endl;
	cout << "minStdDevSDF = " << minStdDevSDF << " , maxStdDevSDF = " << maxStdDevSDF << " , avgStdDevSDF = " << avgStdDevSDF << endl;
	cout << "minStdDevNSDF = " << minStdDevNSDF << " , maxStdDevNSDF = " << maxStdDevNSDF << " , avgStdDevNSDF = " << avgStdDevNSDF << endl;
	cout << "binInterval = " << binInterval;

	vector<vector<int> > vertexIndexBins(numBins);
	Mesh::Vertex_iterator vit = watermarkedMesh->vertices_begin();
	Mesh::Vertex_iterator vit_end = watermarkedMesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		double nsdf = vit->volumeNSDF();
		if (nsdf < minStdDevNSDF || nsdf > maxStdDevNSDF)
			continue;

		int binIntervalIndex = (nsdf - minStdDevNSDF) / binInterval;		
		if (binIntervalIndex > numBins || binIntervalIndex < 0)
			continue;

		if (binIntervalIndex == numBins)
			vertexIndexBins[binIntervalIndex-1].push_back(vit->index());
		else
			vertexIndexBins[binIntervalIndex].push_back(vit->index());
	}

	//////////////////////////////////////////////////////////////////////////
	// Watermark detecting procedure
	//////////////////////////////////////////////////////////////////////////
	cout << "\n\nStart to detect watermark with totalWatermarkBits = " << totalWatermarkBits << endl;

	TNT::Array1D<int> detectedWatermarkSequence = TNT::Array1D<int>(totalWatermarkBits);
	meshDetectigProcedure(watermarkedMesh, binInterval, vertexIndexBins, detectedWatermarkSequence);

	int timerDetect = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Generate pseudo random watermark sequence
	//////////////////////////////////////////////////////////////////////////
	TNT::Array1D<int> originalWatermarkSequence = TNT::Array1D<int>(totalWatermarkBits);
	
	PseudoRandomSequenceGenerator generator;
	generator.generate(seed, originalWatermarkSequence);

	//////////////////////////////////////////////////////////////////////////
	// Statistical analysis
	//////////////////////////////////////////////////////////////////////////
	cout << "\nOriginal Watermark Sequence\n";
	cout << "[";
	for (int watermarkBitIndex = 0; watermarkBitIndex < totalWatermarkBits; watermarkBitIndex++)
	{
		cout << originalWatermarkSequence[watermarkBitIndex];
		if (watermarkBitIndex != (totalWatermarkBits - 1))
			cout << ",";
	}
	cout << "]" << endl;

	cout << "\nDetected Watermark Sequence\n";
	cout << "[";
	for (int watermarkBitIndex = 0; watermarkBitIndex < totalWatermarkBits; watermarkBitIndex++)
	{
		cout << detectedWatermarkSequence[watermarkBitIndex];
		if (watermarkBitIndex != (totalWatermarkBits - 1))
			cout << ",";
	}
	cout << "]" << endl << endl;

	Statistic statistic;

	double bitErrorRate = statistic.calculateBitErrorRate(originalWatermarkSequence, detectedWatermarkSequence);
	double correlationCoefficient = statistic.calcuateCorrelationCoefficient(originalWatermarkSequence, detectedWatermarkSequence);

	cout << "BER = " << bitErrorRate << " (total totalWatermarkBits = " << totalWatermarkBits << ")\n";
	cout << "CC = " << correlationCoefficient << endl;

	//////////////////////////////////////////////////////////////////////////
	// Export log information
	//////////////////////////////////////////////////////////////////////////
	QString exportLogFileName = appObject->fileDir + appObject->fileSeparator + "WD_" + appObject->fileName + "_log.txt";

	ExporterLog exporterLog;
	exporterLog.exportRelationSDFWatermarkDetectingLog(
		exportLogFileName, 
		gridSize, numCones, coneSeparationDegree, coneSeparationRadian, raysInCone, 
		gaussianWeights, removeOutlier, normalizeTypes, logAlpha,
		minSDF, maxSDF, avgSDF, minNSDF, maxNSDF, avgNSDF,
		minStdDevSDF, maxStdDevSDF, avgStdDevSDF, minStdDevNSDF, maxStdDevNSDF, avgStdDevNSDF,
		numBins, binInterval,
		seed, totalWatermarkBits,
		totalVertices, totalFacets, totalEdges,
		vertexIndexBins,
		originalWatermarkSequence, detectedWatermarkSequence, 
		bitErrorRate, correlationCoefficient,
		timerDetect);

	return true;

	/*QTime timer;

	timer.start();

	const double coneSeparationRadian = coneSeparationDegree * DEG2RAD;
	const int numBins = 2 * totalWatermarkBits + 2;

	//////////////////////////////////////////////////////////////////////////
	// Basic model information
	//////////////////////////////////////////////////////////////////////////
	const int totalVertices = watermarkedMesh->size_of_vertices();
	const int totalFacets = watermarkedMesh->size_of_facets();
	const int totalEdges = watermarkedMesh->size_of_halfedges() / 2; 

	//////////////////////////////////////////////////////////////////////////
	// Find minimum SDF and maximum SDF value
	//////////////////////////////////////////////////////////////////////////
	const double minSDF = watermarkedMesh->getMinVolume();
	const double maxSDF = watermarkedMesh->getMaxVolume();
	const double avgSDF = watermarkedMesh->getAvgVolume();
	const double minNSDF = watermarkedMesh->getMinNVolume();
	const double maxNSDF = watermarkedMesh->getMaxNVolume();
	const double avgNSDF = watermarkedMesh->getAvgNVolume();

	const double minStdDevSDF = watermarkedMesh->getMinStdDevVolume();
	const double maxStdDevSDF = watermarkedMesh->getMaxStdDevVolume();
	const double avgStdDevSDF = watermarkedMesh->getAvgStdDevVolume();
	const double minStdDevNSDF = watermarkedMesh->getMinStdDevNVolume();
	const double maxStdDevNSDF = watermarkedMesh->getMaxStdDevNVolume();
	const double avgStdDevNSDF = watermarkedMesh->getAvgStdDevNVolume();

	//////////////////////////////////////////////////////////////////////////
	// Allocate all the vertices into bins
	//////////////////////////////////////////////////////////////////////////
	const double range = maxStdDevNSDF - minStdDevNSDF;
	const double binInterval = range / numBins;
	
	cout.precision(10);
	cout << "minSDF = " << minSDF << " , maxSDF = " << maxSDF << " , avgSDF = " << avgSDF << endl;
	cout << "minNSDF = " << minNSDF << " , maxNSDF = " << maxNSDF << " , avgNSDF = " << avgNSDF << endl;
	cout << "minStdDevSDF = " << minStdDevSDF << " , maxStdDevSDF = " << maxStdDevSDF << " , avgStdDevSDF = " << avgStdDevSDF << endl;
	cout << "minStdDevNSDF = " << minStdDevNSDF << " , maxStdDevNSDF = " << maxStdDevNSDF << " , avgStdDevNSDF = " << avgStdDevNSDF << endl;
	cout << "binInterval = " << binInterval;

	vector<vector<int> > vertexIndexBins(numBins);
	Mesh::Vertex_iterator vit = watermarkedMesh->vertices_begin();
	Mesh::Vertex_iterator vit_end = watermarkedMesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		double nsdf = vit->volumeNSDF();
		if (nsdf < minStdDevNSDF || nsdf > maxStdDevNSDF)
			continue;

		int binIntervalIndex = (nsdf - minStdDevNSDF) / binInterval;		
		if (binIntervalIndex > numBins || binIntervalIndex < 0)
			continue;

		if (binIntervalIndex == numBins)
			vertexIndexBins[binIntervalIndex-1].push_back(vit->index());
		else
			vertexIndexBins[binIntervalIndex].push_back(vit->index());
	}

	//////////////////////////////////////////////////////////////////////////
	// Watermark detecting procedure
	//////////////////////////////////////////////////////////////////////////
	cout << "\n\nStart to detect watermark with totalWatermarkBits = " << totalWatermarkBits << endl;

	TNT::Array1D<int> detectedWatermarkSequence = TNT::Array1D<int>(totalWatermarkBits);
	meshDetectigProcedure(watermarkedMesh, binInterval, vertexIndexBins, detectedWatermarkSequence);

	int timerDetect = timer.elapsed();

	//////////////////////////////////////////////////////////////////////////
	// Generate pseudo random watermark sequence
	//////////////////////////////////////////////////////////////////////////
	TNT::Array1D<int> originalWatermarkSequence = TNT::Array1D<int>(totalWatermarkBits);
	
	PseudoRandomSequenceGenerator generator;
	generator.generate(seed, originalWatermarkSequence);

	//////////////////////////////////////////////////////////////////////////
	// Statistical analysis
	//////////////////////////////////////////////////////////////////////////
	cout << "\nOriginal Watermark Sequence\n";
	cout << "[";
	for (int watermarkBitIndex = 0; watermarkBitIndex < totalWatermarkBits; watermarkBitIndex++)
	{
		cout << originalWatermarkSequence[watermarkBitIndex];
		if (watermarkBitIndex != (totalWatermarkBits - 1))
			cout << ",";
	}
	cout << "]" << endl;

	cout << "\nDetected Watermark Sequence\n";
	cout << "[";
	for (int watermarkBitIndex = 0; watermarkBitIndex < totalWatermarkBits; watermarkBitIndex++)
	{
		cout << detectedWatermarkSequence[watermarkBitIndex];
		if (watermarkBitIndex != (totalWatermarkBits - 1))
			cout << ",";
	}
	cout << "]" << endl << endl;

	Statistic statistic;

	double bitErrorRate = statistic.calculateBitErrorRate(originalWatermarkSequence, detectedWatermarkSequence);
	double correlationCoefficient = statistic.calcuateCorrelationCoefficient(originalWatermarkSequence, detectedWatermarkSequence);

	cout << "BER = " << bitErrorRate << " (total totalWatermarkBits = " << totalWatermarkBits << ")\n";
	cout << "CC = " << correlationCoefficient << endl;

	//////////////////////////////////////////////////////////////////////////
	// Export log information
	//////////////////////////////////////////////////////////////////////////
	QString exportLogFileName = appObject->fileDir + appObject->fileSeparator + "WD_" + appObject->fileName + "_log.txt";

	ExporterLog exporterLog;
	exporterLog.exportRelationSDFWatermarkDetectingLog(
		exportLogFileName, 
		gridSize, numCones, coneSeparationDegree, coneSeparationRadian, raysInCone, 
		gaussianWeights, removeOutlier, normalizeTypes, logAlpha,
		minSDF, maxSDF, avgSDF, minNSDF, maxNSDF, avgNSDF,
		minStdDevSDF, maxStdDevSDF, avgStdDevSDF, minStdDevNSDF, maxStdDevNSDF, avgStdDevNSDF,
		numBins, binInterval,
		seed, totalWatermarkBits,
		totalVertices, totalFacets, totalEdges,
		vertexIndexBins,
		originalWatermarkSequence, detectedWatermarkSequence, 
		bitErrorRate, correlationCoefficient,
		timerDetect);

	return true;*/
}


void SDFWatermarkingAlgorithm::meshDetectigProcedure(
	Mesh* watermarkedMesh, const int binInterval,
	vector<vector<int> > vertexIndexBins, TNT::Array1D<int>& detectedWatermarkSequence)
{
	const int numBins = vertexIndexBins.size();

	int watermarkIndex = 0;
	for (int binIndex = 1; binIndex < (numBins-1); binIndex+=2)
	{
		if((binIndex + 1) == (numBins - 1))
			continue;

		int firstBinIndex = binIndex;
		int secondBinIndex = firstBinIndex + 1;

		int firstBinSize = vertexIndexBins[firstBinIndex].size();
		int secondBinSize = vertexIndexBins[secondBinIndex].size();

		if (firstBinSize <= secondBinSize)
			detectedWatermarkSequence[watermarkIndex] = 1;
		else
			detectedWatermarkSequence[watermarkIndex] = 0;

		watermarkIndex++;
	}
}