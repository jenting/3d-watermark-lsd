#ifndef __APP_MANAGER_H_
#define __APP_MANAGER_H_

#include "stdafx.h"


#include <q3filedialog.h>
#include <QFileDialog>
#include <qfileinfo.h>
#include <qmessagebox.h>
#include <qinputdialog.h>
#include <qstringlist.h>
#include <qfile.h>
#include <qbuffer.h>
#include <q3textstream.h>
#include <qcombobox.h>
#include <q3listbox.h>
#include <qspinbox.h>
#include <qcheckbox.h>
#include <qlayout.h>
#include <qlabel.h>
#include <qpushbutton.h>
#include <qslider.h>
#include <q3hbox.h>
//Added by qt3to4:


#include <Q3HBoxLayout>
#include <QDir>
#include <QFileInfo>
#include <QFileInfoList>


#include "AppManager.h"
#include "AppParams.h"
#include "Attacks.h"
#include "CalculateSDF.h"
#include "CalculateCurvature.h"
#include "ColorSchemeManager.h"
#include "DrawHistogram.h"
#include "Exporter.h"
#include "FeaturePainters.h"
#include "FileBrowser.h"
#include "GreenCoordinate.h"
#include "Importer.h"
#include "MainWindow.h"
#include "MyMatrix.h"
#include "MeanValueCoordinate.h"
#include "Metric.h"
#include "mesh/Parser.h"
#include "PartitionMesh.h"
#include "Watermark.h"
#include "WorldManager.h"


class FacetPainter;
class DebugPainter;
class OverlayPainter;
class SdfOverlayPainter;

class AppParams;
class WorldObject;
class QTextStream;


struct ApplicationParameters
{
	bool debugMode;
	bool multithreaded;
	int threadsNum;
	QString dirRunner;
	QString fileName1;
	QString fileName2;
	bool normalsReversed;

	/// sdf
	int lsd_gridSize;
	int lsd_numberOfCones;
	int lsd_coneSeparation;
	int lsd_raysInCone;
	bool lsd_removeSameDirectionIntersection;
	bool lsd_gaussianWeights;
	bool lsd_removeOutlier;
	NormalizeTypes lsd_normalize;
	int lsd_log_alpha;
	bool lsd_smoothing;
	int  lsd_smoothing_iterations;

	/// partition
    int partition_candidates;
    double partition_threshold;

	/// watermarking
	long watermark_seed;
	int watermark_bits;
	int robustness_factor;
	SearchSpace watermark_searchSpace;
	int watermark_maxCycles;
	int watermark_particleNumbers;

	ApplicationParameters() :
		debugMode(false),
		multithreaded(true),
		threadsNum(8),
		dirRunner(""),
		fileName1(""),
		fileName2(""),
		normalsReversed(false),
		lsd_gridSize(40), 
		lsd_numberOfCones(4),
		lsd_coneSeparation(10),
		lsd_raysInCone(8),
		lsd_removeSameDirectionIntersection(false),
		lsd_gaussianWeights(false),
		lsd_removeOutlier(false),
		lsd_normalize(None),
		lsd_log_alpha(4),
		lsd_smoothing(false),
		lsd_smoothing_iterations(1),
		watermark_seed(123456789),
		watermark_bits(16),
		robustness_factor(4),
		watermark_searchSpace(MIN_DIST),
		watermark_maxCycles(10),
		watermark_particleNumbers(50)
	{}
};

typedef QHash<QString, AppObject*> appObjectsMap_t;


#define FEATURES_FILE_SDF_VERTICES	0x001
#define FEATURES_FILE_SDF_FACETS	0x002

class WorldManager;
class MainWindow;
class PartsSnapViewer;



class FaceColorizer
{
public:
	virtual QColor getColor(Mesh::Facet_handle f) = 0;
};


struct IndexRecord
{
	IndexRecord(int _start = -1, int _count = 0) :start(_start), count(_count) {}
	int start, count;
};


typedef QHash<QString, IndexRecord> TModelIndex;



/**
 *	This is the application manager, it keeps a list of the open meshes
 */
class AppManager : public QObject
{
	Q_OBJECT

public:
	AppManager(MainWindow* parent = NULL, const char* name = NULL);
	~AppManager();
	/** after all objects are constructed and connections are made call this function to initialize everything */
	void Init(WorldManager* worldManager);

	/// currently loaded objects
	appObjectsMap_t m_objects;

//	TFilesAppParams m_filesAppParams;

	/** current selected object id */
	QString m_selectedObject;

	/** dialog for editing application parameters */
	AppParams* m_appParamsDialog;

	/** application parameters */
	ApplicationParameters m_appParameters;

public:

	/// link to world manager
	WorldManager* m_worldManager;

	/** generate a unique id for an object from the filename */
	QString createObjectId(const QString& fileName);

	AppObject* getObject(const QString& id);

	/** calculate initial features on mesh */
	static void doMeshPrecalc(Mesh* mesh);

	/** load features from file for this app object */
	bool loadCage(const QString& fileName, AppObject* appObject, TStubData* part = NULL, DrawMethod dm = DRAW_PART);

	bool loadMeshFeatures(const QString& featuresFileName, AppObject* appObject);

	/** embed/detect watermark */
	bool embedWatermark(AppObject* appObject, const bool showDialog = true);
	bool detectWatermark(AppObject* appObject, const bool showDialog = true);

	/** calculate volume */	
	bool calculateTotalVolume(AppObject* appObject);
	bool calculateShapeDiameterFunction(AppObject* appObject, const bool onVertices, const bool showDialog = true);

	/** calculate curvature for all vertices in mesh */
	bool calculateCurvature(AppObject* appObject);

	/** calcuate Laplace-Beltrami operator */
	// bool laplaceBeltrami(AppObject* appObject);

	bool calculateProjectionDifference(AppObject* appObject);

	/** calculate deformation for different coordinates method */
	bool generateCageMesh(AppObject* appObject);
	bool meanValueCoordinates(AppObject* appObject);
	bool harmonicCoordinates(AppObject* appObject);
	bool greenCoordinates(AppObject* appObject);

	/** partition mesh according to the facets/vertices volume */
    bool partitionFacetsMesh(AppObject* appObject);
	bool partitionVerticesMesh(AppObject* appObject);

	/** mesh attacks */
	static void allAttacks(AppObject* appObject);
	void noiseAddition(AppObject* appObject, const bool preserveBoundaries, const TNT::Array1D<double>& noiseIntensityArray);
	void quantization(AppObject* appObject, const TNT::Array1D<int>& quantizationBitArray);
	void reordering(AppObject* appObject);
	void similarityTransform(AppObject* appObject);
	void smoothing(AppObject* appObject, const double deformFactor, const bool preserveBoundaries, const TNT::Array1D<int>& iterationList);
	void translation(AppObject* appObject);
	void rotation(AppObject* appObject);
	void uniformScaling(AppObject* appObject);
	void nonuniformScaling(AppObject* appObject);

	void simplification(AppObject* appObject, SimplificationType simpType, const TNT::Array1D<double>& simplificationRatioArray);
	void subdivision(AppObject* appObject, SubdivisionType subdivisionType);
	
	/** dir runner */
	void dirRunnerAvgDDegree(AppObject* appObject, const QString& fileDir);
	void dirRunnerMaxDDegree(AppObject* appObject, const QString& fileDir);
	double calculateDeformationDistance(Mesh* originalMesh, Mesh* deformedMesh);
	void calculateMaxDeformationDistance(AppObject* appObject, Mesh* sourceMesh, Mesh* targetMesh);

	/** set object with given ID as selected */
	void selectObject(const QString& id);

	void reverseObjectNormals(AppObject* appObject);

	const ApplicationParameters *getAppParams() { return &m_appParameters; }
	void saveApplicationParameters();
	void loadApplicationParameters();

	QStringList getObjectNames(int& selectedIndex);

	void createProjectionMesh(AppObject* appObject);
	void createProjectionPart(AppObject* appObject, const int partIndex);

public:
	QString selectFileName();
	QStringList selectFileNames();

	/** prompts the user to load a specified mesh */
	bool loadObject();
	bool loadObjects();

	/** 
	  * loads a mesh from a given filename
	  * @param err error tolerance. 1.0 - tolerate anything. 0.3 tolerate less
	  */
	bool loadObject(const QString& fileName, double err, TStubData* part = NULL, DrawMethod dm = DRAW_PART);

	void copyObject(AppObject *from, TStubData* newPart, DrawMethod dm);

	void removeAllObjects();

	AppObject* getSelectedObject();

signals:

	void showParts();

	/** one or more of the objects have changed and probably require redraw */
	void ObjectsChanged();

	/** a new cage object was loaded into the application */
	void CageAdded(AppObject *appObject, TStubData* onlyPart, DrawMethod dm);

	/** a new object was loaded into the application */
	void ObjectAdded(AppObject *appObject, TStubData* onlyPart, DrawMethod dm);

	/** an object was removed from the application */
	void ObjectRemoved(const QString& id);

	/** an object was changed and rendering needs to be redone */
	void ObjectChanged(const QString& id, bool geometryChanged);

	void FacetPainterChanged(int);

	/** send message to status bar */
	void message(const QString& message);

	/** send message to status bar for several seconds */
	void message(const QString& message, int ms);

public slots:

	void openFile(QString filename);

	void saveSnapshot(const QString &fileName);

	void appParamsOk();
	void appParamsApply();

	void OnArrowKeyPressed(int key, Qt::ButtonState s);

	void removeObject(QString id);

	/** load/save an object */
	void OnLoadObjects();
	void OnLoadObjectCage();
	void OnSaveObject();

	/** load/save a file with object features */
	void OnLoadObjectFeatures();
	void OnSaveObjectFeatures();
	void OnDiffObjectFeatures();
	
	/** remove all loaded objects */
	void OnRemoveObjects();

	/** remove selected object */
	void OnRemoveSelectedObject();

	/**	open application parameters window */
	void OnAppParams();

	/** show information about current object */
	void OnSelectedObjectInfo();

	/** user selected a new object */
	void OnObjectSelected(const QString& id);

	/** user selected a specific facet in a specific object */
	void OnFacetSelected(const QString& id, const int facetIndex);

	/** user selected a specific vertex in a specific object */
	void OnVertexSelected(const QString& id, const int vertexIndex);
	
	//////////////////////////////////////////////////////////////////////////
	// Calculate Volume
	//////////////////////////////////////////////////////////////////////////

	/** calculate total volume on selected mesh */
	void OnCalculateTotalVolume();
	/** calculate shape diameter function on selected mesh's vertices/facets */
	void OnCalculateLSDVertices();
	void OnCalculateLSDFacets();
	/** calculates SDF on all loaded mesh's facets */
	void OnCalculateLSDVerticesForAll();
	void OnCalculateLSDFacetsForAll();
	

	//////////////////////////////////////////////////////////////////////////
	// Segmentation
	//////////////////////////////////////////////////////////////////////////
	/** use GMM algorithm to partition on selected mesh based on facets volume */
	void OnPartitionFacetsMesh();
	/** use GMM algorithm to partition on selected mesh based on vertices volume */
	void OnPartitionVerticesMesh();

	//////////////////////////////////////////////////////////////////////////
	// Draw Histogram
	//////////////////////////////////////////////////////////////////////////
	/** draw Facet/Vertex LSD/NLSD value histogram on selected mesh */
	void OnDrawFacetLSDValue();
	void OnDrawFacetNLSDValue();
	void OnDrawVertexLSDValue();
	void OnDrawVertexNLSDValue();

	/** draw Facet/Vertex LSD/NLSD bins histogram on selected mesh */
	void OnDrawFacetLSDBins();
	void OnDrawFacetNLSDBins();
	void OnDrawVertexLSDBins();
	void OnDrawVertexNLSDBins();
	
	/** draw Facet/Vertex LSD/NLSD probability histogram on selected mesh */
	void OnDrawFacetLSDProbability();
	void OnDrawFacetNLSDProbability();
	void OnDrawVertexLSDProbability();
	void OnDrawVertexNLSDProbability();

	//////////////////////////////////////////////////////////////////////////
	// Watermark
	//////////////////////////////////////////////////////////////////////////
	void OnEmbedWatermark();
	void OnDetectWatermark();

	//////////////////////////////////////////////////////////////////////////
	// Pose
	//////////////////////////////////////////////////////////////////////////
	void OnPoseTransfer();
	void OnCalculateDeformationDistance();

	//////////////////////////////////////////////////////////////////////////
	// Metric
	//////////////////////////////////////////////////////////////////////////
	/** calculate signal-to-noise ratio(SNR) */
	void OnCalculateSNR();
	/** calculate peak signal-to-noise ratio(PSNR) */
	void OnCalculatePSNR();
	/** calculate Maximum Root Mean Square Error (MRMS) */
	void OnCalculateMRMS();
	/** calculate Hausdorff distance */
	void OnCalculateHausdorffDistance();


	//////////////////////////////////////////////////////////////////////////
	// Deformation
	//////////////////////////////////////////////////////////////////////////
	/** generate cage mesh for cage-based deformation used */
	void OnGenerateCageMesh();
	/** mean value coordinates deformation */
	void OnMeanValueCoordinates();
	/** harmoic coordinates deformation */
	void OnHarmonicCoordinates();
	/** green coordinates deformatoin */
	void OnGreenCoordinates();

	//////////////////////////////////////////////////////////////////////////
	// Others
	//////////////////////////////////////////////////////////////////////////
	/** reverse normals for selected object */
	void OnReverseNormals();
	/** calculate curvature for all vertices in a mesh */
	void OnCalculateCurvature();

	//////////////////////////////////////////////////////////////////////////
	// Attacks
	//////////////////////////////////////////////////////////////////////////
	/** attack on selected mesh */
	void OnAllAttacks();

	/// geometry ///
	void OnReordering();
	void OnNoiseAddition();
	void OnQuantization();
	void OnSmoothing();
	void OnSimilarityTransform();
	void OnTranslation();
	void OnRotation();
	void OnUniformScaling();
	void OnNonUniformScaling();

	/// topology ///
	void OnSimplificationLindstromturk();
	void OnSimplificationMidpoint();
	void OnSubdivisionCatmullClark();
	void OnSubdivisionDoosabin();
	void OnSubdivisionLoop();
	void OnSubdivisionMidpoint();
	void OnSubdivisionSqrt3();

	//////////////////////////////////////////////////////////////////////////
	// Directory Runner
	//////////////////////////////////////////////////////////////////////////
	/// calculate SDF///
	void OnDirRunnerCalculateLSDVertices();
	void OnDirRunnerCalculateLSDFacets();

	/// watermark ///
	void OnDirRunnerDetectWatermark();

	/// all attacks ///
	void OnDirRunnerAllAttacks();

	/// geometry ///
	void OnDirRunnerReordering();
	void OnDirRunnerNoiseAttack();
	void OnDirRunnerQuantization();
	void OnDirRunnerSmoothing();
	void OnDirRunnerSimilarityTransform();

	/// topology ///
	void OnDirRunnerSimplificationLindstromturk();
	void OnDirRunnerSimplificationMidpoint();
	void OnDirRunnerSubdivisionCatmullClark();
	void OnDirRunnerSubdivisionDoosabin();
	void OnDirRunnerSubdivisionLoop();
	void OnDirRunnerSubdivisionMidpoint();
	void OnDirRunnerSubdivisionSqrt3();

	/// deformation ///
	void OnDirRunnerCalculateDeformationDistance();
	void OnDirRunnerAvgDDegree();
	void OnDirRunnerMaxDDegree();

	/// animation ///
	void OnPlayAnimation();

	/// user approves part structure (from SDF projection) for current part and model
	/// Find correspondence between two objects
	/// Show color table for correspondence results

};

#endif
