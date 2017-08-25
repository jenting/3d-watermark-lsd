/** 
 * @author: Ren-Ting Hsiao
 * @version: 1.1.0
 */ 

// sdf.cpp : Defines the entry point for the application. //

#include "stdafx.h"

#include <QApplication>
#include <QStatusBar>

#include "SdfViewer.h"
#include "MainWindow.h"

#include "AppManager.h"
#include "WorldManager.h"

#define MAX_LOADSTRING 100

class FacetPainter;

bool setupConnections(QApplication& application, AppManager* appManager, WorldManager* worldManager, MainWindow* smf);

QMainWindow *g_main = NULL;
MainWindow *g_sdfmain = NULL;
AppManager *g_appManager = NULL;

LONG WINAPI myExceptionFilter(EXCEPTION_POINTERS* ExceptionInfo)
{
	return EXCEPTION_EXECUTE_HANDLER;
}

void my_invalid_parameter(const wchar_t * expression, const wchar_t * function, const wchar_t * file, unsigned int line, uintptr_t pReserved)
{
	throw QString("invalid");
}

int main(int argc, char** argv)
{
//	SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOOPENFILEERRORBOX);
//	SetUnhandledExceptionFilter(myExceptionFilter);
//	_set_error_mode(_OUT_TO_STDERR);
//	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
//	_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
//	_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
//	_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);

	_set_invalid_parameter_handler(my_invalid_parameter);

	// Read command lines arguments.
	QApplication a(argc,argv);
	QCoreApplication *p = QCoreApplication::instance();

	MainWindow smf(NULL);
	g_main = &smf;
	g_sdfmain = &smf;

//	smf.ui.CorrespondenceToolbar->hide();

	// initialize log object
	ProgSettings ps;
	QString logFile = ps.readEntry("log_file", "e:/three_dimension/watermarking/tmp/log/log.txt"); //E:\three_dimension\watermarking\tmp\log
	int logLevel = ps.readNumEntry("log_level", 8);
	GlobalLog::init(&smf, logFile, logLevel);
	ps.writeEntry("log_file", logFile);
	ps.writeEntry("log_level", logLevel);

	SDFLOG10("===================== Mesh Watermark Starting ========================");
    cout << "===================== Mesh Watermark Starting =========================" << endl;

	//load two main managers
	AppManager* appManager = new AppManager(&smf, "appManager");
	g_appManager = appManager;
	WorldManager* worldManager = new WorldManager(smf.ui.world, &smf, "worldManager");
    //worldManager->selectFacetPainter(1); // default facet painter : SDF Facet Painter

	if (!setupConnections(a, appManager, worldManager, &smf))
	{
		QMessageBox::warning(g_main, "Error on initialization", "Could not setup slot/signal connections");
		return -1;
	}
	//////////////////////////////////////////////////////////////////////////

	//initialize App Manager
	appManager->Init(worldManager);

	//////////////////////////////////////////////////////////////////////////

	smf.show();

	// Set the viewer as the application main widget.
	//application.setMainWidget(smf);

	// Run main loop.
	return a.exec();
}

bool setupConnections(QApplication& application, AppManager* appManager, WorldManager* worldManager, MainWindow* smf)
{
	/** QObject::connect(sender, SIGNAL(signal), receiver, SLOT(slot)) */

	bool success = true;

//	success &= GlobalLog::getLog()->connect((const QObject*)smf->ui.showLogAction,	SIGNAL(activated()), SLOT(ShowLog()));

	/** file */
	success &= appManager->connect((const QObject*)smf->ui.loadMeshAction, SIGNAL(activated()), SLOT(OnLoadObjects()));
	success &= appManager->connect((const QObject*)smf->ui.loadMeshFeaturesAction, SIGNAL(activated()), SLOT(OnLoadObjectFeatures()));
	success &= appManager->connect((const QObject*)smf->ui.loadCageAction, SIGNAL(activated()), SLOT(OnLoadObjectCage()));
	
	success &= appManager->connect((const QObject*)smf->ui.saveMeshAction, SIGNAL(activated()), SLOT(OnSaveObject()));
	success &= appManager->connect((const QObject*)smf->ui.saveMeshFeaturesAction, SIGNAL(activated()), SLOT(OnSaveObjectFeatures()));
	success &= QObject::connect(smf->ui.snapAction, SIGNAL(activated()), smf->ui.world, SLOT(saveManualSnap()));

	success &= appManager->connect((const QObject*)smf->ui.removeSelectedModelAction,	SIGNAL(activated()), SLOT(OnRemoveSelectedObject()));
	success &= appManager->connect((const QObject*)smf->ui.clearModelsAction,	SIGNAL(activated()), SLOT(OnRemoveObjects()));
	success &= appManager->connect((const QObject*)smf->ui.appParametersAction, SIGNAL(activated()), SLOT(OnAppParams()));

	/** tools */
	/// calculate volume
	success &= appManager->connect((const QObject*)smf->ui.calculateLSDFacetsAction, SIGNAL(activated()), SLOT(OnCalculateLSDFacets()));
	success &= appManager->connect((const QObject*)smf->ui.calculateLSDVerticesAction, SIGNAL(activated()), SLOT(OnCalculateLSDVertices()));

	/// histogram
	success &= appManager->connect((const QObject*)smf->ui.facetLSDValueAction, SIGNAL(activated()), SLOT(OnDrawFacetLSDValue()));
	success &= appManager->connect((const QObject*)smf->ui.facetNLSDValueAction, SIGNAL(activated()), SLOT(OnDrawFacetNLSDValue()));
	success &= appManager->connect((const QObject*)smf->ui.vertexLSDValueAction, SIGNAL(activated()), SLOT(OnDrawVertexLSDValue()));
	success &= appManager->connect((const QObject*)smf->ui.vertexNLSDValueAction, SIGNAL(activated()), SLOT(OnDrawVertexNLSDValue()));

	success &= appManager->connect((const QObject*)smf->ui.facetLSDBinsAction, SIGNAL(activated()), SLOT(OnDrawFacetLSDBins()));
	success &= appManager->connect((const QObject*)smf->ui.facetNLSDBinsAction, SIGNAL(activated()), SLOT(OnDrawFacetNLSDBins()));
	success &= appManager->connect((const QObject*)smf->ui.vertexLSDBinsAction, SIGNAL(activated()), SLOT(OnDrawVertexLSDBins()));
	success &= appManager->connect((const QObject*)smf->ui.vertexNLSDBinsAction, SIGNAL(activated()), SLOT(OnDrawVertexNLSDBins()));

	success &= appManager->connect((const QObject*)smf->ui.facetLSDProbabilityAction, SIGNAL(activated()), SLOT(OnDrawFacetLSDProbability()));
	success &= appManager->connect((const QObject*)smf->ui.facetNLSDProbabilityAction, SIGNAL(activated()), SLOT(OnDrawFacetNLSDProbability()));
	success &= appManager->connect((const QObject*)smf->ui.vertexLSDProbabilityAction, SIGNAL(activated()), SLOT(OnDrawVertexLSDProbability()));
	success &= appManager->connect((const QObject*)smf->ui.vertexNLSDProbabilityAction, SIGNAL(activated()), SLOT(OnDrawVertexNLSDProbability()));

	/// watermark
	success &= appManager->connect((const QObject*)smf->ui.embedWatermarkLSDAction, SIGNAL(activated()), SLOT(OnEmbedWatermark()));
	success &= appManager->connect((const QObject*)smf->ui.detectWatermarkLSDAction, SIGNAL(activated()), SLOT(OnDetectWatermark()));

	/// pose transfer
	//success &= appManager->connect((const QObject*)smf->ui.poseTransferAction, SIGNAL(activated()), SLOT(OnPoseTransfer()));
	success &= appManager->connect((const QObject*)smf->ui.calculateDeformationDistanceAction, SIGNAL(activated()), SLOT(OnCalculateDeformationDistance()));

	/// metric
	success &= appManager->connect((const QObject*)smf->ui.signalToNoiseRatioAction, SIGNAL(activated()), SLOT(OnCalculateSNR()));
	success &= appManager->connect((const QObject*)smf->ui.peakSignalToNoiseRatioAction, SIGNAL(activated()), SLOT(OnCalculatePSNR()));
	success &= appManager->connect((const QObject*)smf->ui.maximumRootMeanSquareAction, SIGNAL(activated()), SLOT(OnCalculateMRMS()));
	success &= appManager->connect((const QObject*)smf->ui.hausdorffDistanceAction, SIGNAL(activated()), SLOT(OnCalculateHausdorffDistance()));

	/// deformation
	success &= appManager->connect((const QObject*)smf->ui.generateCageMeshAction, SIGNAL(activated()), SLOT(OnGenerateCageMesh()));
	success &= appManager->connect((const QObject*)smf->ui.mvcAction, SIGNAL(activated()), SLOT(OnMeanValueCoordinates()));
	success &= appManager->connect((const QObject*)smf->ui.hcAction, SIGNAL(activated()), SLOT(OnHarmonicCoordinates()));
	success &= appManager->connect((const QObject*)smf->ui.gcAction, SIGNAL(activated()), SLOT(OnGreenCoordinates()));

	/// others
	success &= appManager->connect((const QObject*)smf->ui.reverseNormalsAction, SIGNAL(activated()), SLOT(OnReverseNormals()));

	/// attacks
	success &= appManager->connect((const QObject*)smf->ui.allAttacksAction, SIGNAL(activated()), SLOT(OnAllAttacks()));

	success &= appManager->connect((const QObject*)smf->ui.noiseAction, SIGNAL(activated()), SLOT(OnNoiseAddition()));
	success &= appManager->connect((const QObject*)smf->ui.quantizationAction, SIGNAL(activated()), SLOT(OnQuantization()));
	success &= appManager->connect((const QObject*)smf->ui.reorderingAction, SIGNAL(activated()), SLOT(OnReordering()));
	success &= appManager->connect((const QObject*)smf->ui.similarityTransformAction, SIGNAL(activated()), SLOT(OnSimilarityTransform()));
	success &= appManager->connect((const QObject*)smf->ui.smoothingAction, SIGNAL(activated()), SLOT(OnSmoothing()));
	success &= appManager->connect((const QObject*)smf->ui.translationAction, SIGNAL(activated()), SLOT(OnTranslation()));
	success &= appManager->connect((const QObject*)smf->ui.rotationAction, SIGNAL(activated()), SLOT(OnRotation()));;
	success &= appManager->connect((const QObject*)smf->ui.uniformScalingAction, SIGNAL(activated()), SLOT(OnUniformScaling()));
	success &= appManager->connect((const QObject*)smf->ui.nonUniformScalingAction, SIGNAL(activated()), SLOT(OnNonUniformScaling()));

	success &= appManager->connect((const QObject*)smf->ui.simplificationLindstromturkAction, SIGNAL(activated()), SLOT(OnSimplificationLindstromturk()));
	success &= appManager->connect((const QObject*)smf->ui.simplificationMidpointAction, SIGNAL(activated()), SLOT(OnSimplificationMidpoint()));
	success &= appManager->connect((const QObject*)smf->ui.subdivisionCatmullClarkAction, SIGNAL(activated()), SLOT(OnSubdivisionCatmullClark()));
	success &= appManager->connect((const QObject*)smf->ui.subdivisionDoosabinAction, SIGNAL(activated()), SLOT(OnSubdivisionDoosabin()));
	success &= appManager->connect((const QObject*)smf->ui.subdivisionLoopAction, SIGNAL(activated()), SLOT(OnSubdivisionLoop()));
	success &= appManager->connect((const QObject*)smf->ui.subdivisionSqrt3Action, SIGNAL(activated()), SLOT(OnSubdivisionSqrt3()));
	success &= appManager->connect((const QObject*)smf->ui.subdivisionMidpointAction, SIGNAL(activated()), SLOT(OnSubdivisionMidpoint()));

	/** dir runner */
	// Calculate LSD //
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerCalculateLSDVerticesAction, SIGNAL(activated()), SLOT(OnDirRunnerCalculateLSDVertices()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerCalculateLSDFacetsAction, SIGNAL(activated()), SLOT(OnDirRunnerCalculateLSDFacets()));

	// Watermark //
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerDetectWatermarkLSDAction, SIGNAL(activated()), SLOT(OnDirRunnerDetectWatermark()));

	// attacks //
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerAllAttacksAction, SIGNAL(activated()), SLOT(OnDirRunnerAllAttacks()));

	success &= appManager->connect((const QObject*)smf->ui.dirRunnerNoiseAction, SIGNAL(activated()), SLOT(OnDirRunnerNoiseAttack()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerQuantizationAction, SIGNAL(activated()), SLOT(OnDirRunnerQuantization()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerReorderingAction, SIGNAL(activated()), SLOT(OnDirRunnerReordering()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSimilarityTransformAction, SIGNAL(activated()), SLOT(OnDirRunnerSimilarityTransform()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSmoothingAction, SIGNAL(activated()), SLOT(OnDirRunnerSmoothing()));

	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSimplificationLindstromturkAction, SIGNAL(activated()), SLOT(OnDirRunnerSimplificationLindstromturk()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSimplificationMidpointAction, SIGNAL(activated()), SLOT(OnDirRunnerSimplificationMidpoint()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSubdivisionCatmullClarkAction, SIGNAL(activated()), SLOT(OnDirRunnerSubdivisionCatmullClark()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSubdivisionDoosabinAction, SIGNAL(activated()), SLOT(OnDirRunnerSubdivisionDoosabin()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSubdivisionLoopAction, SIGNAL(activated()), SLOT(OnDirRunnerSubdivisionLoop()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSubdivisionSqrt3Action, SIGNAL(activated()), SLOT(OnDirRunnerSubdivisionMidpoint()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerSubdivisionMidpointAction, SIGNAL(activated()), SLOT(OnDirRunnerSubdivisionSqrt3()));

	// deformation analysis //
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerCalculateDeformationDistanceAction, SIGNAL(activated()), SLOT(OnDirRunnerCalculateDeformationDistance()));
	success &= appManager->connect((const QObject*)smf->ui.dirRunnerDeformationDegreeAction, SIGNAL(activated()), SLOT(OnDirRunnerMaxDDegree()));

	success &= smf->statusBar()->connect(appManager, SIGNAL(message(const QString&)), SLOT(message(const QString&)));
	success &= smf->statusBar()->connect(appManager, SIGNAL(message(const QString&, int)), 	SLOT(message(const QString&, int)));

	success &= appManager->connect(worldManager, SIGNAL(ObjectSelected(const QString&)), SLOT(OnObjectSelected(const QString&)));
	success &= appManager->connect(worldManager, SIGNAL(FacetSelected(const QString&, const int)), SLOT(OnFacetSelected(const QString&, const int)));
	success &= appManager->connect(worldManager, SIGNAL(VertexSelected(const QString&, const int)), SLOT(OnVertexSelected(const QString&, const int)));
	success &= appManager->connect(worldManager, SIGNAL(doRemoveObject(const QString&)), SLOT(removeObject(const QString&)));


	success &= worldManager->connect(appManager, SIGNAL(ObjectsChanged()), SLOT(invalidate()));
	success &= worldManager->connect(appManager, SIGNAL(ObjectChanged(const QString&, bool)), SLOT(OnChangeObject(const QString&, bool)));
	success &= worldManager->connect(appManager, SIGNAL(FacetPainterChanged(int)), SLOT(OnSelectFacetPainter(int)));

	success &= worldManager->connect(appManager, SIGNAL(CageAdded(AppObject*, TStubData*, DrawMethod)), SLOT(OnAddCage(AppObject*, TStubData*, DrawMethod)));

	success &= worldManager->connect(appManager, SIGNAL(ObjectAdded(AppObject*, TStubData*, DrawMethod)), SLOT(OnAddObject(AppObject*, TStubData*, DrawMethod)));

	success &= worldManager->connect(appManager, SIGNAL(ObjectRemoved(const QString&)),	SLOT(OnRemoveObject(const QString&)));

	success &= worldManager->connect(smf, SIGNAL(functionKeyPressed(int, Qt::ButtonState)), SLOT(OnFunctionKeyPressed(int, Qt::ButtonState)));

	success &= appManager->connect(smf, SIGNAL(arrowKeyPressed(int, Qt::ButtonState)), SLOT(OnArrowKeyPressed(int, Qt::ButtonState)));

	success &= worldManager->connect((const QObject*)smf->ui.renderingParametersAction, SIGNAL(activated()), SLOT(OnEditRenderingParameters()));
	success &= worldManager->connect((const QObject*)smf->ui.glLightsAction, SIGNAL(activated()), SLOT(OnEditGLLights()));

	success &= worldManager->connect(smf, SIGNAL(selectionMode(int)), SLOT(OnSelectionMode(int)));

	success &= worldManager->connect(smf->ui.world, SIGNAL(mydrawNeeded()), SLOT(OnDraw()));
	success &= worldManager->connect(smf->ui.world, SIGNAL(drawWithNamesNeeded()), SLOT(OnSelectDraw()));
	success &= worldManager->connect(smf->ui.world, SIGNAL(selected(const unsigned int, Qt::KeyboardModifiers)), SLOT(OnObjectSelected(const unsigned int, Qt::KeyboardModifiers)));
	success &= worldManager->connect(smf->ui.world, SIGNAL(viewerInitialized()), SLOT(OnViewerInitialized()));
	success &= worldManager->connect(smf->ui.world, SIGNAL(windowResized(int,int)), SLOT(OnWorldSizeChanged(int,int)));

	success &= g_sdfmain->connect(worldManager, SIGNAL(SceneParametersChanged(int,float,float,float,float)), SLOT(ChangeSceneParameters(int,float,float,float,float)));

	success &= smf->statusBar()->connect(worldManager, SIGNAL(message(const QString&)), SLOT(message(const QString&)));
	success &= smf->statusBar()->connect(worldManager, SIGNAL(message(const QString&, int)), SLOT(message(const QString&, int)));

	success &= application.connect((const QObject*)smf->ui.fileExitAction, SIGNAL(activated()), SLOT(quit()));

	return success;
}