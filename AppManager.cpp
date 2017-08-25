/**
 * @ author: Ren-Ting Hsiao
 */

#include "stdafx.h"

#include "AppManager.h"
#include "AppObject.h"
#include "Watermark.h"


typedef CGAL::Simple_cartesian<double> My_kernel;

// ProgSettings "mesh_directory" -> file dir
// AppObject->id() -> file full name

void AppObject::unloadMesh()
{
	delete mesh;
	mesh = NULL;
}

AppManager::AppManager(MainWindow* parent, const char* name) :
QObject(parent, name),
m_worldManager(NULL)
{
	m_appParamsDialog = new AppParams((QWidget*)parent);
	connect(m_appParamsDialog, SIGNAL(pressdOk()), this, SLOT(appParamsOk()));
	connect(m_appParamsDialog, SIGNAL(pressdApply()), this, SLOT(appParamsApply()));


	loadApplicationParameters();

	if (parent->m_browse != NULL)
		connect(parent->m_browse, SIGNAL(openFile(QString)), this, SLOT(openFile(QString)));

}

void AppManager::openFile(QString filename)
{
	removeAllObjects();
	loadObject(filename, 1.0);
}

AppManager::~AppManager()
{
	saveApplicationParameters();

	//removeAllObjects();
}

void AppManager::Init(WorldManager* worldManager)
{
	m_worldManager = worldManager;

}

QString AppManager::createObjectId(const QString& fileName)
{
	int instancesFound = 0;
	QFileInfo fi(fileName);
	QString id = fi.fileName();

	for (appObjectsMap_t::iterator it = m_objects.begin();it!=m_objects.end();it++)
	{
		QFileInfo tempinfo(it.value()->inFileName);
		if (tempinfo.fileName() == id)
			instancesFound++;
	}

	// for the first, without a number
	if (instancesFound > 0)
	{
		int counter = 1;
		QString tempName = id + QString::number(instancesFound+counter);
		while (m_objects.find(tempName.ascii()) != m_objects.end()) {
			counter++;
			tempName = id + QString::number(instancesFound+counter);
		}
		id = tempName;
	}

	return id;
}

void AppManager::doMeshPrecalc(Mesh* mesh)
{
	assert (mesh != NULL);

	mesh->computeNormals();
	mesh->computeBoundingBox();
	mesh->computeTriangleSurfaces();
	mesh->centerOfMass();
	mesh->compute_type();

	//mesh->estimateCurvature();

	//m_mesh->computeVolumeSDF);
	//for (Mesh::Vertex_iterator it = mesh->vertices_begin(); it != mesh->vertices_end(); it++)
	//{
		//it->volumeSDF(FLT_MAX);
		//it->centricity(FLT_MAX);
	//}
}

void AppManager::selectObject(const QString& id)
{
	for (appObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end(); it++)
	{
		bool thisone = (it.key() == id);
		it.value()->meshSelection->setSelected(thisone);
	}

	m_selectedObject = id;
}


#include <CGAL/IO/Polyhedron_iostream.h>

typedef Enriched_polyhedron<Enriched_kernel,Enriched_items>::HalfedgeDS MyHDS;
typedef MyHDS::Vertex::Point MyPoint;

typedef CGAL::Simple_cartesian<double> OKernel;
typedef CGAL::Polyhedron_3<OKernel> OPolyhedron;


/*
bool AppManager::loadObjFile(const QString& fileName, Mesh* mesh)
{
	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly))
		return false;
	QTextStream in(&file);

	CGAL::Polyhedron_incremental_builder_3<MyHDS> builder(mesh->get_hds());

	int vindex = 0, findex = 0, filtered = 0;
	QString line;
	do {
		line = in.readLine();
		if (line.length() == 0)
			continue;
		QStringList flds = line.split(' ', QString::SkipEmptyParts);
		if (line[0] == 'v')
		{
			Mesh::Vertex_handle v = builder.add_vertex(MyPoint(flds[1].todouble(), flds[2].todouble(), flds[3].todouble()));
			v->index() = vindex++;
		}
		else if (line[0] == 'f')
		{
			size_t vertices[3] = { flds[1].toInt() - 1, flds[2].toInt() - 1, flds[3].toInt() - 1};

			if (!builder.test_facet(&(vertices[0]), &(vertices[3])))
			{
				++filtered;
				continue;
			}

			Mesh::Halfedge_handle he = builder.add_facet(&(vertices[0]), &(vertices[3]));
			if (he == NULL || he->facet() == NULL)
				return false;
			he->facet()->index(findex++);
		}

	} while (!line.isNull());

	printf("filtered=%d ", filtered);

	return true;
}*/

void AppManager::copyObject(AppObject *from, TStubData* newPart, DrawMethod dm)
{
	AppObject* appObject = new AppObject;
	appObject->mesh = from->mesh;
	appObject->meshSelection = new MeshSelection;
	appObject->description = from->description;
	appObject->inFileName = from->inFileName;
	appObject->id = createObjectId(appObject->inFileName);
	appObject->hasSdfFacets = from->hasSdfFacets;

	appObject->refCount = from->refCount;
	++(*appObject->refCount);

	m_objects.insert(appObject->id, appObject);
	selectObject(appObject->id);

	emit message("Copied <" + appObject->description + ">");
	emit ObjectAdded(appObject, newPart, dm);
}

QString makeSubinFileName(const QFileInfo& fi, QString sub, bool create)
{
	QDir dir(fi.path());
	if ((!dir.exists(sub)) && create)
		dir.mkdir(sub);
	return fi.path() + "\\" + sub + "\\" + fi.fileName();
}

QString makeSubFilename(const QFileInfo& fi, QString sub, bool create)
{
	QDir dir(fi.path());
	if ((!dir.exists(sub)) && create)
		dir.mkdir(sub);
	return fi.path() + "\\" + sub + "\\" + fi.fileName();
}


bool AppManager::loadObject(const QString& fileName, double err, TStubData* part, DrawMethod dm)
{
	//create app object
	AppObject* appObject = new AppObject;
	appObject->mesh = new Mesh;
	appObject->meshSelection = new MeshSelection;

	QFileInfo fi(fileName);
	/** setting some information about this mesh */
	appObject->description = fi.fileName();
	appObject->id = createObjectId(fileName);

	appObject->fileDir = fi.dirPath(TRUE);
	appObject->sdfDir = "sdfs";
	appObject->fileSeparator = "\\"; ///const information
	appObject->inFileName = fileName; //cout << "infile name : " << qPrintable(appObject->inFileName) << endl;
	appObject->fileName = fi.baseName(); //cout << "file name : " << qPrintable(appObject->fileName) << endl;
	
	//appObject->fileSuffix = "." + fi.completeSuffix();
	appObject->fileSuffix = "." + fi.suffix(); //cout << "file suffix : " << qPrintable(appObject->fileSuffix) << endl;

	//if marked object as normals reversed, do it now
	bool reverseNormals;
	ProgSettings ps;
	reverseNormals = m_appParameters.normalsReversed;
	ps.beginGroup("Reverse Normal Models");
	if (ps.readBoolEntry(appObject->inFileName))
		reverseNormals = !reverseNormals;
	ps.endGroup();

	//actually load the mesh
	bool success = ImporterMesh::loadMesh(fileName, appObject, err);

	if (success)
	{
		//if marked object as normals reversed, do it now
		if (reverseNormals)
		{
			appObject->mesh->inside_out();
			appObject->mesh->reverseNormals();
		}

		// check if to load default volume file		
		QString vsdfName = appObject->fileDir + appObject->fileSeparator + appObject->sdfDir + appObject->fileSeparator + appObject->fileName + "_vsdf.txt"; // plain file with vsdf extension
		QString vnsdfName = appObject->fileDir + appObject->fileSeparator + appObject->sdfDir + appObject->fileSeparator + appObject->fileName + "_vnsdf.txt"; // plain file with vnsdf extension
		QString fsdfName = appObject->fileDir + appObject->fileSeparator + appObject->sdfDir + appObject->fileSeparator + appObject->fileName + "_fsdf.txt"; // plain file with fsdf extension
		QString fnsdfName = appObject->fileDir + appObject->fileSeparator + appObject->sdfDir + appObject->fileSeparator + appObject->fileName + "_fnsdf.txt"; // plain file with fnsdf extension
		
		if (QFile::exists(vsdfName))
			ImporterMeshFeatures::loadSDFText(vsdfName, true, appObject);
		if (QFile::exists(vnsdfName))
			ImporterMeshFeatures::loadNSDFText(vnsdfName, true, appObject);
		if (QFile::exists(fsdfName))
			ImporterMeshFeatures::loadSDFText(fsdfName, false, appObject);
		if (QFile::exists(fnsdfName))
			ImporterMeshFeatures::loadNSDFText(fnsdfName, false, appObject);

		//add to object map
		m_objects.insert(appObject->id, appObject);
	
		selectObject(appObject->id);

		//send notifications
		emit message("Loaded <" + appObject->description + ">");

		printf("Loaded <%s>\n", appObject->description.toAscii().data());

		emit ObjectAdded(appObject, part, dm); //, appObject->id, appObject->description, appObject->mesh, part, dm, appObject->meshSurfaceGraph, appObject->meshSelection, tabact);
	}
	else
	{
		delete appObject;
	}

	return success;
}

bool AppManager::loadObject()
{
	///select file name
	ProgSettings ps;

	QString fileName = Q3FileDialog::getOpenFileName (
		ps.readEntry("mesh_directory", QString::null),
		"Mesh (*.obj *.off)", /*"Mesh (*.ply2 *.simple *.obj *.off)",*/
		NULL,
		"Load Mesh"
		"Choose mesh file" );

	cout << "selected file : " << qPrintable(fileName) << endl;

	bool allSuccess = loadObject(fileName, 1.0);

	if (!allSuccess)
		QMessageBox::warning(g_main, "Error loading", "Could not load mesh");

	return allSuccess;
}

QString AppManager::selectFileName()
{
	///select file names
	ProgSettings ps;

	QString fileName = Q3FileDialog::getOpenFileName (
		ps.readEntry("mesh_directory", QString::null),
		"Mesh (*.obj *.off)", /*"Mesh (*.ply2 *.simple *.obj *.off)",*/
		NULL,
		"Load Mesh"
		"Choose mesh file");

	return fileName;
}

QStringList AppManager::selectFileNames()
{
	///select file names
	ProgSettings ps;

	QStringList fileNames = Q3FileDialog::getOpenFileNames(
		"Mesh(es) (*.obj *.off)",/*"Mesh (*.ply2 *.simple *.obj *.off)",*/
		ps.readEntry("mesh_directory", QString::null),
		NULL,
		"Load Mesh(es)"
		"Choose mesh file(s)");

	return fileNames;
}

bool AppManager::loadObjects()
{
	///select file name
	ProgSettings ps;

	QStringList fileNames = Q3FileDialog::getOpenFileNames(
		"Mesh(es) (*.obj *.off)",/*"Mesh (*.ply2 *.simple *.obj *.off)",*/
		ps.readEntry("mesh_directory", QString::null),
		NULL,
		"Load Mesh(es)"
		"Choose mesh file(s)" );

	bool allSuccess = true;

	QString thisdir;
	if (fileNames.size())
	{
		for (QStringList::iterator it = fileNames.begin(); it != fileNames.end();it++)
		{
			QString fileName = *it;
			QFileInfo fi(fileName);
			thisdir = fi.dirPath(TRUE);
			if (!fi.exists())
			{
				allSuccess = false;
			} 
			else if (!loadObject(fileName, 1.0))
			{
				QMessageBox::warning(g_main, "Error loading", "Could not load mesh");
			}
		}

		if (allSuccess) 
		{
			//save directory for next time
			ps.writeEntry("mesh_directory",thisdir);
		}
	}

	return allSuccess;
}

class VolumeFaceColorizer : public FaceColorizer
{
public:
	VolumeFaceColorizer(WorldManager *wm) { m_cs = ColorSchemeManager::GetScheme(wm->getColorSchemeName()); }
	QColor getColor(Mesh::Facet_handle f) { return m_cs.gradient(f->volumeNSDF()); }
	ColorScheme m_cs;
};

uint qHash(const QColor& v)
{
	return (uint)v.rgb();
}


void AppManager::removeObject(QString id)
{
	appObjectsMap_t::iterator it = m_objects.find(id.ascii());

	if (it != m_objects.end()) {
		AppObject* appObject = it.value();

		if (m_selectedObject == appObject->id)
			selectObject(QString::null);

		delete appObject;

		m_objects.erase(it);
	}

	if (m_selectedObject == QString::null && m_objects.size()) {
		selectObject(m_objects.begin().value()->id);
	}

	emit ObjectRemoved(id);
}

void AppManager::removeAllObjects()
{
	//delete all objects from world
	for (appObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end(); it++)
	{
		AppObject* appObject = it.value();

		QString id = appObject->id;
		delete appObject;

		emit ObjectRemoved(id);
	}

	m_objects.clear();

	selectObject(QString::null);

	emit ObjectsChanged();
}

void AppManager::OnLoadObjects()
{
	loadObjects();
}

void AppManager::OnLoadObjectCage()
{
	//temporary selected object
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	QString filter = appObject->inFileName;
	int beginAt = filter.findRev("\\");
	if (beginAt >= 0)
	{
		int length = filter.find(".", beginAt) - beginAt - 1;
		filter = filter.mid(beginAt+1, length);
		filter.append("*.*");
	}
	else
	{
		filter = "*.*";
	}

	QString dirOnly = appObject->inFileName;
	dirOnly = dirOnly.left(beginAt+1);

	QString fileName;
	if ((fileName =
		Q3FileDialog::getOpenFileName(dirOnly + filter,
		"Cage File (*_cage.obj)",
		(QWidget*)parent(),
		NULL,
		"Open cage file")) != QString::null)
	{
		loadCage(fileName, appObject);
	}
}

void AppManager::OnSaveObject()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	QStringList meshTypes;
	meshTypes.append("obj");
	meshTypes.append("off");

	//featureTypes.append("obj with SDF Facets");

	//featureTypes.append("SDF Vertices");
	//featureTypes.append("SDF Facets");

	//featureTypes.append("NSDF Vertices");
	//featureTypes.append("NSDF Facets");
	
	bool ok;
	QString answer = QInputDialog::getItem(
		"Save Mesh",
		"Select mesh type to save",
		meshTypes, 5, FALSE, &ok,
		(QWidget*)parent(), "meshTypeSelect");

	if (!ok)
		return;

	int meshType = meshTypes.findIndex(answer);

	char* exts[] = { ".obj", ".off" };

	QString ext = exts[meshType];

	/*QString filter = appObject->inFileName;
	int beginAt = filter.lastIndexOf("\\");
	if (beginAt >= 0)
	{
		int length = filter.find(".", beginAt) - beginAt - 1;
		filter = filter.mid(beginAt+1, length);
		filter.append(ext);
	}
	else
	{
		filter = "*." + ext;
	}*/

	//QString dirOnly = appObject->inFileName;
	QString dirOnly = appObject->fileDir + appObject->fileSeparator + appObject->fileName + ext;

	QString fileName;
	if ((fileName = Q3FileDialog::getSaveFileName(dirOnly/* + filter*/, "Mesh Files (*" + ext + ")", (QWidget*)parent(), NULL, "Save mesh file")) != QString::null)
	{
		switch(meshType)
		{
			case 0: ExporterMeshFile::write_obj(fileName, appObject); break;
			case 1: ExporterMeshFile::write_off(fileName, appObject); break;
			/*
			case 2: ExporterMeshFile::write_ply(fileName, appObject); break;
			case 3: ExporterMeshFile::write_VRML1(fileName, appObject); break;
			case 4: ExporterMeshFile::write_VRML2(fileName, appObject); break;
			*/

			//case 1: saveObjWithFaceSDF(fileName, appObject); break;
		}
	}
}

void AppManager::OnLoadObjectFeatures()
{
	//temporary selected object
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	QString filter = appObject->inFileName;
	int beginAt = filter.findRev("\\");
	if (beginAt >= 0)
	{
		int length = filter.find(".", beginAt) - beginAt - 1;
		filter = filter.mid(beginAt+1, length);
		filter.append("*.*");
	}
	else
	{
		filter = "*.*";
	}

	QString dirOnly = appObject->inFileName;
	dirOnly = dirOnly.left(beginAt+1);

	QString fileName;
	if ((fileName =
		Q3FileDialog::getOpenFileName(dirOnly + filter,
		"Feature Files (*.vsdf *.vnsdf *.fsdf *.fnsdf)",
		(QWidget*)parent(),
		NULL,
		"Open SDF file")) != QString::null)
	{
		loadMeshFeatures(fileName, appObject);
	}
}

void AppManager::OnSaveObjectFeatures()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	QStringList featureTypes;

	//featureTypes.append("obj");

	//featureTypes.append("Snapshot");

	//featureTypes.append("Obj with SDF Facets");

	//featureTypes.append("SDF Vertices");
	//featureTypes.append("SDF Facets");
	featureTypes.append("SDF Vertices");
	featureTypes.append("SDF Facets");
	featureTypes.append("SDF Vertices Diff");
	featureTypes.append("SDF Facets Diff");

	//featureTypes.append("NSDF Vertices");
	//featureTypes.append("NSDF Facets");
	featureTypes.append("NSDF Vertices");
	featureTypes.append("NSDF Facets");
	featureTypes.append("NSDF Vertices Diff");
	featureTypes.append("NSDF Facets Diff");
	
	featureTypes.append("VSI Ref Point");
	featureTypes.append("VSI Vertices");
	featureTypes.append("VSI Facets");

	featureTypes.append("Vertex Max Diameter");
	featureTypes.append("Vertex Min Diameter");
	featureTypes.append("Facet Max Diameter");
	featureTypes.append("Facet Min Diameter");

	featureTypes.append("Vertex Max Diameter Volume");
	featureTypes.append("Vertex Min Diameter Volume");
	featureTypes.append("Facet Max Diameter Volume");
	featureTypes.append("Facet Min Diameter Volume");

	featureTypes.append("DDegree");
	featureTypes.append("NDDegree");

	featureTypes.append("Vertex Normal");
	featureTypes.append("Facet Normal");

	bool ok;
	QString answer = QInputDialog::getItem(
		"Save Mesh Features",
		"Select feature to save",
		featureTypes, 5, FALSE, &ok,
		(QWidget*)parent(), "featureFileTypeSelect");

	if (ok)
	{
		int featureType = featureTypes.findIndex(answer);

		//char* exts[] = { ".obj", ".png", ".sdf.obj", ".vsdf", ".sdf", ".vsdf.txt", ".sdf.txt", ".sdf_diff.txt", ".vnsdf", ".nsdf", ".vnsdf.txt", ".nsdf.txt", ".nsdf_diff.txt", ".dd.txt", ".ndd.txt" };
		char* exts[] = { "_vsdf.txt", "_fsdf.txt", "_vsdf_diff.txt", "_fsdf_diff.txt", 
									"_vnsdf.txt", "_fnsdf.txt", "_vnsdf_diff.txt", "_fnsdf_diff.txt", 
									"_dd.txt", "_ndd.txt",
									"_vnormal.txt", "_fnormal.txt"};

		QString ext = exts[featureType];

		QString filter = appObject->inFileName;
		int beginAt = filter.findRev("\\");
		if (beginAt >= 0)
		{
			int length = filter.find(".", beginAt) - beginAt - 1;
			filter = filter.mid(beginAt+1, length);
			filter.append(ext);
		}
		else
		{
			filter = "*." + ext;
		}

		QString dirOnly = appObject->inFileName;

		dirOnly = dirOnly.left(beginAt+1);

		QString fileName;
		if ((fileName = Q3FileDialog::getSaveFileName(dirOnly + filter, "Feature Files (*." + ext + ")", (QWidget*)parent(), NULL, "Save SDF file")) != QString::null)
		{
			switch(featureType)
			{
				ExporterMeshFeature exporterMeshFeature;

				//saveSDFVertices(fileName, appObject); break;
				//saveSDFFacets(fileName, appObject); break;

				case  0: exporterMeshFeature.saveSDFText(fileName, true, appObject); break;
				case  1: exporterMeshFeature.saveSDFText(fileName, false, appObject); break;	
				case  2: exporterMeshFeature.saveSDFDiffText(fileName, true, appObject); break;
				case  3: exporterMeshFeature.saveSDFDiffText(fileName, false, appObject); break;

				case  4: exporterMeshFeature.saveNSDFText(fileName, true, appObject); break;
				case  5: exporterMeshFeature.saveNSDFText(fileName, false, appObject); break;
				case  6: exporterMeshFeature.saveNSDFDiffText(fileName, true, appObject); break;
				case  7: exporterMeshFeature.saveNSDFDiffText(fileName, false, appObject); break;

				case  8: exporterMeshFeature.saveDDegreeText(fileName, appObject); break;
				case  9: exporterMeshFeature.saveNDDegreeText(fileName, appObject); break;

				case 10: exporterMeshFeature.saveNormals(fileName, true, appObject); break;
				case 11: exporterMeshFeature.saveNormals(fileName, false, appObject); break;
			}
		}
	}
}

void AppManager::OnRemoveSelectedObject()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	removeObject(appObject->id);
}

void AppManager::OnRemoveObjects()
{
	if (QMessageBox::question(
		(QWidget*)parent(),
		"Removing all models",
		"Are you sure you want to remove all models?",
		QMessageBox::Yes, QMessageBox::No) == QMessageBox::Yes)
	{
		removeAllObjects();
	}
}

AppObject* AppManager::getSelectedObject()
{
	if (m_selectedObject == QString::null)
		return NULL;

	//	AppObject* appObject;
	appObjectsMap_t::iterator it = m_objects.find(m_selectedObject.ascii());
	if (it != m_objects.end())
	{
		return it.value();
	}
	else
	{
		return NULL;
	}
}

void AppManager::OnDiffObjectFeatures()
{
	//temporary selected object
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	QString filter = appObject->inFileName;
	int beginAt = filter.findRev("\\");
	if (beginAt >= 0) {
		int length = filter.find(".", beginAt) - beginAt - 1;
		filter = filter.mid(beginAt+1, length);
		filter.append("*.*");
	} else {
		filter = "*.*";
	}

	QString dirOnly = appObject->inFileName;
	dirOnly = dirOnly.left(beginAt+1);

	QString fileName;
	if ((fileName =
		Q3FileDialog::getOpenFileName(dirOnly + filter,
		"Feature Files (*.sdf *.fsdf)",
		(QWidget*)parent(),
		NULL,
		"Open SDF file")) != QString::null)
	{
		ImporterMeshFeatures::diffSDFFacets(fileName, appObject);
	}

	emit ObjectChanged(appObject->id, false);
}


bool AppManager::loadMeshFeatures(const QString& featuresFileName, AppObject* appObject)
{
	assert(appObject != NULL);

	QFileInfo fi(featuresFileName);
	if (!fi.exists()) return false;

	bool success = false;

	QString ext = fi.suffix();
	if (ext == "vsdf")
	{
		// vertices SDF
		success = ImporterMeshFeatures::loadSDFVertices(featuresFileName, false, appObject); 
	}
	else if (ext == "vnsdf")
	{
		// vertices NSDF
		success = ImporterMeshFeatures::loadSDFVertices(featuresFileName, true, appObject);
	}
	else if (ext == "fsdf")
	{
		// on faces, facets SDF
		success = ImporterMeshFeatures::loadSDFFacets(featuresFileName, false, appObject);
		appObject->hasSdfFacets = true;
	}
	else if (ext == "fnsdf")
	{
		// on faces, facets NSDF
		success = ImporterMeshFeatures::loadSDFFacets(featuresFileName, true, appObject);
		appObject->hasSdfFacets = true;
	}

	if (success)
	{
		emit message("Loaded features for <" + appObject->description  + "> from <" + fi.fileName() + ">");
		emit ObjectChanged(appObject->id, false);
	}

	return success;
}


bool AppManager::loadCage(const QString& cageFileName, AppObject* appObject, TStubData* part, DrawMethod dm)
{
	assert(appObject != NULL);

	QFileInfo fi(cageFileName);
	if (!fi.exists()) return false;

	/// load cage mesh
	Mesh *cage = new Mesh;
	
	bool success = false;

	// actually load mesh
	success = ImporterMesh::loadMesh(cageFileName, cage, 1.0);

	if (success)
	{
		appObject->cage = cage;

		emit message("Loaded cage for <" + appObject->description + "> from <" + fi.fileName() + ">");
		printf("Loaded cage for <" + appObject->description + "> from <" + fi.fileName() + ">");

		emit ObjectChanged(appObject->id, false);

		emit CageAdded(appObject, part, dm);
	}

	return success;
}


bool AppManager::embedWatermark(AppObject* appObject, const bool showDialog)
{
	Mesh* watermarkedMesh = appObject->mesh;
	
	bool success = false; 

	//////////////////////////////////////////////////////////////////////////
	// Load cover mesh
	//////////////////////////////////////////////////////////////////////////
	QString originalMeshFileName = appObject->fileDir + appObject->fileSeparator + appObject->fileName + appObject->fileSuffix;

	Mesh* originalMesh = new Mesh;

	success = ImporterMesh::loadMesh(originalMeshFileName, originalMesh);

	if (!success)
		return false;
	//////////////////////////////////////////////////////////////////////////

	SDFWatermarkingAlgorithm sdfWatermarkingAlgorithm;

	success = sdfWatermarkingAlgorithm.embed(
			appObject, originalMesh, watermarkedMesh, 
			true, m_appParameters.multithreaded, m_appParameters.threadsNum,
			m_appParameters.watermark_seed, m_appParameters.watermark_bits, m_appParameters.robustness_factor, m_appParameters.watermark_searchSpace,
			m_appParameters.watermark_maxCycles, m_appParameters.watermark_particleNumbers,
			m_appParameters.lsd_gridSize, m_appParameters.lsd_numberOfCones, m_appParameters.lsd_coneSeparation, m_appParameters.lsd_raysInCone, 
			m_appParameters.lsd_removeSameDirectionIntersection, m_appParameters.lsd_gaussianWeights, m_appParameters.lsd_removeOutlier, 
			m_appParameters.lsd_normalize, m_appParameters.lsd_log_alpha,
			m_appParameters.lsd_smoothing, m_appParameters.lsd_smoothing_iterations,
			m_appParameters.debugMode, showDialog);

	delete originalMesh;

	if (!success)
		return false;

	QString s;
	s.sprintf("Embed watermark for <%s>", appObject->id.toAscii());
	SDFLOG8(s);

	emit message("Embed watermark for <" + appObject->description + ">");

	cout << "Embed watermark for <" << qPrintable(appObject->description) << ">" << endl << endl;;

	//emit ObjectChanged(appObject->id, false);

	//emit ObjectsChanged();

	return true;
}


bool AppManager::detectWatermark(AppObject* appObject, const bool showDialog)
{
	Mesh* watermarkedMesh = appObject->mesh;

	SDFWatermarkingAlgorithm sdfWatermarkingAlgorithm;

	bool success = sdfWatermarkingAlgorithm.detect(
			appObject, watermarkedMesh, true, 
			m_appParameters.watermark_seed, m_appParameters.watermark_bits, 
			m_appParameters.lsd_gridSize, m_appParameters.lsd_numberOfCones, m_appParameters.lsd_coneSeparation, m_appParameters.lsd_raysInCone, 
			m_appParameters.lsd_gaussianWeights, m_appParameters.lsd_removeOutlier, m_appParameters.lsd_normalize, m_appParameters.lsd_log_alpha,
			m_appParameters.lsd_smoothing, m_appParameters.lsd_smoothing_iterations,
			m_appParameters.debugMode, showDialog);

	if (!success)
		return false;

	QString s;
	s.sprintf("Detect watermark for <%s>", appObject->id.toAscii());
	SDFLOG8(s);

	emit message("Detect watermark for <" + appObject->description + ">");

	//emit ObjectChanged(appObject->id, false);

	//emit ObjectsChanged();

	return true;
}

bool AppManager::calculateTotalVolume(AppObject* appObject)
{
	Mesh* mesh = appObject->mesh;

	vector<Point_3> points(3);
	int index;
	double totalVolume = 0.0;

	Mesh::Facet_iterator fit = mesh->facets_begin();
	Mesh::Facet_iterator fit_end = mesh->facets_end();
	for (; fit != fit_end; fit++)
	{
		assert(fit->size() == 3);

		index = 0;

		Mesh::Halfedge_around_facet_const_circulator pHalfedge = fit->facet_begin();
		Mesh::Halfedge_around_facet_const_circulator end = pHalfedge;
		CGAL_For_all(pHalfedge, end)
		{
			points[index] = pHalfedge->vertex()->point();

			index++;
		}

		double det = (points[1].x() - points[0].x()) * (points[2].y() - points[0].y()) - (points[2].x() - points[0].x()) * (points[1].y() - points[0].y());
		double localVolume = (points[0].z() + points[1].z() + points[2].z()) * det;

		fit->volume() = localVolume;

		totalVolume += localVolume;
	}

	Mesh::Vertex_iterator vit = mesh->vertices_begin();
	Mesh::Vertex_iterator vit_end = mesh->vertices_end();
	for (; vit != vit_end; vit++)
	{
		double sum = 0.0;		
		int total = 0;
		Mesh::Halfedge_around_vertex_const_circulator pHalfedge = vit->vertex_begin();
		Mesh::Halfedge_around_vertex_const_circulator end = pHalfedge;
		CGAL_For_all(pHalfedge, end);
		{
			if (pHalfedge->opposite()->facet() != NULL)
			{
				sum += pHalfedge->opposite()->facet()->volume();
				total++;
			}
		}

		vit->volume() = sum / total;
	}

	printf("total volume = %g\n", totalVolume);

	bool success = true;

	QString exportFacetFileName = appObject->fileDir + appObject->fileSeparator + "volume_facet_" + appObject->fileName + ".txt"; // save deformed mesh
	cout << "export file to: " << qPrintable(exportFacetFileName) << endl;

	QFile file(exportFacetFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);
		
		Mesh::Facet_iterator fit = mesh->facets_begin();
		Mesh::Facet_iterator fit_end = mesh->facets_end();
		for (; fit != fit_end; fit++)
			ts << fit->volume() << '\n';
		
		file.close();

		success &= true;
	}	
	else
	{
		success &= false;
	}
	
	QString exportVertexFileName = appObject->fileDir + appObject->fileSeparator + "volume_vertex_" + appObject->fileName + ".txt"; // save deformed mesh
	cout << "export file to: " << qPrintable(exportVertexFileName) << endl;

	QFile file1(exportVertexFileName);
	if (file1.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file1);
		ts.setRealNumberPrecision(10);
		
		Mesh::Vertex_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
			ts << vit->volume() << '\n';
		
		file1.close();

		success &= true;
	}
	else
	{
		success &= false;
	}

	return success;
}


bool AppManager::calculateShapeDiameterFunction(
	AppObject* appObject,
	const bool onVertices, 
	const bool showDialog)
{
	//calculate SDF
	CalculateSDF csdf;

	bool success = csdf.go(
		appObject->mesh,
		onVertices,
		m_appParameters.debugMode,
		m_appParameters.multithreaded,
		m_appParameters.threadsNum,
		m_appParameters.lsd_gridSize,
		m_appParameters.lsd_numberOfCones,
		m_appParameters.lsd_coneSeparation,
		m_appParameters.lsd_raysInCone,
		m_appParameters.lsd_removeSameDirectionIntersection,
		m_appParameters.lsd_gaussianWeights,
		m_appParameters.lsd_removeOutlier,
		m_appParameters.lsd_normalize,
		m_appParameters.lsd_log_alpha,
		m_appParameters.lsd_smoothing,
		m_appParameters.lsd_smoothing_iterations,
		showDialog);

	if (!success)
		return false;

	QString s;
	s.sprintf("Calculated SDF for <%s>", appObject->id.toAscii());
	SDFLOG8(s);

	emit message("Calculated SDF for <" + appObject->description + ">");

	emit ObjectChanged(appObject->id, false);

	emit ObjectsChanged();

	return true;
}


bool AppManager::calculateCurvature(AppObject* appObject)
{
	/*CalculateCurvature calculateCurvature;

	//Calculate curvature
	bool success = calculateCurvature.go(
		appObject->mesh,
		m_appParameters.multithreaded);

	if(!success)
		return false;*/

	Mesh* mesh = appObject->mesh;

	Mesh::Vertex_iterator it = mesh->vertices_begin();
	Mesh::Vertex_iterator it_end = mesh->vertices_end();
	for (; it != it_end; it++)
	{
		double cur = mesh->average_curvature_around(it);
		
		it->curvature(cur);

		printf("%f\n", cur);
	}

	return true;
}


bool AppManager::generateCageMesh(AppObject* appObject)
{
	const double enlarge_distance = 0.25;

	Mesh* mesh = appObject->mesh;
	
	mesh->computeNormals();

	for (Mesh::Vertex_iterator ivertex = mesh->vertices_begin(); ivertex != mesh->vertices_end(); ivertex++)
	{
		Mesh::Vertex_handle v = ivertex;

		Point_3 point = v->point();
		Vector_3 vertex_normal = v->normal();

		Vector_3 distance_vector = enlarge_distance * vertex_normal;

		Point_3 cage_point(point.x() + distance_vector.x(), point.y() + distance_vector.y(), point.z() + distance_vector.z());

		v->point() = cage_point;
	}

	/// output file
	QString cageObjFileName = appObject->fileDir + "\\" + appObject->fileName + ".cage.obj";
	mesh->write_obj(cageObjFileName, 1);

	printf("Generate cage mesh done!\n");

	return true;
}


bool AppManager::calculateProjectionDifference(AppObject* appObject)
{
	Mesh* mesh = appObject->mesh;
	
	mesh->computeNormals();

	Mesh::Vertex_iterator vit = mesh->vertices_begin();
	Mesh::Vertex_iterator vit_end = mesh->vertices_end();

	for (; vit != vit_end; vit++)
	{
		Point_3 p = vit->point(); // point
		Vector_3 normal = vit->normal(); // point normal

		// printf("(%f, %f, %f)\n", normal.x(), normal.y(), normal.z());

		//double normal_length = sqrt(normal.squared_length()); // point normal length

		//printf("normal length = %d\n", normal_length);

		Mesh::Halfedge_around_vertex_circulator pHalfEdge = vit->vertex_begin();
		Mesh::Halfedge_around_vertex_circulator end = pHalfEdge;

		double sum_diff_length = 0;
		CGAL_For_all(pHalfEdge, end)
		{
			Point_3  p_neigh = pHalfEdge->opposite()->vertex()->point();

			// find projection point
			double d = normal.x() * p.x() + normal.y() * p.y() + normal.z() * p.z();
			double proj_length = (normal.x() * p_neigh.x() + normal.y() * p_neigh.y() + normal.z() * p_neigh.z() - d)/* / normal_length*/;
			
			Point_3 p_proj(p_neigh.x() - proj_length*normal.x(), p_neigh.y() - proj_length*normal.y(), p_neigh.z() - proj_length*normal.z());

			double diff_length = sqrt(pow(p.x() - p_proj.x(), 2) + pow(p.y() - p_proj.y(), 2) + pow(p.z() - p_proj.z(), 2));

			sum_diff_length += diff_length;
		}

		vit->diff_length(sum_diff_length);
	}

	// calculate difference length done
	QString exportFileName = appObject->fileDir + "\\" + appObject->fileName + "_diff_length.txt"; // save deformed mesh
	cout << "export file to: " << qPrintable(exportFileName) << endl;

	QFile file(exportFileName);
	if (file.open(QIODevice::WriteOnly))
	{
		QTextStream ts(&file);
		ts.setRealNumberPrecision(10);
		
		Mesh::Vertex_iterator vit = mesh->vertices_begin();
		Mesh::Vertex_iterator vit_end = mesh->vertices_end();
		for (; vit != vit_end; vit++)
			ts << vit->diff_length() << '\n';
		
		file.close();

		return true;
	}

	return false;
}


/*bool AppManager::laplaceBeltrami(AppObject* appObject)
{
	Mesh* mesh = appObject->mesh;

	LaplaceBeltrami lb;
	lb.go(mesh);

	return true;
}*/

bool AppManager::meanValueCoordinates(AppObject* appObject)
{
	QString cageObjFileName = appObject->fileDir + "\\" + appObject->fileName + ".cage.obj"; // load cage mesh
	QString deformCageObjFileName = appObject->fileDir + "\\" + appObject->fileName + ".deformcage.obj"; // load deformed cage mesh

	QFileInfo cage_file(cageObjFileName);
	QFileInfo deform_cage_file(deformCageObjFileName);
	if(!cage_file.exists())
	{
		QMessageBox::warning(g_main, "Can not found the file", cageObjFileName);
		return false;
	}
	if(!deform_cage_file.exists())
	{
		QMessageBox::warning(g_main, "Can not found the file", deformCageObjFileName);
		return false;
	}

	printf("file exists!\n");

	bool success;

	Mesh *origin_mesh = appObject->mesh;

	Mesh *cage_mesh = new Mesh;
	success = FileParser::read_Obj(cageObjFileName.toAscii(), cage_mesh, 1.0);
	if (!success)
		return false;

	printf("read cage mesh done!\n");

	
	Mesh *deform_cage_mesh = new Mesh;
	success = FileParser::read_Obj(deformCageObjFileName.toAscii(), deform_cage_mesh, 1.0);
	if (!success)
		return false;

	printf("read deform cage mesh done!\n");

	MeanValueCoordinate mvc;
	success = mvc.go(appObject, origin_mesh, cage_mesh, deform_cage_mesh);
	if (!success)
		return false;

	printf("mean value coordinate done!\n");

	QString deformObjFileName = appObject->fileDir + "\\" + appObject->fileName + ".deform_mvc.obj"; // load cage mesh
	origin_mesh->write_obj(deformObjFileName, 1);

	return true;
}

bool AppManager::harmonicCoordinates(AppObject* appObject)
{
	return true;
}

bool AppManager::greenCoordinates(AppObject* appObject)
{
	QString cageObjFileName = appObject->fileDir + "\\" + appObject->fileName + ".cage.obj"; // load cage mesh
	QString deformCageObjFileName = appObject->fileDir + "\\" + appObject->fileName + ".deformcage.obj"; // load deformed cage mesh

	QFileInfo cage_file(cageObjFileName);
	QFileInfo deform_cage_file(deformCageObjFileName);
	if(!cage_file.exists())
	{
		QMessageBox::warning(g_main, "Can not found the file", cageObjFileName);
		return false;
	}
	if(!deform_cage_file.exists())
	{
		QMessageBox::warning(g_main, "Can not found the file", deformCageObjFileName);
		return false;
	}

	printf("file exists!\n");

	bool success;

	Mesh *origin_mesh = appObject->mesh;

	Mesh *cage_mesh = new Mesh;
	success = FileParser::read_Obj(cageObjFileName.toAscii(), cage_mesh, 1.0);
	if (!success)
		return false;
	
	printf("read cage mesh done!\n");


	Mesh *deform_cage_mesh = new Mesh;
	success = FileParser::read_Obj(deformCageObjFileName.toAscii(), deform_cage_mesh, 1.0);
	if (!success)
		return false;

	printf("read deform cage mesh done!\n");

	printf("calculate green coordinate...\n");

	GreenCoordinate gc;
	success = gc.go(appObject, origin_mesh, cage_mesh, deform_cage_mesh);
	if (!success)
		return false;

	printf("green coordinate done!\n");

	QString deformObjFileName = appObject->fileDir + "\\" + appObject->fileName + ".deform_gc.obj"; // save deformed mesh
	origin_mesh->write_obj(deformObjFileName, 1);

	return true;
}

void AppManager::OnGenerateCageMesh()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	//appObject->mesh->computeTotalVolume();

	//if (!generateCageMesh(appObject))
	//	printf("Failed on generate cage mesh");
}

void AppManager::OnMeanValueCoordinates()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	if (!meanValueCoordinates(appObject))
		printf("Failed on deforming using mean value coordinates\n");
}

void AppManager::OnHarmonicCoordinates()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	if (!harmonicCoordinates(appObject))
		printf("Failed on deforming using mean value coordinates\n");
}

void AppManager::OnGreenCoordinates()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	if (!greenCoordinates(appObject))
		printf("Failed on deforming using mean value coordinates\n");
}

void AppManager::OnCalculateSNR()
{
	QString originalMeshFileName = m_appParameters.fileName1;
	QString deformedMeshFileName = m_appParameters.fileName2;

	Mesh* originalMesh = new Mesh;
	Mesh* deformedMesh = new Mesh;
	
	if (!ImporterMesh::loadMesh(originalMeshFileName, originalMesh))
		return;
	if (!ImporterMesh::loadMesh(deformedMeshFileName, deformedMesh))
		return;

	Quality quality;
	double SNR = quality.measureSNR(originalMesh, deformedMesh);

	cout << "SNR = " << SNR << endl;

	delete originalMesh;
	delete deformedMesh;
}

void AppManager::OnCalculatePSNR()
{
	QString originalMeshFileName = m_appParameters.fileName1;
	QString deformedMeshFileName = m_appParameters.fileName2;

	Mesh* originalMesh = new Mesh;
	Mesh* deformedMesh = new Mesh;
	
	if (!ImporterMesh::loadMesh(originalMeshFileName, originalMesh))
		return;
	if (!ImporterMesh::loadMesh(deformedMeshFileName, deformedMesh))
		return;

	Quality quality;
	double PSNR = quality.measurePSNR(originalMesh, deformedMesh);

	cout << "PSNR = " << PSNR << endl;

	delete originalMesh;
	delete deformedMesh;
}

void AppManager::OnCalculateMRMS()
{
	QString originalMeshFileName = m_appParameters.fileName1;
	QString deformedMeshFileName = m_appParameters.fileName2;

	double MRMS;
	double MRMSwrtBB;
	double Hausdorff;
	double HausdorffwrtBB;

	Quality quality;
	quality.calculateGeoDistances(originalMeshFileName, deformedMeshFileName, MRMS, MRMSwrtBB, Hausdorff, HausdorffwrtBB);

	cout << "MRMS = " << MRMS << endl;
	cout << "MRMSwrtBB = " << MRMSwrtBB << endl;
}

void AppManager::OnCalculateHausdorffDistance()
{
	QString originalMeshFileName = m_appParameters.fileName1;
	QString deformedMeshFileName = m_appParameters.fileName2;

	double MRMS;
	double MRMSwrtBB;
	double Hausdorff;
	double HausdorffwrtBB;

	Quality quality;
	quality.calculateGeoDistances(originalMeshFileName, deformedMeshFileName, MRMS, MRMSwrtBB, Hausdorff, HausdorffwrtBB);

	cout << "Hausdorff = " << Hausdorff << endl;
	cout << "HausdorffwrtBB = " << HausdorffwrtBB << endl;
}

void AppManager::OnCalculateCurvature()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	if (!calculateCurvature(appObject))
		printf("Failed on calculate curvature\n");
}

void AppManager::OnCalculateTotalVolume()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	if (!calculateTotalVolume(appObject))
		printf("Failed on calculate total volume\n");
}

void AppManager::OnCalculateLSDVertices()
{
	try
	{
		AppObject* appObject = getSelectedObject();
		if (!appObject) return;

		if (!calculateShapeDiameterFunction(appObject, true))
			printf("Failed on calculate SDF Vertices\n");
	}
	catch (...)
	{
		QMessageBox::critical(NULL, "SDF Vertices", "Exception");
	}
}

void AppManager::OnCalculateLSDFacets()
{
	try
	{
		AppObject* appObject = getSelectedObject();
		if (!appObject) return;

		if (!calculateShapeDiameterFunction(appObject, false))
			printf("Failed on calculate SDF Facets\n");
	}
	catch (...)
	{
		QMessageBox::critical(NULL, "SDF Facets", "Exception");
	}
}

void AppManager::OnCalculateLSDVerticesForAll()
{
	for (appObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end(); it++)
	{
		AppObject* appObject = it.value();
		if (appObject)
		{
			if (!calculateShapeDiameterFunction(appObject, true))
				printf("Failed on calculate SDF Vertices\n");
		}
	}
}

void AppManager::OnCalculateLSDFacetsForAll()
{
	for (appObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end(); it++)
	{
		AppObject* appObject = it.value();
		if (appObject)
		{
			if (!calculateShapeDiameterFunction(appObject, false))
				printf("Failed on calculate SDF Facets\n");
		}
	}
}


bool AppManager::partitionFacetsMesh(AppObject* appObject)
{
	Mesh* mesh = appObject->mesh;

	PartitionMesh partitionMesh;
	partitionMesh.partition(mesh, false);

	return true;
}


bool AppManager::partitionVerticesMesh(AppObject* appObject)
{
	Mesh* mesh = appObject->mesh;

	PartitionMesh partitionMesh;
	partitionMesh.partition(mesh, true);

	return true;
}


void AppManager::allAttacks(AppObject* appObject)
{
	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

	printf("###Element reordering attack...");
	attack.ElementReordering();
	printf("Done!!!\n");
		
	printf("#Geometry quantization attack...");
	attack.CoordinateQuantization(7);
	attack.CoordinateQuantization(8);
	attack.CoordinateQuantization(9);
	attack.CoordinateQuantization(10);
	attack.CoordinateQuantization(11);
	printf("Done!!!\n");

	printf("#Laplacian smoothing attack...");
	attack.LaplacianSmoothing(0.1, 5, true);
	attack.LaplacianSmoothing(0.1, 10, true);
	attack.LaplacianSmoothing(0.1, 20, true);
	attack.LaplacianSmoothing(0.1, 30, true);
	attack.LaplacianSmoothing(0.1, 40, true);
	printf("Done!!!\n");

	printf("#Noise addition attack...");
	attack.NoiseAdditionUniform(0.0005, true);
	attack.NoiseAdditionUniform(0.001, true);
	attack.NoiseAdditionUniform(0.002, true);
	attack.NoiseAdditionUniform(0.003, true);
	attack.NoiseAdditionUniform(0.004, true);
	printf("Done!!!\n");

	printf("#Random Similarity transformation attack...");
	attack.SimilarityTransformation();
	printf("Done!!!\n");

	printf("#Simplification attack...");
	//attack.Simplification(EDGELENGTHMIDPOINT, 5);
	//attack.Simplification(EDGELENGTHMIDPOINT, 10);
	//attack.Simplification(EDGELENGTHMIDPOINT, 20);
	//attack.Simplification(EDGELENGTHMIDPOINT, 30);
	//attack.Simplification(EDGELENGTHMIDPOINT, 40);
	//attack.Simplification(EDGELENGTHMIDPOINT, 50);
	attack.Simplification(LINDSTROMTURK, 5);
	attack.Simplification(LINDSTROMTURK, 10);
	attack.Simplification(LINDSTROMTURK, 20);
	attack.Simplification(LINDSTROMTURK, 30);
	attack.Simplification(LINDSTROMTURK, 40);
	//attack.Simplification(LINDSTROMTURK, 50);
	printf("Done!!!\n");

	printf("#Subdivision attack...");
	attack.Subdivision(CATMULLCLARK, 1);
	//attack.Subdivision(LOOP, 1);
	//attack.Subdivision(DOOSABIN, 1);
	//attack.Subdivision(SQRT3, 1);
	attack.Subdivision(MIDPOINT, 1);
	printf("Done!!!\n");
}

void AppManager::noiseAddition(AppObject* appObject, const bool preserveBoundaries, const TNT::Array1D<double>& noiseIntensityArray)
{
	printf("\n#Noise addition attack...   ");

	const int size = noiseIntensityArray.dim();
	for(int noiseIndex = 0; noiseIndex < size; noiseIndex++)
	{
		double noise = noiseIntensityArray[noiseIndex];

		Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

		attack.NoiseAdditionUniform(noise, preserveBoundaries);
	}

	printf("done!\n");
}

void AppManager::quantization(AppObject* appObject, const TNT::Array1D<int>& quantizationBitArray)
{
	printf("\n#Geometry quantization attack...   ");

	const int size = quantizationBitArray.dim();
	for(int quantizationIndex = 0; quantizationIndex < size; quantizationIndex++)
	{
		int bitDepth = quantizationBitArray[quantizationIndex];

		Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

		attack.CoordinateQuantization(bitDepth);
	}

	printf("done!\n");
}

void AppManager::reordering(AppObject* appObject)
{
	printf("\n#Element reordering attack...   ");

	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

	attack.ElementReordering();

	printf("done!\n");
}

void AppManager::similarityTransform(AppObject* appObject)
{
	printf("\n#Random Similarity transformation attack...   ");

	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName, appObject->fileDir);

	attack.SimilarityTransformation();

	printf("done!\n");
}

void AppManager::translation(AppObject* appObject)
{
	printf("\n#Translation attack...   ");

	double xTranslation = 1;
	double yTranslation = 1;
	double zTranslation = 1;

	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

	attack.Translation_Save(xTranslation, yTranslation, zTranslation);

	printf("done!\n");
}

void AppManager::rotation(AppObject* appObject)
{
	printf("\n#Rotation attack...   ");

	double xAxis = 1;
	double yAxis = 1;
	double zAxis = 1;
	double angle = 0;

	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

	attack.Rotation_Save(xAxis, yAxis, zAxis, angle);

	printf("done!\n");
}

void AppManager::uniformScaling(AppObject* appObject)
{
	double scalingFactor = 3;

	printf("\n#Uniform scaling attack with scaling factor %f...   ", scalingFactor);

	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

	attack.UniformScaling_Save(scalingFactor);

	printf("done!\n");
}

void AppManager::nonuniformScaling(AppObject* appObject)
{
	double scalingFactorX= 1.2;
	double scalingFactorY = 1.2;
	double scalingFactorZ = 1.2;

	printf("\n#Non-Uniform scaling attack...   ");

	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

	attack.NonUniformScaling_Save(scalingFactorX, scalingFactorY, scalingFactorZ);

	printf("done!\n");
}

void AppManager::smoothing(
	AppObject* appObject, 
	const double deformFactor, const bool preserveBoundaries, const TNT::Array1D<int>& iterationList)
{
	printf("\n#Laplacian smoothing attack...   ");

	const int size = iterationList.dim();
	for(int iterationIndex = 0; iterationIndex < size; iterationIndex++) 
	{
		int iteration = iterationList[iterationIndex];
		
		Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

		attack.LaplacianSmoothing(deformFactor, iteration, preserveBoundaries);
	}

	printf("done!\n");
}

void AppManager::simplification(
	AppObject* appObject, 
	SimplificationType simpType, const TNT::Array1D<double>& simplificationRatioArray)
{
	printf("\n#Simplification attack...   ");

	const int size = simplificationRatioArray.dim();
	for(int simplificationIndex = 0; simplificationIndex < size; simplificationIndex++)
	{
		double ratio = simplificationRatioArray[simplificationIndex];

		Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

		attack.Simplification(simpType, ratio);
	}

	printf("done!\n");
}

void AppManager::subdivision(AppObject* appObject, SubdivisionType subdivisionType)
{
	printf("\n#Subdivision attack...   ");

	int subdivisionIteration = 1;
	//SubdivisionType subdivisionType = MIDPOINT;

	Attack<Enriched_polyhedron<Enriched_kernel, Enriched_items>, Enriched_kernel> attack(appObject->mesh, appObject->fileSuffix, appObject->fileName,  appObject->fileDir);

	attack.Subdivision(subdivisionType, subdivisionIteration);

	printf("done!\n");

	//attack.Subdivision(CATMULLCLARK, 1); // 1 means number of Subdivision iterations
	//attack.Subdivision(LOOP, 1); // 1 means number of Subdivision iterations 
	//attack.Subdivision(DOOSABIN, 1); // 1 means number of Subdivision iterations 
	//attack.Subdivision(SQRT3, 1); // 1 means number of Subdivision iterations 
}



double deformationDistance(double** a, double** b, int order)
{
	MyMatrix myMatrix;
	double** m = myMatrix.matrixSubtraction(a, b, order);
	double f = myMatrix.frobeniusNorm(m, order);

	// release memory
	for(int index = 0; index < order; index++)
		delete [] m[index];
	delete [] m;

	return  f;
}

//////////////////////////////////////////////////////////////////////////////////

#include <QVector>
void AppManager::dirRunnerAvgDDegree(AppObject* appObject, const QString& fileDir)
{
	// 1. calculate orientation matrix of the first frame
	// 2. for 2~n frames
	//     a. calculate orientation matrix of each facet
	//     b. calculate deformation distance of each facet
	//     c. use a facet array to record current max deformation distance of each facet
	// 3. paint this result

	const int order = 3;

	Mesh* sourceMesh = appObject->mesh;

	QString ReferenceFramePath = appObject->inFileName;

	cout << "Reference frame path: " << qPrintable(ReferenceFramePath) << endl;
	
	/// reference frame orientation matrix ///
	int numOfFacets = sourceMesh->size_of_facets();

	QVector<double**> sourceOrientationMatrixs(numOfFacets);

	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
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
		
		double** orientationMatrix = new double*[order];
		for(int i = 0; i < order; i++)
			orientationMatrix[i] = new double[order];

		double** inverseOrientatoinMatrix = new double*[order];
		for(int i = 0; i < order; i++)
			inverseOrientatoinMatrix[i] = new double[order];

		orientationMatrix[0][0] = v1.x();
		orientationMatrix[1][0] = v1.y();
		orientationMatrix[2][0] = v1.z();

		orientationMatrix[0][1] = v2.x();
		orientationMatrix[1][1] = v2.y();
		orientationMatrix[2][1] = v2.z();

		orientationMatrix[0][2] = v3.x();
		orientationMatrix[1][2] = v3.y();
		orientationMatrix[2][2] = v3.z();

		MyMatrix myMatrix;
		myMatrix.matrixInversion(orientationMatrix, order, inverseOrientatoinMatrix);

		int fit_index = fit->index();
		sourceOrientationMatrixs[fit_index] = inverseOrientatoinMatrix;

	    /// release memory ///
		for(int i = 0; i < order; i++)
			delete [] orientationMatrix[i];
		delete [] orientationMatrix;
	}

	printf("Reference frame done!\n");

	Mesh* mesh = new Mesh;

	QVector<double> maxDDegree(numOfFacets);
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
		maxDDegree[facetIndex] = FLT_MIN;

	/// list dir files ///
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for(QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;
		
		if(filePath == ReferenceFramePath) continue;

		cout << qPrintable(filePath) << ":  ";

		bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);

		if(!success) continue;

		// for each facet, calculate its deformation gradient matrix
		// deformation gradient matrix Dt,i = Ot,i * (Or,i)^-1
		// Or,i : reference frame orientation matrix
		// Ot,i : current frame orientation
		// * : matrix product
		QVector<double**> deformationGradientMatrixs(numOfFacets);

		for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
		{
			double** m = new double*[order];
			for(int j = 0; j < order; j++)
				m[j] = new double[order];
			deformationGradientMatrixs[facetIndex] = m;
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

			double** orientationMatrix = new double*[order];
			for(int i = 0; i < order; i++)
				orientationMatrix[i] = new double[order];

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
			double** referenceOrientatoinMatrix = sourceOrientationMatrixs.at(fit_index);

			MyMatrix myMatrix;
			double** deformationGradient = myMatrix.matrixMultiplication(orientationMatrix, referenceOrientatoinMatrix, order);

			deformationGradientMatrixs[fit_index] = deformationGradient;
		}
		
		QVector<double> deforDist(numOfFacets);
		for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
			deforDist[facetIndex] = 0;

		/// sum all the deformation degree of current facet ///
		for (Mesh::Edge_const_iterator eit = mesh->edges_begin(); eit != mesh->edges_end(); eit++)
		{
			if (eit->facet() == NULL || eit->opposite()->facet() == NULL)
				continue;

			Mesh::Facet_const_handle f1 = eit->facet();
			Mesh::Facet_const_handle f2 = eit->opposite()->facet();

			int f_index1 = f1->index();
			int f_index2 = f2->index();

			double** m1 = deformationGradientMatrixs.at(f_index1);
			double** m2 = deformationGradientMatrixs.at(f_index2);

			double dd = deformationDistance(m1, m2, order);

			deforDist[f_index1] += dd;
			deforDist[f_index2] += dd;
		}

		 /// avg deformation degree ///
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

		/// release memory ///
		for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
		{
			double** defor = deformationGradientMatrixs.at(facetIndex);
			for(int j = 0; j < order; j++)
				delete [] defor[j];
			delete [] defor;
		}

		mesh->clear();

		cout << "ok!" << endl;
	}
	
	/// release memory ///
	delete mesh;

	/// save the dd to reference frame mesh ///
	double maxDD = FLT_MIN;
	double minDD = FLT_MAX;
	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		double mdd = maxDDegree[f_index];

		fit->deformationDegree(mdd);

		if(mdd > maxDD) maxDD = mdd;
		if(mdd < minDD) minDD = mdd;
	}

	/// normalize ///
	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		double mdd = fit->deformationDegree();
		double n_mdd = (mdd - minDD) / (maxDD - minDD);
		
		fit->normalizeDeformationDegree(n_mdd);
	}

	/// release memory ///
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
	{
		double** ref = sourceOrientationMatrixs[facetIndex];
		for(int j = 0; j < order; j++)
			delete [] ref[j];
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

	Mesh* sourceMesh = appObject->mesh;

	QString ReferenceFramePath = appObject->inFileName;

	cout << "Reference frame path: " << qPrintable(ReferenceFramePath) << endl;

	/// reference frame orientation matrix ///
	int numOfFacets = sourceMesh->size_of_facets();

	QVector<double**> sourceOrientationMatrixs(numOfFacets);

	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
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
		
		double** orientationMatrix = new double*[order];
		for(int i = 0; i < order; i++)
			orientationMatrix[i] = new double[order];

		double** inverseOrientatoinMatrix = new double*[order];
		for(int i = 0; i < order; i++)
			inverseOrientatoinMatrix[i] = new double[order];

		orientationMatrix[0][0] = v1.x();
		orientationMatrix[1][0] = v1.y();
		orientationMatrix[2][0] = v1.z();

		orientationMatrix[0][1] = v2.x();
		orientationMatrix[1][1] = v2.y();
		orientationMatrix[2][1] = v2.z();

		orientationMatrix[0][2] = v3.x();
		orientationMatrix[1][2] = v3.y();
		orientationMatrix[2][2] = v3.z();

		MyMatrix myMatrix;
		myMatrix.matrixInversion(orientationMatrix, order, inverseOrientatoinMatrix);

		int fit_index = fit->index();
		sourceOrientationMatrixs[fit_index] = inverseOrientatoinMatrix;

	    /// release memory ///
		for(int orderIndex = 0; orderIndex < order; orderIndex++)
			delete [] orientationMatrix[orderIndex];
		delete [] orientationMatrix;
	}

	printf("Reference frame done!\n");

	Mesh* mesh = new Mesh;

	QVector<double> maxDDegree(numOfFacets);
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
		maxDDegree[facetIndex] = FLT_MIN;

	/// list dir files ///
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for(QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		if(filePath == ReferenceFramePath) continue;

		cout << qPrintable(filePath) << ":  ";

		bool success = FileParser::read_Obj(filePath.toAscii(), mesh, 1.0);

		if(!success) continue;

		// for each facet, calculate its deformation gradient matrix
		// deformation gradient matrix Dt,i = Ot,i * (Or,i)^-1
		// Or,i : reference frame orientation matrix
		// Ot,i : current frame orientation
		// * : matrix product
		QVector<double**> deformationGradientMatrixs(numOfFacets);

		for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
		{
			double** m = new double*[order];
			for(int j = 0; j < order; j++)
				m[j] = new double[order];
			deformationGradientMatrixs[facetIndex] = m;
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

			double** orientationMatrix = new double*[order];
			for(int i = 0; i < order; i++)
				orientationMatrix[i] = new double[order];

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
			double** referenceOrientatoinMatrix = sourceOrientationMatrixs.at(fit_index);

			MyMatrix myMatrix;
			double** deformationGradient = myMatrix.matrixMultiplication(orientationMatrix, referenceOrientatoinMatrix, order);

			deformationGradientMatrixs[fit_index] = deformationGradient;
		}

		/// find the largest deformation distance between two adjacent faces ///
		for(Mesh::Facet_iterator fit = mesh->facets_begin(); fit != mesh->facets_end(); fit++)
		{
			double maxDeforDist = FLT_MIN;

			int f_index1 = fit->index();

			Mesh::Halfedge_around_facet_circulator pHalfedge = fit->facet_begin();
			do { 
				if(pHalfedge->opposite()->facet() == NULL)
					continue;
						
				int f_index2 = pHalfedge->opposite()->facet()->index();

				double** m1 = deformationGradientMatrixs.at(f_index1);
				double** m2 = deformationGradientMatrixs.at(f_index2);

				double dd = deformationDistance(m1, m2, order);

				if( dd > maxDeforDist )
					maxDeforDist = dd;
			} while(++pHalfedge != fit->facet_begin());
			
			if(maxDeforDist > maxDDegree[f_index1])
				maxDDegree[f_index1] = maxDeforDist;
		}

		/// release memory ///
		for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
		{
			double** defor = deformationGradientMatrixs.at(facetIndex);
			for(int j = 0; j < order; j++)
				delete [] defor[j];
			delete [] defor;
		}

		mesh->clear();

		cout << "ok!" << endl;
	}
	
	/// release memory ///
	delete mesh;

	/// save the dd to reference frame mesh ///
	double maxDD = FLT_MIN;
	double minDD = FLT_MAX;
	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		double mdd = maxDDegree[f_index];

		fit->deformationDegree(mdd);

		if(mdd > maxDD) maxDD = mdd;
		if(mdd < minDD) minDD = mdd;
	}

	/// normalize ///
	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		double mdd = fit->deformationDegree();
		double n_mdd = (mdd - minDD) / (maxDD - minDD);
		
		fit->normalizeDeformationDegree(n_mdd);
	}

	/// release memory ///
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
	{
		double** ref = sourceOrientationMatrixs[facetIndex];
		for(int j = 0; j < order; j++)
			delete [] ref[j];
		delete [] ref;
	}
}


#include "tnt/jama_lu.h"
#include "tnt/tnt_array2d_utils.h"

double AppManager::calculateDeformationDistance(Mesh* originalMesh, Mesh* deformedMesh)
{
	// 1. calculate orientation matrix of the first frame
	// 2a. calculate orientation matrix of each facet
	// 2b. calculate deformation distance of each facet

	// for each facet, calculate its deformation gradient matrix
	// deformation gradient matrix Dt,i = Ot,i * (Or,i)^-1
	// Or,i : reference frame orientation matrix
	// Ot,i : current frame orientation
	// * : matrix product

	const int order = 3;
	const int numOfFacets = originalMesh->size_of_facets();
	MyMatrix myMatrix;

	QVector<double**> deformationGradientMatrixs(numOfFacets);
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
	{
		double** m = new double*[order];
		for(int j = 0; j < order; j++)
			m[j] = new double[order];
		deformationGradientMatrixs[facetIndex] = m;
	}

	double max = -FLT_MAX;
	Mesh::Facet_iterator original_fit = originalMesh->facets_begin();
	Mesh::Facet_iterator original_fit_end = originalMesh->facets_end();
	Mesh::Facet_iterator deformed_fit = deformedMesh->facets_begin();
	for(; original_fit != original_fit_end; original_fit++, deformed_fit++)
	{
		Mesh::Halfedge_around_facet_const_circulator original_circulator = original_fit->facet_begin();
		Mesh::Halfedge_around_facet_const_circulator deformed_circulator = deformed_fit->facet_begin();

		/// original orientation matrix ///
		Point_3 originalPoint1 = original_circulator->vertex()->point(); ++original_circulator;
		Point_3 originalPoint2 = original_circulator->vertex()->point(); ++original_circulator;
		Point_3 originalPoint3 = original_circulator->vertex()->point();

		Vector_3 originalVector1 = originalPoint2 - originalPoint1;
		Vector_3 originalVector2 = originalPoint3 - originalPoint1;
		Vector_3 originalVector3 = cross_product(originalVector1, originalVector2);
		originalVector3 = originalVector3 / sqrt(originalVector3.squared_length());
		
		double** originalOrientationMatrix = new double*[order];
		double** inverseOriginalOrientatoinMatrix = new double*[order];
		double** deformedOrientationMatrix = new double*[order];
		for(int orderIndex = 0; orderIndex < order; orderIndex++)
		{
			originalOrientationMatrix[orderIndex] = new double[order];
			inverseOriginalOrientatoinMatrix[orderIndex] = new double[order];
			deformedOrientationMatrix[orderIndex] = new double[order];
		}

		originalOrientationMatrix[0][0] = originalVector1.x();
		originalOrientationMatrix[1][0] = originalVector1.y();
		originalOrientationMatrix[2][0] = originalVector1.z();
		originalOrientationMatrix[0][1] = originalVector2.x();
		originalOrientationMatrix[1][1] = originalVector2.y();
		originalOrientationMatrix[2][1] = originalVector2.z();
		originalOrientationMatrix[0][2] = originalVector3.x();
		originalOrientationMatrix[1][2] = originalVector3.y();
		originalOrientationMatrix[2][2] = originalVector3.z();

		myMatrix.matrixInversion(originalOrientationMatrix, order, inverseOriginalOrientatoinMatrix);

		double inverseDD = myMatrix.frobeniusNorm(inverseOriginalOrientatoinMatrix, order);
		cout << "inverseDD = " << inverseDD << endl;
		///////////////////////////////////

		/// deformed orientation matrix ///
		Point_3 deformedPoint1 = deformed_circulator->vertex()->point(); ++deformed_circulator;
		Point_3 deformedPoint2 = deformed_circulator->vertex()->point(); ++deformed_circulator;
		Point_3 deformedPoint3 = deformed_circulator->vertex()->point();

		Vector_3 deformedVector1 = deformedPoint2 - deformedPoint1;
		Vector_3 deformedVector2 = deformedPoint3 - deformedPoint1;
		Vector_3 deformedVector3 = cross_product(deformedVector1, deformedVector2);
		deformedVector3 = deformedVector3 / sqrt(deformedVector3.squared_length());
		
		deformedOrientationMatrix[0][0] = deformedVector1.x();
		deformedOrientationMatrix[1][0] = deformedVector1.y();
		deformedOrientationMatrix[2][0] = deformedVector1.z();
		deformedOrientationMatrix[0][1] = deformedVector2.x();
		deformedOrientationMatrix[1][1] = deformedVector2.y();
		deformedOrientationMatrix[2][1] = deformedVector2.z();
		deformedOrientationMatrix[0][2] = deformedVector3.x();
		deformedOrientationMatrix[1][2] = deformedVector3.y();
		deformedOrientationMatrix[2][2] = deformedVector3.z();
		///////////////////////////////////

		double** deformationGradientMatrix = myMatrix.matrixMultiplication(deformedOrientationMatrix, inverseOriginalOrientatoinMatrix, order);
		deformationGradientMatrixs[original_fit->index()] = deformationGradientMatrix;

		double maxDD = myMatrix.frobeniusNorm(deformationGradientMatrix, order);
		if (maxDD > max)
			max = maxDD;

		/// release memory ///
		for(int orderIndex = 0; orderIndex < order; orderIndex++)
		{
			delete [] originalOrientationMatrix[orderIndex];
			delete [] inverseOriginalOrientatoinMatrix[orderIndex];
			delete [] deformedOrientationMatrix[orderIndex];
		}
		delete [] originalOrientationMatrix;
		delete [] inverseOriginalOrientatoinMatrix;
		delete [] deformedOrientationMatrix;
		//////////////////////
	}

	cout << "max = " << max << endl;

	QVector<double> deforDist(numOfFacets);
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
		deforDist[facetIndex] = 0;

	/// sum all the deformation degree of current facet ///
	for (Mesh::Edge_const_iterator eit = deformedMesh->edges_begin(); eit != deformedMesh->edges_end(); eit++)
	{
		if (eit->facet() == NULL || eit->opposite()->facet() == NULL)
			continue;

		Mesh::Facet_const_handle f1 = eit->facet();
		Mesh::Facet_const_handle f2 = eit->opposite()->facet();

		int f_index1 = f1->index();
		int f_index2 = f2->index();

		double** m1 = deformationGradientMatrixs.at(f_index1);
		double** m2 = deformationGradientMatrixs.at(f_index2);

		/// deformation distance ///
		double dd = deformationDistance(m1, m2, order);
		//cout << dd << endl;
		////////////////////////////

		deforDist[f_index1] += dd;
		deforDist[f_index2] += dd;
	}


	/// avg deformation degree ///
	double maxDeformationDistance = 0.0;
	double minDeformationDistance = FLT_MAX;
	for(Mesh::Facet_iterator fit = deformedMesh->facets_begin(); fit != deformedMesh->facets_end(); fit++)
	{
		int totalNeighboringFacets = 0;
		int f_index = fit->index();

		Mesh::Halfedge_around_facet_const_circulator c = fit->facet_begin();
		do { 
			if (c->facet() == NULL || c->opposite()->facet() == NULL)
				continue;
			totalNeighboringFacets++; 
		} while (++c != fit->facet_begin());

		deforDist[f_index] /= totalNeighboringFacets;

		if (deforDist[f_index] > maxDeformationDistance)
			maxDeformationDistance = deforDist[f_index];
		if (deforDist[f_index] < minDeformationDistance)
			minDeformationDistance = deforDist[f_index];
	}

	
	cout << "maxDeformationDistance = " << maxDeformationDistance << endl;
	cout << "minDeformationDistance = " << minDeformationDistance << endl;

	/// release memory ///
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
	{
		double** defor = deformationGradientMatrixs.at(facetIndex);
		for(int j = 0; j < order; j++)
			delete [] defor[j];
		delete [] defor;
	}

	return maxDeformationDistance;
	/////////////////////
}


void AppManager::calculateMaxDeformationDistance(AppObject* appObject, Mesh* sourceMesh, Mesh* targetMesh)
{
	// 1. calculate orientation matrix of the first frame
	// 2. for 2~n frames
	//     a. calculate orientation matrix of each facet
	//     b. calculate deformation distance of each facet
	//     c. use a facet array to record current max deformation distance of each facet
	// 3. paint this result

	const int order = 3;

	/// reference frame orientation matrix ///
	int numOfFacets = sourceMesh->size_of_facets();

	QVector<double**> sourceOrientationMatrixs(numOfFacets);

	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
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
		
		double** orientationMatrix = new double*[order];
		for(int i = 0; i < order; i++)
			orientationMatrix[i] = new double[order];

		double** inverseOrientatoinMatrix = new double*[order];
		for(int i = 0; i < order; i++)
			inverseOrientatoinMatrix[i] = new double[order];

		orientationMatrix[0][0] = v1.x();
		orientationMatrix[1][0] = v1.y();
		orientationMatrix[2][0] = v1.z();
		orientationMatrix[0][1] = v2.x();
		orientationMatrix[1][1] = v2.y();
		orientationMatrix[2][1] = v2.z();
		orientationMatrix[0][2] = v3.x();
		orientationMatrix[1][2] = v3.y();
		orientationMatrix[2][2] = v3.z();

		MyMatrix myMatrix;
		myMatrix.matrixInversion(orientationMatrix, order, inverseOrientatoinMatrix);

		int fit_index = fit->index();
		sourceOrientationMatrixs[fit_index] = inverseOrientatoinMatrix;

	    /// release memory ///
		for(int orderIndex = 0; orderIndex < order; orderIndex++)
			delete [] orientationMatrix[orderIndex];
		delete [] orientationMatrix;
	}

	printf("Reference frame done!\n");

	QVector<double> maxDDegree(numOfFacets);
	for(int i = 0; i < numOfFacets; i++)
		maxDDegree[i] = FLT_MIN;

	// for each facet, calculate its deformation gradient matrix
	// deformation gradient matrix Dt,i = Ot,i * (Or,i)^-1
	// Or,i : reference frame orientation matrix
	// Ot,i : current frame orientation
	// * : matrix product
	QVector<double**> deformationGradientMatrixs(numOfFacets);
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
	{
		double** m = new double*[order];
		for(int j = 0; j < order; j++)
			m[j] = new double[order];
		deformationGradientMatrixs[facetIndex] = m;
	}

	for(Mesh::Facet_const_iterator fit = targetMesh->facets_begin(); fit != targetMesh->facets_end(); fit++)
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

		double** orientationMatrix = new double*[order];
		for(int i = 0; i < order; i++)
			orientationMatrix[i] = new double[order];

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
		double** referenceOrientatoinMatrix = sourceOrientationMatrixs.at(fit_index);

		MyMatrix myMatrix;
		double** deformationGradient = myMatrix.matrixMultiplication(orientationMatrix, referenceOrientatoinMatrix, order);

		deformationGradientMatrixs[fit_index] = deformationGradient;
	}

	/// find the largest deformation distance between two adjacent faces ///
	for(Mesh::Facet_iterator fit = targetMesh->facets_begin(); fit != targetMesh->facets_end(); fit++)
	{
		double maxDeforDist = FLT_MIN;

		int f_index1 = fit->index();

		Mesh::Halfedge_around_facet_circulator pHalfedge = fit->facet_begin();
		do { 
			if(pHalfedge->opposite()->facet() == NULL)
				continue;
					
			int f_index2 = pHalfedge->opposite()->facet()->index();

			double** m1 = deformationGradientMatrixs.at(f_index1);
			double** m2 = deformationGradientMatrixs.at(f_index2);

			double dd = deformationDistance(m1, m2, order);

			if( dd > maxDeforDist )
				maxDeforDist = dd;
		} while(++pHalfedge != fit->facet_begin());
		
		if(maxDeforDist > maxDDegree[f_index1])
			maxDDegree[f_index1] = maxDeforDist;
	}

	/// release memory ///
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
	{
		double** defor = deformationGradientMatrixs.at(facetIndex);
		for(int j = 0; j < order; j++)
			delete [] defor[j];
		delete [] defor;
	}

	/// save the dd to reference frame mesh ///
	double maxDD = FLT_MIN;
	double minDD = FLT_MAX;
	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		double mdd = maxDDegree[f_index];

		fit->deformationDegree(mdd);

		if(mdd > maxDD) maxDD = mdd;
		if(mdd < minDD) minDD = mdd;
	}

	/// normalize ///
	for(Mesh::Facet_iterator fit = sourceMesh->facets_begin(); fit != sourceMesh->facets_end(); fit++)
	{
		int f_index = fit->index();
		double mdd = fit->deformationDegree();
		double n_mdd = (mdd - minDD) / (maxDD - minDD);
		
		fit->normalizeDeformationDegree(n_mdd);
	}

	/// release memory ///
	for(int facetIndex = 0; facetIndex < numOfFacets; facetIndex++)
	{
		double** ref = sourceOrientationMatrixs[facetIndex];
		for(int j = 0; j < order; j++)
			delete [] ref[j];
		delete [] ref;
	}

}


void AppManager::OnDirRunnerCalculateDeformationDistance()
{
	QString fileDir = m_appParameters.dirRunner;
	QString originalMeshFileName = m_appParameters.fileName1;

	cout << "directory => " << qPrintable(fileDir) << endl;
	cout << "originalMeshFileName => " << qPrintable(originalMeshFileName) << endl;

	QString deformationDistanceFileDir = fileDir + "\\" + "deforDists";

	QDir d(deformationDistanceFileDir);
	if(!d.exists())
		d.mkdir(deformationDistanceFileDir);

	Mesh* originalMesh = new Mesh;
	if (!ImporterMesh::loadMesh(originalMeshFileName, originalMesh))
		return;

	QDir fd(fileDir);
	fd.setFilter(QDir::Files);
	const QFileInfoList list = fd.entryInfoList();
	const int listSize = list.size();

	TNT::Array1D<double> maxDeformationDistanceArray = TNT::Array1D<double>(listSize);
	int listIndex = 0;
	for (QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString deformedMeshFileName = fileDir + "\\" + inFileName;

		Mesh* deformedMesh = new Mesh;
		if (!ImporterMesh::loadMesh(deformedMeshFileName, deformedMesh))
			return;

		double maxDeformationDistance = calculateDeformationDistance(originalMesh, deformedMesh);

		maxDeformationDistanceArray[listIndex] = maxDeformationDistance;

		delete deformedMesh;

		listIndex++;
	}

	QString exportFileName = deformationDistanceFileDir + "\\" + "deformationDistance.txt";

	ExporterArray exporterArray;
	exporterArray.exportArray1D(maxDeformationDistanceArray, exportFileName);
}


void AppManager::OnDirRunnerAvgDDegree()
{
	cout << "OnDirRunnerAvgDDegree" << endl;
	/*
	QString fileDir = "C:/Users/hsiao/Desktop/dance-clone";
	
	QDir d(fileDir);
	if(!d.exists()) { cout << "\nX " << qPrintable(fileDir) << ": not exists!" << endl; return; }

	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	printf("\n#Avg Deformation Degree...\n");

	dirRunnerAvgDDegree(appObject, fileDir);

	printf("#Done!\n");
	*/
}

void AppManager::OnDirRunnerMaxDDegree()
{
	cout << "OnDirRunnerMaxDDegree" << endl;
	/*
	QString fileDir = "C:/Users/hsiao/Desktop/dance-clone";
	
	QDir d(fileDir);
	if(!d.exists()) { cout << "\nX " << qPrintable(fileDir) << ": not exists!" << endl; return; }

	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	printf("\n#Max Deformation Degree...\n");

	dirRunnerMaxDDegree(appObject, fileDir);

	printf("#Done!\n");
	*/
}


void AppManager::OnPlayAnimation()
{
	QString fileDir = m_appParameters.dirRunner;

	cout << "directory => " << qPrintable(fileDir) << endl;

	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for (QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		loadObject(filePath, 1.0);

		AppObject* appObject = getSelectedObject();
		if (!appObject) continue;

		removeObject(appObject->id);
	}
}


void AppManager::OnDirRunnerCalculateLSDVertices()
{
	QString fileDir = m_appParameters.dirRunner;
	QString sdfsFileDir = fileDir + "\\" + "sdfs";

	QDir d(sdfsFileDir);
	if(!d.exists())
		d.mkdir(sdfsFileDir);
		
	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	ExporterMeshFeature exporterMeshFeature;
	const QFileInfoList list = fd.entryInfoList();
	for (QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		loadObject(filePath, 1.0);

		try
		{
			AppObject* appObject = getSelectedObject();
			if (!appObject) continue;

			bool success = calculateShapeDiameterFunction(appObject, true);
			if (!success)
				printf("Failed on calculate SDF Vertices\n");

			QString exportVertexSDFFileName = sdfsFileDir + "\\" + appObject->fileName + "_vsdf.txt";
			QString exportVertexNSDFFileName = sdfsFileDir + "\\" + appObject->fileName + "_vnsdf.txt";

			exporterMeshFeature.saveSDFText(exportVertexSDFFileName, true, appObject);
			exporterMeshFeature.saveNSDFText(exportVertexNSDFFileName, true, appObject);

			removeObject(appObject->id);
		}
		catch (...)
		{
			QMessageBox::critical(NULL, "SDF Vertices", "Exception");
		}
	}
}

void AppManager::OnDirRunnerCalculateLSDFacets()
{
	QString fileDir = m_appParameters.dirRunner;
	QString sdfsFileDir = fileDir + "\\" + "sdfs";

	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	ExporterMeshFeature exporterMeshFeature;
	const QFileInfoList list = fd.entryInfoList();
	for (QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		loadObject(filePath, 1.0);

		try
		{
			AppObject* appObject = getSelectedObject();
			if (!appObject) continue;

			bool success = calculateShapeDiameterFunction(appObject, true);
			if (!success)
				printf("Failed on calculate SDF Vertices\n");

			QString exportFacetSDFFileName = sdfsFileDir + "\\" + appObject->fileName + "_fsdf.txt";
			QString exportFacetNSDFFileName = sdfsFileDir + "\\" + appObject->fileName + "_fnsdf.txt";

			exporterMeshFeature.saveSDFText(exportFacetSDFFileName, false, appObject);
			exporterMeshFeature.saveNSDFText(exportFacetNSDFFileName, false, appObject);

			removeObject(appObject->id);
		}
		catch (...)
		{
			QMessageBox::critical(NULL, "SDF Facets", "Exception");
		}
	}
}

void AppManager::OnDirRunnerDetectWatermark()
{
	QString fileDir = m_appParameters.dirRunner;

	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for (QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString fileName = (*iterator).baseName();

		QString filePath = fileDir + "\\" + inFileName;
		//QString vsdfFeatureFileName = fileDir + "\\" + "sdfs" + "\\" + fileName + "_vsdf.txt";
		//QString vnsdfFeatureFileName = fileDir + "\\" + "sdfs" + "\\" + fileName + "_vnsdf.txt";

		loadObject(filePath, 1.0);
	
		AppObject* appObject = getSelectedObject();
		if (!appObject) continue;
		
		//if (QFile::exists(vsdfFeatureFileName))
		//	ImporterMeshFeatures::loadSDFText(vsdfFeatureFileName, true, appObject);

		//if (QFile::exists(vnsdfFeatureFileName))
		//	ImporterMeshFeatures::loadNSDFText(vnsdfFeatureFileName, true, appObject);

		if (!calculateShapeDiameterFunction(appObject, true))
			cout << "Failed on calculate SDF Vertices\n";

		if (!detectWatermark(appObject, true))
			cout << "Failed on detect watermark\n";

		removeObject(appObject->id);
	}
}

void AppManager::OnDirRunnerAllAttacks()
{
	QString fileDir = m_appParameters.dirRunner;

	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for (QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		loadObject(filePath, 1.0);
	
		AppObject* appObject = getSelectedObject();
		if (!appObject) continue;

		allAttacks(appObject);

		removeObject(appObject->id);
	}
}

void AppManager::OnDirRunnerReordering()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	reordering(appObject);
}

void AppManager::OnDirRunnerNoiseAttack()
{
	AppObject* appObject = getSelectedObject();
	if(!appObject) return;

	const bool preserveBoundaries = true;

	const int size = 5;
	TNT::Array1D<double> noiseIntensityArray = TNT::Array1D<double>(size);
	noiseIntensityArray[0] = 0.0005;
	noiseIntensityArray[1] = 0.001;
	noiseIntensityArray[2] = 0.002;
	noiseIntensityArray[3] = 0.003;
	noiseIntensityArray[4] = 0.004;	

	noiseAddition(appObject, preserveBoundaries, noiseIntensityArray);
}

void AppManager::OnDirRunnerQuantization()
{
	AppObject* appObject = getSelectedObject();
	if(!appObject) return;

	const int size = 5;

	TNT::Array1D<int> quantizationBitArray = TNT::Array1D<int>(size);
	quantizationBitArray[0] = 7;
	quantizationBitArray[1] = 8;
	quantizationBitArray[2] = 9;
	quantizationBitArray[3] = 10;
	quantizationBitArray[4] = 11;

	quantization(appObject, quantizationBitArray);
}

void AppManager::OnDirRunnerSimilarityTransform()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	similarityTransform(appObject);
}

void AppManager::OnDirRunnerSmoothing()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	const double deformFactor = 0.1;
	const bool preserveBoundaries = true;

	const int size = 5;
	TNT::Array1D<int> iterationList = TNT::Array1D<int>(size);
	iterationList[0] = 5;
	iterationList[1] = 10;
	iterationList[2] = 20;
	iterationList[3] = 30;
	iterationList[4] = 40;
	
	smoothing(appObject, deformFactor, preserveBoundaries, iterationList);
}

void AppManager::OnDirRunnerSimplificationLindstromturk()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SimplificationType simpType = LINDSTROMTURK;

	const int size = 5;
	TNT::Array1D<double> simplificationRatioArray = TNT::Array1D<double>(size);
	simplificationRatioArray[0] = 5;
	simplificationRatioArray[1] = 10;
	simplificationRatioArray[2] = 20;
	simplificationRatioArray[3] = 30;
	simplificationRatioArray[4] = 40;

	simplification(appObject, simpType, simplificationRatioArray);
}

void AppManager::OnDirRunnerSimplificationMidpoint()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	SimplificationType simpType = EDGELENGTHMIDPOINT;

	const int size = 5;
	TNT::Array1D<double> simplificationRatioArray = TNT::Array1D<double>(size);
	simplificationRatioArray[0] = 5;
	simplificationRatioArray[1] = 10;
	simplificationRatioArray[2] = 20;
	simplificationRatioArray[3] = 30;
	simplificationRatioArray[4] = 40;

	simplification(appObject, simpType, simplificationRatioArray);
}

void AppManager::OnDirRunnerSubdivisionCatmullClark()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = CATMULLCLARK;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnDirRunnerSubdivisionDoosabin()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = DOOSABIN;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnDirRunnerSubdivisionLoop()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = LOOP;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnDirRunnerSubdivisionMidpoint()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = MIDPOINT;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnDirRunnerSubdivisionSqrt3()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = SQRT3;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnAllAttacks()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	allAttacks(appObject);
}

void AppManager::OnReordering()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	reordering(appObject);
}

void AppManager::OnNoiseAddition()
{
	AppObject* appObject = getSelectedObject();
	if(!appObject) return;

	const bool preserveBoundaries = true;

	const int size = 5;
	TNT::Array1D<double> noiseIntensityArray = TNT::Array1D<double>(size);
	noiseIntensityArray[0] = 0.0005;
	noiseIntensityArray[1] = 0.001;
	noiseIntensityArray[2] = 0.002;
	noiseIntensityArray[3] = 0.003;
	noiseIntensityArray[4] = 0.004;	

	noiseAddition(appObject, preserveBoundaries, noiseIntensityArray);
}

void AppManager::OnQuantization()
{
	AppObject* appObject = getSelectedObject();
	if(!appObject) return;

	const int size = 5;

	TNT::Array1D<int> quantizationBitArray = TNT::Array1D<int>(size);
	quantizationBitArray[0] = 7;
	quantizationBitArray[1] = 8;
	quantizationBitArray[2] = 9;
	quantizationBitArray[3] = 10;
	quantizationBitArray[4] = 11;

	quantization(appObject, quantizationBitArray);
}

void AppManager::OnSimilarityTransform()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	similarityTransform(appObject);
}

void AppManager::OnTranslation()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	translation(appObject);
}

void AppManager::OnRotation()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	rotation(appObject);
}

void AppManager::OnUniformScaling()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	uniformScaling(appObject);
}

void AppManager::OnNonUniformScaling()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	nonuniformScaling(appObject);
}

void AppManager::OnSmoothing()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	const double deformFactor = 0.1;
	const bool preserveBoundaries = true;

	const int size = 5;
	TNT::Array1D<int> iterationList = TNT::Array1D<int>(size);
	iterationList[0] = 5;
	iterationList[1] = 10;
	iterationList[2] = 20;
	iterationList[3] = 30;
	iterationList[4] = 40;

	smoothing(appObject, deformFactor, preserveBoundaries, iterationList);
}

void AppManager::OnSimplificationLindstromturk()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SimplificationType simpType = LINDSTROMTURK;

	const int size = 5;
	TNT::Array1D<double> simplificationRatioArray = TNT::Array1D<double>(size);
	simplificationRatioArray[0] = 5;
	simplificationRatioArray[1] = 10;
	simplificationRatioArray[2] = 20;
	simplificationRatioArray[3] = 30;
	simplificationRatioArray[4] = 40;

	simplification(appObject, simpType, simplificationRatioArray);
}

void AppManager::OnSimplificationMidpoint()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	SimplificationType simpType = EDGELENGTHMIDPOINT;
	
	const int size = 5;
	TNT::Array1D<double> simplificationRatioArray = TNT::Array1D<double>(size);
	simplificationRatioArray[0] = 5;
	simplificationRatioArray[1] = 10;
	simplificationRatioArray[2] = 20;
	simplificationRatioArray[3] = 30;
	simplificationRatioArray[4] = 40;

	simplification(appObject, simpType, simplificationRatioArray);
}

void AppManager::OnSubdivisionCatmullClark()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = CATMULLCLARK;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnSubdivisionDoosabin()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = DOOSABIN;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnSubdivisionLoop()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = LOOP;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnSubdivisionMidpoint()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = MIDPOINT;
	subdivision(appObject, subdivisionType);
}

void AppManager::OnSubdivisionSqrt3()
{
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

	SubdivisionType subdivisionType = SQRT3;
	subdivision(appObject, subdivisionType);
}

//////////////////////////////////////////////////////////////////////////
// Segmentation
//////////////////////////////////////////////////////////////////////////
void AppManager::OnPartitionFacetsMesh()
{
    /*
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

    if(!partitionFacetsMesh(appObject))
		printf("Failed partition mesh");
	*/
}

void AppManager::OnPartitionVerticesMesh()
{
	/*
	AppObject* appObject = getSelectedObject();
    if(!appObject) return;

    if(!partitionVerticesMesh(appObject))
		printf("Failed partition mesh");
	*/
}

//////////////////////////////////////////////////////////////////////////
// Draw Histogram
//////////////////////////////////////////////////////////////////////////
void AppManager::OnDrawFacetLSDValue()
{
    AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFValue(appObject, false, false, m_appParameters.watermark_bits);
}

void AppManager::OnDrawFacetNLSDValue()
{
    AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFValue(appObject, false, true, m_appParameters.watermark_bits);
}

void AppManager::OnDrawVertexLSDValue()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFValue(appObject, true, false, m_appParameters.watermark_bits);
}

void AppManager::OnDrawVertexNLSDValue()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFValue(appObject, true, true, m_appParameters.watermark_bits);
}

void AppManager::OnDrawFacetLSDBins()
{
    AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFBins(appObject, false, false, m_appParameters.watermark_bits, m_appParameters.lsd_log_alpha);
}

void AppManager::OnDrawFacetNLSDBins()
{
    AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFBins(
		appObject, false, true, m_appParameters.watermark_bits, m_appParameters.lsd_log_alpha);
}

void AppManager::OnDrawVertexLSDBins()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFBins(
		appObject, true, false, m_appParameters.watermark_bits, m_appParameters.lsd_log_alpha);
}

void AppManager::OnDrawVertexNLSDBins()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFBins(
		appObject, true, true, m_appParameters.watermark_bits, m_appParameters.lsd_log_alpha);
}

void AppManager::OnDrawFacetLSDProbability()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFProbability(
		appObject, 
		false, 
		false,
		m_appParameters.watermark_bits);
}

void AppManager::OnDrawFacetNLSDProbability()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFProbability(
		appObject, 
		false, 
		true,
		m_appParameters.watermark_bits);
}

void AppManager::OnDrawVertexLSDProbability()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFProbability(
		appObject, 
		true, 
		false,
		m_appParameters.watermark_bits);
}

void AppManager::OnDrawVertexNLSDProbability()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	DrawHistogram dh;
	dh.drawSDFProbability(
		appObject, 
		true, 
		true,
		m_appParameters.watermark_bits);
}

//////////////////////////////////////////////////////////////////////////
// Watermark
//////////////////////////////////////////////////////////////////////////
void AppManager::OnEmbedWatermark()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	if (!embedWatermark(appObject, true))
		cout << "Failed on embed watermark\n";
}

void AppManager::OnDetectWatermark()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	if (!detectWatermark(appObject, true))
		cout << "Failed on detect watermark\n";
}

void AppManager::OnPoseTransfer()
{
	QString fileDir = m_appParameters.dirRunner;
	QString originalMeshFileName1 = m_appParameters.fileName1;
	QString targetMeshFileName1 = m_appParameters.fileName2;

	QDir fd(fileDir);
	fd.setFilter(QDir::Files);

	const QFileInfoList list = fd.entryInfoList();
	for (QFileInfoList::const_iterator iterator = list.begin(); iterator != list.end(); iterator++)
	{
		QString meshFormat = (*iterator).suffix();
		if (meshFormat != "obj") continue;

		QString inFileName = (*iterator).fileName();
		QString filePath = fileDir + "\\" + inFileName;

		if (filePath == originalMeshFileName1)
			continue;

		Mesh* originalMesh1 = new Mesh;
		Mesh* originalMesh2 = new Mesh;
		//Mesh* targetMesh1 = new Mesh;

		cout << qPrintable(filePath) << endl;
		
		QString originalMeshFileName2 = filePath;

		if (!ImporterMesh::loadMesh(originalMeshFileName1, originalMesh1))
			return;
		if (!ImporterMesh::loadMesh(originalMeshFileName2, originalMesh2))
			return;
		//if (!ImporterMesh::loadMesh(targetMeshFileName1, targetMesh1))
		//	return;

		loadObject(targetMeshFileName1, 1.0);

		AppObject* targetAppObject = getSelectedObject();
		Mesh* targetMesh1 = targetAppObject->mesh;

		QString targetMeshFileName2 = targetAppObject->fileDir + "\\" + "Watermark_SDF_Relation_" + inFileName;

		Mesh::Vertex_iterator tar_mesh1_vit = targetMesh1->vertices_begin();
		Mesh::Vertex_iterator ori_mesh2_vit = originalMesh2->vertices_begin();
		for (Mesh::Vertex_iterator ori_mesh1_vit = originalMesh1->vertices_begin(); ori_mesh1_vit != originalMesh1->vertices_end(); ori_mesh1_vit++, ori_mesh2_vit++, tar_mesh1_vit++)
		{
			Point_3 op1 = ori_mesh1_vit->point();
			Point_3 op2 = ori_mesh2_vit->point();
			
			double diff_x = op2.x() - op1.x();
			double diff_y = op2.y() - op1.y();
			double diff_z = op2.z() - op1.z();

			Point_3 tp1 = tar_mesh1_vit->point();
			Point_3 tp2(tp1.x() + diff_x, tp1.y() + diff_y, tp1.z() + diff_z);

			tar_mesh1_vit->point() = tp2;
		}

		originalMesh2->write_obj(targetMeshFileName2);
		
		delete originalMesh1;
		delete originalMesh2;
		//delete targetMesh1;
		removeObject(targetAppObject->id);
	}
}


void AppManager::OnCalculateDeformationDistance()
{
	QString originalMeshFileName = m_appParameters.fileName1;
	QString deformedMeshFileName = m_appParameters.fileName2;

	cout << "originalMeshFileName => " << qPrintable(originalMeshFileName) << endl;
	cout << "deformedMeshFileName => " << qPrintable(deformedMeshFileName) << endl;

	Mesh* originalMesh = new Mesh;
	Mesh* deformedMesh = new Mesh;

	if (!ImporterMesh::loadMesh(originalMeshFileName, originalMesh))
		return;
	if (!ImporterMesh::loadMesh(deformedMeshFileName, deformedMesh))
		return;

	calculateDeformationDistance(originalMesh, deformedMesh);

	delete originalMesh;
	delete deformedMesh;
}


/*bool AppManager::saveSDFVertices(const QString& featuresFileName, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_vertices();
	double* sdf = new double[sdf_size];
	Mesh::Vertex_const_iterator it = appObject->mesh->vertices_begin();
	Mesh::Vertex_const_iterator it_end = appObject->mesh->vertices_end();
	int i=0;
	for (;it != it_end; it++) {
		sdf[i] = it->volumeSDF();
		i++;
	}

	bool ret = FileUtils::save(
		featuresFileName,
		FEATURES_FILE_SDF_VERTICES,
		"SDF Vertices on " + appObject->inFileName,
		(char*)sdf,
		sizeof(double)*sdf_size);

	delete[] sdf;
	return ret;
}*/

/*bool AppManager::saveNSDFVertices(const QString& featuresFileName, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_vertices();
	double* sdf = new double[sdf_size];
	Mesh::Vertex_const_iterator it = appObject->mesh->vertices_begin();
	Mesh::Vertex_const_iterator it_end = appObject->mesh->vertices_end();
	int i=0;
	for (;it != it_end; it++) {
		sdf[i] = it->volumeNSDF();
		i++;
	}

	bool ret = FileUtils::save(
		featuresFileName,
		FEATURES_FILE_SDF_VERTICES,
		"NSDF Vertices on " + appObject->inFileName,
		(char*)sdf,
		sizeof(double)*sdf_size);

	delete[] sdf;
	return ret;
}*/

/*bool AppManager::saveSDFFacets(const QString& featuresFileName, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_facets();
	double* sdf = new double[sdf_size];
	Mesh::Facet_const_iterator it = appObject->mesh->facets_begin();
	Mesh::Facet_const_iterator it_end = appObject->mesh->facets_end();
	int i=0;
	for (;it != it_end; it++) {
		sdf[i] = it->volumeSDF();
		i++;
	}

	bool ret = FileUtils::save(
		featuresFileName,
		FEATURES_FILE_SDF_FACETS,
		"SDF Facets on " + appObject->inFileName,
		(char*)sdf,
		sizeof(double)*sdf_size);

	delete[] sdf;
	return ret;
}*/

/*bool AppManager::saveNSDFFacets(const QString& featuresFileName, AppObject* appObject)
{
	unsigned int sdf_size = appObject->mesh->size_of_facets();
	double* sdf = new double[sdf_size];
	Mesh::Facet_const_iterator it = appObject->mesh->facets_begin();
	Mesh::Facet_const_iterator it_end = appObject->mesh->facets_end();
	int i=0;
	for (;it != it_end; it++) {
		sdf[i] = it->volumeNSDF();
		i++;
	}

	bool ret = FileUtils::save(
		featuresFileName,
		FEATURES_FILE_SDF_FACETS,
		"NSDF Facets on " + appObject->inFileName,
		(char*)sdf,
		sizeof(double)*sdf_size);

	delete[] sdf;
	return ret;
}*/


void AppManager::saveSnapshot(const QString &fileName)
{
	g_sdfmain->ui.world->setSnapshotFormat("PNG");
	g_sdfmain->ui.world->saveSnapshot(fileName, true);
}

void AppManager::reverseObjectNormals(AppObject* appObject)
{
	if (!appObject || !appObject->mesh) return;

	appObject->mesh->reverseNormals();
	appObject->mesh->inside_out();

	//remember that I reversed normals for this object
	ProgSettings ps;
	ps.beginGroup("Reverse Normal Models");

	bool alreadyReversed = ps.readBoolEntry(appObject->inFileName);
	ps.writeEntry(appObject->inFileName, !alreadyReversed);
	ps.endGroup();

	emit ObjectChanged(appObject->id, true);
}

void AppManager::OnReverseNormals()
{
	reverseObjectNormals(getSelectedObject());
}

void AppManager::OnSelectedObjectInfo()
{
	AppObject* appObject = getSelectedObject();
	if (!appObject) return;

	QString info;

	QString generalInfo;
	generalInfo.sprintf("Object %s. #vertices = %d #facets = %d\nSelected facet #%d",
		appObject->inFileName.toAscii(),
		appObject->mesh->size_of_vertices(),
		appObject->mesh->size_of_facets(),
		appObject->meshSelection->selectedFacet());
	info.append(generalInfo);

	QMessageBox::information((QWidget*)parent(), appObject->id + " info", info);
}

void AppManager::OnObjectSelected(const QString& id)
{
	selectObject(id);
}

void AppManager::OnArrowKeyPressed(int key, Qt::ButtonState s)
{
	if (key != Qt::Key_Up)
		return;
}

void AppManager::OnFacetSelected(const QString& id, const int facetIndex)
{
	selectObject(id);
}

void AppManager::OnVertexSelected(const QString& id, const int vertexIndex)
{
	selectObject(id);
}

QStringList AppManager::getObjectNames(int& selectedIndex)
{
	QStringList list;

	selectedIndex = -1;

	appObjectsMap_t::iterator it = m_objects.begin();
	appObjectsMap_t::iterator it_end = m_objects.end();
	int i=0;
	for (;it != it_end; it++) {
		AppObject* appObject = it.value();
		if (appObject->meshSelection->isSelected())
			selectedIndex = i;
		list.append(appObject->id);
		i++;
	}

	return list;
}

AppObject* AppManager::getObject(const QString& id)
{
	appObjectsMap_t::iterator it = m_objects.find(id.ascii());
	if (it == m_objects.end())
		return NULL;
	else
		return it.value();
}