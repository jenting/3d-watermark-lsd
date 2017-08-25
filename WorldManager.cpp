#include "stdafx.h"

#define min min
#define max max
#define _XUTILITY_

#include <qwidget>
#include <Q3DockWindow>

#include "SdfViewer.h"

#include "WorldManager.h"
#include "FileBrowser.h"
#include "AppManager.h"
#include "MainWindow.h"

#include "FeaturePainters.h"

static GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.1f, 1.0f};
static GLfloat light0_diffuse[] =  {1.0f, 1.0f, 1.0f, 1.0f};
static GLfloat light0_position[] = {.5f, .5f, 1.0f, 0.0f};
static GLfloat light0_specular[] =  {1.0f, 1.0f, 1.0f, 1.0f};
static GLfloat specref[]        = { 0.4f, 0.4f, 0.4f, 1.0f };

WorldManager::WorldManager(SdfViewer* viewer, MainWindow* parent, const char* name) :
QObject(parent, name), m_viewer(viewer), m_selectionMode(SELECTION_MODE_OBJECT), m_displayMode(DRAW_MODE_SMOOTH), m_width(0), m_height(0)
{
	m_defaultRenderingParams = new RenderingParams;
	loadDefaultRenderingParams();

	//create default debug painter
	addDebugPainter(new NoneDebugPainter);
	addDebugPainter(new RayCellsPainter);
	addDebugPainter(new SdfFacetDebugPainter);
	addDebugPainter(new SdfVertexDebugPainter);
	addDebugPainter(new DihedralDebugPainter("Dihedral 1", 0));
	//addDebugPainter(new CagePainter);

	//create default facet painter
	addFacetPainter(new SolidFacetPainter);
	addFacetPainter(new SdfFacetsPainter("SDF Facets", true));
	addFacetPainter(new PartitionPainter);

	//addFacetPainter(new SdfVerticesPainter);
	addFacetPainter(new SdfDifferencesFacetsPainter);
	addFacetPainter(new DeformationDegreeFacetsPainter);
	addFacetPainter(new DeformationDegreeDifferencesFacetsPainter);

	m_rdp = new RenderingParamsDialog((QWidget*)parent, "rdp");
	//give the dialog the default rendering params
	m_rdp->OnRenderingParams(m_defaultRenderingParams);
	//when dialog updates setting, be notified
	connect(m_rdp, SIGNAL(renderingParamsChanged(RenderingParams*)), SLOT(OnRenderingPrametersChanged(RenderingParams*)));

	m_cgl = new ConfigureGlLights(parent);
	connect(m_cgl, SIGNAL(lightChanged()), this, SLOT(invalidate()));
}

WorldManager::~WorldManager()
{
	saveDefaultRenderingParams();

	delete m_defaultRenderingParams;

	//removeAllObjects();

	// responsible for killing painters too
	for (int i=0; i<m_facetPainters.size(); i++)
	{
		delete m_facetPainters[i];
	}
	m_facetPainters.clear();

	for (int i=0; i<m_debugPainters.size(); i++)
	{
		delete m_debugPainters[i];
	}
	m_debugPainters.clear();

}

void WorldManager::OnRenderingPrametersChanged(RenderingParams* renderingParams)
{
	if (renderingParams == m_defaultRenderingParams)
	{
		saveDefaultRenderingParams();
		//let all listeners know default parameters have changed
		selectFacetPainter(renderingParams->m_renderModeFacets);
		emit RenderingParametersChanged();
	}
}

void WorldManager::addCage(AppObject *appObject, TStubData* onlyPart, DrawMethod dm)
{
	WorldObject* wo = getSelectedObject();
	if (!wo) return;

	wo->id = appObject->id;
	wo->description = appObject->description;
	wo->cage = appObject->cage;
	wo->cage->centerOfMass();
	wo->cageSelection = appObject->cageSelection;
	wo->renderingParams = m_defaultRenderingParams;
	wo->myAo = appObject;

	QString cageRendererName = "cageRenderer_" + appObject->id;
	wo->cageRenderer = new CageRenderer(this, cageRendererName.toAscii());
	wo->cageRenderer->setRenderingParams(m_defaultRenderingParams);
	wo->cageRenderer->OnFacetPainterChanged(getSelectedFacetPainter());
	wo->cageRenderer->setCage(appObject->cage, NULL, appObject->cageSelection);

	//renderer is notified when rendering params change
	wo->cageRenderer->connect(this, SIGNAL(RenderingParametersChanged()), SLOT(OnObjectChanged()));

	//renderer is notified when painter is changed
	bool success = wo->cageRenderer->connect(this, SIGNAL(FacetPainterSelected(FacetPainter*)), SLOT(OnFacetPainterChanged(FacetPainter*)));

	m_objects.insert(wo->id, wo);
	
	int iindex = 0;

	Point_3 centerMass = wo->cage->centerOfMass();
	float diagonalLength = wo->cage->diagonalLength();
	
	//find index of added object in objects
	int ind=0;
	for (worldObjectsMap_t::iterator it = m_objects.begin();it != m_objects.end();it++)
	{
		if (it.value() == wo)
			break;
		ind++;
	}
	
	int selMode = m_selectionMode;
	m_selectionMode = SELECTION_MODE_OBJECT; // need to be in object mode to select the object
	
	OnObjectSelected(ind, false);
	m_selectionMode = selMode;

	//emit SceneParametersChanged(centerMass.x(), centerMass.y(), centerMass.z(), diagonalLength/2);
	// important: should set scale before object is loaded	
	calculateNewScene(iindex);

	appObject->setWorldObject(wo);

	emit finishedPaint();

	return;
}

void WorldManager::addObject(AppObject* appObject, TStubData* onlyPart, DrawMethod dm)
{
	WorldObject* wo = new WorldObject(this);
	wo->id = appObject->id;
	wo->description = appObject->description;
	wo->mesh = appObject->mesh;
	wo->mesh->centerOfMass();
	wo->meshSelection = appObject->meshSelection;
	wo->renderingParams = m_defaultRenderingParams;
	wo->myAo = appObject;

	QString meshRendererName = "meshRenderer_" + appObject->id;

	wo->meshRenderer = new MeshRenderer(this, meshRendererName.toAscii());
	wo->meshRenderer->setRenderingParams(m_defaultRenderingParams);
	wo->meshRenderer->OnFacetPainterChanged(getSelectedFacetPainter());
	wo->meshRenderer->setMesh(appObject->mesh, NULL, appObject->meshSelection);

	//renderer is notified when rendering params change
	wo->meshRenderer->connect(this, SIGNAL(RenderingParametersChanged()), SLOT(OnObjectChanged()));

	//renderer is notified when painter is changed
	bool success = wo->meshRenderer->connect(this, SIGNAL(FacetPainterSelected(FacetPainter*)), SLOT(OnFacetPainterChanged(FacetPainter*)));

	// scale by this object
	/*if (m_objects.size() == 0)
	{
		m_globalScale = mesh->diagonalLength();
		wo->meshRenderer->scale(1.0);
	}
	else
	{
		wo->meshRenderer->scale(m_globalScale / mesh->diagonalLength());
	}*/

	m_objects.insert(wo->id, wo);
	
	int iindex = 0;

	Point_3 centerMass = wo->mesh->centerOfMass();
	float diagonalLength = wo->mesh->diagonalLength();
	
	//find index of added object in objects
	int ind=0;
	for (worldObjectsMap_t::iterator it = m_objects.begin();it != m_objects.end();it++)
	{
		if (it.value() == wo)
			break;
		ind++;
	}
	
	int selMode = m_selectionMode;
	m_selectionMode = SELECTION_MODE_OBJECT; // need to be in object mode to select the object
	
	OnObjectSelected(ind, false);
	m_selectionMode = selMode;

	//emit SceneParametersChanged(centerMass.x(), centerMass.y(), centerMass.z(), diagonalLength/2);
	// important: should set scale before object is loaded	
	calculateNewScene(iindex);

	appObject->setWorldObject(wo);

	emit finishedPaint();

	return;
}

void WorldManager::invalidate()
{
	m_viewer->updateGL();
}

void WorldManager::calculateNewScene(int tabIndex)
{
	int rows;
	int cols;
	switch (m_objects.size()) {
	case 1: rows = 1; cols = 1; break;
	case 2: rows = 1; cols = 2; break;
	case 3: rows = 2; cols = 2; break;
	case 4: rows = 2; cols = 2; break;
	case 5: rows = 2; cols = 3; break;
	case 6: rows = 2; cols = 3; break;
	case 7: rows = 2; cols = 4; break;
	case 8: rows = 2; cols = 4; break;
	case 9: rows = 3; cols = 3; break;
	case 10: rows = 3; cols = 4; break;
	case 11: rows = 3; cols = 4; break;
	case 12: rows = 3; cols = 4; break;
	case 13: rows = 4; cols = 4; break;
	case 14: rows = 4; cols = 4; break;
	case 15: rows = 4; cols = 4; break;
	case 16: rows = 4; cols = 4; break;
	default: rows = 1; cols = 1; break;
	}

	float centerx = 0.0;
	float centery = 0.0;
	float centerz = 0.0;

	int currentRow = 0;
	int currentCol = 0;
	float cmod = (cols%2==0? 0.5 : 0.0);
	float rmod = (rows%2==0? 0.5 : 0.0);

	for (worldObjectsMap_t::iterator it = m_objects.begin(); it!= m_objects.end(); ++it)
	{
		WorldObject* wo = *it;

		Point_3 pmin = Point_3(wo->mesh->xmin(), wo->mesh->ymin(), wo->mesh->zmin());
		Point_3	pmax = Point_3(wo->mesh->xmax(), wo->mesh->ymax(), wo->mesh->zmax());

		float longest = qMax(pmax.z() - pmin.z(), qMax(pmax.x() - pmin.x(), pmax.y() - pmin.y()));
		float scale = 1.0 / longest;
		wo->meshRenderer->scale(scale);

		float x = (currentCol - cols / 2 + cmod);
		float y = (currentRow - rows / 2 + rmod);
		float z = 0.0;
		wo->meshRenderer->frame()->setPosition(x,y,z);

		//qglviewer::Quaternion q(qglviewer::Vec(1,0,0), -0.77);
		//qglviewer::Quaternion q2(qglviewer::Vec(0,1,0), 0.34);
		//wo->meshRenderer->frame()->setOrientation(q2*q);

		QString temp;
		temp.sprintf("calculateNewScene: <%s> scale=%f, translation=[%f,%f,%f]", wo->id.toAscii(), wo->meshRenderer->scale(), x, y, z);
		SDFLOG6(temp);
		const GLdouble* m = wo->meshRenderer->frame()->matrix();
		temp.sprintf("calculateNewScene: Matrix=[%.2f %.2f %.2f %.2f ; %.2f %.2f %.2f %.2f ; %.2f %.2f %.2f %.2f ; %.2f %.2f %.2f %.2f]",
			m[0],m[1],m[2],m[3], m[4],m[5],m[6],m[7],m[8],m[9],m[10],m[11],m[12],m[13],m[14],m[15]);
		SDFLOG6(temp);

		currentCol++;
		if (currentCol == cols) {
			currentRow++;
			currentCol = 0;
		}
	}

	emit SceneParametersChanged(tabIndex, centerx, centery, centerz, (float) qMax(rows,cols) / 2);
}


void WorldManager::removeObject(const QString& id)
{
	worldObjectsMap_t::iterator it = m_objects.find(id.ascii());

	if (it == m_objects.end())
		return;

	WorldObject* wo = it.value();

	invalidateDebugPainters(wo->mesh);

	if (m_viewer->manipulatedFrame() == wo->meshRenderer->frame())
		m_viewer->setManipulatedFrame(NULL);

	delete wo;

	m_objects.erase(it);

	calculateNewScene(0);
}

void WorldManager::removeAllObjects()
{
	//delete all objects from world
	for (worldObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end(); it++)
	{
		WorldObject* wo = it.value();
		delete wo;
	}

	invalidateDebugPainters();

	calculateNewScene(0);
}

const GLdouble* WorldManager::getObjectMatrix(const QString& id)
{
	WorldObject* wo = findObject(id);

	if (wo) {
		return wo->meshRenderer->frame()->matrix();
	} else {
		return NULL;
	}
}

WorldObject* WorldManager::findObject(const QString& id)
{
	worldObjectsMap_t::iterator it = m_objects.find(id.ascii());
	if (it != m_objects.end()) {
		return it.value();
	} else {
		return NULL;
	}
}

void WorldObject::drawObject(int name, int selectionMode)
{
	glPushMatrix();
	glMultMatrixd(meshRenderer->frame()->matrix());
	float localScale = meshRenderer->scale();
	glScalef(localScale, localScale, localScale);

	// translate to [0,0,0]

	Point_3 center = mesh->computedCenterOfMass();

	glTranslatef(-center.x(), -center.y(), -center.z());

	//QGLViewer::drawAxis();

	if (cageRenderer != NULL)
		cageRenderer->render(selectionMode);

	if (selectionMode == SELECTION_MODE_OBJECT)
	{
		//render all objects
		glPushName(name);
		meshRenderer->render(selectionMode);
		glPopName();
	}
	else if (selectionMode == SELECTION_MODE_FACET) 
	{
		if (meshSelection->isSelected())
		{
			//if selectin facet work only with selected object
			meshRenderer->render(selectionMode);
		}
		else
		{
			glPushName(INT_MAX-name);
			meshRenderer->render(SELECTION_MODE_OBJECT);
			glPopName();
		}
	}
	else if (selectionMode == SELECTION_MODE_NONE)
	{
		DebugPainter *dbgp = mgr->getSelectedDebugPainter();
		//if no selection enabled render everyone
		if ((dbgp != NULL) && selectionMode == SELECTION_MODE_NONE)
		{
			mgr->getSelectedDebugPainter()->paint(mesh, meshSelection, NULL, this);
		}

		bool hasAlpha = (mgr->defaultRenderParams()->m_alpha<255);
		if (hasAlpha)
		{
			//glDepthMask(GL_FALSE);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}

		meshRenderer->render(selectionMode);

		//if no selection enabled render everyone
		if ((dbgp != NULL) && selectionMode == SELECTION_MODE_NONE)
		{
			mgr->getSelectedDebugPainter()->postPaint(mesh, meshSelection, NULL);
		}

		if (hasAlpha)
		{
			//glDepthMask(GL_TRUE);
			glDisable(GL_BLEND);
		}
	}

	glPopMatrix();
}

void WorldObject::drawCage(int name, int selectionMode)
{
	glPushMatrix();
	glMultMatrixd(meshRenderer->frame()->matrix());
	float localScale = meshRenderer->scale();
	glScalef(localScale, localScale, localScale);

	// translate to [0,0,0]

	Point_3 center = mesh->computedCenterOfMass();

	glTranslatef(-center.x(), -center.y(), -center.z());

	//QGLViewer::drawAxis();

	if (selectionMode == SELECTION_MODE_OBJECT)
	{
		//render all objects
		glPushName(name);
		cageRenderer->render(selectionMode);
		glPopName();
	}
	else if (selectionMode == SELECTION_MODE_FACET) 
	{
		if (meshSelection->isSelected())
		{
			//if selectin facet work only with selected object
			cageRenderer->render(selectionMode);
		}
		else
		{
			glPushName(INT_MAX-name);
			cageRenderer->render(SELECTION_MODE_OBJECT);
			glPopName();
		}
	}
	else if (selectionMode == SELECTION_MODE_NONE)
	{
		DebugPainter *dbgp = mgr->getSelectedDebugPainter();

		//if no selection enabled render everyone
		if ((dbgp != NULL) && selectionMode == SELECTION_MODE_NONE)
		{
			mgr->getSelectedDebugPainter()->paint(mesh, meshSelection, NULL, this);
		}

		bool hasAlpha = (mgr->defaultRenderParams()->m_alpha < 255);
		if (hasAlpha)
		{
			//glDepthMask(GL_FALSE);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}

		cageRenderer->render(selectionMode);

		//if no selection enabled render everyone
		if ((dbgp != NULL) && selectionMode == SELECTION_MODE_NONE)
		{
			mgr->getSelectedDebugPainter()->postPaint(mesh, meshSelection, NULL);
		}

		if (hasAlpha)
		{
			//glDepthMask(GL_TRUE);
			glDisable(GL_BLEND);
		}
	}

	glPopMatrix();
}


void WorldManager::drawWorld(int selectionMode)
{
	int i=0;
	for (worldObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end(); it++)
	{
		WorldObject* wo = it.value();
		
		wo->drawObject(i, selectionMode);

		if (wo->cageRenderer !=NULL)
			wo->drawCage(i, selectionMode);
			
		i++;
	}

	//if not drawing for selection and overlay enabled, draw overlay
	if (selectionMode == SELECTION_MODE_NONE && m_defaultRenderingParams->m_overlay)
	{
		drawOverlay();
	}

	glFlush();
}

void WorldManager::drawOverlay()
{
	//setup GL for 2d overlay
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0,1,0,1);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, &viewport[0]);
	glViewport(0,0,m_width >> 3, m_height >> 3);
	for (int i=0; i<m_overlayPainters.size(); i++)
	{
		m_overlayPainters[i]->render();
	}
	glViewport(0.0, 0.0, viewport[2], viewport[3]);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}


//////////////////////////////////////////////////////////////////////////
// Slots
//////////////////////////////////////////////////////////////////////////

void WorldManager::OnDraw()
{
	drawWorld(SELECTION_MODE_NONE);

	//emit RedrawNeeded(); does nothing
}

void WorldManager::OnSelectDraw()
{
	drawWorld(m_selectionMode);
}

/** add a new cage */
void WorldManager::OnAddCage(AppObject *appObject, TStubData* onlyPart, DrawMethod dm)
{
	addCage(appObject, onlyPart, dm);

	emit RedrawNeeded();
}

/** add a new object */
void WorldManager::OnAddObject(AppObject* appObject, TStubData* onlyPart, DrawMethod dm)
{
	addObject(appObject, onlyPart, dm);

	emit RedrawNeeded();
}


/** remove an object */
void WorldManager::OnRemoveObject(const QString& id)
{
	removeObject(id);

	emit RedrawNeeded();
}

void WorldManager::OnChangeObject(const QString& id, bool geometryChanged)
{
	worldObjectsMap_t::iterator it = m_objects.find(id.ascii());

	if (it != m_objects.end())
	{
		WorldObject* wo = it.value();
		invalidateDebugPainters(wo->mesh);
		wo->meshRenderer->OnObjectChanged(geometryChanged);
	}

	emit RedrawNeeded();
}

void WorldManager::OnViewerInitialized()
{
	//init glew
	GLenum err = glewInit();
	if (GLEW_OK == err)
	{
		m_defaultRenderingParams->m_vboSupported = GLEW_ARB_vertex_buffer_object;
	}
	else
	{
		m_defaultRenderingParams->m_vboSupported = false;
	}

	glEnable(GL_COLOR_MATERIAL);
	//glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	glMaterialfv(GL_FRONT, GL_SPECULAR, specref);
	glMateriali(GL_FRONT,GL_SHININESS, 60);

	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

	//glFrontFace(GL_CW);

	//Load default configuration from lights dialog
	m_cgl->LoadConfiguration(0);
	m_cgl->applyLight();

	//TODO: temporary
	glEnable(GL_NORMALIZE);
}

void WorldManager::OnEditRenderingParameters()
{
	m_rdp->onFillFacetPainters(getFacetPainterNames());
	m_rdp->OnFillDebugPainters(getDebugPainterNames());
	m_rdp->OnRenderingParams(m_defaultRenderingParams);
	m_rdp->show();
}

void WorldManager::saveRenderingParams(const QString& prefix, RenderingParams* renderingParams)
{
	ProgSettings ps;

	ps.beginGroup(prefix);

	ps.writeEntry("smoothShading", renderingParams->m_smoothShading);
	ps.writeEntry("culling",renderingParams->m_culling);
	ps.writeEntry("polygonMode",renderingParams->m_polygonMode);
	ps.writeEntry("antialiasing",renderingParams->m_antialiasing);
	ps.writeEntry("superimposeVertices",renderingParams->m_superimposeVertices);
	ps.writeEntry("superimposeEdges",renderingParams->m_superimposeEdges);
	ps.writeEntry("useNormals",renderingParams->m_useNormals);
	ps.writeEntry("smoothNormals", renderingParams->m_smoothNormals);
	ps.writeEntry("lighting",renderingParams->m_lighting);

	ps.writeEntry("facetColor", (int)renderingParams->m_facetColor.rgb());
	ps.writeEntry("edgeColor", (int)renderingParams->m_edgeColor.rgb());
	ps.writeEntry("vertexColor", (int)renderingParams->m_vertexColor.rgb());

	ps.writeEntry("renderModeFacets", renderingParams->m_renderModeFacets);
	ps.writeEntry("renderModeEdges", renderingParams->m_renderModeEdges);

	ps.writeEntry("debugPainter", renderingParams->m_debugPainter);

	ps.writeEntry("overlay", renderingParams->m_overlay);

	ps.writeEntry("vertexThickness", renderingParams->m_vertexThickness);
	ps.writeEntry("edgeThickness", renderingParams->m_edgeThickness);

	ps.writeEntry("colorSchemeName", renderingParams->m_colorSchemeName);

	ps.writeEntry("alpha", renderingParams->m_alpha);

	ps.endGroup();
}

void WorldManager::loadRenderingParams(const QString& prefix, RenderingParams* renderingParams)
{
	ProgSettings ps;

	ps.beginGroup(prefix);

	renderingParams->m_smoothShading = ps.readBoolEntry("smoothShading", true);
	renderingParams->m_culling = ps.readBoolEntry("culling", true);

	renderingParams->m_antialiasing = ps.readBoolEntry("antialiasing", true);
	renderingParams->m_superimposeVertices = ps.readBoolEntry("superimposeVertices", false);
	renderingParams->m_superimposeEdges = ps.readBoolEntry("superimposeEdges", false);
	renderingParams->m_useNormals = ps.readBoolEntry("useNormals", true);
	renderingParams->m_smoothNormals = ps.readBoolEntry("smoothNormals", true);
	renderingParams->m_lighting = ps.readBoolEntry("lighting", true);

	renderingParams->m_polygonMode = ps.readNumEntry("polygonMode", 0x1b02);

	renderingParams->m_facetColor = QColor((QRgb)ps.readNumEntry("facetColor", 0xffffffff));
	renderingParams->m_edgeColor = QColor((QRgb)ps.readNumEntry("edgeColor", 0xfff0d400));
	renderingParams->m_vertexColor = QColor((QRgb)ps.readNumEntry("vertexColor", 0xfff06d6d));

	renderingParams->m_overlay = ps.readBoolEntry("overlay", false);

	renderingParams->m_renderModeFacets = ps.readNumEntry("renderModeFacets", 1);	// SDF
	renderingParams->m_renderModeEdges = ps.readNumEntry("renderModeEdges", 0);

	renderingParams->m_debugPainter = ps.readNumEntry("debugPainter", 0);

	renderingParams->m_vertexThickness = ps.readNumEntry("vertexThickness", 1);
	renderingParams->m_edgeThickness = ps.readNumEntry("edgeThickness", 1);

	renderingParams->m_colorSchemeName = ps.readEntry("colorSchemeName", "default");

	renderingParams->m_alpha = ps.readNumEntry("alpha", 255);

	ps.endGroup();
}

void WorldManager::saveDefaultRenderingParams()
{
	saveRenderingParams("default", m_defaultRenderingParams);
}

void WorldManager::loadDefaultRenderingParams()
{
	loadRenderingParams("default", m_defaultRenderingParams);
}

void WorldManager::OnEditGLLights()
{
	m_cgl->show();
}

int WorldManager::addFacetPainter(FacetPainter* facetPainter)
{
	m_facetPainters.push_back(facetPainter);
	return m_facetPainters.size();
}

QStringList WorldManager::getFacetPainterNames()
{
	QStringList fpnames;
	for (int i=0;i<m_facetPainters.size();i++)
	{
		if (m_facetPainters[i] != NULL)
		{
			fpnames.append(m_facetPainters[i]->name());
		}
	}

	return fpnames;
}

int WorldManager::selectFacetPainter(const int facetPainterIndex)
{
	//int oldFacetPainterIndex = m_defaultRenderingParams->m_renderModeFacets;

	m_defaultRenderingParams->m_renderModeFacets = facetPainterIndex;

	//recalculate colors
	//if (facetPainterIndex != oldFacetPainterIndex)
	emit FacetPainterSelected(m_facetPainters[m_defaultRenderingParams->m_renderModeFacets]);

	//redreaw needed
	m_viewer->updateGL();

	return m_defaultRenderingParams->m_renderModeFacets;
}

int WorldManager::selectFacetPainter(const QString& facetPainterName)
{
	int oldFacetPainterIndex = m_defaultRenderingParams->m_renderModeFacets;

	for (int i=0;i<m_facetPainters.size();i++)
	{
		if (m_facetPainters[i] != NULL && m_facetPainters[i]->name() == facetPainterName) {
			m_defaultRenderingParams->m_renderModeFacets = i;

			//recalculate colors
			if (i != oldFacetPainterIndex)
				emit FacetPainterSelected(m_facetPainters[m_defaultRenderingParams->m_renderModeFacets]);

			return m_defaultRenderingParams->m_renderModeFacets;
		}
	}

	return -1;
}


FacetPainter* WorldManager::getSelectedFacetPainter()
{
	if (!(m_defaultRenderingParams &&
		m_defaultRenderingParams->m_renderModeFacets >=0 &&
		m_defaultRenderingParams->m_renderModeFacets < m_facetPainters.size()))
		return NULL;

	return m_facetPainters[m_defaultRenderingParams->m_renderModeFacets];
}

//////////////////////////////////////////////////////////////////////////

int WorldManager::addDebugPainter(DebugPainter* debugPainter)
{
	m_debugPainters.push_back(debugPainter);
	return m_debugPainters.size();
}

QStringList WorldManager::getDebugPainterNames()
{
	QStringList fpnames;
	for (int i=0;i<m_debugPainters.size();i++)
	{
		if (m_debugPainters[i] != NULL) {
			fpnames.append(m_debugPainters[i]->name());
		}
	}

	return fpnames;
}

int WorldManager::selectDebugPainter(const int debugPainterIndex)
{
	m_defaultRenderingParams->m_debugPainter = debugPainterIndex;

	//recalculate colors
	emit DebugPainterSelected(m_debugPainters[m_defaultRenderingParams->m_debugPainter]);

	return m_defaultRenderingParams->m_debugPainter;
}

int WorldManager::selectDebugPainter(const QString& debugPainterName)
{
	int oldDebugPainterIndex = m_defaultRenderingParams->m_debugPainter;

	for (int i=0;i<m_debugPainters.size();i++)
	{
		if (m_debugPainters[i] != NULL && m_debugPainters[i]->name() == debugPainterName) {
			m_defaultRenderingParams->m_debugPainter = i;

			//recalculate colors
			if (i != oldDebugPainterIndex)
				emit DebugPainterSelected(m_debugPainters[m_defaultRenderingParams->m_debugPainter]);

			return m_defaultRenderingParams->m_debugPainter;
		}
	}

	return -1;
}


int WorldManager::addOverlayPainter(OverlayPainter* overlayPainter)
{
	m_overlayPainters.push_back(overlayPainter);

	return m_overlayPainters.size();
}


DebugPainter* WorldManager::getSelectedDebugPainter()
{
	if (!(m_defaultRenderingParams &&
		m_defaultRenderingParams->m_debugPainter >=0 &&
		m_defaultRenderingParams->m_debugPainter < m_debugPainters.size()))
		return NULL;

	return m_debugPainters[m_defaultRenderingParams->m_debugPainter];
}

WorldObject* WorldManager::getSelectedObject()
{
	worldObjectsMap_t::iterator it = m_objects.begin();
	worldObjectsMap_t::iterator it_end = m_objects.end();
	for (; it != it_end; it++)
	{
		if (it.value()->meshSelection->isSelected())
		{
			return it.value();
		}
	}

	return NULL;
}

// select by either index or id
void WorldManager::objectSelected(int i, QString id)
{
	WorldObject* selwo = NULL;
	int ind = 0;
	for (worldObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end(); ++it)
	{
		WorldObject* wo = it.value();
		if ( ((i != -1) && (ind == i)) ||
			((!id.isEmpty()) && (wo->id == id)) )
		{
			selwo = wo;
		}
		else
		{
			wo->meshSelection->setSelected(false);
		}
		ind++;
	}

	if (selwo == NULL)
	{
		return;
	}
	m_viewer->setManipulatedFrame(selwo->meshRenderer->frame());
	selwo->meshSelection->setSelected(true);

	invalidateDebugPainters(selwo->mesh);

	QString v = QString::number(selwo->mesh->size_of_vertices());
	QString f = QString::number(selwo->mesh->size_of_facets());

	emit message("Selected object <" + selwo->id + ">" + " V = " + v + " , F = " + f);

	emit ObjectSelected(selwo->id);
}


// entry point
void WorldManager::OnObjectSelected(unsigned int i, Qt::KeyboardModifiers modif, int extPart /* =-1 */)
{
	bool multi = false; //(modif == Qt::ShiftModifier);
	bool rotsel = (modif == Qt::ShiftModifier);

	if (m_selectionMode == SELECTION_MODE_OBJECT)
	{
		objectSelected(i);
	}
	else if (m_selectionMode == SELECTION_MODE_FACET)
	{
		WorldObject* wo = getSelectedObject();
		if (!wo) return;

		if (UINT_MAX == i)
		{
			wo->meshSelection->clear();
		}
		else if ((i == UINT_MAX - 1) || (i < wo->mesh->size_of_facets()))
		{
			if (i == UINT_MAX - 1)
				i = wo->meshSelection->selectedFacet();

			invalidateDebugPainters(wo->mesh); // invalidate debug painters

			if (multi)
				wo->meshSelection->selectFacetAdd(i);
			else
				wo->meshSelection->selectFacet(i);

			//OnChangeObject(wo->id, false);
			emit FacetSelected(wo->id, i);

			Mesh::Facet_const_handle f = wo->mesh->findFacet(i);
			if (f != NULL)
			{
				/*if (multi)
					wo->meshSelection->selectPartAdd(f->part());
				else
					wo->meshSelection->selectPart(f->part());*/

				QString temp;
				temp.sprintf("Selected facet <%d> in %s. [SDF=%.3f; NSDF=%.3f; Cluster=%d]", i, wo->description.ascii(), f->volumeSDF(), f->volumeNSDF(), f->cluster());

				emit message(temp);
			}
			else
			{
				emit message("Selected facet <" + QString::number(i) + "> in object " + wo->description);
			}
		}
		else
		{
			objectSelected(INT_MAX - i);
		}
	}
	else if (m_selectionMode == SELECTION_MODE_VERTEX)
	{
		WorldObject* wo = getSelectedObject();
		if (!wo) return;

		if (UINT_MAX == i)
		{
			wo->meshSelection->clear();
		}
		else if ((i == UINT_MAX - 1) || (i < wo->mesh->size_of_vertices()))
		{
			if (i == UINT_MAX - 1)
				i = wo->meshSelection->selectedVertex();

			invalidateDebugPainters(wo->mesh); // invalidate debug painters

			/*if (multi)
				wo->meshSelection->selectVertexAdd(i);
			else
				wo->meshSelection->selectVertex(i);*/

			//OnChangeObject(wo->id, false);
			emit VertexSelected(wo->id, i);

			Mesh::Vertex_const_handle v = wo->mesh->findVertex(i);
			if (v != NULL)
			{
				/*if (multi)
					wo->meshSelection->selectPartAdd(f->part());
				else
					wo->meshSelection->selectPart(f->part());*/

				QString temp;
				temp.sprintf("Selected vertex <%d> in %s. [SDF=%.3f; NSDF=%.3f; Cluster=%d]", i, wo->description.ascii(), v->volumeSDF(), v->volumeNSDF(), v->cluster());
				
				printf("Selected vertex <%d> in %s. [SDF=%.3f; NSDF=%.3f; Cluster=%d]", i, wo->description.ascii(), v->volumeSDF(), v->volumeNSDF(), v->cluster());

				emit message(temp);
			}
			else
			{
				printf("Selected vertex <" + QString::number(i) + "> in object " + wo->description);
				emit message("Selected vertex <" + QString::number(i) + "> in object " + wo->description);
			}
		}
		else
		{
			objectSelected(INT_MAX - i);
		}
	}
}

void WorldManager::OnDrawAxis(bool isDrawAxis)
{
	m_drawAxis = isDrawAxis;
}

void WorldManager::OnSelectionMode(int selectionMode)
{
	m_selectionMode = selectionMode;
}

void WorldManager::OnDisplayMode(int displayMode)
{
	m_displayMode = displayMode;
	
	switch (displayMode)
	{
		case 1: 
			cout << "points\n"; 
			m_defaultRenderingParams->m_polygonMode = 0x1B00;
			m_defaultRenderingParams->m_superimposeVertices = true;
			m_defaultRenderingParams->m_superimposeEdges = false;
			break;
		case 2: 
			cout << "wireframe\n"; 
			m_defaultRenderingParams->m_polygonMode = 0x1B01;
			m_defaultRenderingParams->m_superimposeVertices = false;
			m_defaultRenderingParams->m_superimposeEdges = false;
			break;
		case 3: 
			cout << "flatlines\n"; 
			m_defaultRenderingParams->m_polygonMode = 0x1B02;
			m_defaultRenderingParams->m_superimposeEdges = true;
			break;
		case 4: 
			cout << "flat\n"; 
			m_defaultRenderingParams->m_smoothShading = false;
			m_defaultRenderingParams->m_polygonMode = 0x1B02;
			m_defaultRenderingParams->m_superimposeVertices = false;
			m_defaultRenderingParams->m_superimposeEdges = false;
			break;
		case 5: 
			m_defaultRenderingParams->m_smoothShading = true;
			m_defaultRenderingParams->m_polygonMode = 0x1B02;
			m_defaultRenderingParams->m_superimposeVertices = false;
			m_defaultRenderingParams->m_superimposeEdges = false;
			cout << "smooth\n"; 
			break;
		default: 
			cout << "error\n"; break;
	}

	drawWorld(SELECTION_MODE_NONE);

	//emit RenderingParametersChanged();
	//emit RedrawNeeded();
}

void WorldManager::OnWorldSizeChanged(int width, int height)
{
	m_height = height;
	m_width = width;
}

void WorldManager::OnSelectFacetPainter(int facetPainter)
{
	if (facetPainter<0 || facetPainter>=m_facetPainters.size() || facetPainter == m_defaultRenderingParams->m_renderModeFacets)
		return;

	selectFacetPainter(facetPainter);
	emit RenderingParametersChanged();

	QString temp;
	temp.sprintf("Selected painter %s", m_facetPainters[facetPainter]->name().ascii());
	emit message(temp);

	glDisable(GL_LIGHTING);
	m_viewer->displayMessage(temp);
}

void WorldManager::invalidateDebugPainters(Mesh* mesh)
{
	for (int i=0;i<m_debugPainters.size();i++)
		m_debugPainters[i]->invalidate(mesh);
}

void WorldManager::OnFunctionKeyPressed(int key, Qt::ButtonState state)
{
	if (state == Qt::NoButton)
	{
		OnSelectFacetPainter(key);
	}
	else if (state == Qt::ControlButton)
	{
		//notify debug painter
		DebugPainter* dp = getSelectedDebugPainter();
		if (dp != NULL)
		{
			QString answer = dp->command(key);
			if (answer != QString::null)
			{
				glDisable(GL_LIGHTING);
				m_viewer->displayMessage(answer);
			}
		}
	}
	else if (state == Qt::AltButton)
	{
		FacetPainter* fp = getSelectedFacetPainter();
		if (fp != NULL)
		{
			QString answer = fp->command(key);
			if (answer != QString::null)
			{
				QString s = QString("Facet Painter %1 command %2\n").arg(fp->name()).arg(answer);

				printf(s.toAscii());

				for (worldObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end();it++)
					it.value()->meshRenderer->OnObjectChanged(false);

				glDisable(GL_LIGHTING);
				m_viewer->displayMessage(answer);
			}
		}
	}
}


void WorldManager::redoAllColors()
{
	for (worldObjectsMap_t::iterator it = m_objects.begin(); it != m_objects.end();it++)
		it.value()->meshRenderer->OnObjectChanged(false);
	invalidate();
}