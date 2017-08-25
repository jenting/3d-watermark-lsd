#ifndef __WORLD_MANAGER_H_
#define __WORLD_MANAGER_H_

#include "stdafx.h"

#include <hash_map>
#include <string>

#include "RenderingParams.h"
#include "MeshRenderer.h"

#include "BasicRenderingParams.h"
#include "glLights.h"

class SdfViewer;

class AppObject;
class WorldManager;


/**
 *	Represents an object in the world that needs to be rendered and can be controlled.
 *  Each object should have geometrical data, identification, location in world, rendering
 *  parameters...
 */
struct WorldObject
{
public:
	WorldObject(WorldManager* m) : mgr(m), mesh(NULL), cage(NULL), meshSelection(NULL), cageSelection(NULL), renderingParams(NULL), meshRenderer(NULL), cageRenderer(NULL), myAo(NULL)
	{
	}

	~WorldObject()
	{
		delete meshRenderer;
	}

	void drawObject(int name, int selectionMode);
	void drawCage(int name, int selectionMode);

	WorldManager *mgr;
	AppObject *myAo;

public:
	/** unique string identifying the object */
	QString id;
	/** textual description of the object */
	QString description;
	/** pointer to mesh to be rendered */
	Mesh* mesh;
	/** pointer to cage to be rendered */
	Mesh* cage;
	/** pointer to selection object on mesh */
	MeshSelection* meshSelection;
	/** pointer to selection cage */
	MeshSelection* cageSelection;
	/** pointer to rendering parameters structure */
	RenderingParams* renderingParams;
	/** Object responsible for rendering mesh onscreen */
	//MeshRenderer* renderer;
	MeshRenderer* meshRenderer;
	/** Object responsible for rendering cage onscreen */
	CageRenderer* cageRenderer;

	int *tabIndex; // in which tab I am living.

};


typedef QHash<QString, WorldObject*> worldObjectsMap_t;


/**
 *	This class manages objects in the world, objects can be added to the manager
 * and it is responsible for displaying them, controlling them, reporting events
 * happening to them...
 *
 */
class WorldManager : public QObject
{
	Q_OBJECT
public:

	WorldManager(SdfViewer* viewer, MainWindow* parent = NULL, const char* name = NULL);

	~WorldManager();

	/** hash containing objects in the world */
	worldObjectsMap_t m_objects;

	/** default rendering params assigned to all new objects */
	RenderingParams* m_defaultRenderingParams;

	/** collection of facet painters */
	std::vector<FacetPainter*> m_facetPainters;

	/** collection of debug painters */
	std::vector<DebugPainter*> m_debugPainters;

	/** collection of overlay painters */
	std::vector<OverlayPainter*> m_overlayPainters;

	/** dialog for editing rendering parameters */
	RenderingParamsDialog* m_rdp;

	/** dialog for configuring lighting in scene */
	ConfigureGlLights* m_cgl;

	/** pointer to the viewer */
	SdfViewer* m_viewer;

	/** keep scale */
	float m_globalScale;

	/** to draw axis or not. */
	bool m_drawAxis;

	/** 
	 * When user selects something in the scene, what does he select. 
	 * 0 for none, 1 for objects, 2 for facets, 3 for vertices
	 */
	int m_selectionMode;

	/** 
	 * When user selects something in the draw mode, what does he select.
	 * 0 for points, 1 for wireframe, 2 for flat lines, 3 for flat, 4 for smooth
	 */
	int m_displayMode; 

	int m_width;
	int m_height;

public:

		/** returns pointer to currently selected Debug painter */
	DebugPainter* getSelectedDebugPainter();
	const RenderingParams *defaultRenderParams() const { return m_defaultRenderingParams; }

	/**
	*	adds pointer to a facet painter
	* @param facetPainter pointer to the painter
	* @return number of current painters
	*/
	int addFacetPainter(FacetPainter* facetPainter);
	int addDebugPainter(DebugPainter* debugPainter);
	int addOverlayPainter(OverlayPainter* overlayPainter);

	/** selects facet painter and returns selected index */
	int selectFacetPainter(const int facetPainterIndex);
	/** selects facet painter and returns selected index */
	int selectFacetPainter(const QString& facetPainterName);

protected:
	/** locate object based on id */
	WorldObject* findObject(const QString& id);

	void saveRenderingParams(const QString& prefix, RenderingParams* renderingParams);

	void loadRenderingParams(const QString& prefix, RenderingParams* renderingParams);

	void saveDefaultRenderingParams();

	void loadDefaultRenderingParams();


	/** returns pointer to currently selected facet painter */
	FacetPainter* getSelectedFacetPainter();

	/** selects Debug painter and returns selected index */
	int selectDebugPainter(const int debugPainterIndex);
	/** selects Debug painter and returns selected index */
	int selectDebugPainter(const QString& debugPainterName);

	WorldObject* getSelectedObject();

	// select either by index in the objects array or by id
	void objectSelected(int i, QString id = QString());

	void calculateNewScene(int tab);

	void invalidateDebugPainters(Mesh* mesh = NULL);
	int currentTab();


public:
	int getSelectMode() { return m_selectionMode; }

	int getDisplayMode() { return m_displayMode; }

	/** add new cage to the world, identified by the id and description */
	void addCage(AppObject *appObject, TStubData* onlyPart, DrawMethod dm);

	/** add new object to the world, identified by the id and description */
	void addObject(AppObject* appObject, TStubData* onlyPart, DrawMethod dm);

	/** remove object identified by id */
	void removeObject(const QString& id);

	/** remove all objects */
	void removeAllObjects();

	/** draw all objects in the world */
	void drawWorld(int selectionMode = SELECTION_MODE_NONE);

	/** draw overlays */
	void drawOverlay();

	/** returns list of available facet painters */
	QStringList getFacetPainterNames();

	/** returns list of available debug painters */
	QStringList getDebugPainterNames();

	/**
	 *	returns the matrix of the given object
	 */
	const GLdouble* getObjectMatrix(const QString& id);

	QString getColorSchemeName() {  return m_defaultRenderingParams->m_colorSchemeName; }

	void emitMessage(const QString& s) { emit message(s); }
public slots:

	void redoAllColors();
	void invalidate();

	/** */
	void OnDrawAxis(bool);

	/** trigger drawing */
	void OnDraw();

	/** trigger selection drawing */
	void OnSelectDraw();

	/** add a new cage */
	void OnAddCage(AppObject *appObject, TStubData* onlyPart, DrawMethod dm);

	/** add a new object */
    void OnAddObject(AppObject *appObject, TStubData* onlyPart, DrawMethod dm);

	/** remove an object */
	void OnRemoveObject(const QString& id);

	void OnChangeObject(const QString& id, bool geometryChanged);

	/** viewer was initialized, set it up */
	void OnViewerInitialized();

	/** open dialog to edit rendering parameters */
	void OnEditRenderingParameters();

	/** the rendering params dialog reported a change */
	void OnRenderingPrametersChanged(RenderingParams* renderingParams);

	/** open dialog to edit lights in scene */
	void OnEditGLLights();

	void OnObjectSelected(unsigned int i, Qt::KeyboardModifiers modif, int extPart = -1);

	void OnSelectionMode(int selectionMode);
	void OnDisplayMode(int drawMode);

	void OnWorldSizeChanged(int width, int height);

	void OnSelectFacetPainter(int facetPainter);

	void OnFunctionKeyPressed(int key, Qt::ButtonState state);

signals:
	void finishedPaint();

	/** user selected an object in the world */
	void ObjectSelected(const QString& id);

	/** user selected facet in an object */
	void FacetSelected(const QString& id, const int facetIndex);

	/** user selected vertex in an object */
	void VertexSelected(const QString& id, const int vertexIndex);

	/** global rendering parameters changed */
	void RenderingParametersChanged();

	/** signal the viewer to refresh */
	void RedrawNeeded();

	/** tell the viewer what to show */
	void SceneParametersChanged(int tabIndex, float centerx,float centery, float centerz, float radius);

	void FacetPainterSelected(FacetPainter* facetPainter);

	void DebugPainterSelected(DebugPainter* debugPainter);

	void message(const QString& message);

	void message(const QString& message, int ms);

	void changeFacetPainter(const int index);

	void doRemoveObject(const QString& id);
};

#endif