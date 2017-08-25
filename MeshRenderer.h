#ifndef __MESH_RENDERER_H_
#define __MESH_RENDERER_H_

#include "stdafx.h"

#include <qstringlist.h>
#include "RenderingParams.h"

#include <QGLViewer/qglviewer.h>
#include <QGLViewer/manipulatedFrame.h>

#define SELECTION_MODE_NONE 0x0
#define SELECTION_MODE_OBJECT 0x1
#define SELECTION_MODE_FACET 0x2
#define SELECTION_MODE_VERTEX 0x3

#define DRAW_MODE_POINTS 0x0
#define DRAW_MODE_WIREFRAME 0x1
#define DRAW_MODE_FLAT_LINES 0x2
#define DRAW_MODE_FLAT 0x3
#define DRAW_MODE_SMOOTH 0x4

/**
 *	Abstract class which represents a generic painter, implemented to select colors based on specific features
 */
class Painter {
private:
	
protected:
	
public:
	Painter() {}

	/** name of the painter (for selection purposes) */
	virtual QString name() = 0;
};


class Indexer
{
public:
	virtual ~Indexer() {}
	virtual bool hasNext() = 0;
	virtual int next() = 0;
};

class CountIndexer : public Indexer
{
public:
	CountIndexer(int count) : m_count(count), m_cur(0) {}
	virtual ~CountIndexer() {}
	virtual bool hasNext() { return m_cur < m_count; }
	virtual int next() { Q_ASSERT(m_cur < m_count); return m_cur++; }
private:
	int m_count;
	int m_cur;
};

class ListIndexer : public Indexer
{
public:
	ListIndexer(const TIndexList& list) : m_list(list), m_it(list.constBegin()) {}
	virtual ~ListIndexer() {}
	virtual bool hasNext() { return m_it != m_list.constEnd(); }
	virtual int next() { Q_ASSERT(m_it != m_list.end()); return *(m_it++); }
private:
	const TIndexList& m_list;
	TIndexList::const_iterator m_it;
};


class FacetPainter: public Painter {
public:
	/** if true one color per vertex, o/w one color per vertex per facet */
	virtual bool isSmoothShading() = 0;

	/** what size of buffer must be assigned to this pointer */
	virtual unsigned int requiredBufferSize(Mesh* mesh) = 0;

	/**
		*	Fill given buffer of size requiredBufferSize() with color values
		* @param renderingParams pointer to rendering params with hints as to what to color
		* @param mesh pointer to mesh who holds feature values
		* @param colorBuffer allocated unsigned byte buffer of size requiredBufferSize()
		* @return true if successful
		*/
	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh, TStubData* onlyPart, DrawMethod dm,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer)
	{
		return fillBuffer(renderingParams, mesh, msg, meshSelection, colorBuffer);
	}

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer) { return true; }

	/**
		*	Paints one vertex
		* @param renderingParams pointer to rendering params with hints as to what to color
		* @param pVertex pointer to vertex to paint
		* @param color pointer to location where to place color
		* @return true if successful
		*/
	virtual bool paintVertex(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Vertex_const_handle& pVertex,
		GLubyte* color) 
	{
		return false;
	}

	/**
	*	give debug renderer a command and get back a string signalling the new mode
	*/
	virtual QString command(const int key) {
		return QString::null;
	}

	virtual void setLevel(int level) {}

	/**
	*	Paints one vertex
	* @param renderingParams pointer to rendering params with hints as to what to color
	* @param pFacet pointer to facet in which vertex resides
	* @param pHalfedge halfedge which points to vertex
	* @param color pointer to location where to place color
	* @return true if successful
	*/
	virtual bool paintVertexInFacet(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Facet_const_handle& pFacet,
		const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
		GLubyte* color)
	{
		return false;
	}
};

class EdgePainter : public Painter {
public:
	/** 
	* selects color for edges 		
	* @param pHalfedge const handle to halfedge
	* @param color array of 4 integers in which color is filled
	* @return 0 on failure
	*/
	virtual int ecolor(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
		GLubyte* color) = 0;	

	/**
	 *	checks how many edges needed to render
	 */
	virtual int visibleEdges() const = 0;
};

/**
 *	used for showing debug information
 */
class WorldObject;

class DebugPainter : public Painter {
public:
	/**
	 *	paints relevant debug information, responsible for all drawing operations
	 *  @param mesh pointer to mesh structure
	 */
	virtual void paint(Mesh* mesh, MeshSelection* meshSelection, MeshSurfaceGraph* msg, WorldObject* wo)
	{
		paint(mesh, meshSelection, msg);
	}

	virtual void paint(Mesh* mesh, MeshSelection* meshSelection, MeshSurfaceGraph* msg) {}

	virtual void postPaint(Mesh* mesh, MeshSelection* MeshSelection, MeshSurfaceGraph* msg) {}

	//virtual void paint(Mesh* cage_mesh) {}

	/**
	 *	tells debug painters if certain mesh is invalidated (perhaps deleted), if NULL is passed then
	 * all data is invalidated
	 */
	virtual void invalidate(Mesh* mesh) {}

	/**
	 *	give debug renderer a command and get back a string signalling the new mode
	 */
	virtual QString command(const int key) {
		return QString::null;
	}
};

/**
 *	Do not display debug information
 */ 
class NoneDebugPainter : public DebugPainter {
public:
	virtual QString name() {return "None";}

	virtual void paint(Mesh* mesh, MeshSelection* meshSelection, MeshSurfaceGraph* msg) {
		return;
	}
};

class SolidFacetPainter : public FacetPainter {
	virtual QString name() {return "Solid";}

	virtual bool isSmoothShading() {
		return true;
	}

	virtual unsigned int requiredBufferSize(Mesh* mesh) {
		return mesh->size_of_vertices() * 4;
	}

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer)
	{
		QColor c = renderingParams->m_facetColor;
		GLubyte color[4];
		color[0] = c.red(); color[1] = c.green(); color[2] = c.blue(); color[3] = renderingParams->m_alpha;

		const int vsize = mesh->size_of_vertices();
		int dataIndex = 0;
		for (int i=0; i<vsize; i++) {
			memcpy((void*)&colorBuffer[i * 4], (void*)color, sizeof(GLubyte) * 4);
		}

		return true;
	}

	virtual bool paintVertex(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Vertex_const_handle& pVertex,
		GLubyte* color)
	{
		QColor c = renderingParams->m_facetColor;
		color[0] = c.red(); color[1] = c.green(); color[2] = c.blue(); color[3] = 255;
		return true;
	}

	
	virtual bool paintVertexInFacet(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Facet_const_handle& pFacet,
		const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
		GLubyte* color)
	{
		QColor c = renderingParams->m_facetColor;
		color[0] = c.red(); color[1] = c.green(); color[2] = c.blue(); color[3] = 255;
		return true;
	}		
};

/**
 *	Class used to draw overlay of mesh in a small window, used to show some extra information
 */
class OverlayPainter : public Painter {
public:
	/**
	 *	assume orthogonal 2d surface is setup, just render what you need...
	 */
	virtual void render() = 0;
};

class MeshRenderer : public QObject 
{
	Q_OBJECT
public:
	MeshRenderer(QObject* parent = 0, const char* name = 0);
	~MeshRenderer();

private:
	Mesh* m_mesh;

	MeshSelection* m_meshSelection;

	DrawMethod m_dmet;

	qglviewer::ManipulatedFrame* m_manipulatedFrame;

	FacetPainter* m_selectedFacetPainter;

	std::vector<EdgePainter*> m_edgePainters;
	int m_selectedEdgePainter;

	/** scale for normals */
	GLfloat m_scale;
	GLfloat m_oldScale;

	int m_sizeVertices;
	int m_sizeFacets;
	int m_sizeEdges;

	/** holds coordinates of mesh vertices */
	GLfloat* m_vertexBuffer;

	/** holds color RGB values for each vertex */
	GLubyte* m_vertexColorBuffer;

	/** size of current color buffer */
	unsigned int m_vertexColorBufferSize;
	bool m_vertexColorBufferVBO;

	/** holds normals for each vertex */
	GLfloat* m_vertexNormalBuffer;

	/** holds indices for vertices */
	GLuint* m_vertexIndexesBuffer;

	/** holds indices for mesh edges */
	GLuint* m_edgeIndexesBuffer;

	/** pointer to vertex buffer object */
	GLuint m_vertexBufferPtr;
	/** pointer to vertex normals buffer object */
	GLuint m_vertexNormalBufferPtr;
	/** pointer to vertex indexes buffer object */
	GLuint m_vertexIndexesBufferPtr;
	/** pointer to vertex color buffer object */
	GLuint m_vertexColorBufferPtr;	
	/** id of mesh display list */
	GLuint m_displayListId;
	/** id of edges display list */
	GLuint m_edgesDisplayListId;
	/** id of vertices display list */
	GLuint m_verticesDisplayListId;

	/** holds normals for the facets */	
	GLfloat* m_facetNormalBuffer;
	/** pointer to facet normals for buffer object */
	GLuint m_facetNormalBufferPtr;
	/** duplicated vertices array */
	GLfloat* m_facetVerticesBuffer;
	/** pointer to facet vertices for buffer object */
	GLuint m_facetVerticesBufferPtr;

	/** pointer to struct of rendering parameters */
	RenderingParams* m_renderingParams;

protected:
	/** returns pointer to currently selected facet painter */
	FacetPainter* getSelectedFacetPainter();
	/** returns pointer to currently selected edge painter */
	EdgePainter* getSelectedEdgePainter();

	/** fills vertex coordinates and normals buffer for currently loaded mesh */
	void fillVertexBuffer();

	/** if VBO supported bind the current vertex buffer to the GPU */
	void bindVertexBuffer();
	void unbindVertexBuffer();

	void bindVertexColorBuffer();
	void unbindVertexColorBuffer();

	/** fills vertex color bufferfor this frame of rendering, one color for each vertex */
	void fillVertexColorBuffer();

	/** fills array with indexes of mesh edges */
	void fillEdgeBuffer();

	/** selects color for each vertex in each facet */
	//void fillVertexColorBuffer2();

	/** render using vertex arrays or vertex buffer objects */
	void renderBuffers();

	/** render flat shading using VBO's */
	void renderBuffersFlatShaded();

	/** render using display lists */
	void renderDisplayLists();

	/** render brute force */
	void renderBrute();

	/** render brute force one color per vertex */
	void renderBruteSmooth();

	/** render brute force one color per vertex per facet */
	void renderBruteFlat(const bool selection = false);

	/** render brute force with glnames for the facets */
	void renderBruteWithNames();

	/** do super impose of vertices using a display list */
	void renderVertices();

	/** do super impose of vertices */
	void renderVerticesBrute();

	/** do super impose of edges using a display list*/
	void renderEdges();

	/** do supoer impose of edges */
	void renderEdgesBrute();

	void setGeneralRenderingMode();

	void checkGlErrors();

	void reloadNormals();

public:

	/** clears all memory and buffers, resets gl buffers...*/
	void Cleanup();

	/** 
	 * sets a new mesh in the renderer, clears old mesh and buffers 
	 *
	 * @param mesh pointer to mesh
	 * @param meshSelection information about if the mesh is selected and what within it
	 */
	void setMesh(Mesh* mesh, MeshSurfaceGraph* msg, MeshSelection* meshSelection);
	DrawMethod drawMethod() { return m_dmet; }

	void setCage(Mesh* cage, MeshSurfaceGraph* msg, MeshSelection* meshSelection);

	/** sets rendering parameters for this class */
	void setRenderingParams(RenderingParams* renderingParams);

	const GLfloat& scale() const {return m_scale;}
	//GLfloat& scale() {return m_scale;}
	void scale(const float f);

	/**
	*	adds pointer to an edge painter
	* @param edgePainter pointer to the painter
	* @return number of current painters
	*/
	int addEdgePainter(EdgePainter* edgePainter);	

	/** returns list of available edge painters */
	QStringList getEdgePainterNames();

	/** selects edge painter and returns selected index */
	int selectEdgePainter(const int edgePainterIndex);
	/** selects edge painter and returns selected index */
	int selectEdgePainter(const QString& edgePainterName);

	/** 
	 * renders mesh using current rendering parameters 
	 *
	 * @param selectionMode if other than none draw with glNames so we can select specific parts
	 */
	void render(int selectionMode = SELECTION_MODE_NONE, int part = -1);

	/** return pointer to loaded mesh */
	Mesh* mesh() {return m_mesh;}

	qglviewer::ManipulatedFrame* frame() {return m_manipulatedFrame;}

public slots:
	/** need to rethink how to render class */
	void OnObjectChanged(bool geometryChanged);
	void OnObjectChanged();

	void OnFacetPainterChanged(FacetPainter* facetPainter);
};


class CageRenderer : public QObject 
{
	Q_OBJECT
public:
	CageRenderer(QObject* parent = 0, const char* name = 0);
	~CageRenderer();

private:
	Mesh* m_cage;

	MeshSelection* m_cageSelection;

	DrawMethod m_dmet;

	qglviewer::ManipulatedFrame* m_manipulatedFrame;

	FacetPainter* m_selectedFacetPainter;

	std::vector<EdgePainter*> m_edgePainters;
	int m_selectedEdgePainter;

	/** scale for normals */
	GLfloat m_scale;
	GLfloat m_oldScale;

	int m_sizeVertices;
	int m_sizeFacets;
	int m_sizeEdges;

	/** holds coordinates of mesh vertices */
	GLfloat* m_vertexBuffer;

	/** holds color RGB values for each vertex */
	GLubyte* m_vertexColorBuffer;

	/** size of current color buffer */
	unsigned int m_vertexColorBufferSize;
	bool m_vertexColorBufferVBO;

	/** holds normals for each vertex */
	GLfloat* m_vertexNormalBuffer;

	/** holds indices for vertices */
	GLuint* m_vertexIndexesBuffer;

	/** holds indices for mesh edges */
	GLuint* m_edgeIndexesBuffer;

	/** pointer to vertex buffer object */
	GLuint m_vertexBufferPtr;
	/** pointer to vertex normals buffer object */
	GLuint m_vertexNormalBufferPtr;
	/** pointer to vertex indexes buffer object */
	GLuint m_vertexIndexesBufferPtr;
	/** pointer to vertex color buffer object */
	GLuint m_vertexColorBufferPtr;	
	/** id of mesh display list */
	GLuint m_displayListId;
	/** id of edges display list */
	GLuint m_edgesDisplayListId;
	/** id of vertices display list */
	GLuint m_verticesDisplayListId;

	/** holds normals for the facets */	
	GLfloat* m_facetNormalBuffer;
	/** pointer to facet normals for buffer object */
	GLuint m_facetNormalBufferPtr;
	/** duplicated vertices array */
	GLfloat* m_facetVerticesBuffer;
	/** pointer to facet vertices for buffer object */
	GLuint m_facetVerticesBufferPtr;

	/** pointer to struct of rendering parameters */
	RenderingParams* m_renderingParams;

protected:
	/** returns pointer to currently selected facet painter */
	FacetPainter* getSelectedFacetPainter();
	/** returns pointer to currently selected edge painter */
	EdgePainter* getSelectedEdgePainter();

	/** fills vertex coordinates and normals buffer for currently loaded mesh */
	void fillVertexBuffer();

	/** if VBO supported bind the current vertex buffer to the GPU */
	void bindVertexBuffer();
	void unbindVertexBuffer();

	void bindVertexColorBuffer();
	void unbindVertexColorBuffer();

	/** fills vertex color bufferfor this frame of rendering, one color for each vertex */
	void fillVertexColorBuffer();

	/** fills array with indexes of mesh edges */
	void fillEdgeBuffer();

	/** selects color for each vertex in each facet */
	void fillVertexColorBuffer2();

	/** render using vertex arrays or vertex buffer objects */
	void renderBuffers();

	/** render flat shading using VBO's */
	void renderBuffersFlatShaded();

	/** render using display lists */
	void renderDisplayLists();

	/** render brute force */
	void renderBrute();

	/** render brute force one color per vertex */
	void renderBruteSmooth();

	/** render brute force one color per vertex per facet */
	void renderBruteFlat(const bool selection = false);

	/** render brute force with glnames for the facets */
	void renderBruteWithNames();

	/** do super impose of vertices using a display list */
	void renderVertices();

	/** do super impose of vertices */
	void renderVerticesBrute();

	/** do super impose of edges using a display list*/
	void renderEdges();

	/** do super impose of edges */
	void renderEdgesBrute();

	void setGeneralRenderingMode();

	void checkGlErrors();

	void reloadNormals();

public:

	/** clears all memory and buffers, resets gl buffers...*/
	void Cleanup();

	/** 
	 * sets a new mesh in the renderer, clears old mesh and buffers 
	 *
	 * @param mesh pointer to mesh
	 * @param meshSelection information about if the mesh is selected and what within it
	 */
	void setCage(Mesh* cage, MeshSurfaceGraph* msg, MeshSelection* meshSelection);
	DrawMethod drawMethod() { return m_dmet; }
	

	/** sets rendering parameters for this class */
	void setRenderingParams(RenderingParams* renderingParams);

	const GLfloat& scale() const {return m_scale;}
	//GLfloat& scale() {return m_scale;}
	void scale(const float f);

	/**
	*	adds pointer to an edge painter
	* @param edgePainter pointer to the painter
	* @return number of current painters
	*/
	int addEdgePainter(EdgePainter* edgePainter);	

	/** returns list of available edge painters */
	QStringList getEdgePainterNames();

	/** selects edge painter and returns selected index */
	int selectEdgePainter(const int edgePainterIndex);
	/** selects edge painter and returns selected index */
	int selectEdgePainter(const QString& edgePainterName);

	/** 
	 * renders mesh using current rendering parameters 
	 *
	 * @param selectionMode if other than none draw with glNames so we can select specific parts
	 */
	void render(int selectionMode = SELECTION_MODE_NONE, int part = -1);

	/** return pointer to loaded mesh */
	Mesh* cage() {return m_cage;}

	qglviewer::ManipulatedFrame* frame() {return m_manipulatedFrame;}

public slots:
	/** need to rethink how to render class */
	void OnObjectChanged(bool geometryChanged);
	void OnObjectChanged();

	void OnFacetPainterChanged(FacetPainter* facetPainter);
};


#endif