#include "stdafx.h"

#include "MeshRenderer.h"
//#include "MeshPartition.h"

//using namespace MR;

FacetPainter* MeshRenderer::getSelectedFacetPainter()
{
	return m_selectedFacetPainter;
}

EdgePainter* MeshRenderer::getSelectedEdgePainter()
{
	assert(m_selectedEdgePainter>=0 && m_selectedEdgePainter<m_edgePainters.size());

	return m_edgePainters[m_selectedEdgePainter];
}

MeshRenderer::MeshRenderer(QObject* parent, const char* name) :
	QObject(parent, name),
	m_mesh(NULL),
	m_selectedFacetPainter(NULL),
	m_selectedEdgePainter(-1),
	m_vertexBuffer(NULL),
	m_vertexColorBuffer(NULL),
	m_vertexNormalBuffer(NULL),
	m_facetNormalBuffer(NULL),
	m_facetVerticesBuffer(NULL),
	m_vertexIndexesBuffer(NULL),
	m_edgeIndexesBuffer(NULL),
	m_renderingParams(NULL),
	m_vertexColorBufferVBO(false),
	m_vertexColorBufferSize(0),
	m_displayListId(0),
	m_edgesDisplayListId(0),
	m_manipulatedFrame(NULL),
	m_oldScale(1.0),
	m_scale(1.0)
{
	m_manipulatedFrame = new qglviewer::ManipulatedFrame;
	m_manipulatedFrame->setReferenceFrame(NULL);
}

MeshRenderer::~MeshRenderer()
{
	Cleanup();

	for (int i=0; i<m_edgePainters.size(); i++)
	{
		delete m_edgePainters[i];
	}

	m_edgePainters.clear();

	delete m_manipulatedFrame;
}

void MeshRenderer::Cleanup()
{
	//we can now get rid of the local buffers
	delete[] m_vertexBuffer;
	m_vertexBuffer = NULL;
	delete[] m_vertexNormalBuffer;
	m_vertexNormalBuffer = NULL;
	delete[] m_vertexIndexesBuffer;
	m_vertexIndexesBuffer = NULL;

	delete[] m_facetNormalBuffer;
	m_facetNormalBuffer = NULL;
	delete[] m_facetVerticesBuffer;
	m_facetVerticesBuffer = NULL;

	delete[] m_vertexColorBuffer;
	m_vertexColorBuffer = NULL;

	delete[] m_edgeIndexesBuffer;
	m_edgeIndexesBuffer = NULL;

	//TODO: if created display list remove it
	if (m_displayListId != 0)
	{
		glDeleteLists(m_displayListId, 1);
		m_displayListId = 0;
	}

	if (m_edgesDisplayListId != 0)
	{
		glDeleteLists(m_edgesDisplayListId, 1);
		m_edgesDisplayListId = 0;
	}
}

void MeshRenderer::scale(const float f)
{
	if (f != m_scale)
	{
		m_oldScale = m_scale;
		m_scale = f;
		//reloadNormals();
	}
}


void MeshRenderer::setMesh(Mesh* mesh, MeshSurfaceGraph* msg, MeshSelection* meshSelection)
{
	Cleanup();

	m_mesh = mesh;
	m_meshSelection = meshSelection;

	m_sizeFacets = m_mesh->size_of_facets();
	m_sizeVertices = m_mesh->size_of_vertices();
	m_sizeEdges = m_mesh->size_of_halfedges() / 2;

	//set vertex buffer so each vertex in each facet has 3 coordinates
	m_vertexBuffer = new GLfloat[m_sizeVertices * 3];
	//same for normals
	m_vertexNormalBuffer = new GLfloat[m_sizeVertices * 3];
	//save indexes
	m_vertexIndexesBuffer = new GLuint[m_sizeFacets * 3];
	//save RGBA values for color values
	m_vertexColorBuffer = new GLubyte[m_sizeVertices * 4];

	m_edgeIndexesBuffer = new GLuint[m_mesh->size_of_halfedges()];

	m_facetNormalBuffer = new GLfloat[m_sizeFacets * 9];
	m_facetVerticesBuffer = new GLfloat[m_sizeFacets * 9];

	checkGlErrors();

	//fill vertex and normal buffers
	fillVertexBuffer();

	//fill color information
	fillVertexColorBuffer();

	//fill edge information
	fillEdgeBuffer();

	checkGlErrors();

	//bind vertex and normal buffers as VBO's (vertex buffer objects)
	if (m_renderingParams && m_renderingParams->m_vboSupported)
	{
		bindVertexBuffer();
		//fill initial color information
		bindVertexColorBuffer();
	}
}

void MeshRenderer::reloadNormals()
{
	GLfloat inv_scale = 1.0 / m_scale;

	for (int i=0; i<m_sizeVertices * 3; i++)
	{
		m_vertexNormalBuffer[i] = m_vertexNormalBuffer[i] * m_oldScale * inv_scale;
	}

	for (int i=0; i<m_sizeFacets * 9; i++)
	{
		m_facetNormalBuffer[i] = m_facetNormalBuffer[i] * m_oldScale * inv_scale;
	}

	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexNormalBufferPtr);

	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeVertices * 3, m_vertexNormalBuffer, GL_STATIC_DRAW_ARB);

	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_facetNormalBufferPtr);

	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeFacets * 9, m_facetNormalBuffer, GL_STATIC_DRAW_ARB);
}


void MeshRenderer::setRenderingParams(RenderingParams* renderingParams)
{
	m_renderingParams = renderingParams;

	//m_renderingParams->m_vboSupported = false;
}

int MeshRenderer::addEdgePainter(EdgePainter* edgePainter)
{
	m_edgePainters.push_back(edgePainter);
	return m_edgePainters.size();
}

QStringList MeshRenderer::getEdgePainterNames()
{
	QStringList epnames;
	for (int i=0;i<m_edgePainters.size();i++)
	{
		if (m_edgePainters[i] != NULL) {
			epnames.append(m_edgePainters[i]->name());
		}
	}

	return epnames;
}

int MeshRenderer::selectEdgePainter(const int edgePainterIndex)
{
	m_selectedEdgePainter = edgePainterIndex;

	return m_selectedEdgePainter;
}

int MeshRenderer::selectEdgePainter(const QString& edgePainterName)
{
	for (int i=0;i<m_edgePainters.size();i++)
	{
		if (m_edgePainters[i] != NULL && m_edgePainters[i]->name() == edgePainterName) {
			m_selectedEdgePainter = i;
			return m_selectedEdgePainter;
		}
	}

	return -1;
}

void MeshRenderer::setGeneralRenderingMode()
{
	// shading option
	if(m_renderingParams->m_smoothShading)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	// culling option
	if(m_renderingParams->m_culling)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);

}

void MeshRenderer::render(int selectionMode, int part)
{
	setGeneralRenderingMode();

	//draw model
	glPushMatrix();

	if (selectionMode == SELECTION_MODE_OBJECT)
	{
		//cout << "SELECTION_MODE_OBJECT" << endl;

		glEnable(GL_DEPTH_TEST);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDisable(GL_LIGHTING);
		renderBruteFlat(true);
	}
	else if (selectionMode == SELECTION_MODE_FACET)
	{
		//cout << "SELECTION_MODE_FACET" << endl;

		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDisable(GL_LIGHTING);
		renderBruteWithNames();
	}
	else if (selectionMode == SELECTION_MODE_VERTEX)
	{
		//cout << "SELECTION_MODE_VERTEX" << endl;
	}
	else
	{
		//cout << "SELECTION_MODE_ELSE" << endl;

		//regular rendering

		if (m_renderingParams->m_superimposeVertices)
		{
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_BLEND);
			// enable polygon offset
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(3.0f, 1.0f);
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);

			glPointSize((float)m_renderingParams->m_vertexThickness);
			renderVertices();
		}

		if (m_renderingParams->m_superimposeEdges)
		{
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_BLEND);
			// enable polygon offset
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(3.0f, 1.0f);
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

			glLineWidth((float)m_renderingParams->m_edgeThickness);
			renderEdges();
		}

		// polygon mode (point, line or fill)
		glPolygonMode(GL_FRONT_AND_BACK, m_renderingParams->m_polygonMode);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		// antialiasing
		if(m_renderingParams->m_antialiasing)
		{
			glEnable(GL_LINE_SMOOTH);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
			glLineWidth(1.5f);
		}
		else
		{
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_BLEND);
			glLineWidth(1.0f);
		}

		// lighting option
		if(m_renderingParams->m_lighting)
		{
			m_renderingParams->m_useNormals = true;
			glEnable(GL_LIGHTING);
		}
		else
		{
			m_renderingParams->m_useNormals = false;
			glDisable(GL_LIGHTING);
		}

		glDisable(GL_POLYGON_OFFSET_LINE);

		///////////////////////////////////////////
		//renderBuffersDebug();
		if (getSelectedFacetPainter()->isSmoothShading())
		{
			renderBuffers();
		}
		else
		{
			if (m_renderingParams->m_vboSupported)
				renderBuffersFlatShaded();
			else
				renderDisplayLists();
		}
		///////////////////////////////////////////
	}

	glPopMatrix();
}

void MeshRenderer::fillEdgeBuffer()
{
	assert(m_edgeIndexesBuffer);

	Mesh::Edge_const_iterator hit = m_mesh->edges_begin();
	Mesh::Edge_const_iterator hit_end = m_mesh->edges_end();

	int dataIndex = 0;
	for(; hit != hit_end; hit++)
	{
		m_edgeIndexesBuffer[dataIndex	] = hit->vertex()->index();
		m_edgeIndexesBuffer[dataIndex+1	] = hit->opposite()->vertex()->index();
		dataIndex += 2;
	}
}

void MeshRenderer::renderVertices()
{
	// if list not ready create it
	if (m_verticesDisplayListId == 0)
	{
		m_verticesDisplayListId = glGenLists(1);

		// Render model into the display list
		glNewList(m_verticesDisplayListId, GL_COMPILE);
		{
			renderVerticesBrute();
		}
		// Finished Creating the Display List
		glEndList();
	}

	// call the display list
	glCallList( m_verticesDisplayListId );
}

void MeshRenderer::renderVerticesBrute()
{
	glBegin(GL_POINTS);

	glColor4ub(
		m_renderingParams->m_vertexColor.red(),
		m_renderingParams->m_vertexColor.green(),
		m_renderingParams->m_vertexColor.blue(),
		255);

	const int numVertices = m_mesh->size_of_vertices();
	for (int i=0; i < numVertices; i++)
	{
		glVertex3f(
			m_vertexBuffer[i*3],
			m_vertexBuffer[i*3+1],
			m_vertexBuffer[i*3+2]);
	}

	glEnd();
}

void MeshRenderer::renderEdges()
{
	//if list not ready create it
	if (m_edgesDisplayListId == 0)
	{
		m_edgesDisplayListId = glGenLists(1);

		// Render model into the display list
		glNewList(m_edgesDisplayListId, GL_COMPILE);
		{
			renderEdgesBrute();
		}
		// Finished Creating the Display List
		glEndList();
	}

	//call the display list
	glCallList( m_edgesDisplayListId );
}

void MeshRenderer::renderEdgesBrute()
{
	glBegin(GL_LINES);

	glColor4ub(
		m_renderingParams->m_edgeColor.red(),
		m_renderingParams->m_edgeColor.green(),
		m_renderingParams->m_edgeColor.blue(),
		255);

	const int numEdges = m_mesh->size_of_halfedges() / 2;
	for (int i=0; i < numEdges; i++)
	{
		glVertex3f(
			m_vertexBuffer[m_edgeIndexesBuffer[i*2]*3],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2]*3+1],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2]*3+2]);

		glVertex3f(
			m_vertexBuffer[m_edgeIndexesBuffer[i*2+1]*3],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2+1]*3+1],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2+1]*3+2]);
	}

	glEnd();
}

void MeshRenderer::renderBrute()
{
	if (getSelectedFacetPainter()->isSmoothShading())
		renderBruteSmooth();
	else
		renderBruteFlat();
}

void MeshRenderer::renderBruteSmooth()
{
	glBegin(GL_TRIANGLES);
	for (int i=0; i<m_sizeFacets; i++)
	{
		for (int j = i*3; j< i*3+3; j++)
		{
			glColor4ub(
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4],
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4+1],
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4+2],
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4+3]);

			glVertex3f(
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+1],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+2]);

			glNormal3f(
				m_vertexNormalBuffer[m_vertexIndexesBuffer[j]*3],
				m_vertexNormalBuffer[m_vertexIndexesBuffer[j]*3+1],
				m_vertexNormalBuffer[m_vertexIndexesBuffer[j]*3+2]);
		}
	}
	glEnd();
}

void MeshRenderer::renderBruteWithNames()
{
	int sizeOfFacets = m_sizeFacets;
	for (int i=0; i<sizeOfFacets; i++)
	{
		glPushName(i);
		glBegin(GL_TRIANGLES);
		for (int j = i*3; j< i*3+3; j++)
		{
			glVertex3f(
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+1],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+2]);
			glNormal3f(
				m_facetNormalBuffer[j*3],
				m_facetNormalBuffer[j*3+1],
				m_facetNormalBuffer[j*3+2]);
		}
		glEnd();
		glPopName();
	}

}

void MeshRenderer::renderBruteFlat(const bool selection)
{
	if (m_renderingParams->m_smoothNormals)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	glBegin(GL_TRIANGLES);

	if (selection)
	{
		for (int i=0; i<m_sizeFacets; i++)
		{
			for (int j = i*3; j< i*3+3; j++)
			{
				glVertex3f(
					m_facetVerticesBuffer[j*3],
					m_facetVerticesBuffer[j*3+1],
					m_facetVerticesBuffer[j*3+2]);
				glNormal3f(
					m_facetNormalBuffer[j*3],
					m_facetNormalBuffer[j*3+1],
					m_facetNormalBuffer[j*3+2]);
			}
		}
	}
	else
	{
		for (int i=0; i<m_sizeFacets; i++)
		{
			for (int j = i*3; j< i*3+3; j++)
			{
				glColor4ub(
					m_vertexColorBuffer[j*4],
					m_vertexColorBuffer[j*4+1],
					m_vertexColorBuffer[j*4+2],
					m_vertexColorBuffer[j*4+3]);
				glVertex3f(
					m_facetVerticesBuffer[j*3],
					m_facetVerticesBuffer[j*3+1],
					m_facetVerticesBuffer[j*3+2]);
				glNormal3f(
					m_facetNormalBuffer[j*3],
					m_facetNormalBuffer[j*3+1],
					m_facetNormalBuffer[j*3+2]);
			}
		}
	}

	glEnd();
}

void MeshRenderer::renderDisplayLists()
{
	//if list not ready create it
	if (m_displayListId == 0)
	{
		m_displayListId = glGenLists(1);

		// Render model into the display list
		glNewList(m_displayListId, GL_COMPILE);
		{
			renderBrute();
		}
		// Finished Creating the Display List
		glEndList();
	}

	//call the display list
	glCallList( m_displayListId );
}

void MeshRenderer::checkGlErrors()
{
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR)
	{
		qDebug("GLError = %d", err);
	}
}

void MeshRenderer::renderBuffersFlatShaded()
{
	//assert (m_facetNormalBuffer && m_facetVerticesBuffer && m_vertexColorBuffer);

	if (m_renderingParams->m_smoothNormals)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	if (!m_renderingParams->m_vboSupported) return;

	glEnableClientState( GL_VERTEX_ARRAY );				// Enable Vertex Arrays
	glEnableClientState( GL_NORMAL_ARRAY );				// Enable Normal Arrays
	glEnableClientState( GL_COLOR_ARRAY  );				// Enable Color Arrays

	glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_facetVerticesBufferPtr);
	glVertexPointer( 3, GL_FLOAT, 0, (char *) NULL );		// Set The Vertex Pointer To The Vertex Buffer

	glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_facetNormalBufferPtr);
	glNormalPointer( GL_FLOAT, 0, (char *) NULL );		// Set The normal Pointer To The normals Buffer

	glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexColorBufferPtr);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, (char *) NULL );		// Set The normal Pointer To The color Buffer

	GLint numOfElements = m_sizeFacets * 3;
	//glDrawElements( GL_TRIANGLES, numOfElements, GL_UNSIGNED_INT, (char *)NULL);
	glDrawArrays(GL_TRIANGLES, 0, numOfElements);

	glDisableClientState( GL_VERTEX_ARRAY	);				// Disable Vertex Arrays
	glDisableClientState( GL_NORMAL_ARRAY	);				// Disable Normal Arrays
	glDisableClientState( GL_COLOR_ARRAY	);				// Disable Color Arrays
}

void MeshRenderer::renderBuffers()
{
	assert (m_vertexBuffer && m_vertexNormalBuffer && m_vertexColorBuffer);

	glEnableClientState( GL_VERTEX_ARRAY );				// Enable Vertex Arrays
	glEnableClientState( GL_NORMAL_ARRAY );				// Enable Normal Arrays
	glEnableClientState( GL_COLOR_ARRAY  );				// Enable Color Arrays

	if(m_renderingParams && m_renderingParams->m_vboSupported)
	{
		glEnableClientState( GL_INDEX_ARRAY  );				// Enable Index Arrays

		glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexBufferPtr );
		glVertexPointer( 3, GL_FLOAT, 0, (char *) NULL );		// Set The Vertex Pointer To The Vertex Buffer

		glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexNormalBufferPtr);
		glNormalPointer( GL_FLOAT, 0, (char *) NULL );		// Set The normal Pointer To The normals Buffer

		glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexColorBufferPtr);
		glColorPointer(4, GL_UNSIGNED_BYTE, 0, (char *) NULL );		// Set The normal Pointer To The color Buffer

		glBindBufferARB( GL_ELEMENT_ARRAY_BUFFER_ARB, m_vertexIndexesBufferPtr);
	}
	else
	{
		glVertexPointer( 3, GL_FLOAT, 0, m_vertexBuffer );	// Set The Vertex Pointer To Our Vertex Data
		glNormalPointer(GL_FLOAT, 0, m_vertexNormalBuffer );	// Set The Vertex Pointer To Our TexCoord Data
		glColorPointer(4, GL_UNSIGNED_BYTE, 0, m_vertexColorBuffer); //set color buffer
	}

	GLint numOfElements = m_sizeFacets * 3;

	if (m_renderingParams && m_renderingParams->m_vboSupported)
	{
		glDrawElements( GL_TRIANGLES, numOfElements, GL_UNSIGNED_INT, (char *)NULL);
		glDisableClientState( GL_INDEX_ARRAY);				// Disable Index Arrays
	}
	else
	{
		glDrawElements( GL_TRIANGLES, numOfElements, GL_UNSIGNED_INT, m_vertexIndexesBuffer);
	}

	glDisableClientState( GL_VERTEX_ARRAY	);				// Disable Vertex Arrays
	glDisableClientState( GL_NORMAL_ARRAY	);				// Disable Normal Arrays
	glDisableClientState( GL_COLOR_ARRAY	);				// Disable Color Arrays
}

void MeshRenderer::fillVertexBuffer()
{
	assert(m_mesh && m_vertexBuffer && m_vertexNormalBuffer);

	//I assume buffer is allocated already (when loading new mesh)

	int dataIndex = 0;

	//Go over vertices
	Mesh::Vertex_const_iterator vit = m_mesh->vertices_begin();
	Mesh::Vertex_const_iterator vit_end = m_mesh->vertices_end();

	//GLfloat inv_scale = 1 / m_scale;

	for (;vit != vit_end; vit++)
	{
		Point_3 p1 = vit->point();
		Vector_3 n1 = vit->normal();

		m_vertexBuffer[dataIndex] = p1[0]; 
		m_vertexBuffer[dataIndex+1] = p1[1];
		m_vertexBuffer[dataIndex+2] = p1[2];

		m_vertexNormalBuffer[dataIndex] = n1[0] /* * inv_scale*/;
		m_vertexNormalBuffer[dataIndex+1] = n1[1] /** inv_scale */;
		m_vertexNormalBuffer[dataIndex+2] = n1[2] /** inv_scale */;

		dataIndex += 3;
	}

	//Go over facets (fill indexes and facet normals)
	dataIndex = 0;
	Indexer *indr;
	indr = new CountIndexer(m_mesh->size_of_facets());

	while (indr->hasNext())
	{
		int i = indr->next();
		
		Mesh::Facet_handle fit = m_mesh->findFacet(i);
		Mesh::Halfedge_const_handle c = fit->halfedge();

		//first point of triangle
		Mesh::Vertex_const_handle v1 = c->vertex();
		Mesh::Vertex_const_handle v2 = c->next()->vertex();
		Mesh::Vertex_const_handle v3 = c->next()->next()->vertex();

		m_vertexIndexesBuffer[dataIndex] = v1->index();
		m_vertexIndexesBuffer[dataIndex+1] = v2->index();
		m_vertexIndexesBuffer[dataIndex+2] = v3->index();

		if (m_renderingParams->m_smoothNormals)
		{
			/*m_facetNormalBuffer[dataIndex*3	]	= m_vertexNormalBuffer[v1->index()*3];
			m_facetNormalBuffer[dataIndex*3	+1] = m_vertexNormalBuffer[v1->index()*3 + 1];
			m_facetNormalBuffer[dataIndex*3	+2] = m_vertexNormalBuffer[v1->index()*3 + 2 ];
			m_facetNormalBuffer[dataIndex*3	+3]	= m_vertexNormalBuffer[v2->index()*3];
			m_facetNormalBuffer[dataIndex*3	+4] = m_vertexNormalBuffer[v2->index()*3 + 1];
			m_facetNormalBuffer[dataIndex*3	+5] = m_vertexNormalBuffer[v2->index()*3 + 2];
			m_facetNormalBuffer[dataIndex*3	+6]	= m_vertexNormalBuffer[v3->index()*3];
			m_facetNormalBuffer[dataIndex*3	+7] = m_vertexNormalBuffer[v3->index()*3 + 1];
			m_facetNormalBuffer[dataIndex*3	+8] = m_vertexNormalBuffer[v3->index()*3 + 2];*/
			m_facetNormalBuffer[dataIndex*3	]	= v1->normal().x();
			m_facetNormalBuffer[dataIndex*3	+1] = v1->normal().y();
			m_facetNormalBuffer[dataIndex*3	+2] = v1->normal().z();
			m_facetNormalBuffer[dataIndex*3	+3]	= v2->normal().x();
			m_facetNormalBuffer[dataIndex*3	+4] = v2->normal().y();
			m_facetNormalBuffer[dataIndex*3	+5] = v2->normal().z();
			m_facetNormalBuffer[dataIndex*3	+6]	= v3->normal().x();
			m_facetNormalBuffer[dataIndex*3	+7] = v3->normal().y();
			m_facetNormalBuffer[dataIndex*3	+8] = v3->normal().z();
		}
		else
		{
			m_facetNormalBuffer[dataIndex*3	]	= fit->normal().x();
			m_facetNormalBuffer[dataIndex*3	+1] = fit->normal().y();
			m_facetNormalBuffer[dataIndex*3	+2] = fit->normal().z();
			m_facetNormalBuffer[dataIndex*3	+3]	= fit->normal().x();
			m_facetNormalBuffer[dataIndex*3	+4] = fit->normal().y();
			m_facetNormalBuffer[dataIndex*3	+5] = fit->normal().z();
			m_facetNormalBuffer[dataIndex*3	+6]	= fit->normal().x();
			m_facetNormalBuffer[dataIndex*3	+7] = fit->normal().y();
			m_facetNormalBuffer[dataIndex*3	+8] = fit->normal().z();
		}

		m_facetVerticesBuffer[dataIndex*3	]	= v1->point().x();
		m_facetVerticesBuffer[dataIndex*3	+1] = v1->point().y();
		m_facetVerticesBuffer[dataIndex*3	+2] = v1->point().z();
		m_facetVerticesBuffer[dataIndex*3	+3]	= v2->point().x();
		m_facetVerticesBuffer[dataIndex*3	+4] = v2->point().y();
		m_facetVerticesBuffer[dataIndex*3	+5] = v2->point().z();
		m_facetVerticesBuffer[dataIndex*3	+6]	= v3->point().x();
		m_facetVerticesBuffer[dataIndex*3	+7] = v3->point().y();
		m_facetVerticesBuffer[dataIndex*3	+8] = v3->point().z();
		//}

		dataIndex += 3;
	}

	delete indr;
}

void MeshRenderer::bindVertexColorBuffer()
{
	if (!m_renderingParams || !m_renderingParams->m_vboSupported) return;

	//if not smooth shading we won't use VBO's
	//if (!getSelectedFacetPainter()->isSmoothShading()) return;

	m_vertexColorBufferVBO = true;

	int bufferSize = getSelectedFacetPainter()->requiredBufferSize(m_mesh);

	checkGlErrors();

	//bind color buffer to graphics card
	glGenBuffersARB(1, &m_vertexColorBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexColorBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLubyte) * bufferSize, m_vertexColorBuffer,
		GL_STATIC_DRAW_ARB);


	checkGlErrors();
}

void MeshRenderer::unbindVertexColorBuffer()
{
	if (!m_vertexColorBufferVBO) return;

	glDeleteBuffersARB(1, &m_vertexColorBufferPtr);

	m_vertexColorBufferVBO = false;
}

void MeshRenderer::bindVertexBuffer()
{
	if (!m_renderingParams || !m_renderingParams->m_vboSupported) return;

	checkGlErrors();

	//bind buffer to graphics card
	glGenBuffersARB(1, &m_vertexBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeVertices * 3, m_vertexBuffer,
		GL_STATIC_DRAW_ARB);

	//bind normals buffer to graphics card
	glGenBuffersARB(1, &m_vertexNormalBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexNormalBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeVertices * 3, m_vertexNormalBuffer,
		GL_STATIC_DRAW_ARB);


	//bind indexes buffer to graphics card
	glGenBuffersARB(1, &m_vertexIndexesBufferPtr);
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, m_vertexIndexesBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, sizeof(GLuint) * m_sizeFacets * 3, m_vertexIndexesBuffer,
		GL_STATIC_DRAW_ARB);

	//////////////////////////////////////////////////////////////////////////
	checkGlErrors();

	//bind second buffer of vertices and normals (normals equal)
	glGenBuffersARB(1, &m_facetVerticesBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_facetVerticesBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeFacets * 9, m_facetVerticesBuffer,
		GL_STATIC_DRAW_ARB);

	//bind normals buffer to graphics card
	glGenBuffersARB(1, &m_facetNormalBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_facetNormalBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeFacets * 9, m_facetNormalBuffer,
		GL_STATIC_DRAW_ARB);

	checkGlErrors();

}

void MeshRenderer::unbindVertexBuffer()
{
	glDeleteBuffersARB(1, &m_vertexBufferPtr);
	glDeleteBuffersARB(1, &m_vertexNormalBufferPtr);
	glDeleteBuffersARB(1, &m_vertexIndexesBufferPtr);

	glDeleteBuffersARB(1, &m_facetVerticesBufferPtr);
	glDeleteBuffersARB(1, &m_facetNormalBufferPtr);
}

void MeshRenderer::fillVertexColorBuffer()
{
	assert(m_mesh &&
		m_vertexColorBuffer &&
		m_selectedFacetPainter);

	FacetPainter* fp = getSelectedFacetPainter();

	assert(fp);

	if (m_vertexColorBuffer && m_vertexColorBufferSize != fp->requiredBufferSize(m_mesh))
	{
		//delete the old
		delete[] m_vertexColorBuffer;
		m_vertexColorBuffer = NULL;
		//in with the new
		m_vertexColorBuffer = new GLubyte[fp->requiredBufferSize(m_mesh)];
	}
	else if (!m_vertexColorBuffer)
	{
		//just create a new one
		m_vertexColorBuffer = new GLubyte[fp->requiredBufferSize(m_mesh)];
	}

	fp->fillBuffer(m_renderingParams, m_mesh, NULL, m_dmet, NULL, m_meshSelection, m_vertexColorBuffer);

}

void MeshRenderer::OnObjectChanged()
{
	OnObjectChanged(false);
}

void MeshRenderer::OnObjectChanged(bool geometryChanged)
{
	if (geometryChanged)
	{
		setMesh(m_mesh, NULL, m_meshSelection);
		return;
	}

	if (m_renderingParams && m_renderingParams->m_vboSupported && m_vertexColorBufferVBO)
	{
		unbindVertexColorBuffer();
	}

	fillVertexColorBuffer();

	if (m_renderingParams && m_renderingParams->m_vboSupported)
	{
		bindVertexColorBuffer();
	}

	// display list is no longer relevant, delete it
	if (m_displayListId != 0)
	{
		glDeleteLists(m_displayListId, 1);
		m_displayListId = 0;
	}

	if (m_edgesDisplayListId != 0)
	{
		glDeleteLists(m_edgesDisplayListId, 1);
		m_edgesDisplayListId = 0;
	}
}

void MeshRenderer::OnFacetPainterChanged(FacetPainter* facetPainter)
{
	if (m_selectedFacetPainter == facetPainter) return;

	m_selectedFacetPainter = facetPainter;

	if (m_mesh)
		OnObjectChanged();
}

///////////////////////////////////////////////////////////////////////////////

//#include "CageRenderer.h"

FacetPainter* CageRenderer::getSelectedFacetPainter()
{
	return m_selectedFacetPainter;
}

EdgePainter* CageRenderer::getSelectedEdgePainter()
{
	assert(m_selectedEdgePainter>=0 && m_selectedEdgePainter<m_edgePainters.size());

	return m_edgePainters[m_selectedEdgePainter];
}

CageRenderer::CageRenderer(QObject* parent, const char* name) :
	QObject(parent, name),
	m_cage(NULL),
	m_selectedFacetPainter(NULL),
	m_selectedEdgePainter(-1),
	m_vertexBuffer(NULL),
	m_vertexColorBuffer(NULL),
	m_vertexNormalBuffer(NULL),
	m_facetNormalBuffer(NULL),
	m_facetVerticesBuffer(NULL),
	m_vertexIndexesBuffer(NULL),
	m_edgeIndexesBuffer(NULL),
	m_renderingParams(NULL),
	m_vertexColorBufferVBO(false),
	m_vertexColorBufferSize(0),
	m_displayListId(0),
	m_edgesDisplayListId(0),
	m_manipulatedFrame(NULL),
	m_oldScale(1.0),
	m_scale(1.0)
{
	m_manipulatedFrame = new qglviewer::ManipulatedFrame;
	m_manipulatedFrame->setReferenceFrame(NULL);
}

CageRenderer::~CageRenderer()
{
	Cleanup();

	for (int i=0; i<m_edgePainters.size(); i++)
	{
		delete m_edgePainters[i];
	}

	m_edgePainters.clear();

	delete m_manipulatedFrame;
}

void CageRenderer::Cleanup()
{
	//we can now get rid of the local buffers
	delete[] m_vertexBuffer;
	m_vertexBuffer = NULL;
	delete[] m_vertexNormalBuffer;
	m_vertexNormalBuffer = NULL;
	delete[] m_vertexIndexesBuffer;
	m_vertexIndexesBuffer = NULL;

	delete[] m_facetNormalBuffer;
	m_facetNormalBuffer = NULL;
	delete[] m_facetVerticesBuffer;
	m_facetVerticesBuffer = NULL;

	delete[] m_vertexColorBuffer;
	m_vertexColorBuffer = NULL;

	delete[] m_edgeIndexesBuffer;
	m_edgeIndexesBuffer = NULL;

	//TODO: if created display list remove it
	if (m_displayListId != 0)
	{
		glDeleteLists(m_displayListId, 1);
		m_displayListId = 0;
	}

	if (m_edgesDisplayListId != 0)
	{
		glDeleteLists(m_edgesDisplayListId, 1);
		m_edgesDisplayListId = 0;
	}
}

void CageRenderer::scale(const float f)
{
	if (f != m_scale)
	{
		m_oldScale = m_scale;
		m_scale = f;
		//reloadNormals();
	}
}

void CageRenderer::setCage(Mesh* cage, MeshSurfaceGraph* msg, MeshSelection* cageSelection)
{
	Cleanup();
	
	m_cage = cage;
	m_cageSelection = cageSelection;

	m_sizeFacets = m_cage->size_of_facets();
	m_sizeVertices = m_cage->size_of_vertices();
	m_sizeEdges = m_cage->size_of_halfedges() / 2;

	//set vertex buffer so each vertex in each facet has 3 coordinates
	m_vertexBuffer = new GLfloat[m_sizeVertices * 3];
	//same for normals
	m_vertexNormalBuffer = new GLfloat[m_sizeVertices * 3];
	//save indexes
	m_vertexIndexesBuffer = new GLuint[m_sizeFacets * 3];
	//save RGBA values for color values
	m_vertexColorBuffer = new GLubyte[m_sizeVertices * 4];

	m_edgeIndexesBuffer = new GLuint[m_cage->size_of_halfedges()];

	m_facetNormalBuffer = new GLfloat[m_sizeFacets * 9];
	m_facetVerticesBuffer = new GLfloat[m_sizeFacets * 9];


	checkGlErrors();
	
	//fill vertex and normal buffers
	fillVertexBuffer();

	//fill color information
	fillVertexColorBuffer();

	//fill edge information
	fillEdgeBuffer();

	checkGlErrors();

	//bind vertex and normal buffers as VBO's (vertex buffer objects)
	if (m_renderingParams && m_renderingParams->m_vboSupported)
	{
		bindVertexBuffer();
		//fill initial color information
		bindVertexColorBuffer();
	}
}


void CageRenderer::reloadNormals()
{
	GLfloat inv_scale = 1.0 / m_scale;

	for (int i=0; i<m_sizeVertices * 3; i++)
	{
		m_vertexNormalBuffer[i] = m_vertexNormalBuffer[i] * m_oldScale * inv_scale;
	}

	for (int i=0; i<m_sizeFacets * 9; i++)
	{
		m_facetNormalBuffer[i] = m_facetNormalBuffer[i] * m_oldScale * inv_scale;
	}

	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexNormalBufferPtr);

	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeVertices * 3, m_vertexNormalBuffer, GL_STATIC_DRAW_ARB);

	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_facetNormalBufferPtr);

	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeFacets * 9, m_facetNormalBuffer, GL_STATIC_DRAW_ARB);
}


void CageRenderer::setRenderingParams(RenderingParams* renderingParams)
{
	m_renderingParams = renderingParams;

	//m_renderingParams->m_vboSupported = false;
}

int CageRenderer::addEdgePainter(EdgePainter* edgePainter)
{
	m_edgePainters.push_back(edgePainter);
	return m_edgePainters.size();
}

QStringList CageRenderer::getEdgePainterNames()
{
	QStringList epnames;
	for (int i=0;i<m_edgePainters.size();i++)
	{
		if (m_edgePainters[i] != NULL) {
			epnames.append(m_edgePainters[i]->name());
		}
	}

	return epnames;
}

int CageRenderer::selectEdgePainter(const int edgePainterIndex)
{
	m_selectedEdgePainter = edgePainterIndex;

	return m_selectedEdgePainter;
}

int CageRenderer::selectEdgePainter(const QString& edgePainterName)
{
	for (int i=0;i<m_edgePainters.size();i++)
	{
		if (m_edgePainters[i] != NULL && m_edgePainters[i]->name() == edgePainterName) {
			m_selectedEdgePainter = i;
			return m_selectedEdgePainter;
		}
	}

	return -1;
}

void CageRenderer::setGeneralRenderingMode()
{
	// shading option
	if(m_renderingParams->m_smoothShading)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	// culling option
	if(m_renderingParams->m_culling)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);

}

void CageRenderer::render(int selectionMode, int part)
{
	setGeneralRenderingMode();

	//draw model
	glPushMatrix();

	if (selectionMode == SELECTION_MODE_OBJECT)
	{
		glEnable(GL_DEPTH_TEST);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDisable(GL_LIGHTING);
		renderBruteFlat(true);
	}
	else if (selectionMode == SELECTION_MODE_FACET) 
	{
		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDisable(GL_LIGHTING);
		renderBruteWithNames();
	}
	else
	{
		//regular rendering

		if (m_renderingParams->m_superimposeVertices)
		{
		}

		if(m_renderingParams->m_superimposeEdges)
		{
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_BLEND);
			// enable polygon offset
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(3.0f,1.0f);
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

			glLineWidth((float)m_renderingParams->m_edgeThickness);
			renderEdges();
		}

		// polygon mode (point, line or fill)
		//glPolygonMode(GL_FRONT_AND_BACK, m_renderingParams->m_polygonMode);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		// antialiasing
		if(m_renderingParams->m_antialiasing)
		{
			glEnable(GL_LINE_SMOOTH);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
			glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
			glLineWidth(1.5f);
		}
		else
		{
			glDisable(GL_LINE_SMOOTH);
			glDisable(GL_BLEND);
			glLineWidth(1.0f);
		}

		// lighting option
		if(m_renderingParams->m_lighting)
		{
			m_renderingParams->m_useNormals = true;
			glEnable(GL_LIGHTING);
		}
		else
		{
			m_renderingParams->m_useNormals = false;
			glDisable(GL_LIGHTING);
		}

		glDisable(GL_POLYGON_OFFSET_LINE);

		///////////////////////////////////////////
		//renderBuffersDebug();
		if (getSelectedFacetPainter()->isSmoothShading())
			renderBuffers();
		else
		{
			if (m_renderingParams->m_vboSupported)
				renderBuffersFlatShaded();
			else
				renderDisplayLists();
		}
		///////////////////////////////////////////
	}

	glPopMatrix();
}

void CageRenderer::fillEdgeBuffer()
{
	assert(m_edgeIndexesBuffer);

	Mesh::Edge_const_iterator hit = m_cage->edges_begin();
	Mesh::Edge_const_iterator hit_end = m_cage->edges_end();

	int dataIndex = 0;
	for(; hit != hit_end; hit++)
	{
		m_edgeIndexesBuffer[dataIndex	] = hit->vertex()->index();
		m_edgeIndexesBuffer[dataIndex+1	] = hit->opposite()->vertex()->index();
		dataIndex += 2;
	}
}

void CageRenderer::renderEdges()
{
	//if list not ready create it
	if (m_edgesDisplayListId == 0)
	{
		m_edgesDisplayListId = glGenLists(1);

		// Render model into the display list
		glNewList(m_edgesDisplayListId, GL_COMPILE);
		{
			renderEdgesBrute();
		}
		// Finished Creating the Display List
		glEndList();
	}

	//call the display list
	glCallList( m_edgesDisplayListId );
}

void CageRenderer::renderEdgesBrute()
{
	glBegin(GL_LINES);

	glColor4ub(
		m_renderingParams->m_edgeColor.red(),
		m_renderingParams->m_edgeColor.green(),
		m_renderingParams->m_edgeColor.blue(),
		255);

	const int numEdges = m_cage->size_of_halfedges() / 2;
	for (int i=0; i < numEdges; i++) {
		glVertex3f(
			m_vertexBuffer[m_edgeIndexesBuffer[i*2]*3],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2]*3+1],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2]*3+2]);

		glVertex3f(
			m_vertexBuffer[m_edgeIndexesBuffer[i*2+1]*3],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2+1]*3+1],
			m_vertexBuffer[m_edgeIndexesBuffer[i*2+1]*3+2]);
	}

	glEnd();
}

void CageRenderer::renderBrute()
{
	if (getSelectedFacetPainter()->isSmoothShading())
		renderBruteSmooth();
	else
		renderBruteFlat();
}

void CageRenderer::renderBruteSmooth()
{
	glBegin(GL_TRIANGLES);
	for (int i=0; i<m_sizeFacets; i++) {
		for (int j = i*3; j< i*3+3; j++) {
			glColor4ub(
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4],
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4+1],
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4+2],
				m_vertexColorBuffer[m_vertexIndexesBuffer[j]*4+3]);

			glVertex3f(
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+1],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+2]);

			glNormal3f(
				m_vertexNormalBuffer[m_vertexIndexesBuffer[j]*3],
				m_vertexNormalBuffer[m_vertexIndexesBuffer[j]*3+1],
				m_vertexNormalBuffer[m_vertexIndexesBuffer[j]*3+2]);
		}
	}
	glEnd();
}

void CageRenderer::renderBruteWithNames()
{
	int sizeOfFacets = m_sizeFacets;
	for (int i=0; i<sizeOfFacets; i++) {
		glPushName(i);
		glBegin(GL_TRIANGLES);
		for (int j = i*3; j< i*3+3; j++) {
			glVertex3f(
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+1],
				m_vertexBuffer[m_vertexIndexesBuffer[j]*3+2]);
			glNormal3f(
				m_facetNormalBuffer[j*3],
				m_facetNormalBuffer[j*3+1],
				m_facetNormalBuffer[j*3+2]);
		}
		glEnd();
		glPopName();
	}

}

void CageRenderer::renderBruteFlat(const bool selection)
{
	if (m_renderingParams->m_smoothNormals)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	glBegin(GL_TRIANGLES);

	if (selection) {
		for (int i=0; i<m_sizeFacets; i++) {
			for (int j = i*3; j< i*3+3; j++) {
				glVertex3f(
					m_facetVerticesBuffer[j*3],
					m_facetVerticesBuffer[j*3+1],
					m_facetVerticesBuffer[j*3+2]);
				glNormal3f(
					m_facetNormalBuffer[j*3],
					m_facetNormalBuffer[j*3+1],
					m_facetNormalBuffer[j*3+2]);
			}
		}
	} else {
		for (int i=0; i<m_sizeFacets; i++) {
			for (int j = i*3; j< i*3+3; j++) {
				glColor4ub(
					m_vertexColorBuffer[j*4],
					m_vertexColorBuffer[j*4+1],
					m_vertexColorBuffer[j*4+2],
					m_vertexColorBuffer[j*4+3]);
				glVertex3f(
					m_facetVerticesBuffer[j*3],
					m_facetVerticesBuffer[j*3+1],
					m_facetVerticesBuffer[j*3+2]);
				glNormal3f(
					m_facetNormalBuffer[j*3],
					m_facetNormalBuffer[j*3+1],
					m_facetNormalBuffer[j*3+2]);
			}
		}
	}

	glEnd();
}

void CageRenderer::renderDisplayLists()
{
	//if list not ready create it
	if (m_displayListId == 0)
	{
		m_displayListId = glGenLists(1);

		// Render model into the display list
		glNewList(m_displayListId, GL_COMPILE);
		{
			renderBrute();
		}
		// Finished Creating the Display List
		glEndList();
	}

	//call the display list
	glCallList( m_displayListId );
}

void CageRenderer::checkGlErrors()
{
	GLenum err;
	while ((err = glGetError()) != GL_NO_ERROR)
	{
		qDebug("GLError = %d", err);
	}
}

void CageRenderer::renderBuffersFlatShaded()
{
	//assert (m_facetNormalBuffer && m_facetVerticesBuffer && m_vertexColorBuffer);

	if (m_renderingParams->m_smoothNormals)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	if (!m_renderingParams->m_vboSupported) return;

	glEnableClientState( GL_VERTEX_ARRAY );				// Enable Vertex Arrays
	glEnableClientState( GL_NORMAL_ARRAY );				// Enable Normal Arrays
	glEnableClientState( GL_COLOR_ARRAY  );				// Enable Color Arrays

	glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_facetVerticesBufferPtr);
	glVertexPointer( 3, GL_FLOAT, 0, (char *) NULL );		// Set The Vertex Pointer To The Vertex Buffer

	glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_facetNormalBufferPtr);
	glNormalPointer( GL_FLOAT, 0, (char *) NULL );		// Set The normal Pointer To The normals Buffer

	glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexColorBufferPtr);
	glColorPointer(4, GL_UNSIGNED_BYTE, 0, (char *) NULL );		// Set The normal Pointer To The color Buffer

	GLint numOfElements = m_sizeFacets * 3;
	//glDrawElements( GL_TRIANGLES, numOfElements, GL_UNSIGNED_INT, (char *)NULL);
	glDrawArrays(GL_TRIANGLES, 0, numOfElements);

	glDisableClientState( GL_VERTEX_ARRAY	);				// Disable Vertex Arrays
	glDisableClientState( GL_NORMAL_ARRAY	);				// Disable Normal Arrays
	glDisableClientState( GL_COLOR_ARRAY	);				// Disable Color Arrays
}

void CageRenderer::renderBuffers()
{
	assert (m_vertexBuffer && m_vertexNormalBuffer && m_vertexColorBuffer);


	glEnableClientState( GL_VERTEX_ARRAY );				// Enable Vertex Arrays
	glEnableClientState( GL_NORMAL_ARRAY );				// Enable Normal Arrays
	glEnableClientState( GL_COLOR_ARRAY  );				// Enable Color Arrays

	if(m_renderingParams && m_renderingParams->m_vboSupported)	{
		glEnableClientState( GL_INDEX_ARRAY  );				// Enable Index Arrays

		glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexBufferPtr );
		glVertexPointer( 3, GL_FLOAT, 0, (char *) NULL );		// Set The Vertex Pointer To The Vertex Buffer

		glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexNormalBufferPtr);
		glNormalPointer( GL_FLOAT, 0, (char *) NULL );		// Set The normal Pointer To The normals Buffer

		glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_vertexColorBufferPtr);
		glColorPointer(4, GL_UNSIGNED_BYTE, 0, (char *) NULL );		// Set The normal Pointer To The color Buffer

		glBindBufferARB( GL_ELEMENT_ARRAY_BUFFER_ARB, m_vertexIndexesBufferPtr);
	} else {
		glVertexPointer( 3, GL_FLOAT, 0, m_vertexBuffer );	// Set The Vertex Pointer To Our Vertex Data
		glNormalPointer(GL_FLOAT, 0, m_vertexNormalBuffer );	// Set The Vertex Pointer To Our TexCoord Data
		glColorPointer(4, GL_UNSIGNED_BYTE, 0, m_vertexColorBuffer); //set color buffer
	}

	GLint numOfElements = m_sizeFacets * 3;

	if (m_renderingParams && m_renderingParams->m_vboSupported) {
		glDrawElements( GL_TRIANGLES, numOfElements, GL_UNSIGNED_INT, (char *)NULL);
		glDisableClientState( GL_INDEX_ARRAY);				// Disable Index Arrays
	} else {
		glDrawElements( GL_TRIANGLES, numOfElements, GL_UNSIGNED_INT, m_vertexIndexesBuffer);
	}

	glDisableClientState( GL_VERTEX_ARRAY	);				// Disable Vertex Arrays
	glDisableClientState( GL_NORMAL_ARRAY	);				// Disable Normal Arrays
	glDisableClientState( GL_COLOR_ARRAY	);				// Disable Color Arrays
}

void CageRenderer::fillVertexBuffer()
{
	assert(m_cage && m_vertexBuffer && m_vertexNormalBuffer);

	//I assume buffer is allocated already (when loading new mesh)

	int dataIndex = 0;

	//Go over vertices
	Mesh::Vertex_const_iterator vit = m_cage->vertices_begin();
	Mesh::Vertex_const_iterator vit_end = m_cage->vertices_end();

	//GLfloat inv_scale = 1 / m_scale;

	for (;vit != vit_end; vit++) {
		Point_3 p1 = vit->point();
		Vector_3 n1 = vit->normal();

		m_vertexBuffer[dataIndex] = p1[0]; 
		m_vertexBuffer[dataIndex+1] = p1[1];
		m_vertexBuffer[dataIndex+2] = p1[2];

		m_vertexNormalBuffer[dataIndex] = n1[0] /* * inv_scale*/;
		m_vertexNormalBuffer[dataIndex+1] = n1[1] /** inv_scale */;
		m_vertexNormalBuffer[dataIndex+2] = n1[2] /** inv_scale */;

		dataIndex += 3;
	}

	//Go over facets (fill indexes and facet normals)
	dataIndex = 0;
	Indexer *indr;
	indr = new CountIndexer(m_cage->size_of_facets());

	while (indr->hasNext())
	{
		int i = indr->next();
		
		Mesh::Facet_handle fit = m_cage->findFacet(i);
		Mesh::Halfedge_const_handle c = fit->halfedge();

		//first point of triangle
		Mesh::Vertex_const_handle v1 = c->vertex();
		Mesh::Vertex_const_handle v2 = c->next()->vertex();
		Mesh::Vertex_const_handle v3 = c->next()->next()->vertex();

		m_vertexIndexesBuffer[dataIndex] = v1->index();
		m_vertexIndexesBuffer[dataIndex+1] = v2->index();
		m_vertexIndexesBuffer[dataIndex+2] = v3->index();

		if (m_renderingParams->m_smoothNormals)
		{
			/*m_facetNormalBuffer[dataIndex*3	]	= m_vertexNormalBuffer[v1->index()*3];
			m_facetNormalBuffer[dataIndex*3	+1] = m_vertexNormalBuffer[v1->index()*3 + 1];
			m_facetNormalBuffer[dataIndex*3	+2] = m_vertexNormalBuffer[v1->index()*3 + 2 ];
			m_facetNormalBuffer[dataIndex*3	+3]	= m_vertexNormalBuffer[v2->index()*3];
			m_facetNormalBuffer[dataIndex*3	+4] = m_vertexNormalBuffer[v2->index()*3 + 1];
			m_facetNormalBuffer[dataIndex*3	+5] = m_vertexNormalBuffer[v2->index()*3 + 2];
			m_facetNormalBuffer[dataIndex*3	+6]	= m_vertexNormalBuffer[v3->index()*3];
			m_facetNormalBuffer[dataIndex*3	+7] = m_vertexNormalBuffer[v3->index()*3 + 1];
			m_facetNormalBuffer[dataIndex*3	+8] = m_vertexNormalBuffer[v3->index()*3 + 2];*/
			m_facetNormalBuffer[dataIndex*3	]	= v1->normal().x();
			m_facetNormalBuffer[dataIndex*3	+1] = v1->normal().y();
			m_facetNormalBuffer[dataIndex*3	+2] = v1->normal().z();
			m_facetNormalBuffer[dataIndex*3	+3]	= v2->normal().x();
			m_facetNormalBuffer[dataIndex*3	+4] = v2->normal().y();
			m_facetNormalBuffer[dataIndex*3	+5] = v2->normal().z();
			m_facetNormalBuffer[dataIndex*3	+6]	= v3->normal().x();
			m_facetNormalBuffer[dataIndex*3	+7] = v3->normal().y();
			m_facetNormalBuffer[dataIndex*3	+8] = v3->normal().z();
		} else {
			m_facetNormalBuffer[dataIndex*3	]	= fit->normal().x();
			m_facetNormalBuffer[dataIndex*3	+1] = fit->normal().y();
			m_facetNormalBuffer[dataIndex*3	+2] = fit->normal().z();
			m_facetNormalBuffer[dataIndex*3	+3]	= fit->normal().x();
			m_facetNormalBuffer[dataIndex*3	+4] = fit->normal().y();
			m_facetNormalBuffer[dataIndex*3	+5] = fit->normal().z();
			m_facetNormalBuffer[dataIndex*3	+6]	= fit->normal().x();
			m_facetNormalBuffer[dataIndex*3	+7] = fit->normal().y();
			m_facetNormalBuffer[dataIndex*3	+8] = fit->normal().z();
		}

		m_facetVerticesBuffer[dataIndex*3	]	= v1->point().x();
		m_facetVerticesBuffer[dataIndex*3	+1] = v1->point().y();
		m_facetVerticesBuffer[dataIndex*3	+2] = v1->point().z();
		m_facetVerticesBuffer[dataIndex*3	+3]	= v2->point().x();
		m_facetVerticesBuffer[dataIndex*3	+4] = v2->point().y();
		m_facetVerticesBuffer[dataIndex*3	+5] = v2->point().z();
		m_facetVerticesBuffer[dataIndex*3	+6]	= v3->point().x();
		m_facetVerticesBuffer[dataIndex*3	+7] = v3->point().y();
		m_facetVerticesBuffer[dataIndex*3	+8] = v3->point().z();
		//}

		dataIndex += 3;
	}

	delete indr;
}

void CageRenderer::bindVertexColorBuffer()
{
	if (!m_renderingParams || !m_renderingParams->m_vboSupported) return;

	//if not smooth shading we won't use VBO's
	//if (!getSelectedFacetPainter()->isSmoothShading()) return;

	m_vertexColorBufferVBO = true;

	int bufferSize = getSelectedFacetPainter()->requiredBufferSize(m_cage);

	checkGlErrors();

	//bind color buffer to graphics card
	glGenBuffersARB(1, &m_vertexColorBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexColorBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLubyte) * bufferSize, m_vertexColorBuffer,
		GL_STATIC_DRAW_ARB);


	checkGlErrors();
}

void CageRenderer::unbindVertexColorBuffer()
{
	if (!m_vertexColorBufferVBO) return;

	glDeleteBuffersARB(1, &m_vertexColorBufferPtr);

	m_vertexColorBufferVBO = false;
}

void CageRenderer::bindVertexBuffer()
{
	if (!m_renderingParams || !m_renderingParams->m_vboSupported) return;

	checkGlErrors();

	//bind buffer to graphics card
	glGenBuffersARB(1, &m_vertexBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeVertices * 3, m_vertexBuffer,
		GL_STATIC_DRAW_ARB);

	//bind normals buffer to graphics card
	glGenBuffersARB(1, &m_vertexNormalBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_vertexNormalBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeVertices * 3, m_vertexNormalBuffer,
		GL_STATIC_DRAW_ARB);


	//bind indexes buffer to graphics card
	glGenBuffersARB(1, &m_vertexIndexesBufferPtr);
	glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, m_vertexIndexesBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, sizeof(GLuint) * m_sizeFacets * 3, m_vertexIndexesBuffer,
		GL_STATIC_DRAW_ARB);

	//////////////////////////////////////////////////////////////////////////
	checkGlErrors();

	//bind second buffer of vertices and normals (normals equal)
	glGenBuffersARB(1, &m_facetVerticesBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_facetVerticesBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeFacets * 9, m_facetVerticesBuffer,
		GL_STATIC_DRAW_ARB);

	//bind normals buffer to graphics card
	glGenBuffersARB(1, &m_facetNormalBufferPtr);
	glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_facetNormalBufferPtr);
	//store data into buffer object, which is accelerated gc memory
	//when available
	glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(GLfloat) * m_sizeFacets * 9, m_facetNormalBuffer,
		GL_STATIC_DRAW_ARB);

	checkGlErrors();

}

void CageRenderer::unbindVertexBuffer()
{
	glDeleteBuffersARB(1, &m_vertexBufferPtr);
	glDeleteBuffersARB(1, &m_vertexNormalBufferPtr);
	glDeleteBuffersARB(1, &m_vertexIndexesBufferPtr);

	glDeleteBuffersARB(1, &m_facetVerticesBufferPtr);
	glDeleteBuffersARB(1, &m_facetNormalBufferPtr);
}

void CageRenderer::fillVertexColorBuffer()
{
	assert(m_cage && m_vertexColorBuffer && m_selectedFacetPainter);

	FacetPainter* fp = getSelectedFacetPainter();

	assert(fp);

	if (m_vertexColorBuffer && m_vertexColorBufferSize != fp->requiredBufferSize(m_cage))
	{
		//delete the old
		delete[] m_vertexColorBuffer;
		m_vertexColorBuffer = NULL;
		//in with the new
		m_vertexColorBuffer = new GLubyte[fp->requiredBufferSize(m_cage)];
	}
	else if (!m_vertexColorBuffer)
	{
		//just create a new one
		m_vertexColorBuffer = new GLubyte[fp->requiredBufferSize(m_cage)];
	}

	fp->fillBuffer(m_renderingParams, m_cage, NULL, m_dmet, NULL, m_cageSelection, m_vertexColorBuffer);
}


void CageRenderer::OnObjectChanged()
{
	OnObjectChanged(false);
}

void CageRenderer::OnObjectChanged(bool geometryChanged)
{
	if (geometryChanged)
	{
		setCage(m_cage, NULL, m_cageSelection);
		return;
	}

	if (m_renderingParams && m_renderingParams->m_vboSupported && m_vertexColorBufferVBO)
	{
		unbindVertexColorBuffer();
	}

	fillVertexColorBuffer();

	if (m_renderingParams && m_renderingParams->m_vboSupported)
	{
		bindVertexColorBuffer();
	}

	// display list is no longer relevant, delete it
	if (m_displayListId != 0)
	{
		glDeleteLists(m_displayListId, 1);
		m_displayListId = 0;
	}

	if (m_edgesDisplayListId != 0)
	{
		glDeleteLists(m_edgesDisplayListId, 1);
		m_edgesDisplayListId = 0;
	}
}


void CageRenderer::OnFacetPainterChanged(FacetPainter* facetPainter)
{
	if (m_selectedFacetPainter == facetPainter) return;

	m_selectedFacetPainter = facetPainter;

	if (m_cage)
		OnObjectChanged();
}
