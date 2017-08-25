#ifndef __FEATURE_PAINTERS_H_
#define __FEATURE_PAINTERS_H_

#include "MeshRenderer.h"

class DiscreteColors
{
public:
	DiscreteColors(const QColor &defaultc) :m_default(defaultc) {}

	QColor getColor(int n)
	{
		if ((n >= 0) && (n < size))
			return QColor(table[n]);
		else
			return m_default;
	}
	QColor getColorWrap(int n)
	{
		if (n < 0)
			return m_default;
		return QColor(table[n % size]);
	}

	static QColor getDefColor(int n)
	{
		return defaultCol.getColor(n);
	}
	static QColor getDefColorWrap(int n)
	{
		return defaultCol.getColorWrap(n);
	}

private:
	QColor m_default;
	static unsigned int *table;
	static int size;

	static DiscreteColors defaultCol;

	friend class QColorOrder;
};

extern unsigned int CORRESPONDENCE_COLORS[];
extern int CORRESPONDENCE_COLORS_SIZE;

#include "RayIntersect.h"
class RayCellsPainter : public DebugPainter
{
private:
	RayIntersect m_rayIntersect;
	
public:
	RayCellsPainter() {};

	virtual QString name() { return "Ray Cells"; }

	virtual void paint(Mesh* mesh, MeshSelection* meshSelection = NULL, MeshSurfaceGraph* msg = NULL);
};


//////////////////////////////////////////////////////////////////////////

/**
 *	Renders shape diameter function value as calculated on vertices of the mesh
 */
class SdfVerticesPainter : public FacetPainter 
{
	virtual QString name() {return "SDF Vertices";}

	virtual bool isSmoothShading() { return true; }

	virtual unsigned int requiredBufferSize(Mesh* mesh);

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer);

	virtual bool paintVertex(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Vertex_const_handle& pVertex,
		GLubyte* color);

	virtual bool paintVertexInFacet(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Facet_const_handle& pFacet,
		const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
		GLubyte* color);
};

//////////////////////////////////////////////////////////////////////////

/**
*	Renders shape diameter function value as calculated on vertices of the mesh
*/
class SdfFacetsPainter : public FacetPainter {
private:
	bool m_normalize;
	QString m_name;
	
public:
	SdfFacetsPainter(const QString& name = "SDF Facets", const bool normalize = false) 
		: m_name(name), m_normalize(normalize)
	{		
	}

	virtual QString name() {return m_name;}

	virtual bool isSmoothShading() {
		return false;
	}

	double logit(const double& value, const int& range);

	virtual unsigned int requiredBufferSize(Mesh* mesh);

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer);

	virtual bool paintVertex(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Vertex_const_handle& pVertex,
		GLubyte* color);

	virtual bool paintVertexInFacet(
		const RenderingParams* renderingParams,
		const Enriched_Mesh::Facet_const_handle& pFacet,
		const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
		GLubyte* color);
};

//////////////////////////////////////////////////////////////////////////

class PartitionPainter : public SdfFacetsPainter {
	virtual QString name() {return "Partition";}

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer);
};

//////////////////////////////////////////////////////////////////////////

class SdfDifferencesFacetsPainter : public SdfFacetsPainter {
	virtual QString name() {return "SDF Diff";}

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer);
};

//////////////////////////////////////////////////////////////////////////

class DeformationDegreeFacetsPainter : public SdfFacetsPainter {
	virtual QString name() { return "DDegree"; }

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer);
};

class DeformationDegreeDifferencesFacetsPainter : public SdfFacetsPainter {
	virtual QString name() { return "DDegree Diff"; }

	virtual bool fillBuffer(
		const RenderingParams* renderingParams,
		Mesh* mesh,
		MeshSurfaceGraph* msg,
		MeshSelection* meshSelection,
		GLubyte* colorBuffer);
};

//////////////////////////////////////////////////////////////////////////

/**
*	Renders shape diameter function value max gradient as calculated on facets of the mesh
*/
// class SdfGradientPainter : public FacetPainter {
// private:
// 	bool m_normalize;
// 	QString m_name;
// 
// public:
// 	SdfGradientPainter(const QString& name = "SDF Gradients", const bool normalize = false) 
// 		: m_name(name), m_normalize(normalize)
// 	{		
// 	}
// 
// 	virtual QString name() {return m_name;}
// 
// 	virtual bool isSmoothShading() {
// 		return false;
// 	}
// 
// 	double logit(const double& value, const int& range);
// 
// 	virtual unsigned int requiredBufferSize(Mesh* mesh);
// 
// 	virtual bool fillBuffer(
// 		const RenderingParams* renderingParams,
// 		Mesh* mesh,
// 		MeshSurfaceGraph* msg,
// 		MeshSelection* meshSelection,
// 		GLubyte* colorBuffer);
// 
// 	virtual bool paintVertex(
// 		const RenderingParams* renderingParams,
// 		const Enriched_Mesh::Vertex_const_handle& pVertex,
// 		GLubyte* color);
// 
// 	virtual bool paintVertexInFacet(
// 		const RenderingParams* renderingParams,
// 		const Enriched_Mesh::Facet_const_handle& pFacet,
// 		const Enriched_Mesh::Halfedge_const_handle& pHalfedge,
// 		GLubyte* color);
// };
/*
class CurvatureVertexPainter : public FacetPainter {
private:
	QString m_name;
	int m_type;

	double getval(const Mesh::Vertex_const_handle& v);
public:
	CurvatureVertexPainter(const QString& iname, const int itype) :
	m_name(iname), m_type(itype) 
	{
	}

	virtual QString name() {return m_name;}

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
		GLubyte* colorBuffer);
};*/

//////////////////////////////////////////////////////////////////////////

class SdfFacetDebugPainter : public DebugPainter {

private:

	RayIntersect m_rayIntersect;

	TNT::Array3D<RayIntersectCell*> m_grid;
	triangleVector_t m_triangleVector;
	double m_xmin, m_xmax, m_ymin, m_ymax, m_zmin, m_zmax, m_xspread, m_yspread, m_zspread;

	int m_cones;
	float m_coneSeperationDegree;
	double m_coneSeparationRadian;
	int m_raysInCone;
	int m_gridSize;

public:

	SdfFacetDebugPainter();

	virtual QString name() { return "SDF Facet Rays"; }

	// virtual QString command(const int key);

	virtual void paint(Mesh* mesh, MeshSelection* meshSelection = NULL, MeshSurfaceGraph* msg = NULL);

	void setCones(int cones) { m_cones = cones; }
	void setDegree(float degree) { m_coneSeperationDegree = degree; }
	void setRadian(double radian) { m_coneSeparationRadian = radian; }
	void setRays(int rays) { m_raysInCone = rays; }
	void setGridSize(int gridSize) { m_gridSize = gridSize; }

};

//////////////////////////////////////////////////////////////////////////

class SdfVertexDebugPainter : public DebugPainter {
private:
	int m_cones;
	float m_coneSeperationDegree;
	double m_coneSeparationRadian;
	int m_raysInCone;
	int m_gridSize;

public:
	SdfVertexDebugPainter();

	virtual QString name() { return "SDF Vertex Rays"; }

	virtual void paint(Mesh* mesh, MeshSelection* meshSelection = NULL, MeshSurfaceGraph* msg = NULL);

	void setCones(int cones) { m_cones = cones; }
	void setDegree(float degree) { m_coneSeperationDegree = degree; }
	void setRadian(double radian) { m_coneSeparationRadian = radian; }
	void setRays(int rays) { m_raysInCone = rays; }
	void setGridSize(int gridSize) { m_gridSize = gridSize; }
};

//////////////////////////////////////////////////////////////////////////

class DihedralDebugPainter : public DebugPainter {
private:
	QString m_name;
	int m_type;

	double calcWeight(Mesh::Facet_const_handle& f1, Mesh::Facet_const_handle& f2, Mesh::Halfedge_const_handle& h);
public:

	DihedralDebugPainter(const QString& name, const int type) : m_name(name), m_type(type) {}

	virtual QString name() {return m_name;}

	virtual void paint(Mesh* mesh, MeshSelection* MeshSelection = NULL, MeshSurfaceGraph* msg = NULL);
};


//////////////////////////////////////////////////////////////////////////

/*
class CagePainter : public DebugPainter {
private:
	
public:
	CagePainter();

	virtual QString name() {return "Cage Mesh"; }

	virtual void paint(Mesh* mesh, MeshSelection* MeshSelection = NULL, MeshSurfaceGraph* msg = NULL);
};
*/

#endif