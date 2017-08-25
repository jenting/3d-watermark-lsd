#pragma once

#include <QHash>

#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "Enriched_polyhedron.h"

using namespace std;

struct Face
{
	Face(int *from = NULL, int sz = 0)  : size(sz)
	{
		Q_ASSERT(size < 10);
		for (int i = 0; i < sz; ++i)
			v[i] = from[i];
	}
	
	void add(int i) 
	{ 
		Q_ASSERT(size < 10);
		v[size++] = i; 
	}
	void reset() { size = 0; }
	bool isValid() 
	{ 
		return v[0] >= 0; 
	}
	void setValid(bool b) 
	{ 
		v[0] = b?(qAbs(v[0])):(-qAbs(v[0])); 
	}

	int size;
	int v[10];
	
};

class FaceList
{
public:
	FaceList() :m_readOff(0), m_size(0) {}
	void reserve(quint32 res)
	{
		m_data.reserve(res);
	}
	void add(const Face& f)
	{
		m_data.append(f.size);
		for (int i = 0; i < f.size; ++i)
			m_data.append(f.v[i]);
		++m_size;
	}

	void resetRead() { m_readOff = 0; }

	Face* next()
	{
		Face* p = (Face*)(&m_data[m_readOff]);
		m_readOff += p->size + 1;
		return p;
	}

	bool hasNext()
	{
		return m_readOff < m_data.size();
	}

	quint32 size() { return m_size; }

private:
	QVector<int> m_data;
	int m_readOff;
	int m_size;
};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//typedef QVector<Face> TFaceList;
class MeshBuilder
{
public:
	typedef Mesh::HalfedgeDS::Vertex::Point Point;
	typedef CGAL::Polyhedron_incremental_builder_3<Mesh::HalfedgeDS> Builder;

	MeshBuilder(Mesh *pMesh) : m_vtxIndex(0), m_faceIndex(0), 
		m_builder(pMesh->get_hds(), true), m_mesh(pMesh)
	{
		pMesh->get_hds().clear();
	}

	virtual ~MeshBuilder() {}

	virtual void startSurface(int v, int f, int e)
	{
		initV = v; initF = f, initE = e;
		m_builder.begin_surface(v, f, e);
	}
	virtual void endSurface()
	{
		m_builder.end_surface();
	}
	virtual void setVtxStart(int startIndex) 
	{ 
		m_vtxIndex = startIndex; 
	}
	virtual void addVtx(double x, double y, double z)
	{
		Mesh::Vertex_handle v = m_builder.add_vertex(Point(x,y,z));
		v->index() = m_vtxIndex++;
	}
	virtual bool processFaces(FaceList &faces) = 0;

	class ExposedBuilder : public Builder
	{
	public:
		void emptyEdgeMap()
		{
			for(int i = 0; i < hds.size_of_vertices(); ++i)
				set_vertex_to_edge_map(i, Halfedge_handle());
			m_error = false;
		}
	};


protected:
	enum EAddResult
	{
		ADD_OK, ADD_ILLEGAL, ADD_ERROR
	};

	EAddResult addTriangle(int a, int b, int c);
	EAddResult addFacet(const Face &f);

	int m_vtxIndex, m_faceIndex;
	Builder m_builder;
	Mesh *m_mesh;
	int initV, initF, initE;
};

class SimpleMeshBuilder : public MeshBuilder
{
public:
	SimpleMeshBuilder(Mesh *pMesh) : MeshBuilder(pMesh) {}
	virtual ~SimpleMeshBuilder() {}
	virtual bool processFaces(FaceList &faces);
};

class RepeatingMeshBuilder : public MeshBuilder
{
public:
	RepeatingMeshBuilder(Mesh *pMesh) : MeshBuilder(pMesh) {}
	virtual ~RepeatingMeshBuilder() {}
	virtual bool processFaces(FaceList &faces);
};


class LegalMeshBuilder : public MeshBuilder
{
public:
	virtual ~LegalMeshBuilder() {}
	virtual bool processFaces(FaceList &faces);
};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


class FileReader
{
public:
	FileReader(MeshBuilder* bld) : m_build(bld) {}

	virtual bool read(const char* filename) = 0;

protected:
	MeshBuilder* m_build;
};


class Ply2Reader : public FileReader
{
public:
	Ply2Reader(MeshBuilder *bld) : FileReader(bld) {}
	virtual ~Ply2Reader() {}
	virtual bool read(const char* filename);

private:
	bool read(ifstream& file, bool readFacetSize, bool readEdgesNum);
	bool read_facets(ifstream& file, bool readFacetSize);

	int m_numVertices;
	int m_numTriangles;

};

class ObjReader : public FileReader
{
public:
	ObjReader(MeshBuilder *bld) : FileReader(bld) {}
	virtual ~ObjReader() {}
	virtual bool read(const char* filename);
};
