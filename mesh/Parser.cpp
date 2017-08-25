#include "stdafx.h"

#include "Parser.h"


MeshBuilder::EAddResult MeshBuilder::addTriangle(int a, int b, int c)
{
	int vertices[3] = { a, b, c };
	return addFacet(Face(vertices, 3));
}

MeshBuilder::EAddResult MeshBuilder::addFacet(const Face &f)
{
	if (!m_builder.test_facet(&(f.v[0]), &(f.v[f.size])))
		return ADD_ILLEGAL;

	Mesh::Halfedge_handle he = m_builder.add_facet(&(f.v[0]), &(f.v[f.size]));
	if (he != NULL && he->facet() != NULL)
	{
		he->facet()->index(m_faceIndex++);
	}
	else
	{
		printf("error: %d= ", m_faceIndex);
		for(int i = 0; i < f.size; ++i)
			printf("%d ", f.v[i]);
		printf("\n");
		return ADD_ERROR;
	}
	return ADD_OK;
}

bool SimpleMeshBuilder::processFaces(FaceList &faces)
{
	for(int i = 0; i < faces.size(); ++i)
	{
		faces.resetRead();
		EAddResult ret = addFacet(*faces.next());
		if (ret == ADD_ERROR)
			return false;
	}
	return true;
}

bool RepeatingMeshBuilder::processFaces(FaceList &faces)
{
	int rep = 0;
	bool done = false;
	int screwed = 0;
	while ((rep < 30) && !done)
	{
		bool err = false;
		if (rep != 0)
		{
			Mesh::HDS hds = m_mesh->get_hds();
			hds.edges_clear();
			hds.faces_clear();
			ExposedBuilder* exb = (ExposedBuilder*)&m_builder;
			exb->emptyEdgeMap();
		}
		int i;
		faces.resetRead();
		for(i = 0; i < faces.size(); ++i)
		{
			Face &f = *faces.next();
			if (!f.isValid())
				continue;
			EAddResult ret = addFacet(f);
			if (ret == ADD_ILLEGAL)
			{
				++screwed;
				f.setValid(false);
			}
			else if (ret == ADD_ERROR)
			{
				printf("error at vertex %d, aborting\n", i);
				++screwed;
				f.setValid(false);
				break;
			}
		}
		done = (i == faces.size());
		++rep;
	}

//	printf("filtered=%d ", filtered);

	return done;
}



struct TriAdd
{
	TriAdd(int l = -1, bool* a = NULL) :last(l), toadd(a) {}
	int last;
	bool *toadd;
};

typedef QHash<int, TriAdd> TInCirc;
typedef QHash<int, TInCirc> TriMap;


bool LegalMeshBuilder::processFaces(FaceList &faces)
{
/*	TriMap vtxToTris;
	QVector<bool> tris;
	tris.reserve(m_numTriangles * 2);

	for (int i = 0; i < m_numTriangles; ++i)
	{
		if (vsize == 3)
		{
			tris.push_back(false);
			bool *curt = &tris.back();
			vtxToTris[vertices[0]][vertices[1]] = TriAdd(vertices[2], curt);
			vtxToTris[vertices[1]][vertices[2]] = TriAdd(vertices[0], curt);
			vtxToTris[vertices[2]][vertices[0]] = TriAdd(vertices[1], curt);
		}

	}

	//if (filtered != 0)
	//	printf("filtered=%d\n", filtered);

	for(int v = 0; v < vtxToTris.size(); ++v)
	{
		int w, u = -1;
		TriAdd firstTr;
		TInCirc &localCirc = vtxToTris[v];
		TInCirc::iterator cit = localCirc.begin();
		do {
			w = u;
			TriAdd &tr = cit.value();
			u = tr.last;

			*(tr.toadd) = true;
			localCirc.remove(cit.key());
			cit = localCirc.find(u);
		} while (cit != localCirc.end());


		if (localCirc.size() != 0)
		{
			printf("found double-vertex %d with: ", v);
			foreach(const TriAdd& ti, localCirc)
			{
				*ti.added = true;
				printf("%d ", ti.last);
			}
			printf("\n");
		}
	}
	printf("%d-%d\n", m_index, m_numTriangles);
	return true;


	//return (index == m_numTriangles); // got to the end
*/
	return true;
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


bool Ply2Reader::read(const char* pFilename)
{
	string fname(pFilename);
	string extension = fname.substr(fname.find_last_of("."));

	ifstream file(pFilename);
	if (!file.is_open())
		return false;

	bool readFacetSize = true;
	bool readEdgesNum = false;

	if (extension == ".ply2") 
	{
		readFacetSize = true;
	}
	else if (extension == ".simple")
	{
		readFacetSize = false;
	}
	else if (extension == ".off")
	{
		string offhead;
		file >> offhead;
		readFacetSize = true;
		readEdgesNum = true;
	}

	return read(file, readFacetSize, readEdgesNum);
}


bool Ply2Reader::read(ifstream& file, bool readFacetSize, bool readEdgesNum)
{
	file >> m_numVertices;
	file >> m_numTriangles;
	int edgesNum = 6;
	if (readEdgesNum)
		file >> edgesNum;

	bool success = false;
	int times = 0;

	m_build->startSurface(m_numVertices, m_numTriangles, edgesNum);
	m_build->setVtxStart(0);

	// read the vertices
	float x,y,z;
	for (int i = 0; i < m_numVertices; ++i)
	{
		file >> x;
		file >> y;
		file >> z;
		m_build->addVtx(x, y, z);
	}

	// read faces

	success = read_facets(file, readFacetSize);

	m_build->endSurface();

	return success;
}


// read facets and uv coordinates per halfedge
bool Ply2Reader::read_facets(ifstream& file, bool readFacetSize)
{
	Face f;
	FaceList faces;

	int v, i;
	int fieldCount = 3;

	faces.reserve(m_numTriangles);

	for (i = 0; i < m_numTriangles; ++i)
	{
		f.reset();
		for (int j = 0; j < fieldCount; ++j)
		{
			file >> v;
			if ((j == 0) && readFacetSize)
				fieldCount = v + 1;
			else
				f.add(v);
		}
		faces.add(f);
	}

	return m_build->processFaces(faces);
}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


bool ObjReader::read(const char* filename)
{
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly))
		return false;
	QTextStream in(&file);

	m_build->startSurface(1, 1, 3);
	m_build->setVtxStart(0);

	FaceList faces;

	int filtered = 0;
	QString line;
	do {
		line = in.readLine();
		if (line.length() == 0)
			continue;
		QStringList flds = line.split(' ', QString::SkipEmptyParts);
		if (flds[0] == "v")
		{
			m_build->addVtx(flds[1].toDouble(), flds[2].toDouble(), flds[3].toDouble());
		}
		else if (flds[0] == "f")
		{
			int numpnts = flds.size() - 1;
			if (numpnts > 10)
			{
				printf("Error: got a face with %d vertices\n", numpnts);
				return false;
			}
			Face f;
			for(int i = 1; i < flds.size(); ++i)
			{ // cut away texture and normal information
				int si = flds[i].indexOf('/');
				if (si != -1)
					flds[i].resize(si);
				f.add(flds[i].toInt() - 1);
			}

			faces.add(f);
		}

	} while (!line.isNull());

	m_build->processFaces(faces);

	return true;
}