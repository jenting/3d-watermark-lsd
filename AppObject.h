#ifndef __APPOBJECT_H_INCLUDED__
#define __APPOBJECT_H_INCLUDED__

#undef min
#undef max

class WorldObject;

/**
 *	Represents an object in the application, a mesh and its relevant status's and properties
 */
class AppObject 
{
public:
	AppObject() :  mesh(NULL),  cage(NULL), meshSelection(NULL), 
		 hasSdfFacets(false), refCount(new int(1)), m_myWo(NULL)
	{
		//fileSeparator = "\\";
	}

	void unloadMesh();
		
	~AppObject() 
	{
		--(*refCount);
		if (*refCount == 0)
		{
			delete mesh;
			delete cage;
			delete meshSelection;
			delete refCount;
		}
	}

	void setWorldObject(WorldObject* wo) { m_myWo = wo; }
	
	int numberOfFaces() { return mesh->size_of_facets(); }
	int numberOfVertices() { return mesh->size_of_vertices(); }

public:
	/** unique identifier */
	QString id;
	/** descriptor */
	QString description;
	/** file dir */
	QString fileDir;
	/** file name */
	QString inFileName;
	/** fule sub name */
	QString fileName;
	/** file suffix (default : .obj) */
	QString fileSuffix;
	/** file seperator, in Windows =  "\\" */
	QString fileSeparator;
	/** directory to storage sdf information */
	QString sdfDir;

	/** pointer to actual mesh object */
	Mesh* mesh;

	/** pointer to actual cage mesh object */
	Mesh* cage;

	//MyMesh *omesh;

	/** pointer to selection object on mesh */
	//MeshSelection* selection;
	MeshSelection* meshSelection;
	
	/** pointer to selection cage */
	MeshSelection* cageSelection;

	bool hasSdfFacets;

	WorldObject* m_myWo; 

	int *refCount;

	int maxLevel;
};


#endif // __APPOBJECT_H_INCLUDED__