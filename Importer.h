#ifndef __IMPORTER_H_
#define __IMPORTER_H_

#include "stdafx.h"

#include <QFileInfo>

#include "AppObject.h"
#include "AppManager.h"

class FileParser {

public:

	static bool read_OffPly2(const char* filename, Mesh* mesh, int err = 1.0);
	static bool read_Obj(const char* filename, Mesh* mesh, int err = 1.0);

};


class ImporterMesh {
public:
	/** load a mesh from a file */
	static bool loadMesh(const QString& fileName, AppObject* appObject, double err = 1.0);
	static bool loadMesh(const QString& fileName, Mesh* mesh, double err = 1.0);
};

class ImporterMeshFeatures {
public:
	static bool loadVertexNormals(const QString& featuresFileName, double** vertexNormals);

	static bool loadSDFVertices(const QString& featuresFileName, const bool normalize, AppObject* appObject);
	static bool loadSDFFacets(const QString& featuresFileName, const bool normalize, AppObject* appObject);

	/** load an sdf information file and store the diff in values between current and new file in the volume field */
	static bool diffSDFFacets(const QString& featuresFileName, AppObject* appObject);


	static bool loadSDFText(const QString& importFileName, const bool onVertices, AppObject* appObject);
	//static bool loadSDFDiffText(const QString& featuresFileName, const bool onVertices, AppObject* appObject);

	static bool loadNSDFText(const QString& importFileName, const bool onVertices, AppObject* appObject);
	//static bool saveNSDFDiffText(const QString& featuresFileName, const bool onVertices, AppObject* appObject);
};

#endif