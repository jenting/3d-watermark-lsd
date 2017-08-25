// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef __STDAFX_H_INCLUDED__
#define __STDAFX_H_INCLUDED__

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
// C RunTime Header Files

#pragma warning (disable: 4503) // decorated name length exceeded, name was truncated
#pragma warning (disable: 4018) // signed/unsigned mismatch
#pragma warning (disable: 4244) // conversion loss of data
#pragma warning (disable: 4305) // truncation from to
#pragma warning (disable: 4099) // type name first seen using 'class' now seen using 'struct' (in CGAL)
#pragma warning (disable: 4996) // xxx was declared deprecated

#include <tchar.h>
#include <iostream>
#include <fstream>

#define DEG2RAD 0.01745329f
#define RAD2DEG 57.2957795f

typedef unsigned char byte;


// TODO: reference additional headers your program requires here
#include "GL/glew.h"
#include "gl/glext.h"

#include <vector>
#include <list>
#include <hash_map>
using namespace std;

//CGAL Stuff
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Kd_tree.h>
//#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/aff_transformation_tags.h>

#include <CGAL/squared_distance_3.h>

#include <QtGlobal>
#include <QObject>
#include <QString>
#include <QMessageBox>

typedef QVector<int> TIndexList;

#include "mesh/enriched_polyhedron.h"
#include "mesh/CGALTypes.h"
#include "mesh/Mesh.h"

#include <qsettings.h>
#include <qstringlist.h>

#include "GlobalLog.h"

//#include "mesh/mesh_utils.h"

typedef void TStubData;
typedef void MeshSurfaceGraph;


class MainWindow;
extern MainWindow *g_sdfmain;

#include <QMainWindow.h>
extern QMainWindow *g_main;

typedef QMap<int,int> TPartTable;

//#include <Windows.h>
//#include <gl/GL.h>

enum DrawMethod
{
	DRAW_PART = 0x1,
	DRAW_MODE_COLOR_PART = 0x2,
	DRAW_PART_MIX_COL = 0x5
};


class ProgSettings : public QSettings
{
public:
	/*ProgSettings() : QSettings() {
		setPath("liors.net", "sdf", QSettings::User);
	}*/

	ProgSettings() : QSettings() {
		setPath("hsiao", "MeshWatermark", QSettings::User);
	}

	ProgSettings(const QString& domain, const QString& prog) : QSettings() {
		setPath(domain, prog);
	}

};


class AppManager;
extern AppManager *g_appManager;

#endif
