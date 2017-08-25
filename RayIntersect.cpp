#include "stdafx.h"

#include "RayIntersect.h"

#include <gl/GL.h>
#include <gl/GLU.h>

#include <qfile>
#include <q3textstream>
#include <qdatetime>

//#include "CylLine.h"

void RayIntersectCell::initGeometry(const Point_3& p1, const double xsize, const double ysize, const double zsize)
{
	m_points[0] = p1;
	m_points[1] = Point_3(p1.x(), p1.y(), p1.z() + zsize);
	m_points[2] = Point_3(p1.x(), p1.y() + ysize, p1.z());
	m_points[3] = Point_3(p1.x(), p1.y() + ysize, p1.z() + zsize);
	m_points[4] = Point_3(p1.x() + xsize, p1.y(), p1.z());
	m_points[5] = Point_3(p1.x() + xsize, p1.y(), p1.z() + zsize);
	m_points[6] = Point_3(p1.x() + xsize, p1.y() + ysize, p1.z());
	m_points[7] = Point_3(p1.x() + xsize, p1.y() + ysize, p1.z() + zsize);

	m_planes[0] = Plane_3(m_points[7],m_points[3],m_points[1]);
	m_planes[1] = Plane_3(m_points[1],m_points[3],m_points[2]);
	m_planes[2] = Plane_3(m_points[5],m_points[1],m_points[0]);
	m_planes[3] = Plane_3(m_points[7],m_points[5],m_points[4]);
	m_planes[4] = Plane_3(m_points[4],m_points[0],m_points[2]);
	m_planes[5] = Plane_3(m_points[6],m_points[2],m_points[3]);
}

int RayIntersectCell::intersects(const triangleVector_t& triangleVector, const Ray_3& ray, const bool debugMode)
{
	triangleIndexList_t::const_iterator end_it = m_triangleIndexList.end();
	for (triangleIndexList_t::const_iterator it = m_triangleIndexList.begin(); it != end_it; it++) 
		if (CGAL::do_intersect(*(triangleVector[*it]), ray))
			return *it;

	return -1;
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

RayIntersect::~RayIntersect()
{
	/*
	//delete m_grid
	for (int x = 0; x < grid.dim1(); x++)
		for (int y = 0; y < grid.dim2(); y++)
			for (int z = 0; z < grid.dim3(); z++)
				delete grid[x][y][z];

	for (int i = 0; i < m_triangleVector.size(); i++)
		delete triangleVector[i];

	triangleVector.clear();

	initialized = false;
	*/
}

void RayIntersect::Init(const Mesh& mesh, const int gridSize)
{
	//step 0 - make sure all data structures are empty
	if (m_initialized)
	{
		for (int x=0;x<m_grid.dim1();x++)
			for (int y=0;y<m_grid.dim2();y++)
				for (int z=0;z<m_grid.dim3();z++)
				{
					delete m_grid[x][y][z];
				}

		for (int i=0;i<m_triangleVector.size();i++)
			delete m_triangleVector[i];

		m_triangleVector.clear();
	}

	//step 1 - run over all triangles and add them to triangle vector
	m_triangleVector.reserve(mesh.size_of_facets());

	for (Mesh::Facet_const_iterator f = mesh.facets_begin(); f != mesh.facets_end(); f++)
	{
		CGAL_assertion(f->size() == 3);

		Mesh::Halfedge_const_handle h1 = f->halfedge();
		Mesh::Halfedge_const_handle h2 = h1->next();
		Mesh::Halfedge_const_handle h3 = h2->next();
		Triangle_3_With_Normal* triangle = new Triangle_3_With_Normal(
											h1->vertex()->point(), h2->vertex()->point(),
											h3->vertex()->point(), f->normal());

		m_triangleVector.push_back(triangle);		               
	}
		
	m_xmin = mesh.xmin(); m_xmax = mesh.xmax();
	m_ymin = mesh.ymin(); m_ymax = mesh.ymax();
	m_zmin = mesh.zmin(); m_zmax = mesh.zmax();
	m_xspread = 1 / (m_xmax - m_xmin + 1e-5);
	m_yspread = 1 / (m_ymax - m_ymin + 1e-5);
	m_zspread = 1 / (m_zmax - m_zmin + 1e-5);
	m_diagonal = sqrt(pow(m_xmax-m_xmin,2)+pow(m_ymax-m_ymin,2)+pow(m_zmax-m_zmin,2));

	double xsize = (double) (m_xmax - m_xmin) / gridSize;
	double ysize = (double) (m_ymax - m_ymin) / gridSize;
	double zsize = (double) (m_zmax - m_zmin) / gridSize;

	m_grid = TNT::Array3D<RayIntersectCell*>(gridSize, gridSize, gridSize);
	for (int x=0;x<gridSize;x++)
		for (int y=0;y<gridSize;y++)
			for (int z=0;z<gridSize;z++)
			{
				m_grid[x][y][z] = new RayIntersectCell();
				m_grid[x][y][z]->initGeometry(index2coordinates(x,y,z), xsize, ysize, zsize);
			}

	//step 2 - run over vertices and add their triangle to a grid
	for (Mesh::Vertex_const_iterator v = mesh.vertices_begin(); v != mesh.vertices_end(); v++)
	{
		//get the grid indexes for this vertex
		int xindex,yindex,zindex;
		Point_3 p = v->point();
		coordinates2index(p, xindex, yindex, zindex);
		
		Mesh::Halfedge_around_vertex_const_circulator c = v->vertex_begin();

		int checkCount = 0;
		if (c!= NULL) do {
			if (c->facet() != NULL) {
				int tindex = c->facet()->index();
				m_grid[xindex][yindex][zindex]->addTriangle(tindex);
			}
			if (++checkCount > 1000)
				throw QString("TooMuch!");
		} while (++c != v->vertex_begin());
	}

	//step 3 - if triangle vertices are in 3 different grids add to 4th grid

	//step 4
	m_initialized = true;
}

void RayIntersect::Init(
	Mesh& mesh, TNT::Array3D<RayIntersectCell*>& grid, triangleVector_t& triangleVector,
	double& xmin, double& xmax, double& ymin, double& ymax, double& zmin, double& zmax,
	double& xspread, double& yspread, double& zspread,
	const int gridSize, const bool debugMode)
{
	//step 0 - make sure all data structures are empty
	/*if (initialized)
	{
		for (int x = 0;x < grid.dim1(); x++)		
			for (int y = 0; y < grid.dim2(); y++)
				for (int z = 0; z < grid.dim3(); z++)
					delete grid[x][y][z];

		for (int i = 0; i < triangleVector.size(); i++)
			delete triangleVector[i];

		triangleVector.clear();
	}*/
	
	if (debugMode) 
		cout << "step 0 - make sure all data structures are empty done!" << endl;

	//step 1 - run over all triangles and add them to triangle vector
	triangleVector.reserve(mesh.size_of_facets());

	for (Mesh::Facet_const_iterator f = mesh.facets_begin(); f != mesh.facets_end(); f++)
	{
		assert(f->size == 3);

		Mesh::Halfedge_const_handle h1 = f->halfedge();
		Mesh::Halfedge_const_handle h2 = h1->next();
		Mesh::Halfedge_const_handle h3 = h2->next();

		Triangle_3_With_Normal* triangle = new Triangle_3_With_Normal(h1->vertex()->point(), h2->vertex()->point(), h3->vertex()->point(), f->normal());
		triangleVector.push_back(triangle);
	}
	
	mesh.computeBoundingBox();

	xmin = mesh.xmin(); xmax = mesh.xmax();
	ymin = mesh.ymin(); ymax = mesh.ymax();
	zmin = mesh.zmin(); zmax = mesh.zmax();
	xspread = 1 / (xmax - xmin + 1e-5);
	yspread = 1 / (ymax - ymin + 1e-5);
	zspread = 1 / (zmax - zmin + 1e-5);

	/*
	printf("xmin = %g, xmax = %g\n", xmin, xmax);
	printf("ymin = %g, ymax = %g\n", ymin, ymax);
	printf("zmin = %g, zmax = %g\n", zmin, zmax);
	printf("xspread = %g, yspread = %g, zspread = %g\n", xspread, yspread, zspread);
	*/

	double xsize = (double) (xmax - xmin) / gridSize;
	double ysize = (double) (ymax - ymin) / gridSize;
	double zsize = (double) (zmax - zmin) / gridSize;

	grid = TNT::Array3D<RayIntersectCell*>(gridSize, gridSize, gridSize);
	for (int x = 0; x < gridSize; x++)
		for (int y = 0; y < gridSize; y++)
			for (int z = 0; z < gridSize; z++)
			{
				grid[x][y][z] = new RayIntersectCell();
				grid[x][y][z]->initGeometry(
							index2coordinates(x,y,z,grid.dim1(),grid.dim2(),grid.dim3(), xmin, xmax, ymin, ymax, zmin, zmax), 
							xsize, ysize, zsize);
			}

	if (debugMode)
		cout << "step 1 - run over all triangles and add them to triangle vector done!" << endl;

	//step 2 - run over vertices and add their triangle to a grid
	for (Mesh::Vertex_const_iterator v = mesh.vertices_begin(); v != mesh.vertices_end(); v++)
	{
		//get the grid indexes for this vertex
		int xindex,yindex,zindex;
		Point_3 p = v->point();

		coordinates2index(
			p, 
			xindex, yindex, zindex,
			grid.dim1(), grid.dim2(), grid.dim3(),
			xmin, ymin, zmin,
			xspread, yspread, zspread);

		Mesh::Halfedge_around_vertex_const_circulator c = v->vertex_begin();

		int checkCount = 0;
		if (c!= NULL) do
		{
			if (c->facet() != NULL)
			{
				int tindex = c->facet()->index();
				grid[xindex][yindex][zindex]->addTriangle(tindex);
			}

			if (++checkCount > 1000)
			{
				throw QString("TooMuch!");
			}
		} while (++c != v->vertex_begin());
	}

	if (debugMode) 
		cout << "step 2 - run over vertices and add their triangle to a grid done!" << endl;

	//step 3 - if triangle vertices are in 3 different grids add to 4th grid
	if (debugMode) 
		cout << "step 3 - if triangle vertices are in 3 different grids add to 4th grid done!" << endl;

	//step 4
	m_initialized = true;

	//generateReportFile("e:/three_dimension/watermarking/tmp/log/rayintersect.html"); //E:\three_dimension\watermarking\tmp\log
}


void RayIntersect::Terminate(TNT::Array3D<RayIntersectCell*>& grid, triangleVector_t& triangleVector)
{
	//delete m_grid
	for (int x = 0; x < grid.dim1(); x++)
		for (int y = 0; y < grid.dim2(); y++)
			for (int z = 0; z < grid.dim3(); z++)
				delete grid[x][y][z];

	for (int i = 0; i < triangleVector.size(); i++)
		delete triangleVector[i];

	triangleVector.clear();
}

void RayIntersect::coordinates2index(const Point_3& p, int& x, int& y, int& z) const
{
	x = ((double) (p.x() - m_xmin) * m_xspread * m_grid.dim1());
	x = qMin(x,m_grid.dim1()-1);
	y = ((double) (p.y() - m_ymin) * m_yspread * m_grid.dim2());
	y = qMin(y,m_grid.dim2()-1);
	z = ((double) (p.z() - m_zmin) * m_zspread * m_grid.dim3());
	z = qMin(z,m_grid.dim3()-1);
}

void RayIntersect::coordinates2index(
	const Point_3& p, 
	int& x, int& y, int& z,
	const int& gridDim1, const int& gridDim2, const int& gridDim3, 
	const double& xmin, const double& ymin, const double& zmin,
	const double& xspread, const double& yspread, const double& zspread) const
{
	x = ((double) (p.x() - xmin) * xspread * gridDim1);
	x = qMin(x,gridDim1-1);
	y = ((double) (p.y() - ymin) * yspread * gridDim2);
	y = qMin(y,gridDim2-1);
	z = ((double) (p.z() - zmin) * zspread * gridDim3);
	z = qMin(z,gridDim3-1);
}

double RayIntersect::Query(
	const Ray_3& ray, /*const bool& initialized, */
	const TNT::Array3D<RayIntersectCell*>& grid, const triangleVector_t& triangleVector,
	const double& xmin, const double& ymin, const double& zmin,
	const double& xspread, const double& yspread, const double& zspread,
	const bool& removeSameDirectionIntersection, cellVector_t* debugCellVector, int* debugTriangleIndex, 
	FILE* debugFile, const bool debugMode) const
{
	/*assert(initialized);*/

	int xindex, yindex, zindex;

	if (debugFile) { fprintf(debugFile, "RayIntersect::Query begin\n"); fflush(debugFile); }
	if (debugMode) { cout << "RayIntersect::Query begin" << endl; }

	coordinates2index(
		ray.source(), 
		xindex, yindex, zindex,
		grid.dim1(), grid.dim2(), grid.dim3(),
		xmin, ymin, zmin,
		xspread, yspread, zspread);
	
	if (debugMode)
	{
		printf("grid.x= %d, grid.y= %d, grid.z= %d\n", grid.dim1(), grid.dim2(), grid.dim3());
		printf("ray source=(%g,%g,%g) ; index=(%d,%d,%d)\n", ray.source().x(), ray.source().y(), ray.source().z(), xindex, yindex, zindex);
	}

	if (debugFile) { fprintf(debugFile, "RayIntersect::Query converted coordinates 2 index\n"); fflush(debugFile); }
	if (debugMode) { cout << "RayIntersect::Query converted coordinates 2 index" << endl; }

	if (debugCellVector) debugCellVector->push_back(makeCellIndexVector(xindex,yindex,zindex));

	//look in this cell
	int triangleIndex;
	bool cellLegal = true;
	do {	
		if (debugMode) printf("going to find trignaleIndex... with ray source= (%g, %g, %g), ray direction= (%g, %g, %g)", ray.source().x(), ray.source().y(), ray.source().z(), ray.to_vector().x(), ray.to_vector().y(), ray.to_vector().z());

		triangleIndex = grid[xindex][yindex][zindex]->intersects(triangleVector, ray, debugMode);

		if (debugMode) printf(" Triangle index= %d\n", triangleIndex);

		if (debugFile) { fprintf(debugFile, "RayIntersect::Query triangle index %d\n",triangleIndex); fflush(debugFile); }
		if (debugMode) { printf("RayIntersect::Query triangle index %d\n", triangleIndex); }

		if (triangleIndex < 0)
		{
			//what's the new cell, each cell is neighbor with 26 other cells
			if (!findNextCell(ray, xindex, yindex, zindex, grid, debugFile))
			{
				cellLegal = false;
			}
			else
			{
				if (debugFile) { fprintf(debugFile, "RayIntersect::Query found next cell [%d %d %d]\n", xindex, yindex, zindex); fflush(debugFile); }
				if (debugMode) { cout << "RayIntersect::Query found next cell [" << xindex << " " << yindex << " " << zindex << "]" << endl; }

				//for debug purposes remember the cells we passed through
				if (debugCellVector) debugCellVector->push_back(makeCellIndexVector(xindex,yindex,zindex));
	
				if (xindex < 0 || yindex < 0 || zindex < 0 || xindex > grid.dim1() || yindex > grid.dim2() || zindex > grid.dim3())
				{
					cout << "Cell Index Out Of Range!!!" << endl;
					cellLegal = false;
				}
			}
		}
	} while (triangleIndex < 0 && cellLegal);

	if (debugFile) { fprintf(debugFile, "RayIntersect::Query found triangle"); fflush(debugFile); }
	if (debugMode) { cout << "RayIntersect::Query found triangle" << endl; }

	//triangle was found, now find distance to triangle (and check normal direction)
	if (debugTriangleIndex != NULL) *debugTriangleIndex = triangleIndex;

	if (triangleIndex >= 0 && triangleIndex < triangleVector.size())
	{
		Segment_3 iseg;

		Triangle_3_With_Normal* t = triangleVector[triangleIndex];
		CGAL::Object result = CGAL::intersection(t->supporting_plane(), ray);
		Point_3 ip;
		if (CGAL::assign(ip, result))
		{
			Segment_3 s(ray.source(), ip);
			double distance = sqrt(s.squared_length());

			if (removeSameDirectionIntersection)
			{
				if (t->normal() * ray.to_vector() < 0) 
					return FLT_MAX;
				else
					return distance;
			}
			else
				return distance;
		}
		else
		{
			/*
			if (CGAL::assign(iseg, result))
				cout << "###The segment intersection case" << endl;
			else
				cout << "###No intersection case" << endl;
			*/
		}
	}

	if (debugFile) { fprintf(debugFile, "RayIntersect::Query finished"); fflush(debugFile); }
	if (debugMode) { cout << "RayIntersect::Query finished" << endl; }

	return FLT_MAX;
}


bool RayIntersect::findNextCell(
	const Ray_3& ray, 
	int& x, int& y, int& z, 
	const TNT::Array3D<RayIntersectCell*>& grid,
	FILE* debugFile) const
{
	RayIntersectCell* cell = grid[x][y][z];

	Vector_3 rayVector = ray.to_vector();

	if (debugFile) { fprintf(debugFile, "RayIntersect::findNextCell 1\n"); fflush(debugFile); }

	bool use_pl1 = (cell->planes()[0].orthogonal_vector() * rayVector > 0 ? true : false); //otherwise pl5
	bool use_pl3 = (cell->planes()[2].orthogonal_vector() * rayVector > 0 ? true : false); //o/w pl6
	bool use_pl2 = (cell->planes()[1].orthogonal_vector() * rayVector > 0 ? true : false);// o/w pl4

	if (debugFile) { fprintf(debugFile, "RayIntersect::findNextCell 2\n"); fflush(debugFile); }	

	Plane_3 test_plane1 = (use_pl1 ? cell->planes()[0] : cell->planes()[4]);
	Plane_3 test_plane2 = (use_pl2 ? cell->planes()[1] : cell->planes()[3]);
	Plane_3 test_plane3 = (use_pl3 ? cell->planes()[2] : cell->planes()[5]);	
	
	if (debugFile)
	{
		fprintf(debugFile, "RayIntersect::findNextCell 3\n"); 
		fprintf(debugFile, "Test plane 1: %f %f %f %f\n", test_plane1.a(), test_plane1.b(), test_plane1.c(), test_plane1.d());
		fprintf(debugFile, "Test plane 2: %f %f %f %f\n", test_plane2.a(), test_plane2.b(), test_plane2.c(), test_plane2.d());
		fprintf(debugFile, "Test plane 3: %f %f %f %f\n", test_plane3.a(), test_plane3.b(), test_plane3.c(), test_plane3.d());
		fprintf(debugFile, "Ray: %f %f %f to direction %f %f %f\n", 
			ray.source().x(), ray.source().y(), ray.source().z(),
			ray.to_vector().x(), ray.to_vector().y(), ray.to_vector().z());
		fflush(debugFile);
	}	

	double id2Plane1 = ray2planeDistance(test_plane1, ray, debugFile);

	if (debugFile) { fprintf(debugFile, "RayIntersect::findNextCell 3.1\n"); fflush(debugFile); }

	double id2Plane2 = ray2planeDistance(test_plane2, ray, debugFile);

	if (debugFile) { fprintf(debugFile, "RayIntersect::findNextCell 3.2\n"); fflush(debugFile);	}

	double id2Plane3 = ray2planeDistance(test_plane3, ray, debugFile);

	if (debugFile) { fprintf(debugFile, "RayIntersect::findNextCell 3.3\n"); fflush(debugFile); }	

	if (id2Plane1 < id2Plane2 && id2Plane1 < id2Plane3)
	{
		if (use_pl1) 
			z++; 
		else 
			z--;
	}
	else if (id2Plane2 < id2Plane1 && id2Plane2 < id2Plane3)
	{
		if (use_pl2) 
			x--; 
		else 
			x++;
	}
	else
	{
		if (use_pl3) 
			y--; 
		else 
			y++;
	}

	if (debugFile) { fprintf(debugFile, "RayIntersect::findNextCell 5\n"); fflush(debugFile); }	

	if (x < 0 || y < 0 || z < 0 || x >= grid.dim1() || y >= grid.dim2() || z >= grid.dim3())
		return false;

	return true;
}

double RayIntersect::ray2planeDistance(const Plane_3& plane, const Ray_3& ray, FILE* debugFile) const
{
	double dist = FLT_MAX;

	if (debugFile) { fprintf(debugFile, "RayIntersect::ray2planeDistance begin\n"); fflush(debugFile); }

	CGAL::Object intersection_object = CGAL::intersection(plane, ray);
	Point_3 ip;

	if (debugFile) {fprintf(debugFile, "RayIntersect::ray2planeDistance found intersection\n"); fflush(debugFile);}

	if (CGAL::assign(ip, intersection_object))
	{
		if (debugFile) {fprintf(debugFile, "RayIntersect::ray2planeDistance assign successful\n"); fflush(debugFile);}

		Segment_3 s(ray.source(), ip);
		dist = s.squared_length();

		if (debugFile) {fprintf(debugFile, "RayIntersect::ray2planeDistance found distance %f\n", dist); fflush(debugFile);}
	}

	return dist;
}

Point_3 RayIntersect::index2coordinates(const int x, const int y, const int z) const
{
	Point_3 p(
		m_xmin + ((double) x / m_grid.dim1() * (m_xmax-m_xmin)),
		m_ymin + ((double) y / m_grid.dim2() * (m_ymax-m_ymin)),
		m_zmin + ((double) z / m_grid.dim3() * (m_zmax-m_zmin)));

	return p;
}

Point_3 RayIntersect::index2coordinates(
	const int x, const int y, const int z, 
	const int& gridDim1, const int& gridDim2, const int& gridDim3, 
	const double& xmin, const double& xmax, 
	const double& ymin, const double& ymax, 
	const double& zmin, const double& zmax) const
{
	Point_3 p(
		xmin + ((double) x / gridDim1 * (xmax-xmin)),
		ymin + ((double) y / gridDim2 * (ymax-ymin)),
		zmin + ((double) z / gridDim3 * (zmax-zmin)));
	return p;
}

void RayIntersect::renderSquare(const Point_3& p1, const Point_3& p2, const Point_3& p3, const Point_3& p4) const
{
	//render a cell nicely
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_QUADS);

	glVertex3f(p1.x(), p1.y(), p1.z());
	glVertex3f(p2.x(), p2.y(), p2.z());
	glVertex3f(p3.x(), p3.y(), p3.z());
	glVertex3f(p4.x(), p4.y(), p4.z());

	glEnd();
}

void RayIntersect::renderCell(const std::vector<int>& cellIndexes) const
{
	renderCell(cellIndexes[0], cellIndexes[1], cellIndexes[2]);
}

void RayIntersect::renderCell(const int x, const int y, const int z) const
{
	RayIntersectCell* cell = m_grid[x][y][z];

	renderSquare(cell->points()[7],cell->points()[3],cell->points()[1],cell->points()[5]);
	renderSquare(cell->points()[1],cell->points()[3],cell->points()[2],cell->points()[0]);
	renderSquare(cell->points()[5],cell->points()[1],cell->points()[0],cell->points()[4]);
	renderSquare(cell->points()[7],cell->points()[5],cell->points()[4],cell->points()[6]);
	renderSquare(cell->points()[4],cell->points()[0],cell->points()[2],cell->points()[6]);
	renderSquare(cell->points()[6],cell->points()[2],cell->points()[3],cell->points()[7]);
}

void RayIntersect::renderAllCells() const
{
	for (int x = 0; x < m_grid.dim1(); x++)
		for (int y = 0; y < m_grid.dim2(); y++)
			for (int z = 0; z < m_grid.dim3(); z++) 
			{
				renderCell(x,y,z);
			}
}

/*
void RayIntersect::renderRay(const Ray_3& ray) const
{
	Point_3 p = ray.point(m_diagonal);
	CylLine cylLine(m_diagonal/1000, 24);
	cylLine(
		ray.source().x(),ray.source().y(),ray.source().z(),
		p.x(),p.y(),p.z());	
}

void RayIntersect::renderRayPath(const Ray_3& ray) const
{
	//query the ray and go over all the cells
	cellVector_t debugCellVector;
	int debugTriangleIndex;
	double dist = Query(ray, &debugCellVector, &debugTriangleIndex);
	
	if (debugTriangleIndex>=0 && debugTriangleIndex < m_triangleVector.size()) {
		glColor4d(0.1,0.8,0.1,1.0);
		glBegin(GL_TRIANGLES);
		Triangle_3_With_Normal* t = m_triangleVector[debugTriangleIndex];
		glVertex3f(t->vertex(0).x(),t->vertex(0).y(),t->vertex(0).z());
		glVertex3f(t->vertex(1).x(),t->vertex(1).y(),t->vertex(1).z());
		glVertex3f(t->vertex(2).x(),t->vertex(2).y(),t->vertex(2).z());
		glEnd();
	}
	
	glEnable(GL_LIGHTING);
	glColor4d(0.1,0.1,0.8,1.0);
	renderRay(ray);

	glDepthMask(GL_FALSE);
    glDisable(GL_LIGHTING);
	for (int i=0;i<debugCellVector.size();i++) {
		float perc = (float) (i+1) / debugCellVector.size();
		glColor4d(0.8 * perc, 0.1, 0.1, 0.4);
		renderCell(debugCellVector[i]);
	}

	glDepthMask(GL_TRUE);
	
}
*/

/*
void RayIntersect::generateReportFile(const char* fileName)
{
	QFile file( fileName);
	if ( file.open( QIODevice::WriteOnly ) ) {
		Q3TextStream stream( &file );

		stream << "<html><head/><body>\n";
		stream << "<h1>RayIntersect Report File " << QDateTime::currentDateTime().toString("dd/MM/yyyy hh:mm:ss") << "</h1>\n";
		stream << "Dimensions: " << m_grid.dim1() << "x" << m_grid.dim2() << "x" << m_grid.dim3() << "<br>";
		stream << "Mesh Info: #Triangles=" << m_triangleVector.size() << "<br>\n";

		stream << "<table>";
		stream << "<tr><th>X</th><th>Y</th><th>Z</th><th>Size</th></tr>";
		for (int x=0;x<m_grid.dim1();x++)
			for (int y=0;y<m_grid.dim2();y++)
				for (int z=0;z<m_grid.dim3();z++) 
				{
					if (m_grid[x][y][z]->size()) {
						stream << "<tr><td>" << x << "</td><td>" << y << "</td><td>" << z << "</td><td>";
						stream << m_grid[x][y][z]->size() << "</td></tr>";
					}
				}
		stream << "</table>\n";
		stream << "</body></html>";

		file.close();
	}
}
*/